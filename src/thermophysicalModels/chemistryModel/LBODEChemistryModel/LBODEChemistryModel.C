/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "LBODEChemistryModel.H"
#include "chemistrySolver.H"
#include "reactingMixture.H"
#include "multiComponentMixture.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::LBODEChemistryModel<CompType, ThermoType>::LBODEChemistryModel
(
    const fvMesh& mesh,
    const objectRegistry& obj,
    const word& compTypeName,
    const word& thermoTypeName
)
:
    ODEChemistryModel<CompType, ThermoType>(mesh, obj, compTypeName, thermoTypeName),
    balancer_
    (
        IOdictionary
        (
            IOobject
            (
                "chemistryProperties",
                mesh.time().constant(),
                mesh.db(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        )
    ), 
    mapper_
    (
        IOdictionary
        (
            IOobject
            (
                "chemistryProperties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        ),
        this->thermo().composition()
    ),
    cpuTimes_
    (
        IOobject
        (
            "cellCpuTimes",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        scalar(0.0)
    ),
    refMap_
    (
        IOobject
        (
            "referenceMap",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        scalar(0.0)
    ),
    cpuSolveFile_(),
    Treact_
    (
        IOdictionary
        (
            IOobject
            (
                "chemistryProperties",
                mesh.time().constant(),
                mesh.db(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        ).lookupOrDefault<scalar>("Treact", 0)
    )
{
    if(balancer_.log())
    {
        cpuSolveFile_ = logFile("cpu_solve.out");
        cpuSolveFile_() << "         time" << tab
                        << "  getProblems" << tab  
                        << "  updateState" << tab
                        << "      balance" << tab
                        << "  solveBuffer" << tab
                        << "    unbalance" << tab
                        << "      rank ID" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::LBODEChemistryModel<CompType, ThermoType>::~LBODEChemistryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class CompType, class ThermoType>
template <class DeltaTType>
Foam::scalar Foam::LBODEChemistryModel<CompType, ThermoType>::solve
(
    const DeltaTType& deltaT
)
{
    // CPU time analysis
    clockTime timer;
    scalar t_getProblems(0);
    scalar t_updateState(0);
    scalar t_balance(0);
    scalar t_solveBuffer(0);
    scalar t_unbalance(0);

    // Does nothing in OpenFOAM 8; so skipped here
    //BasicChemistryModel<CompType>::correct();

    if(!this->chemistry_)
    {
        return VGREAT;
    }

    timer.timeIncrement();
    DynamicList<ChemistryProblem> allProblems = getProblems(deltaT);
    t_getProblems = timer.timeIncrement();

    RecvBuffer<ChemistrySolution> incomingSolutions;

    if(balancer_.active())
    {
        timer.timeIncrement();
        balancer_.updateState(allProblems);
        t_updateState = timer.timeIncrement();

        timer.timeIncrement();
        auto guestProblems = balancer_.balance(allProblems);
        auto ownProblems = balancer_.getRemaining(allProblems);
        t_balance = timer.timeIncrement();

        timer.timeIncrement();
        auto ownSolutions = solveList(ownProblems);
        auto guestSolutions = solveBuffer(guestProblems);
        t_solveBuffer = timer.timeIncrement();

        timer.timeIncrement();      
        incomingSolutions = balancer_.unbalance(guestSolutions);
        incomingSolutions.append(ownSolutions);
        t_unbalance = timer.timeIncrement();
    }
    else
    {
        timer.timeIncrement();
        incomingSolutions.append(solveList(allProblems));
        t_solveBuffer = timer.timeIncrement();
    }
        
    if(balancer_.log())
    {
        if(balancer_.active())
        {
            balancer_.printState();
        }
        cpuSolveFile_() << setw(13)
                        << this->time().timeOutputValue()<<tab
                        << setw(13) << t_getProblems<<tab
                        << setw(13) << t_updateState<<tab
                        << setw(13) << t_balance<<tab
                        << setw(13) << t_solveBuffer<<tab
                        << setw(13) << t_unbalance<<tab
                        << setw(13) << Pstream::myProcNo()
                        << endl;
    }

    return updateReactionRates(incomingSolutions);
}


template <class CompType, class ThermoType>
void Foam::LBODEChemistryModel<CompType, ThermoType>::solveSingle
(
    ChemistryProblem& problem, ChemistrySolution& solution
) const
{
    scalar timeLeft = problem.deltaT;
    scalarField c0 = problem.c;

    // Timer begins
    clockTime time;
    time.timeIncrement();

    // Define a const label to pass as the cell index placeholder
    //const label arbitrary = 0;

    // Calculate the chemical source terms
    while(timeLeft > VSMALL)
    {
        scalar dt = timeLeft;
        //this->solve(
        //    problem.pi,
        //    problem.Ti,
        //    problem.c,
        //    arbitrary,
        //    dt,
        //    problem.deltaTChem);
        this->solver().solve
        (
            problem.c,
            problem.Ti,
            problem.pi,
            problem.deltaTChem,
            dt
        );
        timeLeft -= dt;
    }

    solution.c_increment = (problem.c - c0) / problem.deltaT;
    solution.deltaTChem = min(problem.deltaTChem, this->deltaTChemMax_);

    // Timer ends
    solution.cpuTime = time.timeIncrement();

    solution.cellid = problem.cellid;
    solution.rhoi = problem.rhoi;
}


template <class CompType, class ThermoType>
Foam::scalar
Foam::LBODEChemistryModel<CompType, ThermoType>::updateReactionRates
(
    const RecvBuffer<ChemistrySolution>& solutions
)
{
    scalar deltaTMin = VGREAT;

    for(const auto& array : solutions)
    {
        for(const auto& solution : array)
        {

            updateReactionRate(solution, solution.cellid);

            deltaTMin = min(solution.deltaTChem, deltaTMin);
        
            cpuTimes_[solution.cellid] = solution.cpuTime;
        }
    }

    return deltaTMin;
}


template <class CompType, class ThermoType>
Foam::scalar Foam::LBODEChemistryModel<CompType, ThermoType>::solve
(
    const scalarField& deltaT
)
{
    return this->solve<scalarField>(deltaT);
}


template <class CompType, class ThermoType>
Foam::scalar Foam::LBODEChemistryModel<CompType, ThermoType>::solve
(
    const scalar deltaT
)
{
    // Don't allow the time-step to change more than a factor of 2
    return min
    (
        this->solve<UniformField<scalar>>(UniformField<scalar>(deltaT)),
        2 * deltaT
    );
}


template <class CompType, class ThermoType>
Foam::RecvBuffer<Foam::ChemistrySolution>
Foam::LBODEChemistryModel<CompType, ThermoType>::solveBuffer
(
    RecvBuffer<ChemistryProblem>& problems
) const
{
    // allocate the solutions buffer
    RecvBuffer<ChemistrySolution> solutions;
    
    for(auto& p : problems)
    {
        solutions.append(solveList(p));
    }
    return solutions;
}


template <class CompType, class ThermoType>
Foam::DynamicList<Foam::ChemistrySolution>
Foam::LBODEChemistryModel<CompType, ThermoType>::solveList
(
    UList<ChemistryProblem>& problems
) const
{
    List<ChemistrySolution> tsols
    (
        problems.size(),
        ChemistrySolution(this->nSpecie_)
    );
    DynamicList<ChemistrySolution> solutions (tsols.xfer());

    for(label i = 0; i < problems.size(); ++i)
    {
        solveSingle(problems[i], solutions[i]);
    }
    return solutions;
}



template <class CompType, class ThermoType>
template<class DeltaTType>
Foam::DynamicList<Foam::ChemistryProblem>
Foam::LBODEChemistryModel<CompType, ThermoType>::getProblems
(
    const DeltaTType& deltaT
)
{
    const scalarField& T = this->thermo().T();
    const scalarField& p = this->thermo().p();
    tmp<volScalarField> trho(this->thermo().rho());
    const scalarField& rho = trho();
    


    DynamicList<ChemistryProblem> solved_problems;
    DynamicList<ChemistryProblem> mapped_problems;

    solved_problems.resize(p.size(), ChemistryProblem(this->nSpecie_));

    scalarField massFraction(this->nSpecie_);
    scalarField concentration(this->nSpecie_);

    label counter = 0;

    forAll(T, celli)
    {

        if(T[celli] > Treact_)
        {
            for(label i = 0; i < this->nSpecie_; i++)
            {
                concentration[i] = rho[celli] * this->Y_[i][celli] / this->specieThermo_[i].W();
                massFraction[i] = this->Y_[i][celli];
            }
            
            ChemistryProblem problem;
            problem.c = getVariable(concentration, massFraction);        
            problem.Ti = T[celli];
            problem.pi = p[celli];
            problem.rhoi = rho[celli];
            problem.deltaTChem = this->deltaTChem_[celli];
            problem.deltaT = deltaT[celli];
            problem.cpuTime = cpuTimes_[celli];
            problem.cellid = celli;

            // NOTE: Solve everything, side-effects of not mapping reference solutions?
            // This check can only be done based on the concentration as the 
            // reference temperature is not known
            if (mapper_.shouldMap(massFraction))
            {                
                mapped_problems.append(problem);
                refMap_[celli] = 1;
            } else {
                solved_problems[counter] = problem;
                counter++;
                refMap_[celli] = 2;
            }

        }
        else
        {
            for(label i = 0; i < this->nSpecie(); i++)
            {
                this->RR_[i][celli] = 0;
            }
        }

    }

    //the real size is set here
    solved_problems.setSize(counter);

    runtime_assert(solved_problems.size() + mapped_problems.size() == p.size(), "getProblems fails");

    // No mapped problems, solve everything
    this->map(mapped_problems, solved_problems);
    

    return solved_problems;
}

template <class CompType, class ThermoType>
void Foam::LBODEChemistryModel<CompType, ThermoType>::updateReactionRate
(
    const ChemistrySolution& solution, const label& i
)
{
    for(label j = 0; j < this->nSpecie_; j++)
    {
        this->RR_[j][i] = solution.c_increment[j] * this->specieThermo_[j].W();
    }
    this->deltaTChem_[i] = min(solution.deltaTChem, this->deltaTChemMax_);
}

template <class CompType, class ThermoType>
Foam::scalarField Foam::LBODEChemistryModel<CompType, ThermoType>::getVariable
(
    const scalarField& concentration, const scalarField& massFraction
)
{
    return concentration;
}

template <class CompType, class ThermoType>
void Foam::LBODEChemistryModel<CompType, ThermoType>::map
(
    DynamicList<ChemistryProblem>& mapped_problems, 
    DynamicList<ChemistryProblem>& solved_problems
)
{
    if (mapped_problems.size() > 0)
    {

        ChemistryProblem refProblem = mapped_problems[0];
        scalar refTemperature = refProblem.Ti;

        ChemistrySolution refSolution(this->nSpecie_);
        solveSingle(refProblem, refSolution);
        refMap_[refProblem.cellid] = 0;

        for (auto& problem : mapped_problems){

            // Check that the refmap temperature condition is also fullfilled
            if (mapper_.temperatureWithinRange(problem.Ti, refTemperature))
            {
                updateReactionRate(refSolution, problem.cellid);
                cpuTimes_[problem.cellid] = refSolution.cpuTime;
            }
            // Otherwise solve
            else 
            {
                solved_problems.append(problem);
                refMap_[problem.cellid] = 2;
            }
        }


    }
}
// ************************************************************************* //

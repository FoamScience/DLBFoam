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

Class
    Foam::LBODEChemistryModel

Description
    A port of LoadBalancedChemistryModel to FE-4 with Load-Balaning
    functionality and no reference solution mapping.

Author
    Mohammed Elwardi Fadeli (elwardifadeli@gmail.com)

SourceFiles
    LBODEChemistryModelI.H
    LBODEChemistryModel.C

Note
    Depends on ODEChemistryModel from
    https://github.com/FoamScience/fe4-thermodynamics
    and not on the one from Foam-Extend 4. You have been warned

\*---------------------------------------------------------------------------*/

#ifndef LBODEChemistryModel_H
#define LBODEChemistryModel_H

#include "ODEChemistryModel.H"
#include "LoadBalancer.H"
#include "ChemistryProblem.H"
#include "ChemistrySolution.H"
#include "LoadBalancer.H"
#include "OFstream.H"
#include "IOmanip.H"
#include "UniformField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class LBODEChemistryModel Declaration
\*---------------------------------------------------------------------------*/

template<class CompType, class ThermoType>
class LBODEChemistryModel
:
    public ODEChemistryModel<CompType, ThermoType>
{
protected:

    // Protected data

        // Load balancing object
        LoadBalancer balancer_;

        // Field containing chemistry CPU time information    
        volScalarField cpuTimes_;

        // A file to output the balancing stats
        autoPtr<OFstream>        cpuSolveFile_;

    // Protected Member Functions

        //- Disallow default bitwise assignment
        void operator=(const LBODEChemistryModel&);

        //- Get the list of problems to be solved
        template<class DeltaTType>
        DynamicList<ChemistryProblem> getProblems(const DeltaTType& deltaT);

        //- Solve a list of chemistry problems and return a list of solutions
        DynamicList<ChemistrySolution> 
        solveList(UList<ChemistryProblem>& problems) const;

        //- Solve the problem buffer coming from the balancer
        RecvBuffer<ChemistrySolution> 
        solveBuffer(RecvBuffer<ChemistryProblem>& problems) const;

        //- Update the reaction rate of cell i
        virtual void
        updateReactionRate(const ChemistrySolution& solution, const label& i);

        //- Update the reaction rates from a list of solutions
        scalar updateReactionRates(const RecvBuffer<ChemistrySolution>& solutions);

        //- Return the main chemistry variable
        //  (concentration for standard, mass fraction for pyJac)
        virtual scalarField getVariable(const scalarField& concentration, const scalarField& massFraction);

        //- Solve the reaction system for the given time step
        //  of given type and return the characteristic time
        template<class DeltaTType>
        scalar solve(const DeltaTType& deltaT);

        //- Solve the reaction system for the given start time and time
        //  step and return the characteristic time
        virtual scalar solve(const scalar t0, const scalar deltaT) override
        {
            // NOTE: Is there a side-effect of discarding t0?
            return this->solve<UniformField<scalar>>(UniformField<scalar>(deltaT));
        }

public:

    //- Runtime type information
    TypeName("LBODEChemistryModel");


    // Constructors

        //- Construct from components
        LBODEChemistryModel
        (
            const fvMesh& mesh,
            const objectRegistry& obj,
            const word& compTypeName,
            const word& thermoTypeName
        );


    //- Destructor
    virtual ~LBODEChemistryModel();


    // Member Functions

        //- Solve a single ChemistryProblem and put the solution to ChemistrySolution 
        void solveSingle(ChemistryProblem& problem, ChemistrySolution& solution) const;

        // Chemistry model functions (overriding functions in
        // StandardChemistryModel to use the private solve function)

            //- Solve the reaction system for the given time step
            //  and return the characteristic time
            virtual scalar solve(const scalar deltaT);

            //- Solve the reaction system for the given time step
            //  and return the characteristic time
            virtual scalar solve(const scalarField& deltaT);

        //- Create and return a log file of the given name
        inline autoPtr<OFstream> logFile(const word& name) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "LBODEChemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "LBODEChemistryModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

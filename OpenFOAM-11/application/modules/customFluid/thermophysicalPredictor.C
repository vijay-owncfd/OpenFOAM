/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "customFluid.H"
#include "fvcDdt.H"
#include "fvmDiv.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::customFluid::thermophysicalPredictor()
{
    volScalarField& he = thermo_.he();

/*
    // Keeping it here, just for reference
    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + fvm::div(phi, he)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + pressureWork
        (
            he.name() == "e"
          ? fvc::div(phi, p/rho)()
          : -dpdt
        )
      + thermophysicalTransport->divq(he)
     ==
        (
            buoyancy.valid()
          ? fvModels().source(rho, he) + rho*(U & buoyancy->g)
          : fvModels().source(rho, he)
        )
    );
*/

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + fvm::div(phi, he, "div(phi,he)")
      ==
      - thermophysicalTransport->divq(he)
      + fvModels().source(rho, he)
    );

    if(energySolverMode_ == "Static")
    {
        if(he.name() == "h")
        {
            EEqn -= (dpdt + (fvc::grad(p) & U)());
        } 
        else
        {
            EEqn -= (-p * fvc::div(phi, 1/rho)());     // In future, calculate this in the pressure equation
        }  
    }
    else
    {
        EEqn += (fvc::ddt(rho, K) + fvc::div(phi, K));

        if(buoyancy.valid())
        {
            EEqn -= rho*(U & buoyancy->g);
        }
   
        if(he.name() == "h")
        {
            EEqn -= dpdt;
        }
        else
        {
            EEqn -= -fvc::div(phi, p/rho)();
        }
    }

    EEqn.relax();

    fvConstraints().constrain(EEqn);

    EEqn.solve();

    fvConstraints().constrain(he);

    thermo_.correct();
}


// ************************************************************************* //

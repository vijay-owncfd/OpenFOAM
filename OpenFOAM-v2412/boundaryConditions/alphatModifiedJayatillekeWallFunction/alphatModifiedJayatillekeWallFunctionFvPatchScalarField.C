/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2022 OpenCFD Ltd
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

#include "alphatModifiedJayatillekeWallFunctionFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"
#include "nutWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

scalar alphatModifiedJayatillekeWallFunctionFvPatchScalarField::maxExp_ = 50.0;
scalar alphatModifiedJayatillekeWallFunctionFvPatchScalarField::tolerance_ = 0.01;
label alphatModifiedJayatillekeWallFunctionFvPatchScalarField::maxIters_ = 10;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void alphatModifiedJayatillekeWallFunctionFvPatchScalarField::checkType()
{
    if (!isA<wallFvPatch>(patch()))
    {
        FatalErrorInFunction
            << "Patch type for patch " << patch().name() << " must be wall\n"
            << "Current patch type is " << patch().type() << nl
            << exit(FatalError);
    }
}


tmp<scalarField> alphatModifiedJayatillekeWallFunctionFvPatchScalarField::yPlus
(
    const compressible::turbulenceModel& turbModel
) const
{
    const label patchi = patch().index();

    const tmp<volScalarField> tnut = turbModel.nut();
    const volScalarField& nut = tnut();

    if (isA<nutWallFunctionFvPatchScalarField>(nut.boundaryField()[patchi]))
    {
        const auto& nutPf =
            dynamic_cast<const nutWallFunctionFvPatchScalarField&>
            (
                nut.boundaryField()[patchi]
            );

        return nutPf.yPlus();
    }
    else
    {
        const scalarField& y = turbModel.y()[patchi];
        const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];

        return
            y*sqrt(turbModel.nuEff(patchi)*mag(Uw.snGrad()))
           /turbModel.nu(patchi);
    }
}


scalar alphatModifiedJayatillekeWallFunctionFvPatchScalarField::Psmooth
(
    const scalar Prat
) const
{
    return 9.24*(pow(Prat, 0.75) - 1.0)*(1.0 + 0.28*exp(-0.007*Prat));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

alphatModifiedJayatillekeWallFunctionFvPatchScalarField::
alphatModifiedJayatillekeWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    Prt_(0.85),
    kappa_(0.41),
    E_(9.8)
{
    checkType();
}


alphatModifiedJayatillekeWallFunctionFvPatchScalarField::
alphatModifiedJayatillekeWallFunctionFvPatchScalarField
(
    const alphatModifiedJayatillekeWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    Prt_(ptf.Prt_),
    kappa_(ptf.kappa_),
    E_(ptf.E_)
{}


alphatModifiedJayatillekeWallFunctionFvPatchScalarField::
alphatModifiedJayatillekeWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    Prt_(dict.getOrDefault<scalar>("Prt", 0.85)),
    kappa_(dict.getOrDefault<scalar>("kappa", 0.41)),
    E_(dict.getOrDefault<scalar>("E", 9.8))
{
    checkType();
}


alphatModifiedJayatillekeWallFunctionFvPatchScalarField::
alphatModifiedJayatillekeWallFunctionFvPatchScalarField
(
    const alphatModifiedJayatillekeWallFunctionFvPatchScalarField& awfpsf
)
:
    fixedValueFvPatchScalarField(awfpsf),
    Prt_(awfpsf.Prt_),
    kappa_(awfpsf.kappa_),
    E_(awfpsf.E_)
{
    checkType();
}


alphatModifiedJayatillekeWallFunctionFvPatchScalarField::
alphatModifiedJayatillekeWallFunctionFvPatchScalarField
(
    const alphatModifiedJayatillekeWallFunctionFvPatchScalarField& awfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(awfpsf, iF),
    Prt_(awfpsf.Prt_),
    kappa_(awfpsf.kappa_),
    E_(awfpsf.E_)
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void alphatModifiedJayatillekeWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const label patchi = patch().index();

    // Retrieve turbulence properties from model
    const auto& turbModel = db().lookupObject<compressible::turbulenceModel>
    (
        IOobject::groupName
        (
            compressible::turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    const scalarField yPlusp(yPlus(turbModel));

    const scalarField& y = turbModel.y()[patchi];

    const tmp<scalarField> talphaw = turbModel.alpha(patchi);
    const scalarField& alphaw = talphaw();

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const scalarField& rhow = turbModel.rho().boundaryField()[patchi];

    scalarField& alphatw = *this;
    
    // Populate boundary values
    forAll(alphatw, facei)
    {
        // max limiter to prevent floating point exceptions during Tplus computations below
        const scalar yPlus = max(scalar(0.1),yPlusp[facei]);    
        const scalar uTau = yPlus/y[facei]*nuw[facei];

        // Molecular Prandtl number
        const scalar Pr = rhow[facei]*nuw[facei]/alphaw[facei];

        // Molecular-to-turbulent Prandtl number ratio
        const scalar Prat = Pr/Prt_;

        // Thermal sublayer thickness
        const scalar P = Psmooth(Prat);

        // Temperature profile      
        const scalar TplusVis  = Pr*yPlus;
        const scalar TplusLog = Prt_ * (log(E_*yPlus)/kappa_ + P);
           
        // Blended Tplus with fourth-power Harmonic averaging
        const scalar Tplus = pow025(1/(1/pow(TplusVis,4) + 1/pow(TplusLog,4)));
        
        // Evaluate new effective thermal diffusivity
        const scalar alphaEff = rhow[facei]*uTau*y[facei]/(Tplus + VSMALL);

        // Update turbulent thermal diffusivity
        alphatw[facei] = max(scalar(0), alphaEff - alphaw[facei]);

        if (debug)
        {
            Info<< "    uTau           = " << uTau << nl
                << "    Pr             = " << Pr << nl
                << "    Prt            = " << Prt_ << nl
                << "    yPlus          = " << yPlus << nl
                << "    alphaEff       = " << alphaEff << nl
                << "    alphaw         = " << alphaw[facei] << nl
                << "    alphatw        = " << alphatw[facei] << nl
                << endl;
        }
    }

    fixedValueFvPatchField<scalar>::updateCoeffs();
}


void alphatModifiedJayatillekeWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeEntryIfDifferent<scalar>("Prt", 0.85, Prt_);
    os.writeEntryIfDifferent<scalar>("kappa", 0.41, kappa_);
    os.writeEntryIfDifferent<scalar>("E", 9.8, E_);
    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    alphatModifiedJayatillekeWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //

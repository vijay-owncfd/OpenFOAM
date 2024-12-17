/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "nutUKaderWallFunctionFvPatchScalarField.H"
#include "momentumTransportModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalarField> nutUKaderWallFunctionFvPatchScalarField::nut() const
{
    const label patchi = patch().index();

    const momentumTransportModel& turbModel =
        db().lookupType<momentumTransportModel>(internalField().group());

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const scalarField magGradU(mag(Uw.snGrad()));
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    return max
    (
        scalar(0),
        sqr(calcUTau(magGradU))/(magGradU + rootVSmall) - nuw
    );
}


tmp<scalarField> nutUKaderWallFunctionFvPatchScalarField::calcUTau
(
    const scalarField& magGradU
) const
{
    const label patchi = patch().index();

    const momentumTransportModel& turbModel =
        db().lookupType<momentumTransportModel>(internalField().group());

    const scalarField& y = turbModel.y()[patchi];

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const vectorField Urel(Uw.patchInternalField() - Uw);
    const vectorField& nHat = patch().n();

    const vectorField Un((Urel & nHat) * nHat);
    const vectorField Ut(Urel - Un);
    const scalarField magUp(mag(Ut));

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const scalarField& nutw = *this;

    tmp<scalarField> tuTau(new scalarField(patch().size(), 0.0));
    scalarField& uTau = tuTau.ref();

    forAll(uTau, facei)
    {
        scalar ut_n = sqrt((nutw[facei] + nuw[facei])*magGradU[facei]);

        scalar ut_nm1 = magUp[facei]/25.0;  // one guess at u+ = 25
        scalar f_nm1 = funcOfUTau(ut_nm1);
        scalar f_n = 0.0;

        if (ut > rootVSmall)
        {
            int iter = 0;
            scalar err = great;

            do
            {
                f_n = funcOfUTau(ut_n);

                scalar uTauNew = ut_n - f_n * (ut_n - ut_nm1)/(f_n - f_nm1 + 1e-10);

                err = mag((ut_n - uTauNew)/ut_n);
                ut_nm1 = ut_n;
                f_nm1 = f_n;
                ut_n = uTauNew;

            } while (ut_n > rootVSmall && err > 0.01 && ++iter < 10);

            uTau[facei] = max(0.0, ut_n);
        }
    }

    return tuTau;
}

scalar nutUKaderWallFunctionFvPatchScalarField::funcOfUTau
(
    const scalar uTau,
    const scalar y,
    const scalar nu,
    const scalar rho,
    const scalar dpdx
)
{
    scalar u0_Plus = 10;
    scalar y0_Plus =6;

    scalar yPlus_S = 60;
    scalar yPlus = max(y*uTau/nu, 0.01);
    scalar gamma = -0.01 * pow(yPlus, 4) / (1.0 + 5.0 * yPlus);

    scalar alpha = nu * dpdx/max(rho * pow(uTau, 3), 1e-10);
    
    scalar uPlusLam = yPlus * (1.0 + alpha*yPlus/2.0);

    scalar uPlusTurb = 10.0;    // just initialization

    scalar S0 = 1.0 + alpha * y0_Plus;
    scalar Sdash = 1.0;

    if(yPlus < yPlus_S)
    {
        Sdash += alpha*yPlus;
    }
    else
    {
        Sdash += alpha*yPlus_S;
    }
    
    if(pos(Sdash) && pos(S0))
    {
        if(yPlus < yPlus_S)
        {
            Sdash = sqrt(Sdash);
            S0 = sqrt(S0);

            uPlusTurb = u0_Plus + log((mag(Sdash - 1) * (S0 + 1)) / ((Sdash + 1) * mag(S0 - 1))) / kappa_
                        + (Sdash - S0) * 2 / kappa_;
        }
        else
        {
            uPlusTurb = u0_Plus + sqrt(Sdash) * log(yPlus/y0Plus) / kappa_;
        }
    }
    else
    {
        uPlusTurb = log(E_ * yPlus)/kappa_;
    }
    
    scalar uPlus = exp(gamma) * uPlusLam + exp(1.0/gamma) * uPlusTurb; 
    scalar ff = U/uTau - uPlus;
    return ff;
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nutUKaderWallFunctionFvPatchScalarField::
nutUKaderWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nutWallFunctionFvPatchScalarField(p, iF, dict),
    compressible_(false)
{}


nutUKaderWallFunctionFvPatchScalarField::
nutUKaderWallFunctionFvPatchScalarField
(
    const nutUKaderWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    nutWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    compressible_(ptf.compressible_)
{}


nutUKaderWallFunctionFvPatchScalarField::
nutUKaderWallFunctionFvPatchScalarField
(
    const nutUKaderWallFunctionFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(wfpsf, iF),
    compressible_(wfpsf.compressible_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField> nutUKaderWallFunctionFvPatchScalarField::yPlus() const
{
    const label patchi = patch().index();

    const momentumTransportModel& turbModel =
        db().lookupType<momentumTransportModel>(internalField().group());

    const scalarField& y = turbModel.y()[patchi];
    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    return y*calcUTau(mag(Uw.snGrad()))/nuw;
}


void nutUKaderWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    nutUKaderWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

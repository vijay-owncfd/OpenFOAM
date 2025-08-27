/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2019 OpenFOAM Foundation
    Copyright (C) 2017-2024 OpenCFD Ltd.
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

#include "epsilonTwoLayerWallFunctionFvPatchScalarField.H"
#include "nutWallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fvMatrix.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::scalar Foam::epsilonTwoLayerWallFunctionFvPatchScalarField::tolerance_ = 1e-5;

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::epsilonTwoLayerWallFunctionFvPatchScalarField::setMaster()
{
    if (master_ != -1)
    {
        return;
    }

    const auto& epsilon =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = epsilon.boundaryField();

    label master = -1;
    forAll(bf, patchi)
    {
        if (isA<epsilonTwoLayerWallFunctionFvPatchScalarField>(bf[patchi]))
        {
            epsilonTwoLayerWallFunctionFvPatchScalarField& epf = epsilonPatch(patchi);

            if (master == -1)
            {
                master = patchi;
            }

            epf.master() = master;
        }
    }
}


void Foam::epsilonTwoLayerWallFunctionFvPatchScalarField::createAveragingWeights()
{
    const auto& epsilon =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = epsilon.boundaryField();

    const fvMesh& mesh = epsilon.mesh();

    if (initialised_ && !mesh.changing())
    {
        return;
    }

    volScalarField weights
    (
        IOobject
        (
            "weights",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh,
        dimensionedScalar(dimless, Zero)
    );

    DynamicList<label> epsilonPatches(bf.size());
    forAll(bf, patchi)
    {
        if (isA<epsilonTwoLayerWallFunctionFvPatchScalarField>(bf[patchi]))
        {
            epsilonPatches.append(patchi);

            const labelUList& faceCells = bf[patchi].patch().faceCells();
            for (const auto& faceCell : faceCells)
            {
                ++weights[faceCell];
            }
        }
    }

    cornerWeights_.setSize(bf.size());

    for (const auto& patchi : epsilonPatches)
    {
        const fvPatchScalarField& wf = weights.boundaryField()[patchi];
        cornerWeights_[patchi] = 1.0/wf.patchInternalField();
    }

    G_.setSize(internalField().size(), Zero);
    epsilon_.setSize(internalField().size(), Zero);

    initialised_ = true;
}


Foam::epsilonTwoLayerWallFunctionFvPatchScalarField&
Foam::epsilonTwoLayerWallFunctionFvPatchScalarField::epsilonPatch
(
    const label patchi
)
{
    const auto& epsilon =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = epsilon.boundaryField();

    const auto& epf =
        refCast<const epsilonTwoLayerWallFunctionFvPatchScalarField>(bf[patchi]);

    return const_cast<epsilonTwoLayerWallFunctionFvPatchScalarField&>(epf);
}


void Foam::epsilonTwoLayerWallFunctionFvPatchScalarField::calculateTurbulenceFields
(
    const turbulenceModel& turbulence,
    scalarField& G0,
    scalarField& epsilon0
)
{
    // Accumulate all of the G and epsilon contributions
    forAll(cornerWeights_, patchi)
    {
        if (!cornerWeights_[patchi].empty())
        {
            epsilonTwoLayerWallFunctionFvPatchScalarField& epf = epsilonPatch(patchi);

            const List<scalar>& w = cornerWeights_[patchi];

            epf.calculate(turbulence, w, epf.patch(), G0, epsilon0);
        }
    }

    // Apply zero-gradient condition for epsilon
    forAll(cornerWeights_, patchi)
    {
        if (!cornerWeights_[patchi].empty())
        {
            epsilonTwoLayerWallFunctionFvPatchScalarField& epf = epsilonPatch(patchi);

            epf == scalarField(epsilon0, epf.patch().faceCells());
        }
    }
}


void Foam::epsilonTwoLayerWallFunctionFvPatchScalarField::calculate
(
    const turbulenceModel& turbModel,
    const List<scalar>& cornerWeights,
    const fvPatch& patch,
    scalarField& G0,
    scalarField& epsilon0
)
{
    const label patchi = patch.index();

    const scalar Cmu25 = pow025(wallCoeffs_.Cmu());
    const scalar kappa = wallCoeffs_.kappa();
    const scalar yPlusLam = wallCoeffs_.yPlusLam();

    const scalarField& y = turbModel.y()[patchi];

    const labelUList& faceCells = patch.faceCells();

    const tmp<scalarField> tnutw = turbModel.nut(patchi);
    const scalarField& nutw = tnutw();

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const scalarField magGradUw(mag(Uw.snGrad()));

    const volScalarField::Internal& G 
        = db().lookupObject<volScalarField::Internal>(turbModel.GName());

    forAll(nutw, facei)
    {
        const label celli = faceCells[facei];
        const scalar w = cornerWeights[facei];
        
        const scalar uStarVis = sqrt(nuw[facei]*magGradUw[facei]);
        const scalar uStarLog = Cmu25*sqrt(k[celli]);
        const scalar uStar = pow((pow(uStarVis,4)+pow(uStarLog,4)), 0.25);
        const scalar yPlus = uStar*y[facei]/nuw[facei];
        
        // With two-layer treatment, the epsilon value in the near wall region
        // is irrelevant if the first grid point is in the viscous sublayer
        //  yPlus < yPlusLam. epsilon value set here is just a boundary condition
        // for the epsilon transport equation, which gives only the High-Re 
        // epsilon value. Hence scaling the y value
        scalar yScaled = y[facei];
        
        if(yPlus < yPlusLam)
        {
            yScaled *= (yPlusLam/yPlus);
        }

        epsilon0[celli] += w * pow(uStar, 3)/(kappa*yScaled);

        // Unlike epsilon, the tke production term is important and has to be
        // blended between the viscous sublayer and log-law region.
        // k-Production term in the laminar sublayer
        // this is simply nut * (dUidxj + dUjdxi) * dUidxj and hence the same
        // as the standard production term formulation
        const scalar GVis = G[celli];

        // k-Production term in the turbulent layer:  Moulkalled Eqn. 17.104
        const scalar GLog = (nutw[facei] + nuw[facei])*magGradUw[facei] // tauw
                          * uStar / (kappa*y[facei]);

        // k-Production term with Kader blending
        const scalar gamma = -0.01*pow(yPlus,4)/(1+5*yPlus);    
        G0[celli] += w * (GVis*exp(gamma) + GLog*exp(1/gamma));
    }
}


void Foam::epsilonTwoLayerWallFunctionFvPatchScalarField::writeLocalEntries
(
    Ostream& os
) const
{
    wallFunctionBlenders::writeEntries(os);
    wallCoeffs_.writeEntries(os);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::epsilonTwoLayerWallFunctionFvPatchScalarField::
epsilonTwoLayerWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    wallFunctionBlenders(),
    initialised_(false),
    master_(-1),
    wallCoeffs_(),
    G_(),
    epsilon_(),
    cornerWeights_()
{}


Foam::epsilonTwoLayerWallFunctionFvPatchScalarField::
epsilonTwoLayerWallFunctionFvPatchScalarField
(
    const epsilonTwoLayerWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    wallFunctionBlenders(ptf),
    initialised_(false),
    master_(-1),
    wallCoeffs_(ptf.wallCoeffs_),
    G_(),
    epsilon_(),
    cornerWeights_()
{}


Foam::epsilonTwoLayerWallFunctionFvPatchScalarField::
epsilonTwoLayerWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict),
    wallFunctionBlenders(dict, blenderType::STEPWISE, scalar(2)),
    initialised_(false),
    master_(-1),
    wallCoeffs_(dict),
    G_(),
    epsilon_(),
    cornerWeights_()
{
    // Apply zero-gradient condition on start-up
    this->extrapolateInternal();
}


Foam::epsilonTwoLayerWallFunctionFvPatchScalarField::
epsilonTwoLayerWallFunctionFvPatchScalarField
(
    const epsilonTwoLayerWallFunctionFvPatchScalarField& ewfpsf
)
:
    fixedValueFvPatchField<scalar>(ewfpsf),
    wallFunctionBlenders(ewfpsf),
    initialised_(false),
    master_(-1),
    wallCoeffs_(ewfpsf.wallCoeffs_),
    G_(),
    epsilon_(),
    cornerWeights_()
{}


Foam::epsilonTwoLayerWallFunctionFvPatchScalarField::
epsilonTwoLayerWallFunctionFvPatchScalarField
(
    const epsilonTwoLayerWallFunctionFvPatchScalarField& ewfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ewfpsf, iF),
    wallFunctionBlenders(ewfpsf),
    initialised_(false),
    master_(-1),
    wallCoeffs_(ewfpsf.wallCoeffs_),
    G_(),
    epsilon_(),
    cornerWeights_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalarField& Foam::epsilonTwoLayerWallFunctionFvPatchScalarField::G
(
    bool init
)
{
    if (patch().index() == master_)
    {
        if (init)
        {
            G_ = 0.0;
        }

        return G_;
    }

    return epsilonPatch(master_).G();
}


Foam::scalarField& Foam::epsilonTwoLayerWallFunctionFvPatchScalarField::epsilon
(
    bool init
)
{
    if (patch().index() == master_)
    {
        if (init)
        {
            epsilon_ = 0.0;
        }

        return epsilon_;
    }

    return epsilonPatch(master_).epsilon(init);
}


void Foam::epsilonTwoLayerWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const auto& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    setMaster();

    if (patch().index() == master_)
    {
        createAveragingWeights();
        calculateTurbulenceFields(turbModel, G(true), epsilon(true));
    }

    const scalarField& G0 = this->G();
    const scalarField& epsilon0 = this->epsilon();

    typedef DimensionedField<scalar, volMesh> FieldType;

    FieldType& G = db().lookupObjectRef<FieldType>(turbModel.GName());

    FieldType& epsilon = const_cast<FieldType&>(internalField());

    forAll(*this, facei)
    {
        const label celli = patch().faceCells()[facei];

        G[celli] = G0[celli];
        epsilon[celli] = epsilon0[celli];
    }

    fvPatchField<scalar>::updateCoeffs();
}


void Foam::epsilonTwoLayerWallFunctionFvPatchScalarField::updateWeightedCoeffs
(
    const scalarField& weights
)
{
    if (updated())
    {
        return;
    }

    const auto& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    setMaster();

    if (patch().index() == master_)
    {
        createAveragingWeights();
        calculateTurbulenceFields(turbModel, G(true), epsilon(true));
    }

    const scalarField& G0 = this->G();
    const scalarField& epsilon0 = this->epsilon();

    typedef DimensionedField<scalar, volMesh> FieldType;

    FieldType& G = db().lookupObjectRef<FieldType>(turbModel.GName());

    FieldType& epsilon = const_cast<FieldType&>(internalField());

    scalarField& epsilonf = *this;

    // Only set the values if the weights are > tolerance
    forAll(weights, facei)
    {
        const scalar w = weights[facei];

        if (w > tolerance_)
        {
            const label celli = patch().faceCells()[facei];

            G[celli] = (1.0 - w)*G[celli] + w*G0[celli];
            epsilon[celli] = (1.0 - w)*epsilon[celli] + w*epsilon0[celli];
            epsilonf[facei] = epsilon[celli];
        }
    }

    fvPatchField<scalar>::updateCoeffs();
}


void Foam::epsilonTwoLayerWallFunctionFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& matrix
)
{
    if (manipulatedMatrix())
    {
        return;
    }

    matrix.setValues(patch().faceCells(), patchInternalField());

    fvPatchField<scalar>::manipulateMatrix(matrix);
}


void Foam::epsilonTwoLayerWallFunctionFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& matrix,
    const Field<scalar>& weights
)
{
    if (manipulatedMatrix())
    {
        return;
    }

    DynamicList<label> constraintCells(weights.size());
    DynamicList<scalar> constraintValues(weights.size());
    const labelUList& faceCells = patch().faceCells();

    const DimensionedField<scalar, volMesh>& fld = internalField();

    forAll(weights, facei)
    {
        // Only set the values if the weights are > tolerance
        if (weights[facei] > tolerance_)
        {
            const label celli = faceCells[facei];

            constraintCells.append(celli);
            constraintValues.append(fld[celli]);
        }
    }

    if (debug)
    {
        Pout<< "Patch: " << patch().name()
            << ": number of constrained cells = " << constraintCells.size()
            << " out of " << patch().size()
            << endl;
    }

    matrix.setValues(constraintCells, constraintValues);

    fvPatchField<scalar>::manipulateMatrix(matrix);
}


void Foam::epsilonTwoLayerWallFunctionFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);
    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        epsilonTwoLayerWallFunctionFvPatchScalarField
    );
}


// ************************************************************************* //

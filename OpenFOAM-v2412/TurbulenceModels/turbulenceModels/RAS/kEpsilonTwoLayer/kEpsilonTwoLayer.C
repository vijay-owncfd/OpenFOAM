/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

#include "kEpsilonTwoLayer.H"
#include "fvOptions.H"
#include "bound.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void kEpsilonTwoLayer<BasicTurbulenceModel>::correctNut()
{
    this->nut_ = Cmu_*sqr(k_)/epsilon_;

    if(twoLayerTreatment_)
    {
        Rey_ = mag(y_*sqrt(k_.oldTime())/this->nu());
        volScalarField l_mu(y_*ClStar_*(1.0 - exp(-1.0*Rey_/Amu_)));

        lEps_ = y_*ClStar_*(1.0 - exp(-1.0*Rey_/Aeps_));

        dimensionedScalar A = ReyFactor_*ReyStar_/atanh(0.98);
        lambdaEps_ = 0.5*(1.0 + tanh((Rey_ - ReyStar_)/A)); 

        this->nut_ *= lambdaEps_;   // high Re part
        this->nut_ += (1.0-lambdaEps_)*Cmu_*l_mu*sqrt(k_);  // low Re part
    }

    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> kEpsilonTwoLayer<BasicTurbulenceModel>::kSource() const
{
    return tmp<fvScalarMatrix>::New
    (
        k_,
        dimVolume*this->rho_.dimensions()*k_.dimensions()/dimTime
    );
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> kEpsilonTwoLayer<BasicTurbulenceModel>::epsilonSource() const
{
    return tmp<fvScalarMatrix>::New
    (
        epsilon_,
        dimVolume*this->rho_.dimensions()*epsilon_.dimensions()/dimTime
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
kEpsilonTwoLayer<BasicTurbulenceModel>::kEpsilonTwoLayer
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity<RASModel<BasicTurbulenceModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    Cmu_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cmu",
            this->coeffDict_,
            0.09
        )
    ),
    C1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C1",
            this->coeffDict_,
            1.44
        )
    ),
    C2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C2",
            this->coeffDict_,
            1.92
        )
    ),
    C3_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C3",
            this->coeffDict_,
            0
        )
    ),
    sigmak_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "sigmak",
            this->coeffDict_,
            1.0
        )
    ),
    sigmaEps_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "sigmaEps",
            this->coeffDict_,
            1.3
        )
    ),
    twoLayerTreatment_
    (
        Switch::getOrAddToDict
        (
            "twoLayerTreatment",
            this->coeffDict_,
            false
        )
    ),        
    ClStar_(0.41*pow(0.09, -0.75)),
    Amu_(70.0),
    Aeps_(2*ClStar_),
    ReyStar_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "ReyStar",
            this->coeffDict_,
            200.0
        )
    ),
    ReyFactor_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "ReyFactor",
            this->coeffDict_,
            0.05
        )
    ),
    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    epsilon_
    (
        IOobject
        (
            IOobject::groupName("epsilon", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    Rey_
    (
        IOobject
        (
            "Rey",
            this->runTime_.name(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar
        (
            "Rey",
            dimless,
            200.0
        )
    ),
    lEps_
    (
        IOobject
        (
            "lEps",
            this->runTime_.name(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar
        (
            "lEps",
            dimLength,
            0.0
        )
    ),
    lambdaEps_
    (
        IOobject
        (
            "lambdaEps",
            this->runTime_.name(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar
        (
            "lambdaEps",
            dimless,
            1.0
        )
    ),
    y_(wallDist::New(this->mesh_).y())

{
    bound(k_, this->kMin_);
    bound(epsilon_, this->epsilonMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool kEpsilonTwoLayer<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        Cmu_.readIfPresent(this->coeffDict());
        C1_.readIfPresent(this->coeffDict());
        C2_.readIfPresent(this->coeffDict());
        C3_.readIfPresent(this->coeffDict());
        sigmak_.readIfPresent(this->coeffDict());
        sigmaEps_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
void kEpsilonTwoLayer<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    const volScalarField& nut = this->nut_;

    fv::options& fvOptions(fv::options::New(this->mesh_));

    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    const volScalarField::Internal divU
    (
        fvc::div(fvc::absolute(this->phi(), U))().v()
    );

    tmp<volTensorField> tgradU = fvc::grad(U);
    const volScalarField::Internal GbyNu
    (
        IOobject::scopedName(this->type(), "GbyNu"),
        tgradU().v() && devTwoSymm(tgradU().v())
    );
    const volScalarField::Internal G(this->GName(), nut()*GbyNu);
    tgradU.clear();

    // Update epsilon and G at the wall
    epsilon_.boundaryFieldRef().updateCoeffs();
    // Push any changed cell values to coupled neighbours
    epsilon_.boundaryFieldRef().template evaluateCoupled<coupledFvPatch>();

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(alpha, rho, epsilon_)
      + fvm::div(alphaRhoPhi, epsilon_)
      - fvm::laplacian(alpha*rho*DepsilonEff(), epsilon_)
     ==
        C1_*alpha()*rho()*GbyNu*Cmu_*k_()
      - fvm::SuSp(((2.0/3.0)*C1_ - C3_)*alpha()*rho()*divU, epsilon_)
      - fvm::Sp(C2_*alpha()*rho()*epsilon_()/k_(), epsilon_)
      + epsilonSource()
      + fvOptions(alpha, rho, epsilon_)
    );

    epsEqn.ref().relax();
    fvOptions.constrain(epsEqn.ref());
    epsEqn.ref().boundaryManipulate(epsilon_.boundaryFieldRef());
    solve(epsEqn);
    fvOptions.correct(epsilon_);
    bound(epsilon_, this->epsilonMin_);

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
     ==
        alpha()*rho()*G
      - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, k_)
      - (
            twoLayerTreatment_ ?
                (
                    fvm::Sp(alpha()*rho()*lambdaEps_()*epsilon_()/k_(), k_)
                    + fvm::Sp(alpha()*rho()*(1.0-lambdaEps_())*sqrt(k_())/lEps_(), k_)
                )
              : fvm::Sp(alpha()*rho()*epsilon_()/k_(), k_)
        )
      + kSource()
      + fvOptions(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);

    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //

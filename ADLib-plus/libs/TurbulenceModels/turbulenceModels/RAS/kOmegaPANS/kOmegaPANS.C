/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "kOmegaPANS.H"
#include "fvOptions.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void kOmegaPANS<BasicTurbulenceModel>::correctNut()
{
	const dimensionedScalar small("small", this->nut_.dimensions(), VSMALL);

    this->nut_ = max(min(k_/omega_, 100000*this->nu()), small);

    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}

template<class BasicTurbulenceModel>
void kOmegaPANS<BasicTurbulenceModel>::updateFk()
{
	// Evaluate length scales
	const dimensionedScalar smalld("smalld", dimensionSet(0,0,-1,0,0), SMALL);
	const volScalarField lturb = sqrt(k_)/max(betaStar_*omega_*pow(fk_,1.5),smalld);
	Info<< "min/max lturb: " << min(lturb).value() << "/" << max(lturb).value() << endl;

	const volScalarField eta = pow(pow(this->nu(),3)/(betaStar_*omega_*k_), 0.25);
	Info<< "min/max eta  : " << min(eta).value() << "/" << max(eta).value() << endl;

	const scalar twoThird = 2.0/3.0;
	const dimensionedScalar small("small", fk_.dimensions(), VSMALL);

	fk_ = min
	(
	 	max
		(
			pow(delta()/eta, twoThird)
		  * (1.0-pow(eta/delta(),twoThird))
		  * pow(eta/lturb,twoThird)
		  , VSMALL
		),
		1.0
	);

	fw_ = min(1.0/fk_, 1e10);

	Info<< "min/max fk   : " << min(fk_).value() << "/" << max(fk_).value() << endl;
	Info<< "min/max fw   : " << min(fw_).value() << "/" << max(fw_).value() << endl;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
kOmegaPANS<BasicTurbulenceModel>::kOmegaPANS
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

    betaStar_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "betaStar",
            this->coeffDict_,
            0.09
        )
    ),
    beta_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "beta",
            this->coeffDict_,
            0.075
        )
    ),
    gamma_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "gamma",
            this->coeffDict_,
            5.0/9.0
        )
    ),
    alphaK_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaK",
            this->coeffDict_,
            0.5
        )
    ),
    alphaOmega_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaOmega",
            this->coeffDict_,
            0.5
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
    omega_
    (
        IOobject
        (
            IOobject::groupName("omega", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
	delta_
	(
	 	LESdelta::New
		(
		 	IOobject::groupName("delta", alphaRhoPhi.group()),
			*this,
			this->coeffDict_
		)
	),
	fk_
	(
        IOobject
        (
            IOobject::groupName("fk", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
	),
	fw_
	(
        IOobject
        (
            IOobject::groupName("fw", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		1.0/fk_
	)
{
    bound(k_, this->kMin_);
    bound(omega_, this->omegaMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool kOmegaPANS<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        betaStar_.readIfPresent(this->coeffDict());
        beta_.readIfPresent(this->coeffDict());
        gamma_.readIfPresent(this->coeffDict());
        alphaK_.readIfPresent(this->coeffDict());
        alphaOmega_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
void kOmegaPANS<BasicTurbulenceModel>::correct()
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
        tgradU().v() && dev(twoSymm(tgradU().v()))
    );
    const volScalarField::Internal G(this->GName(), nut()*GbyNu);
    tgradU.clear();

    // Update omega and G at the wall
    omega_.boundaryFieldRef().updateCoeffs();

	// Update control parameter
	updateFk();

	// Update model coefficient betaPns
	const volScalarField::Internal betaPns
	(
	 	max
		(
		 	gamma_*betaStar_ + fk_*(beta_-gamma_*betaStar_),
			VSMALL
		)
	);

	/*
	volScalarField DfwDt
	(
	 	fvc::ddt(alpha*rho, fw_) + fvc::div(alphaRhoPhi, fw_)
	);

	Info<< "min/max(DfwDt): " << min(DfwDt).value() << "/" << max(DfwDt).value() << endl;
	*/

    // Turbulence specific dissipation rate equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(alpha, rho, omega_)
      + fvm::div(alphaRhoPhi, omega_)
      - fvm::laplacian(alpha*rho*DomegaEff(), omega_)
     ==
        gamma_*alpha()*rho()*GbyNu
      - fvm::SuSp(((2.0/3.0)*gamma_)*alpha()*rho()*divU, omega_)
      - fvm::Sp(beta_*alpha()*rho()*omega_(), omega_)
      + fvOptions(alpha, rho, omega_)
//	  + fvm::SuSp(DfwDt, omega_)
    );

    omegaEqn.ref().relax();
    fvOptions.constrain(omegaEqn.ref());
    omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());
    solve(omegaEqn);
    fvOptions.correct(omega_);

    bound(omega_, this->omegaMin_);

	volScalarField DfkDt
	(
	 	fvc::ddt(alpha*rho, fk_) + fvc::div(alphaRhoPhi, fk_)
	);

	Info<< "min/max(DfkDt): " << min(DfkDt).value() << "/" << max(DfkDt).value() << endl;

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
     ==
        alpha()*rho()*G
      - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, k_)
      - fvm::Sp(betaStar_*alpha()*rho()*omega_(), k_)
      + fvOptions(alpha, rho, k_)
	  + fvm::SuSp(DfkDt, k_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);

//	k_ = min(k_, 1.5*sqr(mag(U)));
    bound(k_, this->kMin_);

	Info<< "min(omega): " << min(omega_).value() << endl;
	Info<< "min(k)    : " << min(k_).value() << endl;

    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //

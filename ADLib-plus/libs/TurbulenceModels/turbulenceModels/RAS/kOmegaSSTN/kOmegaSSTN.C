/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "kOmegaSSTN.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaSSTN<BasicTurbulenceModel>::F1
(
    const volScalarField& CDkOmega
) const
{
    tmp<volScalarField> CDkOmegaPlus = max
    (
        CDkOmega,
        dimensionedScalar(dimless/sqr(dimTime), 1.0e-21)
    );

    tmp<volScalarField> arg1 = min
    (
        min
        (
            max
            (
                (scalar(1)/this->betaStar_)*sqrt(this->k_)/(this->omega_*this->y_),
                scalar(500)*(this->mu()/this->rho_)/(sqr(this->y_)*this->omega_)
            ),
            (4*this->alphaOmega2_)*this->k_/(CDkOmegaPlus*sqr(this->y_))
        ),
        scalar(10)
    );

    return tanh(pow4(arg1));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaSSTN<BasicTurbulenceModel>::F4
(
    const volTensorField& gradU
) const
{
    tmp<volScalarField> f4 =
		tmp<volScalarField>::New
		(
		 	IOobject
			(
			 	"f4",
				this->runTime_.timeName(),
				this->mesh_,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			this->mesh_,
			dimensionedScalar(dimless, scalar(1))
		);

    if (F4_)
    {
		const scalar CRC(1.4);

    	const volScalarField S2
		(
		 	max
			(
			 	this->S2(gradU), 
				dimensionedScalar("small", dimensionSet(0,0,-2,0,0), 1e-10)
			)
		);
		// Limiting W2 not to exceed S2
    	const volScalarField W2
		(
		 	min(this->W2(gradU), S2)
		);

		const volScalarField WbyS(sqrt(W2)/sqrt(S2));
		const volScalarField Ri(WbyS*(WbyS-scalar(1)));

		/*
		Pout<< "max(WbyS): " << max(WbyS) 
			<< "/min(WbyS): " << min(WbyS)
			<< endl;
		Pout<< "max(Ri): " << max(Ri) 
			<< "/min(Ri): " << min(Ri)
			<< endl;
		*/
		
		f4.ref() /= (scalar(1)+CRC*Ri);
    }

    return f4;
}

template<class BasicTurbulenceModel>
Foam::tmp<Foam::volScalarField> kOmegaSSTN<BasicTurbulenceModel>::W2
(
    const volTensorField& gradU
) const
{
    return 2*magSqr(skew(gradU));
}

template<class BasicTurbulenceModel>
void kOmegaSSTN<BasicTurbulenceModel>::correctNut(const volScalarField& S2)
{
    // Correct the turbulence viscosity
    this->nut_ = this->a1_*this->k_/max(this->a1_*this->omega_, this->b1_*this->F23()*sqrt(S2));
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    // Correct the turbulence thermal diffusivity
    BasicTurbulenceModel::correctNut();
}


template<class BasicTurbulenceModel>
void kOmegaSSTN<BasicTurbulenceModel>::correctNut()
{
    // correctNut(2*magSqr(symm(fvc::grad(this->U_))));
    correctNut(2*magSqr(skew(fvc::grad(this->U_))));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
kOmegaSSTN<BasicTurbulenceModel>::kOmegaSSTN
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
    kOmegaSSTBase<eddyViscosity<RASModel<BasicTurbulenceModel>>>
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
    F4_
    (
        Switch::getOrAddToDict
        (
            "F4",
            this->coeffDict_,
            false
        )
    )
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool kOmegaSSTN<BasicTurbulenceModel>::read()
{
    if (BasicTurbulenceModel::read())
    {
        this->alphaK1_.readIfPresent(this->coeffDict());
        this->alphaK2_.readIfPresent(this->coeffDict());
        this->alphaOmega1_.readIfPresent(this->coeffDict());
        this->alphaOmega2_.readIfPresent(this->coeffDict());
        this->gamma1_.readIfPresent(this->coeffDict());
        this->gamma2_.readIfPresent(this->coeffDict());
        this->beta1_.readIfPresent(this->coeffDict());
        this->beta2_.readIfPresent(this->coeffDict());
        this->betaStar_.readIfPresent(this->coeffDict());
        this->a1_.readIfPresent(this->coeffDict());
        this->b1_.readIfPresent(this->coeffDict());
        this->c1_.readIfPresent(this->coeffDict());
        this->F3_.readIfPresent("F3", this->coeffDict());
        F4_.readIfPresent("F4", this->coeffDict());

        this->setDecayControl(this->coeffDict());

        return true;
    }

    return false;
}

template<class BasicTurbulenceModel>
void kOmegaSSTN<BasicTurbulenceModel>::correct()
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
    volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    BasicTurbulenceModel::correct();

    const volScalarField::Internal divU
    (
        fvc::div(fvc::absolute(this->phi(), U))
    );

    tmp<volTensorField> tgradU = fvc::grad(U);
    const volScalarField S2(this->S2(tgradU()));
    volScalarField::Internal GbyNu0(this->GbyNu0(tgradU(), S2));
    volScalarField::Internal G(this->GName(), nut*GbyNu0);


    // - boundary condition changes a cell value
    // - normally this would be triggered through correctBoundaryConditions
    // - which would do
    //      - fvPatchField::evaluate() which calls
    //      - fvPatchField::updateCoeffs()
    // - however any processor boundary conditions already start sending
    //   at initEvaluate so would send over the old value.
    // - avoid this by explicitly calling updateCoeffs early and then
    //   only doing the boundary conditions that rely on initEvaluate
    //   (currently only coupled ones)

    //- 1. Explicitly swap values on coupled boundary conditions
    // Update omega and G at the wall
    this->omega_.boundaryFieldRef().updateCoeffs();
    // omegaWallFunctions change the cell value! Make sure to push these to
    // coupled neighbours. Note that we want to avoid the re-updateCoeffs
    // of the wallFunctions so make sure to bypass the evaluate on
    // those patches and only do the coupled ones.
    this->omega_.boundaryFieldRef().template evaluateCoupled<coupledFvPatch>();

    ////- 2. Make sure the boundary condition calls updateCoeffs from
    ////     initEvaluate
    ////     (so before any swap is done - requires all coupled bcs to be
    ////      after wall bcs. Unfortunately this conflicts with cyclicACMI)
    //omega_.correctBoundaryConditions();


    const volScalarField CDkOmega
    (
        (2*this->alphaOmega2_)*(fvc::grad(this->k_) & fvc::grad(this->omega_))/this->omega_
    );

    const volScalarField F1(this->F1(CDkOmega));
    const volScalarField F23(this->F23());

    {
        const volScalarField::Internal gamma(this->gamma(F1));
        const volScalarField::Internal beta(this->beta(F1));

        GbyNu0 = this->GbyNu(GbyNu0, this->F23(), S2());

        // Turbulent frequency equation
        tmp<fvScalarMatrix> omegaEqn
        (
            fvm::ddt(alpha, rho, this->omega_)
          + fvm::div(alphaRhoPhi, this->omega_)
          - fvm::laplacian(alpha*rho*this->DomegaEff(F1), this->omega_)
         ==
            alpha()*rho()*gamma*GbyNu0
          - fvm::SuSp((2.0/3.0)*alpha()*rho()*gamma*divU, this->omega_)
          - fvm::Sp(F4(tgradU())()*alpha()*rho()*beta*this->omega_(), this->omega_)
          + fvm::Sp
            (
                alpha()*rho()*(scalar(1) - F1())*CDkOmega()/this->omega_(),
                this->omega_
            )
          + alpha()*rho()*beta*sqr(this->omegaInf_)
          + this->Qsas(S2(), gamma, beta)
          + this->omegaSource()
          + fvOptions(alpha, rho, this->omega_)
        );

        omegaEqn.ref().relax();
        fvOptions.constrain(omegaEqn.ref());
        omegaEqn.ref().boundaryManipulate(this->omega_.boundaryFieldRef());
        solve(omegaEqn);
        fvOptions.correct(this->omega_);
        bound(this->omega_, this->omegaMin_);
    }

    {
        // Turbulent kinetic energy equation
        tmp<fvScalarMatrix> kEqn
        (
            fvm::ddt(alpha, rho, this->k_)
          + fvm::div(alphaRhoPhi, this->k_)
          - fvm::laplacian(alpha*rho*this->DkEff(F1), this->k_)
         ==
            alpha()*rho()*this->Pk(G)
          - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, this->k_)
          - fvm::Sp(alpha()*rho()*this->epsilonByk(F1, tgradU()), this->k_)
          + alpha()*rho()*this->betaStar_*this->omegaInf_*this->kInf_
          + this->kSource()
          + fvOptions(alpha, rho, this->k_)
        );

        tgradU.clear();

        kEqn.ref().relax();
        fvOptions.constrain(kEqn.ref());
        solve(kEqn);
        fvOptions.correct(this->k_);
        bound(this->k_, this->kMin_);
    }

    correctNut(S2);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "LLEASM.H"
#include "fvOptions.H"
#include "bound.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> LLEASM<BasicTurbulenceModel>::F1
(
    const volScalarField& CDkOmega
) const
{
    tmp<volScalarField> CDkOmegaPlus = max
    (
        CDkOmega,
        dimensionedScalar("1.0e-10", dimless/sqr(dimTime), 1.0e-10)
    );

    tmp<volScalarField> arg1 = min
    (
        min
        (
            max
            (
                (scalar(1)/betaStar_)*sqrt(k_)/(omega_*y_),
                scalar(500)*(this->mu()/this->rho_)/(sqr(y_)*omega_)
            ),
            (4*alphaOmega2_)*k_/(CDkOmegaPlus*sqr(y_))
        ),
        scalar(10)
    );

    return tanh(pow4(arg1));
}

template<class BasicTurbulenceModel>
tmp<volScalarField>
LLEASM<BasicTurbulenceModel>::F2() const
{
    tmp<volScalarField> arg2 = min
    (
        max
        (
            (scalar(2)/betaStar_)*sqrt(k_)/(omega_*y_),
            scalar(500)*(this->mu()/this->rho_)/(sqr(y_)*omega_)
        ),
        scalar(100)
    );

    return tanh(sqr(arg2));
}

template<class BasicTurbulenceModel>
tmp<volScalarField>
LLEASM<BasicTurbulenceModel>::F3() const
{
    tmp<volScalarField> arg3 = min
    (
        150*(this->mu()/this->rho_)/(omega_*sqr(y_)),
        scalar(10)
    );

    return 1 - tanh(pow4(arg3));
}

template<class BasicTurbulenceModel>
tmp<volScalarField>
LLEASM<BasicTurbulenceModel>::F23() const
{
    tmp<volScalarField> f23(F2());

    if (F3_)
    {
        f23.ref() *= F3();
    }

    return f23;
}


template<class BasicTurbulenceModel>
void LLEASM<BasicTurbulenceModel>::correctNut
(
    const volScalarField& S2,
    const volScalarField& F2
)
{
    this->nut_ = a1_*k_/max(a1_*omega_, b1_*F2*sqrt(S2));
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}



// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void LLEASM<BasicTurbulenceModel>::correctNut()
{
	correctNut(2*magSqr(symm(fvc::grad(this->U_))), F23());
}

template<class BasicTurbulenceModel>
tmp<volScalarField> LLEASM<BasicTurbulenceModel>::epsilonByk
(
    const volScalarField& F1,
    const volScalarField& F2
) const
{
    return betaStar_*omega_;
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix>
LLEASM<BasicTurbulenceModel>::kSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            k_,
            dimVolume*this->rho_.dimensions()*k_.dimensions()/dimTime
        )
    );
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix>
LLEASM<BasicTurbulenceModel>::omegaSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            omega_,
            dimVolume*this->rho_.dimensions()*omega_.dimensions()/dimTime
        )
    );
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> LLEASM<BasicTurbulenceModel>::Qsas
(
    const volScalarField& S2,
    const volScalarField& gamma,
    const volScalarField& beta
) const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            omega_,
            dimVolume*this->rho_.dimensions()*omega_.dimensions()/dimTime
        )
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
LLEASM<BasicTurbulenceModel>::LLEASM
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
    AnisotropyReynoldsStress<RASModel<BasicTurbulenceModel>>
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

    alphaK1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK1",
            this->coeffDict_,
            0.85
        )
    ),
    alphaK2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK2",
            this->coeffDict_,
            1.0
        )
    ),
    alphaOmega1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega1",
            this->coeffDict_,
            0.5
        )
    ),
    alphaOmega2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega2",
            this->coeffDict_,
            0.856
        )
    ),
    gamma1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma1",
            this->coeffDict_,
            5.0/9.0
        )
    ),
    gamma2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma2",
            this->coeffDict_,
            0.44
        )
    ),
    beta1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta1",
            this->coeffDict_,
            0.075
        )
    ),
    beta2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta2",
            this->coeffDict_,
            0.0828
        )
    ),
    betaStar_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaStar",
            this->coeffDict_,
            0.09
        )
    ),
    a1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "a1",
            this->coeffDict_,
            0.31
        )
    ),
    b1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "b1",
            this->coeffDict_,
            1.0
        )
    ),
    c1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "c1",
            this->coeffDict_,
            10.0
        )
    ),
    F3_
    (
        Switch::lookupOrAddToDict
        (
            "F3",
            this->coeffDict_,
            false
        )
    ),
    c_
    (
        Switch::lookupOrAddToDict
        (
            "c",
            this->coeffDict_,
            false
        )
    ),
	B2_("B2", dimless, (8.0+c_.value())/11.0),
	B3_("B3", dimless, (8.0*c_.value()-2.0)/11.0),
	B4_("B4", dimless, (60.0*c_.value()-4.0)/55.0),
	B5_("B5", dimless, (6.0*c_.value()+4.0)/11.0),

    y_(wallDist::New(this->mesh_).y()),

    k_
    (
        IOobject
        (
            IOobject::groupName("k", U.group()),
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
            IOobject::groupName("omega", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    )
{
    bound(k_, this->kMin_);
    bound(omega_, this->omegaMin_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool LLEASM<BasicTurbulenceModel>::read()
{
    if (AnisotropyReynoldsStress<RASModel<BasicTurbulenceModel>>::read())
    {
        alphaK1_.readIfPresent(this->coeffDict());
        alphaK2_.readIfPresent(this->coeffDict());
        alphaOmega1_.readIfPresent(this->coeffDict());
        alphaOmega2_.readIfPresent(this->coeffDict());
        gamma1_.readIfPresent(this->coeffDict());
        gamma2_.readIfPresent(this->coeffDict());
        beta1_.readIfPresent(this->coeffDict());
        beta2_.readIfPresent(this->coeffDict());
        betaStar_.readIfPresent(this->coeffDict());
        a1_.readIfPresent(this->coeffDict());
        b1_.readIfPresent(this->coeffDict());
        c1_.readIfPresent(this->coeffDict());
        F3_.readIfPresent("F3", this->coeffDict());
        c_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
void LLEASM<BasicTurbulenceModel>::correct()
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

    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));

    tmp<volTensorField> tgradU = fvc::grad(U);

	// Tensorial fields
	volSymmTensorField S(symm(tgradU()));
	volTensorField W(skew(tgradU()));

    volScalarField S2(2*magSqr(S));
    volScalarField GbyNu((tgradU() && dev(twoSymm(tgradU()))));
    volScalarField G(this->GName(), nut*GbyNu);
    tgradU.clear();

    // Update omega and G at the wall
    omega_.boundaryFieldRef().updateCoeffs();

    volScalarField CDkOmega
    (
        (2*alphaOmega2_)*(fvc::grad(k_) & fvc::grad(omega_))/omega_
    );

    volScalarField F1(this->F1(CDkOmega));
    volScalarField F23(this->F23());

    {
        volScalarField gamma(this->gamma(F1));
        volScalarField beta(this->beta(F1));

        // Turbulent frequency equation
        tmp<fvScalarMatrix> omegaEqn
        (
            fvm::ddt(alpha, rho, omega_)
          + fvm::div(alphaRhoPhi, omega_)
          - fvm::laplacian(alpha*rho*DomegaEff(F1), omega_)
         ==
            alpha*rho*gamma
           *min
            (
                GbyNu,
                (c1_/a1_)*betaStar_*omega_*max(a1_*omega_, b1_*F23*sqrt(S2))
            )
          - fvm::SuSp((2.0/3.0)*alpha*rho*gamma*divU, omega_)
          - fvm::Sp(alpha*rho*beta*omega_, omega_)
          - fvm::SuSp
            (
                alpha*rho*(F1 - scalar(1))*CDkOmega/omega_,
                omega_
            )
          + Qsas(S2, gamma, beta)
          + omegaSource()
          + fvOptions(alpha, rho, omega_)
        );

        omegaEqn.ref().relax();
        fvOptions.constrain(omegaEqn.ref());
        omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());
        solve(omegaEqn);
        fvOptions.correct(omega_);
        bound(omega_, this->omegaMin_);
    }

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(F1), k_)
     ==
        min(alpha*rho*G, (c1_*betaStar_)*alpha*rho*k_*omega_)
      - fvm::SuSp((2.0/3.0)*alpha*rho*divU, k_)
      - fvm::Sp(alpha*rho*epsilonByk(F1, F23), k_)
      + kSource()
      + fvOptions(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);

    correctNut(S2, F23);

	// Evaluate the Reynolds stress anisotropy tensor
	tmp<volScalarField> Skk = tr(S)/this->epsilonByk(F1, this->F2());
	S = dev(S)/this->epsilonByk(F1, this->F2());
	W = W/this->epsilonByk(F1, this->F2());

	const scalar L3(B2_.value() + B3_.value() - 1.0);
	const scalar L4(B2_.value() - B3_.value() - 1.0);
	const scalar L2(2./3.*L3 - 0.5*B4_.value());

	// Model coefficient - need to be modified later
	const scalar C10 = 3.0;

	forAll(k_, i)
	{
		const scalar ktol = 1.5*magSqr(0.0001*U[i]);

		if (k_[i] > ktol)
		{
			const scalar L1
			(
				0.5*C10 - 1.0 
			  + (B5_.value() - 2./3.*L3)*Skk()[i] 
			  + G[i]/epsilon().ref()[i]
			);
			scalarSquareMatrix A(6, 0.0);

			symmTensor& Sij = S[i];
			tensor& Wij = W[i];
			const scalar L323(2./3.*L3);
			const scalar L343(2.0*L323);

			scalarList bij(6);
			scalarList RHS(6);

			// Matrix coefficient
			A[0][0] = L1 - L343*Sij.xx();
			A[0][1] = -L323*Sij.xy() - 2.0*L4*Wij.xy();
			A[0][2] = -L343*Sij.xz() - 2.0*L4*Wij.xz();
			A[0][3] = L323*Sij.yy();
			A[0][4] = L343*Sij.yz();
			A[0][5] = L323*Sij.zz();

			A[1][0] = -L3*Sij.xy() - L4*Wij.yx();
			A[1][1] = L1 - L3*(Sij.xx() + Sij.yy());
			A[1][2] = -L3*Sij.yz() - L4*Wij.yz();
			A[1][3] = -L3*Sij.xy() - L4*Wij.xy();
			A[1][4] = -L3*Sij.xz() - L4*Wij.xz();
			A[1][5] = 0.0;

			A[2][0] = -L3*Sij.xz() - L4*Wij.zx();
			A[2][1] = -L3*Sij.yz() - L4*Wij.zy();
			A[2][2] = L1 - L3*(Sij.xx() + Sij.zz());
			A[2][3] = 0.0;
			A[2][4] = -L3*Sij.xy() - L4*Wij.xy();
			A[2][5] = -L3*Sij.xz() - L4*Wij.xz();

			A[3][0] = L323*Sij.xx();
			A[3][1] = -L323*Sij.xy() - 2.0*L4*Wij.yx();
			A[3][2] = L343*Sij.xz();
			A[3][3] = L1 - L343*Sij.yy();
			A[3][4] = -L323*Sij.yz() - 2.0*L4*Wij.yz();
			A[3][5] = L323*Sij.zz();

			A[4][0] = 0.0;
			A[4][1] = -L3*Sij.xz() - L4*Wij.zx();
			A[4][2] = -L3*Sij.xy() - L4*Wij.yx();
			A[4][3] = -L3*Sij.yz() - L4*Wij.zy();
			A[4][4] = L1 - L3*(Sij.yy() + Sij.zz());
			A[4][5] = -L3*Sij.yz() - L4*Wij.yz();

			A[5][0] = -L323*Sij.xx();
			A[5][1] = L343*Sij.xy();
			A[5][2] = -L323*Sij.xz() - 2.0*L4*Wij.zx();
			A[5][3] = L323*Sij.yy();
			A[5][4] = -L323*Sij.yz() - 2.0*L4*Wij.zy();
			A[5][5] = L1 - L343*Sij.zz();

			// RHS
			RHS[0] = L2*Sij.xx();
			RHS[1] = L2*Sij.xy();
			RHS[2] = L2*Sij.xz();
			RHS[3] = L2*Sij.yy();
			RHS[4] = L2*Sij.yz();
			RHS[5] = L2*Sij.zz();

			solve(bij, A, RHS);

			this->B_[i] = 
			symmTensor
			(
				bij[0], bij[1], bij[2], 
						bij[3], bij[4], 
								bij[5]
			);
		}
		else
		{
			this->B_[i] = sphericalTensor(1.0/3.0);
		}
	}

    // Correct wall shear-stresses when applying wall-functions
	this->update(k_);
    this->correctWallShearStress(this->R_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //

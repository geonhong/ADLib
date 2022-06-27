/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
	Released 2004-2011 OpenCFD Ltd. 
    Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "RelativeMotion.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvMatrices.H"
#include "fvcMeshPhi.H"

// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * * //

template<class RhoFieldType>
void Foam::RelativeMotion::makeRelativeRhoFlux
(
 	const RhoFieldType& rho,
	surfaceScalarField& phi
) const
{
	if (!active_)
	{
		return;
	}

	const surfaceVectorField& Cf = mesh_.Cf();
	const surfaceVectorField& Sf = mesh_.Sf();

	// const vector Omega_->value(mesh_.time().timeOutputValue())*axis_;

	const vectorField& Cfi = Cf;
	const vectorField& Sfi = Sf;
	scalarField& phii = phi.primitiveFieldRef();

	// Internal faces
	forAll(phi, i)
	{
		label facei = i;

		vector rvel = Zero;
		vector lvel = Zero; 

		if (rm_.faceType(facei) == 1)
		{
			rvel = rm_.velocity(Cfi[facei]);
		}

		if (lm_.faceType(facei) == 1)
		{
			lvel = lm_.velocity();
		}

		phii[facei] -= rho[facei]*(rvel+lvel) & Sfi[facei];
		
	}

	//forAll(internalFaces_, i)
	//{
	//		label facei = internalFaces_[i];
	//		phii[facei] -= rho[facei]*(Omega ^ (Cfi[facei] - origin_)) & Sfi[facei];
	//}
	
	makeRelativeRhoFlux(rho.boundaryField(), phi.boundaryFieldRef());
}

template<class RhoFieldType>
void Foam::RelativeMotion::makeRelativeRhoFlux
(
 	const RhoFieldType& rho,
	FieldField<fvsPatchField, scalar>& phi
) const
{
	if (!active_)
	{
		return;
	}

    // const surfaceVectorField& Cf = mesh_.Cf();
	// const surfaceVectorField& Sf = mesh_.Sf();
	
	// const vector Omega =	omega_->value(mesh_.time().timeOutputValue())*axis_;
	
	// Included patches
	const polyBoundaryMesh& patches = mesh_.boundaryMesh();

	forAll(patches, patchi)
	{
		makeRelativeRhoFlux(rho[patchi], phi[patchi], patchi);
	}

	/*
	forAll(includedFaces, patchi)
	{
		forAll(includedFaces_[patchi], i)
		{
			label patchFacei = includedFaces_[patchi][i];

			phi[patchi][patchFacei] = 0.0;
		}
	}

	// Excluded patches
	forAll(excludedFaces_, patchi)
	{
		forAll(excludedFaces_[patchi], i)
		{
			label patchFacei = excludedFaces_[patchi][i];	

			phi[patchi][patchFacei] -=
				rho[patchi][patchFacei]
				* (Omega ^ (Cf.boundaryField()[patchi][patchiFacei] - origin_))
				& Sf.boundaryField()[patchi][patchiFacei];
		}
	}
	*/
}

template<class RhoFieldType>
void Foam::RelativeMotion::makeRelativeRhoFlux
(
 	const RhoFieldType& rho,
	Field<scalar>& phi,
	const label patchi
) const 
{
	if (!active_)
	{
		return;
	}

	const surfaceVectorField& Cf = mesh_.Cf();
	const surfaceVectorField& Sf = mesh_.Sf();

	const volVectorField& U = mesh_.thisDb().lookupObject<volVectorField>("U");

	// const vector Omega = Omega_->value(mesh_.time().timeOutputValue())*axis_;
	const polyPatch& pp = mesh_.boundaryMesh()[patchi]; 

	// Skip empty patch, 2020-02-12, Geon-Hong 
	if (pp.type() == "empty")
	{
		return;
	}

	forAll(pp, patchFacei)
	{
		const label facei = pp.start() + patchFacei;

		scalar phii = 0.0;

		if (rm_.faceType(facei) == 1)
		{
			phii -= 0.0;

			// meshPhi is assigned for moving mesh cases, 2020-02-26, Geon-Hong 
			if (mesh_.moving())
			{
				const surfaceScalarField& meshPhi = fvc::meshPhi(U);
				const surfaceScalarField::Boundary& bMeshPhi = meshPhi.boundaryField();

				phi[patchFacei] = bMeshPhi[patchi][patchFacei];

			}
			else
			{
				phi[patchFacei] = 0.0;
			}
		}
		else if (rm_.faceType(facei) == 2)
		{
			phii -= 
				rho[patchFacei]
				*rm_.velocity(Cf.boundaryField()[patchi][patchFacei]) 
				& Sf.boundaryField()[patchi][patchFacei];
		}

		if (lm_.faceType(facei) == 1)
		{
			phii = 0.0;

			// meshPhi is assigned for moving mesh cases, 2020-02-26, Geon-Hong 
			if (mesh_.moving())
			{
				const surfaceScalarField& meshPhi = fvc::meshPhi(U);
				const surfaceScalarField::Boundary& bMeshPhi = meshPhi.boundaryField();

				phi[patchFacei] = bMeshPhi[patchi][patchFacei];
			}
			else 
			{
				phi[patchFacei] = 0.0;
			}
		}
		else if (lm_.faceType(facei) == 2)
		{
			phii -=
				rho[patchFacei]
				* lm_.velocity()
				& Sf.boundaryField()[patchi][patchFacei];
		}

		phi[patchFacei] += phii;
	}

	/*
	// Included patches 
	forAll(includedFaces_[patchi], i)
	{
		label patchFacei = includedFaces_[patchi][i];

		phi[patchFacei] = 0.0;
	}

	//Excluded patches 
	forAll((excludedFaces_[patchi], i)
	{
		label patchFacei = excludedFaces_[patchi][i];

		phi[patchFacei] -= rho[patchFacei]
				* (Omega * (Cf.boundaryField()[patchi][patchFacei] - origin^)) 
				& Sf.boundaryField()[patchi][patchFacei];
	}
	*/
}

template<class RhoFieldType>
void Foam::RelativeMotion::makeAbsoluteRhoFlux
(
 	const RhoFieldType& rho,
	surfaceScalarField& phi
) const 
{
	if(!active_)
	{
		return;
	}

	const surfaceVectorField& Cf = mesh_.Cf();
	const surfaceVectorField& Sf = mesh_.Sf();

	// const vector Omega = omega_->value(mesh_.time().timeOutputValue())*axis_;
	
	const vectorField& Cfi = Cf;
	const vectorField& Sfi = Sf;
	scalarField& phii = phi.primitiveFieldRef();

	// Internal Faces
	forAll(phii, i)
	{
		label facei = i;

		vector rvel = Zero;
		vector lvel = Zero;

		if(rm_.faceType(facei) == 1)
		{
			rvel = rm_.velocity(Cfi[facei]);
		}

		if(lm_.faceType(facei) == 1)
		{
			lvel = lm_.velocity();
		}

		phii[facei] += rho[facei]*(rvel+lvel) & Sfi[facei];
	}

	surfaceScalarField::Boundary& phibf = phi.boundaryFieldRef();

	// Included patches
	const polyBoundaryMesh& patches = mesh_.boundaryMesh();

	forAll(patches, patchi)
	{
		const polyPatch& pp = mesh_.boundaryMesh()[patchi];

		forAll(pp, patchFacei)
		{
			const label facei = pp.start() + patchFacei;

			scalar phii = 0.0;

			if(rm_.faceType(facei) == 1 || rm_.faceType(facei) == 2)
			{
				phii +=
					rho[patchFacei]
					* rm_.velocity(Cf.boundaryField()[patchi][patchFacei])
					& Sf.boundaryField()[patchi][patchFacei];
			}

			if(lm_.faceType(facei) == 1 || lm_.faceType(facei) == 2)
			{
				phii +=
					rho[patchFacei]
					* lm_.velocity()
					& Sf.boundaryField()[patchi][patchFacei];
			}
			
			phibf[patchi][patchFacei] += phii;
		}
	}

	/*
	forAll(includedFaces_, patchi)
	{
		forAll(includedFaces_[patchi], i)
		{
			label patchFacei = includedFaces_[patchi][i];

			phibf[patchi][patchFacei] += 
				rho.boundaryField()[patchi][patchFacei]
				* ( Omega ^ ( Cf.boundaryField()[patchi][patchFacei] - origin_)) 
				& Sf.boundaryField()[patchi][patchFacei];
		}
	}

	// Excluded patches 
	forAll(excludedFaces_, patchi)
	{
		forAll(excludedFaces_[patchi] , i)
		{
			label patchFacei = excludedFaces_[patchi][i];

			phibf[patchi][patchFacei] += 
				rho.boundaryField()[patchi][patchFacei]
				★ (Omega ^ (Cf.boundaryField()[patchi][patchFacei] - origin_)) 
				& Sf.boundaryField()[patchi][patchFacei];
		}
	}
	*/
}

template<class Type>
void Foam::RelativeMotion::zero
(
 	GeometricField<Type, fvsPatchField, surfaceMesh>& phi
) const 
{
	if(!active_)
	{
		return;
	}

	/*
	Field<Type>& phii = phi.primitiveFieldRef();

	forAll(internalFaces_, i)
	{
		phii[internalFaces_[i]] = Zero;
	}
	*/

	phi.primitiveFieldRef() = Zero;

	typename GeometricField<Type, fvsPatchField, surfaceMesh>::Boundary& phibf = phi.boundaryFieldRef();

	const polyBoundaryMesh& patches = mesh_.boundaryMesh();

	forAll(patches, patchi)
	{
		const polyPatch& pp = mesh_.boundaryMesh()[patchi];

		if(pp.type() == "empty")
		{
			continue;
		}

		forAll(pp, patchFacei)
		{
			const label facei = pp.start() + patchFacei;

			scalar phii = phibf[patchi][patchFacei];

			if(rm_.faceType(facei) == 1 || rm_.faceType(facei) == 2)
			{
				phii = Zero;
			}

			if(lm_.faceType(facei) == 1 || lm_.faceType(facei) == 2)
			{
				phii = Zero;
			}

			phibf[patchi][patchFacei] = phii;
		}
	}

	/*
	forAll(includedFaces_, patchi)
	{
		forAll(includedFaces_[patchi], i)
		{
			phibf[patchi][includedFaces_[patchi][i]] = Zero;
		}
	}

	forAll(excludedFaces_, patchi)
	{
		forAll(excludedFaces_[patchi], i)
		{
			phibf[patchi][excludedFaces_[patchi][i]] = Zero;
		}
	}
	*/
}

//- zeroFilter from MRFZoneList 
template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::RelativeMotion::zeroFilter
(
 	const tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>& tphi
) const 
{
	if(!active_)
	{
		return tmp<surfaceScalarField>(tphi, true);
	}

	tmp<surfaceScalarField> zphi
	(
	 	New
		(
		 	tphi,
			"zeroFilter(" + tphi().name() + ')',
			tphi().dimensions(),
			true
		)
	);

	zero(zphi.ref());

	return zphi;
}

// * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * * //
//class linearMotion 
template<class T>
T Foam::RelativeMotion::linearMotion::lastIncrement(const List<T>& L) const 
{
	if(L.size() ==0)
	{
		return Zero;
	}
	else if (L.size() == 1)
	{
		return L.last();
	}
	else 
	{
		return L.last() - L[L.size()-2];
	}
}

// ************************************************************************** //

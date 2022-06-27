/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 AUTHOR,AFFILIATION
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
#include "faceSet.H"
#include "geometricOneField.H"
#include "syncTools.H"
#include "IFstream.H"

#define checkIfCellZoneAvail 	\
	if (rm_.cellZoneID() == -1) \
	{							\
		return;					\
	}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
	defineTypeNameAndDebug(RelativeMotion, 0);
}

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RelativeMotion::RelativeMotion
(
 	const fvMesh& mesh, 
	const dictionary& dict
) 
:
    dict_(dict),
	rm_(mesh, dict.subDict("rotating")), 
	lm_(mesh, dict.subDict("linear")),
	mesh_(mesh),
	active_(dict_.lookupOrDefault("active", true))
{
	Info<< dict_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

//Foam::RelativeMotion::~RelativeMotion()
//{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::RelativeMotion::addCoriolis
(
	const volVectorField& U, 
	volVectorField& ddtU
) const
{
	checkIfCellZoneAvail;

	const labelList& cells = mesh_.cellZones()[rm_.cellZoneID()]; 
	vectorField& ddtUc = ddtU.primitiveFieldRef();
	const vectorField& Uc = U;

	const vector Omega = rm_.Omega();

	forAll(cells, i)
	{
		label celli = cells[i];
		ddtUc[celli] += (Omega ^ Uc[celli]);
	}
}

void Foam::RelativeMotion::addCoriolis
(
 	fvVectorMatrix& UEqn, 
	const bool rhs
) const
{
	checkIfCellZoneAvail;

	const labelList& cells = mesh_.cellZones()[rm_.cellZoneID()];
	const scalarField& V = mesh_.V();
	vectorField& Usource = UEqn.source();
	const vectorField& U = UEqn.psi();

	const vector Omega = rm_.Omega();

	if (rhs)
	{
		forAll(cells, i)
		{
			label celli = cells[i];
			Usource[celli] += V[celli]*(Omega ^ U[celli]);
		}
	}
	else
	{
		forAll(cells, i)
		{
			label celli = cells[i];
			Usource[celli] -= V[celli]*(Omega ^ U[celli]);
		}
	}
}

void Foam::RelativeMotion::addCoriolis
(
 	const volScalarField& rho,
 	fvVectorMatrix& UEqn, 
	const bool rhs
) const
{
	checkIfCellZoneAvail;

	const labelList& cells = mesh_.cellZones()[rm_.cellZoneID()];
	const scalarField& V = mesh_.V();
	vectorField& Usource = UEqn.source();
	const vectorField& U = UEqn.psi();

	const vector Omega = rm_.Omega();

	if (rhs)
	{
		forAll(cells, i)
		{
			label celli = cells[i];
			Usource[celli] += V[celli]*rho[celli]*(Omega ^ U[celli]);
		}
	}
	else
	{
		forAll(cells, i)
		{
			label celli = cells[i];
			Usource[celli] -= V[celli]*rho[celli]*(Omega ^ U[celli]);
		}
	}
}

void Foam::RelativeMotion::makeRelative(volVectorField& U) const
{
	checkIfCellZoneAvail;

	const volVectorField& C = mesh_.C();

	const labelList& cells = mesh_.cellZones() [rm_.cellZoneID()];

	forAll(cells, i)
	{
		label celli = cells[i];
		U[celli] -= rm_.velocity(C[celli]);
	}

	forAll(mesh_.cells(), i)
	{
		label celli = cells[i];
		U[celli] -= lm_.velocity();
	}

	// Make relative on patches
	const polyBoundaryMesh& patches = mesh_.boundaryMesh();
	volVectorField::Boundary& Ubf = U.boundaryFieldRef();

	forAll(patches, patchi)
	{
		const polyPatch& pp = patches[patchi];

		// skip empty patch
		if (pp.type() == "empty")
		{
			continue;
		}

		forAll(pp, patchFacei)
		{
			const label facei = pp.start() + patchFacei;

			vector Ubfi = Zero;

			if (rm_.faceType(facei) == 1)
			{
				Ubfi = Zero;
				Ubf[patchi][patchFacei] = Zero;
			}
			else if (rm_.faceType(facei) == 2) 
			{
				Ubfi -= rm_.velocity(C.boundaryField()[patchi][patchFacei]);
			}

			if (lm_.faceType(facei) == 1)
			{
				Ubfi -= Zero;
				Ubf[patchi][patchFacei] = Zero;
			}
			else if (lm_.faceType(facei) == 2) 
			{
				Ubfi -= lm_.velocity();
			}

			Ubf[patchi][patchFacei] += Ubfi;
		}
	}
}


void Foam::RelativeMotion::makeRelative(surfaceScalarField& phi) const 
{
	makeRelativeRhoFlux(geometricOneField(), phi);
}


void Foam::RelativeMotion::makeRelative(FieldField<fvsPatchField, scalar>& phi) const
{
	makeRelativeRhoFlux(oneFieldField(), phi);
}


void Foam::RelativeMotion::makeRelative(Field<scalar>& phi, const label patchi) const
{
	makeRelativeRhoFlux(oneField(), phi, patchi);
}

void Foam::RelativeMotion::makeRelative
(
 	const surfaceScalarField& rho, 
	surfaceScalarField& phi
) const
{
	makeRelativeRhoFlux(rho, phi) ;
}


void Foam::RelativeMotion::makeAbsolute(volVectorField& U) const
{
	checkIfCellZoneAvail;

	const volVectorField& C = mesh_.C();

	const labelList& cells = mesh_.cellZones() [rm_.cellZoneID()];

	forAll(cells, i)
	{
		label celli = cells[i];
		U[celli] += rm_.velocity(C[celli]);
	}

	forAll(mesh_.cells(), celli)
	{
		U[celli] += lm_.velocity();
	}


	// Included patches
	const polyBoundaryMesh& patches = mesh_.boundaryMesh();
	volVectorField::Boundary& Ubf = U.boundaryFieldRef() ;

	forAll(patches, patchi)
	{
		const polyPatch& pp = patches [patchi] ;

		forAll(pp, patchFacei)
		{
			const label facei = pp.start() + patchFacei;

			vector Ubfi = Zero;

			if (rm_.faceType(facei) == 1)
			{
				Ubfi = rm_.velocity(C.boundaryField() [patchi] [patchFacei]);
				Ubf[patchi][patchFacei] = Zero;
			}
			else if (rm_.faceType(facei) == 2)
			{
				Ubfi = rm_.velocity(C.boundaryField() [patchi] [patchFacei]) ;
			}

			if (lm_.faceType(facei) == 1)
			{
				Ubfi += lm_.velocity();
				Ubf [patchi] [patchFacei] = Zero;
			}
			else if (lm_.faceType(facei) == 2)
			{
				Ubfi += lm_.velocity();
			}

			Ubf[patchi][patchFacei] += Ubfi;
		}
	}
}

void Foam::RelativeMotion::makeAbsolute(surfaceScalarField& phi) const
{
	makeAbsoluteRhoFlux(geometricOneField(), phi);
}

void Foam::RelativeMotion::makeAbsolute
(
	const surfaceScalarField& rho,
	surfaceScalarField& phi
) const
{
	makeAbsoluteRhoFlux(rho, phi);
}

void Foam::RelativeMotion::correctBoundaryVelocity(volVectorField& U) const
{
	if (!active_)
	{
		return;
	}
	// const vector Omega = rm_.Omega();
	// Included patches
	const polyBoundaryMesh& patches = mesh_.boundaryMesh();
	volVectorField::Boundary& Ubf = U.boundaryFieldRef();

	// Additional references and fields for the mesh motion, 2020-02-10, Geon-Hong 
	const pointField& oldPoints = mesh_.oldPoints();

	forAll(patches, patchi)
	{
		const polyPatch& pp = patches[patchi];

		if (pp.type() == "empty" || pp.type() == "processor")
		{
			continue;
		}
		const vectorField& patchC = mesh_.Cf().boundaryField()[patchi];
		vectorField pfld(Ubf[patchi]);
		// Additional references and field for the mesh motion, 
		// 2020-02-10, Geon-Hong 
		// const fvPatch& p = mesh_.boundary()[patchi];
		const fvPatch& p = mesh_.boundary()[patchi] ;
		scalarField phip(p.size());

		if	(mesh_.moving())
		{
			phip =
			(
			 	p.patchField<surfaceScalarField, scalar>(fvc::meshPhi(U))
			);

			// Set velocity on patches due to mesh motion 
			forAll(pp, patchFacei)
			{
				const vector oldFc = pp[patchFacei].centre(oldPoints);
				const scalar deltaT = mesh_.time().deltaTValue();
				const vector Up((pp.faceCentres()[patchFacei] - oldFc)/deltaT);

				const vector n(pp.faceNormals()[patchFacei]);
				const scalar magSf(pp.magFaceAreas()[patchFacei]);
				scalar Un = phip[patchFacei]/(magSf + VSMALL);

				pfld[patchFacei] = Up + n*(Un - (n & Up));
			}
		}

		forAll(pp, patchFacei)
		{
			const label facei = pp.start() + patchFacei;
			vector& pfldi = pfld[patchFacei];
			if (rm_.faceType(facei) == 1 || lm_.faceType(facei) == 1)
			{ 
				// Take account of mesh motion, 2020-02-10, Geon-Hong
				// pfldi = Zero;
				if (mesh_. moving ())
				{
					/*
					const vector oldFc = pp[patchFacei].centre(oldPoints);
					const scalar deltaT = mesh_.time().deltaTValueO;
					const vector Up((pp.faceCentres()[patchFacei] - oldFc》/deltaT);

					const vector n(pp.faceNormals()[patchFacei]);
					const scalar magSf(pp.magFaceAreas()[patchFacei]);
					scalar Un = phip[patchFacei]/(magSf + VSMALL);

					pfldi = Up + n*(Un - (n & Up));
					*/
					// Do not modify the velocity on this patch
				}
				else
				{
					pfldi = Zero;
				}
			}

			if (rm_.faceType(facei) == 1)
			{
				pfldi += rm_.velocity(patchC[patchFacei]);
			}

			if (lm_.faceType(facei) == 1)
			{
				pfldi += lm_.velocity();
			}

			/*
			if (pp.name() == "blade")
			{
				Pout << "Patch " << patchi
					<< ", face " << patchFacei
					<< ", pfldi = " << pfldi
					<< endl;
			}
			*/
		}
		Ubf[patchi] == pfld;
	}
}


//----------------------------------------------------------------------------
//- Functions from MRFZoneList

Foam::tmp<Foam::volVectorField> Foam::RelativeMotion::DDt
(
 	const	volVectorField& U
) const
{
	tmp<volVectorField> tacceleration 
	(
	 	new volVectorField
		(
		 	IOobject
			(
			 	"MRFZoneList:acceleration",
				U.mesh().time().timeName(),
				U.mesh()
			),
			U.mesh(),
			dimensionedVector (U.dimensions() /dimTime, Zero)
		)
	);
	volVectorField& acceleration = tacceleration.ref();
	addCoriolis(U, acceleration);
	return tacceleration;
}

Foam::tmp<Foam::volVectorField> Foam::RelativeMotion::DDt
(
 	const volScalarField& rho, 
	const volVectorField& U
) const
{
	return	rho*DDt(U);
}

//- Return the given absolute boundary flux relative within the
// RelativeMotion region
Foam::tmp<Foam::surfaceScalarField> Foam::RelativeMotion::relative
(
	const tmp<surfaceScalarField>& tphi
) const
{
	if (!active_)
	{
		return tmp<surfaceScalarField>(tphi, true);
	}

	tmp<surfaceScalarField> rphi
	(
		 New
		(
			tphi,
			"relative(" + tphi().name() + ')',
			tphi().dimensions(),
			true
		)
	);

	makeRelative(rphi.ref ());
	tphi.clear();
	return rphi;
}

Foam::tmp<Foam::FieldField<Foam::fvsPatchField, Foam::scalar>>
Foam::RelativeMotion::relative
(
 	const	tmp<FieldField<fvsPatchField, scalar>>& tphi
)	const
{
	if (!active_)
	{
		return tmp<FieldField<fvsPatchField, scalar>>(tphi, true);
	}

	tmp<FieldField<fvsPatchField, scalar>> rphi(New(tphi, true));
	makeRelative(rphi.ref());
	tphi.clear();
	return rphi;
}

Foam::tmp<Foam::Field<Foam::scalar>>
Foam::RelativeMotion::relative
(
 	const	tmp<Field<scalar>>& tphi,
	const	label patchi
)	const
{
	if (!active_)
	{
		return	tmp<Field<scalar>>(tphi, true);
	}
	tmp<Field<scalar>> rphi(New(tphi, true));
	makeRelative(rphi.ref(), patchi);
	tphi.clear();
	return	rphi;
}

//----------------------------------------------------------------------------
//- Print motion information
void Foam::RelativeMotion::printMotionInfo()
{
	if (rm_.active())
	{
		rm_.printMotionInfo();
	}

	if (lm_.active())
	{
		lm_.printMotionInfo();
	}
}

//- Check if initialization is necessary
bool Foam::RelativeMotion::needInitialize(const volVectorField& U)
{
	scalar magU(max(mag(U)).value());

	reduce(magU, sumOp<scalar>());

	if (magU > 1e-10)
	{
		return true;
	}

	return	false;
}

//----------------------------------------------------------------------------
//- Debugging
void Foam::RelativeMotion::checkBoundaryVelocity
(
 	const	word& bName,
	const	string& tag
)	const
{
	const polyBoundaryMesh&	patches = mesh_.boundaryMesh();
	const volVectorField& U(mesh_.objectRegistry::thisDb().lookupObject<volVectorField>("U"));
	const volVectorField::Boundary& Ubf = U.boundaryField();

	forAll(patches, patchi)
	{
		const polyPatch& pp = patches[patchi];

		const vectorField& pfld(Ubf[patchi]);

		if(pp.name()==bName)
		{
			const scalar maxUp(gMax(mag(pfld)));

			Info <<"Patch " << patchi
				<< "(" << pp.name() << ")"
				<< ", " << tag 
				<< ", max(U) = " << maxUp
				<<endl;
		}
	}
}

//---------------------------------------------------------------------------
//- I/O

bool Foam::RelativeMotion::writeData(Ostream& os) const
{
	rm_.writeData(os);
	lm_.writeData(os);

	return true;
}

bool Foam::RelativeMotion::read(const dictionary& dict)
{
	return (rm_.read() && lm_.read());
}

void Foam::RelativeMotion::update()
{
	if(mesh_.topoChanging())
	{
		rm_.update();
		lm_.update();
	}
}

#include "rotatingMotion.C"
#include "linearMotion.C"

//***************************************************************************//



// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //
/*
void Foam::RelativeMotion::operator=(const RelativeMotion& rhs)
{
    if (this == &rhs)
    {
        return;  // Self-assignment is a no-op
    }
}
*/

// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// ************************************************************************* //

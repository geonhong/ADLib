/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Released 2004-2011 OpenCFD Ltd. 
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

#include "movingRMWallVelocityFvPatchVectorField.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::movingRMWallVelocityFvPatchVectorField::
movingRMWallVelocityFvPatchVectorField
(
 	const fvPatch& p,
	const DimensionedField<vector, volMesh>& iF
) 
:
    fixedValueFvPatchVectorField(p, iF),
	U0_(*this),
	latestTimeIndex_(0)
{}

Foam::movingRMWallVelocityFvPatchVectorField::
movingRMWallVelocityFvPatchVectorField
(
 	const fvPatch& p,
	const DimensionedField<vector, volMesh>& iF,
	const dictionary& dict
)
:
	fixedValueFvPatchVectorField(p, iF),
	U0_(*this),
	latestTimeIndex_(0)
{}

Foam::movingRMWallVelocityFvPatchVectorField::
movingRMWallVelocityFvPatchVectorField
(
 	const movingRMWallVelocityFvPatchVectorField& ptf,
	const fvPatch& p,
	const DimensionedField<vector, volMesh>& iF,
	const fvPatchFieldMapper& mapper
)
:
	fixedValueFvPatchVectorField(ptf, p, iF, mapper),
	U0_(*this),
	latestTimeIndex_(0)
{}

Foam::movingRMWallVelocityFvPatchVectorField::
movingRMWallVelocityFvPatchVectorField
(
 	const movingRMWallVelocityFvPatchVectorField& mwvpvf
)
:
	fixedValueFvPatchVectorField(mwvpvf),
	 U0_(*this),
	 latestTimeIndex_(0)
{}

Foam::movingRMWallVelocityFvPatchVectorField::
movingRMWallVelocityFvPatchVectorField
(
 	const movingRMWallVelocityFvPatchVectorField& mwvpvf,
	const DimensionedField<vector, volMesh>& iF
)
:
	fixedValueFvPatchVectorField(mwvpvf, iF),
	 U0_(*this),
	 latestTimeIndex_(0)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::movingRMWallVelocityFvPatchVectorField::updateCoeffs()
{
	if (updated())
	{
		return;
	}

	/*
	const fvMesh& mesh = internalField().mesh();

	if(mesh.moving())
	{
		// Store current velocity field if time has been updated 
		if (latestTimelndex_ != mesh.time().timeIndex())
		{
			U0_.operator=(*this);
			latestTimelndex_ = mesh.time().timeIndex();
		}

		const fvPatch& p = patch();
		const polyPatch& pp = p.patch();
		const pontField& oldPoints = mesh.oldPoints();

		vectorField oldFc(pp.size());

		forAll(oldFc, i)
		{
			oldFc[i] = pp[i].centre(oldPoints);
		}

		const scalar deltaT = mesh.time().deltaTValue();

		const vectorField Up((pp.faceCentres() - oldFc)/deltaT);

		const volVectorField& U =
			static_cast<const volVectorField&>(internalField());

		scalarField phip
		(
			p.patchField<surfaceScalarField, scalar>(fvc::meshPhi(U))
		);

		const vectorField n(p.nf());
		const scalarField& magSf = p.magSf();
		tmp<scalarField> Un = phip/(magSf + VSMALL);

		// vectorField::operator=(U0_ + Up + n*(Un - (n & Up)));
		vectorField::operator=(Up + n * (Un - (n & Up)));
	}
	*/

	fixedValueFvPatchVectorField::updateCoeffs();
}

void Foam::movingRMWallVelocityFvPatchVectorField::write(Ostream& os) const 
{
	fvPatchVectorField::write(os);
	writeEntry("value", os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
	makePatchTypeField
	(
	 	fvPatchVectorField,
		movingRMWallVelocityFvPatchVectorField
	);
}

//***************************************************************************//

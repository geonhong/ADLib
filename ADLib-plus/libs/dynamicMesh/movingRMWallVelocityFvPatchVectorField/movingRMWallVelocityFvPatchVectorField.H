/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
	Released 2004-2011 OpenCFD Ltd.
	Copyright (C) 2011-2016 OpenFOAM Foundation 
	Modified code Copyright (C) 2015 OpenCFD Ltd.
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

Class
    Foam::movingRMWallVelocityFvPatchVectorField

Group
	grpWallBoundaryConditions

Description
	This boundary condition provides a velocity condition for cases with 
	moving walls.

Usage
	Exampes of the boundary condition specification:
	\verbatim 
	<patchName>
	{
		type	movingRMWallVelocity;
		value	uniform (0 0 0); // Initial value
	}
	\endverbatim

See also 
	Foam::fixedValueFvPatchVectorField
	Foam::MRFZone

SourceFiles
    movingRMWallVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef movingRMWallVelocityFvPatchVectorField_H
#define movingRMWallVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class movingRMWallVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class movingRMWallVelocityFvPatchVectorField
:
	public fixedValueFvPatchVectorField
{
	
	// Store previous patch field 
		vectorField U0_;

	// Time index for updating previous patch
		label latestTimeIndex_;

public:

	//- Runtime type information 
	TypeName("movingRMWallVelocity");


	//- Construct from path and internal field 
	movingRMWallVelocityFvPatchVectorField
	(
	 	const fvPatch&,
		const DimensionedField<vector, volMesh>&
	);

	//- Construct from patch and internal field and dictionary 
	movingRMWallVelocityFvPatchVectorField
	(
	 	const fvPatch&,
		const DimensionedField<vector, volMesh>&,
		const dictionary&
	);

	//- Construct by mapping given movingRMWallVelocityFvPatchVectorField
	//	onto a new patch
	movingRMWallVelocityFvPatchVectorField
	(
	 	const movingRMWallVelocityFvPatchVectorField&,
		const fvPatch&,
		const DimensionedField<vector, volMesh>&,
		const fvPatchFieldMapper&
	);

	//- Construct as copy
	movingRMWallVelocityFvPatchVectorField
	(
	 	const movingRMWallVelocityFvPatchVectorField&
	);

	//- Construct and return a clone 
	virtual tmp<fvPatchVectorField> clone() const
	{
		return tmp<fvPatchVectorField>
		(
		 	new movingRMWallVelocityFvPatchVectorField(*this)
		);
	}

	//- Construct as copy setting internal field reference 
	movingRMWallVelocityFvPatchVectorField
	(
	 	const movingRMWallVelocityFvPatchVectorField&,
		const DimensionedField<vector, volMesh>&
	);

	//- Construct and return a clone setting internal field reference 
	virtual tmp<fvPatchVectorField> clone
	(
	 	const DimensionedField<vector, volMesh>& iF
	) const 
	{
		return tmp<fvPatchVectorField>
		(
		 	new movingRMWallVelocityFvPatchVectorField(*this, iF)
		);
	}

	// Member Functions

		//- Update the coefficients associated with the patch field 
		virtual void updateCoeffs();

		//- Write 
		virtual void write(Ostream&) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

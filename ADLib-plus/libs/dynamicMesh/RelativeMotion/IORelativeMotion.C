/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
	Modified code Copyright (C) 2019 OpenCFD Ltd.
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

#include "IORelativeMotion.H"
#include "fvMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::IOobject Foam::IORelativeMotion::createIOobject
(
 	const fvMesh& mesh
) const
{
	IOobject io
	(
	 	"RelativeMotionProperties",
		mesh.time().constant(),
		mesh,
		IOobject::MUST_READ,
		IOobject::NO_WRITE
	);

	if (io.typeHeaderOk<IOdictionary>(true))
	{
		Info << "Creating RelativeMotion zone from " << io.name() << endl;

		io.readOpt() = IOobject::MUST_READ_IF_MODIFIED;
		return io;
	}
	else
	{
		Info << "No RelativeMotion models present" << nl << endl;

		io.readOpt() = IOobject::NO_READ;
		return io;
	}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::IORelativeMotion::IORelativeMotion
(
 	const fvMesh& mesh 
)  
:
    IOdictionary(createIOobject(mesh)),
	RelativeMotion(mesh, *this)
{}

bool Foam::IORelativeMotion::read()
{
	if (regIOobject::read())
	{
		RelativeMotion::read(*this);
		return true;
	}
	return false;
}

//***************************************************************************//

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012 AUTHOR,AFFILIATION
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
    Foam::IORelativeMotion

Description
	List of RelativeMotion zones with IO functionality. RelativeMotion zones are specified by a list 
	of dictionary entries, e.g.

	\verbatim 
	zone1
	{
		cellZone	rotor1;
		active		yes;
		...
	}

	zone2
	{
		cellZone	rotor2;
		active		yes;
		...
	}
	\endverbatim

SourceFiles
    IORelativeMotion.C

\*---------------------------------------------------------------------------*/

#ifndef IORelativeMotion_H
#define IORelativeMotion_H

#include "IOdictionary.H"
#include "RelativeMotion.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class IORelativeMotion Declaration
\*---------------------------------------------------------------------------*/

class IORelativeMotion
:
	public IOdictionary,
	public RelativeMotion
{
private:

	// Private Member Functions
	
		//- Create IO object if dictionary is present 
		IOobject createIOobject(const fvMesh& mesh) const;

		//- No copy construct 
		IORelativeMotion(const IORelativeMotion&) = delete; 

		//- No copy assignment 
		void operator=(const IORelativeMotion&) = delete;

public:

	// Constructors 
		
		//- Construct from mesh 
		IORelativeMotion(const fvMesh& mesh);

		//- Destructor 
		virtual ~IORelativeMotion()
		{}

	// Member Functions

		//- Read dictionary
		virtual bool read();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

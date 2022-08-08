/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::Function1Types::scaledLinearRamp

Description
    Linear ramp function starting from 0 and increasing linearly to 1 from
    \c start over the \c duration and remaining at 1 thereafter.

See also
    Foam::Function1Types::ramp

SourceFiles
    scaledLinearRamp.C

\*---------------------------------------------------------------------------*/

#ifndef scaledLinearRamp_H
#define scaledLinearRamp_H

#include "ramp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{

/*---------------------------------------------------------------------------*\
                           Class scaledLinearRamp Declaration
\*---------------------------------------------------------------------------*/

class scaledLinearRamp
:
    public ramp
{

	// Target value
	// Take scalar into account only since the class ramp inherits
	// Function1<scalar>
	scalar value_;

public:

    // Runtime type information
    TypeName("scaledLinearRamp");


    // Constructors

        //- Construct from entry name and dictionary
        scaledLinearRamp(const word& entryName, const dictionary& dict);


    //- Destructor
    virtual ~scaledLinearRamp() = default;


    // Member Functions

        //- Return value for time t
        virtual inline scalar value(const scalar t) const;

		//- Integrate between two values
		virtual inline scalar integrate(const scalar x1, const scalar x2) const;

		//- Write in dictionary format
		virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "scaledLinearRampI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

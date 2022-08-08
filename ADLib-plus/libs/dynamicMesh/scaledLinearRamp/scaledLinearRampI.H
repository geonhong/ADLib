/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
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

#include "scaledLinearRamp.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::Function1Types::scaledLinearRamp::value
(
    const scalar t
) const
{
    return value_*ramp::linearRamp(t);
}

inline Foam::scalar Foam::Function1Types::scaledLinearRamp::integrate
(
    const scalar x0,
	const scalar x1
) const
{
	scalar ts = max(x0, this->start_);
	scalar te = x1;
	scalar tr = this->start_ + this->duration_;
	scalar ys = value(ts);
	scalar ye = value(x1);

	if (x1 < tr)
	{
		return 0.5*(ye+ys)*(te-ts);
	}
	else
	{
		scalar tramp = tr - ts;
		scalar trest = x1 - tr;

		return 0.5*(value_+ys)*tramp + value_*trest;
	}
}

// ************************************************************************* //

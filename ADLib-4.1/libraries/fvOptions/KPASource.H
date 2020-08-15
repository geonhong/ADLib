/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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

Class
    Foam::KPASource

Description
	KPA source for self-propulsion simulations

	Evaluate the momentum sources from the KPA14 code

Usage
	Example usage:
	\verbatim
	KPASourceCoeffs
	{
		fields		(U);		// names of fields to apply source
		rps			5.4;		// rotational speed of a propeller
		velocity	1.17;		// speed of the objective ship
		density		999.05;		// density of water
		scaleRatio	43.3333;	// scale ratio of the objective ship
		propName	KP505;		// propeller name
		interval	500;		// interval for running KPA14 code

		propAxis				// local coordinate system for
		{						// representing the propeller axis
			origin	(3.5 0 -0.3);
			coordinateRotation
			{
				type	axes;
				e1		(1 0 0);
				e3 		(0 0 1);
			}
		}
	}
	\endverbatim

SourceFiles
    KPASource.C
	KPASourceFwd.H

\*---------------------------------------------------------------------------*/

#ifndef KPASource_H
#define KPASource_H

#include "cellSetOption.H"
#include "coordinateSystem.H"
#include "KPASourceFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                         Class KPASource Declaration
\*---------------------------------------------------------------------------*/

class KPASource
:
    public cellSetOption
{

protected:

	// Protected data

		//- RPS
		scalar rps_;

		//- Speed of a ship in model scale
		scalar VM_;

		//- Density
		scalar rho_;

		//- Scale ratio
		scalar SR_;

		//- Diameter of a propeller
		scalar dia_;

		//- Hub radius
		scalar rh_;

		//- Propeller name and corresponding file names
		word propName_;
		word propGeo_;
		word propMt_;

		//- Update interval
		label interval_;

		//- Momentum source
		vectorField source_;

		//- Coordinate system for the propeller axis
		const coordinateSystem propAxis_;

		//- Points field for extracting velocity
		pointField probes_;

		//- Elements to probe data
		labelList elemList_;

		//- Cells of rotate region in local coordinates
		pointField rotCfld_;

		//- Stencil for interpolating forces estimated from KPA14
		labelListList *stencil_;

		//- The number of updates
		label nUpdate_;

private:

    // Private data

        //- Size of wake array
        const label nR_;
		const label nQ_;

		//- Wake extract data
		wake w_;

		//- KPA data
		KPAData KPA_;

		//- pi
		const scalar& pi;


    // Private Member Functions

		//- Update momentum source terms
		void updateMomentumSources();

		//- Extract velocities at probe points
		void extractWakeVelocity();

		//- Run KPA code
		void runKPA();

		//- Interpolate force data
		vector interpolateData
		(
			const labelList& stencil,
			const label& itrg
		);

		vector interpolateData
		(
			const labelList& stencil,
			const vectorField& fld,
			const point& Ctrg
		);

		//- Check if a given point is located within a cell
		bool pointInCell(const point& p, label celli) const;

		//- Make probe points
		void getProbePoints();

		//- Find element cell for probing data
		label findCell(const labelList&, const point&) const;

		//- Get propeller geometry data from geoFile
		bool getPropGeoData();

        //- Disallow default bitwise copy construct
        KPASource(const KPASource&);

        //- Disallow default bitwise assignment
        void operator=(const KPASource&);


public:

	//- Runtime type information
	TypeName("KPASource");


    // Constructors

        //- Construct from components
        KPASource
		(
			const word& name,
			const word& modelType,
			const dictionary& dict,
			const fvMesh& mesh
		);


    //- Destructor
    virtual ~KPASource()
	{}


    // Member Functions

        // Access

			//- Hub radius
			scalar rh(bool dimless = false) const
			{
				if (dimless)
				{
					return rh_;
				}
				else
				{
					return 0.5*rh_*dia_;
				}
			}

		
		// Add explicit and implicit contributions

			//- Source term to momentum equation
			virtual void addSup
			(
				fvMatrix<vector>& eqn,
				const label fieldi
			);

			//- Source term to compressible momentum equation
			virtual void addSup
			(
				const volScalarField& rho,
				fvMatrix<vector>& eqn,
				const label fieldi
			);


		// IO

			//- Read dictionary
			virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
	#include "KPASourceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

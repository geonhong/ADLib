/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 AUTHOR,AFFILIATION
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
    Foam::RelativeMotion

Description
	RelativeMotion zone definition based on cell zone and parameters
	obtained from a control dictionary constructed from the given stream.

	The rotation of the RelativeMotion region is defined by an origin and axis 
	of rotation and an angular speed.

SourceFiles
    RelativeMotionI.H
    RelativeMotion.C
    RelativeMotionIO.C

\*---------------------------------------------------------------------------*/

#ifndef RelativeMotion_H
#define RelativeMotion_H

#include "dictionary.H"
#include "wordList.H"
#include "labelList.H"
#include "dimensionedScalar.H"
#include "dimensionedVector.H"
#include "volFieldsFwd.H"
#include "surfaceFields.H"
#include "fvMatricesFwd.H"
#include "mapPolyMesh.H"
#include "Function1.H"
#include "autoPtr.H"
#include "forces.H"
//#include "forcesExt.H"
#include "regIOobject.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class fvMesh;


/*---------------------------------------------------------------------------*\
                         Class RelativeMotion Declaration
\*---------------------------------------------------------------------------*/

class RelativeMotion
{
	// Private classes
	
	//- Class for handling rotating motion
	class rotatingMotion
	{
		// Private Data
		
			//- Reference to the mesh database
			const fvMesh& mesh_;

			//- MRF dictionary
			dictionary dict_;

			//- RelativeMotion region active flag
			bool active_;

			//- Name of cell zone
			word cellZoneName_;

			//- Cell zone ID
			label cellZoneID_;

			//- Patch names to exclude
			const wordRes excludedPatchNames_;

			//- Excluded patch label
			labelList excludedPatchLabels_;

			//- Internal faces that are part of RelativeMotion
			labelList faceType_;

			//- Origin of the axis
			const vector origin_;

			//- Axis vector
			vector axis_;

			//- Angular velocity (rad/sec)
			autoPtr<Function1<scalar>> omega_;

		// Private member functions

			//- DIvide faces in frame according to patch
			void setRelativeMotionFaces();

			//- Nocopy construct
			rotatingMotion(const rotatingMotion&) = delete;

			//- No copy assignment
			void operator=(const rotatingMotion&) = delete;

	public:

		// Declare name of the class and its debug switch
		ClassName("rotatingMotion");

		// Constructors

			rotatingMotion(const fvMesh&, const dictionary&);

		// Member functions

			//- Update the mesh corresponding to given map
			void updateMesh(const mapPolyMesh& mpm)
			{
				// Only updates face addressing
				setRelativeMotionFaces();
			}

			//- Return const access to the RelativeMotion active flag
			bool active() const;

			//- Return cellZoneID
			label cellZoneID() const;

			//- Return the current Omega vector
			vector Omega() const;

			//- Return the current velocity
			vector velocity(const point&) const;

			//- Return the face type of the given face
			label faceType(const label&) const;

			//- Update RelativeMotion faces if the mesh topology changes
			void update();

			//- Print rotating motion info
			void printMotionInfo();

		// I-O

			//- Write
			void writeData(Ostream& os) const;

			//- Read RelativeMotion dictionary
			bool read();
	};

	//- Class for handling the linear motion
	class linearMotion
	{

		// Private Data

			//- Reference to the mesh database
			const fvMesh& mesh_;

			//- MRF dictionary
			dictionary dict_;

			const wordRes excludedPatchNames_;

			labelList excludedPatchLabels_;

			//- RelativeMotion region active flag
			bool active_;

			//- Internal faces that are part of RelativeMotion
			labelList faceType_;

			//- forces
			functionObjects::forces fm_;

			//- Velocities
			vector velocity_;		// target velocity (velocity input)
			List<vector> vList_;	// velocity list
			List<scalar> tList_;	// time list
			List<vector> pList_;	// position list

			//- Mass
			scalar mass_;

			//- Ramp time for relaxing the motion at initial stage
			scalar rampTime_;

			//- Clamp time for forcing the linear motion
			//  The motion is released freely after the clamp time
			scalar clampTime_;

			//- Coefficient for accelerating the motion prediction
			scalar motionPredFactor_;
			scalar maxMotionStepSize_;

			//- Target force of trimmed point
			vector targetForce_;

		// Private member functions

			//- Divide faces in frame according to patch
			void setRelativeMotionFaces();

			//- Get the ramp time
			scalar rampTime() const;

			//- Get the latest increment of the given list
			template<class T>
			T lastIncrement(const List<T>& L) const;

			//- Append data list
			void appendData(const scalar&, const vector&);
			void appendData(const scalar&, const vector&, const vector&);

			//- Update motion information
			void updateMotionInfo();

			//- No copy construct
			linearMotion(const linearMotion&) = delete;

			//- No copy assignment
			//	linearMotion(const linearMotion&) = delete;
			void operator=(const linearMotion&) = delete;

	public:

		// Declare name of the class and its debug switch
			ClassName("linearMotion");

		// Constructors
			linearMotion(const fvMesh&, const dictionary&);

		// Member functions

			//- Update the mesh corresponding to given map
			void updateMesh(const mapPolyMesh& mpm)
			{
				// Only updates face addressing
				setRelativeMotionFaces();
			}

			//- Return const access to the RelativeMotion active flag;
			bool active() const;

			//- Return the current velocity
			vector velocity() const;

			//- Return the face type of the given face
			label faceType(const label&) const;

			//- Update RelativeMotion faces if the mesh topology changes
			void update();

			//- Print linear motion info
			void printMotionInfo();

		// I-O

			//- Write
			void writeData(Ostream& os) const;

			//- Read RelativeMotion dictionary
			bool read();
	};

	// Private data
	
		//- Dictionary
		dictionary dict_;

		//- Rotating motion
		rotatingMotion rm_;

		//- Whole domain data (linear motion)
		linearMotion lm_;

		//- Reference to the mesh database
		const fvMesh& mesh_;

		//- RelativeMotion region active flag
		bool active_;


    // Private Member Functions

		//- Make the given absolute mass/vol flux relative within the RelativeMotion region
		template<class RhoFieldType>
		void makeRelativeRhoFlux
		(
		 	const RhoFieldType& rho,
			surfaceScalarField& phi
		) const;

		//- Make the given absolute mass/vol flux relative within the RelativeMotion region
		template<class RhoFieldType>
		void makeRelativeRhoFlux
		(
		 	const RhoFieldType& rho,
			FieldField<fvsPatchField, scalar>& phi
		) const;

		//- Make the given absolute mass/vol flux relative within the RelativeMotion region
		template<class RhoFieldType>
		void makeRelativeRhoFlux
		(
		 	const RhoFieldType& rho,
			Field<scalar>& phi,
			const label patchi
		) const;

		//- Make the given relative mass/vol flux absolute within the RelativeMotion region
		template<class RhoFieldType>
		void makeAbsoluteRhoFlux
		(
		 	const RhoFieldType& rho,
			surfaceScalarField& phi
		) const;

        //- No copy construct
        RelativeMotion(const RelativeMotion&) = delete;

        //- No copy assignment
        void operator=(const RelativeMotion&) = delete;


public:

	// Declare name of the class and its debug switch
	ClassName ("RelativeMotion");


    // Constructors

        //- Construct null
        RelativeMotion();

        //- Construct from fvMesh
        RelativeMotion
		(
		 	const fvMesh& mesh,
			const dictionary& dict
		);

		//- Return clone
		autoPtr<RelativeMotion> clone() const
		{
			NotImplemented;
			return nullptr;
		}


    // Member Functions

		//- Return const access to the RelativeMotion active flag
		inline bool active() const;

		//- Update the mesh corresponding to given map
		void updateMesh(const mapPolyMesh& mpm)
		{
			rm_.updateMesh(mpm);
			lm_.updateMesh(mpm);
		}

		//- Add the Coriolis force contribution to the acceleration field
		void addCoriolis
		(
		 	const volVectorField& U,
			volVectorField& ddtU
		) const;

		//- Add the Coriolis force contribution to the momentum equation
		//  Adds to the lhs of the equation; optionally add to rhs
		void addCoriolis
		(
		 	fvVectorMatrix& UEqn,
			const bool rhs = false
		) const;

		//- Add the Coriolis force contribution to the momentum equation
		//  Adds to the lhs of the equation; optionally add to rhs
		void addCoriolis
		(
		 	const volScalarField& rho,
			fvVectorMatrix& UEqn,
			const bool rhs = false
		) const;

		//- Make the given absolute velocity relative within the RelativeMotion region
		void makeRelative(volVectorField& U) const;

		//- Make the given absolute flux relative within the RelativeMotion region
		void makeRelative(surfaceScalarField& phi) const;

		//- Make the given absolute boundary flux relative
		//  within the RelativeMotion region
		void makeRelative(FieldField<fvsPatchField, scalar>& phi) const;

		//- Make the given absolute patch flux relative
		//  within the RelativeMotion region
		void makeRelative(Field<scalar>& phi, const label patchi) const;

		//- Make the given absolute mass-flux relative within the RelativeMotion region
		void makeRelative
		(
		 	const surfaceScalarField& rho,
			surfaceScalarField& phi
		) const;

		//- Make the given relative velocity absolute within the RelativeMotion region
		void makeAbsolute(volVectorField& U) const;

		//- Make the given relative flux absolute within the RelativeMotion region
		void makeAbsolute(surfaceScalarField& phi) const;

		//- Make the given relative mass-flux absolute within the RelativeMotion region
		void makeAbsolute
		(
		 	const surfaceScalarField& rho,
			surfaceScalarField& phi
		) const;

		//- Correct the boundary velocity for the rotation of the RelativeMotion region
		void correctBoundaryVelocity(volVectorField& U) const;

		//- Zero the RelativeMotion region of the given field
		template<class Type>
		void zero(GeometricField<Type, fvsPatchField, surfaceMesh>& phi) const;

		//- Update RelativeMotion faces if the mesh topology changes
		void update();

		//----------------------------------------------------------------------------
		//- Functions from MRFZoneList

		//- DDt, zeroFilter, and relative methods
		tmp<volVectorField> DDt
		(
		 	const volVectorField& U
		) const;

		tmp<volVectorField> DDt
		(
		 	const volScalarField& rho,
			const volVectorField& U
		) const;

		template<class Type>
		tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> zeroFilter
		(
		 	const tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>& tphi
		) const;

		//- Return the given absolute flux relative within the RelativeMotion region
		tmp<surfaceScalarField> relative
		(
		 	const tmp<surfaceScalarField>& phi
		) const;

		tmp<FieldField<fvsPatchField, scalar>> relative
		(
		 	const tmp<FieldField<fvsPatchField, scalar>>& tphi
		) const;

		tmp<Field<scalar>> relative
		(
		 	const tmp<Field<scalar>>& tphi,
			const label patchi
		) const;

		//----------------------------------------------------------------------------
		//- Print motion information
		void printMotionInfo();

		//- Check if initialization is necessary
		bool needInitialize(const volVectorField&);

		//----------------------------------------------------------------------------
		//- For debugging
		void checkBoundaryVelocity
		(
		 	const word& pName,
			const string& tag = ""
		) const;

	
	// I-O
	
		//- Write
		bool writeData(Ostream& os) const;

		//- Read RelativeMotion dictionary
		bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
	#include "RelativeMotionTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "RelativeMotionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

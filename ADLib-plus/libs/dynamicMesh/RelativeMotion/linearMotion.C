/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
	Released 2004-2011 OpenCFD Ltd. 
    Copyright (C) 2011-2017 OpenFOAM Foundation
	Modified code Copyright (C) 2018 OpenCFD Ltd.
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

#include "mathematicalConstants.H"
#include "gravityMeshObject.H"

// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * * //

void Foam::RelativeMotion::linearMotion::setRelativeMotionFaces()
{
	const polyBoundaryMesh& patches = mesh_.boundaryMesh();

	// Type per face;
	// 0: not in zone
	// 1: moving with frame
	// 2: other 
	// labelList faceType_(mesh_.nFaces(), Zero);
	
	// Determine faces in cell zone 
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// (without constructing cells)
	
	for(label facei = 0; facei<mesh_.nInternalFaces(); facei++)
	{
		faceType_[facei] = 1;
	}

	labelHashSet excludedPatches(excludedPatchLabels_);

	forAll(patches, patchi)
	{
		const polyPatch& pp = patches[patchi];

		if(pp.coupled() || excludedPatches.found(patchi))
		{
			forAll(pp, i)
			{
				label facei = pp.start()+i;

					faceType_[facei] = 2;
			}
		}
		else if (!isA<emptyPolyPatch>(pp))
		{
			forAll(pp, i)
			{
				label facei = pp.start()+i;

					faceType_[facei] = 1;
			}
		}
	}

	// Synchronize the faceType across processor patches 
	syncTools::syncFaceList(mesh_, faceType_, maxEqOp<label>());

	Info << "linearMotion::setRelativeMotionFaces() - completed" << endl;
	return;

	/*
	// Now we have for faceType:
	// 0 : face no in cellZone 
	// 1 : internal face or normal patch face 
	// 2 : coupled patch face or excluded patch face 
	
	// Discriminate internal moving faces 
	movingFaces_.setSize(mesh_.nInternalFaces(), false);
	label nInternal = 0;

	for ( label facei = 0; facei < mesh_.nInternalFaces(); facei++)
	{
		if (faceType_[facei] == 1)
		{
			movingFaces_[facei] = true;
			nInternal++;
		}
	}

	// Discriminate moving faces on patches 
	includedFaces_.setSize(patches.size());
	excludedFaces_.setSize(patches.size());
	forAll(patches, patchi)
	{
		includedFaces_.[patchi].setSize(patches[patchi].size(), false);
		excludedFaces_.[patchi].setSize(patches[patchi].size(), false);
	}

	label nIncluded(0);
	label nExcluded(0);

	forAll(patches, patchi)
	{
		const polyPatch& pp = patches[patchi];

		forAll(pp, patchFacei)
		{
			label facei = pp.start() + patchFacei;

			if (faceType[facei] == 1)
			{
				includedFaces_[patchi][patchFacei] = true;
				nIncluded++;
			}
			else if (faceType[facei] == 2)
			{
				excludedFaces_[patchi][patchFacei] = true;
				nExcluded++;
			}
		}
	}

	reduce(nInternal, sumOp<label>());
	reduce(nIncluded, sumOp<label>());
	reduce(nExcluded, sumOp<label>());

	Info << "linear motion: internal/included/excluded faces:  "
		<< nInternal << "/"
		<< nIncluded << "/"
		<< nExcluded 
		<< endl;

	if (debug)
	{
		faceSet internalFaces (mesh_, "internalFaces", internalFaces_);
		Pout << "Writing " << internalFaces.size()
			<< " internal faces in RelativeMotion zone to faceSet " 
			<< internalFaces.name() << endl;

		internalFaces.write();

		faceSet RelativeMotionFaces(mesh_, "includedFaces", 100);
		forAll(includedFaces_, patchi)
		{
			forAll(includedFaces_[patchi], i)
			{
				label patchFacei = includedFaces_[patchi][i];
				RelativeMotionFaces.insert(patches[patchi].start()+patchFacei);
			}
		}
		Pout << "Writing " << RelativeMotionFaces.size()
			<< " patch faces in RelativeMotion zone to faceSet "
			<< RelativeMotionFaces.name() << endl;
		RelativeMotionFaces.write();

		faceSet excludedFaces(mesh_, "excludedFaces", 100);
		forAll(excludedFaces_, patchi)
		{
			forAll(excludedFaces_[patchi], i)
			{
				label patchFacei = excludedFaces_[patchi][i];
				excludedFaces.inset(patches[patchi].start()+patchFacei);
			}
		}
		Pout << "Writing " << excludedFaces.size()
			<< " faces in RelativeMotion zone with special handling to faceSet " 
			<< excludedFaces.name() << endl;
		excludedFaces.write();
	}
	*/
}

Foam::scalar Foam::RelativeMotion::linearMotion::rampTime() const
{
	if(rampTime_ <= 0)
	{
		return 1.0;
	}

	const scalar currTime = mesh_.time().timeOutputValue();

	const word rampFunc(dict_.lookupOrDefault<word>("rampFunction", "linear"));

	scalar timeRatio = min(1.0, currTime/rampTime_);

	const scalar& PI(Foam::constant::mathematical::pi);

	if(rampFunc == "linear")
	{
		return timeRatio;
	}
	else if (rampFunc == "sin")
	{
		return sin(0.5*PI*timeRatio);
	}
	else if (rampFunc == "tanh")
		return 0.50337*(tanh(5.0*(timeRatio-0.5))+1.0);
	else
	{
		Info << "The function " << rampFunc << " is not supported. " 
			<< " Please check the entity."
			<< endl;

		return 1.0;
	}
}

void Foam::RelativeMotion::linearMotion::appendData
(
 	const scalar& currTime,
	const vector& currVelocity
)
{
	scalar dTime = currTime;

	vector sumV = currVelocity;
	vector currPosition = Zero;

	if (vList_.size() > 0)
	{
		sumV += vList_.last();
		dTime -= tList_.last();
		currPosition = pList_.last();
	}

	currPosition += 0.5*dTime*sumV;

	appendData(currTime, currVelocity, currPosition);
}

void Foam::RelativeMotion::linearMotion::appendData
(
 	const scalar& currTime,
	const vector& currVelocity,
	const vector& currPosition
)
{
	tList_.append(currTime);
	vList_.append(currVelocity);
	pList_.append(currPosition);
}

void Foam::RelativeMotion::linearMotion::updateMotionInfo()
{
	const scalar currTime(mesh_.time().timeOutputValue());

	if (currTime < clampTime_) // forces motion phase 
	{
		vector currVelocity = rampTime() * velocity_;

		if(tList_.size() == 0)
		{
			appendData(currTime, currVelocity);
		}
		else if ( tList_.last() != currTime)
		{
			appendData(currTime, currVelocity);
		}
	}
	else // free motion phase 
	{
		if ( currTime == clampTime_)
		{
			Info << "--- Release the clamp to free the object motion ---"
				<< endl;
		}

		bool updateState(false);

		if(tList_.size() == 0)
		{
			updateState = true;
		}
		else if (tList_.last() != currTime)
		{
			updateState = true;
		}

		// If time index is different, update velocity and position 
		if (updateState)
		{
			fm_.calcForcesMoment();
			vector f(fm_.forceEff());

			f -= targetForce_;

			// Constrain the motion along the velocity vector 
			Switch constrainDir = dict_.lookupOrDefault<Switch>("constrainDir", true);

			if (constrainDir)
			{
				const vector nV(velocity_/mag(velocity_));

				f = (f & nV) * nV;
			}
			else 
			{
				f.z() = 0.0;
			}

			// Sympletic 
			// Previous states
			scalar deltaT0 = lastIncrement(tList_);
			vector qDot0 = vList_.last();
			vector qDdot0 = lastIncrement(vList_)/deltaT0;
			vector q0 = pList_.last();

			// Current states 
			vector qDdot = f/mass_;
			scalar deltaT = currTime - tList_.last();

			// Print the state
			Info<< "estimated force: " << f << " N"
				<< ", accel. = " << qDdot << " m/s2"
				<< endl;

			// Predict new states
			vector qDot = qDot0 + 0.5*deltaT0*qDdot0;
			vector dq = deltaT*qDot;
			qDot += 0.5*deltaT*qDdot;

			const scalar magDq(mag(dq));

			dq = min(magDq, maxMotionStepSize_)*dq/stabilise(magDq, SMALL);

			vector q = q0 + dq;

			appendData(currTime, qDot, q);
		}
	}
}

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * //

Foam::RelativeMotion::linearMotion::linearMotion
(
 	const fvMesh& mesh,
	const dictionary& dict
)
:
	mesh_(mesh),
	dict_(dict),
	excludedPatchNames_
	(
	 	dict_.lookupOrDefault<wordRes>("nonMovingPatches", wordRes())
	),
	active_(dict_.lookupOrDefault("active", true)),
	faceType_(mesh_.nFaces(), Zero),
	fm_
	(
	 	"totalFM",
		mesh_.time(),
		dict_.subDict("fmCoeffs")
	),
	velocity_(dict_.lookup("velocity")),
	//vList_(0, Zero),
	//tList_(0, Zero),
	//pList_(0, Zero),
	rampTime_(dict_.lookupOrDefault("rampTime", 0.0)),
	clampTime_(max(rampTime_, dict_.lookupOrDefault("clampTime", 0.0))),
	motionPredFactor_(dict_.lookupOrDefault("motionPredFactor", 1.0)),
	maxMotionStepSize_(dict_.lookupOrDefault("maxMotionStepSize", 0.1)), 
	targetForce_(dict.lookupOrDefault("targetForce", vector(0.0, 0.0, 0.0)))
{
	if (!active_)
	{
		// Do nothing 
	}
	else 
	{
		const labelHashSet excludedPatchSet
		(
		 	mesh_.boundaryMesh().patchSet(excludedPatchNames_)
		);

		excludedPatchLabels_.setSize(excludedPatchSet.size());

		label i = 0;
		for (const label patchi : excludedPatchSet)
		{
			excludedPatchLabels_[i++] = patchi;
		}

		setRelativeMotionFaces();
	}

	// Estimate the mass of the object 
	scalar g(1.0);
	IOobject io
	(
	 	"g",
		mesh_.time().constant(),
		mesh_,
		IOobject::MUST_READ,
		IOobject::NO_WRITE,
		false
	);

	if (io.typeHeaderOk<uniformDimensionedVectorField>())
	{
		Info << "Reading gravitational acceleration data" << endl;
		uniformDimensionedVectorField gvec(io);
		g = mag(gvec).value();
	}

	Info << "Estimate the initial force exerted on the object" 
		<< " to estimate its mass"
		<< endl;
	fm_.calcForcesMoment();
	//mass_ = fm_.weight(g);
	mass_ = fm_.forceEff().z()/g;

	Info << "RelativeMotion info" << nl
		<< "  - mass of the object: " << mass_ << " kg"
		<< endl;

	// Update the target force 
	//targetForce_.z() = fm_.weight(1.0);
	targetForce_.z() = fm_.forceEff().z();
	Info << " - target force: " << targetForce_ << " N"
		<< endl;
}


// * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * //

bool Foam::RelativeMotion::linearMotion::active() const
{
	return active_;
}

Foam::vector Foam::RelativeMotion::linearMotion::velocity() const 
{
	return vList_.last();
}

Foam::label Foam::RelativeMotion::linearMotion::faceType
(
 	const label& facei
) const
{
	return faceType_[facei];
}

// Print motion info
void Foam::RelativeMotion::linearMotion::printMotionInfo()
{
	Info << "RelativeMotion::linear" ;

	// Update time and velocity vector 
	updateMotionInfo();

	Info << tab
		<< "Time = " << tList_.last() << " s"
		<< "velocity = " << vList_.last() << " m/s"
		<< "at " << pList_.last() << " m"
		<< endl;
}

void Foam::RelativeMotion::linearMotion::writeData(Ostream& os) const
{
	os << nl;

	os.beginBlock("linear");

	os.writeEntry("active", active_);
	os.writeEntry("velocity", velocity_);
	os.writeEntry("rampTime", rampTime_);
	os.writeEntry("clampTime", clampTime_);
	os.writeEntry("motionPredFactor", motionPredFactor_);
	os.writeEntry("maxMotionStepSize", maxMotionStepSize_);

	if(excludedPatchNames_.size())
	{
		os.writeEntry("nonMovingPatches", excludedPatchNames_);
	}
	os.endBlock();
}

bool Foam::RelativeMotion::linearMotion::read()
{
	active_ = dict_.lookupOrDefault("active", true);
	return true;
}

void Foam::RelativeMotion::linearMotion::update()
{
	setRelativeMotionFaces();
}

// ************************************************************************** //

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

// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * * //

void Foam::RelativeMotion::rotatingMotion::setRelativeMotionFaces()
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

	const labelList& own = mesh_.faceOwner();
	const labelList& nei = mesh_.faceNeighbour();

	// Cells in zone 
	boolList zoneCell(mesh_.nCells(), false);

	if (cellZoneID_ != -1)
	{
		const labelList& cellLabels = mesh_.cellZones()[cellZoneID_];
		forAll(cellLabels, i)
		{
			zoneCell[cellLabels[i]] = true;
		}
	}

	label nFaceType1 = 0;
	label nFaceType2 = 0;

	for(label facei = 0; facei<mesh_.nInternalFaces(); facei++)
	{
		if(zoneCell[own[facei]] || zoneCell[nei[facei]])
		{
			faceType_[facei] = 1;
			nFaceType1++;
		}
		else
		{
			faceType_[facei] = 0;
		}
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

				if(zoneCell[own[facei]])
				{
					faceType_[facei] = 2;
					nFaceType2++;
				}
			}
		}
		else if (!isA<emptyPolyPatch>(pp))
		{
			forAll(pp, i)
			{
				label facei = pp.start()+i;

				if(zoneCell[own[facei]])
				{
					faceType_[facei] = 1;
					nFaceType1++;
				}
			}
		}
	}

	// Synchronize the faceType across processor patches 
	syncTools::syncFaceList(mesh_, faceType_, maxEqOp<label>());

	Info << "rotatingMotion::setRelativeMotionFaces() - completed" << endl;
	return; // No need to print the number of faces in rotating zone 

	reduce(nFaceType1, sumOp<label>());
	reduce(nFaceType2, sumOp<label>());

	Info << " Number of faces : type 1 = " << nFaceType1
		<< ", type 2 = " << nFaceType2
		<< endl;
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

	Info << "rotating motion: internal/included/excluded faces:  "
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

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * //

Foam::RelativeMotion::rotatingMotion::rotatingMotion
(
 	const fvMesh& mesh,
	const dictionary& dict
)
:
	mesh_(mesh),
	dict_(dict),
	active_(dict_.lookupOrDefault("active", true)),
	cellZoneName_(dict_.lookup("cellZone")),
	cellZoneID_(-1),
	excludedPatchNames_
	(
	 	dict_.lookupOrDefault<wordRes>("nonRotatingPatches", wordRes())
	),
	faceType_(mesh_.nFaces(), Zero),
	origin_(dict_.lookup("origin")),
	axis_(dict_.lookup("axis")),
	omega_(Function1<scalar>::New("omega", dict_))
{
	if (!active_)
	{
		// Do nothing 
	}
	else 
	{

		cellZoneID_ = mesh_.cellZones().findZoneID(cellZoneName_);

		axis_ = axis_/mag(axis_);

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

		bool cellZoneFound = (cellZoneID_ != -1);
		
		reduce(cellZoneFound, orOp<bool>());

		if (!cellZoneFound)
		{
			FatalErrorInFunction
				<< "cannot find rotatingMotion cellZone " << cellZoneName_
				<< exit(FatalError);
		}

		setRelativeMotionFaces();
	}
}

// * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * //

bool Foam::RelativeMotion::rotatingMotion::active() const
{
	return active_;
}

Foam::label Foam::RelativeMotion::rotatingMotion::cellZoneID() const 
{
	return cellZoneID_;
}

Foam::vector Foam::RelativeMotion::rotatingMotion::Omega() const
{
	return omega_->value(mesh_.time().timeOutputValue())*axis_;
}

Foam::vector Foam::RelativeMotion::rotatingMotion::velocity
(
 	const point& C
) const
{
	return this->Omega() ^ (C-origin_);
}

Foam::label Foam::RelativeMotion::rotatingMotion::faceType
(
 	const label& facei
) const
{
	return faceType_[facei];
}

// Print motion info
void Foam::RelativeMotion::rotatingMotion::printMotionInfo()
{
	Info << "RelativeMotion::rotating" 
		<< tab
		<< "rotating speed " << Omega() << " rad/s"
		<< "at " << origin_
		<< endl;
}

void Foam::RelativeMotion::rotatingMotion::writeData(Ostream& os) const
{
	os << nl;

	os.beginBlock("rotating");

	os.writeEntry("active", active_);
	os.writeEntry("cellZone", cellZoneName_);
	os.writeEntry("origin", origin_);
	os.writeEntry("axis", axis_);
	omega_->writeData(os);

	if(excludedPatchNames_.size())
	{
		os.writeEntry("nonRotatingPatches", excludedPatchNames_);
	}
	os.endBlock();
}

bool Foam::RelativeMotion::rotatingMotion::read()
{
	active_ = dict_.lookupOrDefault("active", true);
	dict_.readEntry("cellZone", cellZoneName_);
	cellZoneID_ = mesh_.cellZones().findZoneID(cellZoneName_);

	return true;
}

void Foam::RelativeMotion::rotatingMotion::update()
{
	setRelativeMotionFaces();
}

// ************************************************************************** //

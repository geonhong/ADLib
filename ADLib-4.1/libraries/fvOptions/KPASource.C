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

\*---------------------------------------------------------------------------*/

#include "KPASource.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
	defineTypeNameAndDebug(KPASource, 0);
	addToRunTimeSelectionTable
	(
	 	option,
		KPASource,
		dictionary
	);

	template<class T>
	class isNotEqOp
	{
	public:
		void operator() (T& x, const T& y) const
		{
			const T unsetVal(-VGREAT*pTraits<T>::one);

			if (x != unsetVal)
			{
				// Keep x
			}
			else
			{
				x = y;
			}
		}
	};
} // End namespace Foam
} // End namespace fv


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::KPASource::updateMomentumSources()
{
	//- Extract and manipulate the wake
	extractWakeVelocity();

	//- Run KPA14
	runKPA();
}

void Foam::fv::KPASource::extractWakeVelocity()
{
	// Sample data
	const vector unsetVal(-VGREAT*pTraits<vector>::one);
	vectorField Uw(probes_.size(), unsetVal);

	const volVectorField& U = mesh_.lookupObject<volVectorField>("U");

	forAll(probes_, i)
	{
		if (elemList_[i] >= 0)
		{
			Uw[i] = U[elemList_[i]];
		}
	}

	Pstream::listCombineGather(Uw, isNotEqOp<vector>());
	Pstream::listCombineScatter(Uw);

	forAll(Uw, i)
	{
		if (Uw[i].x() < -1e10)
		{
			Uw[i] = vector::zero;
		}
	}

	Uw = propAxis_.localVector(Uw);

	// Harmonic analysis
	// This routine is carried out at the master node only
	if (Pstream::master())
	{
		if (KPA_.iV().size() > 0)
		{
			// Remove induced velocity from the extracted wake.
			// Find induced velocity at probe points by using interpolation
			const pointField probesl = propAxis_.localPosition(probes_);
			const labelListList stencilw = KPA_.makeStencil(probesl);

			// Correction factor for the induced velocity
			const scalar err05(-0.044*VM_ + 0.0292);
			const scalar err10(-0.047*VM_ + 0.0899);
			const scalar alpha(err10/(err10-err05));
			const scalar lambda(1.0 - 0.5*alpha);

			forAll(probes_, i)
			{
				if (stencilw[i][0] >= 0)
				{
					vector iV = interpolateData
					(
					 	stencilw[i],
						KPA_.iV(),
						probes_[i]
					);

					// Modify axial velocity only
					Uw[i].x() -= (lambda*iV.x());
				}
			}
		}

		// Limit wake velocity to be always positive along the
		// flow direction to enhance the stability 
		forAll(Uw, i)
		{
			Uw[i].x() = max(Uw[i].x(), 0.2*VM_);
		}

		// Write extracted wake data
		OFstream fwo("extractedWake.plt");
		fwo << "variables = y, z, Wu, Wv, Ww" << nl
			<< "zone t = \"Extracted Wake\""
			<< ", i = " << w_.nQ()
			<< ", j = " << w_.nR()
			<< endl;

		forAll(probes_, i)
		{
			const point& p = probes_[i];
			const vector& U = Uw[i];

			fwo << p.y() << tab << p.z() << tab
				<< U.x() << tab << U.y() << tab << U.z()
				<< endl;
		}

		// Carry out the harmonic analysis
		scalarField& r(w_.r());
		scalarField& q(w_.q());
		const label nR = r.size();
		const label nQ = q.size();
		const label nO = 10;

		scalar An[3][nR][nO];
		scalar Bn[3][nR][nO];

		const scalar dq(q[1]-q[0]);
		const scalar coef(2.0/scalar(nQ)/VM_);

		for (direction dir=0; dir<3; dir++)
		{
			scalarField Uwi(Uw.component(dir));

			for (int iR=0; iR<nR; iR++)
			{
				for (int iO=0; iO<nO; iO++)
				{
					An[dir][iR][iO] = 0.0;
					Bn[dir][iR][iO] = 0.0;

					for (int iQ=0; iQ<nQ; iQ++)
					{
						const label ind = iR*nQ + iQ;
						const scalar arg(dq*scalar(iQ)*scalar(iO));

						An[dir][iR][iO] += coef*Uwi[ind]*Foam::cos(arg);
						Bn[dir][iR][iO] += coef*Uwi[ind]*Foam::sin(arg);
					}
				}

				An[dir][iR][0] *= 0.5;
			}
		}

		// Write the harmonic analysis result to a file
		Info<< "Writing the results of harmonic analysis on wake "
			<< "to file wakeHarmonicAnalysis.wak..." << flush;

		OFstream fout("wakeHarmonicAnalysis.wak");
		List<word> key(3);
		key[0] = "AXIALWAKE";
		key[1] = "RADIALWAKE";
		key[2] = "TANGENTIALWAKE";

		fout<< nR << endl;
		forLoop(iD, 3)
		{
			fout<< key[iD] << endl;

			const scalar R(0.5*dia_);

			writeCoeffs(fout, An);
			writeCoeffs(fout, Bn);
		}

		Info<< "done" << endl;
	}
}

// Run KPA14 code to estimate the momentum sources
// 1. Run KPA14 code at master node
// 2. Broadcast the results to slave nodes
// 3. Interpolate the results to OpenFOAM fields
void Foam::fv::KPASource::runKPA()
{
	Info<< "Run KPA14 to estimate momentum sources"
		<< endl;

	scalarField& Rc(KPA_.R());
	scalarField& Qc(KPA_.Q());
	scalarField& dR(KPA_.dR());
	scalarField& F(KPA_.F());
	scalarField& iV(KPA_.iV());

	Rc.resize(RMAX_KPA);
	dR.resize(RMAX_KPA);
	Qc.resize(QMAX_KPA);

	if (Pstream::master())
	{
		double Qs[RMAX_KPA][QMAX_KPA];
		double mF[3][RMAX_KPA][QMAX_KPA];
		double miV[3] [RMAX_KPA][QMAX_KPA];

		int nBlades(0);		// Read from *.geo file
		int nRadial(0);		// MMR (=20)
		int nTheta(0);		// NSTEPP (=60)
		int nR(0);			// Read from *.wak file
		int nQ(0);			// (=NSTEPP)
		int iSect(1);		// IMT

		// Run KPA14
		bemuf_
		(
		 	propGeo_.c_str(), &rps, &VM_, &SR_, &rho_, &dia_,
			&nBlades, &nRadial, &nTheta, &iSect,
			&nR, &nQ,
			&Rc[0], &dR[0], &Qc[0], Qs,
			mF [0], mF [1], mF [2],
			miV[0], miV[1], miV[2]
		);

		// Shrink size of scalarFields
		Rc.resize(nRadial);
		dR.resize(nRadial);
		Qc.resize(nRadial);

		// Assign results to force and induced velocity fields
		const label fldSize(nRadial*nTheta);
		F.resize(fldSize);
		iV.resize(fldSize);
		
		forLoop(iR, nRadial)
		{
			forLoop(iQ, nTheta)
			{
				const label ind(iR*nTheta + iQ);
				vector& Fi(F[ind]);
				vector& iVi(iV[ind]);

				// Assign force and induced velocity data
				Fi = vector(mF[0][iR][iQ], mF[1][iR][iQ], mF[2][iR][iQ]);
				iVi = vector(miV[0][iR][iQ], miV[1][iR][iQ], miV[2][iR][iQ]);

				// Transform coordinate from r-q to y-z plane
				const scalar arg(Qc[iQ]*pi/180.0);
				const scalar Fy(Fi.y()*Foam::cos(arg) - Fi*z()*Foam::sin(arg));
				const scalar Fz(Fi.y()*Foam::sin(arg) + Fi.z()*Foam::cos(arg));

				Fi.y() = Fy;
				Fi.z() = Fz;
			}
		}

		F = propAxis_.globalVector(F);

		// Estimate force density
		// 1. Estimate area of the disk
		// dimensionalize the resultant data
		scalarField A(nRadial);
		const scalar dq = (Qc[1] - Qc[0])*pi/180.0;
		forAll(Rc, iR)
		{
			Rc[iR] *= 0.5*dia_;
			dR[iR] *= 0.5*dia_;
			A[iR] = dR[iR]*Rc[iR]*dq;
		}

		// 2. Estimate force density
		const scalar dimCoef
		(
		 	0.5*rho_*sqr(VM_*dia_)*0.25*scalar(nBlades)/scalar(nTheta)
		);

		forLoop(iR, nRadial)
		{
			const scalar coef = dimCoef/A[iR];

			forLoop(iQ, nTheta)
			{
				const label ind(iR*nTheta + iQ);
				F[ind] *= coef;
			}
		}
	}

	// Broadcast the variables to slave processors
	KPA_.scatter();

	// Interpolate the estimated forces to obtain momentum sources
	// 1. construct stencil if required
	if (!stencil_)
	{
		vectorField C(cells_.size());
		forAll(C, i)
		{
			C[i] = mesh_.C()[cels_[i]];
		}

		C = propAxis_.localPosition(C);

		stencil_ = new labelListList(KPA_.makeStencil(C, dia_));
	}

	// 2. interpolate forces from KPA results to rotate cellZone
	scalar Fx(0.0);
	label nTrgCells(0);
	scalar trgV(0.0);
	forAll(cells_, i)
	{
		const label& ic(cells_[i]);
		labelList& stencil = stencil_->operator[](i);

		const point& C(mesh_.C()[ic]);
		const point ploc(propAxis_.localPosition(C));
		scalar rloc(mag(vector(0.0, ploc.y(), ploc.z())));

		if (rloc<0.5*dia_)
		{
			nTrgCells++;
			trgV += mesh_.V()[ic];
			source_[i] = interpolateData(stencil, i)/0.05;
		}
		else
		{
			source_[i] = vector::zero;
		}

		Fx -= source_[i].x()*mesh_.V()[ic];
	}

	reduce(Fx, sumOp<scalar>());

	Info<< "Integrated momentum source in x-direction: "
		<< Fx << " N" << nl
		<< "   + number of cells: " << nTrgCells << nl
		<< "   + target cell vol: " << trgV << nl
		<< endl;
}

//- Interpolate the resultant forces to momentum source fields
Foam::vector Foam::fv::KPASource::interpolateData
(
 	const labelList& stencil,
	const label& itrg
)
{
	// Get coordinate data of target cell in local coordinate system,
	// where the local coordinate system represents propeller axis.
	point Ctrg
	(
	 	propAxis_.localPosition(mesh_.C()[cells_[itrg]]);
	);

	return interpolateData(stencil, KPA_.F(), Ctrg);
}

// Interpolate data by means of inverse distance weighting method
Foam::vector Foam::fv::KPASource::interpolateData
(
 	const labelList& stencil,
	const vectorField& fld,
	const point& C
)
{
	vector num(vector::zero);
	scalar den(0.0);

	forAll(stencil, i)
	{
		const label isrc(stencil[i]);

		point Csrc(KPA_.yz(i));

		scalar d(mag(Csrc - vector(0.0, C.y(), C.z())));

		if (d<SMALL)
		{
			return fld[isrc];
		}

		scalar w(1.0/sqr(d));

		num += w*fld(isrc);
		den += w;
	}

	return num/den;
}

//- Get probe points
void Foam::fv::KPASource::getProbePoints()
{
	const label nR(w_.nR());
	const label nQ(w_.nQ());

	probes_.clear();
	probes_.setSize(nR*nQ);

	const scalar rmin(1.005*this->rh());
	const scalar rmax(1.05*0.5*dia_);
	const scalar dr((rmax-rmin)/scalar(nR-1));

	scalarField& r(w_.r());
	r[0] = rmin;
	for (label iR=1; iR<nR; iR++)
	{
		r[iR] = r[iR-1] + dr;
	}

	scalarField& q(w_.q());
	for (label iQ=0; iQ<nQ; iQ++)
	{
		q[iQ] = 2.0*pi*scalar(iQ)/180.0;
	}

	for (label iR=0; iR<nR; iR++)
	{
		for (label iQ=0; iQ<nQ, iQ++)
		{
			const label ind = iR*nQ + iQ;
			probes_[ind] =
				point
				(
					-0.001,
					r[iR]*Foam::cos(q[iQ]),
					r[iR]*Foam::sin(q[iQ])
				);
		}
	}

	// Convert local to global
	tmp<pointField> tpf = propAxis_.globalPosition(probes_);
	probes_ = tpf();

	// Find elements for probing data
	elemList_.clear();
	elemList_.setSize(probes_.size(), -1);

	Info<< "Searching element cells for probing data" << endl;
	time_t startAt = mesh_.time().elapsedClockTime();

	Info<< "> Generating samples for searching element cells...";
	labelList samples;
	forAll(mesh_.C(), i)
	{
		const point& C(propAxis_.localPosition(mesh_.C()[i]));

		if
		(
			(C.x()>-0.1 && C.x()<0.1)
		 && (C.y()>-dia_ && C.y()<dia_)
		 && (C.z()>-dia_ && C.z()<dia_)
		)
		{
			samples.append(i);
		}
	}

	Info<< samples.size() << " cells are selected as samples"
		<< endl;

	forAll(probes_, i)
	{
		const vector& Xb = probes_[i];
		{
			const label cellI = findCell(samples, Xb);
			elemList_[i] = cellI;
		}

		// Print progress bar
		label totalProbes = probes_.size();
		label index = i;

		if ((index+1)%20 == 0 || index==totalProbes-1)
		{
			if (index>21)
			{
				Info<< "\r";
			}Info<< "[";
			const label progress =
				label(scalar((index+1))/scalar(totalProbes)*100.0);

			for (label count=0; count<100; count++)
			{
				if (count == progress-1)
				{
					Info<< ">";
				}
				else if (count < progress)
				{
					Info<< "=";
				}
				else
				{
					Info<< " ";
				}
			}
			Info<< "] " << progress << "%"
				<< " (" << index+1 << "/" << totalProbes << ")"
				<< flush;
		}
	}

	Info<< endl;
	Info<< "It took "
		<< mesh_.time().elapsedClockTime() - startAt
		<< " seconds for searching probe points" << nl
		<< "done." << nl
		<< endl;
}

// Find cell
Foam::label Foam::fv::KPASource::findCell
(
	const labelList& samples,
	const point& X
)
{
	if (samples.size() == 0)
	{
		return -1;
	}

	const vectorField& C(mesh_.C());

	// Find nearest cell
	label itrg(0);
	scalar minProx(magSqr(C[0] - X));

	forAll(samples, i)
	{
		const label currentI(samples[i]);
		scalar prox = magSqr(C[currentI] - X);

		if (prox < minProx)
		{
			itrg = currentI;
			minProx = prox;
		}
	}

	// Check if the point is in the nearest cell
	if (pointInCell(X, itrg))
	{
		return itrg;
	}
	else
	{
		bool cellFound(false);
		label n(0);

		while((!cellFound) && (n<samples.size()))
		{
			if (pointInCell(X, samples[n]))
			{
				cellFound = true;
				itrg = samples[n];
			}
			else
			{
				n++;
			}
		}

		if (cellFound)
		{
			return itrg;
		}
		else
		{
			return -1;
		}
	}
}

// Identical to primitiveMesh::pointInCell
// This is re-written to avoid MPI_Recv error
bool Foam::fv::KPASource::pointInCell(const point& p, label celli) const
{
	const labelList& f = mesh_.cells()[celli];
	const labelList& owner = mesh_.faceOwner();
	const vectorField& Cf = mesh_.faceCentres();
	const vectorField& Sf = mesh_.faceAreas();

	bool inCell = true;

	forAll(f, facei)
	{
		label nFace = f[facei];
		vector proj = p - Cf[nFace];
		vector normal = Sf[nFace];
		if (owner[nFace] != celli)
		{
			normal = -normal;
		}
		inCell = inCell && ((normal & proj) <= 0);
	}

	return inCell;
}

// Get propeller geometry data by reading .geo file
bool Foam::fv::KPASource::getPropGeoData()
{
	propGeo_ = propName_ + ".geo";
	propMt_ = propName_ + ".mt";

	IFstream fin(propGeo_.c_str());

	if (fin.opened())
	{
		Info<< "Read propeller geometry data from " << propGeo_ << endl;
		string line;
		fin.getLine(line);
		fin.getLine(line);
		fin.getLine(line);	// Prop name
		fin.getLine(line);	// Designer
		fin.getLine(line);	// geometry data

		// Get diameter info
		std::stringstream strm(line);
		strm>> dia_;

		// Get hub diameter
		fin.getLine(line);
		fin >> rh_;

		const bool fullScale = 
			coeffs_.lookupOrDefault<bool>("fullScale", false);

		if (fullScale)
		{
			Info<< " >> Simulation will be carried out "
				   "by using full scale model"
				   "(diameter of propeller = " << dia_ << ")"
				<< endl;

			return true;
		}

		if (dia_>1.0)
		{
			WarningIn("KPASource::getPropGeoData()")
				<< "The diameter of given propeller might be given "
				   "in ship scale rather than model scale." << nl
				<< "Given diameter: " << dia_
				<< endl;

			if (dia_ > 1000)
			{
				Info<< " >> Convert dimension from mm to m" << endl;
				dia_ *= 0.001;
			}

			dia_ /= SR_;
		}

		return true;
	}
	else
	{
		WarningIn("KPASource::getPropGeoData()")
			<< "Propeller geometry file " << propGeo_
			<< " is not found."
			   " Please check if appropriate file is provided or "
			   " the file name is properly described in sourceProperties."
			<< endl;

		return false;
	}
}
// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::fv::KPASource::KPASource
(
	const word& name,
	const word& modelType,
	const dictionary& dict,
	const fvMesh& mesh
)
:
	basicSource(name, modelType, dict, mesh),
	nR_(coeffs_.lookupOrDefault<label>("nR", 14)),
	nQ_(coeffs_.lookupOrDefault<label>("nQ", 179)),
	w_(nR_, nQ_),
	source_(cells_.size()),
	propAxis_(coeffs_.subDict("propAxis")),
	stencil_(NULL),
	nUpdate_(0),
	pi(constant::mathematical::pi)
{
	coeffs_.lookup("fieldNames") >> fieldNames_;
	applied_.setSize(fieldNames_.size(), false);

	Info<< "  - Creating KPA zone: " << this->name() << endl;

	read(dict);

	// Get position info of rotate region and localize it
	rotCfld_.clear();
	rotCfld_.setSize(cells_.size());
	forAll(cells_, i)
	{
		rotCfld_[i] = mesh_.C()[cells_[i]];
	}

	tmp<pointField> tCfld = propAxis_.localPosition(rotCfld_);
	rotCfld_ = tCfld();

	Info<< "Center of rotate region: "
		<< Foam::average(rotCfld_) << nl
		<< "Propeller origin: "
		<< propAxis_.origin() << nl
		<< endl;
	
	// Initialize source terms
	// Initial guess for the thrust (KT ~ 0.18)
	const scalar Thrust = 0.18*rho_*sqr(rps_)*sqr(sqr(dia_));
	forAll(source_, i)
	{
		source_[i] = vector(-Thrust, 0.0, 0.0)*mesh_.V()[cells_[i]];
	}

	// Get probe points
	getProbePoints();
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::KPASource>
Foam::KPASource::New()
{
    return autoPtr<KPASource>(new KPASource);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::fv::KPASource::addSup
(
	fvMatrix<vector>& eqn,
	const label fieldI
)
{
	const label& timeIndex = mesh_.time().timeIndex();
	const label runTimeIndex = timeIndex - mesh_.time().startTimeIndex();
	const label nMaxUpdates =
		coeffs_.lookupOrDefault<label>("nMaxUpdates", 100);
	
	const bool initialUpdate =
		coeffs_.lookupOrDefault<bool>("initialUpdate", true);
	
	if
	(
		(timeIndex == 1 && initialUpdate)
	 ||
	 	(
			nUpdate_ < nMaxUpdates
		 &&
		 	(
				timeIndex % interval_ == 0
			 ||
			 	(timeIndex > 1 && runTimeIndex == 1)
			)
		)
	)
	{
		updateMomentumSources();
		nUpdate_++;
		Info<< "Source terms are updated " << nUpdate_
			<< " times out of " << nMaxUpdates
			<< endl;
	}

	vectorField& Usrc = eqn.source();
	// Note that the effect of density is reflected already when
	// the resultant forces of KPA are dimensionalized.
	// Thus the density is not taken into account here.

	forAll(cells_, i)
	{
		const label ic = cells_[i];
		Usrc[ic] += source_[i]*mesh_.V()[ic];
	}
}

void Foam::fv::KPASource::addSup
(
 	const volScalarField& rho,
	fvMatrix<vector>& eqn,
	const label fieldI
)
{
	this->addSup(eqn, fieldI);
}

void Foam::fv::KPASource::writeData(Ostream& os) const
{
	os  << indent << name_ << endl;
	dict_.write(os);
}

bool Foam::fv::KPASource::read(const dictionary& dict)
{
	if (basicSource::read(dict))
	{
		vector velocity;
		coeffs_.readIfPresent("rps", rps_);
		coeffs_.readIfPresent("velocity", velocity);
		coeffs_.readIfPresent("density", rho_);
		coeffs_.readIfPresent("scaleRatio", SR_);
		coeffs_.readIfPresent("propName", propName_);
		coeffs_.readIfPresent("interval", interval_);
		VM_ = velocity.x();

		getPropGeoData();

		Info<< "propName : " << propName_ << endl
			<< "diameter : " << dia_ << endl
			<< "hub radi.: " << this->rh() << endl
			<< "rps      : " << rps_ << endl
			<< "density  : " << rho_ << endl
			<< "velocity : " << VM_ << endl
			<< "scaleRat.: " << SR_ << endl
			<< "interval : " << interval_ << endl
			<< endl;

		return true;
	}
	else
	{
		return false;
	}
}

// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


// ************************************************************************* //

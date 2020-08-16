/*
	Skip the first part of the KPASource.C
*/

// Interpolate the resultant forces to momentum source fields
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
		propAxis_.localPosition(mesh_.C()[cells_[itrg]])
	);

	return interpolationData(stencil, KPA_.F(), Ctrg);
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

		num += w*fld[isrc];
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


// * * * Constructors * * * //

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

// * * * Member Functions * * * //

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

// ***** //
}

/*
	FOAM | OpenFOAM
 */

#include "mathematicalConstants.H"
#include "Tuple2.H"

#define QMAX_KPA 200
#define RMAX_KPA 21

#define forLoop(ind, N) \
	for(int ind = 0; ind < N; ind++)

#define  writeCoeffs(fout, Cn) \
	/* Write radial position info */ \
	forLoop(iR, nR) \
	{ \
		fout<< r[iR]/R << tab; \
	} \
	fout<< endl; \
	\
	/* Write coefficient */ \
	forLoop(iO, nO) \
	{ \
		forLoop(iR, nR) \
		{ \
			fout<< Cn[iD][iR][iO] << tab; \
		} \
		fout<< endl; \
	}

using namespace Foam::constant:: mathematical;

//- Declaration of the function for calling KPA14 code
extern "C"
{
	void bemuf_
	(
		const char* propgeo,
		const double* rps,
		const double* vm,
		const double* SR,
		const double* rho,
		double* dia,
		int* nBlades,
		int* nRadial,
		int* nTheta,
		int* iSect,
		int* nR,
		int* nQ,
		double* Rc,
		double* dR, 
		double* Qline,
		double* (*Qskew)[QMAX_KPA],
		double* (*Fx)[QMAX_KPA],
		double* (*Fy)[QMAX_KPA],
		double* (*Fz)[QMAX_KPA],
		double* (*iVx)[QMAX_KPA],
		double* (*iVy)[QMAX_KPA],
		double* (*iVz)[QMAX_KPA]
	);
}

// * * * //

namespace Foam
{
namespace fv
{

class wake
{
	scalarField r_;
	scalarField q_;

	//- Disallow null constructor
	wake() {}

public:
	
	// Constructor
		
		wake(const label& nR, const label& nQ)
		:
			r_(nR),
			q_(nQ)
		{}
	
	// Access
		
		label nR() const
		{
			return r_.size();
		}

		label nQ() const
		{
			return q_.size();
		}

		scalarField& r()
		{
			return r_;
		}

		scalarField& q()
		{
			return q_;
		}
};

//- Helper class to store and manipulate the results of
//  KPA14 code
class KPAData
{
	scalarField R_;
	scalarField Q_;
	scalarField dR_;
	vectorField F_;
	vectorField iV_;

public:

	KPAData()
	{}

	KPAData(const label& n)
	:
		R_(n),
		Q_(n),
		dR_(n),
		F_(n),
		iV_(n)
	{}

	typedef Tuple2<label, label> labels;

	//- Public memeber functions
		template<class Type>
		void resize(const label& n, Field<Type>& fld)
		{
			fld.clear();
			fld.resize(n);
		}

		label index(const label& iR, const label& iQ) const
		{
			return iR*Q_.size() + iQ;
		}

		const labels index(const label& ind)
		{
			const label iQ = ind % Q_.size();
			const label iR = (ind-iQ)/Q_.size();

			return labels(iR, iQ);
		}

		const point yz(const label& iR, const label& iQ)
		{
			const scalar arg(Q_[iQ]*pi/180.0);
			const scalar y(R_[iR]*Foam::cos(arg));
			const scalar z(R_[iR]*Foam::sin(arg));

			return point(0.0, y, z);
		}

		const point yz(const label& ind)
		{
			const labels iRQ(index(ind));

			return yz(iRQ.first(), iRQ.second());
		}

		void scatter()
		{
			Pstream::scatter(R_);
			Pstream::scatter(Q_);
			Pstream::scatter(dR_);
			Pstream::scatter(F_);
			Pstream::scatter(iV_);
		}

		label nR() const
		{
			return R_.size();
		}

		label nQ() const
		{
			return Q_.size();
		}

		labelListList makeStencil(const vectorField& C)
		{
			return makeStencil(C, 1e10);
		}


		// Make a stencil for interpolation from KPA grid to
		// target grid, represented by vectorField C
		labelListList makeStencil
		(
			const vectorField& C,
			const scalar rMax
		)
		{
			labelListList stencil(C.size());
			const label nR(this->nR());
			const label nQ(this->nQ());

			forAll(C, i)
			{
				// Get local coordinate info of a cell in rotate region
				const point ploc(C[i]);

				scalar rloc(mag(vector(0.0, ploc.y(), ploc.z())));
				scalar qloc(atan2(ploc.z(), ploc.y())*180.0/pi);

				if (qloc < 0.0)
				{
					qloc += 360.0;
				}

				// Search and construct stencils
				if (rloc > 0.5*rMax)
				{
					// out of bound
					stencil[i] = labelList(1, -1);
				}
				else
				{
					// The case of inbound
					label iRtrg(0);

					for (int iR=0; iR<nR-1; iR++)
					{
						scalar chkval((rloc-R_[iR])*(rloc-R_[iR+1]));
						if (chkval<=0.0)
						{
							iRtrg = iR;
							break;
						}
					}

					label iQtrg(0);
					for (int iQ=0; iQ<nQ; iQ++)
					{
						if (iQ == nQ-1)
						{
							iQtrg = iQ;
							break;
						}

						scalar chkval((qloc-Q_[iQ])*(qloc-Q_[iQ+1]));
						if (chkval<=0.0)
						{
							iQtrg = iQ;
							break;
						}
					}

					// Get stencil info
					labelList& stencili = stencil[i];

					label iQnext(iQtrg + 1);
					if (iQnext == nQ)
					{
						iQnext = 0;
					}

					stencili.append(iRtrg*nQ + iQtrg);
					stencili.append(iRtrg*nQ + iQnext);

					if (rloc>=min(R_) && rloc<=max(R_))
					{
						stencili.append((iRtrg+1)*nQ + iQtrg);
						stencili.append((iRtrg+1)*nQ + iQnext);
					}
				}
			}

			return stencil;
		}

	//- Access
		
		scalarField& R()
		{
			return R_;
		}

		scalarField& Q()
		{
			return Q_;
		}

		scalarField& dR()
		{
			return dR_;
		}

		vectorField& F()
		{
			return F_;
		}

		vectorField& iV()
		{
			return iV_;
		}
};

} // End namespace fv
} // End namespace Foam

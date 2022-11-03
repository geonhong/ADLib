#ifndef limitU_H
#define linitU_H

namespace FOAM
{

void limitU(volVectorField& U, const pimpleControl& pimple)
{
	if (pimple.dict().found("limitMagU"))
	{
		dimensionedScalar limitMagU
		(
		 	"limitMagU",
			U.dimensions(),
			readScalar(pimple.dict().lookup("limitMagU"))
		);

		const fvMesh& mesh(U.mesh());

		volScalarField magU(mag(U));
		scalar maxMagU = max(magU).value();

		dimensionedScalar small("small", U.dimensions(), SMALL);


		if (maxMagU > limitMagU.value())
		{
			U.ref() *=
					neg(magU.internalField() - limitMagU)
				+	pos(magU.internalField() - limitMagU) *
					limitMagU/(magU.internalField() + small);
			U.correctBoundaryConditions();

			Info << "mag(U): max: " << maxMagU
				<< " avg: " << magU.weightedAverage(mesh.V()).value()
				<< " ...limiting" << endl;
		}
		else
		{
	//		Info<< endl;
		}
	}
}
}

#endif 

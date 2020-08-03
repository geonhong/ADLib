#ifndef LIMITU_H
#define LIMITU_H

namespace Foam
{

void limitU(volVectorField& U, const solutionControl& control)
{
	const volScalarField magU(mag(U));
	const scalar maxMagU(max(magU).value());

	const scalar limitMagU(readScalar(control.dict().lookup("limitMagU")));

	Info<< "mag(U): max: " << maxMagU
    	<< " avg: " << magU.weightedAverage(U.mesh().V()).value();

	if (limitMagU - maxMagU < 0.0)
	{
		volVectorField::Internal::FieldType& Uref(U.primitiveFieldRef());
		const volScalarField::Internal::FieldType& magUref(magU.primitiveField());
		volScalarField::Internal::FieldType  magUdiff(magUref - limitMagU);

		Uref *= neg(magUdiff) + pos(magUdiff)*limitMagU/(magUref + SMALL);

		U.correctBoundaryConditions();

		scalar maxMagU(max(mag(U)).value());
		reduce(maxMagU, sumOp<scalar>());

		Info<< "...clipping : max(U) = " << maxMagU;
	}

	Info<< endl;
}

} // End namespace Foam

#endif

{
    volScalarField& rDeltaT = trDeltaT.ref();

    scalar rDeltaTSmoothingCoeff
    (
        runTime.controlDict().getOrDefault<scalar>
        (
            "rDeltaTSmoothingCoeff",
            0.02
        )
    );

    // Set the reciprocal time-step from the local Courant number
	maxDeltaT = max(maxDeltaT, maxDeltaT0);
	Info<< "maxDeltaT = " << maxDeltaT << endl;
    rDeltaT.ref() = max
    (
        1/dimensionedScalar("maxDeltaT", dimTime, maxDeltaT),
        fvc::surfaceSum(amaxSf)()()
       /((2*maxCo)*mesh.V())
    );

    // Update tho boundary values of the reciprocal time-step
    rDeltaT.correctBoundaryConditions();

    fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);

	scalar minTscl(gMin(1/rDeltaT.primitiveField()));
	scalar maxTscl(gMax(1/rDeltaT.primitiveField()));

    Info<< "Flow time scale min/max = "
        << minTscl << ", " << maxTscl << endl;

	if (minTscl == maxTscl)
	{
		maxDeltaT *= 1.1;

		// Set the reciprocal time-step from the local Courant number
		Info<< "maxDeltaT = " << maxDeltaT << endl;
		rDeltaT.ref() = max
		(
			1/dimensionedScalar("maxDeltaT", dimTime, maxDeltaT),
			fvc::surfaceSum(amaxSf)()()
		   /((2*maxCo)*mesh.V())
		);

		// Update tho boundary values of the reciprocal time-step
		rDeltaT.correctBoundaryConditions();

		fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);

		scalar minTscl(gMin(1/rDeltaT.primitiveField()));
		scalar maxTscl(gMax(1/rDeltaT.primitiveField()));

		Info<< "Flow time scale min/max = "
			<< minTscl << ", " << maxTscl << endl;
	}

	maxDeltaT0 = maxDeltaT;
}

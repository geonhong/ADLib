    // Momentum predictor

    MRF.correctBoundaryVelocity(U);

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)
      + MRF.DDt(U)
      + turbulence->divDevReff(U)
     ==
        fvOptions(U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

	/*
	// Check UEqn
	const volScalarField& A(UEqn.A()());

	forAll(A, ci) 
	{
		if (mag(A[ci]) <= 1e-10) 
		{
			Pout<< "Cell " << ci
				<< " C = " << mesh.C()[ci]
				<< ", A = " << A[ci]
				<< ", V = " << mesh.V()[ci]
				<< endl;
		}
	}
	*/

    fvOptions.constrain(UEqn);

    if (simple.momentumPredictor())
    {
        solve(UEqn == -fvc::grad(p));

        fvOptions.correct(U);
    }

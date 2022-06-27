    // MRF.correctBoundaryVelocity(U);
	RM.correctBoundaryVelocity(U);

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) 
	  // + fvm::div(rhoPhi, U)
      + RM.DDt(rho, U)
      + turbulence->divDevRhoReff(rho, U)
     ==
        fvOptions(rho, U)
    );

	//-----------------------------------------------------------------------------
	// deferred correction
	Switch deferredCorr
	(
		pimple.dict().lookupOrDefault<Switch>("deferredCorr", "off")
	);

	if (deferredCorr)
	{
		UEqn +=
		(
			fv::gaussConvectionScheme<vector>
			(
				mesh,
				rhoPhi,
				upwind<vector>(mesh, rhoPhi)
			).fvmDiv(rhoPhi, U)
		 ==
			fv::gaussConvectionScheme<vector>
			(
				mesh,
				rhoPhi,
				upwind<vector>(mesh, rhoPhi)
			).fvcDiv(rhoPhi, U)
		  - fvc::div(rhoPhi,U)
		);
	}
	else
	{
		UEqn += fvm::div(rhoPhi, U);
	}
	//-----------------------------------------------------------------------------

	if (pimple.finalIter())
	{
		UEqn.relax(1.0);
	}
	else
	{
		UEqn.relax();
	}

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    mixture.surfaceTensionForce()
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }

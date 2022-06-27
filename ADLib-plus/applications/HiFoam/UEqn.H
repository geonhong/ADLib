RM.correctBoundaryVelocity(U);

//if (mesh.moving())
//{
//		U.correctBoundaryConditions();
//}
//
//RM.checkBoundaryVelocity("blade", "RM.correctBoundaryVelocity(U)");

fvVectorMatrix UEqn
(
 	fvm::ddt(rho, U)
	/*
	+ fvm::div(rhoPhi, U)
	+ fv::gaussConvectionScheme<vector>
		(
		 	mesh,
			rhoPhi,
			upwind<vector>(mesh, rhoPhi)
		).fvmDiv(rhoPhi, U)
		+ MRF.DDt(rho, U)
	*/
		+ RM.DDt(rho, U)
		+ turbulence->divDevRhoReff(rho, U)
		==
		/*
			fv::gaussConvectionScheme<vector>
			(
			 	mesh,
				rhoPhi,
				upwind<vector>(mesh, rhoPhi)
			).fvcDiv(rhoPhi, U)
		- fvc::div(rhoPhi, U)
		*/
		fvOptions(rho, U)
);

// RM.checkBoundaryVelocity("blade", "After discretization");

Switch deferredCorr
(
 	pimple.dict().lookupOrDefault<Switch>("defferredCorr", "off")
);

if (deferredCorr)
{
	Info << "Applying deferred correction for the momentum equation" << endl;
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
		- fvc::div(rhoPhi, U)
	);
}
else
{
	UEqn += fvm::div(rhoPhi, U);
}

if (pimple.finalIter())
{
	UEqn.relax(1.0);
}
else
{
	UEqn.relax(mesh.equationRelaxationFactor(U.name()));
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
					-ghf*fvc::snGrad(rho)
					-fvc::snGrad(p_rgh)
					*mesh.magSf()
				)
			)
		);

	fvOptions.correct(U);
}



{
	if (correctPhi)
	{
		rAU.ref() = 1.0/UEqn.A();
	}
	else
	{
		rAU = 1.0/UEqn.A();
	}

	surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU()));
	volVectorField HbyA(constrainHbyA(rAU()*UEqn.H(), U, p_rgh));

	// Limiting the field value 
	FOAM::limitU(HbyA, pimple);

	surfaceScalarField phiHbyA
	(
	 	"phiHbyA",
		fvc::flux(HbyA)
		// + MRF.zeroFilter(fvc::interpolate(rho*rAU())*fvc::ddtCorr(U, phi, Uf))
		+ RM.zeroFilter(fvc::interpolate(rho*rAU())*fvc::ddtCorr(U, phi, Uf))
	);

	//MRF.makeRelative(phiHbyA);
	RM.makeRelative(phiHbyA);

	//RM.checkBoundaryVelocity("blade", "final");
	
	if (p_rgh.needReference())
	{
		fvc::makeRelative(phiHbyA, U);
		adjustPhi(phiHbyA, U, p_rgh);
		fvc::makeAbsolute(phiHbyA, U);
	}

	surfaceScalarField phig
	(
	 	(
		 	mixture.surfaceTensionForce()
			-ghf*fvc::snGrad(rho)
		)*rAUf*mesh.magSf()
	);

	phiHbyA += phig;

	// Update the pressure BCs to ensure flux consistency 
	// constrainPressure(p_rgh, U, phiHbyA, rAUf, MRF);
	constrainPressure(p_rgh, U, phiHbyA, rAUf, RM);

	// Select solver algorithm
	// pEqnAlgorithm pimple; 		// original pimple algorithm
	// pEqnAlgorithm transimple; 	//transient simple algorithm
	word pEqnAlgorithm = 
		pimple.dict().lookupOrDefault<word>("pEqnAlgorithm", "pimple");

	while (pimple.correctNonOrthogonal())
	{
		fvScalarMatrix p_rghEqn
		(
		 	fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA)
		);

		p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

		if (pEqnAlgorithm == "transimple")
		{
			p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalIter())));
		}
		else
		{
			// Solve the pressure equation using the original PIMPLE algorithm
			// as default 
			p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));
		}

		if (pimple.finalNonOrthogonalIter())
		{
			phi = phiHbyA - p_rghEqn.flux();

			// Original pressure relaxation 
			// p_rgh.relax();

			U = HbyA + rAU()*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf);
			U.correctBoundaryConditions();
			fvOptions.correct(U);

				//FOAM::limitU(U. pimple);
				
				// Relax the pressure explicitly AFTER velocity correction 
			p_rgh.relax();
		}
	}

	#include "continuityErrs.H"

	// Correct Uf if the mesh is moving 
	fvc::correctUf(Uf, U, phi);

	// Make the fluxes relative to the mesh motion 
	fvc::makeRelative(phi, U);

	p == p_rgh + rho*gh;

	if (p_rgh.needReference())
	{
		p += dimensionedScalar
		(
		 	"p",
			p.dimensions(),
			pRefValue - getRefCellValue(p, pRefCell)
		);
		p_rgh = p - rho*gh;
	}

	if (!correctPhi)
	{
		rAU.clear();
	}

	// RM.checkBoundaryVelocity("blade", "fianl");
}




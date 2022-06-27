{
	Info << "Limiting VoF field and fluxes";

	scalar excessiveVol(sum(fvc::volumeIntegrate(max(0.0, alpha1-1))));
	reduce (excessiveVol, sumOp<scalar>());

	Info <<", excessive volume: " << excessiveVol << " m3" << endl;

	// Limit volume fraction
	alpha1 = min(max(0.0, alpha1), 1.001);
	alpha2 = 1-alpha1;

	// Limit face flux
	alphaPhi10 = sign(alphaPhi10)*
		min
		(
		 	mag(alphaPhi10),
			mag(phiCN)
		);
	/*
	 *  Do not limit the VoF flux
	const labelList& own(mesh.owner());
	const labelList& nei(mesh.neighbour());

	forAll(rhoPhi, fi)
	{
		const label cio(own[fi]);
		const label cin(nei[fi]);

		scalar alpha1max = alpha1[cio];

		if (cin != -1)
		{
			alpha1max = max(alpha1[cio], alpha1[cin]);
		}

		scalar rhof = (alpha1max*(rho1-rho2) + rho2).value();

		rhoPhi[fi] = sign(rhoPhi[fi])
			* min
				(
					mag(rhoPhi[fi]).
					rhof*mag(phi[fi])
				);
	}
	*/
}


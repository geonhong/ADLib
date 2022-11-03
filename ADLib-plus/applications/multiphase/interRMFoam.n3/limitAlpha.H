// limitAlpha
{
	Info<< "Limiting VoF field and fluxes";

	scalar excessiveVol(sum(fvc::volumeIntegrate(max(0.0, alpha1-1))));
	reduce(excessiveVol, sumOp<scalar>());

	Info<< ", excessive volume: " << excessiveVol << " m3" << endl;

	// Limit volume fraction
	alpha1 = min( max(0.0, alpha1), 1.001);
	alpha2 = 1 - alpha1;

	// Limit face flux
	alphaPhi10 = sign(alphaPhi10)*
		min
		(
			mag(alphaPhi10),
			mag(phiCN)
		);
}

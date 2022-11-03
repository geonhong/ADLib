// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(U) + fvm::div(phi, U)
  + MRF.DDt(U)
  + turbulence->divDevReff(U)
 ==
    fvOptions(U)
);
fvVectorMatrix& UEqn = tUEqn.ref();

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
    solve(UEqn == -cellMask*fvc::grad(p));

    fvOptions.correct(U);
}

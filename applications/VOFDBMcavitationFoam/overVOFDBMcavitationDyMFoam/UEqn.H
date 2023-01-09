    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      - fvm::Sp(fvc::ddt(rho) + fvc::div(rhoPhi), U)
      + turbulence->divDevRhoReff(rho, U)
      + fvm::Sp(rho*coupling.damping(), U)
      - fvm::Su(cloud.momentumSource(), U)
      ==
        fvOptions(rho, U) 
        + rho*cloud.USource().oldTime() //for transition coupling
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            cellMask*fvc::reconstruct
            (
                (
                    interface.surfaceTensionForce()
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }

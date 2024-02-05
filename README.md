# ADLib
AeroDynamics Laboratory Libraries for the OpenFOAM(R)

This project contains the applications and libraries that have developed
in Aerodynamics Laboratory of Aerospace Engineering Department, KAIST.

## ADLib-plus
ADLib-plus is a ADLib package for ESI version of the OpenFOAM. It contains incompressible, compressible, and multiphase solvers and additional libraries. The libraries in ADLib is as following.

### Solvers
* Incompressible solvers
  - simpleFoam.n3 : simpleFoam with limiting velocity
  - pimpleFoam.n3 : pimpleFoam with limiting velocity and transimple algorithm
* Multiphase solvers
  - interFoam.n3 : two-phase solver with limiting velocity and flux limiter
* Compressible solvers
  - rhoCentralFoam : rhoCentralFoam with limiter on energy equation handler

### Libraries
* Additional function objects
  - nforces : evaluate force on a faceset
* Additional dynamic mesh libraries
  - Relative Motion : relative motion library that can work with MRF
* Turbulence models
  - kOmegaSSTN : kOmegaSST with rotation/curvature correction

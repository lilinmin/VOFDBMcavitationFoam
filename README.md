# VOFDBMcavitationFoam

## Features

- Volume-of-Fluid (VOF) solver with phasechange for cavitating flows coupled with Lagrangian library
- Conversion of small VOF elements to Lagrangian parcels to reduce computational cost
- Support for
  - Particle-particle interaction (collision, coalescence) and secondary breakup
  - Solving RP equation or simplified RP equation for discrete bubbles
  - Adaptive mesh refinement
  - Fully parallelized
- Based on interPhaseChangeFoam

## Compilation

- First, compile the library libVOFDBM in the src folder.
- Afterwards, compile the solver VOFDBMcavitationFoam in the application folder.

## Reference
Linmin Li, Weisen Xu, Bowen Jiang, Xiaojun Li, Zuchao Zhu, A multiscale Eulerian-Lagrangian cavitating flow solver in OpenFOAM, SoftwareX 

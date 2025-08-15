This repository contains code that simulates the physics of toppling dominoes as described in https://zmatt.net/domino-physics/

## Basic simulation
`domino-simulation.py` simulates a line of toppling dominoes according to the assumptions in the blog post (no-slip, etc).
* Run `python3 domino-simulation.py` to run a first simulation, and then gradually inspect the script and edit the variables at the top.

## LMGC90 simulation
For more complex 3D simulations, the `lmgc90` directory contains scripts for running simulations using LMGC90, based on the work by David Cantor and Kajetan Wojtacki at https://github.com/dacg/3D_dominoes .
* Download and build LMGC90 (https://git-xen.lmgc.univ-montp2.fr/lmgc90/lmgc90_user/-/wikis/home)
* Set `PYTHONPATH` to point to the LMGC90 build directory
* Create a new empty directory for the output data
* From the new directory, run (for example) `python3 ../generate.py` and `python3 ../run.py`.
* The output file `DISPLAY/lmgc90.pvd` can be opened in paraview.  In fact you can open it while the simulation is running and interrupt the simulation when time has progressed sufficiently.


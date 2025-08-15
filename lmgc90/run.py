#################################################
#############    Run simulation   ###############
#################################################
# Based on https://github.com/dacg/3D_dominoes
# by David Cantor and Kajetan Wojtacki

# Importing modules
import os
from pylmgc90.chipy import *
import numpy as np

Initialize()                                # Initializing
checkDirectories()                          # Checking/creating mandatory sub-folders
utilities_DisableLogMes()                   # Log Messages

################ Variables #################

# time evolution parameters
dt = 2.0e-4                                 # Time step
max_nb_steps = 50000                        # Steps

theta = 0.5                                 # Theta integrator parameter

deformable = 0                              # Deformable  yes=1, no=0

freq_detect = 1                             # Detection frequency
Rloc_tol = 5.e-2                            # Recycling forces from previous step with tolerance

PRPRx_UseCpCundallDetection(250)            # Detection method. Arg: # of iterations

# nlgs parameters for convergence of the mechanical system
type='Stored_Delassus_Loops         '
norm = 'QM/16'
tol = 1.3e-4
relax = 1.0
gs_it1 = 61
gs_it2 = 101

freq_write   = 1000                         # Writing frequency
freq_display = 1000                         # Display frequency
ref_radius = 2.5e-2                         # Contact forces scale for display

################ Simulation ################
# Set space dimension
SetDimension(3)

# Read and load
#
utilities_logMes('INIT TIME STEPPING')
TimeEvolution_SetTimeStep(dt)
Integrator_InitTheta(theta)
#
utilities_logMes('READ BEHAVIOURS')
ReadBehaviours()
#
utilities_logMes('READ BODIES')
ReadBodies()
#
utilities_logMes('LOAD BEHAVIOURS')
LoadBehaviours()
#
utilities_logMes('READ INI DOF')
ReadIniDof()
#
utilities_logMes('READ DRIVEN DOF')
ReadDrivenDof()
#
utilities_logMes('LOAD TACTORS')
LoadTactors()
#
utilities_logMes('READ INI Vloc Rloc')
ReadIniVlocRloc()

# Paranoid writes
utilities_logMes('WRITE BODIES')
WriteBodies()
utilities_logMes('WRITE BEHAVIOURS')
WriteBehaviours()
utilities_logMes('WRITE DRIVEN DOF')
WriteDrivenDof()


# Open display & postpro
utilities_logMes('DISPLAY & WRITE')
OpenDisplayFiles()
OpenPostproFiles()

# Since constant geometry of bodies, compute elementary mass just once
ComputeMass()

# Time loop
for k in range(1,max_nb_steps+1):

  IncrementStep()
  ComputeFext()                      # Computing external forces

  ComputeBulk()
  ComputeFreeVelocity()              # Bodies free velocities

  SelectProxTactors(freq_detect)     # Contact detection

  RecupRloc(Rloc_tol)                # Recycling forces

  ExSolver(type, norm, tol, relax, gs_it1, gs_it2)      # Solver

  UpdateTactBehav()                  # Updating contact behavior
  StockRloc()                        # Storing forces
  ComputeDof()                       # Updating DOFs
  UpdateStep()

  WriteOutDof(freq_write)            # Write DOFs
  WriteOutVlocRloc(freq_write)       # Write contact information

  WriteDisplayFiles(freq_display,ref_radius)    # Draw display
  WritePostproFiles()
#

# close display & postpro
CloseDisplayFiles()
ClosePostproFiles()

# End of simulation
Finalize()

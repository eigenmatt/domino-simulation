#################################################
####    Building the mechanical systems   #######
#################################################
# Based on https://github.com/dacg/3D_dominoes
# by David Cantor and Kajetan Wojtacki

# Importing modules
import os
import numpy as np
import math
from pylmgc90.pre import *

# Setting the separation between pieces
delta_do = 0.016

# Setting friction
fric_sd = 0.3
fric_dd = 0.3

# Creating folder to store the system
if not os.path.isdir('./DATBOX'):
  os.mkdir('./DATBOX')

# Setting the dimension
dim = 3

# Creating containers for
#   * bodies
bodies = avatars()
#   * materials
mat = materials()
#   * visibility between bodies
svs = see_tables()
#   * contact laws
tacts = tact_behavs()
#   * post-processing commands
post = postpro_commands()

# Creating materials
#   * For the dominoes
plex = material(name='PLEXx', materialType='RIGID', density=1000.)
mat.addMaterial(plex)
#   * For the base surface
tdur = material(name='TDURx', materialType='RIGID', density=1000.)
mat.addMaterial(tdur)

# Creating a mechanical rigid 3D model
mod = model(name='rigid', physics='MECAx', element='Rxx3D', dimension=dim)

# Setting the number of dominoes
nb_particles = 200

# Dimensions of a domino
do_x = 0.0075/2
do_y = do_x*3.2
do_z = do_y*2

# Array of domino vertices
vertices_do = np.array([[ do_x, -do_y, -do_z],
                        [ do_x, -do_y,  do_z],
                        [ do_x,  do_y,  do_z],
                        [ do_x,  do_y, -do_z],
                        [-do_x, -do_y, -do_z],
                        [-do_x, -do_y,  do_z],
                        [-do_x,  do_y,  do_z],
                        [-do_x,  do_y, -do_z]])


# Definition of velocity control for the pushing particle A
def imposeVelocity(t):
  if t <= 0.030:
    return 0.5
  else :
    return  0.0
  #
#

# Creating particle A
body = rigidPolyhedron(model=mod, material=plex, color='BLEUx', generation_type='vertices', vertices=vertices_do)
body.rotate(description='Euler', phi=math.pi/2., theta=math.pi/2., psi=0., center=[0., 0., 0.])
body.translate(dx = -1.3*do_z, dy=0., dz=do_z*2.+1.5*do_x)
body.imposeDrivenDof(component=[2, 3, 4, 5, 6], dofty='vlocy')
body.imposeDrivenDof(component=[1], description='evolution', evolutionFile='veloc_dom.dat', dofty='vlocy')
writeEvolution(f=imposeVelocity, instants=[0.0, 0.030, 0.0301, 9.9999], path='DATBOX/', name='veloc_dom.dat')
bodies += body
#

# Creating the dominoes
for i in range(0,nb_particles,1):
   # A rigid body
   body = rigidPolyhedron(model=mod, material=plex, color='BLEUx', generation_type='vertices', vertices=vertices_do)
   # Translating the body to its corresponding position
   body.translate(dx = delta_do*i, dz=do_z)
   # Adding it to the container
   bodies += body
#

# Setting the dimensions of the base surface
#lx = (nb_particles+10)*delta_do
lx = (nb_particles-1)*delta_do+40*do_x
ly = do_y*15
lz = do_x

# Creating the surface
down =rigidPlan(axe1=0.5*lx, axe2=0.5*ly, axe3=do_x, center=[(nb_particles-1)*delta_do/2, 0., -do_x], model=mod, material=tdur, color='VERTx')

# Freezing DOF
down.imposeDrivenDof(component=[1, 2, 3, 4, 5, 6], dofty='vlocy')

# Adding these bodies to the container
bodies.addAvatar(down)

# Managing the interactions
#  - between dominoes
lprpr=tact_behav(name='iqsc0', law='IQS_CLB', fric=fric_dd)
tacts+=lprpr
#  - between the dominoes and the surface
lprpl=tact_behav(name='iqsc1', law='IQS_CLB', fric=fric_sd)
tacts+=lprpl

# Visibility rules
#  - between dominoes
svprpr = see_table(CorpsCandidat='RBDY3', candidat='POLYR',
   colorCandidat='BLEUx', behav=lprpr, CorpsAntagoniste='RBDY3', 
   antagoniste='POLYR', colorAntagoniste='BLEUx', alert=0.2*do_x)
svs+=svprpr
#  - between the dominoes and the surface
svprpl = see_table(CorpsCandidat='RBDY3', candidat='POLYR',
   colorCandidat='BLEUx', behav=lprpl, CorpsAntagoniste='RBDY3', 
   antagoniste='PLANx', colorAntagoniste='VERTx', alert=0.2*do_x)
svs+=svprpl

post+=postpro_command("SOLVER INFORMATIONS")

# Writing the files. By default gravity is set in the z direction in SI units
writeBodies(bodies, chemin='DATBOX/')
writeBulkBehav(mat, chemin='DATBOX/', dim=dim)
writeTactBehav(tacts, svs, chemin='DATBOX/')
writeDrvDof(bodies, chemin='DATBOX/')
writeDofIni(bodies, chemin='DATBOX/')
writeVlocRlocIni(chemin='DATBOX/')
writePostpro(post, bodies, path='DATBOX/')

#try:
#  visuAvatars(bodies)
#except:
#  pass

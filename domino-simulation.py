#!/usr/bin/python3
#
# Physics simulation of toppling dominoes
# See https://zmatt.net/domino-physics/
# Matthew Chapman <matt@zmatt.net>
#

import numpy as np
import scipy
import matplotlib.pyplot as plt
from matplotlib import patches

# ==========
# parameters
# ==========
# length/height of domino
L=0.048
# width of domino
w=0.024
# thickness/depth of domino
d=0.0075
# density of domino (kg/m^3)
rho=1040
# gravity (ms^-2)
g=9.8
# domino pitch (spacing+thickness)
lmbda=0.016
# co-efficient of dynamic friction between dominoes
mu=0
# initial angular velocity from pushing
omega0=5*np.pi
# time step
dt=0.0001
# number of dominoes
total_dominoes=10
# maximum time steps
max_tsteps=100000
# screen, png, data or none
output='screen'
# display data every N time steps
output_interval=10
# print stats every N dominoes
stats_interval=10
# stepping method (euler or rk4)
step_function='euler'
# force solving force matrix, otherwise shortcuts are taken if no friction
force_solve_forces=False

# =================
# derived variables
# =================
# mass of domino
m=L*w*d*rho
# weight of domino
W=m*g
# distance from pivot to center of mass
rw=np.sqrt(L*L+d*d)/2
# absolute angle from pivot to center of mass, relative to vertical
phi=np.arctan(d/L)
# moment of inertia around pivot
I=m*(L*L+d*d)/3
# angle at which the next domino is hit
collision_angle=np.arcsin((lmbda-d)/L)
# angle at which domino comes to rest
final_angle=np.arccos(d/lmbda)

# =========
# functions
# =========
# calculate lean angle of each domino, assuming each domino is leaning on the next
def calc_theta(ndominoes, theta0):
  theta = np.empty(ndominoes)
  theta[0] = theta0
  for i in range(1,ndominoes):
    theta[i] = theta[i-1] + np.arcsin((lmbda*np.cos(theta[i-1])-d)/L)
  return theta

# calculate angular velocity of each domino (derivative of lean angle) 
def calc_omega(ndominoes, theta, omega0):
  omega = np.empty(ndominoes)
  omega[0] = omega0
  for i in range(1,ndominoes):
    omega[i] = omega[i-1]*(1-(lmbda/L)*np.sin(theta[i-1])/np.cos(theta[i]-theta[i-1]))
  return omega 

# calculate angular acceleration of each domino if acceleration of first domino is known
def calc_alpha(ndominoes, theta, omega, alpha0):
  alpha = np.empty(ndominoes)
  alpha[0] = alpha0
  for i in range(1,ndominoes):
    omega_factor = omega[i]/omega[i-1]
    sinsin = np.sin(theta[i-1])*np.sin(theta[i]-theta[i-1])
    cossq = np.cos(theta[i-1]-theta[i])**2
    alpha[i] = omega_factor*alpha[i-1] - (lmbda*omega[i-1]/L)*(omega[i-1]*np.cos(theta[i])+omega[i]*sinsin)/cossq
  return alpha

# solve for the acceleration of the leading domino, by considering normal forces and torque balances
def solve_alpha_newton(ndominoes, theta0, omega0):
  # possible optimization: we could calculate theta and omega for each domino on the fly in the loop
  # but since we have these functions available...
  theta = calc_theta(ndominoes, theta0)
  omega = calc_omega(ndominoes, theta, omega0)
  normal_coeff = 0
  normal_const = 0
  for i in range(ndominoes):
    # calculate acceleration co-efficients for this domino, relative to the leading domino
    # (alpha = alpha_coeff*alpha0+alpha_const)
    if i == 0:
       alpha_coeff = 1
       alpha_const = 0
    else:
       sin0 = np.sin(theta[i-1])
       cos1 = np.cos(theta[i])
       cos10 = np.cos(theta[i]-theta[i-1])
       sin10 = np.sin(theta[i]-theta[i-1])
       omega_ratio = 1-(lmbda/L)*sin0/cos10
       alpha_ratio = (lmbda/L)*(cos1+sin0*sin10*omega_ratio)/(cos10**2)
       alpha_coeff *= omega_ratio
       alpha_const *= omega_ratio
       alpha_const -= alpha_ratio*omega[i-1]*omega[i-1]

    # calculate normal force co-efficients for this domino, relative to the acceleration of the leading domino
    # (normal = normal_coeff*alpha0+normal_const)
    #
    # effective right torque arm (arm + mu*friction arm)
    # (for the first domino, the r0 arm is irrelevant)
    r0 = 0 if i == 0 else L*(np.cos(theta[i]-theta[i-1])+mu*np.sin(theta[i]-theta[i-1]))
    # effective left torque arm (arm - mu*friction arm)
    # (for the last domino, the r1 arm is irrelevant, we can use any non-zero value)
    r1 = L if i == ndominoes-1 else L*np.cos(theta[i+1])/np.cos(theta[i])-d*np.tan(theta[i])-mu*d
    normal_coeff *= r0/r1
    normal_coeff += (I/r1)*alpha_coeff
    normal_const *= r0/r1
    normal_const += (I/r1)*alpha_const
    normal_const -= (rw/r1)*W*np.sin(theta[i]-phi)
    # possible optimization: terminate loop once normal_const/normal_coeff converges to desired precision

  # now we know that normal=0 at the end of the chain, since there are no more dominoes, so we can solve for alpha0
  alpha0 = -normal_const/normal_coeff
  return alpha0

# solve for the post-collision velocities, by considering normal forces and angular impulse balances
def solve_collision_newton(ndominoes, theta0, omega0):
  # possible optimization: we could calculate theta and omega for each domino on the fly in the loop
  # but since we have these functions available...
  old_theta = calc_theta(ndominoes, theta0)
  old_omega = calc_omega(ndominoes, old_theta, omega0)
  theta = calc_theta(ndominoes+1, 0)
  rel_omega = calc_omega(ndominoes+1, theta, 1)
  impulse_coeff = 0
  impulse_const = 0
  for i in range(ndominoes+1):
    # calculate angular velocity change for this domino, relative to the imparted velocity of the new leading domino (0)
    # (delta_omega = delta_omega_coeff*newomega0+delta_omega_const)
    delta_omega_coeff = rel_omega[i]
    delta_omega_const = 0 if i == 0 else -old_omega[i-1]

    # calculate angular impulse co-efficients for this domino, relative to the imparted velocity of the new leading domino (0)
    # (impulse = impulse_coeff*newomega0+impulse_const)
    #
    # effective right torque arm (arm + mu*friction arm)
    # (for the first domino, the r0 arm is irrelevant)
    r0 = 0 if i == 0 else L*(np.cos(theta[i]-theta[i-1])+mu*np.sin(theta[i]-theta[i-1]))
    # effective left torque arm (arm - mu*friction arm)
    # (for the last domino, the r1 arm is irrelevant, we can use any non-zero value)
    r1 = L if i == ndominoes else L*np.cos(theta[i+1])/np.cos(theta[i])-d*np.tan(theta[i])-mu*d
    impulse_coeff *= r0/r1
    impulse_coeff += (I/r1)*delta_omega_coeff
    impulse_const *= r0/r1
    impulse_const += (I/r1)*delta_omega_const
    # possible optimization: terminate loop once impulse_const/impulse_coeff converges to desired precision

  # now we know that impulse=0 at the end of the chain, since there are no more dominoes, so we can solve for newomega0
  newomega0 = -impulse_const/impulse_coeff
  return newomega0

# solve for the acceleration of the leading domino, using explicit equation from Lagrangian solution
# (should be equivalent to solve_alpha_newton, but friction is not implemented here)
def solve_alpha_lagrange(ndominoes, theta0, omega0):
  assert(mu == 0) # friction not implemented here
  # possible optimization: we could calculate theta for each domino on the fly, but since we have this function available...
  theta = calc_theta(ndominoes, theta0)
  dPE_domega0 = calc_PE_deriv(ndominoes, theta)

  rel_omega = calc_omega(ndominoes, theta, 1)
  rel_omega_sq = rel_omega**2
  Ck = np.cumsum(rel_omega_sq[::-1])[::-1]

  rhs = 0
  for i in range(1,ndominoes):
    sin0 = np.sin(theta[i-1])
    cos1 = np.cos(theta[i])
    cos10 = np.cos(theta[i]-theta[i-1])
    sin10 = np.sin(theta[i]-theta[i-1])
    omega_ratio = 1-(lmbda/L)*sin0/cos10
    alpha_ratio = (lmbda/L)*(cos1+sin0*sin10*omega_ratio)/(cos10**2)
    C = Ck[i]/omega_ratio
    rhs += C*rel_omega[i-1]*alpha_ratio
  rhs *= I*omega0*omega0
  rhs -= dPE_domega0
  alpha0 = rhs/(Ck[0]*I)
  #print(alpha0)
  return alpha0

# solve for the post-collision velocities, assuming no friction
def solve_collision_theoretical(ndominoes, theta0, omega0):
  assert(mu == 0)  # friction not implemented here
  theta = calc_theta(ndominoes+1, 0)
  rel_omega = calc_omega(ndominoes+1, theta, 1)
  newK = sum(rel_omega**2)
  newomega0 = ((newK-1)/newK)*omega0  # =oldK/(oldK+1)
  return newomega0

# evaluate kinetic energy of system
def calc_KE(ndominoes, theta0, omega0):
  theta = calc_theta(ndominoes, theta0)
  omega = calc_omega(ndominoes, theta, omega0)
  return 0.5*I*sum(omega**2)

# evaluate potential energy of system
# equivalent to 0.5*m*g*L*cos(theta)+0.5*m*g*d*sin(theta)
def calc_PE(ndominoes, theta0):
  theta = calc_theta(ndominoes, theta0)
  return m*g*rw*np.sum(np.cos(theta-phi)) + (total_dominoes-ndominoes)*m*g*L/2

# evaluate derivative of potential energy
def calc_PE_deriv(ndominoes, theta):
  rhs = m*g*rw*np.sin(theta[0]-phi)
  derivative = 1
  for i in range(1,ndominoes):
    sin0 = np.sin(theta[i-1])
    cos10 = np.cos(theta[i]-theta[i-1])
    omega_ratio = 1-(lmbda/L)*sin0/cos10
    derivative *= omega_ratio
    rhs += m*g*rw*np.sin(theta[i]-phi)*derivative
  return -rhs

# evaluate normals assuming alpha is known
def calc_normals_forward(ndominoes, theta, alpha):
  normals = np.empty(ndominoes)
  for i in range(ndominoes-1):
    if i == 0:
      t0 = 0
    else:
      r0 = L*(np.cos(theta[i]-theta[i-1])+mu*np.sin(theta[i]-theta[i-1]))
      t0 = normals[i-1]*r0
    tw = -rw*W*np.sin(theta[i]-phi)
    r1 = L*np.cos(theta[i+1])/np.cos(theta[i])-d*np.tan(theta[i])-mu*d
    normals[i] = (t0+tw+I*alpha[i])/r1
  normals[ndominoes-1] = 0
  return normals

# evaluate normals assuming alpha is known (workiing backwards starting from domino N-1)
# should be equivalent to calc_normals_forward, but with many dominoes this may have better numerical stablility
def calc_normals_reverse(ndominoes, theta, alpha):
  normals = np.empty(ndominoes)
  normals[ndominoes-1] = 0
  for i in range(ndominoes-1,0,-1):
    if i == ndominoes-1:
      t1 = 0
    else:
      r1 = L*np.cos(theta[i+1])/np.cos(theta[i])-d*np.tan(theta[i])-mu*d
      t1 = normals[i]*r1
    tw = rw*W*np.sin(theta[i]-phi)
    r0 = L*(np.cos(theta[i]-theta[i-1])+mu*np.sin(theta[i]-theta[i-1]))
    normals[i-1] = (t1+tw-I*alpha[i])/r0
  return normals

# verify correct solution of force equation matrix
def check_residuals(ndominoes, theta, alpha, normals):
  residuals = np.empty(ndominoes)
  for i in range(ndominoes):
    N0 = 0 if i == 0 else normals[i-1]
    N1 = normals[i]
    r0 = 0 if i == 0 else L*(np.cos(theta[i]-theta[i-1])+mu*np.sin(theta[i]-theta[i-1]))
    r1 = 0 if i == ndominoes-1 else L*np.cos(theta[i+1])/np.cos(theta[i])-d*np.tan(theta[i])-mu*d
    residuals[i] = N0*r0-N1*r1-rw*W*np.sin(theta[i]-phi)+I*alpha[i]
  return residuals

# check required co-efficient of friction with floor for each domino
def check_floor_friction(ndominoes, theta0, omega0, alpha0, warn_threshold):
  theta = calc_theta(ndominoes, theta0)
  omega = calc_omega(ndominoes, theta, omega0)
  alpha = calc_alpha(ndominoes, theta, omega, alpha0)
  normals = calc_normals_reverse(ndominoes, theta, alpha)
  assert(not any([N < 0 for N in normals]))
  assert(not any([e > 1e-12 for e in check_residuals(ndominoes, theta, alpha, normals)]))
  friction_coefficients = np.empty(ndominoes)
  for i in range(ndominoes):
    N0x = 0 if i == 0 else normals[i-1]*np.cos(theta[i-1])
    N0y = 0 if i == 0 else normals[i-1]*np.sin(theta[i-1])
    N1x = normals[i]*np.cos(theta[i])
    N1y = normals[i]*np.sin(theta[i])
    Ffloor = N0x - N1x + m*alpha[i]*rw*np.cos(theta[i]-phi)
    Nfloor = N1y - N0y - m*alpha[i]*rw*np.sin(theta[i]-phi) + W
    friction_coefficients[i] = Ffloor/Nfloor # signed
    if abs(friction_coefficients[i]) > warn_threshold:
      print("warning: domino %d of %d requires friction co-efficient %f to not slide" % (i+1, ndominoes, abs(friction_coefficients[i])))
  return friction_coefficients

# one time step using simple Euler stepping
def step_simulation_euler(ndominoes, theta, omega):
  alpha = solve_alpha_f(ndominoes, theta, omega)
  newtheta = theta + omega*dt
  newomega = omega + alpha*dt
  return (newtheta, newomega, alpha)

# one time step using coupled fourth-order Runge-Kutta
# see https://math.stackexchange.com/questions/2023819/using-the-runge-kuttas-method-to-solve-a-2nd-derivative-question
def step_simulation_rk4(ndominoes, theta, omega):
  C0 = omega
  K0 = solve_alpha_f(ndominoes, theta, omega)
  C1 = omega + 0.5*dt*K0
  K1 = solve_alpha_f(ndominoes, theta + 0.5*dt*C0, C1)
  C2 = omega + 0.5*dt*K1
  K2 = solve_alpha_f(ndominoes, theta + 0.5*dt*C1, C2)
  C3 = omega + dt*K2
  K3 = solve_alpha_f(ndominoes, theta + dt*C2, C3)
  thetadot = (C0+2*C1+2*C2+C3)/6
  alpha = (K0+2*K1+2*K2+K3)/6
  newtheta = theta + thetadot*dt
  newomega = omega + alpha*dt
  return (newtheta, newomega, alpha)

# print data at time step
def print_data(tstep, t, ndominoes, theta, omega, alpha):
  print("TIME %f %f %f %f %g %g" % (t, theta, omega, alpha, calc_KE(ndominoes, theta, omega), calc_PE(ndominoes, theta)))

# ========================================
# plot implementation (modify as required)
# ========================================
class Plot:
  def __init__(self, total_dominoes):
    fig, ax = plt.subplots(1, 1, figsize=(6.4,3.6), dpi=300)
    ax.set_aspect('equal')
    ax.set_xlim([-0.02,0.25]) # tweak this as needed, maybe should have a better default
    ax.set_ylim([-0.025,2.1*rw+0.005])
    ax.axhline(-0.0005, zorder=-1, color='k', linewidth=1)
    ax.set_axis_off()
    dominoes = []
    for i in range(total_dominoes):
      domino = ax.add_patch(patches.Rectangle(xy=(i*lmbda, 0), width=-d, height=L, edgecolor='black', facecolor='#ff5f1f80'))
      dominoes.append(domino)
    self.dominoes = dominoes
    fig.tight_layout()

  def update(self, tstep, t, ndominoes, theta0, omega0, alpha):
    theta = calc_theta(ndominoes, theta0)
    omega = calc_omega(ndominoes, theta, omega0)
    vlines = []
    vline_colours = []
    hlines = []
    for i in range(ndominoes):
      j = ndominoes-i-1
      angle = np.degrees(theta[i])
      self.dominoes[j].set_angle(-angle)

    if output == 'png':
      plt.savefig('frame%04d.png' % tstep)
    else:
      plt.draw()
      plt.pause(0.001)

# ======================================================
# assign function pointers (apologies to Python purists)
# ======================================================
if output == 'screen' or output == 'png':
  plot = Plot(total_dominoes)
  output_f = plot.update
elif output == 'data':
  output_f = print_data
else:
  output_f = None

if mu != 0 or force_solve_forces:
  solve_alpha_f = solve_alpha_newton
  solve_collision_f = solve_collision_newton
else:
  solve_alpha_f = solve_alpha_lagrange
  solve_collision_f = solve_collision_theoretical

if step_function == 'rk4':
  step_simulation_f = step_simulation_rk4
else:
  step_simulation_f = step_simulation_euler

# ====================
# main simulation loop
# ====================
def run_simulation():
  ndominoes = 1
  theta = 0
  omega = omega0  # start with some angular velocity from pushing
  write_step = 0
  stats_start_t = 0
  for tstep in range(max_tsteps):
    (newtheta, newomega, alpha) = step_simulation_f(ndominoes, theta, omega)
    if output_f != None and tstep % output_interval == 0:
      output_f(write_step, tstep*dt, ndominoes, theta, omega, alpha)
      write_step += 1
  
    theta = newtheta
    omega = newomega
  
    if theta >= collision_angle:
      if ndominoes == total_dominoes:
        break
      if ndominoes % stats_interval == 0:
        t = tstep*dt
        dps = stats_interval/(t-stats_start_t)
        print("STATS: average collisions per second: %f (%f ms-1)" % (dps, dps*lmbda))
        stats_start_t = t
      newomega = solve_collision_f(ndominoes, theta, omega)
      ndominoes += 1
      theta = 0
      omega = newomega

if  __name__ == '__main__':
  run_simulation()


# Example Rayleigh-Benard Convection
"""
Dedalus script for 2D Boussinesq Rayleigh-Benard convection.

This script uses a Fourier basis in the x direction with periodic boundary
conditions.

Adapted from the 2D Boussinesq Rayleigh-Benard Convection example
as seen in the Dedalus installation.

Equations are non-dimensionalised by the viscous timescale (t in units of viscous time).

Can be run using mpi (on for example, 4 processors) by typing

    mpiexec -n 4 python3 rayleigh_benard.py

instead of simply

    python3 rayleigh_benard.py
"""

import numpy as np
from mpi4py import MPI
import time
import matplotlib.pyplot as plt
import sys

from dedalus import public as de
from dedalus.extras import flow_tools
import pathlib

import logging
logger = logging.getLogger(__name__)

import run_param_file as rpf   # Imports a parameter file "run_param_file.py"

save_direc = "raw_data2/"
pathlib.Path(save_direc).mkdir(parents=True, exist_ok=True)


# Model Parameters
Lx, Lz = rpf.Lx, rpf.Lz
Nx, Nz = rpf.Nx, rpf.Nz
Pr = rpf.Pr
Ra = rpf.Ra

# Create bases and domain
x_basis = de.Fourier('x', Nx, interval=(0, Lx), dealias=3/2)   # Fourier basis in the x
z_basis = de.Chebyshev('z', Nz, interval=(0, Lz), dealias=3/2) # Chebyshev basis in the z
domain = de.Domain([x_basis, z_basis], grid_dtype=np.float64)  # Defining our domain
z = domain.grid(1, scales=1)                                   # accessing the z values

# 2D Boussinesq hydrodynamics
problem = de.IVP(domain, variables=['p', 'T', 'u', 'w', 'Tz', 'uz', 'wz'])
problem.meta['p','T','u','w']['z']['dirichlet'] = True

# Defining model parameters
problem.parameters['Lx'] = Lx
problem.parameters['Lz'] = Lz
problem.parameters['Ra'] = Ra
problem.parameters['Pr'] = Pr
problem.parameters['X'] = Ra/Pr


# Defining d/dz of T, u, and w for reducing our equations to first order
problem.add_equation("Tz - dz(T) = 0")
problem.add_equation("uz - dz(u) = 0")
problem.add_equation("wz - dz(w) = 0")

# mass continuity
problem.add_equation("  dx(u) + wz = 0 ")
# x-component of the momentum equation
problem.add_equation("  dt(u) - dx(dx(u)) - dz(uz) + dx(p) = - ( u*dx(u) + w*uz )  ")
# z-component of the momentum equation
problem.add_equation("  dt(w) - X*T - dz(wz) - dx(dx(w)) + dz(p)  = -( u*dx(w) + w*wz ) ")
# Temperature equation
problem.add_equation("  Pr*dt(T) - dx(dx(T)) - dz(Tz) = -Pr*( u*dx(T) + w*Tz )   ")

# Can change these depending on what we want to do. Just have to obey some mathematical conditions. When dealing with a partial derivate the only way to get an exact solution is to specify exact boundary conditions (at boundary of the domain).
problem.add_bc("left(w) = 0")            # Impermeable bottom boundary
problem.add_bc("right(w) = 0", condition="(nx != 0)")   # Impermeable top boundary
problem.add_bc("right(p) = 0", condition="(nx == 0)")   # Required for equations to be well-posed - see https://bit.ly/2nPVWIg for a related discussion
problem.add_bc("left(uz) = 0")           # Stress-free bottom boundary
problem.add_bc("right(uz) = 0")          # Stress-free top boundary
problem.add_bc("right(T) = 0")           # Fixed temperature at upper boundary
problem.add_bc("left(Tz) = -1")           # Fixed flux at bottom boundary, F = F_cond


# Build solver
solver = problem.build_solver(de.timesteppers.RK222)
logger.info('Solver built')

# Initial conditions
z = domain.grid(1)
T = solver.state['T']
Tz = solver.state['Tz']

# Random perturbations, initialized globally for same results in parallel
gshape = domain.dist.grid_layout.global_shape(scales=1)
slices = domain.dist.grid_layout.slices(scales=1)
rand = np.random.RandomState(seed=42)
noise = rand.standard_normal(gshape)[slices]

# Linear background + perturbations damped at walls
zb, zt = z_basis.interval
pert =  1e-5 * noise * (zt - z) * (z - zb)
T['g'] = pert
T.differentiate('z', out=Tz)

# Initial timestep
dt = rpf.initial_timestep

# Integration parameters --- Note if these are all set to np.inf, simulation will perpetually run.
solver.stop_sim_time = rpf.end_sim_time
solver.stop_wall_time = rpf.end_wall_time
solver.stop_iteration = rpf.end_iterations

# CFL criterion
CFL = flow_tools.CFL(solver, initial_dt=dt, cadence=10, safety=0.5,
                     max_change=1.5, min_change=0.5, max_dt=rpf.max_dt, threshold=0.05)
CFL.add_velocities(('u', 'w'))

# Flow properties
flow = flow_tools.GlobalFlowProperty(solver, cadence=10)
flow.add_property("sqrt(u*u + w*w)", name='Re')

# Saving snapshots
snapshots = solver.evaluator.add_file_handler(save_direc + 'snapshots', sim_dt=rpf.snapshot_freq, max_writes=50)
snapshots.add_system(solver.state)

# Analysis tasks
analysis = solver.evaluator.add_file_handler(save_direc + 'analysis', sim_dt=rpf.analysis_freq, max_writes=5000)
analysis.add_task("integ(T,'x')/Lx", layout='g', name='<T>_x')
analysis.add_task("integ(Tz,'x')/Lx", layout='g', name='<Tz>_x')
'''
The add_task bit uses the same parser as the bit used on the equation. Line 134 says calculate the integral of T in the x direction divided by Lx (makes it a definite integral) so gives the mean temperature in the x direction. 
Line 135 does the integral of Tz (dt/dz) in the x direction divided by Lx (the length of the x direction so divide the definite integral by the length of x to get the mean.)
'''

# Mean Re
analysis.add_task("integ( integ( sqrt(u*u + w*w) , 'x')/Lx, 'z')/Lz", layout='g', name='Re')

# Flux decomposition - Internal energy equation
analysis.add_task("integ(T*w,'x')*Pr/Lx", layout='g', name='L_conv')
analysis.add_task("integ((-1)*Tz, 'x')/Lx", layout='g', name='L_cond')
'''
Flux is calculated by the number density multiplied by the velocity with which you are moving through the density of thing. Energy flux is the energy density multiplied by the velocity. T is the fluctuating temperature field, multiplied by the vertical velocity. 
'''

# Mean KE
analysis.add_task(" integ( (integ(0.5*(u*u + w*w),'x')/Lx), 'z')/Lz", layout='g', name='KE')

# Creating a parameter file
run_parameters = solver.evaluator.add_file_handler(save_direc + 'run_parameters', wall_dt=1e20, max_writes=1)
run_parameters.add_task(Lx, name="Lx")
run_parameters.add_task(Lz, name="Lz")
run_parameters.add_task(Ra, name="Ra")
run_parameters.add_task(Pr, name="Pr")
run_parameters.add_task(Nx, name="Nx")
run_parameters.add_task(Nz, name="Nz")

run_parameters.add_task(rpf.snapshot_freq, name="snap_freq")
run_parameters.add_task(rpf.analysis_freq, name="ana_freq")
run_parameters.add_task(rpf.max_dt,        name="max_dt")

# Main loop
try:
    logger.info('Starting loop')
    start_time = time.time()
    while solver.ok:
        dt = CFL.compute_dt()
        dt = solver.step(dt)

        if (solver.iteration) == 1:
            # Prints various parameters to terminal upon starting the simulation
            logger.info('Parameter values imported form run_param_file.py:')
            logger.info('Lx = {}, Lz = {}; (Resolution of {},{})'.format(Lx, Lz, Nx, Nz))
            logger.info('Ra = {}, Pr = {}'.format(Ra, Pr))
            logger.info('Snapshot files outputted every {}'.format(rpf.snapshot_freq))
            logger.info('Analysis files outputted every {}'.format(rpf.analysis_freq))
            if rpf.end_sim_time != np.inf:
                logger.info('Simulation finishes at sim_time = {}'.format(rpf.end_sim_time))
            elif rpf.end_wall_time != np.inf:
                logger.info('Simulation finishes at wall_time = {}'.format(rpf.end_wall_time))
            elif rpf.end_iterations != np.inf:
                logger.info('Simulation finishes at iteration {}'.format(rpf.end_iterations))
            else:
                logger.info('No clear end point defined. Simulation may run perpetually.')

        if (solver.iteration-1) % 10 == 0:
            # Prints progress information include maximum Reynolds number every 10 iterations
            logger.info('Iteration: %i, Time: %e, dt: %e' %(solver.iteration, solver.sim_time, dt))
            logger.info('Max Re = %f' %flow.max('Re'))

except:
    logger.error('Exception raised, triggering end of main loop.')
    raise
finally:
    # Prints concluding information upon reaching the end of the simulation.
    end_time = time.time()
    logger.info('Iterations: %i' %solver.iteration)
    logger.info('Sim end time: %f' %solver.sim_time)
    logger.info('Run time: %.2f sec' %(end_time-start_time))
    logger.info('Run time: %f cpu-hr' %((end_time-start_time)/60/60*domain.dist.comm_cart.size))

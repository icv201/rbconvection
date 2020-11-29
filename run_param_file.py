"""
Parameter file for use in the Dedalus 2D anelastic convection script.
"""

import numpy as np

Lx, Lz = 4, 1                       # Domain size
Nx, Nz = 256, 64                    # Number of
Pr = 1.                             # Prandtl number
Pm = 1.                             # Magnetic Prandtl number
Ra = 1.6e3                          # Rayleigh number
Np = 0                              # Number of density scale heights
Q = 0                               # Chandrasekhar number
m = 1.5                             # Polytropic index
theta = 1 - np.exp(-Np/m)           # Dimensionaless inverse T scale height

initial_timestep = 1.5e-5                 # Initial timestep
max_dt = 1e-4                         # max dt

snapshot_freq = 1.5e-3              # Frequency snapshot files are outputted
analysis_freq = 1.5e-4              # Frequency analysis files are outputted

end_sim_time = 7.                   # Stop time in simulations units
end_wall_time = np.inf              # Stop time in wall time
end_iterations = np.inf             # Stop time in iterations

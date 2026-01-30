import numpy as np
import pandas as pd


num_particles = 100000
x = np.random.uniform(-5, 5, num_particles)
y = np.random.uniform(-5, 5, num_particles)
z = np.random.uniform(-5, 5, num_particles)

charge = np.random.uniform(-1, 1, num_particles)

np.savetxt("rand_100000.txt", np.column_stack((x, y, z)), fmt="%.6f", delimiter=' ')
np.savetxt("lambda_100000.txt", charge, fmt="%.6f", delimiter=' ')

import numpy as np
import matplotlib.pyplot as plt
import potentials as pot
import energies as ene
import integrator as itg
import initial_conditions as init
import poincare_sections as pcs

h = 0.01
N_inter = 30000
eps = 1e-2

x_0 = 0.0
y_vals = np.linspace(0, 0.6, 50)
E = 0.125  # 1/12

# Initialize lists to store all Poincaré points across different initial conditions
all_pcs_pos_y = []
all_pcs_vel_y = []


for i in range(len(y_vals)):
    y0 = y_vals[i]
    W_part = init.one_part(x_0, y0, 0, 0)
    POT = pot.hh_potential(W_part)[0]
    u_0 = np.sqrt(2 * (E - POT))
    v_0 = 0.0


    W_part[1, 0] = u_0
    W_part[1, 1] = v_0

    # Perform integration
    pos_t, positions = itg.rk4(0, W_part, h, N_inter, pot.hh_evolution)

    # Extract positions and velocities
    pos_x = positions[:, 0, 0]
    pos_y = positions[:, 0, 1]
    vel_x = positions[:, 1, 0]
    vel_y = positions[:, 1, 1]

    # Find Poincaré section points for the current initial condition
    pcs_pos_y, pcs_vel_y = pcs.pcs_find(pos_x, pos_y, vel_x, vel_y)

    # Append the current Poincaré section points to the overall lists
    all_pcs_pos_y.extend(pcs_pos_y)
    all_pcs_vel_y.extend(pcs_vel_y)

# Plotting all Poincaré section points on the same graph
fig, ax = plt.subplots()
ax.scatter(all_pcs_pos_y, all_pcs_vel_y, s=2, color="C4")

# Set labels and legend
ax.set_xlabel("coordinate y")
ax.set_ylabel("velocity v")
ax.set_title("Combined Poincaré Section for Multiple Initial y0 Values")

plt.show()

import numpy as np
import matplotlib.pyplot as plt
import potentials as pot
import integrator as itg
import initial_conditions as init

# Parameters
h = 0.01  # Step size for integrator
N_inter = 500 # Number of iterations for distance calculation
eps = 1e-3  # Perturbation distance
mu_c = 0.6  # Threshold for classifying ergodic vs regular
energy_values = np.array([0.01, 0.02, 0.04, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18])
N_PART = 100

def classify_point(P1, P1_prime):
    pos_t1, positions1 = itg.rk4(0, P1, h, N_inter, pot.hh_evolution)
    pos_t2, positions2 = itg.rk4(0, P1_prime, h, N_inter, pot.hh_evolution)

    # Compute the sum of squared distances (μ)
    mu = 0
    for i in range(N_inter):
        y1, y_dot1 = positions1[i, 0, 1], positions1[i, 1, 1]
        y2, y_dot2 = positions2[i, 0, 1], positions2[i, 1, 1]
        x1, x_dot1 = positions1[i, 0, 0], positions1[i, 1, 0]
        x2, x_dot2 = positions2[i, 0, 0], positions2[i, 1, 0]

        squared_distance = (np.sqrt((x2 - x1)**2 + (y2 - y1) ** 2 + (x_dot2 - x_dot1) ** 2 +(y_dot2 - y_dot1) ** 2))**2
        mu += squared_distance

    return mu

# Analyze for different energy values
relative_areas = []

for E in energy_values:
    # Generate initial conditions for the given energy level
    initial_conditions = init.n_energy_part(pot.hh_potential, N=N_PART, E=E)
    initial_conditions_eps = initial_conditions.copy()
    initial_conditions_eps[0,1] += eps
    n_regular = 0

    # Loop through each initial condition
    for i in range(N_PART):
        P1 = initial_conditions[:, :, i]
        # Perturb the initial point slightly
        P1_prime = initial_conditions_eps[:,:,i]

        # Classify the point based on the evolution
        mu = classify_point(P1, P1_prime)
        np.nan_to_num(mu, copy=False, nan=1.0)

        print(f"the mu value is {mu}")
        if mu < mu_c:
            n_regular += 1

    # Compute the relative area covered by regular points
    relative_area = n_regular / N_PART
    relative_areas.append(relative_area)

# Plot the results
plt.figure()
plt.plot(energy_values, relative_areas, marker='o', color='C0')
plt.xlabel("Energy (E)")
plt.ylabel("Relative Area Covered by Regular Curves")
plt.title("Transition from Ergodic to Stable Regime")
plt.grid()
plt.show()

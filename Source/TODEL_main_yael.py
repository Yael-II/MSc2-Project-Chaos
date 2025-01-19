import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

import potentials as pot
import energies as ene
import integrator as itg
import initial_conditions as init
import poincare_sections as pcs

if "YII_light_1" in plt.style.available: 
    plt.style.use("YII_light_1")
# Loading matplotlib style...
# PARAM
N_grid = 1000
N_inter = int(1e5)
N_part = 200
#all_E = np.array([1/100, 1/12, 1/10, 1/8, 1/6])
all_E = np.linspace(1/100, 1/6, 20)
h = 0.005

mu_c = 1e-4
d_12 = 1e-7

# INIT
W_grid = init.mesh_grid(N_grid, xmin=-1, xmax=1, ymin=-1, ymax=1)
X_grid = W_grid[0]
Y_grid = W_grid[1]
potential = pot.hh_potential(W_grid, position_only=True)

# MAIN
def gen_mu(E):
    pot_valid = np.ma.masked_where(potential > E, potential)
    W_1 = init.n_energy_part(pot.hh_potential, N_part, E)
    W_2 = np.zeros_like(W_1)
    alpha = np.random.uniform(0, 2*np.pi, N_part)
    W_2[0, 0] = W_1[0, 0] + d_12*np.cos(alpha)
    W_2[0, 1] = W_1[0, 1] + d_12*np.sin(alpha)
    W_2[1, 0] = W_1[1, 0]
    W_2[1, 1] = W_1[1, 1]

    t_1, positions_1 = itg.rk4(0, W_1, h, N_inter, pot.hh_evolution)
    x_1 = positions_1[:, 0, 0]
    y_1 = positions_1[:, 0, 1]
    u_1 = positions_1[:, 1, 0]
    v_1 = positions_1[:, 1, 1]
    
    t_2, positions_2 = itg.rk4(0, W_2, h, N_inter, pot.hh_evolution)
    x_2 = positions_2[:, 0, 0]
    y_2 = positions_2[:, 0, 1]
    u_2 = positions_2[:, 1, 0]
    v_2 = positions_2[:, 1, 1]
    dist_sq = (x_2[-25:] - x_1[-25:])**2 \
            + (y_2[-25:] - y_1[-25:])**2 \
            + (u_2[-25:] - u_1[-25:])**2 \
            + (v_2[-25:] - v_1[-25:])**2

    mu = np.sum(dist_sq, axis=0)
    return mu

A = np.zeros_like(all_E)
MU = []
ALL_E = []
for i in range(len(all_E)):
    mu = gen_mu(all_E[i])
    n_total = N_part
    n_ergo = np.sum(mu < mu_c) # count how many times the condition is verified
    A[i] = n_ergo/n_total
    MU.append(mu)
    ALL_E.append([all_E[i]]*len(mu))

MU = np.array(MU)
ALL_E = np.array(ALL_E)

# def lin(x, a, b):
#     return a*x + b

# w_chaos = np.argwhere(A < 1).flatten()
# X = all_E[w_chaos]
# Y = A[w_chaos]

# popt, pcov = curve_fit(lin, X, Y)
# a, b = popt
# da, db = np.sqrt(np.diag(pcov))

fig, ax = plt.subplots(1)
ax.scatter(ALL_E, MU, s=1, 
           color="k", alpha=0.1, label="Data")
ax.scatter(all_E, np.mean(MU, axis=1), s=10, 
           color="C3", marker="*", label="Mean")
ax.scatter(all_E, np.median(MU, axis=1), s=10, 
           color="C2", marker="s", label="Median")
ax.hlines(mu_c, np.min(all_E), np.max(all_E), color="C4")
ax.legend()
ax.set_yscale("log")
ax.set_xlabel("energy E")
ax.set_ylabel("quantity mu")

fig.tight_layout()
# fig.savefig("mu")

fig, ax = plt.subplots(1)
ax.scatter(all_E, A, color="C0", s=5, marker="o", label="Data")
ax.set_xlabel("energy E")
ax.set_ylabel("area A")
ax.legend()
fig.tight_layout()
# fig.savefig("area")
plt.show()



    

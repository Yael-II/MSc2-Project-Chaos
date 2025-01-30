"""
Test: Integrators

Demonstrating Keplerian 2-body orbits using various integrators,
and comparing accuracy and runtime over a range of step sizes.

@ Author: Moussouni, Yaël (MSc student) & Bhat, Junaid Ramzan (MSc student)
@ Institution:  Université de Strasbourg, CNRS, Observatoire astronomique
                de Strasbourg, UMR 7550, F-67000 Strasbourg, France
@ Date: 2025-01-01

Licence:
Order and Chaos in a 2D potential
Copyright (C) 2025 Yaël Moussouni (yael.moussouni@etu.unistra.fr)
                   Bhat, Junaid Ramzan (junaid-ramzan.bhat@etu.unistra.fr)

test_integrators.py
Copyright (C) 2025 Bhat, Junaid Ramzan (junaid-ramzan.bhat@etu.unistra.fr)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see https://www.gnu.org/licenses/.
"""

import numpy as np
import matplotlib.pyplot as plt
import time

import integrator as itg
import initial_conditions as init
import potentials as pot
import energies as ene

from matplotlib.patches import ConnectionPatch

if "YII_1" in plt.style.available: plt.style.use("YII_1")

# ----------------------------
# 1. Setup & global parameters
# ----------------------------
t0 = 0.0
T_final = 8.0
W0 = init.one_part(1, 0, 0, 1)  # [x, y, vx, vy]

h_range = np.logspace(-3.5, -0.1, 25)
h_range = np.append(h_range, 0.001)

# For plotting lines/colors
methods = [
    ("Euler", itg.euler, 'o--', 'C0'),
    ("RK2",   itg.rk2,   's--', 'C2'),
    ("RK4",   itg.rk4,   '^--', 'C3')
]
colors = {'Analytical': 'k', 
          'Euler': 'C0', 
          'RK2': 'C2',
          'RK4': 'C3'}

# Compute machine epsilon
eps = 1.0
while 1.0 + eps/2 > 1.0:
    eps /= 2.0
print(f"Machine epsilon: {eps}")

# Arrays to store final energy errors & times
err_euler, err_rk2, err_rk4 = [], [], []
time_euler, time_rk2, time_rk4 = [], [], []

# ------------------------------------------
# 2. Main loop over step sizes h in h_range
# ------------------------------------------

fig, ax = plt.subplots(1)

for h in h_range:
    ax.cla()
    N = int(T_final / h)

    # Analytical solution
    t_ana, W_ana = itg.kepler_analytical(t0, W0, h, N)
    W_ana_E = np.swapaxes(W_ana, 0, 2)
    W_ana_E = np.swapaxes(W_ana_E, 0, 1)
    E_analytical_final = ene.total(W_ana_E, pot.kepler_potential)

    # Numerical integrators + timing
    all_solutions = {}
    for (label, method, *_), store_err, store_t in zip(
        methods,
        [err_euler, err_rk2, err_rk4],
        [time_euler, time_rk2, time_rk4]
    ):
        start_time = time.time()
        t_num, W_num = itg.integrator_type(t0, W0, h, N, pot.kepler_evolution, method)
        elapsed = time.time() - start_time

        store_t.append(elapsed)
        
        # Final energy error
        W_num_E = np.swapaxes(W_num, 0, 2)
        W_num_E = np.swapaxes(W_num_E, 0, 1)
        E_numerical_final = ene.total(W_num_E, pot.kepler_potential)
        store_err.append(np.max(abs(E_analytical_final - E_numerical_final)))

        all_solutions[label] = W_num

    # Orbit plot (optional, can comment out if too many figures)
    eu_vals  = all_solutions["Euler"]
    rk2_vals = all_solutions["RK2"]
    rk4_vals = all_solutions["RK4"]


    ax.plot(W_ana[:, 0, 0], 
            W_ana[:, 0, 1],
            "-.", 
            color=colors['Analytical'],
            label="Analytical", 
            zorder=4)
    ax.plot(eu_vals[:, 0, 0],  
            eu_vals[:, 0, 1],  
            "-",
            color=colors['Euler'],  
            label="Euler")
    ax.plot(rk2_vals[:, 0, 0],
            rk2_vals[:, 0, 1], 
            "--",
            color=colors['RK2'],    
            label="RK2")
    ax.plot(rk4_vals[:, 0, 0], 
            rk4_vals[:, 0, 1], 
            ":",
            color=colors['RK4'],    
            label="RK4")

    ax.set_title("$h = {:.4f}$".format(h))
    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")
    ax.set_aspect("equal")
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig("Figs/orbit_dt_{:.4f}.pdf".format(h))

    if h == h_range[-1]:
        mosaic = ("AB\n"
                  "AC")
        fig, axs = plt.subplot_mosaic(mosaic)
        axs = list(axs.values())
        for i in [0,1,2]:
            axs[i].plot(W_ana[:, 0, 0], 
                       W_ana[:, 0, 1],
                       "-.", 
                       color=colors['Analytical'],
                       label="Analytical")
            axs[i].plot(eu_vals[:, 0, 0],  
                       eu_vals[:, 0, 1],  
                       "-",
                       color=colors['Euler'],  
                       label="Euler")
            axs[i].plot(rk2_vals[:, 0, 0],
                       rk2_vals[:, 0, 1], 
                       "--",
                       color=colors['RK2'],    
                       label="RK2")
            axs[i].plot(rk4_vals[:, 0, 0], 
                       rk4_vals[:, 0, 1], 
                       ":",
                       color=colors['RK4'],    
                       label="RK4")

            axs[i].set_aspect("equal")

        #fig.suptitle("$\\Var{{t}} = {:.4f}$".format(h))
        axs[0].set_xlabel("$x$")
        axs[0].set_ylabel("$y$")
        axs[0].legend(loc="upper left")

        win_1 = 0.02
        axs[1].set_xlim(0 - win_1, 0 + win_1)
        axs[1].set_ylim(1 - win_1, 1 + win_1)
        #axs[0].indicate_inset_zoom(axs[1], lw=1)
        win_2 = 1e-6
        axs[2].set_xlim(0 - win_2, 0 + win_2)
        axs[2].set_ylim(1 - win_2, 1 + win_2)
        #axs[1].indicate_inset_zoom(axs[2], lw=1)

        ln1 = ConnectionPatch(xyA=(0,1), xyB=(0-win_1,1+win_1), 
                              coordsA="data", coordsB="data",
                              axesA=axs[0], axesB=axs[1], 
                              color="k", lw=1, alpha=0.5)
        ln2 = ConnectionPatch(xyA=(0,1), xyB=(0-win_1,1-win_1), 
                              coordsA="data", coordsB="data",
                              axesA=axs[0], axesB=axs[1], 
                              color="k", lw=1, alpha=0.5)
        fig.add_artist(ln1)
        fig.add_artist(ln2)

        ln3 = ConnectionPatch(xyA=(0,1), xyB=(0-win_2,1+win_2), 
                              coordsA="data", coordsB="data",
                              axesA=axs[1], axesB=axs[2], 
                              color="k", lw=1, alpha=0.5)
        ln4 = ConnectionPatch(xyA=(0,1), xyB=(0+win_2,1+win_2), 
                              coordsA="data", coordsB="data",
                              axesA=axs[1], axesB=axs[2], 
                              color="k", lw=1, alpha=0.5)
        fig.add_artist(ln3)
        fig.add_artist(ln4)
        #fig.tight_layout()
        fig.savefig("Figs/orbit_dt.pdf")



# ---------------------------------------------------------------
# 3. Summary Plots: CPU time and final energy error (Log-Log)
# ---------------------------------------------------------------

# --- Step size vs. CPU Time (Log-Log) ---
fig, ax = plt.subplots()
ax.plot(h_range[:-1], time_euler[:-1], 'o-', color='C0', label="Euler")
ax.plot(h_range[:-1], time_rk2[:-1],   's--', color='C2', label="RK2")
ax.plot(h_range[:-1], time_rk4[:-1],   '^:', color='C3', label="RK4")

ax.set_xscale("log")
ax.set_yscale("log")

ax.set_xlabel("Step size $\\Var{{t}}$")
ax.set_ylabel("CPU Time $t_\\mathrm{CPU}\\axunit{{s}}$")
#ax.minorticks_on()
#ax.grid(True, which="major", linestyle="--", linewidth=0.5, alpha=0.7)
#ax.grid(True, which="minor", linestyle=":", linewidth=0.5, alpha=0.5)
ax.legend(loc="best")
fig.tight_layout()
fig.savefig("Figs/dt_vs_cpu_time_loglog.pdf")

# --- Step size vs. Final Energy Error (Log-Log) ---
fig, ax = plt.subplots()

ax.plot(h_range[:-1], err_euler[:-1], 'o-', color='C0', label="Euler")
ax.plot(h_range[:-1], err_rk2[:-1],   's--', color='C2', label="RK2")
ax.plot(h_range[:-1], err_rk4[:-1],   '^:', color='C3', label="RK4")

ax.set_xscale("log")
ax.set_yscale("log")
# Machine Epsilon line (horizontal)
ax.axhline(eps, color='darkred', ls='-.', 
           label='Machine precision $\\epsilon$')

ax.set_xlabel("Step size $\\Var{{t}}$")
ax.set_ylabel("$\\abs{{E_{\\mathrm{analytical}} - E_{\\mathrm{numerical}}}}$")
#ax.minorticks_on()
#ax.grid(True, which="major", linestyle="--", linewidth=0.5, alpha=0.7)
#ax.grid(True, which="minor", linestyle=":", linewidth=0.5, alpha=0.5)

# Ensure 'Machine Epsilon' is in legend
"""
handles, labels = ax.get_legend_handles_labels()
if 'Machine Epsilon' not in labels:
    import matplotlib.lines as mlines
    h_me = mlines.Line2D([], [], color='darkred', ls='--', label='Machine Epsilon')
    handles.append(h_me)
    labels.append('Machine Epsilon')
ax.legend(handles, labels, loc="best", fontsize=12)
"""
ax.legend()
fig.tight_layout()
fig.savefig("Figs/timestep_vs_final_energy_error_loglog1.pdf")
plt.show()


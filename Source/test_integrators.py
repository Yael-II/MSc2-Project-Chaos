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

import integrator as intg       #  integrator module
import kepler_orbit as kep      # Analytical & numeric Kepler functions
import initial_conditions as ic  # For initial condition

# ----------------------------
# 1. Setup & global parameters
# ----------------------------
t0      = 0.0
T_final = 8.0
y0      = ic.one_part(1, 0, 0, 1)  # [x, y, vx, vy]


h_range = np.logspace(-3.5, -0.1, 25)

# For plotting lines/colors
methods = [
    ("Euler", intg.euler, 'o--', 'red'),
    ("RK2",   intg.rk2,   's--', 'royalblue'),
    ("RK4",   intg.rk4,   '^--', 'limegreen')
]
colors = {'Analytical': 'navy', 'Euler': 'red', 'RK2': 'royalblue', 'RK4': 'limegreen'}

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
for h in h_range:
    N = int(T_final / h)

    # Analytical solution
    t_ana, pos_ana, vel_ana, en_ana = kep.kepler_analytical_orb(y0, h, N)
    E_analytical_final = en_ana[-1]

    # Numerical integrators + timing
    all_solutions = {}
    for (label, method, *_), store_err, store_t in zip(
        methods,
        [err_euler, err_rk2, err_rk4],
        [time_euler, time_rk2, time_rk4]
    ):
        start_time = time.time()
        t_num, y_num = intg.integrator_type(t0, y0, h, N, kep.kepler_orbnum, method)
        elapsed = time.time() - start_time

        store_t.append(elapsed)

        # Final energy error
        en_num = kep.kepler_enrgnum(y_num)
        E_numerical_final = en_num[-1]
        store_err.append(abs(E_analytical_final - E_numerical_final))

        all_solutions[label] = y_num

    # Orbit plot (optional, can comment out if too many figures)
    eu_vals  = all_solutions["Euler"]
    rk2_vals = all_solutions["RK2"]
    rk4_vals = all_solutions["RK4"]

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.plot(pos_ana[:, 0], pos_ana[:, 1], '-.', color=colors['Analytical'],
            label="Analytical", zorder=4)
    ax.plot(eu_vals[:, 0, 0],  eu_vals[:, 0, 1],  color=colors['Euler'],  label="Euler")
    ax.plot(rk2_vals[:, 0, 0], rk2_vals[:, 0, 1], color=colors['RK2'],    label="RK2")
    ax.plot(rk4_vals[:, 0, 0], rk4_vals[:, 0, 1], color=colors['RK4'],    label="RK4")
    ax.set_title(f"Orbit for dt = {h:.4g}")
    ax.set_xlabel("x position(reduced units)")
    ax.set_ylabel("y position(reduced units)")
    ax.legend(loc="best")
    ax.grid(True, which='both', ls='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig(f"orbit_dt_{h:.4g}.png", dpi=300, bbox_inches='tight')
    plt.show()

# ---------------------------------------------------------------
# 3. Summary Plots: CPU time and final energy error (Log-Log)
# ---------------------------------------------------------------

# --- Step size vs. CPU Time (Log-Log) ---
fig, ax = plt.subplots(figsize=(8, 6))
ax.loglog(h_range, time_euler, 'o--', color='red',       label="Euler")
ax.loglog(h_range, time_rk2,   's--', color='royalblue', label="RK2")
ax.loglog(h_range, time_rk4,   '^--', color='limegreen', label="RK4")

ax.set_xlabel(r"Step size $(\Delta t)$", fontsize=12)
ax.set_ylabel("CPU Time (s)", fontsize=12)
ax.set_title(r"$\Delta t$ vs. CPU Time ", fontsize=14)
ax.minorticks_on()
ax.grid(True, which="major", linestyle="--", linewidth=0.5, alpha=0.7)
ax.grid(True, which="minor", linestyle=":", linewidth=0.5, alpha=0.5)
ax.legend(loc="best", fontsize=12)
plt.tight_layout()
plt.savefig("dt_vs_cpu_time_loglog.png", dpi=300, bbox_inches='tight')
plt.show()

# --- Step size vs. Final Energy Error (Log-Log) ---
fig, ax = plt.subplots(figsize=(8, 6))
ax.loglog(h_range, err_euler, 'o--', color='red',       label="Euler")
ax.loglog(h_range, err_rk2,   's--', color='royalblue', label="RK2")
ax.loglog(h_range, err_rk4,   '^--', color='limegreen', label="RK4")

# Machine Epsilon line (horizontal)
ax.axhline(eps, color='darkred', ls='--', label='Machine Epsilon')

ax.set_xlabel(r"Step size $(\Delta t)$", fontsize=12)
ax.set_ylabel(r"$|E_{\mathrm{analytical}} - E_{\mathrm{numerical}}|$", fontsize=12)
ax.set_title(r"$\Delta t$ vs. Final Energy Error ", fontsize=14)
ax.minorticks_on()
ax.grid(True, which="major", linestyle="--", linewidth=0.5, alpha=0.7)
ax.grid(True, which="minor", linestyle=":", linewidth=0.5, alpha=0.5)

# Ensure 'Machine Epsilon' is in legend
handles, labels = ax.get_legend_handles_labels()
if 'Machine Epsilon' not in labels:
    import matplotlib.lines as mlines
    h_me = mlines.Line2D([], [], color='darkred', ls='--', label='Machine Epsilon')
    handles.append(h_me)
    labels.append('Machine Epsilon')
ax.legend(handles, labels, loc="best", fontsize=12)

plt.tight_layout()
plt.savefig("timestep_vs_final_energy_error_loglog1.png", dpi=300, bbox_inches='tight')
plt.show()


#!/usr/bin/env python
"""
Gottwald and Melbourne method

Study chaotic behaviour with the method of Gottwald & Melbourne (2003).

@ Author: Moussouni, Yaël (MSc student) & Bhat, Junaid Ramzan (MSc student)
@ Institution:  Université de Strasbourg, CNRS, Observatoire astronomique
                de Strasbourg, UMR 7550, F-67000 Strasbourg, France
@ Date: 2025-01-01

Licence:
Order and Chaos in a 2D potential
Copyright (C) 2025 Yaël Moussouni (yael.moussouni@etu.unistra.fr)
                   Bhat, Junaid Ramzan (junaid-ramzan.bhat@etu.unistra.fr)

gottwald_melbourne.py
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
from scipy.integrate import cumulative_trapezoid

import potentials as pot
import integrator as itg
import initial_conditions as init

# -----------------------------------------------------------------------------------
# Parameters
# -----------------------------------------------------------------------------------
DEFAULT_N_iter = 30000
DEFAULT_h = 0.01
DEFAULT_c = 1.7

E_all = np.array([1/100, 1/12, 1/10, 1/8, 1/6])

if "YII_1" in plt.style.available: plt.style.use("YII_1")

# -----------------------------------------------------------------------------------
# 1) Integrate Hénon–Heiles and Return x(t), etc.
# -----------------------------------------------------------------------------------
def compute_coordinates(E: float,
                        N_iter: int = DEFAULT_N_iter,
                        h: float = DEFAULT_h,
                        N_part: int = 1) -> tuple:
    """
    Integrate Hénon–Heiles for N_iter steps at step size h for a single
    random initial condition at energy E.
    Returns:
      t_part: array of times of length N_iter
      x_part, y_part, u_part, v_part: arrays of length N_iter each
    """
    # Generate 1 initial condition at energy E
    W_init = init.n_energy_part(pot.hh_potential, N_part, E)
    W0 = W_init[:, :, 0]  # take first particle

    # Perform integration using RK4
    final_t, sol_array = itg.rk4(0.0, W0, h, N_iter, pot.hh_evolution)
    # Reconstruct the time array using step size and number of iterations

    # Extract coordinate arrays
    x_part = sol_array[:, 0, 0]
    y_part = sol_array[:, 0, 1]
    u_part = sol_array[:, 1, 0]
    v_part = sol_array[:, 1, 1]

    return final_t, x_part, y_part, u_part, v_part

# -----------------------------------------------------------------------------------
# 2) Compute theta(t_i) in discrete form
# -----------------------------------------------------------------------------------
def compute_theta_discrete(t_array, x_array, c=DEFAULT_c):
    int_x = cumulative_trapezoid(x_array, t_array, initial=0.0)
    theta_array = c * t_array + int_x
    return theta_array

# -----------------------------------------------------------------------------------
# 3) Compute p(t_i) in discrete form
# -----------------------------------------------------------------------------------
def compute_p_discrete(t_array, x_array, theta_array):
    integrand = x_array * np.cos(theta_array)
    p_array = cumulative_trapezoid(integrand, t_array, initial=0.0)
    return p_array

# -----------------------------------------------------------------------------------
# Main: Compute p(t) for each energy and plot on a single graph
# -----------------------------------------------------------------------------------
if __name__ == "__main__":
    fig, ax = plt.subplots()  # single figure and axis

    # Loop over energies and plot each p(t) on the same axis
    i = 0
    line_styles = ["-", "-", "--", ":", ":"]
    for E in E_all:
        # 1) Integrate Hénon–Heiles for current energy
        t_part, x_part, y_part, u_part, v_part = compute_coordinates(E)

        # 2) Compute theta(t) for current energy
        theta_arr = compute_theta_discrete(t_part, x_part, c=DEFAULT_c)

        # 3) Compute p(t) for current energy
        p_arr = compute_p_discrete(t_part, x_part, theta_arr)

        # Plot p(t) for the current energy on the same axis
        ax.plot(t_part, p_arr, 
                line_styles[i], 
                color="C{}".format(i).replace("C1", "C5"),
                label="$E={:.3f}$".format(E))
        i += 1

    #ax.set_title("$p(t)$ (Melbourne–Gottwald Test)")
    ax.set_xlabel("$t$")
    ax.set_ylabel("$p(t)$")
    ax.legend()

    # Enable both major and minor gridlines
    #ax.grid(which='both', linestyle='--', linewidth=0.5)
    #ax.minorticks_on()

    fig.tight_layout()
    fig.savefig("Figs/gottwald_melbourne.pdf")
    plt.show()


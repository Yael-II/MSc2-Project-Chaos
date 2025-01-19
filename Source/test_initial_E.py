#!/usr/bin/env python
"""
Test: Initial Conditions With a Given Energy

Sample random particles with the same given energy in a valid coordinates range.

@ Author: Moussouni, Yaël (MSc student) & Bhat, Junaid Ramzan (MSc student)
@ Institution:  Université de Strasbourg, CNRS, Observatoire astronomique
                de Strasbourg, UMR 7550, F-67000 Strasbourg, France
@ Date: 2025-01-01

Licence:
Order and Chaos in a 2D potential
Copyright (C) 2025 Yaël Moussouni (yael.moussouni@etu.unistra.fr)
                   Bhat, Junaid Ramzan (junaid-ramzan.bhat@etu.unistra.fr)

test_initial_E.py
Copyright (C) 2025 Yaël Moussouni (yael.moussouni@etu.unistra.fr)
                   Bhat, Junaid Ramzan (junaid-ramzan.bhat@etu.unistra.fr)

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

import potentials as pot
import energies as ene
import integrator as itg
import initial_conditions as init
import poincare_sections as pcs

if "YII_1" in plt.style.available: plt.style.use("YII_1")
# Loading matplotlib style...
# parameters
N_grid = 1000
N_inter = 30000
N_part = 100
E = 1/12
h = 0.01

if __name__ == "__main__":
    # Initial conditions
    W_grid = init.mesh_grid(N_grid, xmin=-1, xmax=1, ymin=-1, ymax=1)
    X_grid = W_grid[0]
    Y_grid = W_grid[1]
    potential = pot.hh_potential(W_grid, position_only=True)
    pot_valid = np.ma.masked_where(potential > E, potential)

    W_all_part = init.n_energy_part(pot.hh_potential, N_part, E)

    # Plot
    fig, ax = plt.subplots(1)

    pcm = ax.pcolormesh(X_grid, Y_grid, pot_valid, vmin = 0)
    sct = ax.scatter(W_all_part[0, 0], W_all_part[0, 1], s=1, color="C3")
    fig.colorbar(pcm, label="potential")
    ax.set_title("$E = {:.2f}$".format(E))
    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")
    ax.set_aspect("equal")

    fig.savefig("Figs/initial_E.png")

    plt.show(block=True)



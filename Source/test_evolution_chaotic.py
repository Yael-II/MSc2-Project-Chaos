#!/usr/bin/env python
"""
Test: Evolution

Evolve a particle with a given energy in a potential and show the result path
followed by the particle in a given time.

@ Author: Moussouni, Yaël (MSc student) & Bhat, Junaid Ramzan (MSc student)
@ Institution:  Université de Strasbourg, CNRS, Observatoire astronomique
                de Strasbourg, UMR 7550, F-67000 Strasbourg, France
@ Date: 2025-01-01

Licence:
Order and Chaos in a 2D potential
Copyright (C) 2025 Yaël Moussouni (yael.moussouni@etu.unistra.fr)
                   Bhat, Junaid Ramzan (junaid-ramzan.bhat@etu.unistra.fr)

test_evolution.py
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

# Constants
N_grid = 5000
h = 0.01
N_inter = 5000
eps = 1e-2

x_0 = 0.3
y_0 = 0.1
E = 1/6

# Main
if __name__ == "__main__":
    W_part = init.one_part(x_0, y_0,0,0)
    POT = pot.hh_potential(W_part)[0]
    u_0 = np.sqrt(2 * (E - POT))
    v_0 = 0.0

    W_part[1,0] = u_0
    W_part[1,1] = v_0

    pos_t, positions = itg.rk4(0, W_part, h, N_inter, pot.hh_evolution)

    pos_x = positions[:,0,0]
    pos_y = positions[:,0,1]
    vel_x = positions[:,1,0]
    vel_y = positions[:,1,1]

    W_grid = init.mesh_grid(N_grid)
    X_grid = W_grid[0]
    Y_grid = W_grid[1]
    potential = pot.hh_potential(W_grid, position_only=True)

    fig, ax = plt.subplots(1)

    pcm = ax.pcolormesh(X_grid, Y_grid, potential)
    line = ax.plot(pos_x, pos_y, color="C3")
    fig.colorbar(pcm, label="potential")

    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")
    ax.set_aspect("equal")

    plt.savefig("Figs/evolution_chaotic.png")

    plt.show(block=True)



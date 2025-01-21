#!/usr/bin/env python
"""
Test: Potential

Draw the Kepler potential and the Hénon--Heils potential

@ Author: Moussouni, Yaël (MSc student) & Bhat, Junaid Ramzan (MSc student)
@ Institution:  Université de Strasbourg, CNRS, Observatoire astronomique
                de Strasbourg, UMR 7550, F-67000 Strasbourg, France
@ Date: 2025-01-01

Licence:
Order and Chaos in a 2D potential
Copyright (C) 2025 Yaël Moussouni (yael.moussouni@etu.unistra.fr)
                   Bhat, Junaid Ramzan (junaid-ramzan.bhat@etu.unistra.fr)

test_potentials.py
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
import initial_conditions as init

if "YII_1" in plt.style.available: plt.style.use("YII_1")
# Loading matplotlib style...

W = init.mesh_grid(300)

def kepler(W):
    """Plots the Kepler potential"""
    X = W[0]
    Y = W[1]
    POT = pot.kepler_potential(W, position_only=True)
    fig, ax = plt.subplots(1)
    ax.set_title("Kepler potential")
    pcm = ax.pcolormesh(X, Y, POT)
    fig.colorbar(pcm, label="potential")
    ax.set_aspect("equal")
    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")

    fig.savefig("Figs/pot_kepler.png")
    return 0

def hh(W):
    """Plots the Hénon--Heils potential"""
    X = W[0]
    Y = W[1]
    POT = pot.hh_potential(W, position_only=True)
    fig, ax = plt.subplots(1)
    ax.set_title("Hénon–Heils potential")
    pcm = ax.pcolormesh(X, Y, POT)
    fig.colorbar(pcm, label="potential")
    ax.set_aspect("equal")
    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")

    fig.savefig("Figs/pot_hh.png")
    return 0

if __name__ == "__main__":
    kepler(W)
    hh(W)

    plt.show(block=True)


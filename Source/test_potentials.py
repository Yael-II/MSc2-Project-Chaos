#!/usr/bin/env python
"""
Test: Potential

Draw the Kepler potential and the Hénon--Heils potential

@ Author: Moussouni, Yaël (MSc student) & Bhat, Junaid Ramzan (MSc student)
@ Institution:  Université de Strasbourg, CNRS, Observatoire astronomique
                de Strasbourg, UMR 7550, F-67000 Strasbourg, France
@ Date: 2024-11-29
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
    return 0

if __name__ == "__main__":
    kepler(W)
    hh(W)

    plt.show(block=True)


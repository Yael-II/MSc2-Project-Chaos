#!/usr/bin/env python
"""
Initial Conditions

Generate initial conditions depending on different criteria

@ Author: Moussouni, Yaël (MSc student) & Bhat, Junaid Ramzan (MSc student)
@ Institution:  Université de Strasbourg, CNRS, Observatoire astronomique
                de Strasbourg, UMR 7550, F-67000 Strasbourg, France
@ Date: 2024-11-29
"""

import numpy as np
import matplotlib.pyplot as plt
POS_MIN = -1
POS_MAX = +1
VEL_MIN = -1
VEL_MAX = +1
N_PART = 1

def mesh_grid(N: int = N_PART,
              xmin: float = POS_MIN, 
              xmax: float = POS_MAX, 
              ymin: float = POS_MIN, 
              ymax: float = POS_MAX) -> np.ndarray:
    """Generates a set of regularly sampled particles with no velocity
    @ params:
        - N: number of particles
        - xmin: minimum value for position x
        - xmax: maximum value for position x
        - ymin: minimum value for position y
        - ymax: maximum value for position y
    @ returns:
        - W: phase-space vector
    """
    X = np.linspace(xmin, xmax, N)
    Y = np.linspace(ymin, ymax, N)
    X,Y = np.meshgrid(X,Y, indexing="ij")
    return np.array([X,Y])

def one_part(x0: float = 0, 
             y0: float = 0, 
             u0: float= 0, 
             v0: float = 0) -> np.ndarray:
    """Generates a particle at position (x0, y0) 
    with velocities (u0,v0)
    @ params:
        - N: number of particles
        - x0: initial position x
        - y0: initial position y
        - u0: initial velocity u
        - v0: initial velocity v
    @ returns:
        - W: phase-space vector
    """
    X = x0
    Y = y0
    U = u0
    V = v0
    return np.array([[X,Y], [U,V]])

def n_energy_part(potential,
                  N: int = N_PART,
                  E: float = 0,
                  xmin: float = -1,
                  xmax: float = +1,
                  ymin: float = -0.5,
                  ymax: float = +1):
    """Generates N particles with an energy E in a potential.
    @ params:
        - potential: gravitational potential
        - N: number of particles
        - E: total energy
        - xmin: minimum value for position x
        - xmax: maximum value for position x
        - ymin: minimum value for position y
        - ymax: maximum value for position y
    @ returns:
        - W: an array of all the positions and velocities.
    """
    X = []
    Y = []
    POT = []
    U = []
    V = []
    while len(X) < N:
        x = np.random.random()*(xmax-xmin)+xmin
        y = np.random.random()*(ymax-ymin)+ymin
        w = np.array([x, y])
        pot = potential(w, position_only=True)[0]
        if pot <= E:
            X.append(x)
            Y.append(y)
            POT.append(pot)
    X = np.array(X)
    Y = np.array(Y)
    POT = np.array(POT)
    U = np.zeros_like(X)
    V = np.zeros_like(Y)
    C = np.sqrt(2 * (E - POT))
    THETA = np.random.random(N)*2*np.pi
    U = C*np.cos(THETA)
    V = C*np.sin(THETA)
    return np.array([[X, Y], [U, V]])

def n_energy_2part(potential,
                   N: int = N_PART,
                   E: float = 0,
                   sep: float = 1e-7,
                   xmin: float = -1,
                   xmax: float = +1,
                   ymin: float = -0.5,
                   ymax: float = +1):
    """Generate a sample of 2N particles with the energy E in a potential in 
    two sets: one "normal" set (see n_energy_part), and a slightly shifted set
    with a separation sep.
    @ params:
        - potential: gravitational potential
        - N: number of particles
        - E: total energy
        - sep: the separation between the two sets
        - xmin: minimum value for position x
        - xmax: maximum value for position x
        - ymin: minimum value for position y
        - ymax: maximum value for position y
    @ returns:
        - (W1, W2): the two arrays of all the positions and velocities for 
        each set.
    """
    W_1 = n_energy_part(potential, N, E)
    W_2 = np.zeros_like(W_1)
    alpha = np.random.uniform(0, 2*np.pi, N)
    W_2[0, 0] = W_1[0, 0] + sep*np.cos(alpha)
    W_2[0, 1] = W_1[0, 1] + sep*np.sin(alpha)
    W_2[1, 0] = W_1[1, 0]
    W_2[1, 1] = W_1[1, 1]
    return (W_1, W_2)






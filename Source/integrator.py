#!/usr/bin/env python
"""
Integrator

Integrate differential equations.

@ Author: Moussouni, Yaël (MSc student) & Bhat, Junaid Ramzan (MSc student)
@ Institution:  Université de Strasbourg, CNRS, Observatoire astronomique
                de Strasbourg, UMR 7550, F-67000 Strasbourg, France
@ Date: 2025-01-01

Licence:
Order and Chaos in a 2D potential
Copyright (C) 2025 Yaël Moussouni (yael.moussouni@etu.unistra.fr)
                   Bhat, Junaid Ramzan (junaid-ramzan.bhat@etu.unistra.fr)

integrator.py
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

def euler(t0: float, 
          W0: np.ndarray, 
          h: float, 
          n: int, 
          func):
    """Euler method adapted for state vector [[x, y], [u, v]]
    @ params
        - t0: initial time value
        - W0: initial state vector [[x, y], [u, v]]
        - h: step size (time step)
        - n: number of steps
        - func: RHS of differential equation
    @returns: 
        - t, W: time and state (solution) arrays
    """
    time = np.zeros(n)
    W = np.zeros((n,) + np.shape(W0))

    t = t0
    w = W0
    for i in range(n):
        k1 = func(t, w)
        w = w + h*k1
        t = t + h

        time[i] = t
        W[i] = w
    return time, W

def rk2(t0: float, 
        W0: np.ndarray, 
        h: float, 
        n: int, 
        func):
    """RK2 method adapted for state vector [[x, y], [u, v]]
    @ params
        - t0: initial time value
        - W0: initial state vector [[x, y], [u, v]]
        - h: step size (time step)
        - n: number of steps
        - func: RHS of differential equation
    @returns: 
        - t, W: time and state (solution) arrays
    """
    time = np.zeros(n)
    W = np.zeros((n,) + np.shape(W0))

    t = t0
    w = W0
    for i in range(n):
        k1 = func(t, w)
        k2 = func(t + h/2, w + h/2*k1)

        w = w + h*k2
        t = t + h

        time[i] = t
        W[i] = w
    return time, W

def rk4(t0: float, 
        W0: np.ndarray, 
        h: float, 
        n: int, 
        func):
    """RK4 method adapted for state vector [[x, y], [u, v]]
    @ params
        - t0: initial time
        - W0: initial state vector [[x, y], [u, v]]
        - h: step size (time step)
        - n: number of steps
        - func: RHS of differential equation
    @returns: 
        - t, W: time and state (solution) arrays
    """
    time = np.zeros(n)
    W = np.zeros((n,) + np.shape(W0))
    # to accommodate the state vector
    t = t0
    w = W0
    for i in range(n):
        k1 = func(t, w)
        k2 = func(t + h/2, w + h/2*k1)
        k3 = func(t + h/2, w + h/2*k2)
        k4 = func(t + h, w + h*k3)

        w = w + h*(k1/6 + k2/3 + k3/3 + k4/6)
        t = t + h

        time[i] = t
        W[i] = w
    return time, W

def integrator_type(t0, W0, h, n, func, integrator):
    return integrator(t0, W0, h, n, func)

def kepler_analytical(t0: float, 
                      W0: np.ndarray, 
                      h: float, 
                      n: int):
    """Computes the evolution from the Kepler potential derivative
    @ params
        - t0: initial time value
        - W0: initial state vector [[x, y], [u, v]]
        - h: step size (time step)
        - n: number of steps
    @returns: 
        - t, W: time and state (solution) arrays  
    """
    X0 = W0[0 ,0]
    Y0 = W0[0, 1]
    U0 = W0[1, 0]
    V0 = W0[1, 1]

    time = np.arange(t0, t0 + n*h, h)
    W = np.zeros((n,) + np.shape(W0))

    R0 = np.sqrt(X0**2 + Y0**2)
    Omega0 = np.sqrt(U0**2 + V0**2)/R0

    X = R0 * np.cos(Omega0 * time)
    Y = R0 * np.sin(Omega0 * time)
    U = -R0 * Omega0 * np.sin(Omega0 * time)
    V = R0 * Omega0 * np.cos(Omega0 * time)

    W = np.array([[X, Y], [U, V]])
    W = np.swapaxes(W, 0, 2)
    W = np.swapaxes(W, 1, 2)
    return time, W

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

def euler(x0, y0, h, n, func): # FIXME cannot be used with vectors
    """DEPRECIATED - DO NOT USE
    Euler method adapted for state vector [[x, y], [u, v]]
    :param x0: initial time value
    :param y0: initial state vector [[x, y], [u, v]]
    :param h: step size (time step)
    :param n: number of steps
    :param func: RHS of differential equation
    :returns: x array, solution array
    """
    x_values = np.zeros(n)
    y_values = np.zeros((n, 2, 2))  #  to accommodate the state vector

    for i in range(n):
        dydt = func(x0, y0)
        y0 = y0 + h * dydt
        x0 = x0 + h

        x_values[i] = x0
        y_values[i, :, :] = y0

    return x_values, y_values

# Updated RK2 integrator
def rk2(x0, y0, h, n, func): # FIXME cannot be used with vectors
    """ DEPRECIATED - DO NOT USE
    RK2 method adapted for state vector [[x, y], [u, v]]
    :param x0: initial time value
    :param y0: initial state vector [[x, y], [u, v]]
    :param h: step size (time step)
    :param n: number of steps
    :param func: RHS of differential equation
    :returns: x array, solution array
    """
    x_values = np.zeros(n)
    y_values = np.zeros((n, 2, 2))  #  to accommodate the state vector

    for i in range(n):
        k1 = func(x0, y0)
        k2 = func(x0 + h / 2., y0 + h / 2. * k1)

        y0 = y0 + h * (k1 / 2. + k2 / 2.)
        x0 = x0 + h

        x_values[i] = x0
        y_values[i, :, :] = y0

    return x_values, y_values

def rk4(t0: float, 
        W0: np.ndarray, 
        h: float, 
        n: int, 
        func):
    """RK4 method adapted for state vector [[x, y], [u, v]]
    @ params
        - x0: initial time value
        - y0: initial state vector [[x, y], [u, v]]
        - h: step size (time step)
        - n: number of steps
        - func: RHS of differential equation
    @returns: 
        - t, W: time and state (solution) arrays, 
    """
    time = np.zeros(n)
    W = np.zeros((n,) + np.shape(W0))
    # to accommodate the state vector
    t = t0
    w = W0
    for i in range(n):
        k1 = func(t, w)
        k2 = func(t + h / 2., w + h / 2. * k1)
        k3 = func(t + h / 2., w + h / 2. * k2)
        k4 = func(t + h, w + h * k3)

        w = w + h * (k1 / 6. + k2 / 3. + k3 / 3. + k4 / 6.)
        t = t + h

        time[i] = t
        W[i] = w

    return t, W


def integrator_type(x0, y0, h, n, func, int_type):
    """DEPRECIATED - DO NOT USE"""
    return int_type(x0, y0, h, n, func)

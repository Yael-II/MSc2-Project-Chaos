#!/usr/bin/env python
"""
Potentials

Functions of the different potentials (and their derivatives for the evolution)

@ Author: Moussouni, Yaël (MSc student) & Bhat, Junaid Ramzan (MSc student)
@ Institution:  Université de Strasbourg, CNRS, Observatoire astronomique
                de Strasbourg, UMR 7550, F-67000 Strasbourg, France
@ Date: 2025-01-01

Licence:
Order and Chaos in a 2D potential
Copyright (C) 2025 Yaël Moussouni (yael.moussouni@etu.unistra.fr)
                   Bhat, Junaid Ramzan (junaid-ramzan.bhat@etu.unistra.fr)

potentials.py
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

MAX_VAL = 1e3

def kepler_potential(W_grid: np.ndarray, 
                     position_only: bool = False) -> np.ndarray:
    """Computes the Kepler potential: V(R) = -G*m1*m2/R 
    (assuming G = 1, m1 = 1, m2 = 1)
    assuming the point mass at (x = 0, y = 0).
    @params:
        - W: Phase-space vector
        - position_only: True if W is np.array([X, Y])
    @returns:
        - computed potential
    """
    if position_only: 
        X = W_grid[0]
        Y = W_grid[1]
    else:
        X = W_grid[0,0]
        Y = W_grid[0,1]
    # If X or Y is not an array (or a list), but rather a scalar, then we
    # create a list of one element so that it can work either way
    if np.ndim(X) == 0: X = np.array([X])
    if np.ndim(Y) == 0: Y = np.array([Y])
    R = np.sqrt(X**2 + Y**2)
    return -1/R

def hh_potential(W_grid: np.ndarray, 
                 position_only=False) -> np.ndarray:
    """Computes the Hénon-Heiles potential.
    @params:
        - W: Phase-space vector
        - position_only: True if W is np.array([X, Y])
    @returns:
        - POT: Potential
    """
    if position_only:
        X = W_grid[0]
        Y = W_grid[1]
    else:
        X = W_grid[0, 0]
        Y = W_grid[0, 1]
        
    # If X or Y is not an array (or a list), but rather a scalar, then we
    # create a list of one element so that it can work either way
    if np.ndim(X) == 0: X = np.array([X])
    if np.ndim(Y) == 0: Y = np.array([Y])

    POT = (X**2 + Y**2 + 2*X**2*Y - 2*Y**3/3)/2
    return POT        

def hh_evolution(t: np.ndarray, W: np.ndarray):
    """Computes the evolution from the HH potential derivative
    @params
        - t: Time (not used)
        - W: Phase space vector
    &returns 
        - dot W: Time derivative of the phase space vector
    """
    X = W[0 ,0]
    Y = W[0, 1]
    U = W[1, 0]
    V = W[1, 1]
    DX = U 
    DY = V
    DU = -(2*X*Y + X)
    DV = -(X**2 - Y**2 + Y)

    return np.array([[DX, DY], [DU, DV]])

#!/usr/bin/env python
"""
Energies

Compute energies (kinetic or total)

@ Author: Moussouni, Yaël (MSc student) & Bhat, Junaid Ramzan (MSc student)
@ Institution:  Université de Strasbourg, CNRS, Observatoire astronomique
                de Strasbourg, UMR 7550, F-67000 Strasbourg, France
@ Date: 2025-01-01

Licence:
Order and Chaos in a 2D potential
Copyright (C) 2025 Yaël Moussouni (yael.moussouni@etu.unistra.fr)
                   Bhat, Junaid Ramzan (junaid-ramzan.bhat@etu.unistra.fr)

energies.py
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

def kinetic(W: np.ndarray) -> np.ndarray:
    """Computes the kinetic energy.
    @param 
        - W: Phase-space vectors
    @returns 
        - T: Kinetic energy
    """
    U = W[1,0]
    V = W[1,1]
    # If U or V is not an array (or a list), but rather a scalar, then we
    # create a list of one element so that it can work either way
    if np.ndim(U) == 0: U = np.array([U])
    if np.ndim(V) == 0: V = np.array([V])
    
    return (U**2 + V**2)/2

def total(W: np.ndarray,
          potential,
          kinetic = kinetic) -> np.ndarray:
    return potential(W) + kinetic(W)

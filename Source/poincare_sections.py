#!/usr/bin/env python
"""
Poincaré Sections

Computes the Poincaré Sections

@ Author: Moussouni, Yaël (MSc student) & Bhat, Junaid Ramzan (MSc student)
@ Institution:  Université de Strasbourg, CNRS, Observatoire astronomique
                de Strasbourg, UMR 7550, F-67000 Strasbourg, France
@ Date: 2025-01-01

Licence:
Order and Chaos in a 2D potential
Copyright (C) 2025 Yaël Moussouni (yael.moussouni@etu.unistra.fr)
                   Bhat, Junaid Ramzan (junaid-ramzan.bhat@etu.unistra.fr)

poincare_sections.py
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
def pcs_find(pos_x, pos_y, vel_x, vel_y):
    """Find Poincaré sections (PCS; x = 0)
    @ params:
        - pos_x: position along the x axis
        - pos_y: position along the y axis
        - pos_x: velocity along the x axis
        - pos_y: velocity along the y axis
    @ returns: (tuple)
        - pcs_pos_y: position of the points in the PCS along the y axis
        - pcs_vel_y: velocity of the points in the PCS along the y axis
    """
    if np.ndim(pos_x) == 1: 
        pos_x = np.array([pos_x])
        pos_y = np.array([pos_y])
        vel_x = np.array([vel_x])
        vel_y = np.array([vel_y])
    pcs_pos_y = []
    pcs_vel_y = []
    for j in range(len(pos_x[0])): # for each particle
        i = 0
        x = pos_x[:,j]
        y = pos_y[:,j]
        u = vel_x[:,j]
        v = vel_y[:,j]
        while i < len(x) - 1:
            if x[i] * x[i+1] < 0:
                y0 = y[i] \
                        + (y[i+1] - y[i])/(x[i+1] - x[i]) \
                        * (0 - x[i])
                v0 = v[i] \
                        + (v[i+1] - v[i])/(x[i+1] - x[i])\
                        * (0 - x[i])
                pcs_pos_y.append(y0)
                pcs_vel_y.append(v0)
            i += 1
    return pcs_pos_y, pcs_vel_y

def pcs_find_legacy(pos_x, pos_y, vel_x, vel_y):
    """DEPRECIATED - DO NOT USE
    Depreciated legacy function that should not be used
    """
    if np.ndim(pos_x) == 1: 
        pos_x = np.array([pos_x])
        pos_y = np.array([pos_y])
        vel_x = np.array([vel_x])
        vel_y = np.array([vel_y])
    pcs_pos_y = []
    pcs_vel_y = []
    for j in range(len(pos_x)): # for each particle
        i = 0
        x = pos_x[j]
        y = pos_y[j]
        u = vel_x[j]
        v = vel_y[j]
        while i < len(x) - 1:
            if x[i] * x[i+1] < 0:
                y0 = y[i] \
                        + (y[i+1] - y[i])/(x[i+1] - x[i]) \
                        * (0 - x[i])
                v0 = v[i] \
                        + (v[i+1] - v[i])/(x[i+1] - x[i])\
                        * (0 - x[i])
                pcs_pos_y.append(y0)
                pcs_vel_y.append(v0)
            i += 1
    return pcs_pos_y, pcs_vel_y

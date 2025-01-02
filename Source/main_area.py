#!/usr/bin/env python
"""
Main: Compute Relative Area

Computes the relative area covered bu the curves for different energies, to 
study ordered and chaotic regimes.

@ Author: Moussouni, Yaël (MSc student) & Bhat, Junaid Ramzan (MSc student)
@ Institution:  Université de Strasbourg, CNRS, Observatoire astronomique
                de Strasbourg, UMR 7550, F-67000 Strasbourg, France
@ Date: 2024-12-13
"""
import numpy as np
from scipy.optimize import curve_fit

import potentials as pot
import energies as ene
import integrator as itg
import initial_conditions as init
import poincare_sections as pcs

OUT_DIR = "./Output/"
FILENAME_PREFIX = "phase_separation_"
EXTENSION = ".csv"
DEFAULT_N_iter = int(1e5)
DEFAULT_N_part = 200
DEFAULT_h = 0.005
E_all = np.linspace(1/100, 1/6, 20)

def compute_mu(E: float,
               N_iter: int = DEFAULT_N_iter,
               N_part: int = DEFAULT_N_part,
               h: float = DEFAULT_h) -> tuple:
    """
    Computes the phase-space squared distances for particles of given energy E.
    @params:
        - E: the total energy of each particles
        - N_iter: the number of iteration
        - N_part: the number of particles
        - h: integration steps
    @returns:
        - mu: phase-space squared distance
    """
    W_1, W_2 = init.n_energy_2part(pot.hh_potential, N_part, E)
    t_1, positions_1 = itg.rk4(0, W_1, h, N_iter, pot.hh_evolution)
    x_1 = positions_1[:, 0, 0]
    y_1 = positions_1[:, 0, 1]
    u_1 = positions_1[:, 1, 0]
    v_1 = positions_1[:, 1, 1]
    
    t_2, positions_2 = itg.rk4(0, W_2, h, N_iter, pot.hh_evolution)
    x_2 = positions_2[:, 0, 0]
    y_2 = positions_2[:, 0, 1]
    u_2 = positions_2[:, 1, 0]
    v_2 = positions_2[:, 1, 1]
    dist_sq = (x_2[-25:] - x_1[-25:])**2 \
            + (y_2[-25:] - y_1[-25:])**2 \
            + (u_2[-25:] - u_1[-25:])**2 \
            + (v_2[-25:] - v_1[-25:])**2
    
    mu = np.sum(dist_sq, axis=0)
    return mu

if __name__ == "__main__":
    mu_all = []
    for i in range(len(E_all)):
        mu = compute_mu(E_all[i])
        filename = OUT_DIR + FILENAME_PREFIX\
                 + str(i) + EXTENSION
        np.savetxt(filename, mu)

#!/usr/bin/env python
"""
Main: Computes Poincaré Sections (Linear Algorithm)

Computes the Poincaré Sections with a linear algorithm 
(i.e. no parallel computing).

@ Author: Moussouni, Yaël (MSc student) & Bhat, Junaid Ramzan (MSc student)
@ Institution:  Université de Strasbourg, CNRS, Observatoire astronomique
                de Strasbourg, UMR 7550, F-67000 Strasbourg, France
@ Date: 2024-11-29
"""
import numpy as np

import potentials as pot
import integrator as itg
import initial_conditions as init
import poincare_sections as pcs

# Parameters
OUT_DIR = "./Output/"
FILENAME_PREFIX = "poincare_sections_linear_"
EXTENSION = ".csv"
DEFAULT_N_iter = 30000
DEFAULT_N_part = 100
DEFAULT_h = 0.01
E_all = np.array([1/100, 1/12, 1/10, 1/8, 1/6])

text_E = ["1/100", "1/12", "1/10", "1/8", "1/6"]

def compute_poincare_sections_linear(E: float,
                                    N_iter: int = DEFAULT_N_iter,
                                    N_part: int = DEFAULT_N_part,
                                    h: float = DEFAULT_h) -> tuple:
    """
    Computes the Poincaré sections for a given energy E.
    @params:
        - E: the total energy of each particles
        - N_iter: the number of iteration
        - N_part: the number of particles
        - h: integration steps
    @returns:
        - y_section, v_section: arrays containing the y and v coordinates of 
          the Poincaré sections
    """
    W_all_part = init.n_energy_part(pot.hh_potential, N_part, E)
    y_section = []
    v_section = []
    for i in range(N_part):
        W_part = W_all_part[:,:,i]
        
        # Perform integration
        t_part, coord_part = itg.rk4(0, W_part, h, N_inter, pot.hh_evolution)

        # Extract positions and velocities
        x_part = coord_part[:, 0, 0]
        y_part = coord_part[:, 0, 1]
        u_part = coord_part[:, 1, 0]
        v_part = coord_part[:, 1, 1]

        # Find Poincaré section points for the current initial condition
        y_pcs, v_pcs = pcs.pcs_find_legacy(x_part, y_part, u_part, v_part)
        # The legacy is important here, the algorithm is the same but the 
        # data format is different...

        # Append the current Poincaré section points to the overall lists
        y_section += y_pcs
        v_section += v_pcs
    return y_section, v_section

if __name__ == "__main__":
    y_section_all = []
    v_section_all = []
    for i in range(len(E_all)):
        y_section, v_section = compute_poincare_sections_linear(E_all[i])
        section = np.array([y_section, v_section])
        filename = OUT_DIR + FILENAME_PREFIX\
                + str(text_E[i][2:]) + EXTENSION
        np.savetxt(filename, section)





    


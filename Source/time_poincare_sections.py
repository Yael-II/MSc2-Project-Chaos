#!/usr/bin/env python
"""
Time: Poincaré Sections (Linear and Parallel)

Computes the running time between the linear and parallel algorithms.

@ Author: Moussouni, Yaël (MSc student) & Bhat, Junaid Ramzan (MSc student)
@ Institution:  Université de Strasbourg, CNRS, Observatoire astronomique
                de Strasbourg, UMR 7550, F-67000 Strasbourg, France
@ Date: 2024-11-29
"""

import time
import numpy as np
import main_poincare_sections_linear as lin
import main_poincare_sections_parallel as par

E_all = np.array([1/100, 1/12, 1/10, 1/8, 1/6])

par_time = []
lin_time = []

print("\033[34m" + "Please wait..." + "\033[0m")
for E in E_all:
    t_0 = time.time()
    par.compute_poincare_sections_numpy(E)
    t_1 = time.time()
    par_time.append(t_1-t_0)

print("\033[34m" + "Still wait..." + "\033[0m")
for E in E_all:
    t_0 = time.time()
    lin.compute_poincare_sections_linear(E)
    t_1 = time.time()
    lin_time.append(t_1-t_0)

print("\033[34m" + "Done!" + "\033[0m")
print("\033[36m" + "=== [ RESULTS ] ===" + "\033[0m")
print("\033[36m" 
      + "- Linear algorithm:   "
      + "\033[0m"
      + "({:07.4f} ± {:.4f}) s".format(np.mean(lin_time), np.std(lin_time))
      + "\033[36m"
      + " per energy iteration"
      + "\033[0m")
print("\033[36m" 
      + "- Parallel algorithm: "
      + "\033[0m"
      + "({:07.4f} ± {:.4f}) s".format(np.mean(par_time), np.std(lin_time))
      + "\033[36m"
      + " per energy iteration"
      + "\033[0m")

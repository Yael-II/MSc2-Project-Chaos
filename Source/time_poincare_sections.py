#!/usr/bin/env python
"""
Time: Poincaré Sections (Linear and Parallel)

Computes the running time between the linear and parallel algorithms.

@ Author: Moussouni, Yaël (MSc student) & Bhat, Junaid Ramzan (MSc student)
@ Institution:  Université de Strasbourg, CNRS, Observatoire astronomique
                de Strasbourg, UMR 7550, F-67000 Strasbourg, France
@ Date: 2025-01-01

Licence:
Order and Chaos in a 2D potential
Copyright (C) 2025 Yaël Moussouni (yael.moussouni@etu.unistra.fr)
                   Bhat, Junaid Ramzan (junaid-ramzan.bhat@etu.unistra.fr)

time_poincare_sections.py
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

#!/usr/bin/env python
"""
Plot: Poincaré Sections (Linear and Parallel)

Plots the Poincaré sections for different energies, computed either with linear
or parallel algorithms.

@ Author: Moussouni, Yaël (MSc student) & Bhat, Junaid Ramzan (MSc student)
@ Institution:  Université de Strasbourg, CNRS, Observatoire astronomique
                de Strasbourg, UMR 7550, F-67000 Strasbourg, France
@ Date: 2025-01-01

Licence:
Order and Chaos in a 2D potential
Copyright (C) 2025 Yaël Moussouni (yael.moussouni@etu.unistra.fr)
                   Bhat, Junaid Ramzan (junaid-ramzan.bhat@etu.unistra.fr)

plot_poincare_sections.py
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
import os
import numpy as np
import matplotlib.pyplot as plt

if "YII_1" in plt.style.available: plt.style.use("YII_1")

OUT_DIR = "./Output/"
FILENAME_PREFIX = "poincare_sections_"
EXTENSION = ".csv"

def plot_poincare_sections(filelist: list, title:str = "") -> int:
    """
    Plot all the Poincaré sections in the file list.
    @params:
        - filelist: the list of files in the output directory, with the format 
        "poincare_sections_{linear, parallel}_[1/E].csv"
        - title: title of the figure
    @returns: 
        - 0.
    """
    orderlist = np.argsort([(int(file
                                 .replace(FILENAME_PREFIX, "")
                                 .replace(EXTENSION, "")
                                 .replace("linear_", "")
                                 .replace("parallel_", ""))) 
                            for file in filelist])
    filelist = np.array(filelist)[orderlist]
    N = len(filelist)
    fig, axs = plt.subplot_mosaic("ABC\nDEF")
    axs = list(axs.values())
    fig.suptitle(title)
    for i in range(N):
        ax = axs[i]
        filename = filelist[i]
        inv_E = (filename
                 .replace(FILENAME_PREFIX, "")
                 .replace(EXTENSION, "")
                 .replace("linear_", "")
                 .replace("parallel_", ""))
        data = np.loadtxt(OUT_DIR + filename)
        y_section = data[0]
        v_section = data[1]
        ax.scatter(y_section, v_section, 
                   s=.1, color="C3", marker=",", alpha=0.5,
                   label="$E = 1/{}$".format(inv_E))
        ax.set_xlabel("$y$")
        ax.set_ylabel("$v$")
        ax.legend(loc="upper right")
    while i < N:
        i += 1
        axs[i].axis('off')
    if "linear" in title: kind = "linear"
    elif "parallel" in title: kind = "parallel"
    else: kind = "error"
    fig.savefig("Figs/pcs_{}.pdf".format(kind))
    return 0
print("\033[32m" 
      + "[P]arallel or [L]inear algorithm result, or [B]oth?" 
      + "\033[0m")
answer = input("\033[32m" + "> " + "\033[0m").upper()

if answer == "P":
    FILENAME_PREFIX += "parallel_"
elif answer == "L":
    FILENAME_PREFIX += "linear_"

filelist = [fname for fname in os.listdir(OUT_DIR) if FILENAME_PREFIX in fname]

if answer in ["L", "B"]:
    filelist_linear = [fname for fname in filelist if "linear_" in fname]
    plot_poincare_sections(filelist_linear, 
                           title=("Poincaré Sections "
                                  "(results from the linear algorithm)"))
if answer in ["P", "B"]:
    filelist_parallel = [fname for fname in filelist if "parallel_" in fname]
    plot_poincare_sections(filelist_parallel, 
                           title=("Poincaré Sections "
                                  "(results from the parallel algorithm)"))

plt.show()

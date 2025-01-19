#!/usr/bin/env python
"""
Plot: Area

Plots areas

@ Author: Moussouni, Yaël (MSc student) & Bhat, Junaid Ramzan (MSc student)
@ Institution:  Université de Strasbourg, CNRS, Observatoire astronomique
                de Strasbourg, UMR 7550, F-67000 Strasbourg, France
@ Date: 2025-01-01

Licence:
Order and Chaos in a 2D potential
Copyright (C) 2025 Yaël Moussouni (yael.moussouni@etu.unistra.fr)
                   Bhat, Junaid Ramzan (junaid-ramzan.bhat@etu.unistra.fr)

plot_area.py
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
FILENAME_PREFIX = "phase_separation_"
EXTENSION = ".csv"

def plot_area(filelist: list, mu_c = 1e-4) -> int:
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
                                 .replace(EXTENSION, "")))
                            for file in filelist])
    filelist = np.array(filelist)[orderlist]
    N = len(filelist)
    E =  np.linspace(1/100, 1/6, N)
    mu = []
    for filename in filelist:
        with open(OUT_DIR + filename) as file:
            data = file.readlines()
            data = [np.float64(d.replace("\n", "")) for d in data]
            mu.append(data)
    mu = np.array(mu)

    fig, ax = plt.subplots(1)
    ax.scatter([], [], s=1, color="k", label="Data")
    for i in range(len(mu)):
        Y = mu[i]
        ax.scatter([E[i]]*len(Y), Y, s=1, color="k", alpha=0.1)
    ax.scatter(E, np.mean(mu, axis=1), s=5, 
               color="C1", marker="o", label="Mean")
    ax.scatter(E, np.median(mu, axis=1), s=5, 
               color="C3", marker="s", label="Median")
    ax.plot(E, [mu_c]*len(E), 
            color="C5", label="Critical value $\\mu_\\mathrm{{c}}$")
    ax.text(0.01, 1e-5, "Regular", va="bottom", ha="left", color="C5")
    ax.text(0.01, 1e-3, "Chaotic", va="top", ha="left", color="C5")
    ax.set_xlabel("Energy $E$")
    ax.set_ylabel("Phase-space squared distance $\\mu$")
    ax.set_yscale("log")
    ax.legend()
    fig.savefig("Figs/mu.pdf")

    fig, ax = plt.subplots(1)
    N_reg = np.count_nonzero(mu < mu_c, axis=1)
    N = np.shape(mu)[1]
    Area = N_reg / N
    ax.scatter(E, Area, s=5, color="C0")
    ax.set_xlabel("Energy $E$")
    ax.set_ylabel("Area $N_\\mathrm{{reg}}/N_\\mathrm{{part}}$")
    fig.savefig("Figs/area.pdf")
    return 0

filelist = [f for f in os.listdir(OUT_DIR) if FILENAME_PREFIX in f]
plot_area(filelist)
plt.show()

# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(1, 720)

def f(x):
    return np.exp(-x) 
def g(x):
    return 0#*2 + np.sin(0.05*x) + np.random.normal(scale=0.125, size=len(x))





Vorticity = np.loadtxt('/Users/Francesco/Sci_Libs/dia_Vorticity_RMS_fct.txt', skiprows=1)
Energy = np.loadtxt('/Users/Francesco/Sci_Libs/dia_Energy_RMS_fct.txt', skiprows=1)



maxTime = 720
varLeft = Vorticity
varRight = Energy




plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Helvetica"],
    "figure.figsize" : (12/2.54, 8/2.54),
    "figure.dpi": 300
    })

title_string = r"$\mathrm{RMS}\left[f^{\prime}\right](t) = " + \
    "\sqrt{V^{-1} \int_{_{V}} ( f^{\prime} ) ^{2}\,\mathrm{d}V}$"




time = np.linspace(0, maxTime, maxTime)

fig, ax1 = plt.subplots()
# Solid
ax1.plot(time, varLeft[0, :], '-', label = r"R27d",linewidth=0.75)
ax1.plot(time, varLeft[1, :], '-', label = r"R3d",linewidth=0.75)
ax1.plot(time, varLeft[2, :], '-', label = r"R3LU",linewidth=0.75)
ax1.set_ylabel("Vorticity")
ax1.set_xlim([0, maxTime])
ax1.set_ylim([0, 4*10**-6])



ax2 = ax1.twinx()
ax2.plot(time, varRight[0, :], ':', label = r"R27d",linewidth=0.75)
ax2.plot(time, varRight[1, :], ':', label = r"R3d",linewidth=0.75)
ax2.plot(time, varRight[2, :], ':', label = r"R3LU",linewidth=0.75)
ax2.set_ylabel("Energy")
ax2.set_ylim([0, 0.25])

plt.title(title_string, fontsize=10)

# plt.savefig("Comparison_RMS.pdf", dpi=fig.dpi, bbox_inches='tight')
# plt.savefig("Comparison_RMS.eps", dpi=fig.dpi, bbox_inches='tight')
plt.show()
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 15:30:52 2017

@author: iniguez
"""
import numpy as np
import os
import DataAnalysisInterSteps as pc2
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib
from matplotlib.colors import from_levels_and_colors

totSimulations = 0
internalHinges = 12
flags = np.empty((6, 6)) # 12 internal hinges, 6 possible flags
flags[:] = np.NaN
ststs = np.empty((6))
ststs[:] = np.NaN
folder_name = "Results/triangular prism/sqp/energy/"

for stepedge in np.arange(6):
    kedge = 5*10**(stepedge/2)
    subfolder_name = "kedge%2.4f/" %(kedge)
    if os.path.exists(folder_name+subfolder_name):
        totSimulations += 1
        flags[stepedge,:], ststs[stepedge] = pc2.ReadandAnalizeFile(folder_name+subfolder_name, plot = False)
        print(folder_name)
            
allflags = np.sum(flags, axis = 1)

#%%
x = 5*10**(np.arange(6)/2)
x2 = 5*10**(np.arange(7)/2)

fig1, axes1 = plt.subplots(nrows=2, ncols=2, figsize=(pc2.cm2inch(35), pc2.cm2inch(20)))

for ax in axes1.flat:
    ax.set_xscale('log')
    pc2.NiceGraph2D(ax, 'kedge', 'Non-Convergence Rate', mincoord = [np.min(x)*10**-0.25, 0], maxcoord = [np.max(x)*10**0.25, 1], divisions = [x,np.nan])
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

axes1[0,0].bar(x, allflags, width=np.diff(x2)*0.3, ec='0.2', color = '#FE9128',align="center")
axes1[0,1].bar(x, flags[:,1], width=np.diff(x2)*0.3, ec='0.2', color = '#FE9128',align="center")
axes1[1,0].bar(x, flags[:,3], width=np.diff(x2)*0.3, ec='0.2', color = '#FE9128',align="center")
axes1[1,1].bar(x, flags[:,5], width=np.diff(x2)*0.3, ec='0.2', color = '#FE9128',align="center")

     
     
fig2, axes2 = plt.subplots(nrows=1, ncols=1, figsize=(pc2.cm2inch(35), pc2.cm2inch(20)))

axes2.set_xscale('log')
pc2.NiceGraph2D(axes2, 'kedge', 'Stable States', mincoord = [np.min(x)*10**-0.25, 0], maxcoord = [np.max(x)*10**0.25, np.max(ststs)+2], divisions = [x,np.nan])
axes2.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

axes2.bar(x, ststs, width=np.diff(x2)*0.3, ec='0.2', color = '#36648B',align="center")

fig1.show()
fig2.show()
fig1.savefig(folder_name + 'ConvergenceKedge.png', transparent = True)
fig2.savefig(folder_name + 'StableStatesKedge.png', transparent = True)
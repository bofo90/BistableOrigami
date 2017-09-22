# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 15:30:52 2017

@author: iniguez
"""
import numpy as np
import os
import plot_csv2 as pc2
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib
from matplotlib.colors import from_levels_and_colors

totSimulations = 0
internalHinges = 9
flags = np.empty((7, 7, 6)) # 12 internal hinges, 7 possible flags
flags[:] = np.NaN
ststs = np.empty((7,7))
ststs[:] = np.NaN
folder_name = "Results/cube/sqp/energy/"

for stephinge in np.arange(7):
    for stepangle in np.arange(7):
        khinge = 0.0005*10**(stephinge)
        kangle = 0.5*10**(stepangle/2)
        subfolder_name = "kangle%2.4f_khinge%2.4f/" %(kangle, khinge)
        if os.path.exists(folder_name+subfolder_name):
            totSimulations += 1
            flags[stephinge,stepangle,:], ststs[stephinge,stepangle] = pc2.ReadandAnalizeFile(folder_name+subfolder_name, plot = False)
            print(folder_name)
            
allflags = np.sum(flags, axis = 2)

#%%
x = 0.0005*10**np.arange(7)
y = 0.5*10**(np.arange(7)/2)
x2 = 0.0005*10**(np.arange(8)-0.5)
y2 = 0.5*10**(np.arange(8)/2-0.25)

fig1, axes1 = plt.subplots(nrows=2, ncols=2, figsize=(pc2.cm2inch(35), pc2.cm2inch(20)))
fig1.subplots_adjust(wspace = 0.23, left = 0.08)

for ax in axes1.flat:
    ax.set_xscale('log')
    ax.set_yscale('log')
    pc2.NiceGraph2D(ax, 'khinge', 'kangle', mincoord = [0.0005*10**(-0.5), 0.5*10**(-0.25)], maxcoord = [500*10**(0.5), 500*10**(0.25)], divisions = [x,y])
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.0E'))

cmap, norm = from_levels_and_colors(np.linspace(0,1,11), cm.summer_r(np.linspace(0, 1, 10)))
allflags = np.ma.masked_where(np.isnan(allflags), allflags)
flags = np.ma.masked_where(np.isnan(flags), flags)

cs1 = axes1[0,0].pcolormesh(x2, y2, allflags.T, cmap=cmap, norm = norm)
axes1[0,1].pcolormesh(x2, y2, flags[:,:,1].T, cmap=cmap, norm = norm)
axes1[1,0].pcolormesh(x2, y2, flags[:,:,3].T, cmap=cmap, norm = norm)
axes1[1,1].pcolormesh(x2, y2, flags[:,:,5].T, cmap=cmap, norm = norm)

for ax in axes1.flat:
    ax.grid(which='major', alpha=0.2, color = '0.2') 
    
cbar = plt.colorbar(cs1, format="%.2f", ax =axes1.ravel().tolist(), fraction=0.046, pad=0.03)
cbar.set_label('Non-Convergence rate', fontsize = 15, color = '0.2')
cbar.ax.tick_params(axis='y',colors='0.2')
cbar.ax.tick_params(axis='x',colors='0.2')
cbar.outline.set_edgecolor('0.2')

fig1.show()
fig1.savefig(folder_name + 'Convergence.png', transparent = True)

#%%

fig2, axes2 = plt.subplots(nrows=1, ncols=1, figsize=(pc2.cm2inch(35), pc2.cm2inch(20)))

axes2.set_xscale('log')
axes2.set_yscale('log')
pc2.NiceGraph2D(axes2, 'khinge', 'kangle', mincoord = [0.0005*10**(-0.5), 0.5*10**(-0.25)], maxcoord = [500*10**(0.5), 500*10**(0.25)], divisions = [x,y])
axes2.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
axes2.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

cmap2, norm2 = from_levels_and_colors(np.arange(np.nanmax(ststs)+1)+1, cm.Set2(np.linspace(0, 1, np.nanmax(ststs))))
ststs = np.ma.masked_where(np.isnan(ststs), ststs)
cs2 = axes2.pcolormesh(x2, y2, ststs.T, cmap=cmap2, norm = norm2)

axes2.grid(which='major', alpha=0.2, color = '0.2') 

cbar2 = plt.colorbar(cs2, format="%d", fraction=0.046, pad=0.03)
cbar2.set_label('Number of Stable States', fontsize = 15, color = '0.2')
cbar2.ax.tick_params(axis='y',colors='0.2')
cbar2.ax.tick_params(axis='x',colors='0.2')
cbar2.outline.set_edgecolor('0.2')

fig2.show()
fig2.savefig(folder_name + 'StableStates.png', transparent = True)
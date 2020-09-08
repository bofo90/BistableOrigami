# -*- coding: utf-8 -*-
"""
Created on Mon May 11 17:22:29 2020

@author: iniguez
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as matl
matl.rcParams['pdf.fonttype'] = 42
matl.rcParams['ps.fonttype'] = 42
matl.rcParams['font.family'] = 'sans-serif'
matl.rcParams['font.sans-serif'] = 'Arial'
matl.rcParams['mathtext.fontset'] = 'cm'

def cm2inch(value):
    return value/2.54

def NiceGraph2D(axes, nameX, nameY, mincoord = [np.NaN, np.NaN], maxcoord = [np.NaN, np.NaN], divisions = [np.NaN, np.NaN],buffer = [0.0, 0.0, 0.0]):
    
    gray = '0.2'
    matl.rcParams.update({'font.size': 9})

    if ~np.isnan(mincoord[0]) and ~np.isnan(maxcoord[0]):
        axes.set_xlim([mincoord[0]-buffer[0], maxcoord[0]+buffer[0]])
        if isinstance(divisions[0], (list, tuple, np.ndarray)):
            if ~np.isnan(divisions[0]).any():
                axes.set_xticks(divisions[0])
        else:
            if ~np.isnan(divisions[0]):
                axes.set_xticks(np.linspace(mincoord[0],maxcoord[0],divisions[0]))
    axes.set_xlabel(nameX,labelpad=0, color = gray)
    
    if ~np.isnan(mincoord[1]) and ~np.isnan(maxcoord[1]):
        axes.set_ylim([mincoord[1]-buffer[1], maxcoord[1]+buffer[1]])
        if isinstance(divisions[1], (list, tuple, np.ndarray)):
            if ~np.isnan(divisions[1]).any():
                axes.set_yticks(divisions[1])
        else:
            if ~np.isnan(divisions[1]):
                axes.set_yticks(np.linspace(mincoord[1],maxcoord[1],divisions[1]))
    axes.set_ylabel(nameY,labelpad=0, color = gray)
   
    axes.xaxis.label.set_color(gray)
    axes.tick_params(axis='x', colors=gray, direction = 'in', width = 0.4)
    axes.yaxis.label.set_color(gray)
    axes.tick_params(axis='y', colors=gray, direction = 'in', width = 0.4)
    axes.tick_params(pad = 2)
    
    axes.tick_params(axis='y', which='minor', colors=gray, direction = 'in', width = 0.4)
    axes.tick_params(axis='x', which='minor', colors=gray, direction = 'in', width = 0.4)
    
    for axis in ['top','bottom','left','right']:
        axes.spines[axis].set_linewidth(0.4)
        axes.spines[axis].set_color(gray)
        
    return

plt.close('all')  
shift = 0.08
das = 0.25
mo = 1-das
#%%
fig1 = plt.figure(figsize=(cm2inch(3.0), cm2inch(2.8)))
ax1 = plt.subplot(111)
fig1.subplots_adjust(top=0.987,
bottom=0.275,
left=0.28,
right=0.985)

NiceGraph2D(ax1, r'$\theta_1/\pi$', r'$\theta_2/\pi$', mincoord = [-1,-1], maxcoord = [1,1], 
            divisions = [3,3], buffer = [0.1, 0.1])

#Branch 1
ax1.plot([0,0],[-1,1], color = '#B3B3B3', linewidth = 1)
#Sub branches 1a and 1b
ax1.plot([-1+shift,1-shift],[1,1], '--', color = '#FFD92F', linewidth = 1)
ax1.plot([-1+shift,1-shift],[-1,-1], '--', color = '#FFD92F', linewidth = 1)
#Branch 2
ax1.plot([-1,1],[0,0], color = '#B3B3B3', linewidth = 1)
#Sub branches 2a and 2b
ax1.plot([1,1],[-1+shift,1-shift], '--', color = '#FFD92F', linewidth = 1)
ax1.plot([-1,-1],[-1+shift,1-shift], '--', color = '#FFD92F', linewidth = 1)

#Dome State
ax1.scatter([das, -das],[das, -das], color = '#66C2A5', marker = 'o', s = 20)
#Saddle State
ax1.scatter([das, -das],[-das, das], color = '#8DA0CB', marker = 'o', s = 20)
#MiuraOriState
ax1.scatter([mo, mo, -mo, -mo],[mo, -mo, mo, -mo], color = '#FFD92F', marker = 'v', s = 20)

fig1.savefig('D:/Documents/Git Programs/nonlinear-bas_Origami/Results/BranchPlot_1.pdf', transparent = True)
fig1.savefig('D:/Documents/Git Programs/nonlinear-bas_Origami/Results/BranchPlot_1.png', transparent = True)

#%%
fig2 = plt.figure(figsize=(cm2inch(4.4), cm2inch(3.1)))
ax2 = plt.subplot(111)
fig2.subplots_adjust(top=0.982,
bottom=0.23,
left=0.280,
right=0.940)

NiceGraph2D(ax2, r'$\theta_1/\pi$', r'$\theta_3/\pi$', mincoord = [-1,-1], maxcoord = [1,1], 
            divisions = [3,3], buffer = [0.1, 0.1])

#Branch 1
ax2.plot([-1,1],[-1,1], color = '#B3B3B3', linewidth = 1, zorder = 10)
#Sub branches 1a and 1b
ax2.scatter([-1,1],[-1,1], marker = 'x', color = '#FFD92F', s = 20)
#Branch 2
ax2.scatter([0],[0], marker = 'o', color = '#B3B3B3', s = 20, zorder = 14)
#Sub branches 2a and 2b
ax2.plot([-1+shift/np.sqrt(2)/2,1+shift/np.sqrt(2)/2],[1+shift/np.sqrt(2)/2,-1+shift/np.sqrt(2)/2], '--', color = '#FFD92F', linewidth = 1, zorder = 11)
ax2.plot([-1-shift/np.sqrt(2)/2,1-shift/np.sqrt(2)/2],[1-shift/np.sqrt(2)/2,-1-shift/np.sqrt(2)/2], '--', color = '#FFD92F', linewidth = 1, zorder = 12)

#Dome State
ax2.scatter([das, -das],[das, -das], color = '#66C2A5', marker = 'o', s = 20, zorder = 15)
# #Saddle State
# ax2.scatter([das, -das],[-das, das], color = '#8DA0CB', marker = 'o', s = 20)
#MiuraOriState
ax2.scatter([mo, mo, -mo, -mo],[mo, -mo, mo, -mo], color = '#FFD92F', marker = 'v', s = 20, zorder = 20)

fig2.savefig('D:/Documents/Git Programs/nonlinear-bas_Origami/Results/BranchPlot_2.pdf', transparent = True)
fig2.savefig('D:/Documents/Git Programs/nonlinear-bas_Origami/Results/BranchPlot_2.png', transparent = True)
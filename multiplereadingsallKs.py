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
import matplotlib.ticker as ticker

totSimulations = 0
internalHinges = 9
flags = np.empty((6, 8, 8, 6)) # 12 internal hinges, 6 possible flags
flags[:] = np.NaN
ststs = np.empty((6, 8, 8))
ststs[:] = np.NaN
folder_name = "Results/triangular prism/sqp/energy/"

uniqueNames = np.empty(0)
uniqueEnergies = np.empty((0,2))
uniqueKs = np.empty((0,3), dtype = int)
simulations = 0

for stephinge in np.arange(6):
    khinge = 0.001*10**(stephinge/2);
    for stepangle in np.arange(8):
        kangle = 0.001*10**(0.5+stepangle/2);
        for stepedge in np.arange(8):
            kedge = 0.001*10**(0.5+stepedge/2)
            subfolder_name = "kh%2.3f_kta%2.3f_ke%2.3f/" %(khinge, kangle, kedge)
            names = np.nan
            energies = np.nan
            if os.path.exists(folder_name+subfolder_name):
                totSimulations += 1
                print(folder_name+subfolder_name)
                flags[stephinge,stepangle,stepedge,:], ststs[stephinge,stepangle,stepedge], names, energies= pc2.ReadandAnalizeFile(folder_name+subfolder_name, plot = False, khinge = khinge, kedge = kedge)
                uniqueNames = np.append(uniqueNames, names)
                uniqueEnergies = np.append(uniqueEnergies, energies, axis = 0)
                for i in np.arange(len(energies)):
                    uniqueKs = np.append(uniqueKs, [[stephinge, stepangle, stepedge]], axis = 0)
                simulations += 1
                
            


#%%
allflags = np.sum(flags, axis = 3)

x = 0.001*10**(np.arange(6)/2)
y = 0.001*10**(0.5+np.arange(8)/2)
x2 = 0.001*10**(np.arange(7)/2-0.25)
y2 =0.001*10**(0.5+np.arange(9)/2-0.25)

fig1, axes1 = plt.subplots(nrows=2, ncols=4, figsize=(pc2.cm2inch(48), pc2.cm2inch(23)))
fig1.subplots_adjust(wspace = 0, hspace = 0.15, left = 0.06, right = 0.98, top = 0.95, bottom = 0.08)

i = 0
for ax, edge in zip(axes1.flat, y):
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_title('kEdge %2.3f' %(edge), fontsize = 15, color = '0.2')
    pc2.NiceGraph2D(ax, 'kHinge', 'kTargetAngle', mincoord = [np.min(x)*10**-0.25, np.min(y)*10**-0.25], maxcoord = [np.max(x)*10**0.25, np.max(y)*10**0.25], divisions = [x,y])
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%2.3f'))
    if i%4 != 0:
        ax.set_yticks([])
        ax.set_ylabel('')
    if i < 4:
        ax.set_xticks([])
        ax.set_xlabel('')
    i += 1


cmap, norm = from_levels_and_colors(np.linspace(0,1,11), cm.gnuplot_r(np.linspace(0, 1, 10)))
allflags = np.ma.masked_where(np.isnan(allflags), allflags)
flags = np.ma.masked_where(np.isnan(flags), flags)

for edge, ax in zip(np.arange(8), axes1.flat):
    cs1 = ax.pcolormesh(x2, y2, allflags[:,:,edge].T, cmap=cmap, norm = norm)
#axes1[0,1].pcolormesh(x2, y2, flags[:,:,1], cmap=cmap, norm = norm)
#axes1[1,0].pcolormesh(x2, y2, flags[:,:,3], cmap=cmap, norm = norm)
#axes1[1,1].pcolormesh(x2, y2, flags[:,:,5], cmap=cmap, norm = norm)
#
cbar = plt.colorbar(cs1, format="%.2f", ax =axes1.ravel().tolist(), fraction=0.05, pad=0.01)
cbar.set_label('Non-Convergence rate', fontsize = 15, color = '0.2')
cbar.ax.tick_params(axis='y',colors='0.2')
cbar.ax.tick_params(axis='x',colors='0.2')
cbar.outline.set_edgecolor('0.2')

fig1.show()
#fig1.savefig(folder_name + 'Convergencekallwo2flags2Flag.png', transparent = True)

#%%

fig2, axes2 = plt.subplots(nrows=2, ncols=4, figsize=(pc2.cm2inch(48), pc2.cm2inch(23)))
fig2.subplots_adjust(wspace = 0, hspace = 0.15, left = 0.06, right = 0.98, top = 0.95, bottom = 0.08)

i = 0
for ax, edge in zip(axes2.flat, y):
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_title('kEdge %2.3f' %(edge), fontsize = 15, color = '0.2')
    pc2.NiceGraph2D(ax, 'kHinge', 'kAngle', mincoord = [np.min(x)*10**-0.25, np.min(y)*10**-0.25], maxcoord = [np.max(x)*10**0.25, np.max(y)*10**0.25], divisions = [x,y])
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%2.3f'))
    if i%4 != 0:
        ax.set_yticks([])
        ax.set_ylabel('')
    if i < 4:
        ax.set_xticks([])
        ax.set_xlabel('')
    i += 1

cmap2, norm2 = from_levels_and_colors(np.arange(22+1)+1, cm.jet(np.linspace(0, 1, 22)))
ststs = np.ma.masked_where(np.isnan(ststs), ststs)

for edge, ax in zip(np.arange(8), axes2.flat):
    cs2 = ax.pcolormesh(x2, y2, ststs[:,:,edge].T, cmap=cmap2, norm = norm2)

cbar2 = plt.colorbar(cs2, format="%d", ax =axes2.ravel().tolist(), fraction=0.05, pad=0.01)
cbar2.set_label('Number of Stable States', fontsize = 15, color = '0.2')
cbar2.ax.tick_params(axis='y',colors='0.2')
cbar2.ax.tick_params(axis='x',colors='0.2')
cbar2.outline.set_edgecolor('0.2')

fig2.show()
#fig2.savefig(folder_name + 'StableStateskallwo2flags.png', transparent = True)

#%%
tolHinge = 0.019
tolEdge = 0.0017
maxHinge = 2.1#np.max(uniqueEnergies[:,0])#120
maxEdge = 0.15#np.max(uniqueEnergies[:,1])#10


fig3, axes3 = plt.subplots(nrows=1, ncols=1, figsize=(pc2.cm2inch(35), pc2.cm2inch(20)))

pc2.NiceGraph2D(axes3, 'Edge Energy',  'Hinge Energy')

cmap = cm.Set2
cmap.set_over('r')
cmap.set_under('0.2')

cs3 = axes3.hist2d(uniqueEnergies[:,1], uniqueEnergies[:,0], bins = [np.arange(0,maxEdge+tolEdge,tolEdge),np.arange(0,maxHinge+tolHinge,tolHinge)],
                   cmap = cm.Set2, vmax = 50, vmin = 3)
axes3.scatter(uniqueEnergies[:,1], uniqueEnergies[:,0], 0.5, '0.5')

cbar3 = plt.colorbar(cs3[3],format="%d", fraction=0.05, pad=0.01, extend='both')
cbar3.set_label('Counts of Simulations', fontsize = 15, color = '0.2')
cbar3.ax.tick_params(axis='y',colors='0.2')
cbar3.ax.tick_params(axis='x',colors='0.2')
cbar3.outline.set_edgecolor('0.2')

goodstst = np.zeros((6,8,8))
goodststName = np.empty((6,8,8,0))
for i in np.arange(cs3[0].shape[0]):
    for j in np.arange(cs3[0].shape[1]):
        if cs3[0][i][j]>=3:
            same = np.where(np.logical_and(np.logical_and(uniqueEnergies[:,1]>=cs3[1][i],
                                                          uniqueEnergies[:,1]<cs3[1][i+1]),
                                           np.logical_and(uniqueEnergies[:,0]>=cs3[2][j], 
                                                          uniqueEnergies[:,0]<cs3[2][j+1])))[0]
            print('Counts', len(same), '\tEnergies', cs3[1][i], cs3[2][j])
#            p = 0
            for stablest in same:
#                if cs3[1][i]<uniqueEnergies[stablest,1] and uniqueEnergies[stablest,1]<cs3[1][i+1] and cs3[2][j]<uniqueEnergies[stablest,0] and uniqueEnergies[stablest,0]< cs3[2][j+1]:
#                    p += 1
#                else:
#                    print(cs3[1][i],uniqueEnergies[stablest,1],cs3[1][i+1])
#                    print(cs3[2][j],uniqueEnergies[stablest,0],cs3[2][j+1])
                        
                goodstst[uniqueKs[stablest][0]][uniqueKs[stablest][1]][uniqueKs[stablest][2]] += 1
                goodststName[uniqueKs[stablest][0]][uniqueKs[stablest][1]][uniqueKs[stablest][2]] = np.append(
                        goodststName[uniqueKs[stablest][0]][uniqueKs[stablest][1]][uniqueKs[stablest][2]], uniqueNames[stablest])
            print(stablest, uniqueNames[stablest], uniqueKs[stablest][0], uniqueKs[stablest][1], uniqueKs[stablest][2])

fig3.show()
#fig3.savefig(folder_name + 'AllStateskallwo2flags.png', transparent = True)

#%%

fig4, axes4 = plt.subplots(nrows=2, ncols=4, figsize=(pc2.cm2inch(48), pc2.cm2inch(23)))
fig4.subplots_adjust(wspace = 0, hspace = 0.15, left = 0.06, right = 0.98, top = 0.95, bottom = 0.08)

i = 0
for ax, edge in zip(axes4.flat, y):
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_title('kEdge %2.3f' %(edge), fontsize = 15, color = '0.2')
    pc2.NiceGraph2D(ax, 'kHinge', 'kAngle', mincoord = [np.min(x)*10**-0.25, np.min(y)*10**-0.25], maxcoord = [np.max(x)*10**0.25, np.max(y)*10**0.25], divisions = [x,y])
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%2.3f'))
    if i%4 != 0:
        ax.set_yticks([])
        ax.set_ylabel('')
    if i < 4:
        ax.set_xticks([])
        ax.set_xlabel('')
    i += 1

cmap4, norm4 = from_levels_and_colors(np.arange(np.max(goodstst)+1)+1, cm.jet(np.linspace(0, 1, np.max(goodstst))))

for edge, ax in zip(np.arange(8), axes4.flat):
    cs4 = ax.pcolormesh(x2, y2, goodstst[:,:,edge].T, cmap=cmap4, norm = norm4)

cbar4 = plt.colorbar(cs4, format="%d", ax =axes4.ravel().tolist(), fraction=0.05, pad=0.01)
cbar4.set_label('Number of Stable States', fontsize = 15, color = '0.2')
cbar4.ax.tick_params(axis='y',colors='0.2')
cbar4.ax.tick_params(axis='x',colors='0.2')
cbar4.outline.set_edgecolor('0.2')

fig4.show()
#fig4.savefig(folder_name + 'StableStatesgoodkallwo2flags.png', transparent = True)
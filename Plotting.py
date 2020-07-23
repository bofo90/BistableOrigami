# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 14:53:45 2020

@author: iniguez
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as matl
from matplotlib.colors import from_levels_and_colors
from mpl_toolkits import mplot3d
import ternary #from https://github.com/marcharper/python-ternary
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

def NiceGraph2Dlog(axes, nameX, nameY, mincoord = [np.NaN, np.NaN], maxcoord = [np.NaN, np.NaN], divisions = [np.NaN, np.NaN],buffer = [0.0, 0.0, 0.0]):
    
    gray = '0.2'
    matl.rcParams.update({'font.size': 9})

    if ~np.isnan(mincoord[0]) and ~np.isnan(maxcoord[0]):
        axes.set_xlim([mincoord[0]*((buffer[0])**(-1)), maxcoord[0]*(buffer[0])])
        if isinstance(divisions[0], (list, tuple, np.ndarray)):
            if ~np.isnan(divisions[0]).any():
                axes.set_xticks(divisions[0])
        else:
            if ~np.isnan(divisions[0]):
                axes.set_xticks(np.logspace(np.log10(mincoord[0]),np.log10(maxcoord[0]),divisions[0]))
#    axes.set_xscale('log')
            
    axes.set_xlabel(nameX)
    if ~np.isnan(mincoord[1]) and ~np.isnan(maxcoord[1]):
        axes.set_ylim([mincoord[1]-buffer[1], maxcoord[1]+buffer[1]])
        if isinstance(divisions[1], (list, tuple, np.ndarray)):
            if ~np.isnan(divisions[1]).any():
                axes.set_yticks(divisions[1])
        else:
            if ~np.isnan(divisions[1]):
                axes.set_yticks(np.linspace(mincoord[1],maxcoord[1],divisions[1]))
                            
    axes.set_ylabel(nameY)
   
    axes.xaxis.label.set_color(gray)
    axes.tick_params(axis='x',colors=gray)
    axes.spines['bottom'].set_color(gray)
    axes.spines['top'].set_color(gray)
    axes.yaxis.label.set_color(gray)
    axes.tick_params(axis='y', colors=gray)
    axes.spines['left'].set_color(gray)
    axes.spines['right'].set_color(gray)
    axes.tick_params(pad = 2)
    
    for axis in ['top','bottom','left','right']:
        axes.spines[axis].set_linewidth(0.4)
        
    return

def NiceTerciaryGraph(ax, name, scale, divisions):
    
    ax.axis('off')
    tax = ternary.TernaryAxesSubplot(ax=ax, scale = scale)
    
    fontsize = 9
    matl.rcParams.update({'font.size': 9})
    gray = '0.2'
    linewidth = 0.4
    mult = scale/divisions
    
    tax.boundary(linewidth=linewidth)
    tax.gridlines(color=gray, multiple=mult)
    
    tax.set_title(name, fontsize=fontsize, color = gray, pad = 15)
    tax.left_axis_label("Angle1", fontsize=fontsize, color = gray, offset = 0.3)
    tax.right_axis_label("Angle2", fontsize=fontsize, color = gray, offset = 0.3)
    tax.bottom_axis_label("Angle3", fontsize=fontsize, color = gray, offset = 0.3)
    
    tax.ticks(axis='lbr', linewidth=linewidth, axes_colors = {'l': gray, 'r':gray, 'b': gray},fontsize = fontsize,
              ticks = [180,150,120,90,60,30,0], clockwise = True, offset = 0.05)
    tax.clear_matplotlib_ticks()
    
    return tax

def getcolormap(allDesigns, colormap):
    
    stst = np.unique(allDesigns)
    cmap = matl.cm.get_cmap(colormap,np.size(stst))
    color = cmap(np.linspace(0,1,np.size(stst)))    
    
    return stst, color

def XYperZ(allDesigns, x, xname, y, yname, z, stst_col, colormap, save = False, Folder_name = '', NameFig = ''):
    
    allDesigns = allDesigns.round(8)
    
    stst, color = getcolormap(allDesigns.iloc[:,stst_col].values, colormap)
    
    variables = np.unique(allDesigns.iloc[:,z])
    
    for i in np.arange(np.size(variables)):
    
        # fig1 = plt.figure(figsize=(cm2inch(4.3), cm2inch(3.1)))
        fig1 = plt.figure(figsize=(cm2inch(8), cm2inch(6)))
        ax1 = plt.subplot(111)
        fig1.subplots_adjust(top=0.982,
        bottom=0.23,
        left=0.280,
        right=0.940)
        
        thisDesBool = allDesigns.iloc[:,z] == variables[i]
        thisDes = allDesigns[thisDesBool]    
        
        NiceGraph2D(ax1, xname, yname)
        
        if yname == r'$E_{tot}$':
            ax1.set_ylim([-0.005, np.max(thisDes.iloc[:,y])+0.005])
        
        if yname == r'$Area$':
            ax1.set_ylim([0, 0.8])
            ax1.axhline(y=0.1, color='r', linestyle='-', linewidth = '0.4')
            
        if xname == r'$\kappa$':
            ax1.set_xscale('log')
            ax1.set_xticks([0.001,0.01,0.1,1])
            ax1.set_xlim([0.0007,1.5])
        
        if xname == r'$\theta_0/\pi$':
            ax1.set_xticks([0,0.5,1])
            ax1.set_xlim([-0.01,1.1])
            
        for k in np.arange(np.size(stst)):
            thisstst = thisDes[thisDes.iloc[:, stst_col] == stst[k]]
            
            ax1.scatter(thisstst.iloc[:,x], thisstst.iloc[:,y], c = matl.colors.rgb2hex(color[k]), s = 8)
        
        fig1.show()
        if save:
            fig1.savefig(Folder_name + '/Images/' + NameFig + '_' + '%.4f' %variables[i] +'.pdf', transparent = True)
            fig1.savefig(Folder_name + '/Images/' + NameFig + '_' + '%.4f' %variables[i] +'.png', transparent = True)
            print(Folder_name + '/Images/' + NameFig + '_' + '%.4f' %variables[i])
    
    return

def XmultYperZ(allDesigns, x, xname, y, yname, z, save = False, Folder_name = '', NameFig = ''):
    
    allDesigns = allDesigns.round(8)
    
    # stst, color = getcolormap(allDesigns.iloc[:,stst_col].values, colormap)
    
    variables = np.unique(allDesigns.iloc[:,z])
    
    for i in np.arange(np.size(variables)):
    
        fig1 = plt.figure(figsize=(cm2inch(4.3), cm2inch(3.1)))
        ax1 = plt.subplot(111)
        fig1.subplots_adjust(top=0.982,
        bottom=0.23,
        left=0.225,
        right=0.940)
        
        thisDesBool = allDesigns.iloc[:,z] == variables[i]
        thisDes = allDesigns[thisDesBool]    
        
        NiceGraph2D(ax1, xname, yname)
                    
        if xname == r'$\kappa$':
            ax1.set_xscale('log')
            ax1.set_xticks([0.001,0.01,0.1,1])
            ax1.set_xlim([0.0007,1.5])
            
        for k in y:
            ax1.scatter(thisDes.iloc[:,x], thisDes.iloc[:,k], s = 8, label = allDesigns.columns[k])
        
        CreateLegend(ax1)
        fig1.show()
               
        if save:
            fig1.savefig(Folder_name + '/Images/' + NameFig + '_' + '%.4f' %variables[i] +'.pdf', transparent = True)
            fig1.savefig(Folder_name + '/Images/' + NameFig + '_' + '%.4f' %variables[i] +'.png', transparent = True)
            print(Folder_name + '/Images/' + NameFig + '_' + '%.4f' %variables[i])
    
    return

def Angles3D(angles, angStSt, colormap, save = False, Folder_name = ''):
        
    stst, color = getcolormap(angStSt, colormap)
    
    stst, inverse = np.unique(angStSt, return_inverse = True)
    
    fig1 = plt.figure(figsize=(cm2inch(10), cm2inch(7)))
    ax1 = plt.subplot(111,projection='3d')
    ax1.set_xlim([-np.pi,np.pi])
    ax1.set_ylim([-np.pi,np.pi])
    ax1.set_zlim([-np.pi,np.pi])
            
    ax1.scatter(angles[:,0],angles[:,1],angles[:,2], c = color[inverse])
        
    for i in np.arange(np.size(stst)):
        ax1.scatter(-400,-400,-400, c = [color[i]], label = stst[i])
    
    CreateLegend(ax1)
    
    if save:
        fig1.savefig(Folder_name + '/Images/Angles3D.pdf', transparent = True)
        fig1.savefig(Folder_name + '/Images/Angles3D.png', transparent = True)
    
    
    return

def CreateColorbar(StableStates, cmap, save = False, Folder_name = '', NameFig = ''):
    
    StableStates = np.round(StableStates,7)
    
    
    stst, nameloc = np.unique(StableStates, return_index = True)
    ststname = StableStates[nameloc]
    cmaptemp = matl.cm.get_cmap(cmap)
    
    cmapfig, normfig = from_levels_and_colors(np.arange(np.size(stst)+1),cmaptemp(np.linspace(0, 1,np.size(stst))))

    fig1 = plt.figure(figsize=(cm2inch(8.8), cm2inch(1.5)))
    fig1.subplots_adjust(top=0.884,
    bottom=0.116,
    left=0.039,
    right=0.961)
    ax1 = plt.subplot(111)
    cbar = plt.colorbar(matl.cm.ScalarMappable(norm=normfig, cmap=cmapfig), ax = ax1, fraction=0.99, pad=0.01, orientation='horizontal')
    cbar.set_ticks(np.arange(np.size(stst))+0.5)
    cbar.ax.set_xticklabels(ststname.astype(int))
    cbar.set_label('Stable State', fontsize = 9, color = '0.2',labelpad = 0)
    cbar.ax.tick_params(colors='0.2', pad=2)
    cbar.outline.set_edgecolor('0.2')
    cbar.outline.set_linewidth(0.4)
    ax1.remove()   
    
    if save: 
        fig1.savefig(Folder_name + '/Images/' + NameFig + '_CB.pdf', transparent = True)
        fig1.savefig(Folder_name + '/Images/' + NameFig + '_CB.png', transparent = True)
    
    return

def CreateContColorbar(Data, Name, cmap, save = False, Folder_name = '', NameFig = ''):
    
    minDat = np.min(Data)
    maxDat = np.max(Data)
    
    cmaptemp = matl.cm.get_cmap(cmap)
    cmapfig, normfig = from_levels_and_colors(np.linspace(minDat, maxDat, 1001),cmaptemp(np.linspace(0, 1,1000)))

    fig1 = plt.figure(figsize=(cm2inch(8.8), cm2inch(1.5)))
    fig1.subplots_adjust(top=0.884,
    bottom=0.116,
    left=0.059,
    right=0.951)
    ax1 = plt.subplot(111)
    cbar = plt.colorbar(matl.cm.ScalarMappable(norm=normfig, cmap=cmapfig), ax = ax1, fraction=0.99, pad=0.01, orientation='horizontal')
    cbar.set_ticks(np.linspace(minDat,maxDat,5))
    # cbar.ax.set_xticklabels(ststname.astype(int))
    cbar.set_label(Name, fontsize = 9, color = '0.2',labelpad = 0)
    cbar.ax.tick_params(colors='0.2', pad=2)
    cbar.outline.set_edgecolor('0.2')
    cbar.outline.set_linewidth(0.4)
    ax1.remove()   
    
    if save: 
        fig1.savefig(Folder_name + '/Images/' + NameFig + '_CB.pdf', transparent = True)
        fig1.savefig(Folder_name + '/Images/' + NameFig + '_CB.png', transparent = True)
    
    return

def CreateLegend(ax, location = 0):
    
    leg = ax.legend(loc = location, fontsize = 7, framealpha = 0.8, edgecolor = 'inherit', fancybox = False, 
               borderpad = 0.3, labelspacing = 0.1, handlelength = 0.4, handletextpad = 0.4)
    plt.setp(leg.get_texts(), color='0.2')
    leg.get_frame().set_linewidth(0.4)

    return

def ConvSim(redDesign,redFlags, colormap, save = False, Folder_name = '', NameFig = ''):
    
    flags = np.unique(redFlags['Flags'])
    flags = flags[flags < 0]
    
    stst = np.unique(redDesign['StableStateAll'])
    
    stddev = np.unique(redDesign['StdDev'])
    
    totsim = np.zeros((np.size(stddev), np.size(flags)+np.size(stst)))
    
    for sd, i in zip(stddev, np.arange(np.size(stddev))):
        j = 0
        for f in flags:
            hereFlag = (redFlags['StdDev'] == sd) & (redFlags['Flags'] == f)
            hereFlag = hereFlag.to_numpy()
            if np.sum(hereFlag)>1:
                print('\nError with uniquenes of flag and stddev\n')
            elif np.sum(hereFlag)== 1:
                totsim[i,j] = redFlags.loc[hereFlag,'amountFlags'].values
            j += 1
        for s in stst:
            hereStSt = (redDesign['StdDev'] == sd) & (redDesign['StableStateAll'] == s)
            hereStSt = hereStSt.to_numpy()
            if np.sum(hereStSt) >= 1:
                totsim[i,j] = np.sum(redDesign.loc[hereStSt,'amountStSt'].values)
            j += 1
    
    allStates = np.concatenate((flags,stst))
    
    cmap = matl.cm.get_cmap(colormap,np.size(allStates))
    color = cmap(np.linspace(0,1,np.size(allStates)))   
    
    fig = plt.figure(figsize=(cm2inch(8), cm2inch(6)))
    ax1 = plt.subplot(111)
    fig.subplots_adjust(top=0.968,
bottom=0.099,
left=0.136,
right=0.982)
    
    NiceGraph2D(ax1, 'Std.Dev.', 'Sim. Amount')
    
    for i in np.arange(np.size(totsim,1)):
        ax1.bar(stddev, totsim[:,i],bottom=np.sum(totsim[:,:i],axis = 1), width = 0.05,
                color = color[i], label = allStates[i], align = 'center')
    
    CreateLegend(ax1, location = 2)
    
    if save: 
        fig.savefig(Folder_name + '/Images/' + NameFig + '_CB.pdf', transparent = True)
        fig.savefig(Folder_name + '/Images/' + NameFig + '_CB.png', transparent = True)
    return

def TotEnergyperZ(allDesigns, x, xname, z, stst_col, colormap, save = False, Folder_name = ''):
    
    allDesigns = allDesigns.round(8)
    
    numhinge = 4
    numedge = 8
    
    stst, color = getcolormap(allDesigns.iloc[:,stst_col].values, colormap)
    
    variables = np.unique(allDesigns.iloc[:,z])
    
    for i in np.arange(np.size(variables)):
    
        fig1 = plt.figure(figsize=(cm2inch(4.3), cm2inch(3.1)))
        ax1 = plt.subplot(111)
        fig1.subplots_adjust(top=0.982,
        bottom=0.23,
        left=0.280,
        right=0.940)
        
        thisDesBool = allDesigns.iloc[:,z] == variables[i]
        thisDes = allDesigns[thisDesBool]    
        
        NiceGraph2D(ax1, xname, r'$E_{tot}$')
        
        Energy_flatSt = thisDes.iloc[:,0]*numhinge/(thisDes.iloc[:,0]*numhinge+numedge)
        Energy_foldSt = (2*thisDes.iloc[:,0])/(thisDes.iloc[:,0]*numhinge+numedge)
        
        ax1.set_ylim([-np.max(Energy_flatSt)*0.3, np.max(Energy_flatSt)*1.07])
            
        if xname == r'$\kappa$':
            ax1.set_xscale('log')
            ax1.set_xticks([0.001,0.01,0.1,1])
            ax1.set_xlim([0.0007,1.5])
            
        ax1.plot(thisDes.iloc[:,x], Energy_flatSt, c = '#36648B', label = 'Flat')
        ax1.plot(thisDes.iloc[:,x], Energy_foldSt, c = '#D75D4A', label = 'Flat')
           
        for k in np.arange(np.size(stst)):
            thisstst = thisDes[thisDes.iloc[:, stst_col] == stst[k]]
            
            ax1.scatter(thisstst.iloc[:,x], thisstst.iloc[:,5], c = matl.colors.rgb2hex(color[k]), s = 8)
            
            
        fig1.show()
        if save:
            fig1.savefig(Folder_name + '/Images/' + 'TotEnergyNorm' + '_' + '%.4f' %variables[i] +'.pdf', transparent = True)
            fig1.savefig(Folder_name + '/Images/' + 'TotEnergyNorm' + '_' + '%.4f' %variables[i] +'.png', transparent = True)
            print(Folder_name + '/Images/' + 'TotEnergyNorm' + '_' + '%.4f' %variables[i])
    
    return

def XYperZwError(allDesigns, x, xname, y, yname, z, stst_col, colormap, yerror, save = False, Folder_name = '', NameFig = ''):
    
    allDesigns = allDesigns.round(8)
    
    stst, color = getcolormap(allDesigns.iloc[:,stst_col].values, colormap)
    
    variables = np.unique(allDesigns.iloc[:,z])
    
    for i in np.arange(np.size(variables)):
    
        fig1 = plt.figure(figsize=(cm2inch(8), cm2inch(6)))
        ax1 = plt.subplot(111)
        fig1.subplots_adjust(top=0.942,
bottom=0.176,
left=0.165,
right=0.957)
        
        thisDesBool = allDesigns.iloc[:,z] == variables[i]
        thisDes = allDesigns[thisDesBool]    
        
        NiceGraph2D(ax1, xname, yname)
        
        if yname == r'$E_{tot}$':
            ax1.set_ylim([-0.005, np.max(thisDes.iloc[:,y])+0.005])
        
        if yname == r'$Area$':
            ax1.set_ylim([0, 0.8])
            ax1.axhline(y=0.1, color='r', linestyle='-', linewidth = '0.4')
            
        if xname == r'$\kappa$':
            ax1.set_xscale('log')
            ax1.set_xticks([0.001,0.01,0.1,1])
            ax1.set_xlim([0.0007,1.5])
        
        if xname == r'$\theta_0/\pi$':
            ax1.set_xticks([0,0.5,1])
            ax1.set_xlim([-0.01,1.1])
            
        for k in np.arange(np.size(stst)):
            thisstst = thisDes[thisDes.iloc[:, stst_col] == stst[k]]
            
            ax1.errorbar(thisstst.iloc[:,x], thisstst.iloc[:,y], yerr = thisstst.iloc[:,yerror], 
                         fmt='o', c = matl.colors.rgb2hex(color[k]), ms = 3, capsize = 3)
        
        fig1.show()
        if save:
            fig1.savefig(Folder_name + '/Images/' + NameFig + '_' + '%.4f' %variables[i] +'.pdf', transparent = True)
            fig1.savefig(Folder_name + '/Images/' + NameFig + '_' + '%.4f' %variables[i] +'.png', transparent = True)
            print(Folder_name + '/Images/' + NameFig + '_' + '%.4f' %variables[i])
    
    return

def XYperZwDoubleError(allDesigns, x, xname, y, yname, z, stst_col, colormap, xerror, yerror, save = False, Folder_name = '', NameFig = ''):
    
    allDesigns = allDesigns.round(8)
    
    stst, color = getcolormap(allDesigns.iloc[:,stst_col].values, colormap)
    
    variables = np.unique(allDesigns.iloc[:,z])
    
    for i in np.arange(np.size(variables)):
    
        fig1 = plt.figure(figsize=(cm2inch(8), cm2inch(6)))
        ax1 = plt.subplot(111)
        fig1.subplots_adjust(top=0.942,
bottom=0.176,
left=0.165,
right=0.957)
        
        thisDesBool = allDesigns.iloc[:,z] == variables[i]
        thisDes = allDesigns[thisDesBool]    
        
        NiceGraph2D(ax1, xname, yname)
        
        if yname == r'$E_{tot}$':
            ax1.set_ylim([-0.005, np.max(thisDes.iloc[:,y])+0.005])
        
        if yname == r'$Area$':
            ax1.set_ylim([0, 0.8])
            ax1.axhline(y=0.1, color='r', linestyle='-', linewidth = '0.4')
            
        if xname == r'$\kappa$':
            ax1.set_xscale('log')
            ax1.set_xticks([0.001,0.01,0.1,1])
            ax1.set_xlim([0.0007,1.5])
        
        if xname == r'$\theta_0/\pi$':
            ax1.set_xticks([0,0.5,1])
            ax1.set_xlim([-0.01,1.1])
            
        for k in np.arange(np.size(stst)):
            thisstst = thisDes[thisDes.iloc[:, stst_col] == stst[k]]
            
            ax1.errorbar(thisstst.iloc[:,x], thisstst.iloc[:,y], xerr = thisstst.iloc[:,xerror],yerr = thisstst.iloc[:,yerror], 
                         fmt='o', c = matl.colors.rgb2hex(color[k]), ms = 3, capsize = 3)
        
        fig1.show()
        if save:
            fig1.savefig(Folder_name + '/Images/' + NameFig + '_' + '%.4f' %variables[i] +'.pdf', transparent = True)
            fig1.savefig(Folder_name + '/Images/' + NameFig + '_' + '%.4f' %variables[i] +'.png', transparent = True)
            print(Folder_name + '/Images/' + NameFig + '_' + '%.4f' %variables[i])
    
    return

def XYperZline(allDesigns, x, xname, y, yname, z, stst_col, colormap, save = False, Folder_name = '', NameFig = ''):
    
    allDesigns = allDesigns.round(8)
    
    stst, color = getcolormap(allDesigns.iloc[:,stst_col].values, colormap)
    
    variables = np.unique(allDesigns.iloc[:,z])
    
    for i in np.arange(np.size(variables)):
    
        fig1 = plt.figure(figsize=(cm2inch(8), cm2inch(6)))
        ax1 = plt.subplot(111)
        fig1.subplots_adjust(top=0.942,
bottom=0.176,
left=0.165,
right=0.957)
        
        thisDesBool = allDesigns.iloc[:,z] == variables[i]
        thisDes = allDesigns[thisDesBool]    
        
        NiceGraph2D(ax1, xname, yname)
        
        if yname == r'$E_{tot}$':
            ax1.set_ylim([-0.005, np.max(thisDes.iloc[:,y])+0.005])
        
        if yname == r'$Area$':
            ax1.set_ylim([0, 0.8])
            ax1.axhline(y=0.1, color='r', linestyle='-', linewidth = '0.4')
            
        if xname == r'$\kappa$':
            ax1.set_xscale('log')
            ax1.set_xticks([0.001,0.01,0.1,1])
            ax1.set_xlim([0.0007,1.5])
        
        if xname == r'$\theta_0/\pi$':
            ax1.set_xticks([0,0.5,1])
            ax1.set_xlim([-0.01,1.1])
            
        for k in np.arange(np.size(stst)):
            thisstst = thisDes[thisDes.iloc[:, stst_col] == stst[k]]
            
            ax1.plot(thisstst.iloc[:,x], thisstst.iloc[:,y], c = matl.colors.rgb2hex(color[k]), linewidth = 0.4)
        
        fig1.show()
        if save:
            fig1.savefig(Folder_name + '/Images/' + NameFig + '_' + '%.4f' %variables[i] +'.pdf', transparent = True)
            fig1.savefig(Folder_name + '/Images/' + NameFig + '_' + '%.4f' %variables[i] +'.png', transparent = True)
            print(Folder_name + '/Images/' + NameFig + '_' + '%.4f' %variables[i])
    
    return

def CurvaturePaper(allDesigns, x, xname, y, yname, z, stst_col, colormap, save = False, Folder_name = '', NameFig = ''):
    
    allDesigns = allDesigns.round(8)
    
    stst, color = getcolormap(allDesigns.iloc[:,stst_col].values, colormap)
    
    variables = np.unique(allDesigns.iloc[:,z])
    
    for i in np.arange(np.size(variables)):
    
        fig1 = plt.figure(figsize=(cm2inch(4.4), cm2inch(3.1)))
        ax1 = plt.subplot(111)
        fig1.subplots_adjust(top=0.982,
        bottom=0.23,
        left=0.280,
        right=0.940)
        
        thisDesBool = allDesigns.iloc[:,z] == variables[i]
        thisDes = allDesigns[thisDesBool]    
        
        NiceGraph2D(ax1, xname, yname)
        
        ax1.set_xticks([0,0.5,1])
        ax1.set_xlim([-0.01,1.01])
            
        for k in np.arange(np.size(stst))[::-1]:
            thisstst = thisDes[thisDes.iloc[:, stst_col] == stst[k]]
            
            ax1.scatter(thisstst.iloc[:,x], thisstst.iloc[:,y], c = matl.colors.rgb2hex(color[k]), s = 8)
        
        fig1.show()
        if save:
            fig1.savefig(Folder_name + '/Images/' + NameFig + '_' + '%.4f' %variables[i] +'.pdf', transparent = True)
            fig1.savefig(Folder_name + '/Images/' + NameFig + '_' + '%.4f' %variables[i] +'.png', transparent = True)
            print(Folder_name + '/Images/' + NameFig + '_' + '%.4f' %variables[i])
    
    return

def StableStatesCounturPlot(allDesigns,x, xname, y, yname, z, color, colorName, colormap, minmax, save = False, Folder_name = '', NameFig = ''):
    
    allDesigns = allDesigns.round(8)
    
    # stst, color = getcolormap(allDesigns.iloc[:,z].values, colormap)
    # cmaptemp = matl.cm.get_cmap(colormap)
    # cmaptemp.set_over('r')
    # cmaptemp.set_under('b')
    stst = np.unique(allDesigns.iloc[:,z])
    
    minCol = minmax[0]#np.min(allDesigns.iloc[:,color])
    maxCol = minmax[1]#np.max(allDesigns.iloc[:,color])

    for i in np.arange(np.size(stst)):
        
        fig1 = plt.figure(figsize=(cm2inch(4.4), cm2inch(4)))
        ax1 = plt.subplot(111)
        fig1.subplots_adjust(top=0.982,
        bottom=0.23,
        left=0.280,
        right=0.940)
        
        thisDesBool = allDesigns.iloc[:,z] == stst[i]
        thisDes = allDesigns[thisDesBool]    
        
        NiceGraph2D(ax1, xname, yname)
        
        ax1.set_xticks([0,0.5,1])
        ax1.set_xlim([-0.01,1.01])
        
        ax1.set_yscale('log')
        ax1.set_yticks([0.001,0.01,0.1,1])
        ax1.set_ylim([0.0007,1.5])
        
        if colorName == r'$Log(Energy)$':
            ax1.scatter(thisDes.iloc[:,x], thisDes.iloc[:,y], marker = 's', cmap = colormap, c = np.log10(thisDes.iloc[:, color]), 
                        vmin = minCol, vmax = maxCol, s = 5)
        else:
            ax1.scatter(thisDes.iloc[:,x], thisDes.iloc[:,y], marker = 's', cmap = colormap, c = thisDes.iloc[:, color], 
                        vmin = minCol, vmax = maxCol, s = 5)
    
        fig1.show()
        if save:
            fig1.savefig(Folder_name + '/Images/' + NameFig + '_' + '%.4f' %stst[i] +'.pdf', transparent = True)
            fig1.savefig(Folder_name + '/Images/' + NameFig + '_' + '%.4f' %stst[i] +'.png', transparent = True)
            print(Folder_name + '/Images/' + NameFig + '_' + '%.4f' %stst[i])
            
    cmaptemp = matl.cm.get_cmap(colormap)
    cmapfig, normfig = from_levels_and_colors(np.linspace(minCol, maxCol, 1000),cmaptemp(np.linspace(0, 1,1001)), extend = 'both')
    
    fig1 = plt.figure(figsize=(cm2inch(8.8), cm2inch(1.5)))
    fig1.subplots_adjust(top=0.884,
    bottom=0.116,
    left=0.059,
    right=0.951)
    ax1 = plt.subplot(111)
    cbar = plt.colorbar(matl.cm.ScalarMappable(norm=normfig, cmap=cmapfig), ax = ax1, fraction=0.99, pad=0.01, orientation='horizontal')
    cbar.set_ticks(np.linspace(minCol,maxCol,5))
    # cbar.ax.set_xticklabels((np.linspace(minCol,maxCol,5)/10).astype(int))
    cbar.set_label(colorName, fontsize = 9, color = '0.2',labelpad = 0)
    cbar.ax.tick_params(colors='0.2', pad=2)
    cbar.outline.set_edgecolor('0.2')
    cbar.outline.set_linewidth(0.4)
    ax1.remove()   
    
    if save: 
        fig1.savefig(Folder_name + '/Images/' + NameFig + '_CB.pdf', transparent = True)
        fig1.savefig(Folder_name + '/Images/' + NameFig + '_CB.png', transparent = True)
    
    return

def StableStatesCounturPlotPaper(allDesigns,x, xname, y, yname, z, color, colorName, colormap, minmax, save = False, Folder_name = '', NameFig = ''):
    
    allDesigns = allDesigns.round(8)
    
    # stst, color = getcolormap(allDesigns.iloc[:,z].values, colormap)
    # cmaptemp = matl.cm.get_cmap(colormap)
    # cmaptemp.set_over('r')
    # cmaptemp.set_under('b')
    stst = np.unique(allDesigns.iloc[:,z])
    
    minCol = minmax[0]#np.min(allDesigns.iloc[:,color])
    maxCol = minmax[1]#np.max(allDesigns.iloc[:,color])

    for i in np.arange(np.size(stst)):
        
        fig1 = plt.figure(figsize=(cm2inch(3.3), cm2inch(3.0)))
        ax1 = plt.subplot(111)
        fig1.subplots_adjust(top=0.962,
        bottom=0.26,
        left=0.305,
        right=0.935)
        
        thisDesBool = allDesigns.iloc[:,z] == stst[i]
        thisDes = allDesigns[thisDesBool]    
        
        NiceGraph2D(ax1, xname, yname)
        
        ax1.set_xticks([0,0.5,1])
        ax1.set_xlim([-0.01,1.01])
        
        ax1.set_yscale('log')
        ax1.set_yticks([0.001,0.01,0.1,1])
        ax1.set_ylim([0.0007,1.5])
            
        ax1.scatter(thisDes.iloc[:,x], thisDes.iloc[:,y], marker = 's', cmap = colormap, c = thisDes.iloc[:, color], 
                    vmin = minCol, vmax = maxCol, s = 2)
    
        fig1.show()
        if save:
            fig1.savefig(Folder_name + '/Images/' + NameFig + '_' + '%.4f' %stst[i] +'.pdf', transparent = True)
            fig1.savefig(Folder_name + '/Images/' + NameFig + '_' + '%.4f' %stst[i] +'.png', transparent = True)
            print(Folder_name + '/Images/' + NameFig + '_' + '%.4f' %stst[i])
            
    cmaptemp = matl.cm.get_cmap(colormap)
    cmapfig, normfig = from_levels_and_colors(np.linspace(minCol, maxCol, 1000),cmaptemp(np.linspace(0, 1,1001)), extend = 'both')
    
    fig1 = plt.figure(figsize=(cm2inch(1.5), cm2inch(2.8)))
    fig1.subplots_adjust(top=0.984,
    bottom=0.021,
    left=0.019,
    right=0.911)
    ax1 = plt.subplot(111)
    cbar = plt.colorbar(matl.cm.ScalarMappable(norm=normfig, cmap=cmapfig), ax = ax1, 
                        fraction=0.99, pad=0.01, orientation='vertical', aspect=15)
    cbar.set_ticks(np.linspace(minCol,maxCol,5))
    # cbar.ax.set_xticklabels(ststname.astype(int))
    cbar.set_label(colorName, fontsize = 9, color = '0.2',labelpad = 0)
    cbar.ax.tick_params(colors='0.2', pad=2)
    cbar.outline.set_edgecolor('0.2')
    cbar.outline.set_linewidth(0.4)
    ax1.remove()   
    
    if save: 
        fig1.savefig(Folder_name + '/Images/' + NameFig + '_CB.pdf', transparent = True)
        fig1.savefig(Folder_name + '/Images/' + NameFig + '_CB.png', transparent = True)
    
    return

def ColorbarPerZKappa(allMat, x, colorbar, z, save = False, Folder_name = '', NameFig = ''):
    
    allMat = allMat.round(8)
    
    z_values = np.unique(allMat.iloc[:,z])
    x_values = np.log10(np.unique(allMat.iloc[:,x]))
    
    # color = ['#66C2A5', '#FFD92F', '#E4BB05', '#B3B3B3', '#339568','#D7704A','#414596'] 
    # color = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', 
    #          '#aec7e8', '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5', '#c49c94', '#f7b6d2', 
    #          '#7f7f7f', '#bcbd22']
    # cmap = matl.cm.get_cmap('Set3',np.size(colorbar))
    # color = cmap(np.linspace(0,1,np.size(colorbar)))   '#D663A9','#e78ac3',  
    color = ['#41A584', '#E26B3C', '#6177AA',  
             '#66c2a5', '#fc8d62', '#8da0cb',             
             '#ffd92f', '#e5c494', '#b3b3b3']
    
    for j in z_values:
        fig = plt.figure(figsize=(cm2inch(8), cm2inch(6)))
        ax1 = plt.subplot(111)
        fig.subplots_adjust(top=0.988,
bottom=0.114,
left=0.154,
right=0.982)
        
        
        
        NiceGraph2D(ax1, r'$Log_{10}\kappa$', 'Sim. Amount', mincoord = [np.min(x_values),30], 
                    maxcoord = [np.max(x_values),1000], divisions = [[-3,-2,-1,0],[0,500,1000]], buffer = [0.2, 30])
        
        # ax1.set_xscale('log')
    
    
        thisMatBool = allMat.iloc[:,z] == j
        thisMat = allMat[thisMatBool]  
        
        for i in np.arange(np.size(colorbar)):
            ax1.bar(np.log10(thisMat.iloc[:,x]), thisMat.iloc[:,colorbar[i]],bottom=np.sum(thisMat.iloc[:,colorbar[0:i]], axis = 1), width = 0.201,
                    color = color[i], label = thisMat.columns[colorbar[i]], align = 'center') #
            
        CreateLegend(ax1)
        
        if save: 
            fig.savefig(Folder_name + '/Images/' + NameFig + '_' + '%.4f' %j + '.pdf', transparent = True)
            fig.savefig(Folder_name + '/Images/' + NameFig + '_' + '%.4f' %j + '.png', transparent = True)
            
    return

def ColorbarPerZ(allMat, x, colorbar, z, save = False, Folder_name = '', NameFig = ''):
    
    allMat = allMat.round(8)
    
    z_values = np.unique(allMat.iloc[:,z])
    x_values = np.unique(allMat.iloc[:,x])
    
    # color = ['#66C2A5', '#FFD92F', '#E4BB05', '#B3B3B3', '#339568','#D7704A','#414596'] 
    # color = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', 
    #          '#aec7e8', '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5', '#c49c94', '#f7b6d2', 
    #          '#7f7f7f', '#bcbd22']
    # cmap = matl.cm.get_cmap('Set3',np.size(colorbar))
    # color = cmap(np.linspace(0,1,np.size(colorbar)))   '#D663A9','#e78ac3',  
    color = ['#41A584', '#E26B3C', '#6177AA',  
             '#66c2a5', '#fc8d62', '#8da0cb',             
             '#ffd92f', '#e5c494', '#b3b3b3']
    
    for j in z_values:
        fig = plt.figure(figsize=(cm2inch(8), cm2inch(6)))
        ax1 = plt.subplot(111)
        fig.subplots_adjust(top=0.988,
bottom=0.114,
left=0.154,
right=0.982)
        
        NiceGraph2D(ax1, 'Material Size', 'Sim. Amount', mincoord = [np.min(x_values),30], 
                    maxcoord = [np.max(x_values),1000], divisions = [[2,5,10,15],[0,500,1000]], buffer = [0.7, 30])
    
        thisMatBool = allMat.iloc[:,z] == j
        thisMat = allMat[thisMatBool]  
        
        for i in np.arange(np.size(colorbar)):
            ax1.bar(thisMat.iloc[:,x], thisMat.iloc[:,colorbar[i]],bottom=np.sum(thisMat.iloc[:,colorbar[0:i]], axis = 1), width = 1,
                    color = color[i], label = thisMat.columns[colorbar[i]], align = 'center') #
            
        CreateLegend(ax1)
        
        if save: 
            fig.savefig(Folder_name + '/Images/' + NameFig + '_' + '%.4f' %j + '.pdf', transparent = True)
            fig.savefig(Folder_name + '/Images/' + NameFig + '_' + '%.4f' %j + '.png', transparent = True)
            
    return

def violinPlot(allDesigns, x, xname, y, yname, z, stst_col, save = False, Folder_name = '', NameFig = ''):
    
    allDesigns = allDesigns.round(8)
    
    color = ['#41A584', '#E26B3C', '#6177AA',  
             '#66c2a5', '#fc8d62', '#8da0cb',             
             '#ffd92f', '#e5c494', '#b3b3b3']
    
    variables = np.unique(allDesigns.iloc[:,z])
    xval = np.unique(allDesigns.iloc[:,x])
    
    for i in np.arange(np.size(variables)):
    
        fig1 = plt.figure(figsize=(cm2inch(8), cm2inch(6)))
        ax1 = plt.subplot(111)
        fig1.subplots_adjust(top=0.987,
bottom=0.111,
left=0.183,
right=0.987,
hspace=0.2,
wspace=0.2)
        
        thisDesBool = allDesigns.iloc[:,z] == variables[i]
        thisDes = allDesigns[thisDesBool]    
        
        NiceGraph2D(ax1, xname, yname)
        
        ax1.set_xlim([1.3,15.7])
        ax1.set_xticks([2,5,10,15])
        
        for k in np.arange(9):
            thisstst = thisDes[thisDes.iloc[:, stst_col] == k+1]
            for j in np.arange(np.size(xval)):
                thisstst2 = thisstst[thisstst.iloc[:, x] == xval[j]]
                if np.size(thisstst2,0)<2:
                    continue
                vio = ax1.violinplot(thisstst2.iloc[:,y], [xval[j]], widths = 0.9, showextrema = False)#, c = matl.colors.rgb2hex(color[k]), s = 8)
                for part in vio['bodies']:
                    part.set_facecolor(color[k])
                    part.set_edgecolor('0.2')
                    part.set_alpha(1)
        
        
        fig1.show()
        if save:
            fig1.savefig(Folder_name + '/Images/' + NameFig + '_' + '%.4f' %variables[i] +'.pdf', transparent = True)
            fig1.savefig(Folder_name + '/Images/' + NameFig + '_' + '%.4f' %variables[i] +'.png', transparent = True)
            print(Folder_name + '/Images/' + NameFig + '_' + '%.4f' %variables[i])
    
    return

def XSizePlot(allDesigns, x, xname, y, yname, z, stst_col, save = False, Folder_name = '', NameFig = ''):
    
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.colors import ListedColormap
    
    allDesigns = allDesigns.round(8)
    
    # color = ['#41A584', '#E26B3C', '#6177AA',  
    #          '#66c2a5', '#fc8d62', '#8da0cb',             
    #          '#ffd92f', '#e5c494', '#b3b3b3']
    
    # colors = ['#A40A20','#AFA80B','#103D73','#61B53D','#D39C47','#65368E']
    # #'#C6435A','#365C8A','#D3CD47' pink to blue to yellow, green to yellow to purple
    # #"#F49CAC", "#860C21", "#ACE693", "#2C7B0B" light pink to dark pinto, light green to dark green
    # nodes = [0.0, 0.2222, 0.4444, 0.5556, 0.7778, 1.0]
    # cmap2 = LinearSegmentedColormap.from_list("mycmap", list(zip(nodes, colors)))
    
    # top = matl.cm.get_cmap('summer', 128)
    # bottom = matl.cm.get_cmap('winter', 128)
    
    # newcolors = np.vstack((top(np.linspace(0, 1, 128)),
    #                        bottom(np.linspace(0, 1, 128))))
    # cmap2 = ListedColormap(newcolors, name='poop')
    cmap2 = 'gist_rainbow'
    
    minLine = 0.5 #np.min(allDesigns.iloc[:,stst_col].values)
    maxLine = 18.5 #np.max(allDesigns.iloc[:,stst_col].values)
    
    variables = np.unique(allDesigns.iloc[:,z])
    
    for i in np.arange(np.size(variables)):
    
        # fig1 = plt.figure(figsize=(cm2inch(4.3), cm2inch(3.1)))
        fig1 = plt.figure(figsize=(cm2inch(8), cm2inch(6)))
        ax1 = plt.subplot(111)
        fig1.subplots_adjust(top=0.987,
bottom=0.111,
left=0.183,
right=0.987,
hspace=0.2,
wspace=0.2)
        
        thisDesBool = allDesigns.iloc[:,z] == variables[i]
        thisDes = allDesigns[thisDesBool]    
        
        NiceGraph2D(ax1, xname, yname)
        
        ax1.set_xlim([1.3,15.7])
        ax1.set_xticks([2,5,10,15])
        
        if xname == r'$\kappa$':
            ax1.set_xscale('log')
            ax1.set_xticks([0.001,0.01,0.1,1])
            ax1.set_xlim([0.0007,1.5])
            
        if (yname == r'$E_{norm}$'):
            ax1.set_yscale('log')
            
        ax1.scatter(thisDes.iloc[:,x]+np.random.rand(np.size(thisDes.iloc[:,x]))*0.5-0.25, thisDes.iloc[:,y],
                    c = thisDes.iloc[:,stst_col],
                    cmap = cmap2, vmin = minLine, vmax = maxLine, s = 15)
        
        fig1.show()
        if save:
            fig1.savefig(Folder_name + '/Images/' + NameFig + '_' + '%.4f' %variables[i] +'.pdf', transparent = True)
            fig1.savefig(Folder_name + '/Images/' + NameFig + '_' + '%.4f' %variables[i] +'.png', transparent = True)
            print(Folder_name + '/Images/' + NameFig + '_' + '%.4f' %variables[i])
    
    fig1 = plt.figure(figsize=(cm2inch(1.5), cm2inch(2.8)))
    fig1.subplots_adjust(top=0.984,
    bottom=0.021,
    left=0.019,
    right=0.911)
    ax1 = plt.subplot(111)
    
    cmaptemp = matl.cm.get_cmap(cmap2)
    # cmaptemp = cmap2
    cmapfig, normfig = from_levels_and_colors(np.linspace(minLine, maxLine, 19),cmaptemp(np.linspace(0, 1,20)), extend = 'both')
    
    cbar = plt.colorbar(matl.cm.ScalarMappable(norm=normfig, cmap=cmapfig), ax = ax1, 
                        fraction=0.99, pad=0.01, orientation='vertical', aspect=15)
    cbar.set_ticks(np.linspace(1,19,7))
    # cbar.ax.set_xticklabels(ststname.astype(int))
    cbar.set_label('Dislocation lines', fontsize = 9, color = '0.2',labelpad = 0)
    cbar.ax.tick_params(colors='0.2', pad=2)
    cbar.outline.set_edgecolor('0.2')
    cbar.outline.set_linewidth(0.4)
    ax1.remove() 
    
    return

def violin_scatter(allDesigns, x, xname, y, yname, stst_col, save = False, Folder_name = '', NameFig = ''):
    
    allDesigns = allDesigns.round(8)
    
    color = ['#41A584', '#E26B3C', '#6177AA',  
             '#66c2a5', '#fc8d62', '#8da0cb']
    
    xval = np.unique(allDesigns.iloc[:,x])
    
    for k in np.arange(3):
    
        fig1 = plt.figure(figsize=(cm2inch(8), cm2inch(6)))
        ax1 = plt.subplot(111)
        fig1.subplots_adjust(top=0.987,
bottom=0.111,
left=0.183,
right=0.987,
hspace=0.2,
wspace=0.2)
        
        NiceGraph2D(ax1, xname, yname)
        
        if (yname == r'$E_{norm}$'):
            ax1.set_yscale('log')
        
        ax1.set_xlim([1.3,15.7])
        ax1.set_xticks([2,5,10,15])
        
        thisstst = allDesigns[allDesigns.iloc[:, stst_col] == k+1+3]
        for j in np.arange(np.size(xval)):
            thisstst2 = thisstst[thisstst.iloc[:, x] == xval[j]]
            if np.size(thisstst2,0)<2:
                continue
            vio = ax1.violinplot(thisstst2.iloc[:,y], [xval[j]], widths = 0.9, showextrema = False)#, c = matl.colors.rgb2hex(color[k]), s = 8)
            for part in vio['bodies']:
                part.set_facecolor(color[k+3])
                part.set_edgecolor('1')
                part.set_alpha(0.7)
          
        thisstst = allDesigns[allDesigns.iloc[:, stst_col] == k+1]
        ax1.scatter(thisstst.iloc[:,x], thisstst.iloc[:,y], c = color[k], s = 8)
        
        
        fig1.show()
        if save:
            fig1.savefig(Folder_name + '/Images/' + NameFig + '_' + '%.4f' %k +'.pdf', transparent = True)
            fig1.savefig(Folder_name + '/Images/' + NameFig + '_' + '%.4f' %k +'.png', transparent = True)
            print(Folder_name + '/Images/' + NameFig + '_' + '%.4f' %k)
    
    return
    
def violin_scatterKappa(allDesigns, x, xname, y, yname, stst_col, save = False, Folder_name = '', NameFig = ''):
    
    allDesigns = allDesigns.round(8)
    
    color = ['#41A584', '#E26B3C', '#6177AA',  
             '#66c2a5', '#fc8d62', '#8da0cb']
    
    xval = np.log10(np.unique(allDesigns.iloc[:,x]))
    
    for k in np.arange(3):
    
        fig1 = plt.figure(figsize=(cm2inch(8), cm2inch(6)))
        ax1 = plt.subplot(111)
        fig1.subplots_adjust(top=0.987,
bottom=0.111,
left=0.183,
right=0.987,
hspace=0.2,
wspace=0.2)
        
        NiceGraph2D(ax1, xname, yname)

        ax1.set_xlim([-3.2,0.2])
        ax1.set_xticks([-3,-2,-1,0])
        
        thisstst = allDesigns[allDesigns.iloc[:, stst_col] == k+1+3]
        for j in np.arange(np.size(xval)):
            thisstst2 = thisstst[thisstst.iloc[:, x] == 10**xval[j]]
            if np.size(thisstst2,0)<2:
                continue
            vio = ax1.violinplot(thisstst2.iloc[:,y], [xval[j]], widths = 0.2, showextrema = False)#, c = matl.colors.rgb2hex(color[k]), s = 8)
            for part in vio['bodies']:
                part.set_facecolor(color[k+3])
                part.set_edgecolor('1')
                part.set_alpha(0.7)
          
        thisstst = allDesigns[allDesigns.iloc[:, stst_col] == k+1]
        ax1.scatter(np.log10(thisstst.iloc[:,x]), thisstst.iloc[:,y], c = color[k], s = 8)
        
        
        fig1.show()
        if save:
            fig1.savefig(Folder_name + '/Images/' + NameFig + '_' + '%.4f' %k +'.pdf', transparent = True)
            fig1.savefig(Folder_name + '/Images/' + NameFig + '_' + '%.4f' %k +'.png', transparent = True)
            print(Folder_name + '/Images/' + NameFig + '_' + '%.4f' %k)
    
    return
    
def XKappawSingleVertPlot(allDesigns, allDesignsSV, x1, x2, xname, y1, y2, yname, stst_col, save = False, Folder_name = '', NameFig = ''):
    
    allDesigns = allDesigns.round(8)
    allDesignsSV = allDesignsSV.round(8)
    
    color = ['#41A584', '#E26B3C', '#6177AA',  
              '#66c2a5', '#fc8d62', '#8da0cb',             
              '#ffd92f', '#e5c494', '#b3b3b3']
    colormap = 'gist_rainbow'
    
    minLine = np.min(allDesigns.iloc[:,stst_col].values) # 0.5 #
    maxLine = np.max(allDesigns.iloc[:,stst_col].values) #18.5 #
    if maxLine < minLine:
        maxLine = minLine + 0.1
    
    SingleVertexStSt = np.array([0, 5, 9])
    
    for i in np.arange(3):
    
        # fig1 = plt.figure(figsize=(cm2inch(4.3), cm2inch(3.1)))
        fig1 = plt.figure(figsize=(cm2inch(8), cm2inch(6)))
        ax1 = plt.subplot(111)
        fig1.subplots_adjust(top=0.987,
bottom=0.111,
left=0.183,
right=0.987,
hspace=0.2,
wspace=0.2)
        
        thisPureBool = allDesigns['StableStateMat'] == i+1
        thisPure = allDesigns[thisPureBool]    
        
        thisNonPureBool = allDesigns['StableStateMat'] == i+1+3
        thisNonPure = allDesigns[thisNonPureBool] 
        
        thisSVBool = allDesignsSV['StableStateAll'] == SingleVertexStSt[i]
        thisSV = allDesignsSV[thisSVBool]
        
        NiceGraph2D(ax1, xname, yname)
        
        if (yname == r'$E_{norm}$'):
            ax1.set_yscale('log')
        
        ax1.set_xscale('log')
        ax1.set_xticks([0.001,0.01,0.1,1])
        ax1.set_xlim([0.0007,1.5])
        
        #+np.random.rand(np.size(thisDes.iloc[:,x]))*0.01
            
        ax1.scatter(thisPure.iloc[:,x1]*10**(np.random.rand(np.size(thisPure.iloc[:,x1]))*0.1), 
                    thisPure.iloc[:,y1], c = color[i], s = 8)
        
        ax1.scatter(thisNonPure.iloc[:,x1]*10**(np.random.rand(np.size(thisNonPure.iloc[:,x1]))*0.1), 
                    thisNonPure.iloc[:,y1], c = thisNonPure.iloc[:,stst_col],
                    cmap = colormap, vmin = minLine, vmax = maxLine, s = 15)
        
        ax1.plot(thisSV.iloc[:,x2], thisSV.iloc[:,y2])
        
        if i == 0:
            thisSVBool = allDesignsSV['StableStateAll'] == 1
            thisSV = allDesignsSV[thisSVBool]
            ax1.plot(thisSV.iloc[:,x2], thisSV.iloc[:,y2])
        
        fig1.show()
        if save:
            fig1.savefig(Folder_name + '/Images/' + NameFig + '_' + '%.4f' %(i+1) +'.pdf', transparent = True)
            fig1.savefig(Folder_name + '/Images/' + NameFig + '_' + '%.4f' %(i+1) +'.png', transparent = True)
            print(Folder_name + '/Images/' + NameFig + '_' + '%.4f' %(i+1))
    
    fig1 = plt.figure(figsize=(cm2inch(1.5), cm2inch(2.8)))
    fig1.subplots_adjust(top=0.984,
    bottom=0.021,
    left=0.019,
    right=0.911)
    ax1 = plt.subplot(111)
    
    cmaptemp = matl.cm.get_cmap(colormap)
    cmapfig, normfig = from_levels_and_colors(np.linspace(minLine, maxLine, 1000),cmaptemp(np.linspace(0, 1,1001)), extend = 'both')
    
    cbar = plt.colorbar(matl.cm.ScalarMappable(norm=normfig, cmap=cmapfig), ax = ax1, 
                        fraction=0.99, pad=0.01, orientation='vertical', aspect=15)
    # cbar.set_ticks(np.linspace(1,19,7))
    cbar.set_ticks(np.linspace(minLine, maxLine, 4))
    # cbar.ax.set_xticklabels(ststname.astype(int))
    cbar.set_label('Dislocation lines', fontsize = 9, color = '0.2',labelpad = 0)
    cbar.ax.tick_params(colors='0.2', pad=2)
    cbar.outline.set_edgecolor('0.2')
    cbar.outline.set_linewidth(0.4)
    ax1.remove() 
    
    if save:
        fig1.savefig(Folder_name + '/Images/' + NameFig + '_CB.pdf', transparent = True)
        fig1.savefig(Folder_name + '/Images/' + NameFig + '_CB.png', transparent = True)
    
    return

def DefectsApperance(allMat, x, xname, save = False, Folder_name = '', NameFig = ''):
    
    allMat = allMat.round(8)
    
    defTypes = [[11,12,13,18], [12,14,15,16,17,18],[13,17,18]]
    defNames = [['GB Twin', 'Interface ds-miura', 'Interface ds-fold', 'Interface all'], 
                ['Interface ds-miura', 'GB Twin', 'GB Slip', 'GB Rotation', 'Interface miura-fold', 'Interface all'], 
                ['Interface ds-fold','Interface miura-fold', 'Interface all']]
    defColors = [[0,1,2,3],[1,4,5,6,7,3],[2,7,3]]
    
    color = ['#41A584', '#E26B3C', '#6177AA',  
              '#66c2a5', '#fc8d62', '#8da0cb',             
              '#ffd92f', '#e5c494', '#b3b3b3']
    
    for i in np.arange(3):
    
        # fig1 = plt.figure(figsize=(cm2inch(4.3), cm2inch(3.1)))
        fig1 = plt.figure(figsize=(cm2inch(8), cm2inch(6)))
        ax1 = plt.subplot(111)
        fig1.subplots_adjust(top=0.987,
bottom=0.111,
left=0.183,
right=0.987,
hspace=0.2,
wspace=0.2)
        
        thisStateBool = allMat['StableStateMat'] == i+1+3
        thisState = allMat[thisStateBool] 
        
        NiceGraph2D(ax1, xname, 'Sim. Amount')
        
        if (xname == r'$\kappa$'):
            ax1.set_xscale('log')
            ax1.set_xticks([0.001,0.01,0.1,1])
            ax1.set_xlim([0.0007,1.5])
        
        for j in np.arange(np.size(defTypes[i])):
            
            ax1.scatter(thisState.iloc[:,x], 
                        thisState.iloc[:,defTypes[i][j]], c = color[defColors[i][j]], s = 8, label = defNames[i][j])
        
        CreateLegend(ax1)  
        
        fig1.show()
        if save:
            fig1.savefig(Folder_name + '/Images/' + NameFig + '_' + '%.4f' %(i+1) +'.pdf', transparent = True)
            fig1.savefig(Folder_name + '/Images/' + NameFig + '_' + '%.4f' %(i+1) +'.png', transparent = True)
            print(Folder_name + '/Images/' + NameFig + '_' + '%.4f' %(i+1))
            
    return

def GrainSize(allMat, x, xname, save = False, Folder_name = '', NameFig = ''):
    
    allMat = allMat.round(8)
    
    # fig1 = plt.figure(figsize=(cm2inch(4.3), cm2inch(3.1)))
    fig1 = plt.figure(figsize=(cm2inch(8), cm2inch(6)))
    ax1 = plt.subplot(111)
    fig1.subplots_adjust(top=0.987,
bottom=0.111,
left=0.183,
right=0.987,
hspace=0.2,
wspace=0.2)
    
    NiceGraph2D(ax1, xname, 'Vertex/Grain')
    
    MatNames = ['Dome-Saddle','All-Miura-ori', 'All-Fold']
    color = ['#66c2a5', '#fc8d62', '#8da0cb']
    
    if (xname == r'$\kappa$'):
        ax1.set_xscale('log')
        ax1.set_xticks([0.001,0.01,0.1,1])
        ax1.set_xlim([0.0007,1.5])
        
        
    for i in np.arange(3):
    
        ax1.scatter(allMat.iloc[:,x], 
                    allMat.iloc[:,12+i], c = color[i], s = 8, label = MatNames[i])
        
    CreateLegend(ax1)  
    
    fig1.show()
    if save:
        fig1.savefig(Folder_name + '/Images/' + NameFig + '.pdf', transparent = True)
        fig1.savefig(Folder_name + '/Images/' + NameFig + '.png', transparent = True)
        print(Folder_name + '/Images/' + NameFig )
            
    return



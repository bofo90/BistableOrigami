# -*- coding: utf-8 -*-
"""
Created on Tue May 21 18:04:51 2019

@author: iniguez
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as matl
from mpl_toolkits import mplot3d
import xarray as xr
import configparser
import os.path
import ternary #from https://github.com/marcharper/python-ternary
import itertools
import copy

def cm2inch(value):
    return value/2.54

def orderAngles(angles, steps, simulations):
    
    finalAngles = np.empty((0,np.size(angles,1)))
    for hinge in np.arange(simulations):
        if np.sum(np.sign(np.around(angles[steps*(hinge+1)-1,:], decimals = 4))) < 0:
            angles[steps*(hinge+1)-1,:] = angles[steps*(hinge+1)-1,:]*(-1)        
        sortAllAngIndex = np.lexsort((angles[steps*(hinge+1)-1,:],angles[steps*hinge,:]))#[0,1,2,3]
        finalAngles = np.append(finalAngles, [angles[steps*(hinge+1)-1,sortAllAngIndex]], axis = 0)
        
    return finalAngles
        
        
def countStableStates(finalAngles, distance, method, plot = False):
    
    import scipy.cluster.hierarchy as hierarch
    from scipy.spatial.distance import pdist
                                       
    Z = hierarch.linkage(finalAngles, method)
    inverse = hierarch.fcluster(Z, distance, criterion='distance')
    c = hierarch.cophenet(Z, pdist(finalAngles))

    if plot:
        
        print('this is the cophenet of the hierarchical linkage', c[0])
        
        plt.figure(0,figsize=(25, 10))
        plt.title('Hierarchical Clustering Dendrogram')
        plt.xlabel('sample index')
        plt.ylabel('distance')
        hierarch.dendrogram(
            Z,
            truncate_mode='lastp',  # show only the last p merged clusters
            p=50,  # show only the last p merged clusters
            leaf_rotation=90.,  # rotates the x axis labels
            leaf_font_size=8.,  # font size for the x axis labels
            show_contracted=True,
        )
        plt.show()

    
    return inverse

def NiceGraph2D(axes, nameX, nameY, mincoord = [np.NaN, np.NaN], maxcoord = [np.NaN, np.NaN], divisions = [np.NaN, np.NaN],
                buffer = [0.0, 0.0, 0.0]):
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
    axes.tick_params(axis='x',colors=gray, width = 0.4)
    axes.spines['bottom'].set_color(gray)
    axes.spines['top'].set_color(gray)
    axes.yaxis.label.set_color(gray)
    axes.tick_params(axis='y', colors=gray, width = 0.4)
    axes.spines['left'].set_color(gray)
    axes.spines['right'].set_color(gray)
    axes.tick_params(pad = 2)
    
    for axis in ['top','bottom','left','right']:
        axes.spines[axis].set_linewidth(0.4)
        
    return

def NiceGraph2Dlog(axes, nameX, nameY, mincoord = [np.NaN, np.NaN], maxcoord = [np.NaN, np.NaN], divisions = [np.NaN, np.NaN],
                buffer = [0.0, 0.0, 0.0]):
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

#%%
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

#%%
def ReadMetadata(file):
    
    if os.path.isfile(file):
        metadata = configparser.RawConfigParser(allow_no_value=True)
        metadata.read(file)
    
        restang = float(metadata.get('options','restang'))
        designang = np.array(metadata.get('options','angDesign').split(),dtype = float)
    else:
        raise FileNotFoundError('No metafile found at the given directory. Changes to the script to put manually the variables are needed\n') 

    return restang, designang

kappas = np.logspace(-3,1,81)
#kappas = np.logspace(-3,1,13)#23

Folder_name = "Results/SingleVertex4/sqp/energy/15-Nov-2019_KandTheta0/DoubleSym"
#Folder_name = "Results/SingleVertex3/sqp/energy/16-Oct-2019_KandTheta0/Angles_120_120"
file_name1 = "/EnergyData.csv" 
file_name2 = "/Hinges.csv"
file_name3 = "/PosStad.csv"
file_name4 = "/Angles.csv"

allDesigns = pd.DataFrame()
allKappasAnalysis = pd.DataFrame()

if not os.path.isdir(Folder_name + '/Images/'):
    os.mkdir(Folder_name + '/Images/')
    os.mkdir(Folder_name + '/Images/SingleDes/')

for subdir in os.listdir(Folder_name):
    
    if subdir == 'Images':
        continue
    
    plt.close('all')
    prevFolder_name = Folder_name + '/' + subdir
    
    allData = pd.DataFrame()
    allAngles = np.empty((0,4))
#    allAngles = np.empty((0,3))
    
    for k in kappas:
        folder_name = prevFolder_name + "/kh%.5f_kta1000.00_ke1.00_kf100.00" %k
    
        dataEnergy = pd.read_csv(folder_name+file_name1)
        simLen = np.size(dataEnergy['Hinge Number'])
        
        dataEnergy['TotalEnergy'] = dataEnergy['EdgeEnergy']+dataEnergy['DiagonalEnergy']+dataEnergy['HingeEnergy']
        dataEnergy['HingeEnergy'] = dataEnergy['HingeEnergy']
        dataEnergy['kappa'] = np.ones(simLen)*k
        
        dataHinges = pd.read_csv(folder_name+file_name2)
        dataEnergy['Curvature'] = dataHinges['Curvature']
    
        dataAngles = np.loadtxt(folder_name+file_name4,skiprows=1, delimiter = ',', dtype = np.float64)
        dataAngles = np.delete(dataAngles, 0, 1)
        dataAnglesOrd = orderAngles(dataAngles, 2, simLen)
        dataEnergy['StableStates'] = np.zeros((simLen,1))
        dataEnergy['StableStates'] = countStableStates(dataAnglesOrd, 10, 'ward')
#        dataEnergy['StableStates'] = countStableStates(dataAnglesOrd, 0.5, 'centroid')
        dataEnergy[['ang1','ang2','ang3','ang4']] = pd.DataFrame(dataAnglesOrd)
#        dataEnergy[['ang1','ang2','ang3']] = pd.DataFrame(dataAnglesOrd)
        
        allAngles = np.append(allAngles, dataAngles, axis = 0)
        allData = allData.append(dataEnergy)
        
    restang, designang = ReadMetadata(folder_name+'/metadata.txt') 
     
    TotSimul = allData.shape[0]
    
    print(allData['Flags'].value_counts())
    exfl = allData.Flags.values.reshape(TotSimul,1)
    flagmask = (exfl !=1) & (exfl !=2)
    flagmask = ~flagmask.any(axis= 1)
    allData['Mask'] = flagmask
    
    allData.set_index(['kappa','Hinge Number'], inplace = True)
    allData.sort_index(inplace=True)
    
    kappasStSt = allData.groupby('kappa').apply(lambda _df: _df.groupby('StableStates')[['EdgeEnergy', 'HingeEnergy', 'TotalEnergy','ang1','ang2','ang3','ang4', 'Curvature']].mean())
#    kappasStSt = allData.groupby('kappa').apply(lambda _df: _df.groupby('StableStates')[['EdgeEnergy', 'HingeEnergy', 'TotalEnergy','ang1','ang2','ang3', 'Curvature']].mean())
    ####Only select stable states that have more than 10% appearance in each simulation
    kappasStSt['amountStSt'] = allData.groupby('kappa').apply(lambda _df: _df.groupby('StableStates')[['DiagonalEnergy']].count())
    more10percent = kappasStSt['amountStSt']>simLen*0.1
    kappasStSt = kappasStSt[more10percent]
    
    kappasStSt = kappasStSt.reset_index(level=0)
    kappasStSt = kappasStSt.reset_index(level=0)#, drop=True)
#    kappasStSt['StableState'] = countStableStates(kappasStSt[['ang1','ang2','ang3','ang4']],1, 'centroid')
#    kappasStSt['StableState'] = countStableStates(kappasStSt[['ang1','ang2','ang3']],0.45, 'centroid')
    
    selection = allData.groupby('kappa').apply(lambda _df: _df.groupby('StableStates').apply(lambda _df2: _df2.sample(1, random_state = 0)))
    selection = selection[more10percent]
    selection = selection.reset_index(level = [0,1], drop = True)
    selection = selection.reset_index(level = [0,1])
    
    kappasStSt = kappasStSt.join(selection['Hinge Number'])
    kappasStSt['LogKappas'] = np.log10(kappasStSt.kappa)
    
    kappasStSt['restang'] = np.ones(np.shape(kappasStSt)[0])*restang
    
    kappasnumStSt = kappasStSt.groupby('kappa')['StableStates'].nunique()
    kappasnumStSt = kappasnumStSt.reset_index()
    kappasnumStSt['restang'] = np.ones(np.shape(kappasnumStSt)[0])*restang

    allKappasAnalysis = allKappasAnalysis.append(kappasnumStSt)
    
    data = kappasStSt.to_xarray()
    
    allDesigns = allDesigns.append(kappasStSt)    
      
allDesigns = allDesigns.reset_index(level=0, drop =True)

#%%
#fig1 = plt.figure(figsize=(cm2inch(8), cm2inch(6)))
#ax1 =plt.subplot(111)
#fig1.subplots_adjust(top=0.967,
#bottom=0.185,
#left=0.215,
#right=0.917)
#
#NiceGraph2D(ax1, r'$Log(\kappa)$', r'$\theta_0/\pi$', mincoord = [np.log10(kappas[0]), allDesigns.restang.min()], 
#            maxcoord = [np.log10(kappas[-1]), allDesigns.restang.max()], divisions = [5,4], buffer = [0.1,0.05])
#ax1.yaxis.set_major_formatter(matl.ticker.FormatStrFormatter('%.2f'))
#
#sep1 = (np.log10(kappas.max())-np.log10(kappas.min()))/np.size(kappas)/2
#sep2 = (allDesigns.restang.max()-allDesigns.restang.min())/(np.size(allDesigns.restang.unique())-1)/2
#
#cs1 = ax1.imshow(allDesigns.TotalEnergy.values.reshape(10,13), extent=[np.log10(kappas[0])-sep1,np.log10(kappas[-1])+sep1,allDesigns.restang.min()-sep2,allDesigns.restang.max()+sep2], 
#                     cmap = matl.cm.nipy_spectral,vmax = allDesigns.TotalEnergy.max(), vmin = 0, aspect = 'auto', origin = 'lower') #nipy_spectral,
##cs1 = ax1.imshow(allDesigns.Curvature.values.reshape(10,13), extent=[np.log10(kappas[0])-sep1,np.log10(kappas[-1])+sep1,allDesigns.restang.min()-sep2,allDesigns.restang.max()+sep2], 
##                     cmap = matl.cm.nipy_spectral,vmax = allDesigns.Curvature.max(), vmin = 0, aspect = 'auto', origin = 'lower') #nipy_spectral,
##cs1 = ax1.imshow(allDesigns.ang1.values.reshape(10,13), extent=[np.log10(kappas[0])-sep1,np.log10(kappas[-1])+sep1,allDesigns.restang.min()-sep2,allDesigns.restang.max()+sep2], 
##                     cmap = matl.cm.nipy_spectral,vmax = allDesigns.ang1.max(), aspect = 'auto', origin = 'lower') #nipy_spectral,
#
#cbar = plt.colorbar(cs1, pad=0.01, format="%.2f")# fraction=0.99,extend = 'max',, orientation='horizontal'
##cbar.set_ticks(np.linspace(0, allDesigns.Curvature.max(), 5))
##cbar.set_label(r'$C$', fontsize = 9, color = '0.2',labelpad = 0)
#cbar.set_ticks(np.linspace(0, allDesigns.TotalEnergy.max(), 5))
#cbar.set_label(r'$E_\mathrm{tot}$', fontsize = 9, color = '0.2',labelpad = 0)
#cbar.ax.tick_params(colors='0.2', pad=2, width=0.4)
#cbar.outline.set_edgecolor('0.2')
#cbar.outline.set_linewidth(0.4)
#
#fig1.show()
##fig1.savefig(Folder_name + '/Images/' + 'Theta0Kappa_Curvature' + '.pdf', transparent = True)
##fig1.savefig(Folder_name + '/Images/' + 'Theta0Kappa_Curvature' + '.png', transparent = True)
#fig1.savefig(Folder_name + '/Images/' + 'Theta0Kappa_Energy' + '.pdf', transparent = True)
#fig1.savefig(Folder_name + '/Images/' + 'Theta0Kappa_Energy' + '.png', transparent = True)

#%%
allDesigns['StableStateAll'] = countStableStates(allDesigns[['ang1','ang2','ang3','ang4']], 1, 'centroid', True)
#allDesigns['StableStateAll'] = countStableStates(allDesigns[['ang1','ang2','ang3']], 0.5, 'centroid')
stst = np.unique(allDesigns['StableStateAll'])
############################################# TO make SS consistence between plots (its done manualy)
newSS = np.array([4,4,3,3,1,1,2,1,2])
invmask_copy = copy.deepcopy(allDesigns.StableStateAll)

for old,new in zip(stst, newSS):
    allDesigns.StableStateAll[invmask_copy == old] = new

delSS = np.array([4,5,6,7,8,9])
stst = np.delete(stst, delSS)
#####################################################################################################

cmap2 = matl.cm.get_cmap('Set2',np.size(kappas))
colors = cmap2(np.linspace(0,1,np.size(kappas)))
#
#plt.scatter(np.log10(allKappasAnalysis['kappa']),allKappasAnalysis['restang'],c = colors[allKappasAnalysis['StableStates'].values])
#
#for state in stst:
#    plt.close('all')
#    
#    fig2 = plt.figure(figsize=(cm2inch(8.6), cm2inch(4.5)))
#    ax1 = plt.subplot(121)
#    ax2 = plt.subplot(122)
#    fig2.subplots_adjust(top=0.885,
#        bottom=0.24,
#        left=0.155,
#        right=0.985,
#        hspace=0.215,
#        wspace=0.5)
#    NiceGraph2D(ax1, 'Kappa', 'RestAngle' , mincoord=[0,0], maxcoord=[180,180], divisions=[7,7], buffer = [5,5])
#    
#    thisstate = allDesigns[allDesigns['StableStateAll'] == state]
#    
#    for i in np.arange(np.size(kappas)):
#        thiskappa = thisstate[thisstate['kappa'] == kappas[i]]
#       
#        if not thiskappa.empty:
#            ax1.scatter(thiskappa['desang1'].values,thiskappa['desang3'].values, c = [colors[i]], marker = 's', s = 4, edgecolor = 'None', zorder = 5-i)
#            ax1.scatter(-thiskappa['desang1'].values+180,-thiskappa['desang3'].values+180, zorder = 5-i, marker = 's',  c = [colors[i]], s = 4, edgecolor = 'None')
#            ax1.scatter(thiskappa['desang3'].values,thiskappa['desang1'].values, c = [colors[i]], zorder = 5-i, marker = 's',  s = 4, edgecolor = 'None')
#            ax1.scatter(-thiskappa['desang3'].values+180,-thiskappa['desang1'].values+180, zorder = 5-i, marker = 's',  c = [colors[i]], s = 4, edgecolor = 'None')
#            ax2.scatter(thiskappa['desang1'].values,thiskappa['desang3'].values, c = [colors[i]], s = 4, edgecolor = 'None', marker = 's', label = str(kappas[i]))
#            ax2.scatter(-thiskappa['desang1'].values+180,-thiskappa['desang3'].values+180, c = [colors[i]], marker = 's', s = 4, edgecolor = 'None')
#            ax2.scatter(thiskappa['desang3'].values,thiskappa['desang1'].values, c = [colors[i]], marker = 's', s = 4, edgecolor = 'None')
#            ax2.scatter(-thiskappa['desang3'].values+180,-thiskappa['desang1'].values+180, c = [colors[i]], marker = 's', s = 4, edgecolor = 'None')
#        else:
#            ax2.scatter([],[],c = [colors[i]], s = 4, edgecolor = 'None', marker = 's', label = str(kappas[i]) )
#            
#            
#    leg = ax2.legend(loc = 5, fontsize = 7, framealpha = 0.8, edgecolor = 'inherit', fancybox = False)
#    plt.setp(leg.get_texts(), color='0.2')
#    leg.get_frame().set_linewidth(0.4)
#
#
#    fig2.show()
#    fig2.savefig(Folder_name + '/Images/' + 'DesignSpaceStStLastK' + str(state) + '.pdf', transparent = True)
#    fig2.savefig(Folder_name + '/Images/' + 'DesignSpaceStStLastK' + str(state) + '.png', transparent = True)
##    fig2.savefig(Folder_name + '/Images/' + 'DesignSpaceStStFirstK' + str(state) + '.pdf', transparent = True)
##    fig2.savefig(Folder_name + '/Images/' + 'DesignSpaceStStFirstK' + str(state) + '.png', transparent = True)

    #%%
restangles = allKappasAnalysis.restang.drop_duplicates().values

#maxTotEn = allDesigns['TotalEnergy'].max()#quantile(0.95)
#maxHinEn = allDesigns['HingeEnergy'].max()#quantile(0.95)
#maxStrEn = allDesigns['EdgeEnergy'].max()#quantile(0.95)
#maxCurv = allDesigns['Curvature'].max()#quantile(0.95)
#minCurv = allDesigns['Curvature'].min()

allDesigns['StableStatesAng'] = 0

maxTotEn = allDesigns['TotalEnergy'].max()#quantile(0.95)
maxCurv = allDesigns['Curvature'].max()#quantile(0.95)
minCurv = allDesigns['Curvature'].min()
    
fig1 = plt.figure(figsize=(cm2inch(8), cm2inch(6)))
ax1 = plt.subplot(111)
fig1.subplots_adjust(top=0.987,
bottom=0.18,
left=0.227,
right=0.982)

markers = np.array(['o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X'])
markers = ['o','^','s']
cmap = matl.cm.get_cmap('Set2',np.size(stst))
colors = cmap(np.linspace(0,1,np.size(stst)))
#    NiceGraph2D(ax1,  r'$Log(\kappa)$', r'$E_\mathrm{tot}$', mincoord=[np.log10(kappas[0]),0], maxcoord=[np.log10(kappas[-1]),maxTotEn], divisions=[5, 3], buffer=[0.1, 0.01])
#NiceGraph2D(ax1, 'Kappa', 'Hinge Energy', mincoord=[kappas[0],0], maxcoord=[kappas[-1],maxHinEn], divisions=[5, 3], buffer=[0, 0.001])
#NiceGraph2D(ax1, 'Kappa', 'Edge Energy', mincoord=[kappas[0],0], maxcoord=[kappas[-1],maxStrEn], divisions=[5, 3], buffer=[0, 0.005])
NiceGraph2D(ax1, r'$Log(\kappa)$', r'$C$', mincoord=[np.log10(kappas[0]),-3], maxcoord=[np.log10(kappas[-1]),3], divisions=[5, 5], buffer=[0.1, 1])
ax1.yaxis.set_major_formatter(matl.ticker.FormatStrFormatter('%.2f'))

for i in np.arange(np.size(restangles)):
    thisangBool = allDesigns['restang'] == restangles[i]
#    angles = allDesigns[['ang1','ang2','ang3','ang4']]
#    allDesigns.loc[thisangBool,'StableStatesAng'] = countStableStates(angles[thisangBool], 2, 'ward', False)
    
    thisang = allDesigns[thisangBool]
   
#    plt.close('all')

#    ax1.scatter(np.log10(thisang['kappa']), thisang['TotalEnergy'], c = colors[thisang['StableStateAll']-1])
#    ax1.scatter(thisang['kappa'], thisang['EdgeEnergy'], c = colors[thisang['StableStateAll']-1])
#    ax1.scatter(thisstate['kappa'], thisstate['EdgeEnergy'], c = [j])
    ax1.scatter(np.log10(thisang['kappa']), thisang['Curvature'], c = colors[thisang['StableStateAll'].values.astype('int')-1],
                s = 5, marker = markers[i])

#leg = ax1.legend(loc = 2, fontsize = 7, framealpha = 0.8, edgecolor = 'inherit', fancybox = False) 
##           borderpad = 0.3, labelspacing = 0.1, handlelength = 0.4, handletextpad = 0.4)
#plt.setp(leg.get_texts(), color='0.2')
#leg.get_frame().set_linewidth(0.4)
    
fig1.show()
#    fig1.savefig(Folder_name + '/Images/Energy_Restang_' + str(i.astype(float))+'.pdf', transparent = True)
#    fig1.savefig(Folder_name + '/Images/Energy_Restang_' + str(i.astype(float))+ '.png', transparent = True)
fig1.savefig(Folder_name + '/Images/All_Restang_' + str(i.astype(float))+'.pdf', transparent = True)
fig1.savefig(Folder_name + '/Images/All_Restang_' + str(i.astype(float))+ '.png', transparent = True)
#%%
    
    
for d in designs:
    plt.close('all')
    
    thisdesign = allDesigns[(allDesigns[['desang1','desang2','desang3','desang4']] == d).all(axis = 1)]
    
    fig1,axes  = plt.subplots(1,4, figsize=(cm2inch(20), cm2inch(5.5)))
    fig1.subplots_adjust(top=0.98,
    bottom=0.215,
    left=0.085,
    right=0.995,
    hspace=0.25,
    wspace=0.495)
    
    cmap = matl.cm.get_cmap('Set2',np.size(stst))
    colors = cmap(np.linspace(0,1,np.size(stst)))
    
    for i, ax in enumerate(axes.flat):
        ax.set_xscale('log')
        NiceGraph2Dlog(ax, 'Kappa', 'Angle '+ str(i+1) +' [rad/pi]', mincoord=[kappas[0],-0.5], maxcoord=[kappas[-1],0.5], divisions=[3, 3], buffer=[2, 0.1])
        
        for s in stst:
            thisstate = thisdesign[thisdesign['StableStateAll']== s]
            ax.scatter(thisstate['kappa'], thisstate.iloc[:,5+i]/np.pi, c = colors[thisstate['StableStateAll']-1], label = s)
            ax.scatter(thisstate['kappa'], -thisstate.iloc[:,5+i]/np.pi, c = colors[thisstate['StableStateAll']-1], marker ='^')
    
    leg = axes[-1].legend(loc = 0, fontsize = 7, framealpha = 0.8, edgecolor = 'inherit', fancybox = False) 
    #           borderpad = 0.3, labelspacing = 0.1, handlelength = 0.4, handletextpad = 0.4)
    plt.setp(leg.get_texts(), color='0.2')
    leg.get_frame().set_linewidth(0.4)
    
    fig1.show()
    fig1.savefig(Folder_name + '/Images/SingleDes/Angles_' + str(d[0].astype(int))+ '_' + str(d[3].astype(int)) + '.pdf', transparent = True)
    fig1.savefig(Folder_name + '/Images/SingleDes/Angles_' + str(d[0].astype(int))+ '_' + str(d[3].astype(int)) + '.png', transparent = True)



#%%
#fig3 = plt.figure(figsize=(cm2inch(8), cm2inch(6)))
#ax3 = plt.subplot(111)
#NiceGraph2D(ax3, 'Desang1', 'Desang2', mincoord=[0,0], maxcoord=[180,180], divisions=[6,6],buffer=[5,5])
#
#onekappa = allKappasAnalysis[allKappasAnalysis['kappa'] == kappas[0]]
#onekappa['desang0'] = np.ones(np.shape(onekappa['desang1']))*180
#
#
#ax3.scatter(onekappa['desang1'],onekappa['desang3'], 
#            c = onekappa[['desang1','desang3','desang0','desang0']].values/180)
#ax3.scatter(onekappa['desang3'],onekappa['desang1'], 
#            c = onekappa[['desang3','desang1','desang0','desang0']].values/180)
#onekappa['desang0'] = onekappa['desang0']*0
#ax3.scatter(180-onekappa['desang1'],180-onekappa['desang3'], 
#            c = 1-onekappa[['desang1','desang3','desang0','desang0']].values/180)
#ax3.scatter(180-onekappa['desang3'],180-onekappa['desang1'], 
#            c = 1-onekappa[['desang3','desang1','desang0','desang0']].values/180)



fig3 = plt.figure(figsize=(cm2inch(17.8), cm2inch(7)))
ax3 = plt.subplot(111,projection='3d')
ax3.set_xlim([-np.pi,np.pi])
ax3.set_ylim([-np.pi,np.pi])
ax3.set_zlim([-np.pi,np.pi])
    
#cmap2 = matl.cm.get_cmap('Set2',np.size(kappas))
#colors = cmap2(np.linspace(0,1,np.size(kappas)))    
cmap2 = matl.cm.get_cmap('Set3',np.size(stst))
colors = cmap2(np.linspace(0,1,np.size(stst)))

order = [5,6,7,8]

for kappa in kappas:
    
#    fig3 = plt.figure(figsize=(cm2inch(8), cm2inch(6)))
#    ax3 = plt.subplot(111,projection='3d')
#    
#    ax3.set_xlim([-np.pi,np.pi])
#    ax3.set_ylim([-np.pi,np.pi])
#    ax3.set_zlim([-np.pi,np.pi])
    
    thisstate = allDesigns[allDesigns['kappa'] == kappa]
#    thisstate['desang0'] = np.zeros(np.shape(thisstate['desang1']))

    if not thisstate.empty:
        
#        for order in np.array(list(itertools.permutations([5,6,7,8],3)))[[0,9,16,18,5,7,14,23]]:
#            thisstate['desang0'] = thisstate['desang0']+180
            ax3.scatter(thisstate.iloc[:,order[0]].values,thisstate.iloc[:,order[1]].values,thisstate.iloc[:,order[2]].values, 
#                        c = (thisstate.iloc[:,[order[0]+7,order[2]+7,17,17]].values)/180)
#                        c = colors[thisstate['LogKappas'].astype(int).values+3])
                        c = colors[thisstate['StableStateAll']-1])
#            ax3.scatter(thisstate.iloc[:,order[2]].values,thisstate.iloc[:,order[3]].values,thisstate.iloc[:,order[0]].values, 
#                        c = (thisstate.iloc[:,[order[2]+7,order[0]+7,17,17]].values)/180)
##                        c = colors[thisstate['LogKappas'].astype(int).values+3])
##                        c = colors[thisstate['StableStateAll']-1])
#            thisstate['desang0'] = thisstate['desang0']*0
#            ax3.scatter(thisstate.iloc[:,order[0]].values,thisstate.iloc[:,order[3]].values,thisstate.iloc[:,order[2]].values, 
#                        c = 1-(thisstate.iloc[:,[order[0]+7,order[2]+7,17,17]].values)/180)
##                        c = colors[thisstate['LogKappas'].astype(int).values+3])
##                        c = colors[thisstate['StableStateAll']-1])
#            ax3.scatter(thisstate.iloc[:,order[2]].values,thisstate.iloc[:,order[1]].values,thisstate.iloc[:,order[0]].values, 
#                        c = 1-(thisstate.iloc[:,[order[2]+7,order[0]+7,17,17]].values)/180)
##                        c = colors[thisstate['LogKappas'].astype(int).values+3])
##                        c = colors[thisstate['StableStateAll']-1])

#%%
allDesigns[['kappa','Hinge Number','StableStateAll','restang']].to_csv(Folder_name + '/Images/InfoforImages.csv', index = False)

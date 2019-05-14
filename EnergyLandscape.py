import matplotlib as matl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import from_levels_and_colors
import numpy as np
import configparser
import os.path
import scipy.cluster.hierarchy as hierarch
from scipy.spatial.distance import pdist
import copy
import matplotlib.animation as animation
matl.rcParams['pdf.fonttype'] = 42
matl.rcParams['ps.fonttype'] = 42
matl.rcParams['font.family'] = 'sans-serif'
matl.rcParams['font.sans-serif'] = 'Arial'
matl.rcParams['mathtext.fontset'] = 'cm'


def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False
    
def cm2inch(value):
    return value/2.54

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
    axes.xaxis.labelpad = -5
    if ~np.isnan(mincoord[1]) and ~np.isnan(maxcoord[1]):
        axes.set_ylim([mincoord[1]-buffer[1], maxcoord[1]+buffer[1]])
        if isinstance(divisions[1], (list, tuple, np.ndarray)):
            if ~np.isnan(divisions[1]).any():
                axes.set_yticks(divisions[1])
        else:
            if ~np.isnan(divisions[1]):
                axes.set_yticks(np.linspace(mincoord[1],maxcoord[1],divisions[1]))
    axes.set_ylabel(nameY)
    axes.yaxis.labelpad = -5
   
    axes.xaxis.label.set_color(gray)
    axes.tick_params(axis='x',colors=gray, width=0.4)
    axes.spines['bottom'].set_color(gray)
    axes.spines['top'].set_color(gray)
    axes.yaxis.label.set_color(gray)
    axes.tick_params(axis='y', colors=gray, width=0.4)
    axes.spines['left'].set_color(gray)
    axes.spines['right'].set_color(gray)
    axes.tick_params(pad = 2)
    
    for axis in ['top','bottom','left','right']:
        axes.spines[axis].set_linewidth(0.4)
    
    return

def NiceGraph3D(axes, nameX, nameY, nameZ, mincoord = [np.NaN, np.NaN, np.NaN], maxcoord = [np.NaN, np.NaN, np.NaN],
                divisions = [np.NaN, np.NaN, np.NaN], buffer = [0.0, 0.0, 0.0]):
    gray = '0.2'
    matl.rcParams.update({'font.size': 9})

    if ~np.isnan(mincoord[0]) and ~np.isnan(maxcoord[0]):
        axes.set_xlim3d([mincoord[0]-buffer[0], maxcoord[0]+buffer[0]])
        if ~np.isnan(divisions[0]):
            axes.set_xticks(np.linspace(mincoord[0],maxcoord[0],divisions[0]))
    axes.set_xlabel(nameX)
#    axes.xaxis.labelpad = 20
    if ~np.isnan(mincoord[1]) and ~np.isnan(maxcoord[1]):
        axes.set_ylim3d([mincoord[1]-buffer[1], maxcoord[1]+buffer[1]])
        if ~np.isnan(divisions[1]):
            axes.set_yticks(np.linspace(mincoord[1],maxcoord[1],divisions[1]))
    axes.set_ylabel(nameY)
#    axes.yaxis.labelpad = 20
    if ~np.isnan(mincoord[2]) and ~np.isnan(maxcoord[2]):
        axes.set_zlim3d([mincoord[2]-buffer[2], maxcoord[2]+buffer[2]])
        if ~np.isnan(divisions[2]):
            axes.set_zticks(np.linspace(mincoord[2],maxcoord[2],divisions[2]))
    axes.set_zlabel(nameZ)
#    axes.zaxis.labelpad = 15
    
    axes.w_xaxis.set_pane_color((0,0,0,0))
    axes.w_yaxis.set_pane_color((0,0,0,0))
    axes.w_zaxis.set_pane_color((0,0,0,0))
    
    axes.w_xaxis.line.set_color(gray)
    axes.w_yaxis.line.set_color(gray)
    axes.w_zaxis.line.set_color(gray)
    
    axes.w_xaxis.label.set_color(gray)
    axes.tick_params(axis='x',colors=gray, width=0.4)
    axes.w_yaxis.label.set_color(gray)
    axes.tick_params(axis='y', colors=gray, width=0.4)
    axes.w_zaxis.label.set_color(gray)
    axes.tick_params(axis='z', colors=gray, width=0.4)    
    return

folder_name = "D:/Documents/Git Programs/nonlinear-bas/Results/SquareTiling/sqp/energy/09-May-2019_LayersUC_Analysis/kh0.000_kta100.000_ke10.000_kf100.000"
#inverted = False
release = False
plt.close('all')
#%%
#######################################################################################################################
##################### Reading Files
#######################################################################################################################

file_name1 = "/EnergyData.csv" 
file_name2 = "/Hinges.csv"
file_name3 = "/PosStad.csv"
file_name4 = "/Angles.csv"

dataEnergy = np.loadtxt(folder_name+file_name1,skiprows=1, delimiter = ',', unpack = True)
hingeNum = dataEnergy[0,:].astype(int)
eEdge = dataEnergy[1,:]
eDiag = dataEnergy[2,:]
eFace = dataEnergy[3,:]
eHinge = dataEnergy[4,:]
eTAngle = dataEnergy[5,:]
exfl = dataEnergy[6,:].astype(int)

closingAngl1 = np.loadtxt(folder_name+file_name2,skiprows=1, delimiter = ',', unpack = True, usecols = [3])
closingAngl2 = np.loadtxt(folder_name+file_name2,skiprows=1, delimiter = ',', unpack = True, usecols = [4])
#angleNum = np.loadtxt(folder_name+file_name2,skiprows=1, delimiter = ',', usecols = [3,4])

dataPos = np.loadtxt(folder_name+file_name3,skiprows=1, delimiter = ',', unpack = True)
lowerRadius = dataPos[1,:]
upperRadius = dataPos[2,:]

#dataAngles = np.loadtxt(folder_name+file_name4,skiprows=1, delimiter = ',', dtype = np.float64)
#dataAngles = np.delete(dataAngles, 0, 1)

TotSimul = np.size(closingAngl1)
#IterPerSimul = np.int(np.size(dataAngles,0)/TotSimul)
IterPerSimulEnergy = np.int(np.size(dataEnergy,1)/TotSimul)
#NumberAngles = np.size(dataAngles,1)

metadataFile = '/metadata.txt'

if os.path.isfile(folder_name+metadataFile):
    metadata = configparser.RawConfigParser(allow_no_value=True)
    metadata.read(folder_name+metadataFile)

    totalnumberHinges = int(metadata.get('extUnitCell', 'Hinges'))
    totalnumberEdges = int(metadata.get('extUnitCell', 'Edges'))
    totalnumberDiag = int(metadata.get('extUnitCell', 'Diag'))
    khinge = float(metadata.get('options','kHinge'))
    kedge = float(metadata.get('options','kEdge'))
    kdiag = float(metadata.get('options','kDiag'))
    kface = float(metadata.get('options','kFace'))
    alpha = metadata.get('extUnitCell', 'alpha')
    alpha = np.array([int(s) for s in alpha.split('  ')])
    beta = metadata.get('extUnitCell', 'beta')
    beta = np.array([int(s) for s in beta.split('  ')])
else:
    raise FileNotFoundError('No metafile found at the given directory. Changes to the script to put manually the variables are needed\n') 

#eEdge = eEdge/kedge/totalnumberEdges
#eDiag = eDiag/kdiag/totalnumberDiag
#eHinge = eHinge/khinge/totalnumberHinges

eTotal= eEdge + eDiag + eHinge
#eTotal = eTotal/totMaxEnergy*100

if release:
    shift = 2
else:
    shift = 1

#%%
###Check if a folding pattern didnt converge
exfl = exfl.reshape(TotSimul,IterPerSimulEnergy)
flagmask = np.logical_and(np.logical_and(exfl !=2, exfl !=0),np.logical_and(exfl !=2, exfl !=1))
flagmask = flagmask.any(axis= 1)
if flagmask.any():
    print('Error: There was at least one non convergent fold pattern.\n')

##Do the energy landscape plot
divitheta1 = len(np.unique(closingAngl1))
divitheta2 = len(np.unique(closingAngl2))

tickstheta1 = 4
tickstheta2 = 5

sortAngl = np.lexsort((closingAngl2, closingAngl1))
flagmask_n = np.logical_not(flagmask[sortAngl])
closingAngl1 = closingAngl1[sortAngl]
closingAngl2 = closingAngl2[sortAngl]

#angleNum = angleNum[sortAngl[::-1],:]
#initial_pos = np.logical_and(angleNum[:,0] == 1,angleNum[:,1] == 1)

theta1 = closingAngl1.reshape((divitheta1,divitheta2))
theta2 = closingAngl2.reshape((divitheta1,divitheta2))

totEnergysort = eTotal[IterPerSimulEnergy-shift::IterPerSimulEnergy]
#totEnergysort = eDiag[IterPerSimulEnergy-shift::IterPerSimulEnergy]+eEdge[IterPerSimulEnergy-shift::IterPerSimulEnergy]
#totEnergysort = eEdge[IterPerSimulEnergy-shift::IterPerSimulEnergy]
#totEnergysort = eDiag[IterPerSimulEnergy-shift::IterPerSimulEnergy]
#totEnergysort = eHinge[IterPerSimulEnergy-shift::IterPerSimulEnergy]
#totEnergysort = eTAngle[IterPerSimulEnergy-shift::IterPerSimulEnergy]

#totEnergysort_n = totEnergysort
totEnergysort = np.ma.masked_array(totEnergysort, mask=flagmask)
totEnergysort = totEnergysort[sortAngl]
totEnergysort[totEnergysort == 0] = np.nan

totEnergyMat = totEnergysort.reshape((divitheta1,divitheta2))


sep1 = (closingAngl1[-1]-closingAngl1[0])/(divitheta1-1)/2
sep2 = (closingAngl2[-1]-closingAngl2[0])/(divitheta2-1)/2

fig1 = plt.figure(0,figsize=(cm2inch(8.7), cm2inch(7)))
fig1.subplots_adjust(top=0.99,
bottom=0.115,
left=0.130,
right=0.89)
ax1 = plt.subplot(111)

NiceGraph2D(ax1, r'$Unit Cells$', r'$Layers$', 
            mincoord = [theta2[0,0],theta1[0,0]], maxcoord = [theta2[0,-1],theta1[-1,0]],
            divisions = [tickstheta1, tickstheta2], buffer = [sep1, sep2])

maxEnergy = 116#2.4  #np.max(totEnergyMat) #0.05#16
cs1 = ax1.imshow(totEnergyMat, extent=[theta2[0,0]-sep2,theta2[0,-1]+sep2,theta1[0,0]-sep1,theta1[-1,0]+sep1], 
                     cmap = cm.nipy_spectral, aspect = 'auto',vmax = maxEnergy, origin = 'lower') #nipy_spectral

#SStheta1 = -dataAngles[1+shift::IterPerSimul,beta-1]/np.pi
#SStheta1 = np.mean(SStheta1,axis = 1)
#SStheta1 = SStheta1[sortAngl]
#SStheta2 = dataAngles[1+shift::IterPerSimul,alpha-1]/np.pi
#SStheta2 = np.mean(SStheta2,axis = 1)
#SStheta2 = SStheta2[sortAngl]

ax1.set_xticklabels(["$%d$" %theta2[0,0], "", "", r"$%d$" %theta2[0,-1]])
ax1.set_yticklabels(["$%d$" %theta1[0,0], "", "", "", r"$%d$" %theta1[-1,0]])

cbar = plt.colorbar(cs1,  ax = ax1, fraction=0.05, pad=0.01, extend = 'max', format="%.1f") 
cbar.set_ticks(np.linspace(0, maxEnergy, 4))
cbar.set_label('Energy', fontsize = 9, color = '0.2',labelpad = 0)
cbar.ax.tick_params(colors='0.2', pad=2, width=0.4)
cbar.outline.set_edgecolor('0.2')
cbar.outline.set_linewidth(0.4)

fig1b = plt.figure(1,figsize=(cm2inch(8.7), cm2inch(7)))
ax1b = fig1b.add_subplot(111, projection='3d')
NiceGraph3D(ax1b, r'$Unit Cells$', r'$Layers$', 'log(Energy)',
            mincoord = [theta2[0,0],theta1[0,0],np.nan], maxcoord = [theta2[0,-1],theta1[-1,0],np.nan],
            divisions = [tickstheta1, tickstheta2, np.nan], buffer = [sep1, sep2, np.nan])

#ax1b.set_zlim(0, 0.000005)
ax1b.plot_wireframe(theta2, theta1, np.log10(totEnergyMat), rstride=1, cstride=1)
ax1b.set_xticklabels(["$%d$" %theta2[0,0], "", "", r"$%d$" %theta2[0,-1]])
ax1b.set_yticklabels(["$%d$" %theta1[0,0], "", "", "", r"$%d$" %theta1[-1,0]])
#    
#def update_lines(num):
#    ax1b.view_init(30, num)
#    return fig1b,
#
#line_ani = animation.FuncAnimation(fig1b, update_lines, 360, interval=10, blit=True)
#
#line_ani.save(folder_name + '/EnergyLand.mp4', fps=20, extra_args=['-vcodec', 'libx264'], bitrate = -1, dpi = 300)

######################################################################
#Analysis for stable states

#finalAngles = np.empty((0,np.size(dataAngles,1)))
#for hinge in sortAngl:
#    sortAllAngIndex = np.lexsort((dataAngles[IterPerSimul*(hinge+1)-1,:],dataAngles[IterPerSimul*hinge,:]))
#    finalAngles = np.append(finalAngles, [dataAngles[IterPerSimul*(hinge+1)-1,sortAllAngIndex]], axis = 0)
#
#eStretch = eEdge + eDiag
#totEStretchsort = eStretch[IterPerSimulEnergy-shift::IterPerSimulEnergy]
#totEStretchsort = totEStretchsort[sortAngl]
#
#Z = hierarch.linkage(finalAngles[flagmask_n,:], 'centroid')
#inverse_masked = hierarch.fcluster(Z, 1, criterion='distance')
#c = hierarch.cophenet(Z, pdist(finalAngles[flagmask_n,:]))
#print('this is the cophenet of the hierarchical linkage', c[0])
#
##    plt.figure(figsize=(25, 10))
##    plt.title('Hierarchical Clustering Dendrogram')
##    plt.xlabel('sample index')
##    plt.ylabel('distance')
##    hierarch.dendrogram(
##        Z,
##        truncate_mode='lastp',  # show only the last p merged clusters
##        p=50,  # show only the last p merged clusters
##        leaf_rotation=90.,  # rotates the x axis labels
##        leaf_font_size=8.,  # font size for the x axis labels
##        show_contracted=True,
##    )
##    plt.show()
#
#AllSimul = np.arange(TotSimul)
#AllSimul = AllSimul[flagmask_n]
#SS, SSpos = np.unique(inverse_masked, return_index = True)
#SSpos = AllSimul[SSpos]

######################################################################################################



#for i in SSpos:
#    print('Angles: ', angleNum[i,:], '\tStretching Energy: ', totEStretchsort[i], '\tTotal Energy: ', totEnergysort[i])
#
#Simul = sortAngl
#Simul = Simul[flagmask_n]
#inverse = np.zeros(TotSimul)
#inverse[Simul] = inverse_masked
#inverse = inverse[sortAngl]
#inverse = np.ma.masked_array(inverse, mask=flagmask[sortAngl])
#
#if inverted:
#    stableStateMat = inverse.reshape((divitheta1,divitheta2))
#else:
#    stableStateMat = np.rot90(np.rot90(inverse.reshape((divitheta1,divitheta2)).T))
#
#fig2 = plt.figure(2,figsize=(cm2inch(8.7), cm2inch(7)))
#fig2.subplots_adjust(top=0.99,
#bottom=0.115,
#left=0.130,
#right=0.89)
#ax2 = plt.subplot(111)
#
#NiceGraph2D(ax2, r'$\alpha$', r'$\beta$', 
#            mincoord = [0,0], maxcoord = [0.5,0.5],divisions = [tickstheta1, tickstheta2], buffer = [sep1, sep2])
#
#cmap2, norm2 = from_levels_and_colors(np.linspace(0,12,13),#np.int(np.max(inverse)),np.int(np.max(inverse)+1)),
#                                      cm.Set3(np.linspace(0, 1,12)))# np.int(np.max(inverse))))) #gist_rainbow #Set3
#
#cs2 = ax2.imshow(stableStateMat, extent=[theta2[0,0]-sep2,theta2[0,-1]+sep2,theta1[0,0]-sep1,theta1[-1,0]+sep1], 
#                     cmap = cmap2, aspect = 'auto', origin = 'lower', vmax=12)
#
#ax2.scatter(SStheta1[SSpos], SStheta2[SSpos], c = inverse[SSpos], cmap = cmap2, s = 100, marker = '*', edgecolor = 'k', lw = 0.2, vmax = 12)              
#
#
#ax2.set_xticklabels(["$0$", "", "", "", r"$\pi/2$"])
#ax2.set_yticklabels(["$0$", "", "", "", r"$\pi/2$"])
#
#
#cbar2 = plt.colorbar(cs2,  ax = ax2, fraction=0.05)#, extend = 'max'
#cbar2.set_ticks([])
#cbar2.ax.tick_params(colors='0.2', pad=0)
#cbar2.outline.set_edgecolor('0.2')
#cbar2.outline.set_linewidth(0.4)

#############################################################################
#Plot of real angles

#realtheta1 = dataAngles[2::IterPerSimul,alpha-1]/np.pi
#realtheta1 = np.mean(realtheta1,axis = 1)
#realtheta1 = realtheta1[sortAngl]
#realtheta2 = -dataAngles[2::IterPerSimul,beta-1]/np.pi
#realtheta2 = np.mean(realtheta2,axis = 1)
#realtheta2 = realtheta2[sortAngl]
#
#totEnergysort = eTAngle[IterPerSimulEnergy-shift::IterPerSimulEnergy]
#totEnergysort = np.ma.masked_array(totEnergysort, mask=flagmask)
#totEnergysort = totEnergysort[sortAngl]
#
#maxTA = 0.12 #np.max(totEnergysort)
#
#fig3 = plt.figure(4,figsize=(cm2inch(8.7), cm2inch(7)))
#fig3.subplots_adjust(top=0.99,
#bottom=0.115,
#left=0.130,
#right=0.89)
#ax3 = plt.subplot(111)
#NiceGraph2D(ax3, r'$\alpha$', r'$\beta$',
#            mincoord = [theta2[0,0],theta1[0,0]], maxcoord = [theta2[0,-1],theta1[-1,0]],
#            divisions = [tickstheta1, tickstheta2], buffer = [sep1, sep2])
#
#cmap3, norm3 = from_levels_and_colors(np.linspace(0,maxTA,1000), cm.rainbow(np.linspace(0, 1, 1000)), extend = 'max')
#
#cs3 = ax3.scatter(realtheta1, realtheta2, c = totEnergysort, cmap = cmap3, vmax = maxTA, s = 20, marker = 's') #150 #360
#
#ax3.set_xticklabels(["$%.2f\pi$" %theta2[0,0], "", "", "", r"$%.2f\pi$" %theta2[0,-1]])
#ax3.set_yticklabels(["$%.2f\pi$" %theta1[0,0], "", "", "", r"$%.2f\pi$" %theta1[-1,0]])
#
#cbar3 = plt.colorbar(cs3, ax = ax3, fraction=0.05, pad=0.01, extend = 'max', format="%.2f")
#cbar3.set_ticks(np.linspace(0, maxTA, 5))
#cbar3.set_label('Energy', fontsize = 11, color = '0.2')
#cbar3.ax.tick_params(colors='0.2', width=0.4)
#cbar3.outline.set_edgecolor('0.2')
#cbar3.outline.set_linewidth(0.4)

################################################################################

#lowerRadius[np.isnan(lowerRadius)]=10**10
lowerRadius = np.array(lowerRadius[sortAngl])
#upperRadius[np.isnan(upperRadius)]=10**10
upperRadius = np.array(upperRadius[sortAngl])
maxRadius = 0
minRadius = -2

#if inverted:
RadiusMat = upperRadius.reshape((divitheta1,divitheta2))
#else:
#    RadiusMat = np.rot90(np.rot90(lowerRadius.reshape((divitheta1,divitheta2)).T))

CurvatureMat = 1/RadiusMat

fig4 = plt.figure(5,figsize=(cm2inch(8.7), cm2inch(7)))
fig4.subplots_adjust(top=0.970,
bottom=0.115,
left=0.130,
right=0.875)
ax4 = plt.subplot(111)
#ax5 = plt.subplot(122)
NiceGraph2D(ax4, r'$Unit Cells$', r'$Layers$', 
            mincoord = [theta2[0,0],theta1[0,0]], maxcoord = [theta2[0,-1],theta1[-1,0]],
            divisions = [tickstheta1, tickstheta2], buffer = [sep1, sep2])

#cmap4, norm4 = from_levels_and_colors(np.linspace(minRadius,maxRadius,1000), cm.rainbow_r(np.linspace(0, 1, 1000)), extend = 'min')

#cs4 = ax4.scatter(realtheta1, realtheta2, c = np.log10(upperRadius), cmap = cmap4, vmax = maxRadius, vmin = minRadius, s = 50, marker = 's') #150 #360
cs4 = ax4.imshow(np.log10(CurvatureMat), extent=[theta2[0,0]-sep2,theta2[0,-1]+sep2,theta1[0,0]-sep1,theta1[-1,0]+sep1], 
                     cmap = cm.jet, aspect = 'auto',vmax = maxRadius, vmin = minRadius, origin = 'lower') #nipy_spectral



ax4.set_xticklabels(["$%d$" %theta2[0,0], "", "", r"$%d$" %theta2[0,-1]])
ax4.set_yticklabels(["$%d$" %theta1[0,0], "", "", "", r"$%d$" %theta1[-1,0]])

cbar4 = plt.colorbar(cs4, ax = ax4, fraction=0.05, pad=0.01, extend = 'min', format=r"$10^{%.1f}$")
cbar4.set_ticks(np.linspace(minRadius,maxRadius,5))
cbar4.set_label('Energy', fontsize = 11, color = '0.2')
cbar4.ax.tick_params(colors='0.2', width=0.4)
cbar4.outline.set_edgecolor('0.2')
cbar4.outline.set_linewidth(0.4)

#angleMat = np.zeros((divitheta1,divitheta2,2))
#angleMat[:,:,0] = angleNum[:,0].reshape((divitheta1,divitheta2))
#angleMat[:,:,1] = angleNum[:,1].reshape((divitheta1,divitheta2))
#angleMat[:,:,0] = np.rot90(np.rot90(angleNum[:,0].reshape((divitheta1,divitheta2)).T))
#angleMat[:,:,1] = np.rot90(np.rot90(angleNum[:,1].reshape((divitheta1,divitheta2)).T))









#############################################################################
fig1.show()
fig1.savefig(folder_name + '/EnergyLand.pdf', transparent = True, pad_inches=0, dpi=300)
fig1.savefig(folder_name + '/EnergyLand.png', transparent = True, pad_inches=0, dpi=300)
#fig1.savefig(folder_name + '/EnergyAllEdges.png', transparent = True)
#fig1.savefig(folder_name + '/EnergyEdge.png', transparent = True)
#fig1.savefig(folder_name + '/EnergyDiag.png', transparent = True)
#fig1.savefig(folder_name + '/EnergyHinge.png', transparent = True)
#fig1.savefig(folder_name + '/EnergyTA.png', transparent = True)

#fig2.show()
#fig2.savefig(folder_name + '/StableStates.pdf', transparent = True, pad_inches=0, dpi=300)
#fig2.savefig(folder_name + '/StableStates.png', transparent = True, pad_inches=0, dpi=300)

#fig3.show()
#fig3.savefig(folder_name + '/RealAngles-EnergyTA.pdf', transparent = True, pad_inches=0, dpi=300)
#fig3.savefig(folder_name + '/RealAngles-EnergyTA.png', transparent = True, pad_inches=0, dpi=300)

fig4.show()
fig4.savefig(folder_name + '/Curvature.pdf', transparent = True, pad_inches=0, dpi=300)
fig4.savefig(folder_name + '/Curvature.png', transparent = True, pad_inches=0, dpi=300)
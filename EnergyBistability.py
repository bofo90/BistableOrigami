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

def NiceGraph2DCenter(axes, nameX, nameY, mincoord = [np.NaN, np.NaN], maxcoord = [np.NaN, np.NaN], divisions = [np.NaN, np.NaN],
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
    axes.set_xlabel(nameX,labelpad=-3, color = gray)
    
    if ~np.isnan(mincoord[1]) and ~np.isnan(maxcoord[1]):
        axes.set_ylim([mincoord[1]-buffer[1], maxcoord[1]+buffer[1]])
        if isinstance(divisions[1], (list, tuple, np.ndarray)):
            if ~np.isnan(divisions[1]).any():
                axes.set_yticks(divisions[1])
        else:
            if ~np.isnan(divisions[1]):
                axes.set_yticks(np.linspace(mincoord[1],maxcoord[1],divisions[1]))
    axes.set_ylabel(nameY,labelpad=-3, color = gray,rotation=0)
   
    axes.xaxis.label.set_color(gray)
    axes.tick_params(axis='x', colors=gray, direction = 'out', width = 0.4)
    axes.yaxis.label.set_color(gray)
    axes.tick_params(axis='y', colors=gray, direction = 'inout', width = 0.4, length = 15)
    axes.tick_params(pad = 2)
    
    axes.tick_params(axis='y', which='minor', colors=gray, direction = 'in', width = 0.4)
    axes.tick_params(axis='x', which='minor', colors=gray, direction = 'in', width = 0.4)
    
    for axis in ['top','bottom','left','right']:
        axes.spines[axis].set_linewidth(0.4)
        axes.spines[axis].set_color(gray)
    
    axes.spines['left'].set_position(('data', 0.0))
    axes.spines['bottom'].set_position(('data', 0.0))
    axes.spines['right'].set_color('none')
    axes.spines['top'].set_color('none')    
    return

def force(x, kappa, restang):
    x = x*np.pi
    restang = restang*np.pi
    Theta0 = np.ones(np.shape(x))*restang
    ThetaVar = x-Theta0
    
    F = 4*kappa/(restang**4)*(0.25*ThetaVar**4+ThetaVar**3.*Theta0+ThetaVar**2.*Theta0**2)
    
    return F

def force2(x, kappa, restang):
    x = x*np.pi
    restang = restang*np.pi
    Theta0 = np.ones(np.shape(x))*restang
    
    F = kappa/Theta0**4*(x**2-Theta0**2)**2
    return F

plt.close('all')

angle = np.linspace(-1.5,1.5,100000)
restang = 0.5
restangles = np.array([0.25,0.5,0.75])
kappa = 1

fig1 = plt.figure(figsize=(cm2inch(3.8), cm2inch(3.3)))
ax1 = plt.subplot(111)
fig1.subplots_adjust(top=0.995,
bottom=0.22,
left=0.020,
right=0.985)

NiceGraph2DCenter(ax1, r'$\theta$', r'$E_\mathregular{h}$', mincoord = [-1.6*np.max(restangles),0], maxcoord = [1.6*np.max(restangles),1.5*kappa], 
            divisions = [[-1,-0.5,0,0.5,1],[kappa]], buffer = [0, 0])

ax1.set_xticklabels([r"-180°",r"-90°",r"0°",r"90°",r"180°"])
ax1.set_yticklabels([r"$\kappa$"])
ax1.xaxis.set_label_coords(0.97, 0.13)
ax1.yaxis.set_label_coords(0.57, 0.85)

cmap = matl.cm.get_cmap('summer_r',np.size(restangles)+1)
colors = cmap(np.linspace(0,1,np.size(restangles)+1))

labels = [r"45°", r"90°",r"135°"]

for restang, i in zip(restangles, np.arange(np.size(restangles))):
# for kappa in [0.01,0.1, 1]:

    ax1.plot(angle, force(angle, kappa, restang), color = colors[i+1], label = labels[i], linewidth = 0.7)

leg = ax1.legend(loc = 1, fontsize = 7, framealpha = 1, edgecolor = 'inherit', fancybox = False, 
                 labelspacing = 0.1, handlelength = 0.5, handletextpad = 0.5, borderaxespad = 0, borderpad = 0.3)
plt.setp(leg.get_texts(), color='0.2')
leg.get_frame().set_linewidth(0)

fig1.show()
fig1.savefig('D:/Documents/Git Programs/nonlinear-bas_Origami/Results/Bistability.pdf', transparent = True)
fig1.savefig('D:/Documents/Git Programs/nonlinear-bas_Origami/Results/Bistability.png', transparent = True)


#%%

fig2 = plt.figure(figsize=(cm2inch(7), cm2inch(5)))
ax2 = plt.subplot(111)
fig2.subplots_adjust(top=0.98,
bottom=0.145,
left=0.015,
right=0.98)

NiceGraph2DCenter(ax2, r'$\theta$', r'$E_\mathregular{H}$', mincoord = [-1.6*np.max(restangles),0], maxcoord = [1.6*np.max(restangles),1.5*kappa], 
            divisions = [[-0.5,0,0.5],[kappa]], buffer = [0, 0])

ax2.set_xticklabels([r"$-\theta_0$",r"$0$",r"$\theta_0$"])
ax2.set_yticklabels([r"$k_\mathregular{H}$"])
ax2.xaxis.set_label_coords(0.97, 0.13)
ax2.yaxis.set_label_coords(0.57, 0.85)

i = 1
restang = restangles[i]
    
ax2.plot(angle, force(angle, kappa, restang), color = '#36648B', linewidth = 1)

fig2.savefig('D:/Documents/Git Programs/nonlinear-bas_Origami/Results/BistabilitySingle.pdf', transparent = True)
fig2.savefig('D:/Documents/Git Programs/nonlinear-bas_Origami/Results/BistabilitySingle.png', transparent = True)

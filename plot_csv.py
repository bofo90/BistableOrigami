import matplotlib as matl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pylab as P
import csv
import glob
import ntpath
import sys

plt.close('all')

def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False
    
def cm2inch(value):
    return value/2.54

def NiceGraph2Dwticks(axes, X, Y, mincoord, maxcoord, xdivisions, ydivisions, buffery):
    gray = '0.2'
    matl.rcParams.update({'font.size': 15})
    bufferx = 0.1

    axes.set_xlim([mincoord[0]-bufferx, maxcoord[0]+bufferx])
    axes.set_xticks(np.linspace(mincoord[0],maxcoord[0],xdivisions))
    axes.set_xlabel(X)
    axes.set_ylim([mincoord[1]-buffery, maxcoord[1]+buffery])
    axes.set_yticks(np.linspace(mincoord[1],maxcoord[1],ydivisions))
    axes.set_ylabel(Y)
   
    axes.xaxis.label.set_color(gray)
    axes.tick_params(axis='x',colors=gray)
    axes.spines['bottom'].set_color(gray)
    axes.spines['top'].set_color(gray)
    axes.yaxis.label.set_color(gray)
    axes.tick_params(axis='y', colors=gray)
    axes.spines['left'].set_color(gray)
    axes.spines['right'].set_color(gray)
    
    axes.grid( 'on', axis='y', color ='0.2', alpha = 0.3)

    return

def NiceGraph2D(axes, X, Y):
    gray = '0.2'
    matl.rcParams.update({'font.size': 15})
    
    axes.set_xlabel(X)
    axes.set_ylabel(Y)
   
    axes.xaxis.label.set_color(gray)
    axes.tick_params(axis='x',colors=gray)
    axes.spines['bottom'].set_color(gray)
    axes.spines['top'].set_color(gray)
    axes.yaxis.label.set_color(gray)
    axes.tick_params(axis='y', colors=gray)
    axes.spines['left'].set_color(gray)
    axes.spines['right'].set_color(gray)
    return

idealAngle = [];
realAngle = [];

enEdge = [];
enFace = [];
enHinge = [];
enTotal = [];

kAngleConst = [];
errorAngle = []; 

datapoints = 80
nruns = 8

dirname = "DataEnergy\\"
file_name = "khinge.csv"

with open(dirname+file_name,'r',newline='') as ffit: 
    readers = csv.reader(ffit,delimiter=',');
    for row in readers:
        if row:
            if isfloat(row[0].encode('ascii','ignore')):
                idealAngle.append(np.float(row[0]));
                realAngle.append(np.float(row[1]));
                enEdge.append(np.float(row[2]));
                enFace.append(np.float(row[3]));
                enHinge.append(np.float(row[4]));
                enTotal.append(np.float(row[5]));
                kAngleConst.append(np.float(row[6]))
                errorAngle.append(np.float(row[7]))


#Make all data as float
idealAngle = np.array(idealAngle);
realAngle = np.array(realAngle);

enEdge = np.array(enEdge);
enFace = np.array(enFace);
enHinge = np.array(enHinge);
enTotal = np.array(enTotal);

kAngleConst = np.array(kAngleConst);
errorAngle = np.array(errorAngle);

enAtAngle = []
normEnHinge = []
enHingekHinge = []
kvector = []
selection = 59

for angle in np.arange(nruns):
    enAtAngle.append(np.float(enHinge[angle*datapoints+selection]))
    kvector.append(np.float(kAngleConst[angle*datapoints+selection]))
    if np.float(kAngleConst[angle*datapoints+selection]) != 0:
        enHingekHinge.append(np.float(enHinge[angle*datapoints+selection])/np.float(kAngleConst[angle*datapoints+selection]))
    else:
        enHingekHinge.append(np.float(0.0))

for run in np.arange(nruns):
    for angle in np.arange(datapoints):
        if np.float(enAtAngle[run]) == 0:
            normEnHinge.append(np.float(0.0))
        else:
            normEnHinge.append(np.float(enHinge[run*datapoints+angle])/np.float(enAtAngle[run]))


fig = plt.figure(0,figsize=(cm2inch(35), cm2inch(20)))
ax1 = plt.subplot(121)
ax2 = plt.subplot(122)

offset1 = 0.1
offset2 = 0.003
offset3 = 0.15

NiceGraph2Dwticks(ax1, 'Real Angle [rad]', 'Energy [a.u.] (with offset)', [-3.1416, 0.],[0,1.2], 5, 13, 0.01)
NiceGraph2Dwticks(ax2, 'Real Angle [rad]', 'ErrorAngle [rad] (with offset of %.3f rad)' %offset2, [-3.1416, 0.0],[0,0.024], 5,9,0.001)

fig2 = plt.figure(1,figsize=(cm2inch(35), cm2inch(20)))
ax3 = plt.subplot(121)
ax4 = plt.subplot(122)

NiceGraph2Dwticks(ax3, 'Real Angle [rad]', 'Normalized Hinge Energy with offset', [-3.1416, 0.],[0,2.10], 5, 15, 0.01)
NiceGraph2D(ax4, 'kHinge/kEdge', 'EHinge/kHinge at -pi rad')


colors = cm.Set1(np.linspace(0, 1, nruns))#['#36438B','#CB5F3D','#C2C93C']

for krun,c in zip(np.arange(nruns),colors):
    ax1.plot(realAngle[(nruns-(krun+1))*datapoints:(nruns-krun)*datapoints],enTotal[(nruns-(krun+1))*datapoints:(nruns-krun)*datapoints]+offset1*(nruns-1-krun), '.',c=c)
    ax1.plot(realAngle[(nruns-(krun+1))*datapoints:(nruns-krun)*datapoints],enTotal[(nruns-(krun+1))*datapoints:(nruns-krun)*datapoints]+offset1*(nruns-1-krun), '-',c=c)
    ax1.plot(realAngle[(nruns-(krun+1))*datapoints:(nruns-krun)*datapoints],enHinge[(nruns-(krun+1))*datapoints:(nruns-krun)*datapoints]+offset1*(nruns-1-krun), '-',c=c, alpha = 0.5)
    ax1.plot(realAngle[(nruns-(krun+1))*datapoints:(nruns-krun)*datapoints],enEdge[(nruns-(krun+1))*datapoints:(nruns-krun)*datapoints]+offset1*(nruns-1-krun), '-',c=c, alpha = 0.5)
    ax2.plot(realAngle[(nruns-(krun+1))*datapoints:(nruns-krun)*datapoints],errorAngle[(nruns-(krun+1))*datapoints:(nruns-krun)*datapoints]+offset2*(nruns-1-krun), '.',label = kAngleConst[(nruns-(krun+1))*datapoints],c=c)
    ax2.plot(realAngle[(nruns-(krun+1))*datapoints:(nruns-krun)*datapoints],errorAngle[(nruns-(krun+1))*datapoints:(nruns-krun)*datapoints]+offset2*(nruns-1-krun), '-',c=c)
    if krun != 0:
        ax3.plot(realAngle[krun*datapoints:(krun+1)*datapoints],normEnHinge[krun*datapoints:(krun+1)*datapoints]+offset3*(nruns-1-krun), '.',label = kAngleConst[krun*datapoints], c = colors[-krun-1])
        ax3.plot(realAngle[krun*datapoints:(krun+1)*datapoints],normEnHinge[krun*datapoints:(krun+1)*datapoints]+offset3*(nruns-1-krun), '-', c = colors[-krun-1])

ax4.plot(kvector[1:-1], enHingekHinge[1:-1], '.', c='#36438B')
ax4.plot(kvector[1:-1], enHingekHinge[1:-1], '-', c='#36438B')

leg = ax2.legend(fontsize = 15, frameon = True, framealpha=0.5)
leg2 = ax3.legend(fontsize = 15, frameon = True, framealpha=0.5)

fig.tight_layout()
fig2.tight_layout()
fig.show()
fig2.show()
fig.savefig(dirname+file_name[:-4]+'.png', transparent = True)
fig2.savefig(dirname+file_name[:-4]+'2.png', transparent = True)
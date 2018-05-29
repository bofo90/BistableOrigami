# -*- coding: utf-8 -*-
"""
Created on Tue May 29 10:21:17 2018

@author: iniguez
"""

import EnergyLandscape as el
import matplotlib as matl
import matplotlib.pyplot as plt
import numpy as np

folder_name1 = "Results/triangular prism/active-set/energy/17-May-2018_EnergyAllAngles_Kdep_8_3\kh0.075_kta100.000_ke3.000_kf100.000"
folder_name2 = "Results/triangular prism/active-set/energy/03-Apr-2018_EnergyAllAngles_8_3\kh0.001_kta100.000_ke3.000_kf100.000"
inverted = True
maxEnergy = 4.5
[angle, ek_1] = el.EnergyLandscape(folder_name1, inverted, maxEnergy, plot = False)
[angle, ek_2] = el.EnergyLandscape(folder_name2, inverted, maxEnergy, plot = False)
ek_0 = np.zeros(np.size(ek_1))+7
ek_0[13:] = 50000

fig1 = plt.figure(3,figsize=(el.cm2inch(24.1), el.cm2inch(20)))
ax1 = plt.subplot(111)
el.NiceGraph2D(ax1, 'Hinge 3 [$\pi$ rad]', 'Energy/$k_H$ [a.u.]', [0,5],[1,2000])
ax1.xaxis.set_major_formatter(matl.ticker.FormatStrFormatter('%.2g'))

ax1.plot(angle,ek_0, label = r'$0$')
ax1.plot(angle,ek_2/0.001, label = '~0.001')
ax1.plot(angle,ek_1/0.075, label = '~0.01')


ax1.set_yscale('log')
ax1.yaxis.set_major_formatter(matl.ticker.FormatStrFormatter('%g'))

lgnd = ax1.legend(loc = 2)

#lgnd = fig1.legend(fig1, loc=8, mode="expand", ncol= 5)
#for handle in lgnd.legendHandles:
#    handle.set_sizes([500])
for label in lgnd.get_texts():
    label.set_fontsize('11')
    label.set_color('0.2')
lgnd.get_frame().set_alpha(0)

fig1.tight_layout()
fig1.show()
fig1.savefig("Results/triangular prism/active-set/energy" + '/EnergyHingeKs.pdf', transparent = True)

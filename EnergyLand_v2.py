# -*- coding: utf-8 -*-
"""
Created on Tue May 21 18:04:51 2019

@author: iniguez
"""

import pandas as pd
import numpy as np
import configparser
import os.path

folder_name = "Results/SingleVertex3/active-set/energy/21-May-2019_temp/kh0.000_kta100.000_ke1.000_kf100.000"

file_name1 = "/EnergyData.csv" 
file_name2 = "/Hinges.csv"
file_name3 = "/PosStad.csv"
file_name4 = "/Angles.csv"

dataEnergy = pd.read_csv(folder_name+file_name1)
dataEnergy['TotalEnergy'] = dataEnergy['EdgeEnergy']+dataEnergy['DiagonalEnergy']+dataEnergy['HingeEnergy']

dataVar = pd.read_csv(folder_name+file_name2)

dataAngles = np.loadtxt(folder_name+file_name4,skiprows=1, delimiter = ',', dtype = np.float64)
dataAngles = np.delete(dataAngles, 0, 1)

TotSimul = dataVar.shape[0]
IterPerSimulEnergy = np.int(dataEnergy.shape[0]/TotSimul)

print(dataEnergy['Flags'].value_counts())
exfl = dataEnergy.Flags.values.reshape(TotSimul,IterPerSimulEnergy)
flagmask = (exfl !=1) & (exfl !=5) & (exfl !=5)
flagmask = ~flagmask.any(axis= 1)

dataVar['Mask'] = flagmask
dataVar = dataVar.join(dataEnergy[['Hinge Number','TotalEnergy']][IterPerSimulEnergy-2::IterPerSimulEnergy].set_index('Hinge Number'), on = 'HingeNumber')
dataVar.set_index(['Kappa','TargetAngle'], inplace = True)
dataVar.sort_index(inplace=True)

dataEnergy.set_index('Hinge Number', inplace = True)
dataEnergy.sort_index(inplace=True)


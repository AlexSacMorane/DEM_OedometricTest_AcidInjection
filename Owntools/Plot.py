# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This file contains the different functions to plot used in the simulation.
"""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import os

#-------------------------------------------------------------------------------

def Plot_PSD(namefile, dict_sample):
    """
    Plot the PSD.

        Input :
            a file name (a str)
            a sample dictionnary (a dict)
        Output :
            Nothing, but a png file is generated (a file)
    """
    plt.figure(1,figsize=(16,9))
    plt.plot(dict_sample['L_radius'], dict_sample['L_cumulative_mass'])
    plt.title('Particle Size Distribution')
    plt.xlabel('Radius (µm)')
    plt.ylabel('Cumulative mass retained (kg)')
    plt.savefig(namefile)
    plt.close(1)

#-------------------------------------------------------------------------------

def Plot_DEM_tracker(dict_tracker):
    '''
    Plot the trackers of the DEM step.

        Input :
            a tracker dictionnary (a dict)
        Output :
            Nothing but a .png file is generated (file)
    '''
    #look for the name of the new plot
    template_name = 'Debug/DEM_tracker/PFDEM_'
    j = 1
    plotpath = Path(template_name+str(j)+'.png')
    while plotpath.exists():
        j = j + 1
        plotpath = Path(template_name+str(j)+'.png')
    name = template_name+str(j)+'.png'


    plt.figure(1,figsize=(16,9))

    plt.subplot(221)
    plt.plot(dict_tracker['Ecin'])
    plt.title('Kinetic energy')

    plt.subplot(222)
    plt.plot(dict_tracker['y_box_max_DEM'])
    plt.title('Upper wall position')

    plt.subplot(223)
    plt.plot(dict_tracker['Force_on_upper_wall'])
    plt.title('Force on the upper wall')

    plt.savefig(name)
    plt.close(1)

#-------------------------------------------------------------------------------

def Plot_yboxmax(dict_tracker):
    '''
    Plot the evolution of the upper wall position.

        Input :
            a tracker dictionnary (a dict)
        Output :
            Nothing but a .png file is generated (file)
    '''
    plt.figure(1,figsize=(16,9))
    plt.plot(dict_tracker['L_t'], dict_tracker['L_y_box_max'], marker ='x')
    plt.title('Upper wall position')
    plt.savefig('Debug/Evolution_yboxmax.png')
    plt.close(1)

#-------------------------------------------------------------------------------

def Plot_porosity(dict_tracker):
    '''
    Plot the evolution of the sample porosity (= grain surface / box surface).

        Input :
            a tracker dictionnary (a dict)
        Output :
            Nothing but a .png file is generated (file)
    '''
    plt.figure(1,figsize=(16,9))
    plt.plot(dict_tracker['L_t'], dict_tracker['L_porosity'], marker ='x')
    plt.title('Evolution of the sample porosity')
    plt.savefig('Debug/Evolution_porosity.png')
    plt.close(1)

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

def Plot_mass(namefile, dict_tracker):
    """
    Plot the evolution of the mass distribution.

        Intput :
            a file name (a str)
            a tracker dictionnary (a dict)
        Output :
            Nothing, but a png file is generated (a file)
    """
    plt.figure(1,figsize=(16,9))
    plt.plot(dict_tracker['L_mass'], label ='Grains mass')
    plt.plot(dict_tracker['L_mass_dissolved'], label ='Mass dissolved')
    plt.xlabel('Dissolution iteration (-)')
    plt.ylabel('Mass (kg)')
    plt.legend()
    plt.savefig(namefile)
    plt.close(1)

#-------------------------------------------------------------------------------

def Plot_custom(namefile, L_data_x, L_data_y, L_label_name, L_label):
    """
    Plot a custom curve.

        Intput :
            a file name (a str)
            a list of data x (a list)
            a list of data y (a list)
            a list of label (a list)
            a list of label (a list)
        Output :
            Nothing, but a png file is generated (a file)
    """
    plt.figure(1,figsize=(16,9))
    plt.plot(L_data_x, L_data_y)
    if 'title' in L_label_name:
        plt.title(L_label[L_label_name.index('title')])
    if 'xlabel' in L_label_name:
        plt.xlabel(L_label[L_label_name.index('xlabel')])
    if 'ylabel' in L_label_name:
        plt.ylabel(L_label[L_label_name.index('ylabel')])
    plt.savefig(namefile)
    plt.close(1)

#-------------------------------------------------------------------------------

def Plot_DEM_trackers(namefile, Force_tracker, Ecin_tracker, Zmax_tracker, s_top_tracker):
    """
    Plot trackers from DEM during loading.

        Input :
            a namefile (a str)
            four lists of tracker (lists)
        Output :
            Nothing, but a .png file is generated (a file)
    """
    plt.figure(1,figsize=(16,9))

    plt.subplot(221)
    plt.title('Force applied on particles (µN)')
    plt.plot(Force_tracker)

    plt.subplot(222)
    plt.title('Kinetic energy of particles (10-12 J)')
    plt.plot(Ecin_tracker)

    plt.subplot(223)
    plt.title('Upper wall position (µm)')
    plt.plot(Zmax_tracker)

    plt.subplot(224)
    plt.title('Vertical stress on top wall (µN)')
    plt.plot(s_top_tracker)

    plt.savefig(namefile)
    plt.close(1)

#-------------------------------------------------------------------------------

def Plot_ir_trackers(namefile, Ecin_tracker, Ratio_Displacement_MeanRadius_tracker, n_contact_tracker):
    """
    Plot trackers from DEM during the radius expansion.

        Input :
            a namefile (a str)
            trackers (lists)
        Output :
            Nothing, but a .png file is generated (a file)
    """
    plt.figure(1,figsize=(16,9))

    plt.subplot(131)
    plt.title('Kinetic energy of particles (10-12 J)')
    plt.plot(Ecin_tracker)

    plt.subplot(132)
    plt.title('Mean Displacement / Mean Radius (-)')
    plt.plot(Ratio_Displacement_MeanRadius_tracker)

    plt.subplot(133)
    plt.title('Number of contacts (-)')
    plt.plot(n_contact_tracker)

    plt.savefig(namefile)
    plt.close(1)

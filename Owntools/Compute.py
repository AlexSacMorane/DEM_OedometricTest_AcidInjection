# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This file contains the different functions used to compute parameters or variables in the simulation.
"""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

import numpy as np
import math

#-------------------------------------------------------------------------------

def Compute_PSD(dict_sample):
    """
    Compute the PSD of the sample.

        Input :
            a sample dictionnary (a dict)
        Output :
            Nothing but the sample dictionnary gets new attributes describing the PSD
    """
    #extract data
    L_mass = []
    L_radius = []
    for grain in dict_sample['L_g']:
        L_mass.append(grain.mass)
        L_radius.append(grain.radius)

    #compute cumulative mass + sort list of radius
    L_radius_sorted = []
    L_cumulative_mass = []
    while L_radius != []:
        i_r_min = L_radius.index(min(L_radius))
        L_radius_sorted.append(L_radius[i_r_min])
        if L_cumulative_mass == []:
            L_cumulative_mass.append(L_mass[i_r_min])
        else :
            L_cumulative_mass.append(L_cumulative_mass[-1] + L_mass[i_r_min])
        L_radius.pop(i_r_min)
        L_mass.pop(i_r_min)

    #update dictionnaries
    dict_sample['L_radius'] = L_radius
    dict_sample['L_cumulative_mass'] = L_cumulative_mass

#-------------------------------------------------------------------------------

def Compute_compacity(dict_sample):
    '''
    Compute the compacity (grain volume / box volume) of the sample.

        Input :
            a sample dictionnary (a dict)
        Output :
            Nothing but the dictionnary gets an updated value for the compacity (a float)
    '''
    Vg = 0
    for grain in dict_sample['L_g']:
        Vg = Vg + grain.volume
    Vb = math.pi*dict_sample['D_oedo']**2/4*(dict_sample['z_box_max']-dict_sample['z_box_min'])

    #update element in dict
    dict_sample['compacity'] = Vg/Vb

#-------------------------------------------------------------------------------

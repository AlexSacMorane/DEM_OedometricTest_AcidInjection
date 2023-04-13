# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This file contains the different functions used in the simulation.
"""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

import os
import math
import numpy as np
from pathlib import Path

#-------------------------------------------------------------------------------

def error_on_zmax_f(dz,overlap_L,k_L,Force_target) :
    """
    Compute the function f to control the upper wall. It is the difference between the force applied and the target value.

        Input :
            an increment of the upper wall position (a float)
            a list of overlap for contact between grain and upper wall (a list)
            a list of spring for contact between grain and upper wall (a list)
            a confinement force (a float)
        Output :
            the difference between the force applied and the confinement (a float)
    """
    f = Force_target
    for i in range(len(overlap_L)):
        f = f - k_L[i]*(max(overlap_L[i]-dz,0))**(3/2)
    return f

#-------------------------------------------------------------------------------

def error_on_zmax_df(dz,overlap_L,k_L) :
    """
    Compute the derivative function df to control the upper wall (error_on_zmax_f()).

        Input :
            an increment of the upper wall position (a float)
            a list of overlap for contact between grain and upper wall (a list)
            a list of spring for contact between grain and upper wall (a list)
        Output :
            the derivative of error_on_zmax_f() (a float)
    """
    df = 0
    for i in range(len(overlap_L)):
        df = df + 3/2*k_L[i]*(max(overlap_L[i]-dz,0))**(1/2)
    return df

#-------------------------------------------------------------------------------

def Control_z_max_NR(dict_sample, dict_sollicitation):
    """
    Control the upper wall to apply force.

    A Newton-Raphson method is applied to verify the confinement.
        Input :
            a sample dictionnary (a dict)
            a sollicitations dictionnary (a dict)
        Output :
            Nothing, but the sample dictionnary is updated, concerning the upper wall position and force applied (two floats)
    """
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
    #load data needed
    Force_target = dict_sollicitation['Vertical_Confinement_Force']
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

    F = 0
    overlap_L = []
    k_L = []
    for contact in dict_sample['L_contact_gw']:
        if contact.nature == 'gwy_max':
            F = F + contact.Fwg_n
            overlap_L.append(contact.overlap)
            k_L.append(contact.k)
            #compute force applied, save contact overlap and spring

    if overlap_L != []:
        i_NR = 0
        dz = 0
        ite_criteria = True
        if -0.01*Force_target<error_on_zmax_f(dz,overlap_L,k_L,Force_target) and error_on_zmax_f(dz,overlap_L,k_L,Force_target)<0.01*Force_target:
            ite_criteria = False
        while ite_criteria :
            i_NR = i_NR + 1
            dz = dz - error_on_zmax_f(dz,overlap_L,k_L,Force_target)/error_on_zmax_df(dz,overlap_L,k_L)
            if i_NR > 100:
                ite_criteria = False
            if -0.01*Force_target<error_on_zmax_f(dz,overlap_L,k_L,Force_target) and error_on_zmax_f(dz,overlap_L,k_L,Force_target)<0.01*Force_target:
                ite_criteria = False
        dict_sample['z_box_max'] = dict_sample['z_box_max'] + dz

    else :
        #if there is no contact with the upper wall, the wall is reset
        dict_sample['z_box_max'] = Reset_z_max(dict_sample['L_g'],Force_target)

    #Update dict
    dict_sollicitation['Force_on_upper_wall'] = F

#-------------------------------------------------------------------------------

def Reset_z_max(L_g,Force):
    """
    The upper wall is located as a single contact verify the target value.

        Input :
            the list of grains (a list)
            the confinement force (a float)
        Output :
            the upper wall position (a float)
    """
    z_max = None
    id_grain_max = None
    for id_grain in range(len(L_g)):
        grain = L_g[id_grain]
        z_max_grain = grain.center[2] + grain.radius

        if z_max != None and z_max_grain > z_max:
            z_max = z_max_grain
            id_grain_max = id_grain
        elif z_max == None:
            z_max = z_max_grain
            id_grain_max = id_grain

    factor = 5
    k = factor*4/3*L_g[id_grain_max].y/(1-L_g[id_grain_max].nu*L_g[id_grain_max].nu)*math.sqrt(L_g[id_grain_max].radius)
    z_max = z_max - (Force/k)**(2/3)

    return z_max

#-------------------------------------------------------------------------------

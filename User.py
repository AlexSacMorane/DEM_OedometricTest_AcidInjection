# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This is the file where the user can change the different parameters for the simulation.
"""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

import math
import numpy as np
from pathlib import Path

#Own functions and classes
import Grain

#-------------------------------------------------------------------------------
#User
#-------------------------------------------------------------------------------

def All_parameters():
    '''
    This function is called in main.py to have all the parameters needed in the simulation.

        Input :
            Nothing
        Output :
            an algorithm dictionnary (a dictionnary)
            a material dictionnary (a dictionnary)
            a sample dictionnary (a dictionnary)
            a sollicitation dictionnary (a dictionnary)
    '''
    #---------------------------------------------------------------------------
    #Geometry parameters

    N_grain = 500 #total number of grains

    #a normal law is assumed for the PSD
    R_50 = 350 #µm expectation
    sigma_psd = 25 #µm standard deviation
    #95% are in [R_50-2*sigma_psd, R_50+2*sigma_psd]
    #99,7% is in [R_50-3*sigma_psd, R_50+3*sigma_psd]

    #write dict
    dict_geometry = {
    'N_grain' : N_grain,
    'R_50' : R_50,
    'sigma_psd' : sigma_psd
    }

    #---------------------------------------------------------------------------
    #Sample parameters

    #Box définition
    z_box_min = 0 #µm
    D_oedo = (dict_geometry['N_grain']*11)**(1/3)*dict_geometry['R_50'] #µm the diameter of the oedometer

    dict_sample = {
    'z_box_min' : z_box_min,
    'D_oedo' : D_oedo
    }

    #---------------------------------------------------------------------------
    #Material parameters

    #DEM
    Y = 70*(10**9)*(10**6)*(10**(-12)) #Young Modulus µN/µm2
    nu = 0.3 #Poisson's ratio
    rho = 2500*10**(-6*3) #density kg/µm3
    mu_friction_gg = 0.5 #grain-grain
    mu_friction_gw = 0.5 #grain-wall
    coeff_restitution = 0.2 #1 is perfect elastic

    dict_material = {
    'Y' : Y,
    'nu' : nu,
    'rho' : rho,
    'mu_friction_gg' : mu_friction_gg,
    'mu_friction_gw' : mu_friction_gw,
    'coeff_restitution' : coeff_restitution
    }

    #---------------------------------------------------------------------------
    #Algorithm parameters

    #DEM parameters
    dt_DEM_crit = math.pi*R_50/(0.16*nu+0.88)*math.sqrt(rho*(2+2*nu)/Y) #s critical time step from O'Sullivan 2011
    dt_DEM = dt_DEM_crit/8 #s time step during DEM simulation
    factor_neighborhood = 2.5 #margin to detect a grain into a neighborhood
    i_update_neighborhoods = 50 #the frequency of the update of the neighborhood of the grains and the walls
    #Stop criteria of the DEM
    i_DEM_stop = 1500 #maximum iteration for one DEM simulation
    Ecin_ratio = 0.0002
    n_window_stop = 100
    dz_box_max_stop = 35

    #List of plot to do
    Debug = True #plot configuration before and after DEM simulation
    Debug_DEM = False #plot configuration inside DEM
    i_print_plot = 300 #frenquency of the print and plot (if Debug_DEM) in DEM step
    #
    L_flag_plot = []

    #Save the simulation
    SaveData = False #Save data or not
    foldername = 'Data_Oedo_Acid' #name of the folder where data are saved
    template = 'Run' #template of the name of the simulation
    if SaveData :
        i_run = 1
        folderpath = Path('../'+foldername+'/'+template+'_'+str(i_run))
        while folderpath.exists():
            i_run = i_run + 1
            folderpath = Path('../'+foldername+'/'+template+'_'+str(i_run))
        namefile = template+'_'+str(i_run)
    else :
        namefile = template

    dict_algorithm = {
    'Debug' : Debug,
    'Debug_DEM' : Debug_DEM,
    'i_print_plot' : i_print_plot,
    'factor_neighborhood' : factor_neighborhood,
    'SaveData' : SaveData,
    'namefile' : namefile,
    'dt_DEM' : dt_DEM,
    'i_update_neighborhoods': i_update_neighborhoods,
    'i_DEM_stop' : i_DEM_stop,
    'Ecin_ratio' : Ecin_ratio,
    'n_window_stop' : n_window_stop,
    'dz_box_max_stop' : dz_box_max_stop,
    'foldername' : foldername,
    'L_flag_plot' : L_flag_plot,
    }

    #---------------------------------------------------------------------------
    #External sollicitation parameters

    gravity = 0
    Vertical_Confinement_Linear_Force = Y*2*R_50/100 #µN/µm used to compute the Vertical_Confinement_Force
    Vertical_Confinement_Force = Vertical_Confinement_Linear_Force*(D_oedo) #µN

    dict_sollicitation = {
    'gravity' : gravity,
    'Vertical_Confinement_Force' : Vertical_Confinement_Force
    }

    #---------------------------------------------------------------------------
    #Initial condition parameters

    n_generation = 3 #number of grains generation
    factor_ymax_box = 3 #margin to generate grains
    N_test_max = 5000 # maximum number of tries to generate a grain without overlap
    i_DEM_stop_IC = 3000 #stop criteria for DEM during IC
    Debug_DEM_IC = True #plot configuration inside DEM during IC
    i_print_plot_IC = 300 #frenquency of the print and plot (if Debug_DEM_IC) for IC
    dt_DEM_IC = dt_DEM_crit/6 #s time step during IC
    Ecin_ratio_IC = 0.001
    factor_neighborhood_IC = 1.5 #margin to detect a grain into a neighborhood
    i_update_neighborhoods_gen = 50 #the frequency of the update of the neighborhood of the grains and the walls during IC generations
    i_update_neighborhoods_com = 200 #the frequency of the update of the neighborhood of the grains and the walls during IC combination

    #write dict
    dict_ic = {
    'n_generation' : n_generation,
    'i_update_neighborhoods_gen': i_update_neighborhoods_gen,
    'i_update_neighborhoods_com': i_update_neighborhoods_com,
    'factor_ymax_box' : factor_ymax_box,
    'i_DEM_stop_IC' : i_DEM_stop_IC,
    'Debug_DEM' : Debug_DEM_IC,
    'dt_DEM_IC' : dt_DEM_IC,
    'Ecin_ratio_IC' : Ecin_ratio_IC,
    'i_print_plot_IC' : i_print_plot_IC,
    'factor_neighborhood_IC' : factor_neighborhood_IC,
    'N_test_max' : N_test_max
    }

    #---------------------------------------------------------------------------

    return dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitation

#-------------------------------------------------------------------------------

def Criteria_StopSimulation(dict_algorithm):
    '''
    Define a stop criteria for the PFDEM simulation.

    The simulation stops if the number of iterations reaches a maximum value.

        Input :
            an algorithm dictionnary (a dictionnary)
        Output :
            The result depends on the fact if the stop criteria is reached or not (a bool)
    '''
    Criteria_Verified = False
    if dict_algorithm['i_PFDEM'] >= dict_algorithm['n_t_PFDEM']:
        Criteria_Verified = True
    return Criteria_Verified

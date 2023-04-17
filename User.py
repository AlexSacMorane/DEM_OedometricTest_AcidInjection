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

    N_grain = 1000 #total number of grains

    #a normal law is assumed for the PSDs

    #PSD 1
    R_50_1 = 300 #µm expectation
    sigma_psd_1 = 50 #µm standard deviation
    #95% are in [R_50-2*sigma_psd, R_50+2*sigma_psd]
    #99,7% is in [R_50-3*sigma_psd, R_50+3*sigma_psd]

    #PSD 2
    R_50_2 = 80 #µm expectation
    sigma_psd_2 = 10 #µm standard deviation
    #95% are in [R_50-2*sigma_psd, R_50+2*sigma_psd]
    #99,7% is in [R_50-3*sigma_psd, R_50+3*sigma_psd]

    #distribution between PSD 1 and PSD 2
    mass_ratio_1_1and2 = 1

    #recompute r_mean with 1 + 2
    number_ratio_1_1and2 = (mass_ratio_1_1and2*R_50_2**3)/(mass_ratio_1_1and2*R_50_2**3+(1-mass_ratio_1_1and2)*R_50_1**3)
    R_50 = number_ratio_1_1and2*R_50_1 + (1-number_ratio_1_1and2)*R_50_2

    #target for compacity e = Vg/Vb
    e_target = 0.7 #for IC with grain radius

    #write dict
    dict_geometry = {
    'N_grain' : N_grain,
    'R_50_1' : R_50_1, #not used for the moment
    'sigma_psd_1' : sigma_psd_1,
    'R_50_2' : R_50_2,
    'sigma_psd_2' : sigma_psd_2,
    'mass_ratio_1_1and2' : mass_ratio_1_1and2,
    'R_50' : R_50,
    'e_target' : e_target
    }

    #---------------------------------------------------------------------------
    #Sample parameters

    #Box définition
    z_box_min = 0 #µm
    D_oedo = (dict_geometry['N_grain']*14)**(1/3)*dict_geometry['R_50'] #µm the diameter of the oedometer

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
    mu_friction_gw = 0 #grain-wall
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

    #stop criteria
    n_dissolution = 50 #number of dissolution increment
    f_mass0_dissolved_mas = 0.5 #maximum of the initial mass dissolved

    #DEM parameters
    dt_DEM_crit = math.pi*dict_geometry['R_50']/(0.16*nu+0.88)*math.sqrt(rho*(2+2*nu)/Y) #s critical time step from O'Sullivan 2011
    dt_DEM = dt_DEM_crit/6 #s time step during DEM simulation
    factor_neighborhood = 1.9 #margin to detect a grain into a neighborhood
    i_update_neighborhoods = 200 #the frequency of the update of the neighborhood of the grains and the walls
    #Stop criteria of the DEM
    i_DEM_stop = 4000 #maximum iteration for one DEM simulation
    Ecin_ratio = 0.0001
    n_window = 100

    #List of plot to do
    Debug_DEM = True #plot configuration
    i_print_plot = 400 #frenquency of the print and plot (if Debug_DEM) in DEM step

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
    'n_dissolution' : n_dissolution,
    'f_mass0_dissolved_mas' : f_mass0_dissolved_mas,
    'dt_DEM' : dt_DEM,
    'factor_neighborhood' : factor_neighborhood,
    'i_update_neighborhoods': i_update_neighborhoods,
    'i_DEM_stop' : i_DEM_stop,
    'Ecin_ratio' : Ecin_ratio,
    'n_window' : n_window,
    'Debug_DEM' : Debug_DEM,
    'i_print_plot' : i_print_plot,
    'SaveData' : SaveData,
    'foldername' : foldername,
    'namefile' : namefile,
    }

    #---------------------------------------------------------------------------
    #External sollicitation parameters

    gravity = 0 #µm/s2
    Vertical_Confinement_Surface_Force = 100*10**-3 #µN/µm2 used to compute the Vertical_Confinement_Force
    Vertical_Confinement_Force = Vertical_Confinement_Surface_Force*(math.pi*D_oedo**2/4) #µN
    f_R50_0_dissolved = 0.005 #fraction of the initial mean radius dissolved

    dict_sollicitation = {
    'gravity' : gravity,
    'Vertical_Confinement_Force' : Vertical_Confinement_Force,
    'f_R50_0_dissolved' : f_R50_0_dissolved
    }

    #---------------------------------------------------------------------------
    #Initial condition parameters

    #method to generate ic
    method_ic = 'Deposition' #IncreaseRadius or Deposition

    n_generation = 5 #number of grains generation
    factor_zmax_box = 2 #margin to generate grains (Deposition)
    N_test_max = 5000 # maximum number of tries to generate a grain without overlap
    i_DEM_stop_IC = 5000 #stop criteria for DEM during IC
    Debug_DEM_IC = True #plot configuration inside DEM during IC
    i_print_plot_IC = 300 #frenquency of the print and plot (if Debug_DEM_IC) for IC
    dt_DEM_IC = dt_DEM_crit/5 #s time step during IC
    n_window = 100 #window to detect the steady-state
    Ecin_ratio_IC = 0.0005 #criteria on kinetic energy to detect the steady-state
    factor_neighborhood_IC = 1.5 #margin to detect a grain into a neighborhood
    i_update_neighborhoods_gen = 50 #the frequency of the update of the neighborhood of the grains and the walls during IC generations
    i_update_neighborhoods_com = 200 #the frequency of the update of the neighborhood of the grains and the walls during IC combination
    gravity = 100*N_grain/n_generation*(factor_zmax_box)*dict_geometry['R_50']**3/D_oedo**2/i_DEM_stop_IC**2/dt_DEM_IC**2 #apply only in the ic phase with one generation µm/s2

    #write dict
    dict_ic = {
    'method_ic' : method_ic,
    'n_generation' : n_generation,
    'factor_zmax_box' : factor_zmax_box,
    'N_test_max' : N_test_max,
    'i_DEM_stop_IC' : i_DEM_stop_IC,
    'Debug_DEM' : Debug_DEM_IC,
    'i_print_plot_IC' : i_print_plot_IC,
    'dt_DEM_IC' : dt_DEM_IC,
    'n_window' : n_window,
    'Ecin_ratio_IC' : Ecin_ratio_IC,
    'factor_neighborhood_IC' : factor_neighborhood_IC,
    'i_update_neighborhoods_gen': i_update_neighborhoods_gen,
    'i_update_neighborhoods_com': i_update_neighborhoods_com,
    'gravity' : gravity
    }

    #---------------------------------------------------------------------------

    return dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitation

#-------------------------------------------------------------------------------

def Criteria_StopSimulation(dict_algorithm, dict_tracker):
    '''
    Define a stop criteria for the PFDEM simulation.

    The simulation stops if the number of iterations reaches a maximum value.

        Input :
            an algorithm dictionnary (a dictionnary)
            a tracker dictionnary (a dictionnary)
        Output :
            The result depends on the fact if the stop criteria is reached or not (a bool)
    '''
    Criteria_Verified = False
    if dict_algorithm['i_dissolution'] >= dict_algorithm['n_dissolution']:
        Criteria_Verified = True
    if dict_tracker['L_mass_dissolved'][-1]/dict_tracker['L_mass'][0] >= dict_algorithm['f_mass0_dissolved_mas']:
        Criteria_Verified = True
    return Criteria_Verified

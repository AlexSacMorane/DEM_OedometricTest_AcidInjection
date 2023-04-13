# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This is the main file.
"""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

from pathlib import Path
import numpy as np
import os
import shutil
import math

#Own functions and classes
import Grain
import Contact_gg
import Contact_gw
import Owntools
import Owntools.Compute
import Owntools.Plot
import Owntools.Save
import Owntools.Write
import Create_IC
import Create_IC.Grain_ic
import Create_IC.Contact_gg_ic
import Create_IC.Contact_gw_ic
import User
import Report

#-------------------------------------------------------------------------------

def open_simulation():
    """
    Open a simulation.

        Input :
            Nothing
        Output :
            an algorithm dictionnary (a dict)
            a geometry dictionnary (a dict)
            an initial condition dictionnary (a dict)
            a material dictionnary (a dict)
            a sample dictionnary (a dict)
            a sollicitation dictionnary (a dict)
            a simulation report (a report)
    """
    #Plan simulation
    if Path('Debug').exists():
        shutil.rmtree('Debug')
    os.mkdir('Debug')
    if Path('Dicts').exists():
        shutil.rmtree('Dicts')
    os.mkdir('Dicts')

    #create a simulation report
    simulation_report = Report.Report('Debug/Report')

    #general parameters
    dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitation = User.All_parameters()
    if dict_algorithm['SaveData']:
        if not Path('../'+dict_algorithm['foldername']).exists():
            os.mkdir('../'+dict_algorithm['foldername'])
        #tempo save of the user file
        shutil.copy('User.py','../'+dict_algorithm['foldername']+'/User_'+dict_algorithm['namefile']+'_tempo.txt')

    #prepare plot
    if dict_algorithm['Debug_DEM'] or dict_ic['Debug_DEM']:
        os.mkdir('Debug/Configuration')
        if dict_ic['Debug_DEM'] :
            os.mkdir('Debug/Configuration/Init')
        if dict_algorithm['Debug_DEM'] :
            os.mkdir('Debug/Configuration/Main')

    return dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitation, simulation_report

#-------------------------------------------------------------------------------

def create_ic(dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitation, simulation_report):
    """
    Create an initial condition.

        Input :
            a geometry dictionnary (a dict)
            an initial condition dictionnary (a dict)
            a material dictionnary (a dict)
            a sample dictionnary (a dict)
            a sollicitation dictionnary (a dict)
            a simulation report (a Report)
        Output :
            Nothing, but the dictionnaries are updated
    """
    simulation_report.tic_tempo()
    #create the initial configuration
    Create_IC.LG_tempo(dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitation, simulation_report)

    #Convert the tempo grain to real grain
    Create_IC.From_LG_tempo_to_usable(dict_ic, dict_sample)

    #compute and plot the PSD
    Owntools.Compute.Compute_PSD(dict_sample)
    Owntools.Plot.Plot_PSD('Debug/PSD.png', dict_sample)

    #save data
    Owntools.Save.save_dicts_ic('Dicts/save_ic', dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitation, simulation_report)

    simulation_report.tac_tempo('Initialisation')

#-------------------------------------------------------------------------------

def main_simulation(dict_algorithm, dict_geometry, dict_material, dict_sample, dict_sollicitation, dict_tracker, simulation_report):
    """
    It is the main simulation.

    Once a steady-state is detected, dissolution of the material occurs. It is repeated until the stop criteria is reached.

        Input :
            an algorithm dictionnary (a dict)
            a geometry dictionnary (a dict)
            a material dictionnary (a dict)
            a sample dictionnary (a dict)
            a sollicitation dictionnary (a dict)
            a tracker dictionnary (a dict)
            a simulation report (a report)
        Output :
            Nothing, but dictionnaries are updated
    """
    # Preparation and add elements in dicts
    dict_algorithm['i_dissolution'] = 0
    dict_algorithm['i_DEM'] = 0
    dict_sample['L_contact_gw'] = []
    dict_sample['L_contact_gw_ij'] = []
    dict_sample['id_contact_gw'] = 0
    dict_sample['L_contact'] = []
    dict_sample['L_contact_ij'] = []
    dict_sample['id_contact'] = 0

    while not User.Criteria_StopSimulation(dict_algorithm, dict_tracker):
        simulation_report.tic()
        dict_algorithm['i_dissolution'] = dict_algorithm['i_dissolution'] + 1

        dissolve_material(dict_geometry, dict_sample, dict_sollicitation, dict_tracker, simulation_report)

        simulation_report.write_and_print('\nStep '+str(dict_algorithm['i_dissolution'])+'/'+str(dict_algorithm['n_dissolution'])+' : percentage of initial mass dissolved '+str(round(dict_tracker['L_mass_dissolved'][-1]/dict_tracker['L_mass'][0],2))+'/'+str(dict_algorithm['f_mass0_dissolved_mas'])+'\n\n',\
                                          '\nStep '+str(dict_algorithm['i_dissolution'])+'/'+str(dict_algorithm['n_dissolution'])+' : percentage of initial mass dissolved '+str(round(dict_tracker['L_mass_dissolved'][-1]/dict_tracker['L_mass'][0],2))+'/'+str(dict_algorithm['f_mass0_dissolved_mas'])+'\n')

        simulation_report.tac('New dissolution step')

#-------------------------------------------------------------------------------

def dissolve_material(dict_geometry, dict_sample, dict_sollicitation, dict_tracker, simulation_report):
    """
    Dissolve the granular material.

        Input :
            a geometry dictionnary (a dict)
            a sample dictionnary (a dict)
            a sollicitation dictionnary (a dict)
            a tracker dictionnary (a dict)
            a simulation report (a report)
        Output :
            Nothing, but dictionnaries are updated
    """
    for i_grain in range(len(dict_sample['L_g'])):
        grain = dict_sample['L_g'][i_grain]
        #dissolution
        grain.radius = grain.radius - dict_geometry['R_50']*dict_sollicitation['f_R50_0_dissolved']
        if grain.radius > 0 :
            #update properties
            grain.volume = 4/3*math.pi*grain.radius**3
            grain.mass = grain.rho*grain.volume
            self.inertia = 2/5*self.mass*grain.radius**2
        else : #no grain anymore
            grain.radius = grain.radius + dict_geometry['R_50']*dict_sollicitation['f_R50_0_dissolved']
            #if grain is deleted, be carefull with :
            #   - delete contact / contact_gw
            #   - update L_contact_ij / L_contact_gw_ij as indentation change
            #for the moment, just undissolve the grain

    #compute mass
    Owntools.Compute.Compute_mass(dict_sample)
    dict_tracker['L_mass'].append(dict_sample['grains_mass'])
    dict_tracker['L_mass_dissolved'].append(dict_tracker['L_mass_dissolved'][-1] + (dict_sample['L_mass'][-2]-dict_sample['L_mass'][-1]))

#-------------------------------------------------------------------------------

def DEM_loading(dict_algorithm, dict_geometry, dict_material, dict_sample, dict_sollicitation, simulation_report):
    """
    Loading the granular system to find a steady state.

        Input :
            an algorithm dictionnary (a dict)
            a geometry dictionnary (a dict)
            a material dictionnary (a dict)
            a smaple dictionnary (a dict)
            a sollicitations dictionnary (a dict)
            a simulation report (a report)
        Output :
            Nothing, but sample dictionnary is updated
    """

    i_DEM_0 = dict_algorithm['i_DEM']
    DEM_loop_statut = True

    #trackers and stop conditions
    Force_tracker = []
    Force_stop = 0
    Ecin_tracker = []
    Ecin_stop = 0
    Zmax_tracker = []
    for grain in dict_sample['L_g']:
        Force_stop = Force_stop + 0.5*grain.mass*dict_sollicitation['gravity']
        Ecin_stop = Ecin_stop + 0.5*grain.mass*(dict_algorithm['Ecin_ratio']*grain.radius/dict_algorithm['dt_DEM'])**2

    while DEM_loop_statut :

        dict_algorithm['i_DEM'] = dict_algorithm['i_DEM'] + 1

        #Contact detection
        if (dict_algorithm['i_DEM']-i_DEM_0-1) % dict_algorithm['i_update_neighborhoods']  == 0:
            Contact_gg.Update_Neighborhoods(dict_algorithm, dict_sample)
        Contact_gg.Grains_contact_Neighborhoods(dict_material, dict_sample)

        # Detection of contacts between grain and walls
        if (dict_algorithm['i_DEM']-i_DEM_0-1) % dict_algorithm['i_update_neighborhoods']  == 0:
            wall_neighborhood = Contact_gw.Update_wall_Neighborhoods(dict_sample['L_g'], dict_algorithm['factor_neighborhood'], dict_sample['D_oedo'], dict_sample['z_box_min'], dict_sample['z_box_max'])
        Contact_gw.Grains_Wall_contact_Neighborhood(wall_neighborhood, dict_material, dict_sample)

        #Sollicitation computation
        for grain in dict_sample['L_g']:
             grain.init_F_control(dict_sollicitation['gravity'])
        for contact in dict_sample['L_contact']+dict_sample['L_contact_gw']:
            contact.normal()
            contact.tangential(dict_algorithm['dt_DEM'])

        #Move grains
        for grain in dict_sample['L_g']:
            grain.euler_semi_implicite(dict_algorithm['dt_DEM'], 10*dict_algorithm['Ecin_ratio'])

        #check if some grains are outside of the study box
        L_ig_to_delete = []
        for id_grain in range(len(dict_sample['L_g'])):
            if np.linalg.norm([dict_sample['L_g'][id_grain].center[0], dict_sample['L_g'][id_grain].center[1]]) > dict_sample['D_oedo']/2:
                L_ig_to_delete.append(id_grain)
            elif dict_sample['L_g'][id_grain].center[2] < dict_sample['z_box_min'] :
                L_ig_to_delete.append(id_grain)
            elif dict_sample['L_g'][id_grain].center[2] > dict_sample['z_box_max'] :
                L_ig_to_delete.append(id_grain)
        L_ig_to_delete.reverse()
        for id_grain in L_ig_to_delete:
            simulation_report.write_and_print('Grain '+str(dict_sample['L_g'][id_grain].id)+' has been deleted because it is out of the box\n','Grain '+str(dict_sample['L_g'][id_grain].id)+' has been deleted because it is out of the box')
            dict_sample['L_g'].pop(id_grain)
            #if grain is deleted, be carefull with :
            #   - delete contact / contact_gw
            #   - update L_contact_ij / L_contact_gw_ij as indentation change
            #for the moment it is not done

        #Control the z_max to have the pressure target
        Reset, result, Fv = Owntools.Control_z_max_NR(dict_sample['z_box_max'], dict_sollicitation['Vertical_Confinement_Force'], dict_sample['L_contact_gw'], dict_sample['L_g'])
        if Reset :
            dict_sample['z_box_max'] = result
        else :
            dict_sample['z_box_max'] = dict_sample['z_box_max'] + result

        #Tracker
        F = Owntools.Compute.Compute_F_total(dict_sample['L_g'])
        Ecin = Owntools.Compute.Compute_E_cin_total(dict_sample['L_g'])
        Force_tracker.append(F)
        Ecin_tracker.append(Ecin)
        Zmax_tracker.append(dict_sample['z_box_max'])

        if dict_algorithm['i_DEM'] % dict_algorithm['i_print_plot'] ==0:
            if dict_sollicitation['gravity'] > 0 :
                print('i_DEM',str(dict_algorithm['i_DEM'])+'/'+str(dict_algorithm['i_DEM_stop']+i_DEM_0),'and Ecin',int(100*Ecin/Ecin_stop),'% and Force',int(100*F/Force_stop),'% and Confinement',int(100*Fv/dict_sollicitation['Vertical_Confinement_Force']),'%')
            else :
                print('i_DEM',str(dict_algorithm['i_DEM'])+'/'+str(dict_algorithm['i_DEM_stop']+i_DEM_0),'and Ecin',int(100*Ecin/Ecin_stop),'% and Confinement',int(100*Fv/dict_sollicitation['Vertical_Confinement_Force']),'%')
            if dict_algorithm['Debug_DEM'] :
                Owntools.Write.Write_vtk('Debug/Configuration/Main/configuration_'+str(dict_algorithm['i_DEM'])+'.vtk', dict_sample['L_g'])

        #Check stop conditions for DEM
        if dict_algorithm['i_DEM'] >= dict_algorithm['i_DEM_stop'] + i_DEM_0:
             DEM_loop_statut = False
        if dict_sollicitation['gravity'] > 0:
            if Ecin < Ecin_stop and F < Force_stop and (0.95*dict_sollicitation['Vertical_Confinement_Force']<Fv and Fv<1.05*dict_sollicitation['Vertical_Confinement_Force']):
                  DEM_loop_statut = False
        else:
            if Ecin < Ecin_stop and dict_algorithm['i_DEM'] >= dict_algorithm['i_DEM_stop']*0.1 + i_DEM_0 and (0.95*dict_sollicitation['Vertical_Confinement_Force']<Fv and Fv<1.05*dict_sollicitation['Vertical_Confinement_Force']):
                DEM_loop_statut = False
        if dict_sample['L_g'] == []:
            DEM_loop_statut = False

#-------------------------------------------------------------------------------

def close_main(dict_algorithm, dict_material, dict_sample, dict_sollicitation, dict_tracker, simulation_report):
    '''
    Close the PFDEM.

        Input :
            an algorithm dictionnary (a dict)
            a material dictionnary (a dict)
            a sample dictionnary (a dict)
            a sollicitation dictionnary (a dict)
            a tracker dictionnary (a dict)
            a simulation report (a Report)
        Output :
            Nothing but the dictionnaries and the report are updated
    '''
    #make movie of the different configuration
    if 'Movie' in dict_algorithm['L_flag_plot'] and 'Config' in dict_algorithm['L_flag_plot']:
        Owntools.Plot.Plot_mp4('Debug/Configuration/Configuration_','Debug/Configuration.mp4')
    #make movie of the init-current shape
    if 'Movie' in dict_algorithm['L_flag_plot'] and 'Init_Current_Shape' in dict_algorithm['L_flag_plot']:
        Owntools.Plot.Plot_mp4('Debug/Comparison_Init_Current/Init_Current_Shape_','Debug/Init_Current_Shape.mp4')

    simulation_report.end()

    #final save
    if dict_algorithm['SaveData']:

        Owntools.Save.save_dicts_final(dict_algorithm, dict_material, dict_sample, dict_sollicitation, dict_tracker, simulation_report)
        name_actual_folder = os.path.dirname(os.path.realpath(__file__))
        shutil.copytree(name_actual_folder, '../'+dict_algorithm['foldername']+'/'+dict_algorithm['namefile'])
        os.remove('../'+dict_algorithm['foldername']+'/User_'+dict_algorithm['namefile']+'_tempo.txt')
        os.remove('../'+dict_algorithm['foldername']+'/Report_'+dict_algorithm['namefile']+'_tempo.txt')

#-------------------------------------------------------------------------------

if '__main__' == __name__:

    dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitation, simulation_report = open_simulation()

    create_ic(dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitation, simulation_report)

    #trackers
    Owntools.Compute.Compute_mass(dict_sample)
    Owntools.Compute.Compute_compacity(dict_sample)
    dict_tracker = {
    'L_mass' : [dict_sample['grains_mass']],
    'L_mass_dissolved' : [0],
    'L_z_box_max' : [dict_sample['z_box_max']],
    'L_compacity' : [dict_sample['compacity']],
    'L_k0' : [] # add initial k0
    }

    main_simulation(dict_algorithm, dict_geometry, dict_material, dict_sample, dict_sollicitation, dict_tracker, simulation_report)

    if False :
        close_main(dict_algorithm, dict_material, dict_sample, dict_sollicitation, dict_tracker, simulation_report)

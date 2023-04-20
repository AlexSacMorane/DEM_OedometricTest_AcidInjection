# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This file contains functions used in the simulation.
"""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

import numpy as np
import random
import math
import matplotlib.pyplot as plt

#Own
import Create_IC_IncreaseRadius.Grain_ic
import Create_IC_IncreaseRadius.Contact_gg_ic
import Create_IC_IncreaseRadius.Contact_gw_ic
import Grain
import Owntools.Write
import Owntools.Plot
import Owntools.Save

#-------------------------------------------------------------------------------
#Function
#-------------------------------------------------------------------------------

def LG_tempo(dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitation, simulation_report):
    """
    Create an initial condition

        Input :
            an algorithm dictionnary (a dict)
            a geometry dictionnary (a dict)
            an initial condition dictionnary (a dict)
            a material dictionnary (a dict)
            a sample dictionnary (a dict)
            a sollicitations dictionnary (a dict)
            a simulation report (a report)
        Output :
            Nothing, but dictionnaries are updated
    """
    #Creation of grains
    #grains generation is decomposed in several steps (creation of grain then settlement)
    dict_ic['i_DEM_IC']  = 0
    dict_ic['L_L_g_tempo'] = []
    dict_ic['last_id'] = 0

    #compute dz from e_target and PSD
    dz = dict_sample['D_oedo']*0.6/dict_ic['n_generation']
    dict_sample['z_box_min_ic'] = dict_sample['z_box_min']
    dict_sample['z_box_max'] = dict_sample['z_box_min_ic'] + dz

    #---------------------------------------------------------------------------

    for i_generation in range(1,dict_ic['n_generation']+1) :

        print(f'Generation {i_generation} of grains')

        #add elements in dicts
        dict_ic['L_g_tempo'] = []
        dict_ic['i_generation'] = i_generation

        #create grains
        Create_grains(dict_ic, dict_geometry, dict_sample, dict_material, simulation_report)

        #increase radius
        Increase_radius(dict_ic, dict_material, dict_sample, simulation_report)

        #update element in dict
        dict_sample['z_box_min_ic'] = dict_sample['z_box_max']
        dict_sample['z_box_max'] = dict_sample['z_box_min_ic'] + dz
        dict_ic['L_L_g_tempo'].append(dict_ic['L_g_tempo'].copy())

    simulation_report.tac_tempo('Radius expansion')

    #---------------------------------------------------------------------------

    print('\nCombine generations of grains')

    dict_ic['L_g_tempo'] = []
    for L_g_tempo in dict_ic['L_L_g_tempo']:
        for g_tempo in L_g_tempo:
            dict_ic['L_g_tempo'].append(g_tempo)

    #save
    Owntools.Save.save_dicts_ic('Dicts/save_ic_before_calming_down', dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitation, simulation_report)

    #grains calm down
    print('\nCalm down grains')
    dict_ic['i_generation'] = dict_ic['n_generation']+1
    DEM_loading(dict_ic, dict_geometry, dict_material, dict_sample, dict_sollicitation, False, simulation_report)

    #save
    Owntools.Save.save_dicts_ic('Dicts/save_ic_before_loading', dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitation, simulation_report)

    #apply load
    print('\nApply load')
    dict_ic['i_generation'] = dict_ic['n_generation']+2
    DEM_loading(dict_ic, dict_geometry, dict_material, dict_sample, dict_sollicitation, True, simulation_report)

    simulation_report.write_and_print(str(len(dict_ic['L_g_tempo']))+'/'+str(dict_geometry['N_grain'])+' grains have been created\n','\n'+str(len(dict_ic['L_g_tempo']))+' / '+str(dict_geometry['N_grain'])+' grains have been created')
    simulation_report.write_and_print('H/D = '+str(round((dict_sample['z_box_max']-dict_sample['z_box_min'])/dict_sample['D_oedo'],2))+'\n','H/D = '+str(round((dict_sample['z_box_max']-dict_sample['z_box_min'])/dict_sample['D_oedo'],2)))

#-------------------------------------------------------------------------------

def Create_grains(dict_ic, dict_geometry, dict_sample, dict_material, simulation_report):
    """
    Generate the grains.

    A position is tried. This position must not create overlap with already creared temporary grain. If there is no overlap, a temporary grai nis created.

        Input :
            an initial condition dictionnary (a dict)
            a geometry dictionnary (a dict)
            a sample dictionnary (a dict)
            a material dictionnary (a dict)
            a generation id (a int)
            a simulation report (a report)
        Output :
            Nothing, but initial configuration dictionnary is updated
    """
    #Parameters for the method
    n_created = 0
    for L_g_tempo in dict_ic['L_L_g_tempo']:
        n_created = n_created + len(L_g_tempo)
    number_ratio_1_1and2 = (dict_geometry['mass_ratio_1_1and2']*dict_geometry['R_50_2']**3)/(dict_geometry['mass_ratio_1_1and2']*dict_geometry['R_50_2']**3+(1-dict_geometry['mass_ratio_1_1and2'])*dict_geometry['R_50_1']**3)

    for i in range(n_created, int(dict_geometry['N_grain']*dict_ic['i_generation']/dict_ic['n_generation'])):
        if random.uniform(0,1) < number_ratio_1_1and2 :
            radius = np.random.normal(dict_geometry['R_50_1'], dict_geometry['sigma_psd_1'])
        else :
            radius = np.random.normal(dict_geometry['R_50_2'], dict_geometry['sigma_psd_2'])
        if radius > 0:
            r_to_center = random.uniform(0,dict_sample['D_oedo']/2-1.1*radius)
            angle = random.uniform(0, 2*math.pi)
            center = np.array([r_to_center*math.cos(angle),\
                               r_to_center*math.sin(angle),\
                               random.uniform(dict_sample['z_box_min_ic']+1.1*radius, dict_sample['z_box_max'])-1.1*radius])
            g_tempo = Grain_ic.Grain_Tempo(dict_ic['last_id']+1, center, radius, dict_material)
            dict_ic['L_g_tempo'].append(g_tempo)
            dict_ic['last_id'] = dict_ic['last_id'] + 1
        else :
            simulation_report.write_and_print('Grain '+str(dict_ic['last_id']+1)+' has not been created as r <= 0\n','Grain '+str(dict_ic['last_id']+1)+' has not been created as r <= 0')

#-------------------------------------------------------------------------------

def Increase_radius(dict_ic, dict_material, dict_sample, simulation_report):
    """
    Increase the radius of the particle in several steps.

        Input :
            an initial condition dictionnary (a dict)
            a material dictionnary (a dict)
            a sample dictionnary (a dict)
            a simulation_report (a report)
        Output :
            Nothing, but dictionnaries are updated
    """
    L_Ecin_tracker = []
    L_Ratio_Displacement_MeanRadius_tracker = []
    L_n_contact_tracker = []

    for i in range(1,dict_ic['n_step_increase_radius']+1):
        print('Radius increase '+str(i)+'/'+str(dict_ic['n_step_increase_radius']))

        #increae radius and update properties
        L_radius = []
        for grain in dict_ic['L_g_tempo']:
            grain.v = np.array([0, 0, 0])
            grain.w = np.array([0, 0, 0])
            grain.radius = grain.radius_potential*i/dict_ic['n_step_increase_radius']
            grain.volume = 4/3*math.pi*grain.radius**3
            grain.mass = grain.rho*grain.volume
            grain.inertia = 2/5*grain.mass*grain.radius**2
            L_radius.append(grain.radius)

        #update dt_DEM
        dt_DEM_crit = math.pi*min(L_radius)/(0.16*dict_material['nu']+0.88)*math.sqrt(dict_material['rho']*(2+2*dict_material['nu'])/dict_material['Y']) #s critical time step from O'Sullivan 2011
        dt_DEM = dt_DEM_crit/dict_ic['ratio_dt_DEM_crit_dt_DEM'] #s time step during DEM simulation
        dict_ic['dt_DEM_IC'] = dt_DEM

        #Initialisation
        i_DEM_0 = dict_ic['i_DEM_IC']
        DEM_loop_statut = True
        dict_ic['L_contact'] = []
        dict_ic['L_contact_ij'] = []
        dict_ic['L_contact_gw'] = []
        dict_ic['L_contact_gw_ij'] = []
        dict_ic['id_contact'] = 0

        #trackers
        Ecin_tracker = []
        Ratio_Displacement_MeanRadius_tracker = []
        n_contact_tracker = []

        while DEM_loop_statut :

            dict_ic['i_DEM_IC'] = dict_ic['i_DEM_IC'] + 1

            #trackers
            Ecin_tracker.append(E_cin_total(dict_ic['L_g_tempo'])/len(dict_ic['L_g_tempo']))
            L_Ecin_tracker.append(Ecin_tracker[-1])
            Ratio_Displacement_MeanRadius_tracker.append(Mean_v(dict_ic['L_g_tempo'])*dict_ic['dt_DEM_IC']/np.mean(L_radius))
            L_Ratio_Displacement_MeanRadius_tracker.append(Ratio_Displacement_MeanRadius_tracker[-1])
            n_contact_tracker.append(len(dict_ic['L_contact'])+len(dict_ic['L_contact_gw']))
            L_n_contact_tracker.append(n_contact_tracker[-1])

            #Contact detection
            if (dict_ic['i_DEM_IC']-i_DEM_0-1) % dict_ic['i_update_neighborhoods_ir']  == 0:
                Contact_gg_ic.Update_Neighborhoods(dict_ic, dict_ic['factor_neighborhood_ir'])
            Contact_gg_ic.Grains_contact_Neighborhoods(dict_ic,dict_material)

            # Detection of contacts between grain and walls
            if (dict_ic['i_DEM_IC']-i_DEM_0-1) % dict_ic['i_update_neighborhoods_ir']  == 0:
                wall_neighborhood = Contact_gw_ic.Update_wall_Neighborhoods(dict_ic['L_g_tempo'], dict_ic['factor_neighborhood_ir'], dict_sample['D_oedo'], dict_sample['z_box_min_ic'], dict_sample['z_box_max'])
            Contact_gw_ic.Grains_Wall_contact_Neighborhood(wall_neighborhood, dict_sample['D_oedo'], dict_sample['z_box_min_ic'], dict_sample['z_box_max'], dict_ic, dict_material)

            #Sollicitation computation
            for grain in dict_ic['L_g_tempo']:
                 grain.init_F_control(0)
            for contact in dict_ic['L_contact']+dict_ic['L_contact_gw']:
                contact.normal()
                contact.tangential(dict_ic['dt_DEM_IC'])

            #Move grains
            for grain in dict_ic['L_g_tempo']:
                grain.euler_semi_implicite(dict_ic['dt_DEM_IC'])

            #check if some grains are outside of the study box
            L_ig_to_delete = []
            for id_grain in range(len(dict_ic['L_g_tempo'])):
                if np.linalg.norm([dict_ic['L_g_tempo'][id_grain].center[0], dict_ic['L_g_tempo'][id_grain].center[1]]) > dict_sample['D_oedo']/2:
                    L_ig_to_delete.append(id_grain)
                elif dict_ic['L_g_tempo'][id_grain].center[2] < dict_sample['z_box_min_ic'] :
                    L_ig_to_delete.append(id_grain)
                elif dict_ic['L_g_tempo'][id_grain].center[2] > dict_sample['z_box_max'] :
                    L_ig_to_delete.append(id_grain)
            L_ig_to_delete.reverse()
            for id_grain in L_ig_to_delete:
                simulation_report.write_and_print('Grain '+str(dict_ic['L_g_tempo'][id_grain].id)+' has been deleted because it is out of the box\n','Grain '+str(dict_ic['L_g_tempo'][id_grain].id)+' has been deleted because it is out of the box')
                dict_ic['L_g_tempo'].pop(id_grain)

            if dict_ic['i_DEM_IC'] % dict_ic['i_print_plot_IC'] == 0:
                print('\ti_DEM',str(dict_ic['i_DEM_IC'])+'/'+str(dict_ic['i_DEM_stop_ir']+i_DEM_0),':\n',\
                      '\t\tMean displacement', str(int(1000*Ratio_Displacement_MeanRadius_tracker[-1]))+'‰ of mean radius')
                if max(Ecin_tracker) != 0:
                    print('\t\tKinetic energy',str(int(100*Ecin_tracker[-1]/max(Ecin_tracker)))+'% of max reached')
                if dict_ic['Debug_DEM'] :
                    Owntools.Write.Write_grains_vtk('Debug/Configuration/Init/grains_'+str(dict_ic['i_DEM_IC'])+'.vtk', dict_ic['L_g_tempo'])
                    Owntools.Write.Write_box_vtk('Debug/Configuration/Init/box_'+str(dict_ic['i_DEM_IC'])+'.vtk', dict_sample)
                    Owntools.Plot.Plot_ir_trackers('Debug/Trackers/Init/Radius_Expansion.png', L_Ecin_tracker, L_Ratio_Displacement_MeanRadius_tracker, L_n_contact_tracker)

            #Check stop conditions for DEM
            if dict_ic['i_DEM_IC'] >= dict_ic['i_DEM_stop_ir'] + i_DEM_0:
                 DEM_loop_statut = False
            if dict_ic['L_contact'] == [] and dict_ic['L_contact_gw'] == []:
                DEM_loop_statut = False
            if dict_ic['i_DEM_IC'] > 0.1*dict_ic['i_DEM_stop_ir'] + i_DEM_0:
                if Ecin_tracker[-1] < dict_ic['ratio_Ecin_maxEcin_ir']*max(Ecin_tracker):
                    DEM_loop_statut = False
            if dict_ic['L_g_tempo'] == []:
                DEM_loop_statut = False

#-------------------------------------------------------------------------------

def DEM_loading(dict_ic, dict_geometry, dict_material, dict_sample, dict_sollicitation, control_z_max, simulation_report):
    """
    Loading the granular system.

        Input :
            an initial condition dictionnary (a dict)
            a geometry dictionnary (a dict)
            a material dictionnary (a dict)
            a smaple dictionnary (a dict)
            a sollicitations dictionnary (a dict)
            a simulation report (a report)
        Output :
            Nothing, but initial condition dictionnary is updated
    """
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
    #load data needed
    if dict_ic['i_generation'] >= dict_ic['n_generation']+1:
        z_min = dict_sample['z_box_min']
    else :
        z_min = dict_sample['z_box_min_ic']
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

    i_DEM_0 = dict_ic['i_DEM_IC']
    DEM_loop_statut = True

    #Initialisation
    dict_ic['L_contact'] = []
    dict_ic['L_contact_ij'] = []
    dict_ic['L_contact_gw'] = []
    dict_ic['L_contact_gw_ij'] = []
    dict_ic['id_contact'] = 0

    #update dt_DEM
    L_radius = []
    for grain in dict_ic['L_g_tempo']:
        L_radius.append(grain.radius)
    dt_DEM_crit = math.pi*min(L_radius)/(0.16*dict_material['nu']+0.88)*math.sqrt(dict_material['rho']*(2+2*dict_material['nu'])/dict_material['Y']) #s critical time step from O'Sullivan 2011
    dt_DEM = dt_DEM_crit/dict_ic['ratio_dt_DEM_crit_dt_DEM'] #s time step during DEM simulation
    dict_ic['dt_DEM_IC'] = dt_DEM

    #trackers and stop conditions
    dict_ic['Force_tracker'] = []
    Force_stop = 0
    dict_ic['Ecin_tracker'] = []
    Ecin_stop = 0
    Ratio_Displacement_MeanRadius_tracker = []
    dict_ic['Zmax_tracker'] = []
    dict_ic['s_top_tracker']= []
    k0_tracker = []
    k0_mean_tracker = []
    for grain in dict_ic['L_g_tempo']:
        Force_stop = Force_stop + 0.5*grain.mass*dict_sollicitation['gravity']
        Ecin_stop = Ecin_stop + 0.5*grain.mass*(dict_ic['ratio_meanDisplacement_meanRadius_load']*grain.radius/dict_ic['dt_DEM_IC'])**2

    while DEM_loop_statut :

        dict_ic['i_DEM_IC'] = dict_ic['i_DEM_IC'] + 1

        #Contact detection
        if (dict_ic['i_DEM_IC']-i_DEM_0-1) % dict_ic['i_update_neighborhoods_load']  == 0:
            Contact_gg_ic.Update_Neighborhoods(dict_ic, dict_ic['factor_neighborhood_load'])
        Contact_gg_ic.Grains_contact_Neighborhoods(dict_ic,dict_material)

        # Detection of contacts between grain and walls
        if (dict_ic['i_DEM_IC']-i_DEM_0-1) % dict_ic['i_update_neighborhoods_load']  == 0:
            wall_neighborhood = Contact_gw_ic.Update_wall_Neighborhoods(dict_ic['L_g_tempo'], dict_ic['factor_neighborhood_load'], dict_sample['D_oedo'], z_min, dict_sample['z_box_max'])
        Contact_gw_ic.Grains_Wall_contact_Neighborhood(wall_neighborhood, dict_sample['D_oedo'], z_min, dict_sample['z_box_max'], dict_ic, dict_material)

        #Sollicitation computation
        for grain in dict_ic['L_g_tempo']:
             grain.init_F_control(dict_sollicitation['gravity'])
        for contact in dict_ic['L_contact']+dict_ic['L_contact_gw']:
            contact.normal()
            contact.tangential(dict_ic['dt_DEM_IC'])

        #Move grains
        for grain in dict_ic['L_g_tempo']:
            grain.euler_semi_implicite(dict_ic['dt_DEM_IC'])

        #check if some grains are outside of the study box
        L_ig_to_delete = []
        for id_grain in range(len(dict_ic['L_g_tempo'])):
            if np.linalg.norm([dict_ic['L_g_tempo'][id_grain].center[0], dict_ic['L_g_tempo'][id_grain].center[1]]) > dict_sample['D_oedo']/2:
                L_ig_to_delete.append(id_grain)
            elif dict_ic['L_g_tempo'][id_grain].center[2] < z_min :
                L_ig_to_delete.append(id_grain)
            elif dict_ic['L_g_tempo'][id_grain].center[2] > dict_sample['z_box_max'] :
                L_ig_to_delete.append(id_grain)
        L_ig_to_delete.reverse()
        for id_grain in L_ig_to_delete:
            simulation_report.write_and_print('Grain '+str(dict_ic['L_g_tempo'][id_grain].id)+' has been deleted because it is out of the box\n','Grain '+str(dict_ic['L_g_tempo'][id_grain].id)+' has been deleted because it is out of the box')
            dict_ic['L_g_tempo'].pop(id_grain)

        #Control the z_max to have the pressure target
        #and compute k0
        Force_on_upper_wall = 0
        Force_on_lateral_wall = 0
        for contact in dict_ic['L_contact_gw']:
            if contact.nature == 'gwz_max':
                Force_on_upper_wall = Force_on_upper_wall + contact.Fwg_n
            if contact.nature == 'gwlat':
                Force_on_lateral_wall = Force_on_lateral_wall + contact.Fwg_n
        dz_max = dict_sollicitation['kp_wall']*(Force_on_upper_wall - dict_sollicitation['Vertical_Confinement_Force'])
        if abs(dz_max) > dict_ic['factor_neighborhood_load']*min(L_radius)/dict_ic['i_update_neighborhoods_load']:
            dz_max = np.sign(dz_max)*dict_ic['factor_neighborhood_load']*min(L_radius)/dict_ic['i_update_neighborhoods_load']
        if not control_z_max:
            dz_max = 0
        dict_sample['z_box_max'] = dict_sample['z_box_max'] + dz_max

        #Tracker
        F = F_total(dict_ic['L_g_tempo'])
        Ecin = E_cin_total(dict_ic['L_g_tempo'])
        dict_ic['Force_tracker'].append(F)
        dict_ic['Ecin_tracker'].append(Ecin)
        Ratio_Displacement_MeanRadius_tracker.append(Mean_v(dict_ic['L_g_tempo'])*dict_ic['dt_DEM_IC']/np.mean(L_radius))
        dict_ic['Zmax_tracker'].append(dict_sample['z_box_max'])
        dict_ic['s_top_tracker'].append(Force_on_upper_wall/(math.pi*dict_sample['D_oedo']**2/4))
        k0_tracker.append(dict_sample['D_oedo']/4/(dict_sample['z_box_max']-dict_sample['z_box_min'])*Force_on_lateral_wall/Force_on_upper_wall)
        if len(k0_tracker) >= dict_ic['n_window'] :
            k0_mean_tracker.append(np.mean(k0_tracker[-dict_ic['n_window']:]))

        if dict_ic['i_DEM_IC'] % dict_ic['i_print_plot_IC'] ==0:
            print('\ti_DEM',str(dict_ic['i_DEM_IC'])+'/'+str(dict_ic['i_DEM_stop_load']+i_DEM_0),':\n',\
                  '\t\tMean displacement', str(int(1000*Ratio_Displacement_MeanRadius_tracker[-1]))+'‰ of mean radius\n',\
                  '\t\tConfinement', str(int(100*dict_ic['s_top_tracker'][-1]/dict_sollicitation['Vertical_Confinement_Surface_Force']))+'% of the target value')
            if max(dict_ic['Ecin_tracker']) != 0:
                print('\t\tKinetic energy',str(int(100*dict_ic['Ecin_tracker'][-1]/max(dict_ic['Ecin_tracker'])))+'% of max reached')
            if dict_ic['Debug_DEM'] :
                if not control_z_max :
                    Owntools.Plot.Plot_DEM_trackers('Debug/Trackers/Init/DEM_trackers_init_calming_down.png', dict_ic['Force_tracker'], dict_ic['Ecin_tracker'], Ratio_Displacement_MeanRadius_tracker, dict_ic['Zmax_tracker'], dict_ic['s_top_tracker'], k0_tracker, k0_mean_tracker)
                else :
                    Owntools.Plot.Plot_DEM_trackers('Debug/Trackers/Init/DEM_trackers_init_loading.png', dict_ic['Force_tracker'], dict_ic['Ecin_tracker'], Ratio_Displacement_MeanRadius_tracker, dict_ic['Zmax_tracker'], dict_ic['s_top_tracker'], k0_tracker, k0_mean_tracker)
                Owntools.Write.Write_grains_vtk('Debug/Configuration/Init/grains_'+str(dict_ic['i_DEM_IC'])+'.vtk', dict_ic['L_g_tempo'])
                Owntools.Write.Write_box_vtk('Debug/Configuration/Init/box_'+str(dict_ic['i_DEM_IC'])+'.vtk', dict_sample)

        #Check stop conditions for DEM
        if dict_ic['i_DEM_IC'] >= dict_ic['i_DEM_stop_load'] + i_DEM_0:
             DEM_loop_statut = False
        if dict_sollicitation['gravity'] > 0:
            if Ecin < Ecin_stop and F < Force_stop :
                window_s_top = dict_ic['s_top_tracker'][-dict_ic['n_window']:]
                if (0.95*dict_sollicitation['Vertical_Confinement_Surface_Force']<min(window_s_top) and max(window_s_top)<1.05*dict_sollicitation['Vertical_Confinement_Surface_Force']):
                    DEM_loop_statut = False
        else:
            if Ecin < Ecin_stop and dict_ic['i_DEM_IC'] >= dict_ic['i_DEM_stop_load']*0.1 + i_DEM_0 :
                window_s_top = dict_ic['s_top_tracker'][-dict_ic['n_window']:]
                window_k0_top = k0_tracker[-dict_ic['n_window']:]
                if control_z_max :
                    if (0.95*dict_sollicitation['Vertical_Confinement_Surface_Force']<min(window_s_top) and max(window_s_top)<1.05*dict_sollicitation['Vertical_Confinement_Surface_Force']) and \
                       (max(window_k0_top) - min(window_k0_top) < dict_ic['dk0_window']):
                        DEM_loop_statut = False
                else :
                    DEM_loop_statut = False
        if dict_ic['L_g_tempo'] == []:
            DEM_loop_statut = False


#-------------------------------------------------------------------------------

def E_cin_total(L_g):
    """
    Compute total kinetic energy.

        Input :
            a list of temporary grains (a list)
        Output :
            the total kinetic energy (a float)
    """
    Ecin = 0
    for grain in L_g:
        Ecin = Ecin + 1/2*grain.mass*np.dot(grain.v,grain.v)
    return Ecin

#-------------------------------------------------------------------------------

def Mean_v(L_g):
    """
    Compute the mean particle velocity.

        Input :
            a list of temporary grains (a list)
        Output :
            a mean velocity (a float)
    """
    v_mean = 0
    for grain in L_g:
        v_mean = v_mean + np.linalg.norm(grain.v)
    return v_mean/len(L_g)

#-------------------------------------------------------------------------------

def F_total(L_g):
    """
    Compute total force applied on grains in the sample.

        Input :
            a list of temporary grains (a list)
        Output :
            the total force applied (a float)
    """
    F = 0
    for grain in L_g:
        F = F + np.linalg.norm([grain.fx, grain.fy, grain.fz])
    return F

#-------------------------------------------------------------------------------

def From_LG_tempo_to_usable(dict_ic, dict_sample):
    """
    Create a real grain from a temporary grain.

        Input :
            an initial condition dictionnary (a dict)
            a sample dictionnary (a dict)
        Output :
            Nothing, but the sample dictionnary is updated with the list of real grains
    """
    L_g = []
    for grain_tempo in dict_ic['L_g_tempo']:
        #create real grain
        L_g.append(Grain.Grain(grain_tempo))
    #Add element in dict
    dict_sample['L_g'] = L_g

#-------------------------------------------------------------------------------

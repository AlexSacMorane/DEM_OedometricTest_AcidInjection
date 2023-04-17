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
import Create_IC.Grain_ic
import Create_IC.Contact_gg_ic
import Create_IC.Contact_gw_ic
import Grain
import Owntools.Write
import Owntools.Plot

#-------------------------------------------------------------------------------
#Function
#-------------------------------------------------------------------------------

def LG_tempo(dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitation, simulation_report):
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
    #define the z_max for the grains generation
    dz_creation = dict_geometry['N_grain']/dict_ic['n_generation']*dict_ic['factor_zmax_box']*7*dict_geometry['R_50']**3/dict_sample['D_oedo']**2
    dict_sample['dz_creation'] = dz_creation

    #Creation of grains
    #grains generation is decomposed in several steps (creation of grain then settlement)
    dict_ic['i_DEM_IC']  = 0
    dict_ic['L_L_g_tempo'] = []
    dict_sample['z_box_min_ic'] = dict_sample['z_box_min']
    dict_ic['last_id'] = 0

    #---------------------------------------------------------------------------

    for i_generation in range(1,dict_ic['n_generation']+1) :

        print(f'Generation {i_generation} of grains')

        #add elements in dicts
        dict_ic['L_g_tempo'] = []
        dict_ic['i_generation'] = i_generation

        #create not overlaping grains
        Create_grains(dict_ic, dict_geometry, dict_sample, dict_material, simulation_report)

        #DEM to find the steady-state configuration after loading
        #find the maximum y (center+radius)
        z_max = dict_sample['z_box_min_ic']
        for grain in dict_ic['L_g_tempo']:
            if grain.center[2] + grain.radius > z_max:
                z_max = grain.center[2] + grain.radius
        #add element in dict
        dict_sample['z_box_max'] = z_max

        DEM_loading(dict_ic, dict_geometry, dict_material, dict_sample, dict_sollicitation, simulation_report)

        #update element in dict
        dict_sample['z_box_min_ic'] = dict_sample['z_box_max']

    #---------------------------------------------------------------------------

    print('Combine generations of grains')

    dict_ic['i_generation'] = dict_ic['n_generation']+1

    dict_ic['L_g_tempo'] = []
    for L_g_tempo in dict_ic['L_L_g_tempo']:
        for g_tempo in L_g_tempo:
            dict_ic['L_g_tempo'].append(g_tempo)

    DEM_loading(dict_ic, dict_geometry, dict_material, dict_sample, dict_sollicitation, simulation_report)

    simulation_report.write_and_print(str(len(dict_ic['L_g_tempo']))+'/'+str(dict_geometry['N_grain'])+' grains have been created\n','\n'+str(len(dict_ic['L_g_tempo']))+' / '+str(dict_geometry['N_grain'])+' grains have been created')
    simulation_report.write_and_print('H/D = '+str(round((dict_sample['z_box_max']-dict_sample['z_box_min'])/dict_sample['D_oedo'],2))+'\n','H/D = '+str(round((dict_sample['z_box_max']-dict_sample['z_box_min'])/dict_sample['D_oedo'],2)))

#-------------------------------------------------------------------------------

def DEM_loading(dict_ic, dict_geometry, dict_material, dict_sample, dict_sollicitation, simulation_report):
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
    if dict_ic['i_generation'] == dict_ic['n_generation']+1 :
        i_update_neighborhoods = dict_ic['i_update_neighborhoods_com']
        z_min = dict_sample['z_box_min']
        gravity = dict_ic['gravity']
    else :
        i_update_neighborhoods = dict_ic['i_update_neighborhoods_gen']
        z_min = dict_sample['z_box_min_ic']
        gravity = dict_ic['gravity']
    #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

    i_DEM_0 = dict_ic['i_DEM_IC']
    DEM_loop_statut = True

    #Initialisation
    dict_ic['L_contact'] = []
    dict_ic['L_contact_ij'] = []
    dict_ic['L_contact_gw'] = []
    dict_ic['L_contact_gw_ij'] = []
    dict_ic['id_contact'] = 0

    #trackers and stop conditions
    Force_tracker = []
    Force_stop = 0
    Ecin_tracker = []
    Ecin_stop = 0
    Zmax_tracker = []
    F_top_tracker= []
    s_top_tracker= []
    for grain in dict_ic['L_g_tempo']:
        Force_stop = Force_stop + 0.5*grain.mass*dict_sollicitation['gravity']
        Ecin_stop = Ecin_stop + 0.5*grain.mass*(dict_ic['Ecin_ratio_IC']*grain.radius/dict_ic['dt_DEM_IC'])**2

    while DEM_loop_statut :

        dict_ic['i_DEM_IC'] = dict_ic['i_DEM_IC'] + 1

        #Contact detection
        if (dict_ic['i_DEM_IC']-i_DEM_0-1) % i_update_neighborhoods  == 0:
            Contact_gg_ic.Update_Neighborhoods(dict_ic)
        Contact_gg_ic.Grains_contact_Neighborhoods(dict_ic,dict_material)

        # Detection of contacts between grain and walls
        if (dict_ic['i_DEM_IC']-i_DEM_0-1) % i_update_neighborhoods  == 0:
            wall_neighborhood = Contact_gw_ic.Update_wall_Neighborhoods(dict_ic['L_g_tempo'], dict_ic['factor_neighborhood_IC'], dict_sample['D_oedo'], z_min, dict_sample['z_box_max'])
        Contact_gw_ic.Grains_Wall_contact_Neighborhood(wall_neighborhood, dict_sample['D_oedo'], z_min, dict_sample['z_box_max'], dict_ic, dict_material)

        #Sollicitation computation
        for grain in dict_ic['L_g_tempo']:
             grain.init_F_control(gravity)
        for contact in dict_ic['L_contact']+dict_ic['L_contact_gw']:
            contact.normal()
            contact.tangential(dict_ic['dt_DEM_IC'])

        #Move grains
        for grain in dict_ic['L_g_tempo']:
            grain.euler_semi_implicite(dict_ic['dt_DEM_IC'], 10*dict_ic['Ecin_ratio_IC'])

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
        Reset, result, Fv = Control_z_max_NR(dict_sample['z_box_max'], dict_sollicitation['Vertical_Confinement_Force'], dict_ic['L_contact_gw'], dict_ic['L_g_tempo'])
        if Reset :
            dict_sample['z_box_max'] = result
        else :
            dict_sample['z_box_max'] = dict_sample['z_box_max'] + result

        #Tracker
        F = F_total(dict_ic['L_g_tempo'])
        Ecin = E_cin_total(dict_ic['L_g_tempo'])
        Force_tracker.append(F)
        Ecin_tracker.append(Ecin)
        Zmax_tracker.append(dict_sample['z_box_max'])
        F_top_tracker.append(Fv)
        s_top_tracker.append(Fv/(math.pi*dict_sample['D_oedo']**2/4))

        if dict_ic['i_DEM_IC'] % dict_ic['i_print_plot_IC'] ==0:
            if dict_sollicitation['gravity'] > 0 :
                print('i_DEM',str(dict_ic['i_DEM_IC'])+'/'+str(dict_ic['i_DEM_stop_IC']+i_DEM_0),'and Ecin',int(100*Ecin/Ecin_stop),'% and Force',int(100*F/Force_stop),'% and Confinement',int(100*Fv/dict_sollicitation['Vertical_Confinement_Force']),'%')
            else :
                print('i_DEM',str(dict_ic['i_DEM_IC'])+'/'+str(dict_ic['i_DEM_stop_IC']+i_DEM_0),'and Ecin',int(100*Ecin/Ecin_stop),'% and Confinement',int(100*Fv/dict_sollicitation['Vertical_Confinement_Force']),'%')
            if dict_ic['Debug_DEM'] :
                Owntools.Plot.Plot_DEM_trackers('Debug/Configuration/DEM_trackers_init_'+str(dict_ic['i_generation'])+'.png', Force_tracker, Ecin_tracker, Zmax_tracker, s_top_tracker)
                Owntools.Write.Write_grains_vtk('Debug/Configuration/Init/grains_'+str(dict_ic['i_DEM_IC'])+'.vtk', dict_ic['L_g_tempo'])
                Owntools.Write.Write_box_vtk('Debug/Configuration/Init/box_'+str(dict_ic['i_DEM_IC'])+'.vtk', dict_sample)

        #Check stop conditions for DEM
        if dict_ic['i_DEM_IC'] >= dict_ic['i_DEM_stop_IC'] + i_DEM_0:
             DEM_loop_statut = False
        if dict_sollicitation['gravity'] > 0:
            if Ecin < Ecin_stop and F < Force_stop :
                window_F_top = F_top_tracker[-dict_ic['n_window']:]
                if (0.95*dict_sollicitation['Vertical_Confinement_Force']<min(window_F_top) and max(window_F_top)<1.05*dict_sollicitation['Vertical_Confinement_Force']):
                    DEM_loop_statut = False
        else:
            if Ecin < Ecin_stop and dict_ic['i_DEM_IC'] >= dict_ic['i_DEM_stop_IC']*0.5 + i_DEM_0 :
                window_F_top = F_top_tracker[-dict_ic['n_window']:]
                if (0.95*dict_sollicitation['Vertical_Confinement_Force']<min(window_F_top) and max(window_F_top)<1.05*dict_sollicitation['Vertical_Confinement_Force']):
                    DEM_loop_statut = False
        if dict_ic['L_g_tempo'] == []:
            DEM_loop_statut = False

    #Update dict
    dict_ic['L_L_g_tempo'].append(dict_ic['L_g_tempo'].copy())

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
        #radius = np.random.normal(dict_geometry['R_50'], dict_geometry['sigma_psd'])
        if radius > 0:
            i_test = 0
            grain_created = False
            while (not grain_created) and i_test < dict_ic['N_test_max']:
                i_test = i_test + 1
                r_to_center = random.uniform(0,dict_sample['D_oedo']/2-1.1*radius)
                angle = random.uniform(0, 2*math.pi)
                center = np.array([r_to_center*math.cos(angle),\
                                   r_to_center*math.sin(angle),\
                                   random.uniform(dict_sample['z_box_min_ic']+1.1*radius, dict_sample['z_box_min_ic'] + dict_sample['dz_creation'])])
                g_tempo = Grain_ic.Grain_Tempo(dict_ic['last_id']+1, center, radius, dict_material)
                grain_created = True
                for grain in dict_ic['L_g_tempo']:
                    if Contact_gg_ic.Grains_contact_f(g_tempo, grain):
                        grain_created = False
            if i_test == dict_ic['N_test_max'] and not grain_created:
                simulation_report.write_and_print('Grain '+str(dict_ic['last_id']+1)+' has not been created after '+str(i_test)+' tries\n','Grain '+str(dict_ic['last_id']+1)+' has not been created after '+str(i_test)+' tries')
            else :
                dict_ic['L_g_tempo'].append(g_tempo)
                dict_ic['last_id'] = dict_ic['last_id'] + 1
        else :
            simulation_report.write_and_print('Grain '+str(dict_ic['last_id']+1)+' has not been created as r <= 0\n','Grain '+str(dict_ic['last_id']+1)+' has not been created as r <= 0')

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

def Control_z_max_NR(z_max, Force_target, L_contact_gw, L_g):
    """
    Control the upper wall to apply force.

    A Newton-Raphson method is applied to verify the confinement.
        Input :
            a coordinate of the upper wall (a float)
            a confinement value (a float)
            a list of contact grain - wall (a list)
            a list of temporary grain (a list)
        Output :
            the coordinate of the upper wall (a float)
            a force applied on the upper wall before control (a float)
    """
    F = 0
    overlap_L = []
    k_L = []
    for contact in L_contact_gw:
        if contact.nature == 'gwz_max':
            F = F + contact.Fwg_n
            overlap_L.append(contact.overlap)
            k_L.append(contact.k)
            #compute force applied, save contact overlap and spring

    if overlap_L != []:
        i_NR = 0
        dz = 0
        ite_criteria = True
        #control the upper wall
        if -0.01*Force_target<error_on_zmax_f(dz,overlap_L,k_L,Force_target) and error_on_zmax_f(dz,overlap_L,k_L,Force_target)<0.01*Force_target:
            ite_criteria = False
        while ite_criteria :
            i_NR = i_NR + 1
            dz = dz - error_on_zmax_f(dz,overlap_L,k_L,Force_target)/error_on_zmax_df(dz,overlap_L,k_L)
            if i_NR > 100: #Maximum try
                ite_criteria = False
            if -0.01*Force_target<error_on_zmax_f(dz,overlap_L,k_L,Force_target) and error_on_zmax_f(dz,overlap_L,k_L,Force_target)<0.01*Force_target:
                ite_criteria = False
        #z_max = z_max + dz
        return False, dz, F

    else :
        #if there is no contact with the upper wall, the wall is reset
        z_max = Reset_z_max(L_g, Force_target)
        return True, z_max, F

#-------------------------------------------------------------------------------

def error_on_zmax_f(dz,overlap_L,k_L,Force_target) :
    """
    Compute the function f to control the upper wall. It is the difference between the force applied and the target value.

        Input :
            an increment of the upper wall position (a float)
            a list of overlap for contact between temporary grain and upper wall (a list)
            a list of spring for contact between temporary grain and upper wall (a list)
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
            a list of overlap for contact between temporary grain and upper wall (a list)
            a list of spring for contact between temporary grain and upper wall (a list)
        Output :
            the derivative of error_on_zmax_f() (a float)
    """
    df = 0
    for i in range(len(overlap_L)):
        df = df + 3/2*k_L[i]*(max(overlap_L[i]-dz,0))**(1/2)
    return df

#-------------------------------------------------------------------------------

def Reset_z_max(L_g,Force):
    """
    The upper wall is located as a single contact verify the target value.

        Input :
            the list of temporary grains (a list)
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

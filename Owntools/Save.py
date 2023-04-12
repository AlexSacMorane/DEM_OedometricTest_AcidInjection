# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This file contains the different functions used to save in the simulation.
"""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

import pickle

#-------------------------------------------------------------------------------

def save_dicts_ic(namesave, dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitation, simulation_report):
    '''
    Save dictionnaries.

        Input :
            a file name (a string)
            an algorithm dictionnary (a dictionnary)
            a geometry dictionnary (a dictionnary)
            an initial condition dictionnary (a dictionnary)
            a material dictionnary (a dictionnary)
            a sample dictionnary (a dictionnary)
            a sollicitation dictionnary (a dictionnary)
            a simulation report (a report)
        Output :
            Nothing but a save file is generated (a file)
    '''
    outfile = open(namesave, 'wb')
    dict_save = {}
    dict_save['algorithm'] = dict_algorithm
    dict_save['geometry'] = dict_geometry
    dict_save['ic'] = dict_ic
    dict_save['material'] = dict_material
    dict_save['sample'] = dict_sample
    dict_save['sollicitation'] = dict_sollicitation
    dict_save['report'] = simulation_report
    pickle.dump(dict_save,outfile)
    outfile.close()

#-------------------------------------------------------------------------------

def save_dicts(namesave, dict_algorithm, dict_geometry, dict_material, dict_sample, dict_sollicitation, dict_tracker, simulation_report):
    '''
    Save dictionnaries.

        Input :
            a file name (a string)
            an algorithm dictionnary (a dictionnary)
            a geometry dictionnary (a dictionnary)
            a material dictionnary (a dictionnary)
            a sample dictionnary (a dictionnary)
            a sollicitation dictionnary (a dictionnary)
            a tracker dictionnary (a dictionnary)
            a simulation report (a report)
        Output :
            Nothing but a save file is generated (a file)
    '''
    outfile = open(namesave, 'wb')
    dict_save = {}
    dict_save['algorithm'] = dict_algorithm
    dict_save['geometry'] = dict_geometry
    dict_save['material'] = dict_material
    dict_save['sample'] = dict_sample
    dict_save['sollicitation'] = dict_sollicitation
    dict_save['tracker'] = dict_tracker
    dict_save['report'] = simulation_report
    pickle.dump(dict_save,outfile)
    outfile.close()

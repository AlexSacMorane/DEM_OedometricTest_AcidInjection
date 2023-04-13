# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This is the file to restart a simulation after a crash.
There is a save at the end of each dissolution iteration.
"""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

from pathlib import Path
from datetime import datetime
import numpy as np
import os
import shutil
import math
import pickle

#Own functions and classes
import Grain
import User
import Report
import main
import Owntools.Compute

#-------------------------------------------------------------------------------
#User
#-------------------------------------------------------------------------------

name_to_load = 'Dicts/save_ic'

#-------------------------------------------------------------------------------
#load data
#-------------------------------------------------------------------------------

toload = open(name_to_load,'rb')
dict_save = pickle.load(toload,encoding = 'bytes')
toload.close()
dict_algorithm = dict_save['algorithm']
dict_geometry = dict_save['geometry']
dict_material = dict_save['material']
dict_sample = dict_save['sample']
dict_sollicitation = dict_save['sollicitation']
if name_to_load == 'Dicts/save_tempo':
    dict_tracker = dict_save['tracker']
if name_to_load == 'Dicts/save_ic':
    dict_ic = dict_save['ic']
simulation_report = dict_save['report']

#-------------------------------------------------------------------------------
#Plan the simulation
#-------------------------------------------------------------------------------

simulation_report.write('\nA crash occurs...\n\n')

#-------------------------------------------------------------------------------
#main
#-------------------------------------------------------------------------------

if name_to_load == 'Dicts/save_ic' :

    #trackers
    Owntools.Compute.Compute_mass(dict_sample)
    Owntools.Compute.Compute_compacity(dict_sample)
    Owntools.Compute.Compute_k0_ic(dict_ic, dict_sample)
    dict_tracker = {
    'L_mass' : [dict_sample['grains_mass']],
    'L_mass_dissolved' : [0],
    'L_z_box_max' : [dict_sample['z_box_max']],
    'L_compacity' : [dict_sample['compacity']],
    'L_k0' : [dict_sample['k0']]
    }

main.main_simulation(dict_algorithm, dict_geometry, dict_material, dict_sample, dict_sollicitation, dict_tracker, simulation_report)
close_main(dict_algorithm, dict_geometry, dict_material, dict_sample, dict_sollicitation, dict_tracker, simulation_report)


#-------------------------------------------------------------------------------
#close simulation
#-------------------------------------------------------------------------------

main.close_main(dict_algorithm, dict_material, dict_sample, dict_sollicitation, dict_tracker, simulation_report)

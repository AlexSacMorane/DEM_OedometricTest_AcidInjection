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

#-------------------------------------------------------------------------------
#User
#-------------------------------------------------------------------------------

name_to_load = 'Dicts/save_tempo'

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
dict_tracker = dict_save['tracker']
simulation_report = dict_save['report']

#-------------------------------------------------------------------------------
#Plan the simulation
#-------------------------------------------------------------------------------

simulation_report.write('\nA crash occurs...\n\n')

#-------------------------------------------------------------------------------
#main
#-------------------------------------------------------------------------------

main.main_simulation(dict_algorithm, dict_geometry, dict_material, dict_sample, dict_sollicitation, dict_tracker, simulation_report)

#-------------------------------------------------------------------------------
#close simulation
#-------------------------------------------------------------------------------

main.close_main(dict_algorithm, dict_material, dict_sample, dict_sollicitation, dict_tracker, simulation_report)

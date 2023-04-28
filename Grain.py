# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

The goal of this file is to define a new class.
The new class is about the grains
"""

#-------------------------------------------------------------------------------
#Libs
#-------------------------------------------------------------------------------

import math
import numpy as np

#-------------------------------------------------------------------------------
#Class
#-------------------------------------------------------------------------------

class Grain:

    #---------------------------------------------------------------------------

    def __init__(self, grain_tempo):
        '''
        Defining a sphere grain

            Input :
                a tempo grain (a grain_ic)
            Output :
                a grain (a grain)
        '''
        self.id = grain_tempo.id
        self.center = grain_tempo.center.copy()
        self.theta_x = 0
        self.theta_y = 0
        self.theta_z = 0
        self.radius = grain_tempo.radius
        self.volume = grain_tempo.volume
        self.rho = grain_tempo.rho
        self.mass = grain_tempo.mass
        self.inertia = grain_tempo.inertia
        self.y = grain_tempo.y
        self.nu = grain_tempo.nu
        self.g = grain_tempo.g #shear modulus
        self.fx = 0
        self.fy = 0
        self.fz = 0
        self.v = np.array([0, 0, 0])
        self.mx = 0
        self.my = 0
        self.mz = 0
        self.w = np.array([0, 0, 0])

    #-------------------------------------------------------------------------------

    def euler_semi_implicite(self, dt_DEM, v_limit):
        """
        Move the grain following a semi implicit euler scheme.

            Input :
                itself (a grain)
                a time step (a float)
            Output :
                Nothing, but the grain is moved
        """
        #translation
        a_i = np.array([self.fx, self.fy, self.fz])/self.mass
        self.v = self.v + a_i*dt_DEM
        self.center = self.center + self.v*dt_DEM

        #rotation
        dw_i = np.array([self.mx, self.my, self.mz])/self.inertia
        self.w = self.w + dw_i*dt_DEM
        self.theta_x = self.theta_x + self.w[0]*dt_DEM
        self.theta_y = self.theta_y + self.w[1]*dt_DEM
        self.theta_z = self.theta_z + self.w[2]*dt_DEM

    #-------------------------------------------------------------------------------

    def init_F_control(self, gravity):
        """
        Initialize the force applied to the grain.

        A gravity of g is applied.

            Input :
                itself (a grain)
                a gravity (a float)
            Ouput :
                Nothing, but the force applied on the grain is initialized
        """
        self.fx = 0
        self.fy = 0
        self.fz = -gravity*self.mass
        self.f = np.array([self.fx, self.fy, self.fz])
        self.mx = 0
        self.my = 0
        self.mz = 0
        self.m = np.array([self.mx, self.my, self.mz])

    #-------------------------------------------------------------------------------

    def add_F(self, F, p_application):
        """
        Add a force to the grain.

            Input :
                itself (a grain)
                a force applied (a 1 x 3 numpy array)
                a application point (a 1 x 3 numpy array)
            Output :
                Nothing, but attributes are updated (three floats)
        """
        self.fx = self.fx + F[0]
        self.fy = self.fy + F[1]
        self.fz = self.fz + F[2]
        v1 = np.array(p_application) - np.array(self.center)
        v2 = np.array(F)
        self.mx = self.mx + np.cross(v1,v2)[0]
        self.my = self.my + np.cross(v1,v2)[1]
        self.mz = self.mz + np.cross(v1,v2)[2]

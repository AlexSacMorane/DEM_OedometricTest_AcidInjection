# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This file contains functions used in the simulation.
"""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

import math
import numpy as np

#-------------------------------------------------------------------------------
#Class
#-------------------------------------------------------------------------------

class Grain_Tempo:
  """
  A temporary grain used to generated an initial condition.
  """

#-------------------------------------------------------------------------------

  def __init__(self, ID, Center, Radius, dict_material):
    """Defining the grain.

        Input :
            itself (a grain_tempo)
            an id (a int)
            a center coordinate (a 1 x 2 numpy array)
            a radius (a float)
            a material dictionnary (a dict)
            a grain type, disk or square (a float)
        Output :
            Nothing, but a temporary grain is generated (a grain_tempo)
    """
    self.id = ID
    self.center = np.array(Center)
    self.theta_x = 0
    self.theta_y = 0
    self.theta_z = 0
    self.radius = 0
    self.radius_potential = Radius
    self.volume = 4/3*math.pi*self.radius**3
    self.rho = dict_material['rho']
    self.mass = self.rho*self.volume
    self.inertia = 2/5*self.mass*self.radius**2
    self.y = dict_material['Y']
    self.nu = dict_material['nu']
    self.g = dict_material['Y']/2/(1+dict_material['nu']) #shear modulus
    self.fx = 0
    self.fy = 0
    self.fz = 0
    self.v = np.array([0, 0, 0])
    self.mx = 0
    self.my = 0
    self.mz = 0
    self.w = np.array([0, 0, 0])

#-------------------------------------------------------------------------------

  def add_F(self, F, p_application):
      """
      Add a force to the grain.

        Input :
            itself (a grain_tempo)
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

#-------------------------------------------------------------------------------

  def init_F_control(self,g):
      """
      Initialize the force applied to the grain.

      A gravity is assumed.

        Input :
            itself (a grain_tempo)
            a gravity value (a float)
        Output :
            Nothing, but attributes concerning the force applied are initialized (floats)
      """
      self.fx = 0
      self.fy = 0
      self.fz = -g*self.mass
      self.mx = 0
      self.my = 0
      self.mz = 0

#-------------------------------------------------------------------------------

  def euler_semi_implicite(self, dt_DEM):
    """
    Move the grain following a semi implicit euler scheme.

        Input :
            itself (a grain_tempo)
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

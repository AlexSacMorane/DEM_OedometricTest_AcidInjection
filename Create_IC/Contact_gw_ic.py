# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This file contains functions used in the simulation.
"""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

import numpy  as np
import math

#Own
import Create_IC.Grain_ic

#-------------------------------------------------------------------------------
#Class
#-------------------------------------------------------------------------------

class Contact_gw_Tempo:
  """
  A temporary contact grain - wall used to generated an initial condition.
  """

#-------------------------------------------------------------------------------

  def __init__(self, ID, G, dict_material, Nature, Nwg, Overlap):
    """
    Defining the contact grain-wall.

        Input :
             itself (a contact_gw_tempo)
             an id (a int)
             a grain (a grain_tempo)
             a material dictionnary (a dict)
             the nature of the wall (a string)
             the normal of the contact (a 3x1 numpy array)
             an overlap (a float)
         Output :
             a contact grain - wall is generated
    """
    self.id = ID
    self.g = G
    self.nature = Nature
    self.overlap = Overlap
    self.nwg = Nwg
    factor = 5 #factor just to increase the stiffness
    self.k = factor*4/3*self.g.y/(1-self.g.nu*self.g.nu)*math.sqrt(self.g.radius) #Hertz law
    self.coeff_restitution = dict_material['coeff_restitution']
    self.mu = 0
    self.kt = 0
    self.tangential_old_statut = False
    self.overlap_tangential = np.array([0, 0, 0])

#-------------------------------------------------------------------------------

  def update(self, new_overlap, new_nwg):
    '''
    Update the overlap of a contact already created.

        Input :
            itself (a contact_gw_tempo)
            an overlap (a float)
            a normal (a 3 x 1 numpy array)
        Output :
            Nothing, but the attribut concerning the overlap is updated (a float)
    '''
    self.overlap = new_overlap
    self.nwg = new_nwg

#-------------------------------------------------------------------------------

  def normal(self):
    """
    Compute the normal reaction of a contact grain-wall.

    Here a pontual spring is considered

        Input :
            itself (a contact_gw_tempo)
        Output :
            Nothing, but attributes are updated
    """
    #conditions "if" are defined and same for each wall nature
    if self.nature == 'gwlat':
        #unlinear stiffness
        Fwg_n = self.k*self.overlap**(3/2)
        Fwg = Fwg_n*self.nwg
        self.Fwg_n = Fwg_n
        self.g.add_F(Fwg, self.g.center - self.g.radius*self.nwg)
        #damping
        gamma = -math.log(self.coeff_restitution)/math.sqrt(math.pi**2+math.log(self.coeff_restitution)**2)
        mass_eq = self.g.mass
        eta = 2 * gamma * math.sqrt(mass_eq*self.k)
        Fwg_damp_n = -np.dot(self.g.v,self.nwg)*eta
        Fwg_damp = Fwg_damp_n*self.nwg
        self.Fwg_damp_n = Fwg_damp_n
        self.g.add_F(Fwg_damp, self.g.center - self.g.radius*self.nwg)

    elif self.nature == 'gwz_min':
        #unlinear stiffness
        Fwg_n = self.k*self.overlap**(3/2)
        Fwg = Fwg_n*self.nwg
        self.Fwg_n = Fwg_n
        self.g.add_F(Fwg, self.g.center - self.g.radius*self.nwg)
        #damping
        gamma = -math.log(self.coeff_restitution)/math.sqrt(math.pi**2+math.log(self.coeff_restitution)**2)
        mass_eq = self.g.mass
        eta = 2 * gamma * math.sqrt(mass_eq*self.k)
        Fwg_damp_n = -np.dot(self.g.v,self.nwg)*eta
        Fwg_damp = Fwg_damp_n*self.nwg
        self.Fwg_damp_n = Fwg_damp_n
        self.g.add_F(Fwg_damp, self.g.center - self.g.radius*self.nwg)

    elif self.nature == 'gwz_max':
        #unlinear stiffness
        Fwg_n = self.k*self.overlap**(3/2)
        Fwg = Fwg_n*self.nwg
        self.Fwg_n = Fwg_n
        self.g.add_F(Fwg, self.g.center - self.g.radius*self.nwg)
        #damping
        Fwg_damp_n = 0
        self.Fwg_damp_n = Fwg_damp_n

#-------------------------------------------------------------------------------

  def tangential(self, dt_DEM):
   """
   Compute the tangential reaction of a contact grain-wall.

   Here a pontual spring is considered.

        Input :
            itself (a contact_gw_tempo)
            a time step (a float)
        Output :
            Nothing, but attributes are updated
   """
   if self.overlap > 0 and self.mu > 0:

       if self.tangential_old_statut:
           #if a reaction has been already computed
           #need to project the tangential overlap into the new tangential plane
           self.overlap_tangential = self.overlap_tangential - np.dot(self.overlap_tangential, self.pc_normal)*self.pc_normal
       else:
           self.tangential_old_statut = True

       #compute the tangential overlap
       r = self.g.radius - self.overlap
       v_g = self.g.v + r*np.cross(self.pc_normal, self.g.w)
       v_tangential = v_g - np.dot(v_g, self.pc_normal)*self.pc_normal
       self.overlap_tangential = self.overlap_tangential + v_tangential*dt_DEM

       #compute the tangential force
       F_g_w_t = - self.kt*self.overlap_tangential
       if np.linalg.norm(F_g_w_t) > abs(self.Fwg_n*self.mu) or kt == 0:
           F_g_w_t = abs(self.Fwg_n*self.mu) * F_g_w_t/np.linalg.norm(F_g_w_t)
           self.ft = 0
       else :
           self.ft = np.linalg.norm(F_g_w_t)
       self.g.add_F( F_g_w_t, self.g.center + self.g.radius*self.pc_normal)

   #no contact finally
   else :
       tangential_old_statut = False
       self.overlap_tangential = np.array([0, 0, 0])
       self.ft = 0
       self.ft_damp = 0

#-------------------------------------------------------------------------------
#Function
#-------------------------------------------------------------------------------

def Update_wall_Neighborhoods(L_g_tempo, factor_neighborhood_IC, d_oedo, z_min, z_max):
    """
    Determine a neighborhood for wall.

    This function is called every x time step. The grain - wall contact is determined by Grains_Wall_contact_Neighborhood().
    A factor determines the size of the neighborhood window.

        Input :
            a list of temporary grains (a list)
            a factor to determine the neighborhood window (a float)
            the coordinates of the walls (floats)
        Output :
            a list of temporary grains in the neighborhood of the walls (a list)
    """
    wall_neighborhood = []
    for grain in L_g_tempo:

        r_to_center = np.linalg.norm([grain.center[0], grain.center[1]])
        p_z_min = grain.center[2] - grain.radius
        p_z_max = grain.center[2] + grain.radius

        #grain-wall lateral wall
        if d_oedo/2 - (r_to_center+grain.radius) < factor_neighborhood_IC*grain.radius :
            wall_neighborhood.append(grain)
        #grain-wall z_min
        if p_z_min - z_min < factor_neighborhood_IC*grain.radius :
            wall_neighborhood.append(grain)
        #grain-wall z_max
        if z_max - p_z_max < factor_neighborhood_IC*grain.radius :
            wall_neighborhood.append(grain)

    return wall_neighborhood

#-------------------------------------------------------------------------------

def Grains_Wall_contact_Neighborhood(wall_neighborhood, d_oedo, z_box_min, z_box_max, dict_ic, dict_material):
  """
  Detect contact grain in the neighborhood of the wall and the wall.

  The neighborhood is updated with Update_wall_Neighborhoods(). An iteration over the grains in the wall neighborhood is done. A comparison is done with the coordinates of the wall to determine if there is a contact.

        Input :
            a walls neighborhood (a list)
            the coordinates of the walls (floats)
            an initial condition dictionnary (a dict)
            a material dictionnary (a dict)
        Output :
            Nothing, but the initial condition dictionnary is updated with the contact grain - walls.
  """
  #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
  #load data needed
  L_ij_contact_gw = dict_ic['L_contact_gw_ij']
  L_contact_gw = dict_ic['L_contact_gw']
  id_contact = dict_ic['id_contact']
  #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-

  for grain in wall_neighborhood:

      r_to_center = np.linalg.norm([grain.center[0], grain.center[1]])
      p_z_min = grain.center[2] - grain.radius
      p_z_max = grain.center[2] + grain.radius

      #grain-wall lateral wall
      if d_oedo/2 < (r_to_center+grain.radius) and (grain.id,-1) not in L_ij_contact_gw:
          overlap = (r_to_center+grain.radius) - d_oedo/2
          nwg = np.array([-grain.center[0], -grain.center[1], 0])/np.linalg.norm(np.array([-grain.center[0], -grain.center[1], 0]))
          L_contact_gw.append(Contact_gw_Tempo(id_contact, grain, dict_material, 'gwlat', nwg, overlap))
          L_ij_contact_gw.append((grain.id,-1))
          id_contact = id_contact + 1
      elif d_oedo/2 < (r_to_center+grain.radius) and (grain.id,-1) in L_ij_contact_gw:
          overlap = (r_to_center+grain.radius) - d_oedo/2
          nwg = np.array([-grain.center[0], -grain.center[1], 0])/np.linalg.norm(np.array([-grain.center[0], -grain.center[1], 0]))
          L_contact_gw[L_ij_contact_gw.index((grain.id,-1))].update(overlap, nwg)
      elif d_oedo/2 > (r_to_center+grain.radius) and (grain.id,-1) in L_ij_contact_gw:
          i_contact = L_ij_contact_gw.index((grain.id,-1))
          L_contact_gw.pop(i_contact)
          L_ij_contact_gw.pop(i_contact)
      #grain-wall z_min
      if p_z_min < z_box_min and (grain.id,-2) not in L_ij_contact_gw:
          overlap = z_box_min - p_z_min
          nwg = np.array([0, 0, 1])
          L_contact_gw.append(Contact_gw_Tempo(id_contact, grain, dict_material, 'gwz_min', nwg, overlap))
          L_ij_contact_gw.append((grain.id,-2))
          id_contact = id_contact + 1
      elif p_z_min < z_box_min and (grain.id,-2) in L_ij_contact_gw:
          overlap = z_box_min - p_z_min
          nwg = np.array([0, 0, 1])
          L_contact_gw[L_ij_contact_gw.index((grain.id,-2))].update(overlap, nwg)
      elif p_z_min > z_box_min and (grain.id,-2) in L_ij_contact_gw:
          i_contact = L_ij_contact_gw.index((grain.id,-2))
          L_contact_gw.pop(i_contact)
          L_ij_contact_gw.pop(i_contact)
      #grain-wall z_max
      if z_box_max < p_z_max and (grain.id,-3) not in L_ij_contact_gw:
          overlap = p_z_max - z_box_max
          nwg = np.array([0, 0, -1])
          L_contact_gw.append(Contact_gw_Tempo(id_contact, grain, dict_material, 'gwz_max', nwg, overlap))
          L_ij_contact_gw.append((grain.id,-3))
          id_contact = id_contact + 1
      elif z_box_max < p_z_max and (grain.id,-3) in L_ij_contact_gw:
          overlap = p_z_max - z_box_max
          nwg = np.array([0, 0, -1])
          L_contact_gw[L_ij_contact_gw.index((grain.id,-3))].update(overlap, nwg)
      elif z_box_max > p_z_max and (grain.id,-3) in L_ij_contact_gw:
          i_contact = L_ij_contact_gw.index((grain.id,-3))
          L_contact_gw.pop(i_contact)
          L_ij_contact_gw.pop(i_contact)

      #Update dict
      dict_ic['L_contact_gw_ij'] = L_ij_contact_gw
      dict_ic['L_contact_gw'] = L_contact_gw
      dict_ic['id_contact'] = id_contact

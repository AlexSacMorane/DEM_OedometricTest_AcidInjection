# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

The goal of this file is to define a new class.
The new class is about the contact between two grains.
"""

#-------------------------------------------------------------------------------
#Librairies
#-------------------------------------------------------------------------------

import numpy as np
import math

#Own
import Grain

#-------------------------------------------------------------------------------
#Class
#-------------------------------------------------------------------------------

class Contact_gg:
  """
  A contact grain - grain used.
  """

#-------------------------------------------------------------------------------

  def __init__(self, ID, G1, G2, dict_material):
    """
    Defining the contact grain-grain.

        Input :
            itself (a contact grain - grain)
            an id (a int)
            two grains (two grains)
            a material dictionnary (a dict)
        Output :
            Nothing, but the contact grain - grain is generated
    """
    self.id = ID
    self.g1 = G1
    self.g2 = G2
    Y_eq = 1/((1-self.g1.nu*self.g1.nu)/self.g1.y+(1-self.g2.nu*self.g2.nu)/self.g2.y)
    R_eq = 1/(1/self.g1.radius+1/self.g2.radius)
    k = 4/3*Y_eq*math.sqrt(R_eq)
    self.k = k
    self.coeff_restitution = dict_material['coeff_restitution']
    self.mu = dict_material['mu_friction_gg']
    self.tangential_old_statut = False
    self.overlap_tangential = np.array([0, 0, 0])

#-------------------------------------------------------------------------------

  def init_contact(self,L_g):
    """
    Initialize the contact.

    The grains are updated, the tangetial reaction is set to 0 and a boolean attribute is set to False (new contact grain - grain).

        Input :
            itself (a contact grain - grain)
            a list of grains (a list)
        Output :
            Nothing, but attributes are updated (two grains, two floats and a boolean)
    """
    self.g1 = L_g[self.g1.id]
    self.g2 = L_g[self.g2.id]
    self.tangential_old_statut = False
    self.overlap_tangential = np.array([0, 0, 0])

#-------------------------------------------------------------------------------

  def normal(self):
    """
    Compute the normal reaction of a contact grain-grain.

    Here a pontual spring is considered.

        Input :
            itself (a contact_tempo)
        Output :
            Nothing, but attributes are updated
    """
    # Compute the normal and tangential planes
    PC_normal = (self.g2.center-self.g1.center)/np.linalg.norm(self.g2.center-self.g1.center)
    self.pc_normal = PC_normal #n12

    # Compute the overlap
    overlap = self.g1.radius + self.g2.radius - np.linalg.norm(self.g1.center - self.g2.center)
    self.overlap_normal = overlap

    # Compute the reaction
    if overlap > 0:

        #Spring term
        F_2_1_n = - self.k * overlap**(3/2)  #unlinear spring
        F_2_1 = F_2_1_n * PC_normal
        self.F_2_1_n = F_2_1_n
        self.g1.add_F( F_2_1, self.g1.center + self.g1.radius*self.pc_normal)
        self.g2.add_F(-F_2_1, self.g2.center - self.g2.radius*self.pc_normal)

        #Damping term
        gamma = -math.log(self.coeff_restitution)/math.sqrt(math.pi**2+math.log(self.coeff_restitution)**2)
        mass_eq = self.g1.mass*self.g2.mass/(self.g1.mass+self.g2.mass)
        eta = 2 * gamma * math.sqrt(mass_eq*self.k)
        F_2_1_damp_n = np.dot(self.g2.v - self.g1.v, PC_normal)*eta
        F_2_1_damp = F_2_1_damp_n *PC_normal
        self.F_2_1_damp = F_2_1_damp_n
        self.g1.add_F( F_2_1_damp, self.g1.center + self.g1.radius*self.pc_normal)
        self.g2.add_F(-F_2_1_damp, self.g2.center - self.g2.radius*self.pc_normal)

    #no contact finally
    else :
        self.F_2_1_n = 0
        self.F_2_1_damp = 0

#-------------------------------------------------------------------------------

  def tangential(self, dt_DEM):
    """
    Compute the tangential reaction of a contact grain-grain.

    Here a pontual spring is considered

        Input :
            itself (a contact_tempo)
            a time step (a float)
        Output :
            Nothing, but attributes are updated
    """
    if self.overlap_normal > 0 and self.mu > 0:

        if self.tangential_old_statut:
          #if a reaction has been already computed
          #need to project the tangential overlap into the new tangential plane
          self.overlap_tangential = self.overlap_tangential - np.dot(self.overlap_tangential, self.pc_normal)*self.pc_normal
        else:
          self.tangential_old_statut = True

        #spring law
        G_eq = 1/((1-self.g1.nu)/self.g1.g+(1-self.g2.nu)/self.g2.g)
        R_eq = 1/(1/self.g1.radius+1/self.g2.radius)
        kt0 = 8 * G_eq *math.sqrt(R_eq*abs(self.overlap_normal))
        kt = kt0*math.sqrt(max(1-2/3*kt0*abs(self.overlap_tangential)/self.mu/abs(self.F_2_1_n),0)) #not linear spring

        #damping law
        gamma = -math.log(self.coeff_restitution)/math.sqrt(math.pi**2+math.log(self.coeff_restitution)**2)
        mass_eq = self.g1.mass*self.g2.mass/(self.g1.mass+self.g2.mass)
        eta = 2 * gamma * math.sqrt(mass_eq*kt)/2 #damping tangential term is corrected from the damping normal term

        #compute the tangential overlap
        r1 = self.g1.radius - self.overlap_normal/2
        r2 = self.g2.radius - self.overlap_normal/2
        v_12 = self.g1.v - self.g2.v + r1*np.cross(self.pc_normal, self.g1.w) + r2*np.cross(self.pc_normal, self.g2.w)
        v_tangential = v_12 - np.dot(v_12, self.pc_normal)*self.pc_normal
        self.overlap_tangential = self.overlap_tangential + v_tangential*dt_DEM

        #compute the tangential force
        F_1_2_t = - kt*self.overlap_tangential - eta*v_tangential
        if np.linalg.norm(F_1_2_t) > abs(self.F_2_1_n*self.mu) or kt == 0:
            F_1_2_t = abs(self.F_2_1_n*self.mu) * F_1_2_t/np.linalg.norm(F_1_2_t)
            self.ft = 0
            self.ft_damp = 0
        else :
            self.ft = np.linalg.norm(-kt*self.overlap_tangential)
            self.ft_damp = np.linalg.norm(-eta*v_tangential)
        self.g1.add_F( F_1_2_t, self.g1.center + self.g1.radius*self.pc_normal)
        self.g2.add_F(-F_1_2_t, self.g2.center - self.g2.radius*self.pc_normal)

    #no contact finally
    else :
        tangential_old_statut = False
        self.overlap_tangential = np.array([0, 0, 0])
        self.ft = 0
        self.ft_damp = 0

#-------------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------------

def Grains_contact_f(g1, g2):
  """
  Detect the contact grain-grain.

    Input :
        two grains (two grains)
    Output :
        a Boolean, True if there is contact between the two grains (a Boolean)
  """
  return np.linalg.norm(g1.center-g2.center) < g1.radius+g2.radius

#-------------------------------------------------------------------------------

def Update_Neighborhoods(dict_algorithm, dict_sample):
    """
    Determine a neighborhood for each grain.

    This function is called every x time step. The contact is determined by Grains_contact_Neighborhoods().
    Notice that if there is a potential contact between grain_i and grain_j, grain_i is not in the neighborhood of grain_j.
    Whereas grain_j is in the neighborhood of grain_i. With i_grain < j_grain.

        Input :
            an algorithm dictionnary (a dict)
            a sample dictionnary (a dict)
        Output :
            Nothing, but the neighborhood of the temporary grains is updated
    """
    for i_grain in range(len(dict_sample['L_g'])-1) :
        neighborhood = []
        for j_grain in range(i_grain+1,len(dict_sample['L_g'])):
            if np.linalg.norm(dict_sample['L_g'][i_grain].center-dict_sample['L_g'][j_grain].center) < dict_algorithm['factor_neighborhood_IC']*(dict_sample['L_g'][i_grain].radius+dict_sample['L_g'][j_grain].radius):
                neighborhood.append(dict_sample['L_g'][j_grain])
        dict_sample['L_g'][i_grain].neighborhood = neighborhood


#-------------------------------------------------------------------------------

def Grains_contact_Neighborhoods(dict_material, dict_sample):
    """
    Detect contact between a grain and grains from its neighborhood.

    The neighborhood is updated with Update_Neighborhoods().

        Input :
            a material dictionnary (a dict)
            a sample dictionnary (a dict)
        Output :
            Nothing, but the sample dictionnary is updated with grain - grain contacts
    """
    for i_grain in range(len(dict_sample['L_g'])-1) :
        grain_i = dict_sample['L_g'][i_grain]
        for neighbor in dict_sample['L_g'][i_grain].neighborhood:
            j_grain = neighbor.id
            grain_j = neighbor
            if Grains_contact_f(grain_i, grain_j):
                if (i_grain,j_grain) not in dict_sample['L_contact_ij']:  #contact not detected previously
                   #creation of contact
                   dict_sample['L_contact_ij'].append((i_grain,j_grain))
                   dict_sample['L_contact'].append(Contact(dict_sample['id_contact'], grain_i, grain_j, dict_material))
                   dict_sample['id_contact'] = dict_sample['id_contact'] + 1

            else :
                if (i_grain,j_grain) in dict_sample['L_contact_ij'] : #contact detected previously is not anymore
                       dict_sample['L_contact'].pop(dict_sample['L_contact_ij'].index((i_grain,j_grain)))
                       dict_sample['L_contact_ij'].remove((i_grain,j_grain))

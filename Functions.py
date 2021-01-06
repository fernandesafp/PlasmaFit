# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 18:11:27 2019

Integral calculations and the corresponding intensities

@author: André Fernandes - afp.fernandes@campus.fct.unl.pt
"""

import numpy as np
import math

kt_c = 1e3
kt_h = 2e4

m_e = 9.1093835611e-31
mec2 = 5.10998946131e5
e_j = 1.602176620898e-19
c_Mw = 2.67618617422916e16
c_nMw = 8.70366940390332e28

x_i = np.asarray([0.0705399, 0.372127, 0.916582, 1.70731, 2.7492, 4.04893,
                  5.61517, 7.45902, 9.59439, 12.0388, 14.8143, 17.9489,
                  21.4788, 25.4517, 29.9326, 35.0134, 40.8331, 47.62, 55.8108,
                  66.5244])
w_i = np.asarray([0.168747, 0.291254, 0.266686, 0.166002, 0.0748261, 0.0249644,
                  0.00620255, 0.00114496, 0.000155742, 1.54014e-05,
                  1.08649e-06, 5.33012e-08, 1.75798e-09, 3.7255e-11,
                  4.76753e-13, 3.37284e-15, 1.15501e-17, 1.53952e-20,
                  5.28644e-24, 1.65646e-28])

#Code that calculates the integral of Maxwellian and non-Maxwellian distributions of the plasma with hot and cold temperatures
#------------------------------------------------------------------------------
def integral(hw_min, cross_section_c, cross_section_h, fraction):

  E_i = x_i*kt_c + hw_min
  f_Mw = np.sum(w_i*E_i*np.sqrt(e_j)*(np.sqrt(E_i+2*mec2)/(E_i+mec2))*
                cross_section_c)*c_Mw*np.exp(-hw_min/kt_c)

  E_i = x_i*kt_h + hw_min
  f_nMw = np.sum(w_i*(E_i*e_j)**1.5*(1+E_i/(2*mec2))**1.5*cross_section_h
                 )*np.sqrt(2/m_e)*c_nMw*kt_h*e_j*np.exp(-hw_min/kt_h)

  n_rate = fraction*f_Mw + (1-fraction)*f_nMw

  return (n_rate)
#------------------------------------------------------------------------------


#Code that calculates the statistical probability (2J+1 over the sum of all 2J+1 states)
#------------------------------------------------------------------------------
def nCr(n,r):
    f = math.factorial
    return f(n) // f(r) // f(n-r)

def orb(orbital):
  order = ['s', 'p', 'd', 'f']

  l = orbital[0]
  n = (4*order.index(l) + 2)                                                   #Number of electrons that fit in the orbital (2n+1)*2 (spin up/down)
  r = int(orbital[1])
  j = nCr(n, r)                                                                #Combination of electronic configurations possible in that state
  return j

def stat_prob(lsj_i, state):
  if len(lsj_i) == 5:
    j_i = 2*int(lsj_i[-2]) + 1
  else:
    fraq = int(lsj_i[3]) / int(lsj_i[5])
    j_i = 2*fraq + 1

  n = int((len(state)-1)/4)                                                    #Know how many orbitals there are based on the length of s (following the standard with a space in between)
  orbital = []
  for i in range(n + 1):
    orbital.append(state[1 + 4*i] + state[2 + 4*i])
  sum_j = 1
  for i in range(len(orbital)):
    sum_j *= orb(orbital[i])
  return j_i/sum_j
#------------------------------------------------------------------------------


#Code can only check up to 10 electron configuration for ground states. It provides a true or false whether or not it is possible for a certain state to come from a ground state through a certain process
#------------------------------------------------------------------------------
ground_states_list = ['1s1','1s2','1s2 2s1','1s2 2s2','1s2 2s2 2p1',
                      '1s2 2s2 2p2','1s2 2s2 2p3','1s2 2s2 2p4','1s2 2s2 2p5',
                      '1s2 2s2 2p6']
def ground_state(state_exc, n):
  '''Check if ground state phase is possible with n ionizations (0 is excitation)

  All processes add one electron to the K-shell, if the first orbital gets
  added an electron too much, then it's impossible to have it come from a
  ground state. Also the first orbital must be 1s1 at least. States such as
  1s2 2p1 are never considered i.e. adding a 2s orbital because it also is
  impossible to come from a ground state
  '''

  s = list(state_exc)

  s[2] = str(int(s[2])+1)                                                      #Adds an electron to the 1s orbital
  if int(s[2]) > 2 or s[1] != 's': return False                                #If we have something like 1s3 or not start with s orbital then it is impossible

  if n == 0:                                                                   #Through excitation K
    s[-1] = str(int(s[-1])-1)                                                  #Removes an electron in the last position
    if s[-1] == '0':                                                           #If there is no electrons, it removes the orbital
      s = s[:-4]
    s = ''.join(s)
    return s in ground_states_list                                             #It does not try to remove electron in the middle, because that immediately makes it impossible to be ground state

  if n == 1:                                                                   #Through ionization K
    s = ''.join(s)
    return s in ground_states_list                                             #Was already added to K-shell

  if n == 2:                                                                   #Through double ionization KL
    s[6] = str(int(s[6])+1)                                                    #Considers adding in the orbital next to the first one
    s = ''.join(s)
    if s in ground_states_list: return True
    elif len(s) == 11:                                                         #If there is a third orbital it adds one there instead
      s = list(s)
      s[6] = str(int(s[6])-1)
      s[10] = str(int(s[10])+1)                                                #If there is more, it adds one in the last position
      s = ''.join(s)
      if s in ground_states_list: return True
    else:                                                                      #There is not a third orbital, then one is added
      s = list(s)
      s[6] = str(int(s[6])-1)                                                  #Remove it to transfer for new orbital
      s.append(' ')
      s.append('2')
      s.append('p')
      s.append('1')                                                            #Adds orbital 2p with one electron
      s = ''.join(s)
      return s in ground_states_list

  if n == 3:                                                                   #Through triple ionization KLL
    s[6] = str(int(s[6])+2)                                                    #Adds two electrons in the following orbital
    s = ''.join(s)
    if s in ground_states_list: return True
    elif len(s) == 11:                                                         #If there's a third orbital transfer one from the 6th index to the 10th
      s = list(s)
      s[6] = str(int(s[6])-1)
      s[10] = str(int(s[10])+1)                                                #Transfers one electron in the middle to the last
      s = ''.join(s)
      if s in ground_states_list: return True
      else:                                                                    #Tries once more if not true
        s = list(s)
        s[6] = str(int(s[6])-1)
        s[10] = str(int(s[10])+1)                                              #Transfers once more to the last orbital
        s = ''.join(s)
        return s in ground_states_list
    else:                                                                      #If no third orbital, add 2p
      s = list(s)
      s[6] = str(int(s[6])-1)
      s.append(' ')
      s.append('2')
      s.append('p')
      s.append('1')
      s = ''.join(s)
      if s in ground_states_list: return True
      else:                                                                    #Transfer one more electron to 2p
        s = list(s)
        s[6] = str(int(s[6])-1)
        s[10] = str(int(s[10])+1)
        s = ''.join(s)
        if s in ground_states_list: return True
        else:                                                                  #Retry and think that maybe we have 1s2 2px and add 1s2 2s2 2px
          s = list(state_exc)                                                  #Reset
          s.insert(3, '2')
          s.insert(3, 's')
          s.insert(3, '2')                                                     #Add orbital 2s in the middle with two electrons
          s.insert(3, ' ')
          s = ''.join(s)
          return s in ground_states_list

  return False                                                                 #If nothing is ever returned as true, then false
#------------------------------------------------------------------------------


#Code that provides the voigt profile (fraction = 1 Gaussian, = 0 Lorentzian)
#------------------------------------------------------------------------------
def voigt(A, x, u, width_gauss, width_lorentz, fraction):
  sigma_gauss = width_gauss/(2*np.sqrt(2*np.log(2)))
  gauss = np.exp(-.5*((x-u)/sigma_gauss)**2)/(sigma_gauss*np.sqrt(2*np.pi))

  lorentz = .5*width_lorentz/(np.pi*((x-u)**2+(.5*width_lorentz)**2))

  return A*(fraction*gauss+(1-fraction)*lorentz)
#------------------------------------------------------------------------------
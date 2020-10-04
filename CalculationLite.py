# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 02:20:25 2019

Calculates Spectrum charge-state density for LM fit

@author: André Fernandes - afp.fernandes@campus.fct.unl.pt
"""

import numpy as np

from Functions import integral, stat_prob, ground_state, voigt
from scipy.interpolate import interp1d                                         #To do the transfer function fit


#Code that provides the final plot
#------------------------------------------------------------------------------
def spectrum_plot_lite(data, counts_exp_err, x, y0, hw_trans, cs_K_exc,
                       cs_K_ion, y, cs_KL_ion, cs_KLL_ion, csd, trans_file,
                       fraction_voigt, width_gauss, width_loren, y0_method,
                       y0_m, y0_b, hw0, fraction_mw, kexc, kion, klion, kllion,
                       metastates):

  if trans_file != '':                                                         #If a transfer function was loaded
    def tf(x):
      data_tf = np.genfromtxt(trans_file, delimiter = ',')
      xtf = data_tf[:,0]
      ytf = data_tf[:,1]
      inter = interp1d(xtf, ytf/max(ytf), kind = 'cubic')                      #Transfer function fit normalized
      return inter(x)
  else:                                                                        #Else, no change in the TF
    def tf(x):
      return np.ones(len(x))

  #Search cross sectional database and calculation
  #----------------------------------------------------------------------------
  #Empty profiles
  voigt_K_exc_profile = np.zeros(len(x))
  voigt_K_ion_profile = np.zeros(len(x))
  voigt_KL_ion_profile = np.zeros(len(x))
  voigt_KLL_ion_profile = np.zeros(len(x))

  counts_K_exc = []
  hw_K_exc = []
  counts_K_ion = []
  hw_K_ion = []
  counts_KL_ion = []
  hw_KL_ion = []
  counts_KLL_ion = []
  hw_KLL_ion = []

  for i in range(len(hw_trans)):                                               #In all the energy transition it will try to attribute cs electronic combinations that will result in the initial position of the decay and then present the radiative decay yield
    #K-SHELL EXCITATION
    #--------------------------------------------------------------------------
    if kexc:
      if ground_state(hw_trans[i,5],0) == 1:                                   #If it is possible to have the energy state to be reached from a ground state position, then it will proceed, else, it will ask if it is a metastate
        for j in range(len(cs_K_exc)):
          if hw_trans[i,5] == cs_K_exc[j,2] and hw_trans[i,7] == cs_K_exc[j,4]:
            n_rate = integral(cs_K_exc[j,0], cs_K_exc[j,6], cs_K_exc[j,7],
                              fraction_mw)
            hw = hw_trans[i,0]
            ty = hw_trans[i,4]
            try: nq = csd.get(str(cs_K_exc[j,5])) * 10**csd.get('power')       #Ion might not exist in dictionary, so it will only try
            except: nq = 0

            counts_K_exc.append(n_rate * hw * ty * nq)                         #Append counts = n_rate * energy * transition yield * charge-state density
            hw_K_exc.append(i)                                                 #Append energy

      elif metastates:                                                         #Else wonder if there's a metastate
        for j in range(len(cs_K_exc)):
          if (cs_K_exc[j,-1] == hw_trans[i,0] and
              hw_trans[i,5] == cs_K_exc[j,2] and
              hw_trans[i,7] == cs_K_exc[j,4]):
            n_rate = integral(cs_K_exc[j,0], cs_K_exc[j,6], cs_K_exc[j,7],
                              fraction_mw)
            hw = hw_trans[i,0]
            ty = hw_trans[i,4]
            try: nq = csd.get(str(cs_K_exc[j,5])+'m') * 10**csd.get('power')
            except: nq = 0

            counts_K_exc.append(n_rate * hw * ty * nq)
            hw_K_exc.append(i)
    #--------------------------------------------------------------------------


    #K-SHELL IONIZATION
    #--------------------------------------------------------------------------
    if kion:                                                                   #Same rationale will be applied for the following processes
      if ground_state(hw_trans[i,5],1) == 1:
        for j in range(len(cs_K_ion)):
          if hw_trans[i,5] == cs_K_ion[j,2]:
            n_rate = integral(cs_K_ion[j,0], cs_K_ion[j,4], cs_K_ion[j,5],
                              fraction_mw)
            hw = hw_trans[i,0]
            ty = hw_trans[i,4]
            try: nq = csd.get(str(cs_K_ion[j,3] - 1)) * 10**csd.get('power')
            except: nq = 0
            sp = stat_prob(hw_trans[i,7], hw_trans[i,5])

            counts_K_ion.append(n_rate * hw * ty * nq * sp)
            hw_K_ion.append(i)

      elif metastates:                                                         #Else wonder if there's a metastate
        for j in range(len(cs_K_ion)):
          if (cs_K_ion[j,-1] == hw_trans[i,0] and
              hw_trans[i,5] == cs_K_ion[j,2]):
            n_rate = integral(cs_K_ion[j,0], cs_K_ion[j,4], cs_K_ion[j,5],
                              fraction_mw)
            hw = hw_trans[i,0]
            ty = hw_trans[i,4]
            try: nq = (csd.get(str(cs_K_ion[j,3] - 1)+'m') * 10**
                       csd.get('power'))
            except: nq = 0
            sp = stat_prob(hw_trans[i,7], hw_trans[i,5])

            counts_K_ion.append(n_rate * hw * ty * nq * sp)
            hw_K_ion.append(i)
    #--------------------------------------------------------------------------


    #KL SHELLS DOUBLE IONIZATION
    #--------------------------------------------------------------------------
    if klion:
      if ground_state(hw_trans[i,5],2) == 1:
        for j in range(len(cs_KL_ion)):
          if hw_trans[i,5] == cs_KL_ion[j,2]:
            n_rate = integral(cs_KL_ion[j,0], cs_KL_ion[j,4], cs_KL_ion[j,5],
                              fraction_mw)
            hw = hw_trans[i,0]
            ty = hw_trans[i,4]
            try: nq = csd.get(str(cs_KL_ion[j,3] - 2)) * 10**csd.get('power')
            except: nq = 0
            sp = stat_prob(hw_trans[i,7], hw_trans[i,5])

            counts_KL_ion.append(n_rate * hw * ty * nq * sp)
            hw_KL_ion.append(i)

      elif metastates:                                                         #Else wonder if there's a metastate
        for j in range(len(cs_KL_ion)):
          if (cs_KL_ion[j,-1] == hw_trans[i,0] and
              hw_trans[i,5] == cs_KL_ion[j,2]):
            n_rate = integral(cs_KL_ion[j,0], cs_KL_ion[j,4], cs_KL_ion[j,5],
                              fraction_mw)
            hw = hw_trans[i,0]
            ty = hw_trans[i,4]
            try: nq = (csd.get(str(cs_KL_ion[j,3] - 2)+'m') * 10**
                       csd.get('power'))
            except: nq = 0
            sp = stat_prob(hw_trans[i,7], hw_trans[i,5])

            counts_KL_ion.append(n_rate * hw * ty * nq * sp)
            hw_KL_ion.append(i)
    #--------------------------------------------------------------------------


    #KLL SHELLS TRIPLE IONIZATION
    #--------------------------------------------------------------------------
    if kllion:
      if ground_state(hw_trans[i,5],3) == 1:
        for j in range(len(cs_KLL_ion)):
          if hw_trans[i,5] == cs_KLL_ion[j,2]:
            n_rate = integral(cs_KLL_ion[j,0], cs_KLL_ion[j,4],
                              cs_KLL_ion[j,5], fraction_mw)
            hw = hw_trans[i,0]
            ty = hw_trans[i,4]
            try: nq = csd.get(str(cs_KLL_ion[j,3] - 3)) * 10**csd.get('power')
            except: nq = 0
            sp = stat_prob(hw_trans[i,7], hw_trans[i,5])

            counts_KLL_ion.append(n_rate * hw * ty * nq * sp)
            hw_KLL_ion.append(i)

      elif metastates:                                                         #Else wonder if there's a metastate
        for j in range(len(cs_KLL_ion)):
          if (cs_KLL_ion[j,-1] == hw_trans[i,0] and
              hw_trans[i,5] == cs_KLL_ion[j,2]):
            n_rate = integral(cs_KLL_ion[j,0], cs_KLL_ion[j,4],
                              cs_KLL_ion[j,5], fraction_mw)
            hw = hw_trans[i,0]
            ty = hw_trans[i,4]
            try: nq = (csd.get(str(cs_KLL_ion[j,3] - 3)+'m') * 10**
                       csd.get('power'))
            except: nq = 0
            sp = stat_prob(hw_trans[i,7], hw_trans[i,5])

            counts_KLL_ion.append(n_rate * hw * ty * nq * sp)
            hw_KLL_ion.append(i)
    #--------------------------------------------------------------------------

  #Simulation
  #----------------------------------------------------------------------------
  if kexc:
    for i in range(len(hw_K_exc)):                                             #Profile making for K shell excitation
      voigt_K_exc_profile += voigt(counts_K_exc[i], x,
                                   hw_trans[hw_K_exc[i], 0] + hw0,
                                   width_gauss, width_loren, fraction_voigt)
  if kion:
    for i in range(len(hw_K_ion)):                                             #Profile making for K shell ionization
      voigt_K_ion_profile += voigt(counts_K_ion[i], x,
                                   hw_trans[hw_K_ion[i], 0] + hw0,
                                   width_gauss, width_loren, fraction_voigt)
  if klion:
    for i in range(len(hw_KL_ion)):                                            #Profile making for KL shells ionization
      voigt_KL_ion_profile += voigt(counts_KL_ion[i], x,
                                    hw_trans[hw_KL_ion[i], 0] + hw0,
                                    width_gauss, width_loren, fraction_voigt)
  if kllion:
    for i in range(len(hw_KLL_ion)):                                           #Profile making for KLL shells ionization
      voigt_KLL_ion_profile += voigt(counts_KLL_ion[i], x,
                                     hw_trans[hw_KLL_ion[i], 0] + hw0,
                                     width_gauss, width_loren,
                                     fraction_voigt)

  voigt_profile = (voigt_K_exc_profile + voigt_K_ion_profile
                   + voigt_KL_ion_profile + voigt_KLL_ion_profile)
  #----------------------------------------------------------------------------

  #Apply transfer function and normalization
  #----------------------------------------------------------------------------
  #TF, profile, y0
  voigt_profile *= tf(x)                                                       #Profile follows the intensity dist.
  #----------------------------------------------------------------------------

  #Chi2 calculation and plotting
  #----------------------------------------------------------------------------
  if y0_method == 'Chi':
    for i in range(len(data[0])):                                              #Restrict exp data to domain
      if data[0][i] >= min(x):
        exp_i = i
        break
    for j in range(len(data[0])):
      if data[0][-j-1] <= max(x):
        exp_f = len(data[0]) - j
        break
    data[0] = data[0][exp_i:exp_f]
    data[1] = data[1][exp_i:exp_f]
    counts_exp_err = counts_exp_err[exp_i:exp_f]

    chi2_max = np.int(2e8)                                                     #High number. Will change to inf later.
    for j in range(len(y)):                                                    #Minimizes chi2 based on different offsets y0
      chi2 = np.int(0)
      spectrum_csd = voigt_profile + y[j]
      for i in range(len(data[0])):
        dif = abs(x - data[0][i])
        n = np.argmin(dif)
        chi2 += (spectrum_csd[n] - data[1][i])**2/(counts_exp_err[i]**2
                *(len(data[0]) - 1))
      if chi2 < chi2_max:
        chi2_max = chi2
        y0 = y[j]
      else: break                                                              #Makes computation quicker but can cause local minimi
    spectrum_csd = voigt_profile + y0                                          #Profile w/ noise, y0 can either come from the one calculated
  else:
    spectrum_csd = voigt_profile + y0                                          #Profile w/ noise y0 set by the user, chi2 calculation proceeds
  #----------------------------------------------------------------------------

  return spectrum_csd
#------------------------------------------------------------------------------
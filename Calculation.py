# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 01:51:21 2019

Calculates and plots spectrum with CSD guess

@author: André Fernandes - afp.fernandes@campus.fct.unl.pt
"""

import numpy as np
from tkinter import messagebox

from Functions import integral, stat_prob, ground_state, voigt
from scipy.interpolate import interp1d                                         #To do the transfer function fit

#Code that provides the final plot
#------------------------------------------------------------------------------
def spectrum_plot(hw_trans, cs_K_exc, cs_K_ion, cs_KL_ion, cs_KLL_ion,
                  spectrum, fig, csd, csv_file, trans_file, stick_sim, scale,
                  param_changed, csd_changed, tf_changed, fraction_voigt,
                  width_gauss, width_loren, step, y0_method, y0_m, y0_b,
                  norm_method, norm_x, norm_y, hw0, hw_min, hw_max,
                  fraction_mw, kexc, kion, klion, kllion, contributions,
                  metastates):

  spectrum.clear()                                                             #New plot
  spectrum.set_xlabel('Energy (eV)')
  spectrum.set_ylabel('Counts/s')

  if csv_file != '': #and stick_sim == 'Sim.':                                 #The experimental data will not be loaded on stickplot mode (only sim.)
    exp_data = True
    data_exp = np.genfromtxt(csv_file, delimiter = ',')
    hw_exp = data_exp[:,0]
    counts_exp = data_exp[:,1]
    spectrum.plot(hw_exp, counts_exp, 'ko', linewidth = .75, markersize = 3,
                  label = 'Experimental data', mfc='none')

    hw_exp_middle = (max(hw_exp) + min(hw_exp))/2                              #Plotting domain
    hw_exp_size = max(hw_exp) - min(hw_exp)
    x_min = hw_exp_middle - hw_exp_size*1.1/2
    x_max = hw_exp_middle + hw_exp_size*1.1/2

    try:                                                                       #Attempt to find errorbars as they're not mandatory
      hw_exp_err = data_exp[:,2]
      counts_exp_err = data_exp[:,3]
      spectrum.errorbar(hw_exp, counts_exp, yerr = counts_exp_err,
                        xerr = hw_exp_err, fmt = 'none', ecolor='k',
                        elinewidth = 0.5, capsize = 1, zorder=1)
    except Exception as ex:
      print('Could not load count error bars. Will be shown as zero. ' + ex)
      counts_exp_err = np.ones(len(hw_exp))                                    #For chi2 display without errors
  else: exp_data = False

  x = np.arange(hw_min, hw_max, step)                                          #Domain arranged
  y0 = np.array(y0_m * x + y0_b)                                               #Offset as an array to simply add to profile

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

  #Domain search in the database
  try:
    for i in range(len(hw_trans)):
      if hw_trans[i,0] >= hw_min:
        trans_i = i
        break
    for j in range(len(hw_trans)):
      if hw_trans[-j-1,0] <= hw_max:
        trans_f = len(hw_trans)-j
        break
    hw_trans = hw_trans[trans_i:trans_f]
  except Exception as ex:                                                      #If it cannot load the database it will warn the user
    print('Could not load the database. Exception: ' + ex)
    messagebox.showerror('Transition energies database', 'The transition ener'+
                         'gies database was incorrectly loaded. The calculati'+
                         'on will stop.')
    return

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
            try:
              nq = csd.get(str(cs_K_exc[j,5])+'m') * 10**csd.get('power')
            except Exception as ex:
              print('Could not obtain CSD. ' + ex)
              nq = 0

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

  if stick_sim == 'Sim.':
    #Simulation
    #--------------------------------------------------------------------------
    if kexc:
      for i in range(len(hw_K_exc)):                                           #Profile making for K shell excitation
        voigt_K_exc_profile += voigt(counts_K_exc[i], x,
                                     hw_trans[hw_K_exc[i], 0] + hw0,
                                     width_gauss, width_loren, fraction_voigt)
    if kion:
      for i in range(len(hw_K_ion)):                                           #Profile making for K shell ionization
        voigt_K_ion_profile += voigt(counts_K_ion[i], x,
                                     hw_trans[hw_K_ion[i], 0] + hw0,
                                     width_gauss, width_loren, fraction_voigt)
    if klion:
      for i in range(len(hw_KL_ion)):                                          #Profile making for KL shells ionization
        voigt_KL_ion_profile += voigt(counts_KL_ion[i], x,
                                      hw_trans[hw_KL_ion[i], 0] + hw0,
                                      width_gauss, width_loren, fraction_voigt)
    if kllion:
      for i in range(len(hw_KLL_ion)):                                         #Profile making for KLL shells ionization
        voigt_KLL_ion_profile += voigt(counts_KLL_ion[i], x,
                                       hw_trans[hw_KLL_ion[i], 0] + hw0,
                                       width_gauss, width_loren,
                                       fraction_voigt)

    voigt_profile = (voigt_K_exc_profile + voigt_K_ion_profile
                     + voigt_KL_ion_profile + voigt_KLL_ion_profile)
    #--------------------------------------------------------------------------

    #Apply transfer function and normalization
    #--------------------------------------------------------------------------
    #TF, profile, y0
    voigt_profile *= tf(x)                                                     #Profile follows the intensity dist. (normalized)

    #If in the norm_method the user chose no normalization, nothing is done
    if norm_method == 'Data': #and exp_data:                                     #If will look for the spot on x to force the theoretical line on the experimental line
      dif_exp = abs(hw_exp - norm_x)
      n_exp = np.argmin(dif_exp)
      dif_prof = abs(x - norm_x)
      n_prof = np.argmin(dif_prof)
      normalization = (counts_exp[n_exp]-y0[n_prof])/voigt_profile[n_prof]
      voigt_profile *= normalization                                           #Multiply itself by the normalization factor
    elif norm_method == 'Coor':                                                #It doesn't need exp_data, it will force theoretical line on (x,y) position
      dif_prof = abs(x - norm_x)
      n_prof = np.argmin(dif_prof)
      voigt_profile *= norm_y/voigt_profile[n_prof]
    #--------------------------------------------------------------------------

    #Chi2 calculation and plotting
    #--------------------------------------------------------------------------
    if exp_data:
      for i in range(len(hw_exp)):                                             #Restrict exp data to domain
        if hw_exp[i] >= min(x):
          exp_i = i
          break
      for j in range(len(hw_exp)):
        if hw_exp[-j-1] <= max(x):
          exp_f = len(hw_exp) - j
          break
      hw_exp = hw_exp[exp_i:exp_f]
      counts_exp = counts_exp[exp_i:exp_f]
      counts_exp_err = counts_exp_err[exp_i:exp_f]

      if y0_method == 'Chi':
        chi2_max = np.int(2e8)                                                 #High number. Will change to inf later.
        y = np.arange(min(counts_exp), np.average(counts_exp), step)
        for j in range(len(y)):                                                #Minimizes chi2 based on different offsets y0
          chi2 = np.int(0)
          spectrum_csd = voigt_profile + y[j]
          for i in range(len(hw_exp)):
            dif = abs(x - hw_exp[i])
            n = np.argmin(dif)
            chi2 += (spectrum_csd[n] - counts_exp[i])**2/(counts_exp_err[i]**2
                    *(len(hw_exp) - 1))
          if chi2 < chi2_max:
            chi2_max = chi2
            y0 = y[j]
          else: break                                                          #Makes computation quicker but can cause local minimi
        chi2 = chi2_max
        spectrum_csd = voigt_profile + y0                                      #Profile w/ noise, y0 can either come from the one calculated

        spectrum.plot(x, spectrum_csd, 'r',                                    #Plotting the spectrum and displaying the y0 and chi2
                      label = 'Theoretical spectrum, $y_0 = {}$, $\chi^2 = {}$'
                      .format(round(y0,2),round(chi2,2)))


      else:
        spectrum_csd = voigt_profile + y0                                      #Profile w/ noise y0 set by the user, chi2 calculation proceeds
        chi2 = np.int(0)
        for i in range(len(hw_exp)):
          dif = abs(x - hw_exp[i])
          n = np.argmin(dif)
          chi2 += (spectrum_csd[n] - counts_exp[i])**2/(counts_exp_err[i]**2   #If there's no exp error, it will just divide by 1
                  *(len(hw_exp) - 1))
        if y0_method == 'Set':                                                 #If it was an offset input with no m value then it displays only b
          spectrum.plot(x, spectrum_csd, 'r', label =
                        'Theoretical spectrum, $y_0= {}$, $\chi^2 = {}$'
                        .format(round(y0_b,2),round(chi2,2)))
        else:                                                                  #Else, if it was with the two points, then it shows y=mx+b
          spectrum.plot(x, spectrum_csd, 'r',
                        label = 'Theoretical spectrum, '+
                        '$y_0 = {}x + {}$, $\chi^2 = {}$'
                        .format(round(y0_m,4),round(y0_b,2), round(chi2,2)))

      spectrum.set_ylim(min(min(counts_exp), min(spectrum_csd)),               #No matter the options it sets the same plotting limits
                        max(max(spectrum_csd), max(counts_exp))*1.1)
    #--------------------------------------------------------------------------

    #Plotting without chi2
    #--------------------------------------------------------------------------
    else:                                                                      #No experimental data
      spectrum_csd = voigt_profile + y0                                        #Profile w/ noise
      if y0_method == 'Set':
        spectrum.plot(x, spectrum_csd, 'r', label =
                      'Theoretical spectrum, $y_0= {}$'.format(round(y0_b,2)))
      else:                                                                    #Only option y0_method == 'TwoPoints'
        spectrum.plot(x, spectrum_csd, 'r',
                      label = 'Theoretical spectrum, $y_0 = {}x + {}$'
                      .format(round(y0_m,2),round(y0_b,2)))

      spectrum.set_ylim(min(spectrum_csd)*.9, max(spectrum_csd)*1.1)
    #--------------------------------------------------------------------------

    #Plotting the contributions of each process (user request)
    #--------------------------------------------------------------------------
    if contributions:
      if kllion:
        spectrum_KLL_ion = voigt_KLL_ion_profile*tf(x) + y0
        spectrum.plot(x, spectrum_KLL_ion, 'y--',
                      label = 'KLL shells triple ionization contribution')
      if klion:
        spectrum_KL_ion = voigt_KL_ion_profile*tf(x) + y0
        spectrum.plot(x, spectrum_KL_ion, 'g--',
                      label = 'KL shells double ionization contribution')
      if kion:
        spectrum_K_ion = voigt_K_ion_profile*tf(x) + y0
        spectrum.plot(x, spectrum_K_ion, 'b--',
                      label = 'K shell ionization contribution')
      if kexc:
        spectrum_K_exc = voigt_K_exc_profile*tf(x) + y0
        spectrum.plot(x, spectrum_K_exc, 'c--',
                      label = 'K shell excitation contribution')

  else:                                                                        #Comes from the if simulation statement (other option is stick)
    #Stick plotting
    #--------------------------------------------------------------------------
    if contributions:
      if kexc:
        markerline = (
            spectrum.stem(hw_trans[hw_K_exc, 0], counts_K_exc, bottom=None,
                          basefmt=' ', markerfmt='co',linefmt='c-',
                          label = 'K shell excitation contribution'))
        markerline.set_markerfacecolor('none')
      if kion:
        markerline = (
            spectrum.stem(hw_trans[hw_K_ion, 0], counts_K_ion, bottom=None,
                          basefmt=' ', markerfmt='bo',linefmt='b-',
                          label = 'K shell ionization contribution'))
        markerline.set_markerfacecolor('none')
      if klion:
        markerline = (
            spectrum.stem(hw_trans[hw_KL_ion, 0], counts_KL_ion, bottom=None,
                          basefmt=' ', markerfmt='go',linefmt='g-',
                          label = 'KL shells double ionization contribution'))
        markerline.set_markerfacecolor('none')
      if kllion:
        markerline = (
            spectrum.stem(hw_trans[hw_KLL_ion, 0], counts_KLL_ion, bottom=None,
                          basefmt=' ', markerfmt='yo',linefmt='y-',
                          label = 'KLL shells triple ionization contribution'))
        markerline.set_markerfacecolor('none')

    #Stemplot with total consideration
    #--------------------------------------------------------------------------
    else:
      hw_stem = []
      counts_stem = []
      added = False
      if kexc:                                                                 #If kexc is considered, it adds without checking
        for i in range(len(hw_K_exc)):
          hw_stem.append(hw_K_exc[i])                                          #Adds index of the energy in the hw_trans database
          counts_stem.append(counts_K_exc[i])
      if kion:
        for i in range(len(hw_K_ion)):
          added = False
          for j in range(len(hw_stem)):                                        #For the first energy value kw_K_ion, it searches if there's already one in the hw_stem
            if hw_K_ion[i] == hw_stem[j]:                                      #Once one is found, it adds to the counts and breaks out the loop
              counts_stem[j] += counts_K_ion[i]
              added = True                                                     #Now it informs that something was added
              break                                                            #It stops looping
          if not added:                                                        #If the previous loop went on without being added, it appends new value
            hw_stem.append(hw_K_ion[i])
            counts_stem.append(counts_K_ion[i])
      if klion:                                                                #Same thought process for the following considerations
        for i in range(len(hw_KL_ion)):
          added = False
          for j in range(len(hw_stem)):
            if hw_KL_ion[i] == hw_stem[j]:
              counts_stem[j] += counts_KL_ion[i]
              added = True
              break
          if not added:
            hw_stem.append(hw_KL_ion[i])
            counts_stem.append(counts_KL_ion[i])
      if kllion:
        for i in range(len(hw_KLL_ion)):
          added = False
          for j in range(len(hw_stem)):
            if hw_KLL_ion[i] == hw_stem[j]:
              counts_stem[j] += counts_KLL_ion[i]
              added = True
              break
          if not added:
            hw_stem.append(hw_KLL_ion[i])
            counts_stem.append(counts_KLL_ion[i])

      for i in range(len(hw_stem)):                                            #Converts index values to the energies in the transition database
        hw_stem[i] = hw_trans[hw_stem[i], 0]

      markerline = (
          spectrum.stem(hw_stem, counts_stem, bottom=None, basefmt=' ',
                        markerfmt='ro',linefmt='r-', label = 'Stem plot'))
      markerline.set_markerfacecolor('none')

    #--------------------------------------------------------------------------

  #Setting limits and plotting details
  #----------------------------------------------------------------------------
  if not exp_data:
    hw_exp_middle = (hw_min + hw_max)/2
    hw_exp_size = hw_max - hw_min
    x_min = hw_exp_middle - hw_exp_size*1.1/2
    x_max = hw_exp_middle + hw_exp_size*1.1/2

  spectrum.legend(loc='upper right')
  spectrum.set_xlim(x_min, x_max)

  if scale == 'Log': spectrum.set_yscale('log')                                #If user input log
  fig.canvas.draw()
  #----------------------------------------------------------------------------
  return spectrum_csd                                                          #Returns for Levenberg-Marquardt fit purposes
#------------------------------------------------------------------------------
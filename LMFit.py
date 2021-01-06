# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 18:11:27 2019

Plots a spectrum with LM fit

@author: André Fernandes - afp.fernandes@campus.fct.unl.pt
"""

import numpy as np
from lmfit import minimize, Parameters, report_fit
from Calculation import spectrum_plot
from CalculationLite import spectrum_plot_lite
from tkinter import messagebox

def LM_spectrum_plot(hw_trans, cs_K_exc, cs_K_ion, cs_KL_ion, cs_KLL_ion,
                     spectrum, fig, csd, csv_file, trans_file, stick_sim,
                     scale, param_changed, csd_changed, tf_changed,
                     fraction_voigt, width_gauss, width_loren, step, y0_method,
                     y0_m, y0_b, norm_method, norm_x, norm_y, hw0, hw_min,
                     hw_max, fraction_mw, kexc, kion, klion, kllion,
                     contributions, metastates, fixed):

  #Starting inputs
  ions = list(csd.keys())[1:]                                                  #Removes the 'power' at the start
  csds = list(csd.values())[1:]
  fixeds = list(fixed.values())
  #----------------------------------------------------------------------------
  def LM_res(params, x, data, counts_exp_err):
    for i in ions:
      csd[i] = params['nq' + i]

    #Uses lite as it is quicker to compute
    spectrum_csd = spectrum_plot_lite(data, counts_exp_err, data[0], y0,
                                      hw_trans, cs_K_exc, cs_K_ion, y,
                                      cs_KL_ion, cs_KLL_ion, csd, trans_file,
                                      fraction_voigt, width_gauss, width_loren,
                                      y0_method, y0_m, y0_b, hw0, fraction_mw,
                                      kexc, kion, klion, kllion, metastates)

    result = (data[1] - spectrum_csd)/counts_exp_err

    return result
  #----------------------------------------------------------------------------

  params = Parameters()
  for i in range(len(ions)):
    #params.add('nq' + ions[i], value = csds[i], min = csds[i] * .5,
     #          max = csds[i] * 2, vary = not fixeds[i])                       #Added limit of 10 times the cds guess value or 10% because sometimes it would go to really high values
    params.add('nq' + ions[i], value = csds[i], min = 0,#csds[i] * .5,         #Removed the max limit because restrictions are bad for computation time
               vary = not fixeds[i])                                           #Needs to be added an expression to have a gaussian distribution within the parameters

  #----------------------------------------------------------------------------
  x = np.arange(hw_min, hw_max, step)                                          #Domain arranged

  data_exp = np.genfromtxt(csv_file, delimiter = ',')
  hw_exp = data_exp[:,0]
  counts_exp = data_exp[:,1]
  try:
    counts_exp_err = data_exp[:,3]
  except Exception as ex:
    print('Could not load count error bars. Will be shown as the square root of experimental counts. ' + ex)
    counts_exp_err = np.sqrt(counts_exp)

  for i in range(len(hw_exp)):                                                 #Restrict exp data to domain
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
  data = [hw_exp, counts_exp]

  #----------------------------------------------------------------------------
  #In the case of using SPlotLite
  y0 = np.array(y0_m * hw_exp + y0_b)                                          #Offset as an array to simply add to profile
  y = np.arange(min(data[1]), np.average(data[1]), step)

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
    print('Could not load the transition energies database. ' + ex)
    messagebox.showerror('Transition energies database', 'The transition ener'+
                         'gies database was incorrectly loaded. The calculati'+
                         'on will stop.')
    return
  #----------------------------------------------------------------------------
  #----------------------------------------------------------------------------

  out = minimize(LM_res, params, method = 'leastsq', args =(x, data,
                                                            counts_exp_err))
  report_fit(out)                                                              #Analyze in console the outcome
  values = out.params.valuesdict()

  for i in ions:
    csd[i] = values['nq' + i]


  spectrum_plot(hw_trans, cs_K_exc, cs_K_ion, cs_KL_ion, cs_KLL_ion, spectrum,
                fig, csd, csv_file, trans_file, stick_sim, scale,
                param_changed, csd_changed, tf_changed, fraction_voigt,
                width_gauss, width_loren, step, y0_method, y0_m, y0_b,
                norm_method, norm_x, norm_y, hw0, hw_min, hw_max, fraction_mw,
                kexc, kion, klion, kllion, contributions, metastates)
  return csd
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 18:11:27 2019

Function to load the database

@author: André Fernandes - afp.fernandes@campus.fct.unl.pt
"""
import os
import pandas as pd                                                            #To deal with the csv files with dataframe
from ast import literal_eval
from tkinter import messagebox

#This code loads the database
#------------------------------------------------------------------------------
def dbloader(path1, Z, kexc, kion, klion, kllion, total):
  #Checks if the database files exists w/ transition energies & cross sections
  hw_trans_exists = os.path.isfile(str(path1+str(Z)+'/'+str(Z)+
                                       '-transitions.csv'))
  cs_K_exc_exists = os.path.isfile(str(path1+str(Z)+'/'+str(Z)+
                                       '-cs_K_exc.csv'))
  cs_K_ion_exists = os.path.isfile(str(path1+str(Z)+'/'+str(Z)+
                                       '-cs_K_ion.csv'))
  cs_KL_ion_exists = os.path.isfile(str(path1+str(Z)+'/'+str(Z)+
                                        '-cs_KL_ion.csv'))
  cs_KLL_ion_exists = os.path.isfile(str(path1+str(Z)+'/'+str(Z)+
                                         '-cs_KLL_ion.csv'))

  cs_K_exc = []
  cs_K_ion = []
  cs_KL_ion = []
  cs_KLL_ion = []

  if hw_trans_exists:
    hw_trans = pd.read_csv(str(path1+str(Z)+'/'+str(Z)+'-transitions.csv'),
                           sep=',', header=None).values
  else:
    messagebox.showerror('Transition energies database missing or incorrect',
                         'The transition energy database for this element is '+
                         'not present or incorrectly put. Make sure it is nam'+
                         'ed {}-transitions.csv in the {} folder.'.format(Z,Z))
    return
  if cs_K_exc_exists:
    cs_K_exc = pd.read_csv(str(path1+str(Z)+'/'+str(Z)+'-cs_K_exc.csv'),
                           sep=', ', header=None, engine = 'python')           #engine = 'python' gets rid of a warning
    cs_K_exc[[6,7]] = cs_K_exc[[6,7]].applymap(literal_eval)                   #Converts the columns into lists
    cs_K_exc = cs_K_exc.values                                                 #Converts dataframe into an array (due to how the program was first made)
    kexc = True                                                                #K-shell excitation
  else:
    messagebox.showwarning('No K shell excitation cross section database',
                           'The cross section values for the K shell excitati'+
                           'on database are either not present or incorrectly'+
                           ' put. Make sure it is named '+
                           '{}-cs_K_exc.csv in the {} folder.'.format(Z,Z))
    kexc = False
    total = False

  if cs_K_ion_exists:
    cs_K_ion = pd.read_csv(str(path1+str(Z)+'/'+str(Z)+'-cs_K_ion.csv'),
                           sep=', ', header=None, engine = 'python')
    cs_K_ion[[4,5]] = cs_K_ion[[4,5]].applymap(literal_eval)
    cs_K_ion = cs_K_ion.values
    kion = True                                                                #K-shell ionization
  else:
    messagebox.showwarning('No K shell ionization cross section database',
                           'The cross section values for the K shell ionizati'+
                           'on database are either not present or incorrectly'+
                           ' put. Make sure it is named '+
                           '{}-cs_K_ion.csv in the {} folder.'.format(Z,Z))
    kion = False
    total = False

  if cs_KL_ion_exists:
    cs_KL_ion = pd.read_csv(str(path1+str(Z)+'/'+str(Z)+'-cs_KL_ion.csv'),
                            sep=', ', header=None, engine = 'python')
    cs_KL_ion[[4,5]] = cs_KL_ion[[4,5]].applymap(literal_eval)
    cs_KL_ion = cs_KL_ion.values
    klion = True                                                               #KL-shell ionization
  else:
    messagebox.showwarning('No KL shells double ionization cross section data'+
                           'base', 'The cross section values for the KL shell'+
                           's double ionization database are either not prese'+
                           'nt or incorrectly put. Make sure it is named '+
                           '{}-cs_KL_ion.csv in the {} folder.'.format(Z,Z))
    klion = False
    total = False

  if cs_KLL_ion_exists:
    cs_KLL_ion = pd.read_csv(str(path1+str(Z)+'/'+str(Z)+'-cs_KLL_ion.csv'),
                             sep=', ', header=None, engine = 'python')
    cs_KLL_ion[[4,5]] = cs_KLL_ion[[4,5]].applymap(literal_eval)
    cs_KLL_ion = cs_KLL_ion.values
    kllion = True                                                              #KLL-shell ionization
  else:
    messagebox.showwarning('No KLL shells triple ionization cross section dat'+
                           'abase', 'The cross section values for the KLL she'+
                           'lls triple ionization database are either not pre'+
                           'sent or incorrectly put. Make sure it is named '+
                           '{}-cs_KLL_ion.csv in the {} folder.'.format(Z,Z))
    kllion = False
    total = False
  return(kexc, kion, klion, kllion, total, hw_trans, cs_K_exc, cs_K_ion,
         cs_KL_ion, cs_KLL_ion, cs_K_exc_exists, cs_K_ion_exists,
         cs_KL_ion_exists, cs_KLL_ion_exists)
#------------------------------------------------------------------------------
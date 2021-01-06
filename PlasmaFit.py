# -*- coding: utf-8 -*-
"""
Created on Sun Aug  4 20:52:44 2019

Plasma analysis main program

@author: André Fernandes - fernandesafp@gmail.com;
                           afp.fernandes@campus.fct.unl.pt
"""
import csv, os, sys
from numpy import genfromtxt, ones, asarray
from tkinter import (Tk, ttk, IntVar, DoubleVar, BooleanVar, StringVar,
                     Toplevel, Frame, messagebox, OptionMenu, Checkbutton,
                     Label, Entry, LEFT, RIGHT, TOP, BOTTOM, X, YES, END, BOTH)
from tkinter.filedialog import askopenfilename
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                               NavigationToolbar2Tk)

#André Fernandes 2018/2019
from Load import dbloader                                                      #To load the databases
from Calculation import spectrum_plot                                          #To use functions for the plotting
from LMFit import LM_spectrum_plot                                             #To perform the LM fit
from mpl_toolkits.axes_grid1.inset_locator import inset_axes                   #To have a subplot in the corner of the main plot

#Global variable for the mouse click coordinates
coords = []
#Global variables for the CSD guess
#------------------------------------------------------------------------------
i=2
rows = []
power = []

csd_standard = {'power': 16, '11': 8, '12': 4, '13': 2, '14': 1, '15': .19,    #Arbitrary values to input for the user
                '16': .021, '13m': 1, '14m': 1}
fixed_standard = {'11': 0, '12': 0, '13': 0, '14': 0, '15': 0, '16': 0,
                  '13m': 1, '14m': 1}

csd = csd_standard
fixed = fixed_standard
#------------------------------------------------------------------------------

if __name__ == "__main__":
  Z = 18
  path1=os.path.dirname(sys.argv[0])+'/'

  sim = Tk()
  sim.title('Plasma diagnostics')
  #sim.iconbitmap(r'atom.ico')                                                 #Optional .ico (does not work in my Ubuntu)

  kexc_var = BooleanVar()
  kion_var = BooleanVar()
  klion_var = BooleanVar()
  kllion_var = BooleanVar()
  total_var = BooleanVar()
  total_var.set(True)                                                          #Consider all

  #Database loader
  db = dbloader(path1, Z, kexc_var.get(), kion_var.get(), klion_var.get(),
                kllion_var.get(), total_var.get())
  if db == None: sys.exit()                                                    #Does not open without a transition database
  #Database assignment
  kexc_var.set(db[0])
  kion_var.set(db[1])
  klion_var.set(db[2])
  kllion_var.set(db[3])
  total_var.set(db[4])
  (hw_trans, cs_K_exc, cs_K_ion, cs_KL_ion, cs_KLL_ion, cs_K_exc_exists,
   cs_K_ion_exists, cs_KL_ion_exists, cs_KLL_ion_exists) = db[5:]

  tf_var = StringVar()
  tf_var.set('')                                                               #Transfer function file location
  load_var = StringVar()
  load_var.set('')                                                             #Load file location

  metastates_var = BooleanVar()
  metastates_var.set(True)
  contributions_var = BooleanVar()
  contributions_var.set(False)                                                 #Ask to plot each process' contribution
  plot_done = BooleanVar()                                                     #To know if there is something on that was calculated
  plot_done.set(False)
  param_changed = BooleanVar()                                                 #Know if parameters were changed/saved by the user
  param_changed.set(False)
  csd_changed = BooleanVar()                                                   #Know if charge-state densities were changed/saved by the user
  csd_changed.set(False)
  tf_changed = BooleanVar()                                                    #Know if user added a transfer function
  tf_changed.set(False)

  #Default parameters
  param_voigt = DoubleVar()
  param_voigt.set(.56)
  param_gauss = DoubleVar()
  param_gauss.set(.62)
  param_loren = DoubleVar()
  param_loren.set(.19)
  param_step = DoubleVar()
  param_step.set(.01)
  param_hwoff = DoubleVar()
  param_hwoff.set(0)
  param_min = DoubleVar()
  param_min.set(3087.94)
  param_max = DoubleVar()
  param_max.set(3118.73)
  param_mw = DoubleVar()
  param_mw.set(.99)

  #Default offset and normalization values
  param_m = DoubleVar()
  param_m.set(0)
  param_b = DoubleVar()
  param_b.set(.7)
  param_norm_x = DoubleVar()
  param_norm_x.set(3104.161)
  param_norm_y = DoubleVar()
  param_norm_y.set(6.213744)
  check_off = StringVar()
  check_off.set('Set')
  check_norm = StringVar()
  check_norm.set('None')

  frame = Frame(sim)
  frame.grid(row=0,column=0, sticky='nsew')                                    #Test

  f = Figure(figsize=(10,5), dpi=125)                                          #Spectrum canvas
  spectrum = f.add_subplot(111)                                                #Spectrum placement
  spectrum.set_xlabel('Energy (eV)')
  spectrum.set_ylabel('Counts/s')
  ax_tf = inset_axes(spectrum, width = '12.5%', height = '12.5%', loc = 2)     #Transfer function canvas
  ax_tf.axis('off')

  f.tight_layout()                                                             #To have a better look

  figure_frame = Frame(sim)                                                    #Plot frame
  figure_frame.grid(row=0,column=0,rowspan=2,columnspan=2, sticky='nesw')      #Plot frame
  canvas = FigureCanvasTkAgg(f, master=figure_frame)

  canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

  toolbar_frame = Frame(sim)
  toolbar_frame.grid(row=4,column=0,columnspan=3, sticky='w') 
  #toolbar = NavigationToolbar2TkAgg( canvas, toolbar_frame )
  NavigationToolbar2Tk( canvas, toolbar_frame )

  def LM_plotting():                                                           #Goes on to plot with Levenberg-Marquardt algorithm
    if load_var.get() != '' and choice_var.get() == 'Sim.':                    #Only plots with simulation and with a file to compare
      progress_var.set(25)
      if not param_changed.get():
        messagebox.showinfo('Default parameters.', 'Parameters were not set, '+
                            'using default.')
        param_changed.set(True)
        progress_var.set(30)

      if not csd_changed.get():
        messagebox.showinfo('Default charge state density.', 'Charge-state de'+
                            'nsities were not set, using default based on def'+
                            'ault database.')
        csd_changed.set(True)
        progress_var.set(40)
      if not tf_changed.get():
        messagebox.showinfo('Linear transfer function.', 'No transfer functio'+
                            'n loaded, using a flat intensity distribution. T'+
                            'o load one, click "Intensity dist."')
        tf_changed.set(True)
        progress_var.set(50)
      #if check_norm.get() != 'None':                                           #Normalization is not good for lmfit
       #messagebox.showinfo('New parameter.', 'No normalization must be set f'+
        #                   'or Levenberg-Marquardt fit.')
        #progress_var.set(60)
        #check_norm.set('None')
      progress_var.set(75)
      global csd
      csd = LM_spectrum_plot(hw_trans, cs_K_exc, cs_K_ion, cs_KL_ion,          #Goes on to plot and save charge-state densities
                             cs_KLL_ion, spectrum, f, csd, load_var.get(),
                             tf_var.get(), choice_var.get(), scale_var.get(),
                             param_changed.get(), csd_changed.get(),
                             tf_changed.get(), param_voigt.get(),
                             param_gauss.get(), param_loren.get(),
                             param_step.get(), check_off.get(), param_m.get(),
                             param_b.get(), check_norm.get(),
                             param_norm_x.get(), param_norm_y.get(),
                             param_hwoff.get(), param_min.get(),
                             param_max.get(), param_mw.get(), kexc_var.get(),
                             kion_var.get(), klion_var.get(), kllion_var.get(),
                             contributions_var.get(), metastates_var.get(),
                             fixed)
      progress_var.set(100)
      plot_done.set(True)
      messagebox.showinfo('LM fit done.','Levenberg-Marquardt fit done, re-op'+
                          'en the CSD window to see the new charge-state dens'+
                          'ities.')
      progress_var.set(0)
    elif load_var.get() == '':                                                 #If there is no file
      messagebox.showerror('No experimental data.','It is impossible to apply'+
                           ' the Levenberg-Marquardt algorithm without experi'+
                           'mental data.')
    else:                                                                      #Else choice_var as stickplot, it forces to simulation
      progress_var.set(5)
      choice_var.set('Sim.')
      LM_plotting()

  def plotting():                                                              #Goes on to calculate and plot
    progress_var.set(25)

    if not param_changed.get():
      messagebox.showinfo('Default parameters.', 'Parameters were not set, us'+
                          'ing default.')
      param_changed.set(True)
      progress_var.set(30)

    if not csd_changed.get():
      messagebox.showinfo('Default charge state density.', 'Charge state dens'+
                          'ities were not set, using default based on default'+
                          ' database.')
      csd_changed.set(True)
      progress_var.set(40)

    if not tf_changed.get():
      messagebox.showinfo('Linear transfer function.', 'No transfer function '+
                          'loaded, using a flat intensity distribution. To lo'+
                          'ad one, click "Intensity dist."')
      tf_changed.set(True)
      progress_var.set(50)
    progress_var.set(75)
    spectrum_plot(hw_trans, cs_K_exc, cs_K_ion, cs_KL_ion, cs_KLL_ion,
                  spectrum, f, csd, load_var.get(), tf_var.get(),
                  choice_var.get(), scale_var.get(), param_changed.get(),
                  csd_changed.get(), tf_changed.get(), param_voigt.get(),
                  param_gauss.get(), param_loren.get(), param_step.get(),
                  check_off.get(), param_m.get(), param_b.get(),
                  check_norm.get(), param_norm_x.get(), param_norm_y.get(),
                  param_hwoff.get(), param_min.get(), param_max.get(),
                  param_mw.get(), kexc_var.get(), kion_var.get(),
                  klion_var.get(), kllion_var.get(), contributions_var.get(),
                  metastates_var.get())
    progress_var.set(100)
    plot_done.set(True)
    progress_var.set(0)

  def _quit():
    sim.quit()                                                                 #Stops mainloop
    sim.destroy()                                                              #This is necessary on Windows to prevent
                                                                               #Fatal Python Error: PyEval_RestoreThread: NULL tstate

  def load_data():                                                             #Gets the path and plots the experimental data
    fname = askopenfilename(filetypes=(('Energy, intensity, energy error, int'+
                                        'ensity error', '*.txt;*.csv'),
      ('All files', '*.*')))
    if fname != '':                                                            #If a file as chosen, then it proceeds
      if load_var.get() != '':                                                 #If there is a file loaded already if clears
        spectrum.clear()
        spectrum.set_xlabel('Energy (eV)')
        spectrum.set_ylabel('Counts/s')

      if plot_done.get():                                                      #If there's a calculation, it runs again trying to show with load file
        try:
          load_var.set(fname)
          choice_var.set('Sim.')                                               #It plots data if it isn't in stickplot mode
          plotting()
        except Exception as ex:
          print('Could not open and plot file. ' + ex)
          messagebox.showerror('Invalid file.', 'It was not possible to open '+
                               'and plot {}.'.format(fname))
          return

      else:
        try:                                                                   #Tries to load experimental data
          data = genfromtxt(fname, delimiter=',')
          hw_exp = data[:,0]
          counts_exp = data[:,1]
          spectrum.plot(hw_exp, counts_exp, 'ko', linewidth = .75,
                        markersize = 3, label ='Experimental data', mfc='none',
                        zorder=1)
        except Exception as ex:
          print('Could not open and plot file. ' + ex)
          messagebox.showerror('Invalid file.', 'It was not possible to open '+
                               'and plot {}.'.format(fname))
          return

        try:                                                                   #Attempt to find errorbars as they're not mandatory
          hw_exp_err = data[:,2]
          counts_exp_err = data[:,3]
          spectrum.errorbar(hw_exp, counts_exp, yerr = counts_exp_err,
                            xerr = hw_exp_err, fmt = 'none', ecolor='k',
                            elinewidth = 0.5, capsize = 1, zorder=1)
        except Exception as ex:
          print('Could not load count error bars. They will not be shown.' + ex)
          counts_exp_err = ones(len(hw_exp))                                   #For chi2 display without errors

        load_var.set(fname)

        hw_exp_middle = (max(hw_exp) + min(hw_exp))/2                          #Axis parameters
        hw_exp_size = max(hw_exp) - min(hw_exp)
        x_min = hw_exp_middle - hw_exp_size*1.1/2
        x_max = hw_exp_middle + hw_exp_size*1.1/2

        spectrum.set_xlim(x_min, x_max)
        if scale_var.get() == 'Log': spectrum.set_yscale('log')
        spectrum.set_ylim(min(counts_exp), max(counts_exp)*1.1)
        spectrum.legend()

        f.canvas.draw()                                                        #Shows the plot

  def ask_plot():                                                              #If there is a plot, it redoes
    if plot_done.get(): plotting()

  def load_tf():                                                               #Gets the path of the transfer function file
    fname = askopenfilename(filetypes=(('Energy, intensity', '*.txt;*.csv'),
                                       ('All files', '*.*')))
    if fname != '':                                                            #If something was loaded, then it plots
      try:
        if tf_var.get() != '':
          ax_tf.clear()
          ax_tf.axis('off')
        data = genfromtxt(fname, delimiter=',')
        hw_tf = data[:,0]
        counts_tf = data[:,1]
        ax_tf.plot(hw_tf, counts_tf, 'ko', markersize = .5)                    #It plots in the corner of the canvas
      except:
        messagebox.showerror('Invalid file.', 'It was not possible to open an'+
                             'd plot {}.'.format(fname))
        return

      tf_var.set(fname)
      tf_changed.set(True)
      ask_plot()

      hw_tf_middle = (max(hw_tf) + min(hw_tf))/2                               #For title placement
      ax_tf.text(hw_tf_middle, max(counts_tf), 'Transfer function',
                 fontsize = 7, ha = 'center')

      f.canvas.draw()                                                          #Shows the plot

  def off_two_grabber():
    check_off.set('TwoPoints')                                                 #This forces equation to set as default in case user doesn't choose two points
    if plot_done.get() or load_var.get() != '':                                #If there is plot or exp data, you can click to get inclination
      try:
        for t in spectrum.texts: t.set_visible(False)                          #Removes x marks to prepare for new offset
        spectrum.grid(b=True, which='both', axis='both')
        f.canvas.draw()

        def onclick(event):                                                    #On click it saves the coordinates
          ix, iy = float(event.xdata), float(event.ydata)
          global coords
          coords.append((ix,iy))

          spectrum.text(ix, iy, 'x', ha='center', va='center')                 #Draws an 'x'
          f.canvas.draw()

          if len(coords) == 2:                                                 #Once it saves two clicks, disconnects the mouse
            f.canvas.mpl_disconnect(cid)
            spectrum.grid(False)                                               #Gets rid of grids and saves coordenates
            spectrum.grid(False, which='minor')
            f.canvas.draw()
            x1 = coords[0][0]
            y1 = coords[0][1]
            x2 = coords[1][0]
            y2 = coords[1][1]
            param_m.set((y2 - y1)/(x2 - x1))
            param_b.set(y2 - param_m.get() * x2)
            for t in spectrum.texts: t.set_visible(False)                      #Gets rid of the 'x's
            coords = []
            ask_plot()

        cid = f.canvas.mpl_connect('button_press_event', onclick)
      except Exception as ex:
        print('Exception thrown whilst obtaining/applying the offset.' + ex)
        messagebox.showerror('Error', 'Something went wrong with the offset.')

    else:                                                                      #If there is no plot or no data, it opens a new window
      fields = ['m', 'b']
      defaults = [0, .7]

      def quit_window():
        root.destroy()

      def fetch(entries):
        val = []
        for entry in entries:
          try:
            val.append(float(entry.get()))
          except:
            messagebox.showerror('Invalid parameter', 'One input is not a num'+
                                 'ber. Parameters will reset.')
            redo_form(entries)
        param_m.set(val[0])
        param_b.set(val[1])
        ask_plot()

      def make_form(root, fields):
        entries = []
        row = Frame(root)
        lab = Label(row, width=36, text='With no plot, choose y=mx+b paramete'+
                    'rs.', anchor='w')
        row.pack(side=TOP, fill=X, padx=5, pady=5)
        lab.pack(side=LEFT)

        row = Frame(root)
        lab = Label(row, width=5, text=fields[0], anchor='w')
        ent = Entry(row)
        ent.insert(0, defaults[0])
        row.pack(side=TOP, fill=X, padx=5, pady=5)
        lab.pack(side=LEFT)
        ent.pack(side=RIGHT, expand=YES, fill=X)
        entries.append(ent)

        row = Frame(root)
        lab = Label(row, width=5, text=fields[1], anchor='w')
        ent = Entry(row)
        ent.insert(0, defaults[1])
        row.pack(side=TOP, fill=X, padx=5, pady=5)
        lab.pack(side=LEFT)
        ent.pack(side=RIGHT, expand=YES, fill=X)
        entries.append(ent)

        return entries

      def redo_form(entries):
        i = 0
        for entry in entries:
          entry.delete(0, END)
          entry.insert(0, defaults[i])
          i += 1
        fetch(entries)

      #Offset submain
      #------------------------------------------------------------------------
      root = Toplevel()
      root.title('Offset points')
      #root.iconbitmap(r'favico.ico')
      ents = make_form(root, fields)
      b1 = ttk.Button(root, text='Save', command=(lambda e=ents: fetch(e)))
      b1.pack(side=LEFT, padx=5, pady=5)
      root.bind('<Return>', lambda events, e=ents: fetch(e))                   #If presses enter, it saves ('presses' b1 button)
      b2 = ttk.Button(root, text='Reset', command=(lambda e=ents:redo_form(e)))
      b2.pack(side=LEFT, padx=5, pady=5)
      b3 = ttk.Button(root, text='Quit', command=quit_window)
      b3.pack(side=LEFT, padx=5, pady=5)

  def norm_grabber():
    check_norm.set('Coor')                                                     #This forces equation to set as default in case user doesn't choose two points
    if plot_done.get() or load_var.get() != '':
      try:
        for t in spectrum.texts: t.set_visible(False)                          #Removes x marks to prepare for new offset
        spectrum.grid(b=True, which='both', axis='both')
        f.canvas.draw()

        def onclick(event):
          ix, iy = float(event.xdata), float(event.ydata)
          spectrum.text(ix, iy, 'x', ha='center', va='center')
          f.canvas.draw()

          global coords
          coords.append((ix,iy))
          param_norm_x.set(coords[0][0])
          param_norm_y.set(coords[0][1])

          f.canvas.mpl_disconnect(cid)
          spectrum.grid(b=False)
          spectrum.grid(False, which='minor')
          f.canvas.draw()
          for t in spectrum.texts: t.set_visible(False)
          coords = []
          ask_plot()

        cid = f.canvas.mpl_connect('button_press_event', onclick)
      except Exception as ex:
        print('Exception thrown whilst obtaining/applying normalization.' + ex)
        messagebox.showerror('Error', 'Something went wrong with the normaliz'+
                             'ation.')

    else:                                                                      #If no plot or exp data, it opens a window for input
      fields = ['x', 'y']
      defaults = [3104.161, 6.213744]

      def quit_window():
        root.destroy()

      def fetch(entries):
        val = []
        for entry in entries:
          try:
            val.append(float(entry.get()))
          except:
            messagebox.showerror('Invalid parameter', 'One input is not a num'+
                                 'ber. Parameters will reset.')
            redo_form(entries)
        param_norm_x.set(val[0])
        param_norm_y.set(val[1])
        ask_plot()

      def make_form(root, fields):
        entries = []
        row = Frame(root)
        lab = Label(row, width=40, text='With no plot, choose coordinates to '+
                    'normalize.', anchor='w')
        row.pack(side=TOP, fill=X, padx=5, pady=5)
        lab.pack(side=LEFT)

        row = Frame(root)
        lab = Label(row, width=5, text=fields[0], anchor='w')
        ent = Entry(row)
        ent.insert(0, defaults[0])
        row.pack(side=TOP, fill=X, padx=5, pady=5)
        lab.pack(side=LEFT)
        ent.pack(side=RIGHT, expand=YES, fill=X)
        entries.append(ent)

        row = Frame(root)
        lab = Label(row, width=5, text=fields[1], anchor='w')
        ent = Entry(row)
        ent.insert(0, defaults[1])
        row.pack(side=TOP, fill=X, padx=5, pady=5)
        lab.pack(side=LEFT)
        ent.pack(side=RIGHT, expand=YES, fill=X)
        entries.append(ent)

        return entries

      def redo_form(entries):
        i = 0
        for entry in entries:
          entry.delete(0, END)
          entry.insert(0, defaults[i])
          i += 1
        fetch(entries)

      #Normalizer submain
      #------------------------------------------------------------------------
      root = Toplevel()
      root.title('Normalizer point')
      #root.iconbitmap(r'favico.ico')
      ents = make_form(root, fields)
      b1 = ttk.Button(root, text='Save', command=(lambda e=ents: fetch(e)))
      b1.pack(side=LEFT, padx=5, pady=5)
      root.bind('<Return>', lambda events, e=ents: fetch(e))                   #If presses enter, it saves ('presses' b1 button)
      b2 = ttk.Button(root, text='Reset', command=(lambda e=ents:redo_form(e)))
      b2.pack(side=LEFT, padx=5, pady=5)
      b3 = ttk.Button(root, text='Quit', command=quit_window)
      b3.pack(side=LEFT, padx=5, pady=5)

  def plot_parameters():
    def quit_window():
      root.destroy()

    fields = ('Voigt fraction', 'Gaussian width', 'Lorentzian width',
              'Energy axis step', 'Energy offset', 'Min. Energy',
              'Max. Energy', 'Maxwellian dist.', 'Background noise',
              'Normalization')
    defaults = asarray([.56, .62, .19, .01, 0, 3087.94, 3118.73, .99, .7,
                        3104.161])

    def fetch(entries):
      val = []
      for entry in entries:
        try:
          val.append(float(entry.get()))
        except:
          messagebox.showerror('Invalid parameter', 'One input is not a numbe'+
                               'r. Parameters will reset.')
          redo_form(entries)

      param_changed.set(True)

      if val[0] >= 0 and val[0] <= 1:
        param_voigt.set(val[0])
      else:
        messagebox.showwarning('Invalid parameter.', 'Voigt fraction must be '+
                               'between 0 and 1. 1 being a fully Gaussian dis'+
                               'tribution and 0 a fully Lorentzian one. Value'+
                               ' was not changed.')
      if val[1] >= 0:
        param_gauss.set(val[1])
      else:
        messagebox.showwarning('Invalid parameter.', 'Gaussian width must be '+
                               'greater than 0. Value was not changed.')
      if val[2] >= 0:
        param_loren.set(val[2])
      else:
        messagebox.showwarning('Invalid parameter.', 'Lorentzian width must b'+
                               'e greater than 0. Value was not changed.')
      if val[3] > 0:
        param_step.set(val[3])
      else:
        messagebox.showwarning('Invalid parameter.', 'Step must be greater th'+
                               'an 0. Value was not changed.')
      param_hwoff.set(val[4])
      if val[5] < val[6]:
        param_min.set(val[5])
        param_max.set(val[6])
      else:
        messagebox.showwarning('Invalid parameters.', 'The energy interval is'+
                               ' not valid. The minimum value is greater than'+
                               ' the maximum value. Values were not changed.')
      if val[7] >= 0 and val[7] <= 1:
        param_mw.set(val[7])
      else:
        messagebox.showwarning('Invalid parameter.', 'The Maxwellian vs. non-'+
                               'Maxwellian distribution fraction must be betw'+
                               'een 0 and 1. 1 being a fully Maxwellian distr'+
                               'ibution and 0 a fully non-Maxwellian one. Val'+
                               'ue was not changed.')
      if check_off.get() == 'Set':
        param_m.set(0)
        param_b.set(ent_off.get())
      if check_norm.get() == 'Data':
        param_norm_x.set(ent_norm.get())
      ask_plot()

    def off_set():
      if check_off.get() == 'Set':
        ent_off.configure(state='normal')
        param_m.set(0)
        param_b.set(ent_off.get())
        ask_plot()
      elif check_off.get() == 'Chi':
        if load_var.get() == '':
          messagebox.showwarning('Warning: No data', 'Missing data to plot wi'+
                                 'th offset chi2 minimization.')
          check_off.set('Set')
          off_set()
        else:
          ent_off.configure(state='disabled')
          ask_plot()
      else:                                                                    #Else two point
        ent_off.configure(state='disabled')
        off_two_grabber()

    def norm():
      if check_norm.get() == 'None':
        ent_norm.configure(state='disabled')
        ask_plot()
      elif check_norm.get() == 'Data':
        if load_var.get() == '':
          messagebox.showwarning('Warning: No data', 'Missing data plot to no'+
                                 'rmalize from.')
          check_norm.set('None')
          norm()
        else:
          ent_norm.configure(state='normal')
          param_norm_x.set(ent_norm.get())
          ask_plot()
      else:                                                                    #Else coordinate
        ent_norm.configure(state='disabled')
        norm_grabber()

    def make_form(root, fields):
      entries = []
      i = 0
      for field in fields[:8]:
        row = Frame(root)
        lab = Label(row, width=15, text=field, anchor='w')
        ent = Entry(row)
        ent.insert(0, defaults[i])
        i += 1
        row.pack(side=TOP, fill=X, padx=5, pady=5)
        lab.pack(side=LEFT)
        ent.pack(side=RIGHT, expand=YES, fill=X)
        entries.append(ent)
      return entries

    def redo_form(entries):
      i = 0
      for entry in entries:
        if i < 8:
          entry.delete(0, END)
          entry.insert(0, defaults[i])
          i += 1
        elif i == 8:
          ent_off.delete(0, END)
          ent_off.insert(0, defaults[i])
          i += 1
        elif i == 9:
          ent_norm.delete(0, END)
          ent_norm.insert(0, defaults[i])
      fetch(entries)

    #Parameters main
    #--------------------------------------------------------------------------
    root = Toplevel()
    root.title('Parameters')
    #root.iconbitmap(r'favico.ico')
    ents = make_form(root, fields)

    row = Frame(root)                                                          #For the offset parameter
    lab = Label(row, width=15, text=fields[8], anchor='w')
    row.pack(side=TOP, fill=X, padx=5, pady=5)
    lab.pack(side=LEFT)
    ent_off = Entry(row, width = 8)
    ent_off.insert(0, defaults[8])
    c1 = ttk.Radiobutton(row, text='Fixed', variable=check_off, value='Set',
                         command=off_set)
    c2 = ttk.Radiobutton(row, text='Min. \u03c7^2', variable=check_off,
                         value='Chi', command=off_set)
    c3 = ttk.Radiobutton(row, text='Two points', variable=check_off,
                         value='TwoPoints', command=off_set)
    c3.pack(side=RIGHT, fill=X)
    c2.pack(side=RIGHT, fill=X)
    ent_off.pack(side=RIGHT)
    c1.pack(side=RIGHT, fill=X)
    ents.append(ent_off)

    row = Frame(root)                                                          #For the normalization parameter
    lab = Label(row, width=15, text=fields[9], anchor='w')
    row.pack(side=TOP, fill=X, padx=5, pady=5)
    lab.pack(side=LEFT)
    ent_norm = Entry(row, width = 8)
    ent_norm.insert(0, defaults[9])
    d1 = ttk.Radiobutton(row, text='None', variable=check_norm, value='None',
                         command=norm)
    d2 = ttk.Radiobutton(row, text='Data point',variable=check_norm,
                         value='Data', command=norm)
    d3 = ttk.Radiobutton(row, text='Coordinate',variable=check_norm,
                         value='Coor', command=norm)
    d3.pack(side=RIGHT, fill=X)
    ent_norm.pack(side=RIGHT)
    d2.pack(side=RIGHT, fill=X)
    d1.pack(side=RIGHT, fill=X)
    ents.append(ent_norm)
    norm()                                                                     #So that window checks proper radiobuttons

    b1 = ttk.Button(root, text='Save',
                    command=(lambda e=ents: fetch(e)))
    b1.pack(side=LEFT, padx=5, pady=5)
    root.bind('<Return>', lambda events, e=ents: fetch(e))                     #If presses enter, it saves ('presses' b1 button)
    b2 = ttk.Button(root, text='Reset',
                    command=(lambda e=ents:redo_form(e)))
    b2.pack(side=LEFT, padx=5, pady=5)
    b3 = ttk.Button(root, text='Quit', command=quit_window)
    b3.pack(side=LEFT, padx=5, pady=5)

    root.mainloop()
    #--------------------------------------------------------------------------

  def csd_guess():                                                             #Maybe insert parameter here to tell to write LM fit
    def quit_window():
      root.destroy()

    def fetch():
      global csd, fixed
      try:
        csd = {}                                                               #Reset the csd to plot
        fixed = {}
        val = []
        for row in rows:
          j = 0
          for column in row:
            if j == 0:
              val.append(str(column.get()))                                    #Add ion
            if j == 1:
              if float(column.get()) < 0:
                messagebox.showwarning('Invalid parameter.', 'Charge-state de'+
                                       'nsity must be a positive number. Ther'+
                                       'e will be a reset.')
                redo_form()
              else:
                val.append(float(column.get()))                                #Add csd
            if j == 2:
              val.append(column.val.get())                                     #Add fixed parameter for LM fit
            j += 1
        for p in power: csd['power'] = float(p.get())
      except:
        messagebox.showerror('Invalid parameter', 'One input is not a number.'+
                             ' Parameters will reset.')
        redo_form()
        return                                                                 #Solved a loop messagebox error
      csd_changed.set(True)

      for j in enumerate(rows):
        csd[val[3*j]] = val[3*j+1]
        fixed[val[3*j]] = val[3*j+2]
      ask_plot()

    def add_row():
      global i 
      i=i+1
      items = []
      for j in range(0,2):                                                     #Columns
        b = Entry(root, width=10)
        items.append(b)
        b.grid(row=i, column=j)
      var = IntVar()
      c = Checkbutton(root, variable = var, width=10, anchor='w')
      c.val = var
      items.append(c)
      c.grid(row = i, column = 2)

      rows.append(items)

    def delete_row():
      for rowno, row in reversed(list(enumerate(rows))):
        for i in row:
          i.destroy()
        rows.pop(rowno)
        break

    def redo_form():
      global csd, csd_standard
      csd = csd_standard
      for rowno, row in reversed(list(enumerate(rows))):
        for i in row:
          i.destroy()
        rows.pop(rowno)
      make_form(root)
      fetch()

    def make_form(root):
      global csd, i, fixed, fixed_standard, power, rows
      ions = list(csd.keys())[1:]                                              #To ignore power
      csds = list(csd.values())[1:]
      fixeds = list(fixed_standard.values())
      rows = []
      power = []
      for k in enumerate(ions):
        i=i+1
        items = []
        for j in range(0,2):                                                   #Columns
          b = Entry(root, width=10)
          if j == 0: b.insert(j, ions[k])
          else: b.insert(j, csds[k])
          items.append(b)
          b.grid(row=i, column=j)

        var = IntVar()
        var.set(fixeds[k])                                                     #To already have checked buttons
        c = Checkbutton(root, variable = var, width=10, anchor='w')
        c.val = var
        items.append(c)
        c.grid(row = i, column = 2)
        rows.append(items)
      e3 = Entry(root, width=5)
      e3.insert(0, csd.get('power'))
      e3.grid(row = 1, column = 2, padx = 5, pady = 5, sticky = 'e')
      power.append(e3)
      e4 = Label(root, text = 'Power 1e+', width=8)
      e4.grid(row = 1, column = 2, padx = 5, pady = 5, sticky = 'w')

    #CSD guess main
    #--------------------------------------------------------------------------
    #global csd_standard, i, ions, csds, fixed, rows, power
    root = Toplevel()
    root.title('Charge-state distribution')
    #root.iconbitmap(r'favico.ico')

    make_form(root)

    b1 = ttk.Button(root, text='Save', command=fetch)
    b1.grid(row = 0, column = 0)
    root.bind('<Return>', lambda events: fetch())                              #If presses enter, it saves ('presses' b1 button)
    b2 = ttk.Button(root, text='Reset', command=redo_form)
    b2.grid(row = 0, column = 1)
    b3 = ttk.Button(root, text='Quit', command=quit_window)
    b3.grid(row = 0, column = 2)
    bt = ttk.Button(root, text='Add row', command=add_row)
    bt.grid(row =1, column=0)
    dl = ttk.Button(root, text='Delete row', command=delete_row)
    dl.grid(row =1, column=1)

    e0 = Label(root, text = 'Ion')
    e0.grid(row = 2, column = 0)
    e1 = Label(root, text= 'CSD (ion/cm3)')
    e1.grid(row = 2, column = 1)
    e2 = Label(root, width=15, text='LM fixed', anchor='w')
    e2.grid(row = 2, column = 2)

    root.mainloop()

  def clear():
    spectrum.clear()
    spectrum.set_xlabel('Energy (eV)')
    spectrum.set_ylabel('Counts/s')
    ax_tf.clear()
    ax_tf.axis('off')
    f.canvas.draw()
    load_var.set('')
    tf_var.set('')
    plot_done.set(False)
    param_changed.set(False)
    csd_changed.set(False)
    global csd, fixed
    csd = csd_standard
    fixed = fixed_standard
    tf_changed.set(False)
    scale_var.set('Linear')
    check_norm.set('None')
    check_off.set('Set')                                                       #Known bug that leaves the box disabled if parameters window is open

  def _help():
    messagebox.showinfo('Under construction.','Help section is in the making.'+
                        " Questions may be answered in André Fernandes' MSc t"+
                        'hesis from 2019')

  buttons_frame = Frame(sim)                                                   #Transition buttons frame
  buttons_frame.grid(row=2,column=0,columnspan=1, sticky='w')
  buttons_frame2 = Frame(sim)                                                  #Choice buttons and calculate/clear buttons
  buttons_frame2.grid(row=0,column=2, sticky='ne')
  buttons_frame3 = Frame(sim)                                                  #Quit, open file and export buttons
  buttons_frame3.grid(row=1,column=2, sticky='se')
  buttons_frame4 = Frame(sim)                                                  #Unused frames
  buttons_frame4.grid(row=2,column=1,columnspan=2, sticky='ne')
  buttons_frame5 = Frame(sim)                                                  #Progress bar frame
  buttons_frame5.grid(row=3,column=0,columnspan=3, sticky='we')

  progress_var = DoubleVar()                                                   #Progress bar float

  #--------------------------------------------------------------------------
  #TRANSITION BUTTONS
  def check_all():                                                             #Checks every button if total, unchecks if not
    if total_var.get():
      if (cs_K_exc_exists and cs_K_ion_exists and cs_KL_ion_exists and
          cs_KLL_ion_exists):
        kexc_var.set(True)
        kion_var.set(True)
        klion_var.set(True)
        kllion_var.set(True)
      else:
        messagebox.showwarning('Warning: Cannot consider total', 'Missing cro'+
                               'ss section database for total consideration.')
        total_var.set(False)
    else:
      kexc_var.set(False)
      kion_var.set(False)
      klion_var.set(False)
      kllion_var.set(False)
    ask_plot()

  def check_total():                                                           #As soon something changes, it checks whether or not total should be checked
    if not (kexc_var.get() and kion_var.get() and klion_var.get() and
            kllion_var.get()):
      total_var.set(False)
    elif (kexc_var.get() and kion_var.get() and klion_var.get() and
          kllion_var.get()):
      total_var.set(True)

    if kexc_var.get() and not cs_K_exc_exists:                                 #This checks if it is possible to consider process
      messagebox.showwarning('Warning: Cannot consider process', 'Missing K s'+
                             'hell excitation cross section database.')
      kexc_var.set(False)
    if kion_var.get() and not cs_K_ion_exists:
      messagebox.showwarning('Warning: Cannot consider process', 'Missing K s'+
                             'hell ionization cross section database.')
      kion_var.set(False)
    if klion_var.get() and not cs_KL_ion_exists:
      messagebox.showwarning('Warning: Cannot consider process', 'Missing KL '+
                             'shells double ionization cross section database'+
                             '.')
      klion_var.set(False)
    if kllion_var.get() and not cs_KLL_ion_exists:
      messagebox.showwarning('Warning: Cannot consider process', 'Missing KLL'+
                             ' shells triple ionization cross section databas'+
                             'e.')
      kllion_var.set(False)
    ask_plot()

  def scale(stype):
    if load_var.get() != '' or plot_done.get():
      if stype == 'Log': spectrum.set_yscale('log')
      else: spectrum.set_yscale('linear')
      f.canvas.draw()

  M1 = Checkbutton(buttons_frame, text = u'K-shell exc.', variable = kexc_var,
                   command = check_total, onvalue = True, offvalue = False)
  M2 = Checkbutton(buttons_frame, text = u'K-shell ion.', variable = kion_var,
                   command = check_total, onvalue = True, offvalue = False)
  M3 = Checkbutton(buttons_frame, text = u'KL-shell ion.',variable = klion_var,
                   command = check_total, onvalue = True, offvalue = False)
  M4 = Checkbutton(buttons_frame, text = u'KLL-shell ion.',
                   variable = kllion_var, command = check_total,
                   onvalue = True, offvalue = False)
  M5 = Checkbutton(buttons_frame, text = 'Total', variable = total_var,
                   command = check_all, onvalue = True, offvalue = False)
  M6 = Checkbutton(buttons_frame, text = 'Metastable states',
                   variable = metastates_var, command = ask_plot,
                   onvalue = True, offvalue = False)
  M7 = Checkbutton(buttons_frame, text = 'Processes contributions',
                   variable = contributions_var, command = ask_plot,
                   onvalue = True, offvalue = False)
  M1.pack(side=LEFT)
  M2.pack(side=LEFT)
  M3.pack(side=LEFT)
  M4.pack(side=LEFT)
  M5.pack(side=LEFT)
  M6.pack(side=LEFT)
  M7.pack(side=LEFT)
  #----------------------------------------------------------------------------

  #----------------------------------------------------------------------------
  def write_to_xls():
    path='data.csv'
    print('Plasma diagnostics')
    with open(path, 'w', newline='') as csvfile:
      w1 = csv.writer(csvfile, delimiter=',',quotechar='|',
                      quoting=csv.QUOTE_MINIMAL)
    with open(path, 'a', newline='') as csvfile2:
      w2 = csv.writer(csvfile2, delimiter=',',quotechar='|',
                      quoting=csv.QUOTE_MINIMAL)
    return w1, w2
  #----------------------------------------------------------------------------

  ttk.Style().configure('red/black.TButton', foreground='red',
           background='black')                                                 #Calculate button
  ttk.Style().configure('parameters.TButton', foreground='black',
           background='black')                                                 #Parameters button
  ttk.Style().configure('clear.TButton', foreground='green',
           background='black')                                                 #Clear plot button
  ttk.Button(master=buttons_frame3, text='Quit', command=_quit
             ).pack(side=BOTTOM)
  ttk.Button(master=buttons_frame3, text='Help', command=_help
             ).pack(side=BOTTOM)
  ttk.Button(master=buttons_frame3, text='Export fit',
             command=write_to_xls()).pack(side=BOTTOM)                  #To do
  ttk.Button(master=buttons_frame3, text='Intensity dist.', command=load_tf
             ).pack(side=BOTTOM)
  ttk.Button(master=buttons_frame3, text='Load exp. file', command=load_data
             ).pack(side=BOTTOM)

  ttk.Button(master=buttons_frame2, text='Calculate',
             command=plotting(),style='red/black.TButton'
             ).pack(side=BOTTOM)
  ttk.Button(master=buttons_frame2, text='LM fit',
             command=LM_plotting(), style= 'red/black.TButton'
             ).pack(side=BOTTOM)

  ttk.Button(master=buttons_frame2, text='Clear plot', command = clear,
             style = 'clear.TButton').pack(side=BOTTOM)
  ttk.Button(master=buttons_frame2, text='Parameters', command =
             plot_parameters, style = 'parameters.TButton').pack(side=BOTTOM)
  ttk.Button(master=buttons_frame2, text='CSD', command = csd_guess,
             style = 'parameters.TButton').pack(side=BOTTOM)

  scale_var = StringVar(value='Linear')
  OptionMenu(buttons_frame2, scale_var, 'Linear', 'Log', command = scale
             ).pack(side=BOTTOM)

  choice_var = StringVar(value='Sim.')
  OptionMenu(buttons_frame2, choice_var, 'Sim.', 'Stick',
             command = lambda choice_var: ask_plot()).pack(side=BOTTOM)

  progressbar = ttk.Progressbar(buttons_frame5, variable=progress_var,
                                maximum=100)
  progressbar.pack(fill=X, expand=1)

  sim.mainloop()
# PlasmaFit
PlasmaFit was was developed during my Master's dissertation in Engineering Physics. The application was done specifically for highly charged ion plasma diagnostics for astrophysical and energy applications.

Plasma diagnostics are crucial for projects like the International Thermonuclear Experimental Reactor (ITER), the world’s largest tokamak, which is being built in the south of France. These diagnostics demand theoretical and experimental studies in order to understand the origin of spectral emissions observed in the plasma. From the balance between the creation and decay of excited states, one can infer on the ionic abundance within the plasma and hence on their quality. Thus, electron-impact ionization and excitation, which require cross section values for any creation process, need to be evaluated for a large number of states and for a wide energy range. Typically, the values are determined computationally with models such as the distorted wave Born approximation (DWBA) and, due to the simplicity of the approach and the large amount of atomic data needed for such codes, the modified relativistic binary encounter Bethe (MRBEB). With these, and the transition energies from the excited states, it is possible to determine the charge-state distribution within the plasma. With the ion structure information, we can determine, for example, the ion temperature and impurities from wall contamination in the plasma. The methodology is also relevant in the field of astrophysics, wherein theoretical calculations make it possible to know the characteristics of distant plasma bodies. This work presents a code which can load x-ray experimental spectra and experimental transfer functions for irregular x-ray detection. The user can input several parameters and charge-state densities for the ions in order to present a simulated spectrum. A [Levenberg-Marquardt algorithm](https://lmfit.github.io/lmfit-py/installation.html) was implemented in order to approximate the ion densities to the experimental data.

The link for the dissertation can be found [here](https://run.unl.pt/handle/10362/89460).

The final version inclues a periodic table and databases that contain elements other than Ar (Z=18), however, these were not done by me and therefore not in this repository.

# Images
Main window GUI:

![Main window.](/Images/main_window.png?raw=true)

Secondary windows with adjustable parameters:

<img src="/Images/sec_windows.png" alt="Secondary windows" width="750"/>

Simulated spectrum fitted over experimental data:

![Best fit spectrum](/Images/best_fit_spectrum.png?raw=true)

# PlasmaFit
PlasmaFit was was developed during my Master's dissertation in Engineering Physics. The application was done specifically for highly charged ion plasma diagnostics for astrophysical and energy applications.

Plasma diagnostics are crucial for projects like the International Thermonuclear Experimental Reactor (ITER), the world’s largest tokamak, which is being built in the south of France. These diagnostics demand theoretical and experimental studies in order to understand the origin of spectral emissions observed in the plasma. From the balance between the creation and decay of excited states, one can infer on the ionic abundance within the plasma and hence on their quality. Thus, electron-impact ionization and excitation, which require cross section values for any creation process, need to be evaluated for a large number of states and for a wide energy range. Typically, the values are determined computationally with models such as the distorted wave Born approximation (DWBA) and, due to the simplicity of the approach and the large amount of atomic data needed for such codes, the modified relativistic binary encounter Bethe (MRBEB). With these, and the transition energies from the excited states, it is possible to determine the charge-state distribution within the plasma. With the ion structure information, we can determine, for example, the ion temperature and impurities from wall contamination in the plasma. The methodology is also relevant in the field of astrophysics, wherein theoretical calculations make it possible to know the characteristics of distant plasma bodies. This work presents a code which can load x-ray experimental spectra and experimental transfer functions for irregular x-ray detection. The user can input several parameters and charge-state densities for the ions in order to present a simulated spectrum. A [Levenberg-Marquardt algorithm](https://lmfit.github.io/lmfit-py/installation.html) was implemented in order to approximate the ion densities to the experimental data.

The link for the dissertation can be found [here](https://run.unl.pt/handle/10362/89460).

The final version inclues a periodic table and databases that contain elements other than Ar (Z=18), however, these were not done by me and therefore not in this repository.

# GUI's
Main GUI from PlasmaFit.py, highlighting four sections with the possible inputs from the user (See Figure 3.5 from my dissertation).

![Main GUI from PlastmaFit.py.](/screenshots/main_window.png?raw=true)


Charge-state distribution GUI from PlasmaFit.py, showing the CSD the user can input. The “Add row” and “Delete row” buttons allow the user to insert or remove an ion. Any box checked in the “LM fixed” column fixes the parameter during the Levenberg-Marquardt calculation.

<img src="/screenshots/cds_window.png" alt="Parameters GUI's from PlasmaFit.py." width="350"/>


Parameters GUI’s from PlasmaFit.py, showing the parameters the user can input. The two windows on the right side are shown in case no plotting is done, allowing the user to manually insert the parameters for the background noise and normalization. With a plot shown, the user can simply click on the coordinates instead.

<img src="/screenshots/parameters_window.png" alt="Parameters GUI's from PlasmaFit.py." width="750"/>


Final simulated spectrum over the experimental data from the ion source in Paris (SIMPA). A slope to the background noise was added by using the “Two points” parameter. After clicking in two different coordinates, a linear equation is given as offset for the background noise.

![Final spectrum](/screenshots/best_fit_spectrum.png?raw=true)

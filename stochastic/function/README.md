# Functions

This folder contains all the python files needed to run the stochastic simulations used in the article "Stochastic dynamics at the back of gene drive propagation: drive eradication, wild-type recolonisation or chasing". 

`drive_simu_1D.py` simulates the stochastic model presented in the article, in one dimension. It corresponds to the global simulations with drive and wild-type individuals. The data are saved in a folder named `stoch_not_save` (which as its name suggests, is not stored in the Github folder because it is too heavy) and used in `main.py`. It also produces some snapshops of the traveling wave and a kymograph.

`drive_simu_2D.py`:exactly the same as `drive_simu_1D.py`, but in two spatial dimensions and it does not produce a kymograph.

`gw_space.py` contains all the functions to simulate the Galton-Watson process introduced in the article (with migration and an exponential initial condition). The dynamics of this process are given by the wild-type population considering that the proportion of drive is one (approximation of the condition at the back of the wild-type wave). These functions are used in `main.py`.

`L_speed_eigen_values.py` contains functions that provide some analytically values based on the exponential profile of the pulled wave.

`main.py` brings all together the data and the functions introduced before, and plot the figures of the article (all expect the snapshots of the traveling wave in 1D and 2D, the kymograph, and the illustrations made with inkscape). It uses the data computed with `drive_simu_1D.py`, and the functions from `gw_space.py` and `L_speed_eigen_values.py`.
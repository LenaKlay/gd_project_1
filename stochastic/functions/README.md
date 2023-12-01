# Functions

This folder contains all the code needed to run the stochastic simulations used in the soon submitted article "Stochastic dynamics at the back of gene drive propagation: drive eradication, wild-type recolonisation or chasing". 

`global_simu.py` simulates the stochastic model presented in the article. It corresponds to the global simulations with drive and wild-type individuals. The datas are saved in a folder named `stoch_not_save` which, as its name suggests, is not stored in the Github folder because it is too heavy.

`galton_watson_simu.py` contains the functions that simulate a Galton-Watson process in space with migration. The dynamics is given by the wild-type population considering that the proportion of drive is one (approximation of the condition at the back of the wild-type wave).

`L_speed_eigen_values.py` contains functions that provide some analytically values based on the exponential profile of the pulled wave.

`figure.py` plot the figures of the article. It uses the datas from `galton_watson_simu.py` and the functions in `galton_watson_simu.py` and `L_speed_eigen_values.py` to do so.
# energy-loss-simulation
How particles loss energy when their decay rate depend on their energy. 
A numerical simulation of the analytical work presented in http://inspirehep.net/record/1276809

energyloss.py contains the code of both the simulation and the visualizations. 
It consits of two different classes, namely, EnergyLoss (the central class) and 
gamma_gen(rv_continuous), that implements a probability distribution.

sc_energyloss.py is a script to run energyloss.py from the command line:  
$ python sc_energyloss.py [constant--fluctuating] n_particles 
where n_particles is an integer

visualizations.ipynb is a notebook with the main visualizations for both 
constant and fluctuaing energy loss per decay.

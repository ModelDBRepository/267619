This model of the Lateral IntraParietal area (LIP) and the Frontal Eye Fields (FEF) was described in the eLife article "Interacting rhythms enhance sensitivity of target detection in a fronto-parietal computational model of visual attention" by Amélie Aussel, Ian C. Fiebelkorn, Sabine Kastner, Nancy J. Kopell and Benjamin R. Pittman-Polletta.   (https://doi.org/10.7554/eLife.67684)


1) Requirements :
This simulator was developed using Python 3.9 with Brian2 version 2.4.2, under a Linux environment. Other important packages are cython (0.29.30) and joblib (1.1.0).
The complete specifcations used for the environment can be found in the file :
environment_LIPFEF.yml
A similar virtual environment can therefore be created by using the command : conda env
create -f environment_LIPFEF.yml


2) User interface : 
One file is provided to set up and start simulations easily : user_interface_single_simulation.py.
Executing this file opens a Tkinter window with different tabs enabling to choose
which brain region(s) to model as well as some parameters to be used in the simulation.
When a simulation is run, a folder is created with a name of the form
« results_date_time », which contains simulations results along with a « parameter.txt »
file with the parameter choices.


3) Parallel processing :
The file parallel_simulations.py can be used to setup multiple simulations and run them
in parallel (using the joblib library).
A user interface is not provided at the moment. Instead, the choice of parameters of the model should be entered manually in the .py file.
The values of the target presentation times and the inhibition time constants to be used should be entered as lists. The complete set of simulations will then be the cartesian product of all the list of parameter values provided. Other simulation parameters are to be entered as single values and will be fixed upon all simulations (details are available as comments next to each parameter name in the file).


4) Other files in the simulator :
The "model_files" folder contains the necessary files to run the model.
The "cells" subfolder contains the files defining the different neuron types used in the model. (Please note that SOM cells in the model are sometimes referred to as "SI" cells for "Slow Interneurons")
- FEF visuomotor_module.py : defines the FEF visuomotor module
- FEF_visual_module.py : defines the FEF visual module
- FEF_full.py : puts together the FEF visuomotor and FEF visual module defined above and adds decision cells.
- LIP_superficial_layer.py : defines the superficial layer of LIP
- LIP_beta1.py : puts together the superficial layer of LIP with deep IB cells. (together, these neuron populations produce the beta1 rhythm observed in LIP during the poor theta phase)
- LIP_full.py : puts together the LIP beta1 module defined above with input layer RS and FS cells as well as deep layer SOM cells to obtain the full LIP model.
- FEF_and_LIP_single_simulation.py : puts together the FEF and LIP modules defined above. 




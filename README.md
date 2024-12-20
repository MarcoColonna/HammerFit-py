Version 0.0.1

HammerFit-py is a python interface to access and read Hammer reweighted histograms and include the Hammer's degrees of freedom (Wilson Coefficients and Form Factor parameters) in the pythonic fitting interfaces (i.e. iMinuit) for Semileptonic studies.

The definition of the model is all contained in a json file where the paths to the files containing the histograms (.dat files for Hammer-reweighted histograms, .root files for non-Hammer-reweighted histograms) are defined together with all the properties the module requires to perform fits and plot.

In the next future the module will be fully documented.

TO INSTALL LOCALLY HammerFit.py do:

git clone https://github.com/MarcoColonna/HammerFit-py.git

cd HammerFit-py

git checkout v0.0.1

pip install .
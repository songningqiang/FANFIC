The main code is main.cc which outputs four columns in the output file: alpha_e, alpha_mu, alpha_tau, chi2. The output file can be specified by user. It can be compiled by

g++ main.cc prototype.cc -o main

Several inputs in main.cc need to be specified before running the code:

Line 9: mass ordering

Line 11: the octant of theta_23, this only makes a difference if DUNE or HYPERK is used

Line 18: list of experiments to be used, instructions can be found in the comments above

Line 19: output file name

prototype.h is the header file containing class declaration, and prototype.cc include the definition of functions.

oscillationparams is a class containing all oscillation parameters and their uncertainties. sigmaplus is the standard deviation above the best-fit, and sigmaminus is the value below. By default, normal ordering is assumed but inverted ordering can also be specified explicitly. Nufit 5.0 with SK included is set as default.

Additional experiments can be added as the derived class of oscillationexperiment by overwriting the virtual function setosc. Inside this function, the standard deviation is supposed to be set to the smallest number expected. In this way, the order of the experiments won't matter. However, if several experiments are specified each with a different t23-dcp chi2 file, only the last chi2 file will be used. HyperK has not been implemented very properly but it can still be used now to account for the chi2 file from HyperK.

flavorregion is a class to calculate the final flavor composition at the earth. This is done with the member function evolveflavor. New physics in neutrino propagation can be introduced by overwriting this class.

likelihood is a class to calculate the total chi2. If no chi2 file served, the member function chisq will be called assuming gaussian. If served, chisqfromdata will be called instead.

prior is a class to randomly generate initial flavor composition and oscillation parameters. If necessary, priors other than flat can be implemented later.

More comments can be found in the code.




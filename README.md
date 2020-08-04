prototype.cc is a prototype to set the oscillation parameters.

oscillationparams is a class containing all oscillation parameters and their uncertainties. *sigmaplus is the standard deviation above the best-fit, and *sigmaminus is the value below. By default, normal ordering is assumed but inverted ordering can also be specified explicitly. Nufit 5.0 with SK included is set as default.

Additional experiments can be added as the derived class of oscillationexperiment by overwriting the virtual function setosc. Inside this function, the standard deviation is supposed to be set to the smallest number expected. In this way, the order of the experiments won't matter.
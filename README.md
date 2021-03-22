# FlavorRegion: A Monte Carlo Code to Find the Composition of Astrophysical Neutrinos

Authors: Ningqiang Song, Shirley Li, Mauricio Bustamante, Aaron Vincent, Carlos A. Arguelles

Reference: arXiv: 2012.12893

## Introduction and getting started

This is a Monte Carlo code to find the flavor composition of astrophysical neutrinos at Earth from the flavor composition at the source with neutrino oscillations and new physics. The oscillation parameters are constrained by the combined analysis of NuFIT or future neutrino oscillation experiments, including JUNO, DUNE and HyperK. The Earth flavor compositions can be constrained by current IceCube data or future neutrino telescopes, including IceCube-Gen2, P-ONE, KM3NeT and TAMBO. We also consider beyond Standard Model scenarios, including non-unitarity neutrino mixing and neutrino invisible decay. Some example codes are given below.

**main.cc**: find the flavor compositions at Earth from certain source flavor compositions assuming the oscillation parameters constrained by NuFIT or future oscillation experiments. To compile and run it, simply 
```
make
```
and
```
./main
```
Prior to run, please specify the number of samples, NuFIT version, the oscillation experiments to use (can be empty), the mass ordering, the value of \delta_{CP} and the octant of \theta_{23}, the initial flavor composition and the name of the output file. See the comments in the code for details. To enable parallelization with openmp, please uncomment the line
```
#pragma omp parallel for 
```
The code will draw samples according to the priors on source flavor compositions and neutrino oscillation parameters., and then calculate the Earth flavor composition and the correpsonding \chi^2. The output file is in the format 
f_{e,\oplus}, f_{\mu,\oplus}, f_{\tau,\oplus}, \chi^2.

**nonunitarity.cc**: find the flavor compositions at Earth from certain source flavor compositions assuming nonunitarity  or future oscillation experiments. To compile and run it, simply
```
make nonunitarity
```
and
```
./nonunitarity
```

**infersourcecomposition.cc**: infer the posterior of f_e in the source flavor composition (f_e:1-f_e:0) with the oscillation parameters constrained by NuFIT or future oscillation experiments and the Earth flavor compositions constrained by IceCube or future neutrino telescopes

**infersourcefraction.cc**: infer the posterior of the fraction of different production meachanisms in the source (k_\pi, k_\mu, k_n) (pion decay, damped muon, neutrond decay) with the oscillation parameters constrained by NuFIT or future oscillation experiments and the Earth flavor compositions constrained by IceCube or future neutrino telescopes.

**infersourcefraction_kpikmu.cc**: infer the posterior of the fraction of different production meachanisms in the source (k_\pi, k_\mu, 0) (pion decay, damped muon, neutrond decay, k_\pi+k_\mu=1), with the oscillation parameters constrained by NuFIT or future oscillation experiments and the Earth flavor compositions constrained by IceCube or future neutrino telescopes.

**neutrinodecay.cc**: infer the posterior of the neutrino decay rate m/\tau by integrating neutrino sources at different redshifts and by assuming certain source flavor compositions with the oscillation parameters constrained by NuFIT or future oscillation experiments and the Earth flavor compositions constrained by IceCube or future neutrino telescopes.

**neutrinodecay_masseigenstates.cc**: find the flavor compositions at Earth from certain mass eigenstates at Earth because of neutrino decay f_{\alpha,\oplus} = \sum_i |U_{\alpha i}|^2*f_i.

**neutrinodecay_kpigaussian.cc: infer the posterior of the neutrino decay rate m/\tau by integrating neutrino sources at different redshifts and by assuming gaussian source flavor compositions with the oscillation parameters constrained by NuFIT or future oscillation experiments and the Earth flavor compositions constrained by IceCube or future neutrino telescopes.

The main code is main.cc which outputs four columns in the output file: alpha_e, alpha_mu, alpha_tau, chi2. The output file can be specified by user. It can be compiled by

## Main modules in prototype.cc:

oscillationparams is a class containing all oscillation parameters and their uncertainties. sigmaplus is the standard deviation above the best-fit, and sigmaminus is the value below. By default, normal ordering is assumed but inverted ordering can also be specified explicitly. Nufit 5.0 with SK included is set as default.

Additional experiments can be added as the derived class of oscillationexperiment by overwriting the virtual function setosc. Inside this function, the standard deviation is supposed to be set to the smallest number expected. In this way, the order of the experiments won't matter. However, if several experiments are specified each with a different t23-dcp chi2 file, only the last chi2 file will be used. HyperK has not been implemented very properly but it can still be used now to account for the chi2 file from HyperK.

flavorregion is a class to calculate the final flavor composition at the earth. This is done with the member function evolveflavor. New physics in neutrino propagation can be introduced by overwriting this class.

likelihood is a class to calculate the total chi2. If no chi2 file served, the member function chisq will be called assuming gaussian. If served, chisqfromdata will be called instead.

prior is a class to randomly generate initial flavor composition and oscillation parameters. If necessary, priors other than flat can be implemented later.

More comments can be found in the code.




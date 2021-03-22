# FlavorRegion: A Monte Carlo Code to Find the Composition of Astrophysical Neutrinos

Authors: Ningqiang Song, Shirley Li, Mauricio Bustamante, Aaron Vincent, Carlos A. Arguelles

Reference: arXiv: [2012.12893](https://arxiv.org/pdf/2012.12893.pdf)

## Introduction and getting started

This is a Monte Carlo code to find the flavor composition of astrophysical neutrinos at Earth from the flavor composition at the source with neutrino oscillations and new physics. The oscillation parameters are constrained by the combined analysis of NuFIT or future neutrino oscillation experiments, including JUNO, DUNE and HyperK. The Earth flavor compositions can be constrained by current IceCube data or future neutrino telescopes, including IceCube-Gen2, P-ONE, KM3NeT and TAMBO. We also consider beyond Standard Model scenarios, including non-unitarity neutrino mixing and neutrino invisible decay. Some example codes are given below.

**Prerequisite**:  boost library [https://www.boost.org/](https://www.boost.org/)

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
f_{e,\oplus}, f_{\mu,\oplus}, f_{\tau,\oplus}, \chi^2. The posteriors can be obtained with kernel density estimate.

**nonunitarity.cc**: find the flavor compositions at Earth from certain source flavor compositions assuming nonunitarity  or future oscillation experiments. To compile and run it, simply
```
make nonunitarity
```
and
```
./nonunitarity
```
Other codes can also be compiled in a similar way by referring to the Makefile.

**infersourcecomposition.cc**: infer the posterior of f_e in the source flavor composition (f_e:1-f_e:0) with the oscillation parameters constrained by NuFIT or future oscillation experiments and the Earth flavor compositions constrained by IceCube or future neutrino telescopes

**infersourcefraction.cc**: infer the posterior of the fraction of different production meachanisms in the source (k_\pi, k_\mu, k_n) (pion decay, damped muon, neutrond decay) with the oscillation parameters constrained by NuFIT or future oscillation experiments and the Earth flavor compositions constrained by IceCube or future neutrino telescopes.

**infersourcefraction_kpikmu.cc**: infer the posterior of the fraction of different production meachanisms in the source (k_\pi, k_\mu, 0) (pion decay, damped muon, neutrond decay, k_\pi+k_\mu=1), with the oscillation parameters constrained by NuFIT or future oscillation experiments and the Earth flavor compositions constrained by IceCube or future neutrino telescopes.

**neutrinodecay.cc**: infer the posterior of the neutrino decay rate m/\tau by integrating neutrino sources at different redshifts and by assuming certain source flavor compositions with the oscillation parameters constrained by NuFIT or future oscillation experiments and the Earth flavor compositions constrained by IceCube or future neutrino telescopes.

**neutrinodecay_masseigenstates.cc**: find the flavor compositions at Earth from certain mass eigenstates at Earth because of neutrino decay f_{\alpha,\oplus} = \sum_i |U_{\alpha i}|^2*f_i.

**neutrinodecay_kpigaussian.cc**: infer the posterior of the neutrino decay rate m/\tau by integrating neutrino sources at different redshifts and by assuming gaussian source flavor compositions with the oscillation parameters constrained by NuFIT or future oscillation experiments and the Earth flavor compositions constrained by IceCube or future neutrino telescopes.

The main code is main.cc which outputs four columns in the output file: alpha_e, alpha_mu, alpha_tau, chi2. The output file can be specified by user. It can be compiled by

## Main classes:

**oscillationparams**:  set oscillation parameters and their uncertainties. sigmaplus is the standard deviation above the best-fit, and sigmaminus is the value below.

**oscillationexperiment**: base class to set neutrino oscillation parameters from neutrino oscillation experiment.

**JUNO**: set oscillation parameters with JUNO experiment.

**DUNE**: set oscillation parameters with DUNE experiment.

**HYPERK**: set oscillation parameters with HYPERK experiment.

**NUFIT**: set \sin^2\theta_{23} and \delta_{CP} with NUFIT, this is to include their correlations.

**flavorregion**: evolve neutrino flavors from the source to Earth.

**nonunitflavorregion**: set the oscillation matrix in case it is non-unitary and evolve neutrino flavors from the source to Earth. The corresponding \chi^2 should be given in files.

**neutrino decay**: evolve neutrino flavors from the source to Earth assuming neutrino decay.

**prior**: set the prior of parameters. prior.randomInitialFlavor() return a random flavor composition which sums up to 1. This is equivalent to a symmetric Dirichlet distribution on a 3-simplex, with concentration parameters a = 1. prior.random2d() returns 2d symmetric Dirichlet-distributed numbers. prior.flatPrior returns uniformly distributed numbers.

**gaussianprior**: set gaussian prior.

**likelihood**: find the \chi^2 for oscillation parameters. Oscillation experiments must be specified before using this class.

**likelihood_ice**: find the \chi^2 for Earth flavor compositions.


## Results:

Plots of the previous runs can be found in the "results" folder. plot_arrays.pdf summarizes part of the results, and more figures can be found in the "figures" folder.

## Using the code:

Feel free to use, modify or distribute the code. If you use the code in your publication, please cite the paper 
[https://arxiv.org/pdf/2012.12893.pdf](https://arxiv.org/pdf/2012.12893.pdf).



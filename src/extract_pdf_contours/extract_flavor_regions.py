#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@nbi.ku.dk"


"""
extract_flavor_regions.py:
    Routines to find the confidence regions of flavor composition using pre-
    sampled data.  For plotting, see plot_flavor_regions.py.

Created: 2020/09/15
Last modified: 2021/04/20
"""


import numpy as np
# from numpy import *
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d
from scipy.integrate import dblquad
from scipy.optimize import root_scalar
import json
import os
from pylab import *
from matplotlib import *
import pandas as pd


def read_data_raw(fname):

    # Pandas dataframe
    # Each element is a 4D list: [fe, fmu, ftau, chi2]
    df = pd.read_csv(fname, delim_whitespace=True, dtype=float,
            header=None, comment='#')

    # We extract and return individual columns as numpy arrays
    data_raw = [df[df.columns[i]].to_numpy(dtype=float, copy=False) \
            for i in range(len(df.columns))]

    data = np.column_stack((data_raw[0],data_raw[1])) # We only need [fe, fmu]
    weights = np.exp(-0.5*data_raw[3]) # The weights are log(-chi^2/2) values

    return data, weights


def read_data_raw_old(fname):

    # Each element of data_raw is a 4D list: [fe, fmu, ftau, chi2]
    data_raw = np.genfromtxt(fname, max_rows=100000000)

    data = data_raw[:,0:2]  # We only need fe, fmu
    weights = data_raw[:,3] # The weights are the chi2 values
    weights = [exp(-chi2/2.0) for chi2 in weights]

    return data, weights


def save_pdf(fname, fe, fmu, pdf):
   
    data = \
    { 
        'f_e_Earth': list(fe),
		'f_mu_Earth': list(fmu),
		'pdf': list(pdf)
    }

    with open(fname, 'w') as file:
        json.dump(data, file, indent=4, ensure_ascii=False)
    file.close()

    return


def read_pdf(fname):

    with open(fname) as file:
        data = json.load(file)

    fe = data['f_e_Earth']
    fmu = data['f_mu_Earth']
    pdf = data['pdf']

    return fe, fmu, pdf


def save_contours(fname, contours, contour_labels):

    contours = np.array(contours)
    data = {}
    # data['best_fit'] = {}
    # data['best fit']['f_e_Earth'] = bf[0] 
    # data['best fit']['f_mu_Earth'] = bf[1]
    for i in range(len(contours)):
        contour_fe = contours[i][0]
        contour_fmu = contours[i][1]
        data[contour_labels[i]] = {}
        data[contour_labels[i]]['f_e_Earth'] = list(contour_fe)
        data[contour_labels[i]]['f_mu_Earth'] = list(contour_fmu)
   
    with open(fname, 'w') as file:
        json.dump(data, file, indent=4, ensure_ascii=False)
    file.close()

    return


def filenames_no_ext(path):

    fnames = [os.path.splitext(f)[0] for f in os.listdir(path) 
        if (os.path.isfile(os.path.join(path, f)) and 
        	os.path.splitext(f)[1] == '.txt')]
    
    return fnames


def build_2d_kde(data, weights=None, kernel='gaussian', bw=0.01, norm=np.inf, 
	grid_points=100):

    # Each element of data is a 2D list of coordinates: [fe, fmu]
    # Each element of weights is a scalar weight, chi2

    grid_data_1d = linspace(0.0, 1.0, grid_points)
    fe = list(grid_data_1d)
    fmu = list(grid_data_1d)
    grid_data = np.stack(np.meshgrid(fmu, fe), -1).reshape(-1, 2)
    grid_data[:, [0, 1]] = grid_data[:, [1, 0]] # Swap indices

    kde = FFTKDE(kernel=kernel, bw=bw, norm=norm)
    pdf = kde.fit(data, weights=weights).evaluate(grid_data)

    return fe, fmu, pdf


def build_2d_kde_from_raw_data_file(fname, kernel='gaussian', bw=0.01, 
	norm=np.inf, grid_points=100, epszero=1.e-15):

    data, weights = read_data_raw(fname)

    fe, fmu, pdf = build_2d_kde(data, weights=weights, kernel=kernel, bw=bw, 
    	norm=norm, grid_points=grid_points)

    # Turn all small numbers to zero
    pdf[np.abs(pdf) < epszero] = 1.e-200

    return fe, fmu, pdf


def build_pdf(data, weights=None, grid_points=100):

    # Each element of data is a 2D list of coordinates: [fe, fmu]
    # Each element of weights is a scalar weight, chi2
    
    fe_data = data[:,0]
    fmu_data = data[:,1]

    pdf, fe, fmu = np.histogram2d(fe_data, fmu_data, bins=grid_points, 
        range=[[0.0, 1.0],[0.0, 1.0]], density=True, weights=weights)

    # histogram2d returns fe and fmu as 1D arrays of size grid_points (these are
    # the bin edges), and pdf as 2D array of size 
    # (grid_points-1)*(grid_points-1).  We take the central value of each bin of
    # fe and fmu in order for the arrays of fe and fmu that are returned to be 
    # also of size (grid_points-1).

    fe_central = [0.5*(fe[i]+fe[i+1]) for i in range(len(fe)-1)]
    fmu_central = [0.5*(fmu[i]+fmu[i+1]) for i in range(len(fmu)-1)]

    # Return the flattened pdf (to save it in json)
    pdf = pdf.flatten()
    # print((1.0/grid_points)**2*sum(pdf)) # == 1 (normalized)

    return fe_central, fmu_central, pdf
       

def build_pdf_from_raw_data_file(fname, grid_points=100, use_weights=False):

    data, weights = read_data_raw(fname)

    if use_weights:
        fe, fmu, pdf = build_pdf(data, grid_points=grid_points, weights=weights)
    else:
        fe, fmu, pdf = build_pdf(data, grid_points=grid_points, weights=None)

    return fe, fmu, pdf


def build_pdf_interp(fe, fmu, pdf_2d, kind='linear'):

    # fe and fmu are 1D lists
    # pdf_2d is a 2D list with dimensions (len(fe), len(fmu))

    pdf_interp = interp2d(fe, fmu, pdf_2d, kind=kind, bounds_error=False)

    return pdf_interp


def build_pdf_interp_from_data_file(fname, kind='linear'):

    # Also return the minimum and maximum values of the PDF in the file; we
    # will need them later when extracting contours

    fe, fmu, pdf = read_pdf(fname)
    # pdf_min = min(pdf)
    # pdf_max = max(pdf)
    pdf_2d = np.array(pdf).reshape((len(fe), len(fmu)))
    pdf_interp = build_pdf_interp(fe, fmu, pdf_2d, kind=kind)

    return pdf_interp #, pdf_min, pdf_max


def extract_pdf_contour(fe, fmu, pdf_2d, pdf_value, save_two_contours=False):

    # Extract the iso-contour of pdf(fe, fmu) == pdf_value
    
    # To do this, we use matplotlib to plot the contour, and then we extract
    # the vertices.  This leaves all internal smoothening to matplotlib, so we
    # do not have to do it ourselves.  

    # fe and fmu are 1D lists
    # pdf is a 2D list with dimensions (len(fe), len(fmu))

    cs = plt.contour(fe, fmu, pdf_2d, [pdf_value])

    if not save_two_contours:

        # Find the biggest contour (spurious ones are small)
        # Use this to find flavor regions
        index_sel = 0
        num_points = 0    
        for i in range(len(cs.collections[0].get_paths())):
            if (len(cs.collections[0].get_paths()[i]) > num_points):
                index_sel = i
                num_points = len(cs.collections[0].get_paths()[i])

        points = cs.collections[0].get_paths()[index_sel]
        v = points.vertices
        fe_ct = v[:,0]  # fe values on the contour
        fmu_ct = v[:,1] # fmu values on the contour

    else:

        # Save the two biggest contours
        # Use this to find regions of kpi, kmu, kn
        fe_ct = []
        fmu_ct = []
        # print(len(cs.collections[0].get_paths()))
        # quit()
        pathlengths = [[i, len(cs.collections[0].get_paths()[i])] \
            for i in range(len(cs.collections[0].get_paths()))]
        # print(pathlengths)
        ordered_pathlengths = sorted(pathlengths, key=lambda x: x[1])
        # print(ordered_pathlengths)
        # quit()
        for i in [ordered_pathlengths[-1][0], ordered_pathlengths[-2][0]]:
            points = cs.collections[0].get_paths()[i]
            v = points.vertices
            fe_ct += list(v[:,0])  # fe values on the contour
            fmu_ct += list(v[:,1]) # fmu values on the contour
        fe_ct += [fe_ct[0]]
        fmu_ct += [fmu_ct[0]]

    return fe_ct, fmu_ct


def split_pdf_contour(fe_ct, fmu_ct, npoints_branch=100, kind='linear'):

    # To integrate pdf_interp within an iso-contour, we split the contour into
    # a lower branch and an upper branch of fmu, and then interpolate so that 
    # both branches are evaluated at the same values of fe

    points = list([list(x) for x in zip(fe_ct, fmu_ct)])
    points.sort() # Sort in ascending order of fe

    # Both branches have the same common starting point
    lower_branch = [points[0]]
    upper_branch = [points[0]]
    fmu_ct_pivot = points[0][1] # First value of fmu

    if (len(points) > 0):
        # Go from the minimum to the maximum value of fe in the contour.  If a 
        # point has an fmu value smaller than the fmu pivot, then it belongs to
        # the lower branch; otherwise, it belongs to the upper branch.
        for p in points[1:]:
            if (p[1] <= fmu_ct_pivot):
                lower_branch.append(p)
            else:
                upper_branch.append(p)

    # Now evaluate both branches at the same values of fe.  They already match
    # at the smallest value of fe; we make them match at the highest value of 
    # fe, too, so that the contour closes.

    fe_min = points[0][0]
    fe_max = max(lower_branch[-1][0], upper_branch[-1][0])
    fe_branch = np.linspace(fe_min, fe_max, npoints_branch)

    if (lower_branch[-1][0] > upper_branch[-1][0]):
        upper_branch.append(lower_branch[-1])
    elif (upper_branch[-1][0] > lower_branch[-1][0]):
        lower_branch.append(upper_branch[-1])

    lower_branch = np.array(lower_branch)
    upper_branch = np.array(upper_branch)
    fmu_lower_branch_interp = interp1d(lower_branch[:,0], lower_branch[:,1],
        kind=kind, bounds_error=False, assume_sorted=True)
    fmu_upper_branch_interp = interp1d(upper_branch[:,0], upper_branch[:,1],
        kind=kind, bounds_error=False, assume_sorted=True)

    fmu_lower_branch = [fmu_lower_branch_interp(x) for x in fe_branch]
    fmu_upper_branch = [fmu_upper_branch_interp(x) for x in fe_branch]

    return fe_branch, fmu_lower_branch, fmu_upper_branch


def integrate_pdf_within_contour(fe, fmu, pdf_2d, pdf_values, 
    npoints_branch=100, kind='linear', epsabs=1.e-4, epsrel=1.e-4,
    use_interpolation=False, save_two_contours=False):

    # Integrate pdf_interp within the iso-contour of 
    # pdf_interp(fe, fmu) == pdf_value

    # fe, fmu are 1D lists
    # pdf is a 1D list with size len(fe)*len(fmu)

    integrals = []

    if not use_interpolation:

        # Do not use pdf_interp, use the pre-computed pdf_2d directly

        pdf = pdf_2d.flatten()
        area_bin = (fmu[1]-fmu[0])*(fe[1]-fe[0])

        for pdf_value in pdf_values:
            # print('  pdf_value = '+str(pdf_value))
            integral = area_bin*sum(x for x in pdf if x >= pdf_value)
            integrals.append(integral)
            # print('    integral = '+str(integral))
 
    else:

        # Build pdf 2D interpolating function, pdf_interp(fe, fmu)
        pdf_interp = build_pdf_interp(fe, fmu, pdf_2d, kind=kind)

        # # Test normalization of pdf_interp --- it correctly outputs 1
        # result, error = dblquad( \
        #     lambda fmu, fe: pdf_interp(fe, fmu), 
        #     0.0, 
        #     1.0,
        #     lambda fe: 0.0,
        #     lambda fe: 1.0-fe,
        #     epsabs=epsabs, epsrel=epsrel)
        # print(result)
        # quit()

        for pdf_value in pdf_values:
    
            # print('  pdf_value = '+str(pdf_value))

            # Extract the coordinates of the iso-contour
            fe_ct, fmu_ct = extract_pdf_contour(fe, fmu, pdf_2d, pdf_value,
                save_two_contours=save_two_contours)

            # print('    len(fe_ct) = '+str(len(fe_ct))+', len(fmu_ct) = '
                # +str(len(fmu_ct)))

            # fig = plt.figure(figsize=[9,9])
            # ax = fig.add_subplot(1,1,1)
            # ax.plot(fe_ct, fmu_ct)
            # plt.show()

            # Split the iso-contour into a lower branch of fmu and an upper branch 
            # of fmu; then build interpolating functions out of them
            fe_branch, fmu_lower_branch, fmu_upper_branch \
               = split_pdf_contour(fe_ct, fmu_ct, npoints_branch=npoints_branch, 
                     kind=kind)
            fmu_lower_branch_interp = interp1d(fe_branch, fmu_lower_branch,
                kind=kind, bounds_error=False, assume_sorted=True)
            fmu_upper_branch_interp = interp1d(fe_branch, fmu_upper_branch,
                kind=kind, bounds_error=False, assume_sorted=True)
            fe_branch_min = fe_branch[0]
            fe_branch_max = fe_branch[-1]

            # print(fe_branch)
            # fig = plt.figure(figsize=[9,9])
            # ax = fig.add_subplot(1,1,1)
            # ax.plot(fe_branch, fmu_lower_branch)
            # ax.plot(fe_branch, fmu_upper_branch)
            # plt.show()

            # The parameter order in dblquad is counter-intuitive:
            # scipy.integrate.dblquad(func, a, b, gfun, hfun) returns
            # the double (definite) integral of func(y, x) from x = a..b and 
            # y = gfun(x)..hfun(x).
            result, error = dblquad( \
                lambda fmu, fe: pdf_interp(fe, fmu) if (fmu <= 1.0-fe) else 0.0, 
                fe_branch_min, 
                fe_branch_max,
                lambda fe: fmu_lower_branch_interp(fe),
                lambda fe: fmu_upper_branch_interp(fe),
                epsabs=epsabs, epsrel=epsrel)

            # print('    integral = '+str(result))

            integrals.append(result)

    return integrals


def find_pdf_value_for_given_cl(fe, fmu, pdf_2d, cls, num_pdf_test_values=100,
    npoints_branch=100, kind='linear', epsabs=1.e-4, epsrel=1.e-4,
    root_xtol=1.e-4, root_rtol=1.e-4, use_interpolation=False,
    save_two_contours=False):

    # For a given value of integral (e.g., 0.90, 0.95), find pdf_value such that
    # integral of pdf_interp within the iso-contour pdf_interp(fe, fmu) == 
    # pdf_value equals cl (credibility level)

    # Do not compute the iso-contour for pdf_max, since none will be found
    # pdf_min = np.min(pdf_2d)
    # pdf_max = np.max(pdf_2d)
    # pdf_values = np.linspace(pdf_min, pdf_max, num_pdf_test_values)[::-1][1:]

    pdf = pdf_2d.flatten()
    pdf[np.abs(pdf) == 0.0] = 1.e-200
    log10_pdf = log10(pdf)
    log10_pdf = [x for x in log10_pdf if x > -200.0]
    log10_pdf_min = np.min(log10_pdf)
    log10_pdf_max = np.max(log10_pdf)
    pdf_min = 10.**log10_pdf_min
    pdf_max = 10.**log10_pdf_max
    # print(pdf_min, pdf_max)
    log10_pdf_values = np.linspace(log10_pdf_min, log10_pdf_max, 
        num_pdf_test_values)[::-1][1:]
    pdf_values = [10.**x for x in log10_pdf_values]

    # First compute the integral of the pdf within iso-contours defined by 
    # pdf_interp(fe,fmu) = pdf_value, varying pdf_value between pdf_min and
    # pdf_max
    integrals = integrate_pdf_within_contour(fe, fmu, pdf_2d, pdf_values, 
        npoints_branch=npoints_branch, kind=kind, epsabs=epsabs, epsrel=epsrel,
        use_interpolation=use_interpolation, 
        save_two_contours=save_two_contours)
    # print(integrals)

    # Now build an interpolating function of the integral as a function of 
    # pdf_value, so that we can find for what values of pdf_value the integral
    # equals the requested cls
    integral_interp = interp1d(pdf_values, integrals, 
        kind=kind, bounds_error=False, assume_sorted=False)
    # print(pdf_values)
    # print(integral_interp(30.))

    roots = []
    for cl in cls:
        print('  cl = '+str(cl))
        pdf_lo = pdf_values[0]
        pdf_hi = pdf_values[-1]
        for i in range(len(integrals)-1):
            if (cl >= integrals[i]) and (cl < integrals[i+1]):
                pdf_lo = pdf_values[i+1]
                pdf_hi = pdf_values[i]
                break
        # print(integral_interp(pdf_lo)-cl, integral_interp(pdf_hi)-cl)
        sol = root_scalar(lambda pdf_value: integral_interp(pdf_value)-cl, 
            bracket=[pdf_lo, pdf_hi], 
            # bracket=[pdf_min, pdf_max], 
            # x0=0.5*(pdf_max-pdf_min),
            # x1=0.6*(pdf_max-pdf_min),
            xtol=root_xtol, rtol=root_rtol)
        roots.append(sol.root)

        # print(pdf_lo, pdf_hi)
        # print(sol.root)

    return roots


def find_cl_contours(fe, fmu, pdf_2d, cls, num_pdf_test_values=100,
    npoints_branch=100, kind='linear', epsabs=1.e-4, epsrel=1.e-4,
    root_xtol=1.e-4, root_rtol=1.e-4, use_interpolation=False,
    save_two_contours=False):

    # Values of pdf at which to extract the iso-contours of 
    # pdf_interp(fe,fmu) = pdf_value such that the integrals of pdf_interp
    # within these iso-contours equal the values in cls
    pdf_values = find_pdf_value_for_given_cl(fe, fmu, pdf_2d, cls, 
        num_pdf_test_values=num_pdf_test_values, npoints_branch=npoints_branch, 
        kind=kind, epsabs=epsabs, epsrel=epsrel, root_xtol=root_xtol, 
        root_rtol=root_rtol, use_interpolation=use_interpolation,
        save_two_contours=save_two_contours)
    # print(pdf_values)

    # Now extract the contours.  Each element in contours contains one contour,
    # as a list of points [fe, fmu].
    contours = [extract_pdf_contour(fe, fmu, pdf_2d, pdf_value, 
        save_two_contours=save_two_contours) 
        for pdf_value in pdf_values]

    return contours


def find_cl_contours_from_pdf_file(fname, cls, num_pdf_test_values=100,
    npoints_branch=100, kind='linear', epsabs=1.e-4, epsrel=1.e-4,
    root_xtol=1.e-4, root_rtol=1.e-4, use_interpolation=False,
    save_two_contours=False):

    fe, fmu, pdf = read_pdf(fname)
    pdf_2d = np.array(pdf).reshape((len(fe), len(fmu)))

    contours = find_cl_contours(fe, fmu, pdf_2d, cls, 
        num_pdf_test_values=num_pdf_test_values,
        npoints_branch=npoints_branch, kind=kind, 
        epsabs=epsabs, epsrel=epsrel,
        root_xtol=root_xtol, root_rtol=root_rtol,
        use_interpolation=use_interpolation,
        save_two_contours=save_two_contours)

    return contours


"""
**Note: Before running the code below, generate the files of flavor ratios vs.
chi2 (or fractions of subcomponents vs. chi2)**

Run the code below to generate the PDFs for each pre-computed file of flavor 
ratios vs. chi2.  Only need to run this once.

Comment the cases that you do not wish to generate.  
"""

def main():

    # flavor/std/var_delta/main/all_fS: 
    # Flavor regions, standard mixing, varying delta_CP, NuFit 4.0, varying all f_S
    PATH_CHI2 = '../results/chi2/flavor/std/var_delta/main/all_fS/'
    PATH_PDF = '../results/pdf/flavor/std/var_delta/main/all_fS/'
    PATH_CONTOURS = '../results/contours/flavor/std/var_delta/main/all_fS/'

    # flavor/std/var_delta/main/fixed_fS: 
    # Flavor regions, standard mixing, varying delta_CP, NuFit 4.0, fixed f_S
    PATH_CHI2 = '../results/chi2/flavor/std/var_delta/main/fixed_fS/'
    PATH_PDF = '../results/pdf/flavor/std/var_delta/main/fixed_fS/'
    PATH_CONTOURS = '../results/contours/flavor/std/var_delta/main/fixed_fS/'

    # flavor/std/fixed_delta/main/all_fS: 
    # Flavor regions, standard mixing, fixed delta_CP, NuFit 4.0, varying all f_S
    PATH_CHI2 = '../results/chi2/flavor/std/fixed_delta/main/all_fS/'
    PATH_PDF = '../results/pdf/flavor/std/fixed_delta/main/all_fS/'
    PATH_CONTOURS = '../results/contours/flavor/std/fixed_delta/main/all_fS/'

    # flavor/std/fixed_delta/main/fixed_fS: 
    # Flavor regions, standard mixing, fixed delta_CP, NuFit 4.0, fixed f_S
    PATH_CHI2 = '../results/chi2/flavor/std/fixed_delta/main/all_fS/'
    PATH_PDF = '../results/pdf/flavor/std/fixed_delta/main/all_fS/'
    PATH_CONTOURS = '../results/contours/flavor/std/fixed_delta/main/all_fS/'

    # flavor/std/var_delta/nufit_versions/all_fS: 
    # Flavor regions, standard mixing, varying delta_CP, previous NuFit versions, varying all f_S
    PATH_CHI2 = '../results/chi2/flavor/std/var_delta/nufit_versions/all_fS/'
    PATH_PDF = '../results/pdf/flavor/std/var_delta/nufit_versions/all_fS/'
    PATH_CONTOURS = '../results/contours/flavor/std/var_delta/nufit_versions/all_fS/'

    # flavor/std/var_delta/nufit_versions/fixed_fS: 
    # Flavor regions, standard mixing, varying delta_CP, previous NuFit versions, fixed f_S
    PATH_CHI2 = '../results/chi2/flavor/std/var_delta/nufit_versions/fixed_fS/'
    PATH_PDF = '../results/pdf/flavor/std/var_delta/nufit_versions/fixed_fS/'
    PATH_CONTOURS = '../results/contours/flavor/std/var_delta/nufit_versions/fixed_fS/'

    # flavor/decay/var_delta/main/all_fS: 
    # Flavor regions, neutrino decay, varying delta_CP, NuFit 4.0, varying all f_S
    PATH_CHI2 = '../results/chi2/flavor/decay/var_delta/main/all_fS/'
    PATH_PDF = '../results/pdf/flavor/decay/var_delta/main/all_fS/'
    PATH_CONTOURS = '../results/contours/flavor/decay/var_delta/main/all_fS/'

    # flavor/decay/var_delta/main/fixed_fS: 
    # Flavor regions, neutrino decay, varying delta_CP, NuFit 4.0, fixed f_S
    PATH_CHI2 = '../results/chi2/flavor/decay/var_delta/main/fixed_fS/'
    PATH_PDF = '../results/pdf/flavor/decay/var_delta/main/fixed_fS/'
    PATH_CONTOURS = '../results/contours/flavor/decay/var_delta/main/fixed_fS/'

    # flavor/nonunitary/var_delta/main/all_fS: 
    # Flavor regions, non-unitarity, varying delta_CP, NuFit 4.0, varying all f_S
    PATH_CHI2 = '../results/chi2/flavor/nonunitary/var_delta/main/all_fS/'
    PATH_PDF = '../results/pdf/flavor/nonunitary/var_delta/main/all_fS/'
    PATH_CONTOURS = '../results/contours/flavor/nonunitary/var_delta/main/all_fS/'

    # flavor/nonunitary/var_delta/main/fixed_fS: 
    # Flavor regions, non-unitarity, varying delta_CP, NuFit 4.0, fixed f_S
    PATH_CHI2 = '../results/chi2/flavor/nonunitary/var_delta/main/fixed_fS/'
    PATH_PDF = '../results/pdf/flavor/nonunitary/var_delta/main/fixed_fS/'
    PATH_CONTOURS = '../results/contours/flavor/nonunitary/var_delta/main/fixed_fS/'

    # subcomponents: Fractions of sources with different production mechanism
    PATH_CHI2 = os.getcwd()+'../results/chi2/subcomponents/2D/'
    PATH_PDF = os.getcwd()+'../results/pdf/subcomponents/2D/'
    PATH_CONTOURS = os.getcwd()+'../results/contours/subcomponents/2D/'

    fnames = filenames_no_ext(PATH_CHI2) # Keep the same file names throughout

    print("Generating the PDFs from pre-computed chi2 files... ")
    if not os.path.exists(PATH_PDF):
        os.mkdir(PATH_PDF)
    for fname in fnames:
        print(fname+'_pdf.json... ', end='')
        fe, fmu, pdf = build_pdf_from_raw_data_file(PATH_CHI2+fname+'.txt', 
            grid_points=250, use_weights=True)
        save_pdf(PATH_PDF+fname+'_pdf.json', fe, fmu, pdf)
        print('Done')

    print("Generating the flavor contours from pre-computed PDF files... ")
    if not os.path.exists(PATH_CONTOURS):
        os.mkdir(PATH_CONTOURS)
    cls = [0.68, 0.90, 0.95, 0.997] # Credible levels of the contours to compute
    cls_labels = ['0.68', '0.90', '0.95', '0.997'] # Labels to save in the json file
    for fname in fnames:
        fname_in = PATH_PDF+fname+'_pdf.json'
        fname_out = PATH_CONTOURS+fname+'_contours.json' 
        print(fname+'... ')
        contours = find_cl_contours_from_pdf_file(fname_in, cls, 
            num_pdf_test_values=200, 
            npoints_branch=200, 
            kind='linear', 
            epsabs=1.e-4, epsrel=1.e-4,
            root_xtol=1.e-6, root_rtol=1.e-6,
            save_two_contours=True)
        save_contours(fname_out, contours, cls_labels)

    return


if __name__ == '__main__':
    main()
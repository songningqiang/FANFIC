#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@nbi.ku.dk"


"""
extract_1d_subcomponents.py:
    Routines to find the confidence regions of f_pi, f_mu, f_nu using pre-
    sampled data.  For plotting, see plot_1d_subcomponents.py.

Created: 2020/12/10
Last modified: 2021/04/20
"""

import numpy as np
# from numpy import *
from scipy.interpolate import interp1d
# from scipy.interpolate import interp2d
# from scipy.integrate import dblquad
from scipy.integrate import quad
from scipy.optimize import root_scalar
import json
import os
from pylab import *
from matplotlib import *
from scipy.signal import savgol_filter

from root_solver import *


def read_data_raw(fname):

    # Each element of data_raw is a 3D list: [fpi, fmu, chi2]
    data_raw = np.genfromtxt(fname)

    data = data_raw[:,0]  # We only need fpi (because fmu = 1-fpi)
    weights = data_raw[:,2] # The weights are the chi2 values

    return data, weights


def save_pdf(fname, fe, pdf):
   
    data = \
    { 
        'f_pi': list(fe),
		'pdf': list(pdf)
    }

    with open(fname, 'w') as file:
        json.dump(data, file, indent=4, ensure_ascii=False)
    file.close()

    return


def read_pdf(fname):

    with open(fname) as file:
        data = json.load(file)

    fe = data['f_pi']
    pdf = data['pdf']

    return fe, pdf
    # return fe, fmu, pdf


def filenames_no_ext(path, ext='txt'):

    fnames = [os.path.splitext(f)[0] for f in os.listdir(path) 
        if (os.path.isfile(os.path.join(path, f)) and 
            os.path.splitext(f)[1] == '.'+ext)]
    
    return fnames


def build_pdf(data, weights=None, grid_points=100):

    # Each element of data is a scalar, fpi
    # Each element of weights is a scalar weight, chi2
    
    fpi_data = data

    pdf, fpi = np.histogram(fpi_data, bins=grid_points, 
        range=[0.0, 1.0], density=True, weights=weights)

    # histogram returns fe as 1D array of size grid_points (these are
    # the bin edges), and pdf as 1D array of size 
    # (grid_points-1).  We take the central value of each bin of
    # fpi in order for the array of fe that is returned to be 
    # also of size (grid_points-1).

    fpi_central = [0.5*(fpi[i]+fpi[i+1]) for i in range(len(fpi)-1)]

    # Return the flattened pdf (to save it in json)
    pdf = pdf.flatten()
    # print((1.0/grid_points)**2*sum(pdf)) # == 1 (normalized)

    return fpi_central, pdf
       

def build_pdf_from_raw_data_file(fname, grid_points=100, use_weights=False):

    data, weights = read_data_raw(fname)
    weights = [exp(-chi2/2.0) for chi2 in weights]

    if use_weights:
        fpi, pdf = build_pdf(data, grid_points=grid_points, weights=weights)
    else:
        fpi, pdf = build_pdf(data, grid_points=grid_points, weights=None)

    return fpi, pdf


def build_pdf_interp(fpi, pdf, kind='linear', fill_value=None):

    # fpi and pdf_1d are 1D lists with the same length

    pdf_interp = interp1d(fpi, pdf, kind=kind, bounds_error=False,
        fill_value=fill_value)

    return pdf_interp


def build_pdf_interp_from_data_file(fname, kind='linear', fill_value=None):

    # Also return the minimum and maximum values of the PDF in the file; we
    # will need them later when extracting contours

    fpi, pdf = read_pdf(fname)
    pdf = savgol_filter(pdf, 21, 3)
    pdf_interp = build_pdf_interp(fpi, pdf, kind=kind, fill_value=fill_value)

    return pdf_interp #, pdf_min, pdf_max


def find_credible_interval(fname, lst_cl=[0.68], f_norm=1.0):

    def Find_PDF_Intercepts(interp_pdf, pdf_value, pdf_min, pdf_max, \
        f_eS_min, f_eS_max, eps=1.e-3):

        # The requested pdf_value does not intersect the pdf
        if ((pdf_value < pdf_min) or (pdf_value > pdf_max)):
            return -1

        if (pdf_value == pdf_min):
            return [f_eS_min]

        if (pdf_value == pdf_max):
            return [f_eS_max]

        # Find all the roots
        lst_roots = roots(lambda x: interp_pdf(x)-pdf_value, 0.0, 1.0,
                    eps=eps) 

        return lst_roots


    def Area_PDF(interp_pdf, pdf_value, pdf_min, pdf_max, f_eS_min, f_eS_max,
        eps=1.e-3):

        lst_roots = Find_PDF_Intercepts(interp_pdf, pdf_value, \
                    pdf_min, pdf_max, f_eS_min, f_eS_max)

        if (lst_roots != -1):

            num_roots = len(lst_roots)

            if (num_roots == 1):

                if (f_eS_min > f_eS_max):
                    f_eS_low = 0.0
                    f_eS_high = lst_roots[0]
                elif (f_eS_min < f_eS_max):
                    f_eS_low = lst_roots[0]
                    f_eS_high = 1.0
                area_pdf = quad(lambda f_eS: interp_pdf(f_eS), \
                            f_eS_low, f_eS_high, epsabs=eps, epsrel=eps)[0]

            else:

                # If there is an odd number of roots, we need to add an extra
                # fake root
                if (num_roots%2 == 1):
                    if (f_eS_min > f_eS_max):
                        lst_roots = [0.0] + lst_roots
                    elif (f_eS_min < f_eS_max):
                        lst_roots = lst_roots + [1.0]

                # Now the number of roots is even, calculate the area between
                # every pair of roots
                area_pdf = 0.0
                for i in range(0, len(lst_roots)-1, 2):
                    f_eS_low = lst_roots[i]
                    f_eS_high = lst_roots[i+1]
                    area_pdf += quad(lambda f_eS: interp_pdf(f_eS), \
                            f_eS_low, f_eS_high, epsabs=eps, epsrel=eps)[0]
        else:

            area_pdf = 0.0

        return area_pdf

    # Read in the pre-computed pdf of f_eS
    fe, pdf = read_pdf(fname)
    pdf_min = min(pdf)
    pdf_max = max(pdf)
    fe_pdf_max = fe[pdf.index(pdf_max)] # Value where the pdf is maximum

    # Interpolating function built from pre-compute pdf data
    interp_pdf = build_pdf_interp_from_data_file(fname, kind='linear',
        fill_value='0.0')

    """
    # Test the normalization
    norm = quad(lambda f_eS: interp_pdf(f_eS), 0.0, 1.0,
        epsabs=1.e-6, epsrel=1.e-6)[0]
    print(norm)
    """

    pdf_min = min(pdf)
    pdf_max = max(pdf)
    lst_merge = [[fe[i], pdf[i]] for i in range(len(fe))]
    lst_merge = sorted(lst_merge, key=lambda x: x[1])
    f_eS_max = lst_merge[-1][0] 
    f_eS_min = lst_merge[0][0]

    # Now find the interval of f_eS for the requested confidence levels
    lst_results = []
    for cl in lst_cl:
        # print(cl)
        pdf_value_requested = roots( \
            lambda pdf_value: Area_PDF(interp_pdf, pdf_value, 
                pdf_min, 
                pdf_max, 
                f_eS_min, 
                f_eS_max, eps=1.e-3)-cl, \
            pdf_min,
            pdf_max,
            eps=1.e-3)[0] 
        # print(pdf_value_requested)
        lst_roots_requested = \
            Find_PDF_Intercepts(interp_pdf, pdf_value_requested, 
                pdf_min, 
                pdf_max, 
                f_eS_min, 
                f_eS_max, 
                eps=1.e-3)
        # print(lst_roots_requested)
        # print(Area_PDF(interp_pdf, pdf_value_requested, pdf_min, pdf_max, \
        #     f_eS_min, f_eS_max))
        # print()
        if (len(lst_roots_requested) == 1):
            if (f_eS_min > f_eS_max):
                lst_roots_requested = [0.0] + lst_roots_requested
            elif (f_eS_min < f_eS_max):
                lst_roots_requested = lst_roots_requested + [1.0]
        lst_result = [cl, fe_pdf_max, lst_roots_requested]
        lst_results.append(lst_result)

    return lst_results


"""
**Note: Before running the code below, generate the files of the subcomponent
fraction kpi vs. chi2.**

Run the code below to generate the PDFs for each pre-computed file of kpi vs. 
chi2.  Only need to run this once.
"""

def main():

    PATH_CHI2 = '../results/chi2/subcomponents/1D/'
    PATH_PDF = '../results/pdf/subcomponents/1D/'

    if not os.path.exists(PATH_PDF):
        os.mkdir(PATH_PDF)
    fnames = filenames_no_ext(PATH_CHI2) # Keep the same file names throughout

    print("Generating the PDFs of kpi from pre-computed chi2 files... ")
    for fname in fnames:
        print(fname+'_pdf.json... ', end='')
        fe, pdf = build_pdf_from_raw_data_file(PATH_CHI2+fname+'.txt', 
            grid_points=300, use_weights=True)
        save_pdf(PATH_PDF+fname+'_pdf.json', fe, pdf)
        print('Done')

    print("Finding the credible intervals of kpi... ")
    for fname in fnames:
        print(fname)
        for cl in [0.68, 0.95, 0.997]:
            interval = find_credible_interval(PATH_PDF+fname+'.json', 
                lst_cl=[cl], f_norm=1.0)
            print(interval)

    return


if __name__ == '__main__':
    main()

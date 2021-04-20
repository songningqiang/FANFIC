#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@nbi.ku.dk"


"""
plot_1d_subcomponents.py:
    Contains routines to plot the f_pi PDF

Created: 2020/12/10
Last modified: 2021/04/20
"""


import sys
sys.path.insert(0, "../extract_pdf_contours")

import os
import gc

from numpy import *
import numpy as np

from pylab import *
from matplotlib import *
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.signal import savgol_filter

import json

from extract_1d_subcomponents import filenames_no_ext


def Plot_PDF_fpi(fname_in, fname_out, label, color, zorder, 
    output_format='pdf'):

    # Open the plot and format it
    mpl.rcParams['xtick.labelsize']=25
    mpl.rcParams['ytick.labelsize']=25
    mpl.rcParams['legend.fontsize']=16
    mpl.rcParams['legend.borderpad']=0.8#0.4
    mpl.rcParams['legend.borderaxespad']=1.0
    mpl.rcParams['legend.handlelength']=0.7#1.6
    mpl.rcParams['axes.labelpad']=10
    mpl.rcParams['ps.fonttype']=42
    mpl.rcParams['pdf.fonttype']=42
    mpl.rcParams['xtick.top']=True
    mpl.rcParams['xtick.bottom']=True
    mpl.rcParams['ytick.left']=True
    mpl.rcParams['ytick.right']=True
    mpl.rcParams['xtick.major.size']=10
    mpl.rcParams['xtick.minor.size']=5
    mpl.rcParams['ytick.major.size']=10
    mpl.rcParams['ytick.minor.size']=5
    mpl.rcParams['xtick.direction']='in'
    mpl.rcParams['ytick.direction']='in'
    mpl.rcParams['savefig.dpi']=300
    mpl.rcParams['savefig.bbox']='tight'
    mpl.rcParams['axes.labelsize']=25

    fig = plt.figure(figsize=[9,9])
    ax = fig.add_subplot(1,1,1)

    x_min = 0.0
    x_max = 1.0
    y_min = 0.0
    y_max = 17.0

    # First plot the shaded regions
    for file_index in range(len(fname_in)):

        if not (file_index < 5): continue

        with open(fname_in[file_index]+'.json') as file:
            data = json.load(file)
        fe = data['f_pi']
        pdf = data['pdf']
        # print(sum(pdf)*(fe[1]-fe[0])) # == 1
        if file_index != 10:
            pdf = savgol_filter(pdf, 21, 3)

        ax.fill_between(fe, pdf, y2=0, ec=color[file_index], 
            fc=color[file_index], zorder=zorder[file_index])

    # Now plot the borders
    for file_index in range(len(fname_in)):

        if not (file_index < 5): continue

        with open(fname_in[file_index]+'.json') as file:
            data = json.load(file)
        fe = data['f_pi']
        pdf = data['pdf']
        # print(sum(pdf)*(fe[1]-fe[0])) # == 1
        if file_index != 10:
            pdf = savgol_filter(pdf, 21, 3)

        ax.plot(fe, pdf, c=color[file_index], lw=2.0, label=label[file_index],
            zorder=max(zorder))

    ax.set_xlabel(r'Pion-decay fraction, $k_\pi$',
        fontsize=25)
    ax.set_ylabel(r'Probability density $\mathcal{P}(k_\pi)$',
        fontsize=25)

    ax.set_xlim([x_min, x_max])
    ax.set_ylim([y_min, y_max])

    ax.xaxis.set_major_locator(MultipleLocator(0.1))
    ax.xaxis.set_minor_locator(MultipleLocator(0.02))
    ax.yaxis.set_major_locator(MultipleLocator(2.0))
    ax.yaxis.set_minor_locator(MultipleLocator(0.5))
    ax.yaxis.set_major_locator(ticker.MaxNLocator(prune='lower'))

    leg = ax.legend(loc='upper left', ncol=1)
    for line in leg.get_lines():
        line.set_linewidth(15.0)

    pylab.savefig(fname_out+'.'+output_format, bbox_inches='tight', dpi=300)


def main():
    
    PATH_FIG = '../results/figures/subcomponents/1D/'
    fname_out = PATH_FIG+'fpi_recon_1D_pdf_NO_all'
    if not os.path.exists(PATH_FIG):
         os.mkdir(PATH_FIG)

    PATH_PDF = '../results/pdf/subcomponents/1D/'
    fname_in = [
        PATH_PDF+'output_NUFIT5.0_NO_IceCube2015_1D_pdf',
        PATH_PDF+'output_2020_NO_1D_pdf',
        # PATH_PDF+'output_2028_NO_1D_pdf', 
        PATH_PDF+'output_2040_NO_1D_pdf',
        PATH_PDF+'output_2040_NO_comb_1D_pdf'
    ]

    label = [
        r'2020 (measured):'+'\n'+r'IC (ApJ 1, 98) $\otimes$ NuFit 5.0'
            +'\n'+r'$k_\pi = 0.37_{-0.37}^{+0.59}~([0.00, 1.00])$',
        r'2020 (projected):'+'\n'+r' IC 8 yr $\otimes$ NuFit 5.0'
            +'\n'+r'$k_\pi = 0.88_{-0.31}^{+0.12}~([0.00, 1.00])$',
        # r'2028 (projected):'+'\n'+r' IC 15 yr $\otimes$ (NuFit 5.0+JUNO)',
        r'2040 (projected):'+'\n'+r' (IC 15 yr+IC-Gen2 10 yr) $\otimes$'+'\n'+r'(NuFit 5.0+JUNO+DUNE+HK)'
            +'\n'+r'$k_\pi = 0.99_{-0.06}^{+0.00}~([0.73, 1.00])$',
        r'2040 (projected):'+'\n'+r' (Combined $\nu$ telescopes) $\otimes$'+'\n'+r'(NuFit 5.0+JUNO+DUNE+HK)'
            +'\n'+r'$k_\pi = 0.99_{-0.05}^{+0.00}~([0.81, 1.00])$',
    ]

    color = [
        '#adcbe3',
        '#63ace5',
        # '#4b86b4',
        '#2a4d69',
        '#FFC914'
    ]

    zorder = [
        1.0,
        1.0,
        # 1.0,
        1.0,
        1.0
    ]

    Plot_PDF_fpi(fname_in, fname_out, label, color, zorder, output_format='pdf')
    Plot_PDF_fpi(fname_in, fname_out, label, color, zorder, output_format='png')

    return


if __name__ == '__main__':
    main()
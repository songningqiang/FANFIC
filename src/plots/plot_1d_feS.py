#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@nbi.ku.dk"


"""
plot_1d_feS.py:
    Contains routines to plot the f_e,S PDF

Created: 2020/10/27
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

from extract_1d_feS import filenames_no_ext


def Plot_PDF_feS(fname_in, fname_out, label, color, zorder, 
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
    y_max = 27.0 #17.0

    # First plot the shaded regions
    for file_index in range(len(fname_in)):

        with open(fname_in[file_index]+'.json') as file:
            data = json.load(file)
        fe = data['f_e_S']
        pdf = data['pdf']
        # print(sum(pdf)*(fe[1]-fe[0])) # == 1
        if file_index != 3:
            pdf = savgol_filter(pdf, 21, 3)

        ax.fill_between(fe, pdf, y2=0, ec=color[file_index], 
            fc=color[file_index], zorder=zorder[file_index])

    # Now plot the borders
    for file_index in range(len(fname_in)):

        with open(fname_in[file_index]+'.json') as file:
            data = json.load(file)
        fe = data['f_e_S']
        pdf = data['pdf']
        # print(sum(pdf)*(fe[1]-fe[0])) # == 1
        if file_index != 4:
            pdf = savgol_filter(pdf, 21, 3)

        ax.plot(fe, pdf, c=color[file_index], lw=2.0, label=label[file_index],
            zorder=max(zorder))

    file_index = 3
    with open(fname_in[file_index]+'.json') as file:
        data = json.load(file)
        fe = data['f_e_S']
        pdf = data['pdf']
        pdf = savgol_filter(pdf, 21, 3)
        ax.plot(fe, pdf, c=color[file_index], lw=2.0, label='',
            zorder=max(zorder))

    # Muon-damped
    ax.plot([0.004, 0.004], [0.0, y_max], c='#ff7417', ls='--', lw=2.0,
        zorder=max(zorder)+1)
    ax.annotate( r'$\mu$-damped', 
        xy = (0.02,7.0),
        xycoords='data', color='#ff7417', 
        fontsize=18, horizontalalignment='left', rotation=90.,
        zorder=max(zorder)+1 )
    # Pion decay
    ax.plot([1./3., 1./3.], [0.0, y_max], c='salmon',
        ls='--', lw=2.0,
        zorder=max(zorder)+1)
    ax.annotate( r'$\pi$ decay (assumed real value for $\geq 2020$)', 
        xy = (0.26,7.0), 
        xycoords='data', c='salmon', 
        fontsize=18, horizontalalignment='left', rotation=90.,
        zorder=max(zorder)+1 )
    # Neutron decay
    ax.plot([0.997, 0.997], [0.0, y_max], c='#00a572', ls='--', lw=2.0,
        zorder=max(zorder)+1)
    ax.annotate( r'$n$ decay', 
        xy = (0.95,7.0), 
        xycoords='data', color='#00a572', 
        fontsize=18, horizontalalignment='left', rotation=90.,
        zorder=max(zorder)+1 )
 
    ax.annotate( r'Fixed $f_{\tau,{\rm S}} = 0$', 
        xy = (0.95,0.05), 
        xycoords='axes fraction', color='k', 
        fontsize=18, horizontalalignment='right', rotation=0.,
        zorder=max(zorder)+1 )

    ax.set_xlabel(r'Average $\nu_e$ fraction at sources, $f_{e, {\rm S}}$',
        fontsize=25)
    ax.set_ylabel(r'Probability density $\mathcal{P}(f_{e, {\rm S}})$',
        fontsize=25)

    ax.set_xlim([x_min, x_max])
    ax.set_ylim([y_min, y_max])

    ax.xaxis.set_major_locator(MultipleLocator(0.1))
    ax.xaxis.set_minor_locator(MultipleLocator(0.02))
    ax.yaxis.set_major_locator(MultipleLocator(2.0))
    ax.yaxis.set_minor_locator(MultipleLocator(0.5))
    ax.yaxis.set_major_locator(ticker.MaxNLocator(prune='lower'))

    leg = ax.legend(loc='upper right', ncol=1)
    for line in leg.get_lines():
        line.set_linewidth(15.0)

    pylab.savefig(fname_out+'.'+output_format, bbox_inches='tight', dpi=300)


def main():

    PATH_FIG = '../results/figures/flavor/std/1D/'
    fname_out = PATH_FIG+'feS_recon_pdf_NO'
    if not os.path.exists(PATH_FIG):
        os.mkdir(PATH_FIG)

    PATH_PDF = '../results/pdf/flavor/std/1D/'
    fname_in = [
        PATH_PDF+'output_NUFIT5.0_NO_IceCube2015_pdf',
        PATH_PDF+'output_2020_NO_pdf',
        # PATH_PDF+'output_2028_NO_pdf', 
        PATH_PDF+'output_2040_NO_pdf',
        PATH_PDF+'output_2040_NO_comb_pdf'
    ]

    label = [
        r'2020 (measured):'+'\n'+r'IC (ApJ 1, 98) $\otimes$ NuFit 5.0'
            +'\n'+r'$f_{e, {\rm S}} = 0.05_{-0.05}^{+0.71}~([0.00, 1.00])$',
        r'2020 (projected):'+'\n'+r' IC 8 yr $\otimes$ NuFit 5.0'
            +'\n'+r'$f_{e, {\rm S}} = 0.31_{-0.13}^{+0.08}~([0.00, 0.52])$',
        # r'2028 (projected):'+'\n'+r' IC 15 yr $\otimes$ (NuFit 5.0+JUNO)',
        r'2040 (projected):'+'\n'+r' (IC 15 yr+IC-Gen2 10 yr) $\otimes$'+'\n'+r'(NuFit 5.0+JUNO+DUNE+HK)'
            +'\n'+r'$f_{e, {\rm S}} = 0.33 \pm 0.02~([0.26, 0.40])$',
        r'2040 (projected):'+'\n'+r' (Combined $\nu$ telescopes) $\otimes$'+'\n'+r'(NuFit 5.0+JUNO+DUNE+HK)'
            +'\n'+r'$f_{e, {\rm S}} = 0.33_{-0.01}^{+0.02}~([0.28, 0.38])$',
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

    Plot_PDF_feS(fname_in, fname_out, label, color, zorder, output_format='pdf')
    Plot_PDF_feS(fname_in, fname_out, label, color, zorder, output_format='png')

    return


if __name__ == '__main__':
    main()
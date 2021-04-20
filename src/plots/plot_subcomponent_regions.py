#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@nbi.ku.dk"


"""
plot_subcomponent_regions.py:
    Contains routines to plot the fractions of sources with different production
    mechanisms.  Adapted from plot_flavor_regions.py.

Created: 2020/11/23
Last modified: 2020/11/23
"""


from global_tools import *
import os
import gc

from numpy import *
import numpy as np

from pylab import *
from matplotlib import *
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patheffects as path_effects
import ternary

import json

from scipy.signal import savgol_filter

from extract_flavor_regions import filenames_no_ext


def plot_ternary_subcomponent_fractions_compare_years(fnames_in, fname_out, 
    plot_label, output_format='png'):

    # Open the plot and format it
    mpl.rcParams['xtick.labelsize']=26
    mpl.rcParams['ytick.labelsize']=26
    mpl.rcParams['legend.fontsize']=23
    mpl.rcParams['legend.borderpad']=0.4
    mpl.rcParams['axes.labelpad']=10
    mpl.rcParams['ps.fonttype']=42
    mpl.rcParams['pdf.fonttype']=42

    # if show_benchmark_astro_labels:
    #     colors_gen = ['#2a4d69', '#4b86b4', '#adcbe3', '#D3E4F4', '#63ace5']
    #     colors_010 = ['#964000', '#ff7417', '#ffbf00']
    #     colors_100 = ['#043927', '#00a572', '#d0f0c0']
    #     colors_120 = ['#420d09', '#b80f0a', '#fa8072']
    #     colors = [colors_010, colors_100, colors_120]
    # elif show_benchmark_decay_labels:
    #     colors_gen = ['#2a4d69', '#4b86b4', '#adcbe3', '#D3E4F4', '#63ace5']
    #     colors_010 = [None, '#7E1F86', '#B9A0CF'] # nu_1
    #     # colors_100 = [None, '#8C6D69', '#BFACAA'] # nu_2
    #     # colors_100 = [None, '#DAB558', '#F4E9CD'] # nu_2
    #     # colors_100 = [None, '#40858C', '#7EBDC3'] # nu_2
    #     # colors_100 = [None, '#188B7B', '#1EA896'] # nu_2
    #     colors_100 = [None, '#F5D07A', '#FCF2D9'] # nu_2
    #     colors_120 = [None, '#05668D', '#91C4F2'] # nu_3
    #     colors = [colors_010, colors_100, colors_120]

    # Boundary and gridlines
    fig, tax = ternary.figure(scale=1.0)
    tax.ax.axis("off")
    fig.set_facecolor('w')

    # Draw boundary and gridlines
    tax.boundary(linewidth=1.0)
    tax.gridlines(color="gray", multiple=0.1, linewidth=0.5, ls='-', alpha=0.5)

    fontsize = 15
    label_left = r'Fraction of neutron-decay sources, $k_n$'
    label_right = r'Fraction of muon-damped sources, $k_\mu$'
    label_bottom = r'Fraction of pion-decay sources, $k_\pi$'
    tax.left_axis_label(label_left, fontsize=fontsize, offset=0.16)
    tax.right_axis_label(label_right, fontsize=fontsize, offset=0.16)
    tax.bottom_axis_label(label_bottom, fontsize=fontsize, offset=0.08)

    # Set ticks blr
    tax.ticks(axis='blr', linewidth=0.5, multiple=0.1, offset=0.022, \
        clockwise=False, tick_formats="%.1f")

    ###########################################################################
    # Standard Model regions
    ###########################################################################

    color = [
        '#adcbe3',
        '#63ace5',
        '#4b86b4',
        '#2a4d69',
    ]

    tax.plot([[0.4, 0.3], [0.4, 0.2]], color=color[0])

    for file_index in range(len(fnames_in)):

        # Plot the region for all-fS in the background
        with open(fnames_in[file_index]+'.json') as file:
            data = json.load(file)
        kpi = data['0.90']['f_e_Earth']
        kmu = data['0.90']['f_mu_Earth']
        kn = [1.0-kpi[i]-kmu[i] for i in range(len(kpi))]
        kpi_rot = [ 0.5 * (2.*kpi[i]+kmu[i]) / (kpi[i]+kmu[i]+kn[i]) 
            for i in range(len(kpi)) ]
        kmu_rot = [ sqrt(3.)/2 * (kmu[i]) / (kpi[i]+kmu[i]+kn[i]) 
            for i in range(len(kpi)) ]
        kn_rot = [1.0-kpi_rot[i]-kmu_rot[i] for i in range(len(kpi))]
        tax.ax.fill(kpi_rot, kmu_rot, color=color[file_index])
        # tax.plot([[kpi[i],kmu[i]] for i in range(len(kpi))], 
            # color=color[file_index])


    ###########################################################################
    # Annotations
    ###########################################################################

    tax.ax.annotate( plot_label, \
        xy = (-0.03, 0.845), xycoords='data', color='k', \
        fontsize=10, horizontalalignment='left' )


    # tax.ax.scatter([-0.01], [0.80], \
    #     marker='s', color=color_2020_0, edgecolor='k', \
    #     s=60, linewidths=0.1, zorder=2, alpha=1.)
    # tax.ax.scatter([0.03], [0.80], \
    #     marker='s', color=color_2020_1, edgecolor='k', \
    #     s=60, linewidths=0.1, zorder=2, alpha=1.)
    # tax.ax.scatter([0.07], [0.80], \
    #     marker='s', color=color_2020_2, edgecolor='k', \
    #     s=60, linewidths=0.1, zorder=2, alpha=1.)

    tax.ax.annotate( r'2015', \
        xy = (0.10,0.787), xycoords='data', color='k', \
        fontsize=10, horizontalalignment='left' )
    tax.ax.annotate( r'2020', \
        xy = (0.10,0.737), xycoords='data', color='k', \
        fontsize=10, horizontalalignment='left' )
    tax.ax.annotate( r'2028',
        xy = (0.12,0.695), xycoords='data', color='k', \
        fontsize=10, horizontalalignment='left' )
    tax.ax.annotate( r'2040',
        xy = (0.12,0.66), xycoords='data', color='k', \
        fontsize=10, horizontalalignment='left' )

    tax.clear_matplotlib_ticks()
    tax._redraw_labels()
    ternary.plt.tight_layout()

    ternary.plt.savefig(fname_out+'.'+output_format, dpi=200)

    print("Saved "+fname_out+'.'+output_format)

    plt.close()

    return


def plot_ternary_subcomponent_fractions_single_year(fname_in, fname_out, 
    plot_label, output_format='png'):

    # Open the plot and format it
    mpl.rcParams['xtick.labelsize']=26
    mpl.rcParams['ytick.labelsize']=26
    mpl.rcParams['legend.fontsize']=23
    mpl.rcParams['legend.borderpad']=0.4
    mpl.rcParams['axes.labelpad']=10
    mpl.rcParams['ps.fonttype']=42
    mpl.rcParams['pdf.fonttype']=42

    # if show_benchmark_astro_labels:
    #     colors_gen = ['#2a4d69', '#4b86b4', '#adcbe3', '#D3E4F4', '#63ace5']
    #     colors_010 = ['#964000', '#ff7417', '#ffbf00']
    #     colors_100 = ['#043927', '#00a572', '#d0f0c0']
    #     colors_120 = ['#420d09', '#b80f0a', '#fa8072']
    #     colors = [colors_010, colors_100, colors_120]
    # elif show_benchmark_decay_labels:
    #     colors_gen = ['#2a4d69', '#4b86b4', '#adcbe3', '#D3E4F4', '#63ace5']
    #     colors_010 = [None, '#7E1F86', '#B9A0CF'] # nu_1
    #     # colors_100 = [None, '#8C6D69', '#BFACAA'] # nu_2
    #     # colors_100 = [None, '#DAB558', '#F4E9CD'] # nu_2
    #     # colors_100 = [None, '#40858C', '#7EBDC3'] # nu_2
    #     # colors_100 = [None, '#188B7B', '#1EA896'] # nu_2
    #     colors_100 = [None, '#F5D07A', '#FCF2D9'] # nu_2
    #     colors_120 = [None, '#05668D', '#91C4F2'] # nu_3
    #     colors = [colors_010, colors_100, colors_120]

    # Boundary and gridlines
    fig, tax = ternary.figure(scale=1.0)
    tax.ax.axis("off")
    fig.set_facecolor('w')

    # Draw boundary and gridlines
    tax.boundary(linewidth=1.0)
    tax.gridlines(color="gray", multiple=0.1, linewidth=0.5, ls='-', alpha=0.5)

    fontsize = 15
    # label_left = r'Fraction of $n$-decay sources, $k_n$'
    # label_right = r'Fraction of $\mu$-damped sources, $k_\mu$'
    # label_bottom = r'Fraction of $\pi$-decay sources, $k_\pi$'
    label_left = r'$n$-decay fraction, $k_n$'
    label_right = r'$\mu$-damped fraction, $k_\mu$'
    label_bottom = r'$\pi$-decay fraction, $k_\pi$'
    tax.left_axis_label(label_left, fontsize=fontsize, offset=0.16)
    tax.right_axis_label(label_right, fontsize=fontsize, offset=0.16)
    tax.bottom_axis_label(label_bottom, fontsize=fontsize, offset=0.08)

    # Set ticks blr
    tax.ticks(axis='blr', linewidth=0.5, multiple=0.1, offset=0.022, \
        clockwise=False, tick_formats="%.1f")

    ###########################################################################
    # Plot pdf
    ###########################################################################

    # Each point: [kpi, kmu, kn, posterior]
    with open(fname_in[0]+'.json') as file:
        data = json.load(file)

    kpi = data['f_e_Earth']
    kmu = data['f_mu_Earth']
    kn = [1.0-kpi[i]-kmu[i] for i in range(len(kpi))]
    posterior = data['pdf']
    # posterior_max = max(posterior)
    # posterior_min = min(posterior)
    posterior_max = log10(max(posterior))
    posterior_set = set(posterior)
    posterior_set.remove(min(posterior_set))
    print(min(posterior_set))
    posterior_min = log10(min(posterior_set))

    kpi_mesh, kmu_mesh = np.meshgrid(kpi, kmu)
    posterior = np.reshape(posterior, (len(kpi), len(kmu)))

    cmap = matplotlib.cm.get_cmap('Reds')

    counter = 0 
    for i in range(len(kpi_mesh)):
        for j in range(len(kmu_mesh)):
            kpi_sel = kpi_mesh[i][j]
            kmu_sel = kmu_mesh[i][j]
            posterior_sel = posterior[i][j]
            # print(color_fraction)
            if ((posterior_sel != 0.0) and (kpi_sel + kmu_sel <= 1) \
                and (counter%1 == 0)):
                color_fraction = (log10(posterior_sel)-posterior_min)/(posterior_max-posterior_min)
                # color_fraction = (posterior_sel-posterior_min)/(posterior_max-posterior_min)
                # Produces a flipped image is [kpi, kmu] is used instead
                tax.scatter([[kmu_sel, kpi_sel]], color=cmap(color_fraction),
                    s=0.3)
                counter = 0
            counter += 1

    ###########################################################################
    # Plot contours
    ###########################################################################

    # color = [
    #     '#adcbe3',
    #     '#63ace5',
    #     '#4b86b4',
    #     '#2a4d69',
    # ]

    color = [
        '#adcbe3',
        '#63ace5',
        # '#4b86b4',
        '#2a4d69',
    ]

    with open(fname_in[1]+'.json') as file:
        data = json.load(file)

    for j, cl in enumerate(['0.68', '0.90', '0.997']):

        # Plot the region for all-fS in the background
        kpi = data[cl]['f_e_Earth'] 
        # kpi = savgol_filter(kpi, 91, 5)
        # kpi = savgol_filter(kpi, 3, 1)
        kmu = data[cl]['f_mu_Earth']
        # kmu = savgol_filter(kmu, 91, 5)
        # kmu = savgol_filter(kmu, 3, 1)
        kn = [1.0-kpi[i]-kmu[i] for i in range(len(kpi))]
        kpi_rot = [ 0.5 * (2.*kpi[i]+kmu[i]) / (kpi[i]+kmu[i]+kn[i]) 
            for i in range(len(kpi)) ]
        kmu_rot = [ sqrt(3.)/2 * (kmu[i]) / (kpi[i]+kmu[i]+kn[i]) 
            for i in range(len(kpi)) ]
        kn_rot = [1.0-kpi_rot[i]-kmu_rot[i] for i in range(len(kpi))]
        # tax.ax.fill(kpi_rot, kmu_rot, color=color[j])
        # Produces a flipped image is [kpi, kmu] is used instead
        tax.plot([[kmu[i], kpi[i]] for i in range(len(kpi))], 
            color=color[j])


    ###########################################################################
    # Annotations
    ###########################################################################

    tax.ax.scatter([1.0], [0.0], \
        marker='*', color='gold', edgecolor='k', \
        s=70, linewidths=0.3, zorder=10, alpha=1.)
    tax.ax.scatter([0.13], [0.04], \
        marker='*', color='gold', edgecolor='k', \
        s=70, linewidths=0.3, zorder=10, alpha=1.)
    tax.ax.annotate( r'Assumed real value', 
        xy = (0.15,0.03), xycoords='data', color='k', 
        fontsize=10, horizontalalignment='left' )

    tax.ax.annotate( plot_label, \
        xy = (-0.03, 0.845), xycoords='data', color='k', \
        fontsize=10, horizontalalignment='left' )

    tax.ax.scatter([-0.01], [0.80], \
        marker='s', color=color[0], edgecolor='k', \
        s=60, linewidths=0.1, zorder=2, alpha=1.)
    tax.ax.annotate( r'68\%~C.R.', \
        xy = (0.02,0.787), xycoords='data', color='k', \
        fontsize=10, horizontalalignment='left' )

    tax.ax.scatter([-0.01], [0.75], \
        marker='s', color=color[1], edgecolor='k', \
        s=60, linewidths=0.1, zorder=2, alpha=1.)
    tax.ax.annotate( r'90\%~C.R.', \
        xy = (0.02,0.737), xycoords='data', color='k', \
        fontsize=10, horizontalalignment='left' )

    tax.ax.scatter([-0.01], [0.70], \
        marker='s', color=color[2], edgecolor='k', \
        s=60, linewidths=0.1, zorder=2, alpha=1.)
    tax.ax.annotate( r'99.7\%~C.R.', \
        xy = (0.02,0.687), xycoords='data', color='k', \
        fontsize=10, horizontalalignment='left' )

    # Add the colorbar
    ax1 = fig.add_axes([0.70, 0.87, 0.25, 0.03])
    norm = mpl.colors.Normalize(vmin=posterior_min,\
                                vmax=posterior_max)
    cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap, norm=norm, \
        orientation='horizontal',
        extend='min',
        ticks=[-4.0, -3.0, -2.0, -1.0, 0.0, 1.0])
    cb1.set_label(r'${\rm Log}_{10}(\mathcal{P}(k_i))$', size=10, labelpad=5)
    cb1.ax.tick_params(labelsize=8)
    cb1.ax.set_xticklabels([r'-4', r'-3', r'-2', r'-1', r'0', r'1'])

    tax.clear_matplotlib_ticks()
    tax._redraw_labels()
    ternary.plt.tight_layout()

    ternary.plt.savefig(fname_out+'.'+output_format, dpi=200)

    print("Saved "+fname_out+'.'+output_format)

    plt.close()

    return


# """
PATH = os.getcwd()+'/flavor_region_output_subcomponent_post/3D/'

# fname_in = PATH+'output_NUFIT5.0_NO_IceCube2015_contours'
# fname_in = PATH+'output_2020_NO_contours'
# fname_in = PATH+'output_2028_NO_contours'
# fname_in = PATH+'output_2040_NO_contours'

# fname_in = [
#     PATH+'output_NUFIT5.0_NO_IceCube2015_pdf',
#     PATH+'output_NUFIT5.0_NO_IceCube2015_contours'
# ]
# fname_out = PATH+'output_subcomponent_contours_2015'

# fname_in = [
#     PATH+'output_2020_NO_pdf',
#     PATH+'output_2020_NO_contours'
# ]
# fname_out = PATH+'output_subcomponent_contours_2020'

fname_in = [
    PATH+'output_2040_NO_comb_pdf',
    PATH+'output_2040_NO_comb_contours'
]
fname_out = PATH+'output_subcomponent_contours_2040_comb'
plot_label = r'JUNO + DUNE + HK'+'\n'+r'Combined $\nu$ telescopes'

plot_ternary_subcomponent_fractions_single_year(fname_in, fname_out, 
    plot_label, output_format='pdf')
plot_ternary_subcomponent_fractions_single_year(fname_in, fname_out, 
    plot_label, output_format='png')

quit()
# """

"""
PATH = os.getcwd()+'/flavor_region_output_subcomponent_post/'

fnames_in = [
    PATH+'output_NUFIT5.0_NO_IceCube2015_contours',
    PATH+'output_2020_NO_contours',
    PATH+'output_2028_NO_contours',
    PATH+'output_2040_NO_contours'
]

fname_out = PATH+'output_subcomponent_contours_year_compare'

plot_label = r'NO, all regions 68\% Cr.L.'
# r'NO, upper $\theta_{23}$ octant'+'\n'+r'All regions 99.7\% Cr.L.']

for i, fname in enumerate(fnames_in):
    print(i)

    plot_ternary_subcomponent_fractions_compare_years( \
        fnames_in, fname_out, plot_label, output_format='pdf')

quit()
"""

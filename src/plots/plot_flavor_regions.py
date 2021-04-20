#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@nbi.ku.dk"


"""
plot_flavor_regions.py:
    Contains routines to plot flavor ratios at Earth

Created: 2020/09/16
Last modified: 2021/04/20
"""


import sys
sys.path.insert(0, "../extract_pdf_contours")

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

from extract_flavor_regions import filenames_no_ext


def plot_ternary_flavor_ratios_all_fS_v0(fname, plot_label, order='no', 
    theta23='upper', output_format='png'):

    # Open the plot and format it
    mpl.rcParams['xtick.labelsize']=26
    mpl.rcParams['ytick.labelsize']=26
    mpl.rcParams['legend.fontsize']=23
    mpl.rcParams['legend.borderpad']=0.4
    mpl.rcParams['axes.labelpad']=10
    mpl.rcParams['ps.fonttype']=42
    mpl.rcParams['pdf.fonttype']=42

    colors = ['#2a4d69', '#4b86b4', '#adcbe3', '#D3E4F4', '#63ace5']

    # Boundary and gridlines
    fig, tax = ternary.figure(scale=1.0)
    tax.ax.axis("off")
    fig.set_facecolor('w')

    # Draw boundary and gridlines
    tax.boundary(linewidth=1.0)
    tax.gridlines(color="gray", multiple=0.1, linewidth=0.5, ls='-', alpha=0.5)

    fontsize = 15
    label_left = r'Fraction of $\nu_\tau$, $f_{\tau, \oplus}$'
    label_right = r'Fraction of $\nu_\mu$, $f_{\mu, \oplus}$'
    label_bottom = r'Fraction of $\nu_e$, $f_{e, \oplus}$'
    tax.left_axis_label(label_left, fontsize=fontsize, offset=0.16)
    tax.right_axis_label(label_right, fontsize=fontsize, offset=0.16)
    tax.bottom_axis_label(label_bottom, fontsize=fontsize, offset=0.08)

    # Set ticks blr
    tax.ticks(axis='blr', linewidth=0.5, multiple=0.1, offset=0.022, \
        clockwise=False, tick_formats="%.1f")

    ###########################################################################
    # LIV-like region
    ###########################################################################

    # x = [0.0, 1.0, 0.0]
    # y = [0.0, 0.0, 1.0]
    # z = [1.0-x[i]-y[i] for i in range(len(x))]
    # x_rot = [ 0.5 * (2.*x[i]+y[i]) / (x[i]+y[i]+z[i]) for i in range(len(x)) ]
    # y_rot = [ sqrt(3.)/2 * (y[i]) / (x[i]+y[i]+z[i]) for i in range(len(x)) ]
    # tax.ax.fill(x_rot, y_rot, color=lst_colors[4], alpha=0.4)


    ###########################################################################
    # Decay-like region
    ###########################################################################

    # filename_in = os.getcwd()+"/flavor_ratios_earth_np_decay_3s.dat"
    # contour_sm = Read_Data_File(filename_in)
    # x = contour_sm[1]
    # y = contour_sm[0]
    # z = [1.0-x[i]-y[i] for i in range(len(x))]
    # x_rot = [ 0.5 * (2.*x[i]+y[i]) / (x[i]+y[i]+z[i]) for i in range(len(x)) ]
    # y_rot = [ sqrt(3.)/2 * (y[i]) / (x[i]+y[i]+z[i]) for i in range(len(x)) ]
    # tax.ax.fill(x_rot, y_rot, color=lst_colors[1])


    ###########################################################################
    # Standard Model region
    ###########################################################################

    """
    x = [0.0, 0.1, 0.1 ,0.0]
    y = [0.0, 0.0, 0.1, 0.1]
    z = [1.0-x[i]-y[i] for i in range(len(x))]
    tax.ax.fill(x,y)
    x_rot = [ 0.5 * (2.*x[i]+y[i]) / (x[i]+y[i]+z[i]) for i in range(len(x)) ]
    y_rot = [ sqrt(3.)/2 * (y[i]) / (x[i]+y[i]+z[i]) for i in range(len(x)) ]
    tax.ax.fill(x_rot, y_rot)
    """

    # filename_in = os.getcwd()+"/flavor_ratios_earth_full_var_std_NH_3s_outline.dat"
    # contour_sm = Read_Data_File(filename_in)
    # x = contour_sm[1]
    # y = contour_sm[0]
    # z = [1.0-x[i]-y[i] for i in range(len(x))]
    # x_rot = [ 0.5 * (2.*x[i]+y[i]) / (x[i]+y[i]+z[i]) for i in range(len(x)) ]
    # y_rot = [ sqrt(3.)/2 * (y[i]) / (x[i]+y[i]+z[i]) for i in range(len(x)) ]
    # tax.ax.fill(x_rot, y_rot, color=lst_colors[0])


    ###########################################################################
    # Standard Model regions
    ###########################################################################

    with open(fname+'.json') as file:
        data = json.load(file)

    for j, cl in enumerate(list(data.keys())[::-1]):
        fmu = data[cl]['f_e_Earth']
        fe = data[cl]['f_mu_Earth']
        ftau = [1.0-fe[i]-fmu[i] for i in range(len(fe))]
        fe_rot = [ 0.5 * (2.*fe[i]+fmu[i]) / (fe[i]+fmu[i]+ftau[i]) 
            for i in range(len(fe)) ]
        fmu_rot = [ sqrt(3.)/2 * (fmu[i]) / (fe[i]+fmu[i]+ftau[i]) 
            for i in range(len(fe)) ]
        ftau_rot = [1.0-fe_rot[i]-fmu_rot[i] for i in range(len(fe))]
        # tax.ax.fill(fe_rot, fmu_rot, color=colors[j])
        tax.ax.fill(fe_rot, fmu_rot, color=colors[j])


    ###########################################################################
    # Draw IceCube 2015 contours
    ###########################################################################

    """
    # lst_filenames_ic =  [ \
    #     os.getcwd()+"/flavor_ratios_earth_ic_2015_bf.dat", \
    #     os.getcwd()+"/flavor_ratios_earth_ic_2015_1s.dat", \
    #     # os.getcwd()+"/flavor_ratios_earth_ic_2015_2s_a.dat", \
    #     # os.getcwd()+"/flavor_ratios_earth_ic_2015_2s_b.dat", \
    #     os.getcwd()+"/flavor_ratios_earth_ic_2015_3s_a.dat", \
    #     os.getcwd()+"/flavor_ratios_earth_ic_2015_3s_b.dat"]

    lst_filenames_ic =  [ \
        os.getcwd()+"/icgen2_icrc2017_pion_68cl.dat",
        os.getcwd()+"/icgen2_icrc2017_pion_90cl.dat",
        os.getcwd()+"/icgen2_icrc2017_muon_68cl.dat",
        os.getcwd()+"/icgen2_icrc2017_muon_90cl.dat"]

    lst_color = ['gold', 'gold', 'w', 'w'] #'y'
    lst_ls = ['--', '-', '--', '-']

    # Each row: f_eE, f_mE
    for file_index in range(len(lst_filenames_ic)):

        # Read data file with IceCube flavor ratios
        lst_ratios_earth = Read_Data_File(lst_filenames_ic[file_index])
        # Recast as [f_eE, f_mE, f_tE]
        lst_ratios_earth = \
            [[lst_ratios_earth[1][i], lst_ratios_earth[2][i],\
            lst_ratios_earth[0][i]] \
            for i in range(len(lst_ratios_earth[0]))]

        tax.plot(lst_ratios_earth, color=lst_color[file_index],
                ls=lst_ls[file_index], linewidth=2.5, zorder=1)

    """

    # IceCube, 8 yr
    with open('icecube_8yr.json') as json_file:
        data = json.load(json_file)
        fraction_nue = data['data']['nue_fraction']
        fraction_numu = data['data']['numu_fraction']
        fraction_nutau = \
            [1.0-fraction_nue[i]-fraction_numu[i] \
            for i in range(len(fraction_nue))]
        confidence_level = data['data']['confidence_level']
        cs = tax.ax.contour(fraction_nue, fraction_numu, confidence_level,
                levels=[68, ], colors=None, linestyles=None, alpha=0.0)
        vertices = cs.collections[0].get_paths()[0].vertices
        x = vertices[:,0]
        y = vertices[:,1]
        z = [1.0-x[i]-y[i] for i in range(len(vertices))]
        ratios_earth = [[x[i], y[i], z[i]] for i in range(len(x))]
        tax.plot(ratios_earth, color='tomato',
                ls='-', linewidth=1.5, zorder=1)

    # IceCube, 15 yr
    with open('icecube_15yr.json') as json_file:
        data = json.load(json_file)
        fraction_nue = data['data']['nue_fraction']
        fraction_numu = data['data']['numu_fraction']
        fraction_nutau = \
            [1.0-fraction_nue[i]-fraction_numu[i] \
            for i in range(len(fraction_nue))]
        confidence_level = data['data']['confidence_level']
        cs = tax.ax.contour(fraction_nue, fraction_numu, confidence_level,
                levels=[68, ], colors=None, linestyles=None, alpha=0.0)
        vertices = cs.collections[0].get_paths()[0].vertices
        x = vertices[:,0]
        y = vertices[:,1]
        z = [1.0-x[i]-y[i] for i in range(len(vertices))]
        ratios_earth = [[x[i], y[i], z[i]] for i in range(len(x))]
        tax.plot(ratios_earth, color='tan',
                ls='-', linewidth=1.5, zorder=1)

    # IceCube-Gen2 10-yr + IceCube 15-yr
    with open('icecube+gen2-inice.json') as json_file:
        data = json.load(json_file)
        fraction_nue = data['data']['nue_fraction']
        fraction_numu = data['data']['numu_fraction']
        fraction_nutau = \
            [1.0-fraction_nue[i]-fraction_numu[i] \
            for i in range(len(fraction_nue))]
        confidence_level = data['data']['confidence_level']
        cs = tax.ax.contour(fraction_nue, fraction_numu, confidence_level,
                levels=[68, ], colors=None, linestyles=None, alpha=0.0)
        vertices = cs.collections[0].get_paths()[0].vertices
        x = vertices[:,0]
        y = vertices[:,1]
        z = [1.0-x[i]-y[i] for i in range(len(vertices))]
        ratios_earth = [[x[i], y[i], z[i]] for i in range(len(x))]
        tax.plot(ratios_earth, color='gold',
                ls='-', linewidth=1.5, zorder=1)

    ###########################################################################
    # Annotations
    ###########################################################################

    tax.ax.annotate( plot_label, \
        xy = (-0.03, 0.845), xycoords='data', color='k', \
        fontsize=10, horizontalalignment='left' )

    tax.ax.scatter([-0.01], [0.80], \
        marker='s', color=colors[3], edgecolor='k', \
        s=60, linewidths=0.1, zorder=2, alpha=1.)
    tax.ax.annotate( r'68\%~C.L.', \
        xy = (0.02,0.787), xycoords='data', color='k', \
        fontsize=10, horizontalalignment='left' )

    tax.ax.scatter([-0.01], [0.75], \
        marker='s', color=colors[2], edgecolor=colors[2], \
        s=60, linewidths=0.1, zorder=2, alpha=1.)
    tax.ax.annotate( r'90\%~C.L.', \
        xy = (0.02,0.737), xycoords='data', color='k', \
        fontsize=10, horizontalalignment='left' )

    tax.ax.scatter([-0.01], [0.70], \
        marker='s', color=colors[1], edgecolor=colors[1], \
        s=60, linewidths=0.1, zorder=2, alpha=1.)
    tax.ax.annotate( r'95\%~C.L.', \
        xy = (0.02,0.687), xycoords='data', color='k', \
        fontsize=10, horizontalalignment='left' )

    tax.ax.scatter([-0.01], [0.65], \
        marker='s', color=colors[0], edgecolor=colors[0], \
        s=60, linewidths=0.1, zorder=2, alpha=1.)
    tax.ax.annotate( r'99.7\%~C.L.', \
        xy = (0.02,0.637), xycoords='data', color='k', \
        fontsize=10, horizontalalignment='left' )

    # tax.ax.scatter([-0.01], [0.872], \
    #     marker='s', color=colors[0], edgecolor=colors[0], \
    #     s=60, linewidths=0.1, zorder=2, alpha=1.)
    # tax.ax.annotate( r'Standard Model ($\nu$SM)', \
    #     xy = (0.02,0.86), xycoords='data', color='k', \
    #     fontsize=10, horizontalalignment='left' )

    # tax.ax.scatter([-0.01], [0.79], \
    #     marker='s', color=colors[1], edgecolor=colors[1], \
    #     s=60, linewidths=0.1, zorder=2, alpha=1.)
    # tax.ax.annotate( r'New physics:', \
    #     xy = (0.02,0.80), xycoords='data', color='k', \
    #     fontsize=10, horizontalalignment='left' )
    # tax.ax.annotate( r'$\nu$ decay-like', \
    #     xy = (0.02,0.755), xycoords='data', color='k', \
    #     fontsize=10, horizontalalignment='left' )

    # tax.ax.scatter([-0.01], [0.685], \
    #     marker='s', color=colors[4], edgecolor=colors[4], \
    #     s=60, linewidths=0.1, zorder=2, alpha=0.4)
    # tax.ax.annotate( r'New physics:', \
    #     xy = (0.02,0.695), xycoords='data', color='k', \
    #     fontsize=10, horizontalalignment='left' )
    # tax.ax.annotate( r'Lorentz violation', \
    #     xy = (0.02,0.65), xycoords='data', color='k', \
    #     fontsize=10, horizontalalignment='left' )

    tax.plot([[0.03, 0.14], [0.07, 0.14]], color='tomato',
                ls='-', linewidth=1.5, zorder=1)
    tax.ax.annotate( r'IceCube 8 yr (68\%)', \
        xy = (0.15,0.11), xycoords='data', color='k', \
        fontsize=10, horizontalalignment='left' )
    tax.plot([[0.05, 0.095], [0.09, 0.095]], color='tan',
                ls='-', linewidth=1.5, zorder=1)
    tax.ax.annotate( r'IceCube 15 yr (68\%)', \
        xy = (0.15,0.07), xycoords='data', color='k', \
        fontsize=10, horizontalalignment='left' )
    tax.plot([[0.073, 0.05], [0.113, 0.05]], color='gold',
                ls='-', linewidth=1.5, zorder=1)
    tax.ax.annotate( r'IceCube 15 yr + Gen2 10 yr (68\%)', \
        xy = (0.15,0.03), xycoords='data', color='k', \
        fontsize=10, horizontalalignment='left' )


    # tax.ax.annotate( r'$1\sigma$', \
    #     xy = (0.58,0.56), xycoords='data', color='w', rotation=15., \
    #     fontsize=14, horizontalalignment='left',
    #     path_effects=[path_effects.withSimplePatchShadow(offset=(1,-1),
    #                 shadow_rgbFace='k', alpha=0.8)] )
    # tax.ax.annotate( r'$3\sigma$', \
    #     xy = (0.51,0.70), xycoords='data', color='w', rotation=15., \
    #     fontsize=14, horizontalalignment='left',
    #     path_effects=[path_effects.withSimplePatchShadow(offset=(1,-1),
    #                 shadow_rgbFace='k', alpha=0.8)] )


    # text = plt.text(0.092, 0.10, 'IceCube 2015', rotation=5.,
    #             fontdict={'color': 'w'},
    #             path_effects=[path_effects.withSimplePatchShadow(offset=(1,-1),
    #                 shadow_rgbFace='k', alpha=0.8)])


    ####################################################################
    # Add regions of selected compositions
    ####################################################################

    if (order.lower() == 'no'):
        if (theta23 == 'upper'):
            flavor_ratios_pion \
                = [0.2982716010150718, 0.359320837204981, 0.34240756177994724]
            flavor_ratios_muon \
                = [0.17128846521772986, 0.4533370231986065, 0.37537451158366353]
            flavor_ratios_neutron \
                = [0.5522378726097557, 0.17128846521772986, 0.2764736621725147]
        elif (theta23 == 'lower'):
            flavor_ratios_pion \
                = [0.32800656074697265, 0.34529257725052065, 0.3267008620025065]
            flavor_ratios_muon \
                = [0.2158909048155812, 0.40999341346799034, 0.3741156817164281]
            flavor_ratios_neutron \
                = [0.5522378726097557, 0.2158909048155812, 0.23187122257466325]
        elif (theta23 == 'max'):
            flavor_ratios_pion \
                = [0.31656024191908616, 0.3474455438823001, 0.3359942141986137]
            flavor_ratios_muon \
                = [0.19872142657375141, 0.4218076025365745, 0.37947097088967396]
            flavor_ratios_neutron \
                = [0.5522378726097557, 0.19872142657375141, 0.24904070081649313]
    elif (order.lower() == 'io'):
        if (theta23 == 'upper'):
            flavor_ratios_pion \
                = [0.3180471687709684, 0.3474505778780882, 0.3345022533509433]
            flavor_ratios_muon \
                = [0.2010801776404412, 0.4206357779969117, 0.378284044362647]
            flavor_ratios_neutron \
                = [0.5519811510320229, 0.2010801776404412, 0.24693867132753594]
        elif (theta23 == 'lower'):
            flavor_ratios_pion \
                = [0.34834738484344396, 0.3277980374633891, 0.32385457769316683]
            flavor_ratios_muon \
                = [0.2465305017491545, 0.36843180532050646, 0.38503769293033896]
            flavor_ratios_neutron \
                = [0.5519811510320229, 0.2465305017491545, 0.20148834721882264]
        elif (theta23 == 'max'):
            flavor_ratios_pion \
                = [0.33698743757165983, 0.3317785331561326, 0.33123402927220735]
            flavor_ratios_muon \
                = [0.22949058084147833, 0.3829225093134597, 0.3875869098450617]
            flavor_ratios_neutron \
                = [0.5519811510320229, 0.22949058084147833, 0.21852826812649873]

    lst_flavor_ratios_bf = [flavor_ratios_pion, flavor_ratios_muon, \
                            flavor_ratios_neutron]
    lst_color = ['salmon', 'orange', 'limegreen']
    lst_markers = ['o', 's', '^']

    for i in range(len(lst_flavor_ratios_bf)):
        # flavor_ratios = [lst_flavor_ratios_bf[i][0], lst_flavor_ratios_bf[i][2],
            # lst_flavor_ratios_bf[i][1]] 
        tax.scatter([lst_flavor_ratios_bf[i]], \
        # tax.scatter([flavor_ratios], \
            marker=lst_markers[i], color=lst_color[i], \
            edgecolor=lst_color[i], s=40, linewidths=0.2, zorder=4, \
            alpha=1.)

    # Annotate
    tax.scatter([[0.20, 0.98, 0.]], \
        marker='o', color='salmon', edgecolor='salmon', \
        s=40, linewidths=0.1, zorder=2, alpha=1.)
    tax.ax.annotate( r'$\pi$ decay: $\left( 1:2:0 \right)_{\rm S}$', \
        xy = (0.72,0.83), xycoords='data', color='k', \
        fontsize=10, horizontalalignment='left' )
    tax.scatter([[0.235, 0.91, 0.]], \
        marker='s', color='orange', edgecolor='orange', \
        s=40, linewidths=0.1, zorder=2, alpha=1.)
    tax.ax.annotate( r'$\mu$-damped: $\left( 0:1:0 \right)_{\rm S}$', \
        xy = (0.72,0.77), xycoords='data', color='k', \
        fontsize=10, horizontalalignment='left' )
    tax.scatter([[0.27, 0.843, 0.]], \
        marker='^', color='limegreen', edgecolor='limegreen', \
        s=40, linewidths=0.1, zorder=2, alpha=1.)
    tax.ax.annotate( r'$n$ decay: $\left( 1:0:0 \right)_{\rm S}$', \
        xy = (0.72,0.71), xycoords='data', color='k', \
        fontsize=10, horizontalalignment='left' )


    tax.clear_matplotlib_ticks()
    tax._redraw_labels()
    ternary.plt.tight_layout()

    ternary.plt.savefig(fname, dpi=200)
    print("Saved "+fname)

    plt.close()

    return


def plot_ternary_flavor_ratios_all_fS_v1(fname, plot_label, order='no', 
    theta23='upper', year='2020', fixed_deltacp=None, nufit_version=None,
    show_benchmarks=True, show_benchmark_astro_labels=True,
    show_benchmark_decay_labels=False, output_format='png'):

    # Open the plot and format it
    mpl.rcParams['xtick.labelsize']=26
    mpl.rcParams['ytick.labelsize']=26
    mpl.rcParams['legend.fontsize']=23
    mpl.rcParams['legend.borderpad']=0.4
    mpl.rcParams['axes.labelpad']=10
    mpl.rcParams['ps.fonttype']=42
    mpl.rcParams['pdf.fonttype']=42

    colors = ['#2a4d69', '#4b86b4', '#adcbe3', '#D3E4F4', '#63ace5']

    # Boundary and gridlines
    fig, tax = ternary.figure(scale=1.0)
    tax.ax.axis("off")
    fig.set_facecolor('w')

    # Draw boundary and gridlines
    tax.boundary(linewidth=1.0)
    tax.gridlines(color="gray", multiple=0.1, linewidth=0.5, ls='-', alpha=0.5)

    fontsize = 15
    label_left = r'Fraction of $\nu_\tau$, $f_{\tau, \oplus}$'
    label_right = r'Fraction of $\nu_\mu$, $f_{\mu, \oplus}$'
    label_bottom = r'Fraction of $\nu_e$, $f_{e, \oplus}$'
    tax.left_axis_label(label_left, fontsize=fontsize, offset=0.16)
    tax.right_axis_label(label_right, fontsize=fontsize, offset=0.16)
    tax.bottom_axis_label(label_bottom, fontsize=fontsize, offset=0.08)

    # Set ticks blr
    tax.ticks(axis='blr', linewidth=0.5, multiple=0.1, offset=0.022, \
        clockwise=False, tick_formats="%.1f")

    ###########################################################################
    # Standard Model regions
    ###########################################################################

    with open(fname+'.json') as file:
        data = json.load(file)

    for j, cl in enumerate(['0.68', '0.95', '0.997'][::-1]): #enumerate(list(data.keys())[::-1]):
        fmu = data[cl]['f_e_Earth']
        fe = data[cl]['f_mu_Earth']
        ftau = [1.0-fe[i]-fmu[i] for i in range(len(fe))]
        fe_rot = [ 0.5 * (2.*fe[i]+fmu[i]) / (fe[i]+fmu[i]+ftau[i]) 
            for i in range(len(fe)) ]
        fmu_rot = [ sqrt(3.)/2 * (fmu[i]) / (fe[i]+fmu[i]+ftau[i]) 
            for i in range(len(fe)) ]
        ftau_rot = [1.0-fe_rot[i]-fmu_rot[i] for i in range(len(fe))]
        # tax.ax.fill(fe_rot, fmu_rot, color=colors[j])
        tax.ax.fill(fe_rot, fmu_rot, color=colors[j])


    ###########################################################################
    # IceCube contours
    ###########################################################################

    if (year == '2020' or year is None):
        # IceCube, 8 yr
        # file_ic_contour = 'icecube_8yr.json'
        file_ic_contour = 'icecube_8yr_recenter.json'
        ic_contour_label = r'IceCube 8 yr (68\%, 95\%, 99.7\% C.R.)'
    elif (year == '2030'):
        # IceCube, 15 yr
        # file_ic_contour = 'icecube_15yr.json'
        file_ic_contour = 'icecube_15yr_recenter.json'
        ic_contour_label = r'IceCube 15 yr (68\%, 95\%, 99.7\% C.R.)'
    elif (year == '2040'):
        # IceCube-Gen2 10-yr + IceCube 15-yr
        # file_ic_contour = 'icecube+gen2-inice.json'
        file_ic_contour = 'icecube+gen2-inice_recenter.json'
        ic_contour_label = r'IceCube 15 yr + Gen2 10 yr (68\%, 95\%, 99.7\% C.R.)'
    elif (year == '2040_comb'):
        # All neutrino telescopes combined
        file_ic_contour = 'combineexp_recenter.json'
        ic_contour_label = r'Combined $\nu$ telescopes (68\%, 95\%, 99.7\% C.R.)'

    tax.plot([[0.073, 0.05], [0.113, 0.05]], color='C0',
                ls='-', linewidth=1.5, zorder=1)
    tax.ax.annotate( ic_contour_label, 
        xy = (0.15,0.03), xycoords='data', color='k', 
        fontsize=10, horizontalalignment='left' )

    with open(file_ic_contour) as json_file:
        data = json.load(json_file)
        fraction_nue = data['data']['nue_fraction']
        fraction_numu = data['data']['numu_fraction']
        fraction_nutau = \
            [1.0-fraction_nue[i]-fraction_numu[i] \
            for i in range(len(fraction_nue))]
        confidence_level = data['data']['confidence_level']
        cs = tax.ax.contour(fraction_nue, fraction_numu, confidence_level,
                levels=[68., 95., 99.7], colors=None, linestyles=None, alpha=0.0)
        for j in range(len(cs.collections)):
            vertices = cs.collections[j].get_paths()[0].vertices
            x = vertices[:,0]
            y = vertices[:,1]
            z = [1.0-x[i]-y[i] for i in range(len(vertices))]
            ratios_earth = [[x[i], y[i], z[i]] for i in range(len(x))]
            tax.plot(ratios_earth, color='C0',
                    ls='-', linewidth=1., zorder=5)

    ###########################################################################
    # Annotations
    ###########################################################################

    tax.ax.annotate( plot_label, \
        xy = (-0.03, 0.845), xycoords='data', color='k', \
        fontsize=10, horizontalalignment='left' )

    tax.ax.scatter([-0.01], [0.80], \
        marker='s', color=colors[2], edgecolor='k', \
        s=60, linewidths=0.1, zorder=2, alpha=1.)
    tax.ax.annotate( r'68\%~C.R.', \
        xy = (0.02,0.787), xycoords='data', color='k', \
        fontsize=10, horizontalalignment='left' )

    tax.ax.scatter([-0.01], [0.75], \
        marker='s', color=colors[1], edgecolor=colors[2], \
        s=60, linewidths=0.1, zorder=2, alpha=1.)
    tax.ax.annotate( r'95\%~C.R.', \
        xy = (0.02,0.737), xycoords='data', color='k', \
        fontsize=10, horizontalalignment='left' )

    tax.ax.scatter([-0.01], [0.70], \
        marker='s', color=colors[0], edgecolor=colors[1], \
        s=60, linewidths=0.1, zorder=2, alpha=1.)
    tax.ax.annotate( r'99.7\%~C.R.', \
        xy = (0.02,0.687), xycoords='data', color='k', \
        fontsize=10, horizontalalignment='left' )

    ####################################################################
    # Add regions of selected compositions
    ####################################################################

    if (nufit_version is None):

        if (order.lower() == 'no'):
            if (theta23 == 'upper'):
                if (fixed_deltacp is None):
                    # Best-fit flavor ratios for benchmark points
                    flavor_ratios_pion \
                        = [0.2982716010150718, 0.359320837204981, 0.34240756177994724]
                    flavor_ratios_muon \
                        = [0.17128846521772986, 0.4533370231986065, 0.37537451158366353]
                    flavor_ratios_neutron \
                        = [0.5522378726097557, 0.17128846521772986, 0.2764736621725147]
                    # Best-fit flavor content of nu_1, nu_2, nu_3
                    flavor_cont_nu1 \
                        = [0.6808662622895465, 0.07368823020070894, 0.24544550750974464]
                    flavor_cont_nu2 \
                        = [0.29692740578832083, 0.3659953978327661, 0.33707719637891315]
                    flavor_cont_nu3 \
                        = [0.022206331922132744, 0.560316371966525, 0.41747729611134216]
                else:
                    if (fixed_deltacp == 0.0):
                        flavor_ratios_pion \
                            = [0.3322161727912954, 0.3480562782080484, 0.319727549000656]
                        flavor_ratios_muon \
                            = [0.22220532288206532, 0.41098175587104, 0.3668129212468944]
                        flavor_ratios_neutron \
                            = [0.5522378726097557, 0.22220532288206532, 0.22555680450817917]
                    elif (fixed_deltacp == 0.5*np.pi):
                        flavor_ratios_pion \
                            = [0.3148647971219528, 0.34795485221718947, 0.33718035066085766]
                        flavor_ratios_muon \
                            = [0.1961782593780514, 0.4238431486367586, 0.3799785919851899]
                        flavor_ratios_neutron \
                            = [0.5522378726097557, 0.1961782593780514, 0.2515838680121932]
                    elif (fixed_deltacp == np.pi):
                        flavor_ratios_pion \
                            = [0.3148647971219528, 0.34795485221718947, 0.33718035066085766]
                        flavor_ratios_muon \
                            = [0.1961782593780514, 0.4238431486367586, 0.3799785919851899]
                        flavor_ratios_neutron \
                            = [0.5522378726097557, 0.1961782593780514, 0.2515838680121932]
                    else:
                        print("Invalid value of fixed_deltacp")
                        quit()
            elif (theta23 == 'lower'):
                if (fixed_deltacp is None):
                    flavor_ratios_pion \
                        = [0.32800656074697265, 0.34529257725052065, 0.3267008620025065]
                    flavor_ratios_muon \
                        = [0.2158909048155812, 0.40999341346799034, 0.3741156817164281]
                    flavor_ratios_neutron \
                        = [0.5522378726097557, 0.2158909048155812, 0.23187122257466325]
                else:
                    print("Cannot choose fixed_deltacp != None for theta = 'lower'")
                    quit()
            elif (theta23 == 'max'):
                if (fixed_deltacp is None):
                    flavor_ratios_pion \
                        = [0.31656024191908616, 0.3474455438823001, 0.3359942141986137]
                    flavor_ratios_muon \
                        = [0.19872142657375141, 0.4218076025365745, 0.37947097088967396]
                    flavor_ratios_neutron \
                        = [0.5522378726097557, 0.19872142657375141, 0.24904070081649313]
                else:
                    print("Cannot choose fixed_deltacp != None for theta = 'max'")
                    quit()
        elif (order.lower() == 'io'):
            if (theta23 == 'upper'):
                if (fixed_deltacp is None):
                    # Best-fit flavor ratios for benchmark points
                    flavor_ratios_pion \
                        = [0.3180471687709684, 0.3474505778780882, 0.3345022533509433]
                    flavor_ratios_muon \
                        = [0.2010801776404412, 0.4206357779969117, 0.378284044362647]
                    flavor_ratios_neutron \
                        = [0.5519811510320229, 0.2010801776404412, 0.24693867132753594]
                    # Best-fit flavor content of nu_1, nu_2, nu_3
                    flavor_cont_nu1 \
                        = [0.6806215638904127, 0.15226374101164764, 0.16711469509793975]
                    flavor_cont_nu2 \
                        = [0.2970466806723904, 0.2858044688729501, 0.41714885045465944]
                    flavor_cont_nu3 \
                        = [0.02233175543719699, 0.5619317901154023, 0.41573645444740065]
                else:
                    if (fixed_deltacp == 0.0):
                        flavor_ratios_pion \
                            = [0.33181153325437795, 0.34861177147789035, 0.3195766952677316]
                        flavor_ratios_muon \
                            = [0.22172672436555554, 0.41205429503405777, 0.3662189806003866]
                        flavor_ratios_neutron \
                            = [0.5519811510320229, 0.22172672436555554, 0.22629212460242162]
                    elif (fixed_deltacp == 1.0*np.pi):
                        flavor_ratios_pion \
                            = [0.29705684749618344, 0.3605546628791725, 0.34238848962464397]
                        flavor_ratios_muon \
                            = [0.16959469572826372, 0.4560346464546269, 0.3743706578171092]
                        flavor_ratios_neutron \
                            = [0.5519811510320229, 0.16959469572826372, 0.2784241532397134]
                    elif (fixed_deltacp == 1.5*np.pi):
                        flavor_ratios_pion \
                            = [0.3144341903752807, 0.34842595462849185, 0.3371398549962274]
                        flavor_ratios_muon \
                            = [0.19566071004690963, 0.424808576919283, 0.37953071303380737]
                        flavor_ratios_neutron \
                            = [0.5519811510320229, 0.19566071004690963, 0.25235813892106745]
                    else:
                        print("Invalid value of fixed_deltacp")
                        quit()
            elif (theta23 == 'lower'):
                if (fixed_deltacp is None):
                    flavor_ratios_pion \
                        = [0.34834738484344396, 0.3277980374633891, 0.32385457769316683]
                    flavor_ratios_muon \
                        = [0.2465305017491545, 0.36843180532050646, 0.38503769293033896]
                    flavor_ratios_neutron \
                        = [0.5519811510320229, 0.2465305017491545, 0.20148834721882264]
                else:
                    print("Cannot choose fixed_deltacp != None for theta = 'lower'")
                    quit()
            elif (theta23 == 'max'):
                if (fixed_deltacp is None):
                    flavor_ratios_pion \
                        = [0.33698743757165983, 0.3317785331561326, 0.33123402927220735]
                    flavor_ratios_muon \
                        = [0.22949058084147833, 0.3829225093134597, 0.3875869098450617]
                    flavor_ratios_neutron \
                        = [0.5519811510320229, 0.22949058084147833, 0.21852826812649873]
                else:
                    print("Cannot choose fixed_deltacp != None for theta = 'max'")
                    quit()

    else: # nufit_version is not None

        # if (nufit_version == 'nufit_v1.0'):
        #     if (order.lower() == 'no'):
        #         flavor_ratios_pion \
        #             = xxx
        #         flavor_ratios_muon \
        #             = xxx
        #         flavor_ratios_neutron \
        #             = xxx
        #     elif (order.lower() == 'io'):
        #         flavor_ratios_pion \
        #             = xxx
        #         flavor_ratios_muon \
        #             = xxx
        #         flavor_ratios_neutron \
        #             = xxx
        # if (nufit_version == 'nufit_v1.1'):
        #     if (order.lower() == 'no'):
        #         flavor_ratios_pion \
        #             = xxx
        #         flavor_ratios_muon \
        #             = xxx
        #         flavor_ratios_neutron \
        #             = xxx
        #     elif (order.lower() == 'io'):
        #         flavor_ratios_pion \
        #             = xxx
        #         flavor_ratios_muon \
        #             = xxx
        #         flavor_ratios_neutron \
        #             = xxx
        # if (nufit_version == 'nufit_v1.2'):
        #     if (order.lower() == 'no'):
        #         flavor_ratios_pion \
        #             = xxx
        #         flavor_ratios_muon \
        #             = xxx
        #         flavor_ratios_neutron \
        #             = xxx
        #     elif (order.lower() == 'io'):
        #         flavor_ratios_pion \
        #             = xxx
        #         flavor_ratios_muon \
        #             = xxx
        #         flavor_ratios_neutron \
        #             = xxx
        # if (nufit_version == 'nufit_v1.3'):
        #     if (order.lower() == 'no'):
        #         flavor_ratios_pion \
        #             = xxx
        #         flavor_ratios_muon \
        #             = xxx
        #         flavor_ratios_neutron \
        #             = xxx
        #     elif (order.lower() == 'io'):
        #         flavor_ratios_pion \
        #             = xxx
        #         flavor_ratios_muon \
        #             = xxx
        #         flavor_ratios_neutron \
        #             = xxx
        if (nufit_version == 'nufit_v2.0'):
            if (order.lower() == 'no'):
                flavor_ratios_pion \
                    = [0.3556613043722851, 0.3255837945733705, 0.318754901054344]
                flavor_ratios_muon \
                    = [0.25727620733858775, 0.3597375881907619, 0.38298620447064996]
                flavor_ratios_neutron \
                    = [0.5524314984396799, 0.25727620733858775, 0.19029229422173227]
            elif (order.lower() == 'io'):
                flavor_ratios_pion \
                    = [0.3085776926666689, 0.3516915412977889, 0.3397307660355422]
                flavor_ratios_muon \
                    = [0.18670502760224328, 0.43418479814556177, 0.37911017425219495]
                flavor_ratios_neutron \
                    = [0.5523230227955201, 0.18670502760224328, 0.2609719496022368]
        if (nufit_version == 'nufit_v2.1'):
            if (order.lower() == 'no'):
                flavor_ratios_pion \
                    = [0.31499629123871653, 0.3481219497568094, 0.3368817590044737]
                flavor_ratios_muon \
                    = [0.19770980806911492, 0.42332802060065666, 0.37896217133022814]
                flavor_ratios_neutron \
                    = [0.5495692575779199, 0.19770980806911492, 0.252720934352965]
            elif (order.lower() == 'io'):
                flavor_ratios_pion \
                    = [0.3090660431370772, 0.35143623465391693, 0.33949772220900554]
                flavor_ratios_muon \
                    = [0.1890301412593759, 0.4326392813511875, 0.3783305773894363]
                flavor_ratios_neutron \
                    = [0.5491378468924799, 0.1890301412593759, 0.26183201184814403]
        if (nufit_version == 'nufit_v2.2'):
            if (order.lower() == 'no'):
                flavor_ratios_pion \
                    = [0.3541798378010116, 0.3264941108322114, 0.31932605136677694]
                flavor_ratios_muon \
                    = [0.2564473535897559, 0.3615174894534392, 0.382035156956805]
                flavor_ratios_neutron \
                    = [0.5496448062235231, 0.2564473535897559, 0.19390784018672091]
            elif (order.lower() == 'io'):
                flavor_ratios_pion \
                    = [0.31156306080426477, 0.3512832938335091, 0.3371536453622259]
                flavor_ratios_muon \
                    = [0.1925869393553972, 0.4306314710725651, 0.37678158957203756]
                flavor_ratios_neutron \
                    = [0.5495153037019999, 0.1925869393553972, 0.25789775694260275]
        if (nufit_version == 'nufit_v3.0'):
            if (order.lower() == 'no'):
                flavor_ratios_pion \
                    = [0.34569828207852565, 0.33114156111530313, 0.32316015680617094]
                flavor_ratios_muon \
                    = [0.24300229079762686, 0.37521119627414135, 0.3817865129282315]
                flavor_ratios_neutron \
                    = [0.5510902646403233, 0.24300229079762686, 0.20590744456204982]
            elif (order.lower() == 'io'):
                flavor_ratios_pion \
                    = [0.31324903482344135, 0.35147907508403603, 0.3352718900925226]
                flavor_ratios_muon \
                    = [0.19439875631305437, 0.4300192344695269, 0.3755820092174188]
                flavor_ratios_neutron \
                    = [0.5509495918442153, 0.19439875631305437, 0.2546516518427304]
        if (nufit_version == 'nufit_v3.1'):
            if (order.lower() == 'no'):
                flavor_ratios_pion \
                    = [0.3053534737546712, 0.3527746876435402, 0.34187183860178844]
                flavor_ratios_muon \
                    = [0.18301214319568435, 0.43765595986746814, 0.37933189693684743]
                flavor_ratios_neutron \
                    = [0.5500361348726449, 0.18301214319568435, 0.2669517219316705]
            elif (order.lower() == 'io'):
                flavor_ratios_pion \
                    = [0.3182553364828055, 0.34681283459069917, 0.3349318289264954]
                flavor_ratios_muon \
                    = [0.20245670395710264, 0.4189908999074975, 0.37855239613540004]
                flavor_ratios_neutron \
                    = [0.5498526015342113, 0.20245670395710264, 0.24769069450868617]
        if (nufit_version == 'nufit_v3.2'):
            if (order.lower() == 'no'):
                flavor_ratios_pion \
                    = [0.3135479127848382, 0.346028878345814, 0.34042320886934785]
                flavor_ratios_muon \
                    = [0.1953631853698009, 0.4213617248338205, 0.38327508979637864]
                flavor_ratios_neutron \
                    = [0.5499173676149128, 0.1953631853698009, 0.2547194470152863]
            elif (order.lower() == 'io'):
                flavor_ratios_pion \
                    = [0.32198676464987386, 0.3423370406991942, 0.3356761946509318]
                flavor_ratios_muon \
                    = [0.2081347790102386, 0.409438171543672, 0.38242704944608913]
                flavor_ratios_neutron \
                    = [0.5496907359291444, 0.2081347790102386, 0.24217448506061715]
        if (nufit_version == 'nufit_v4.0'):
            if (order.lower() == 'no'):
                flavor_ratios_pion \
                    = [0.2988473594059442, 0.3589894623100733, 0.3421631782839824]
                flavor_ratios_muon \
                    = [0.17459388557291625, 0.45118725067865184, 0.3742188637484317]
                flavor_ratios_neutron \
                    = [0.5473543070720001, 0.17459388557291625, 0.2780518073550838]
            elif (order.lower() == 'io'):
                flavor_ratios_pion \
                    = [0.3153279935881048, 0.34979331257813706, 0.334878693833758]
                flavor_ratios_muon \
                    = [0.19943830128706735, 0.424970818223672, 0.3755908804892607]
                flavor_ratios_neutron \
                    = [0.5471073781901797, 0.19943830128706735, 0.2534543205227526]
        if (nufit_version == 'nufit_v4.1'):
            if (order.lower() == 'no'):
                flavor_ratios_pion \
                    = [0.30438521493764115, 0.3536866096566163, 0.3419281754057424]
                flavor_ratios_muon \
                    = [0.18288455868137182, 0.4390876351442386, 0.3780278061743895]
                flavor_ratios_neutron \
                    = [0.5473865274501799, 0.18288455868137182, 0.2697289138684482]
            elif (order.lower() == 'io'):
                flavor_ratios_pion \
                    = [0.3202665854332483, 0.3448720154249937, 0.33486139914175794]
                flavor_ratios_muon \
                    = [0.20682472295246254, 0.41389566166125935, 0.3792796153862782]
                flavor_ratios_neutron \
                    = [0.5471503103948199, 0.20682472295246254, 0.24602496665271745]
        if (nufit_version == 'nufit_v5.0'):
            if (order.lower() == 'no'):
                flavor_ratios_pion \
                    = [0.29829637092300665, 0.3593068266653908, 0.34239680241160236]
                flavor_ratios_muon \
                    = [0.17144024550093245, 0.45324011724762003, 0.37531963725144746]
                flavor_ratios_neutron \
                    = [0.5520086217671551, 0.17144024550093245, 0.2765511327319123]
            elif (order.lower() == 'io'):
                flavor_ratios_pion \
                    = [0.3179873235848132, 0.34751722386393324, 0.33449545255125346]
                flavor_ratios_muon \
                    = [0.20107959603042933, 0.42073603778068525, 0.37818436618888535]
                flavor_ratios_neutron \
                    = [0.5518027786935809, 0.20107959603042933, 0.24711762527598982]

    if show_benchmark_astro_labels:
    
        lst_flavor_ratios_bf = [flavor_ratios_pion, flavor_ratios_muon, \
            flavor_ratios_neutron]

        lst_color = ['salmon', 'orange', 'limegreen']
        lst_markers = ['o', 's', '^']

        for i in range(len(lst_flavor_ratios_bf)):
            tax.scatter([lst_flavor_ratios_bf[i]], \
                marker=lst_markers[i], 
                color=lst_color[i], \
                edgecolor='k', #lst_color[i], 
                s=20, #40 
                linewidths=0.5, 
                zorder=4, #4
                alpha=1.)

        # Annotate
        tax.scatter([[0.20, 0.98, 0.]], \
            marker='o', color='salmon', edgecolor='k', \
            s=40, linewidths=0.5, zorder=2, alpha=1.)
        tax.ax.annotate( r'$\pi$ decay: $\left( 1:2:0 \right)_{\rm S}$', \
            xy = (0.72,0.83), xycoords='data', color='k', \
            fontsize=10, horizontalalignment='left' )
        tax.scatter([[0.235, 0.91, 0.]], \
            marker='s', color='orange', edgecolor='k', \
            s=40, linewidths=0.5, zorder=2, alpha=1.)
        tax.ax.annotate( r'$\mu$-damped: $\left( 0:1:0 \right)_{\rm S}$', \
            xy = (0.72,0.77), xycoords='data', color='k', \
            fontsize=10, horizontalalignment='left' )
        tax.scatter([[0.267, 0.843, 0.]], \
            marker='^', color='limegreen', edgecolor='k', \
            s=40, linewidths=0.5, zorder=2, alpha=1.)
        tax.ax.annotate( r'$n$ decay: $\left( 1:0:0 \right)_{\rm S}$', \
            xy = (0.715,0.71), xycoords='data', color='k', \
            fontsize=10, horizontalalignment='left' )

    elif show_benchmark_decay_labels:
    
        lst_flavor_content_bf = [flavor_cont_nu1, flavor_cont_nu2,
            flavor_cont_nu3]
        lst_color = ['salmon', 'limegreen', 'orange']
        lst_markers = ['o', '^', 's']

        for i in range(len(lst_flavor_content_bf)):
            tax.scatter([lst_flavor_content_bf[i]], \
                marker=lst_markers[i], 
                color=lst_color[i], \
                edgecolor='k', #lst_color[i], 
                s=20, #40 
                linewidths=0.5, 
                zorder=4, #4
                alpha=1.)

        # Annotate
        tax.scatter([[0.20, 0.98, 0.]], \
            marker='o', color='salmon', edgecolor='k', \
            s=40, linewidths=0.5, zorder=2, alpha=1.)
        tax.ax.annotate( r'$\nu_1$', \
            xy = (0.72,0.83), xycoords='data', color='k', \
            fontsize=10, horizontalalignment='left' )
        tax.scatter([[0.235, 0.91, 0.]], \
            marker='s', color='orange', edgecolor='k', \
            s=40, linewidths=0.5, zorder=2, alpha=1.)
        tax.ax.annotate( r'$\nu_2$', \
            xy = (0.72,0.77), xycoords='data', color='k', \
            fontsize=10, horizontalalignment='left' )
        tax.scatter([[0.267, 0.843, 0.]], \
            marker='^', color='limegreen', edgecolor='k', \
            s=40, linewidths=0.5, zorder=2, alpha=1.)
        tax.ax.annotate( r'$\nu_3$', \
            xy = (0.715,0.71), xycoords='data', color='k', \
            fontsize=10, horizontalalignment='left' )

    tax.clear_matplotlib_ticks()
    tax._redraw_labels()
    ternary.plt.tight_layout()

    ternary.plt.savefig(fname+'.'+output_format, dpi=300)
    print("Saved "+fname+'.'+output_format)

    plt.close()

    return


def plot_ternary_flavor_ratios_fixed_fS_separate(fname, plot_label, fS='120', 
    order='no', 
    theta23='upper', output_format='png'):

    # Open the plot and format it
    mpl.rcParams['xtick.labelsize']=26
    mpl.rcParams['ytick.labelsize']=26
    mpl.rcParams['legend.fontsize']=23
    mpl.rcParams['legend.borderpad']=0.4
    mpl.rcParams['axes.labelpad']=10
    mpl.rcParams['ps.fonttype']=42
    mpl.rcParams['pdf.fonttype']=42

    colors = ['#2a4d69', '#4b86b4', '#adcbe3', '#D3E4F4', '#63ace5']

    # Boundary and gridlines
    fig, tax = ternary.figure(scale=1.0)
    tax.ax.axis("off")
    fig.set_facecolor('w')

    # Draw boundary and gridlines
    tax.boundary(linewidth=1.0)
    tax.gridlines(color="gray", multiple=0.1, linewidth=0.5, ls='-', alpha=0.5)

    fontsize = 15
    label_left = r'Fraction of $\nu_\tau$, $f_{\tau, \oplus}$'
    label_right = r'Fraction of $\nu_\mu$, $f_{\mu, \oplus}$'
    label_bottom = r'Fraction of $\nu_e$, $f_{e, \oplus}$'
    tax.left_axis_label(label_left, fontsize=fontsize, offset=0.16)
    tax.right_axis_label(label_right, fontsize=fontsize, offset=0.16)
    tax.bottom_axis_label(label_bottom, fontsize=fontsize, offset=0.08)

    # Set ticks blr
    tax.ticks(axis='blr', linewidth=0.5, multiple=0.1, offset=0.022, \
        clockwise=False, tick_formats="%.1f")

    ###########################################################################
    # LIV-like region
    ###########################################################################

    # x = [0.0, 1.0, 0.0]
    # y = [0.0, 0.0, 1.0]
    # z = [1.0-x[i]-y[i] for i in range(len(x))]
    # x_rot = [ 0.5 * (2.*x[i]+y[i]) / (x[i]+y[i]+z[i]) for i in range(len(x)) ]
    # y_rot = [ sqrt(3.)/2 * (y[i]) / (x[i]+y[i]+z[i]) for i in range(len(x)) ]
    # tax.ax.fill(x_rot, y_rot, color=lst_colors[4], alpha=0.4)


    ###########################################################################
    # Decay-like region
    ###########################################################################

    # filename_in = os.getcwd()+"/flavor_ratios_earth_np_decay_3s.dat"
    # contour_sm = Read_Data_File(filename_in)
    # x = contour_sm[1]
    # y = contour_sm[0]
    # z = [1.0-x[i]-y[i] for i in range(len(x))]
    # x_rot = [ 0.5 * (2.*x[i]+y[i]) / (x[i]+y[i]+z[i]) for i in range(len(x)) ]
    # y_rot = [ sqrt(3.)/2 * (y[i]) / (x[i]+y[i]+z[i]) for i in range(len(x)) ]
    # tax.ax.fill(x_rot, y_rot, color=lst_colors[1])


    ###########################################################################
    # Standard Model region
    ###########################################################################

    """
    x = [0.0, 0.1, 0.1 ,0.0]
    y = [0.0, 0.0, 0.1, 0.1]
    z = [1.0-x[i]-y[i] for i in range(len(x))]
    tax.ax.fill(x,y)
    x_rot = [ 0.5 * (2.*x[i]+y[i]) / (x[i]+y[i]+z[i]) for i in range(len(x)) ]
    y_rot = [ sqrt(3.)/2 * (y[i]) / (x[i]+y[i]+z[i]) for i in range(len(x)) ]
    tax.ax.fill(x_rot, y_rot)
    """

    # filename_in = os.getcwd()+"/flavor_ratios_earth_full_var_std_NH_3s_outline.dat"
    # contour_sm = Read_Data_File(filename_in)
    # x = contour_sm[1]
    # y = contour_sm[0]
    # z = [1.0-x[i]-y[i] for i in range(len(x))]
    # x_rot = [ 0.5 * (2.*x[i]+y[i]) / (x[i]+y[i]+z[i]) for i in range(len(x)) ]
    # y_rot = [ sqrt(3.)/2 * (y[i]) / (x[i]+y[i]+z[i]) for i in range(len(x)) ]
    # tax.ax.fill(x_rot, y_rot, color=lst_colors[0])


    ###########################################################################
    # Standard Model regions
    ###########################################################################

    with open(fname+'.json') as file:
        data = json.load(file)

    for j, cl in enumerate(list(data.keys())[::-1]):
        fmu = data[cl]['f_e_Earth']
        fe = data[cl]['f_mu_Earth']
        ftau = [1.0-fe[i]-fmu[i] for i in range(len(fe))]
        fe_rot = [ 0.5 * (2.*fe[i]+fmu[i]) / (fe[i]+fmu[i]+ftau[i]) 
            for i in range(len(fe)) ]
        fmu_rot = [ sqrt(3.)/2 * (fmu[i]) / (fe[i]+fmu[i]+ftau[i]) 
            for i in range(len(fe)) ]
        ftau_rot = [1.0-fe_rot[i]-fmu_rot[i] for i in range(len(fe))]
        # tax.ax.fill(fe_rot, fmu_rot, color=colors[j])
        tax.ax.fill(fe_rot, fmu_rot, color=colors[j])


    ###########################################################################
    # Draw IceCube 2015 contours
    ###########################################################################

    """
    # lst_filenames_ic =  [ \
    #     os.getcwd()+"/flavor_ratios_earth_ic_2015_bf.dat", \
    #     os.getcwd()+"/flavor_ratios_earth_ic_2015_1s.dat", \
    #     # os.getcwd()+"/flavor_ratios_earth_ic_2015_2s_a.dat", \
    #     # os.getcwd()+"/flavor_ratios_earth_ic_2015_2s_b.dat", \
    #     os.getcwd()+"/flavor_ratios_earth_ic_2015_3s_a.dat", \
    #     os.getcwd()+"/flavor_ratios_earth_ic_2015_3s_b.dat"]

    lst_filenames_ic =  [ \
        os.getcwd()+"/icgen2_icrc2017_pion_68cl.dat",
        os.getcwd()+"/icgen2_icrc2017_pion_90cl.dat",
        os.getcwd()+"/icgen2_icrc2017_muon_68cl.dat",
        os.getcwd()+"/icgen2_icrc2017_muon_90cl.dat"]

    lst_color = ['gold', 'gold', 'w', 'w'] #'y'
    lst_ls = ['--', '-', '--', '-']

    # Each row: f_eE, f_mE
    for file_index in range(len(lst_filenames_ic)):

        # Read data file with IceCube flavor ratios
        lst_ratios_earth = Read_Data_File(lst_filenames_ic[file_index])
        # Recast as [f_eE, f_mE, f_tE]
        lst_ratios_earth = \
            [[lst_ratios_earth[1][i], lst_ratios_earth[2][i],\
            lst_ratios_earth[0][i]] \
            for i in range(len(lst_ratios_earth[0]))]

        tax.plot(lst_ratios_earth, color=lst_color[file_index],
                ls=lst_ls[file_index], linewidth=2.5, zorder=1)

    """

    # IceCube, 8 yr
    with open('icecube_8yr.json') as json_file:
        data = json.load(json_file)
        fraction_nue = data['data']['nue_fraction']
        fraction_numu = data['data']['numu_fraction']
        fraction_nutau = \
            [1.0-fraction_nue[i]-fraction_numu[i] \
            for i in range(len(fraction_nue))]
        confidence_level = data['data']['confidence_level']
        cs = tax.ax.contour(fraction_nue, fraction_numu, confidence_level,
                levels=[68, ], colors=None, linestyles=None, alpha=0.0)
        vertices = cs.collections[0].get_paths()[0].vertices
        x = vertices[:,0]
        y = vertices[:,1]
        z = [1.0-x[i]-y[i] for i in range(len(vertices))]
        ratios_earth = [[x[i], y[i], z[i]] for i in range(len(x))]
        tax.plot(ratios_earth, color='tomato',
                ls='-', linewidth=1.5, zorder=1)

    # IceCube, 15 yr
    with open('icecube_15yr.json') as json_file:
        data = json.load(json_file)
        fraction_nue = data['data']['nue_fraction']
        fraction_numu = data['data']['numu_fraction']
        fraction_nutau = \
            [1.0-fraction_nue[i]-fraction_numu[i] \
            for i in range(len(fraction_nue))]
        confidence_level = data['data']['confidence_level']
        cs = tax.ax.contour(fraction_nue, fraction_numu, confidence_level,
                levels=[68, ], colors=None, linestyles=None, alpha=0.0)
        vertices = cs.collections[0].get_paths()[0].vertices
        x = vertices[:,0]
        y = vertices[:,1]
        z = [1.0-x[i]-y[i] for i in range(len(vertices))]
        ratios_earth = [[x[i], y[i], z[i]] for i in range(len(x))]
        tax.plot(ratios_earth, color='tan',
                ls='-', linewidth=1.5, zorder=1)

    # IceCube-Gen2 10-yr + IceCube 15-yr
    with open('icecube+gen2-inice.json') as json_file:
        data = json.load(json_file)
        fraction_nue = data['data']['nue_fraction']
        fraction_numu = data['data']['numu_fraction']
        fraction_nutau = \
            [1.0-fraction_nue[i]-fraction_numu[i] \
            for i in range(len(fraction_nue))]
        confidence_level = data['data']['confidence_level']
        cs = tax.ax.contour(fraction_nue, fraction_numu, confidence_level,
                levels=[68, ], colors=None, linestyles=None, alpha=0.0)
        vertices = cs.collections[0].get_paths()[0].vertices
        x = vertices[:,0]
        y = vertices[:,1]
        z = [1.0-x[i]-y[i] for i in range(len(vertices))]
        ratios_earth = [[x[i], y[i], z[i]] for i in range(len(x))]
        tax.plot(ratios_earth, color='gold',
                ls='-', linewidth=1.5, zorder=1)

    ###########################################################################
    # Annotations
    ###########################################################################

    tax.ax.annotate( plot_label, \
        xy = (-0.03, 0.845), xycoords='data', color='k', \
        fontsize=10, horizontalalignment='left' )

    tax.ax.scatter([-0.01], [0.80], \
        marker='s', color=colors[3], edgecolor='k', \
        s=60, linewidths=0.1, zorder=2, alpha=1.)
    tax.ax.annotate( r'68\%~C.L.', \
        xy = (0.02,0.787), xycoords='data', color='k', \
        fontsize=10, horizontalalignment='left' )

    tax.ax.scatter([-0.01], [0.75], \
        marker='s', color=colors[2], edgecolor=colors[2], \
        s=60, linewidths=0.1, zorder=2, alpha=1.)
    tax.ax.annotate( r'90\%~C.L.', \
        xy = (0.02,0.737), xycoords='data', color='k', \
        fontsize=10, horizontalalignment='left' )

    tax.ax.scatter([-0.01], [0.70], \
        marker='s', color=colors[1], edgecolor=colors[1], \
        s=60, linewidths=0.1, zorder=2, alpha=1.)
    tax.ax.annotate( r'95\%~C.L.', \
        xy = (0.02,0.687), xycoords='data', color='k', \
        fontsize=10, horizontalalignment='left' )

    tax.ax.scatter([-0.01], [0.65], \
        marker='s', color=colors[0], edgecolor=colors[0], \
        s=60, linewidths=0.1, zorder=2, alpha=1.)
    tax.ax.annotate( r'99.7\%~C.L.', \
        xy = (0.02,0.637), xycoords='data', color='k', \
        fontsize=10, horizontalalignment='left' )

    # tax.ax.scatter([-0.01], [0.872], \
    #     marker='s', color=colors[0], edgecolor=colors[0], \
    #     s=60, linewidths=0.1, zorder=2, alpha=1.)
    # tax.ax.annotate( r'Standard Model ($\nu$SM)', \
    #     xy = (0.02,0.86), xycoords='data', color='k', \
    #     fontsize=10, horizontalalignment='left' )

    # tax.ax.scatter([-0.01], [0.79], \
    #     marker='s', color=colors[1], edgecolor=colors[1], \
    #     s=60, linewidths=0.1, zorder=2, alpha=1.)
    # tax.ax.annotate( r'New physics:', \
    #     xy = (0.02,0.80), xycoords='data', color='k', \
    #     fontsize=10, horizontalalignment='left' )
    # tax.ax.annotate( r'$\nu$ decay-like', \
    #     xy = (0.02,0.755), xycoords='data', color='k', \
    #     fontsize=10, horizontalalignment='left' )

    # tax.ax.scatter([-0.01], [0.685], \
    #     marker='s', color=colors[4], edgecolor=colors[4], \
    #     s=60, linewidths=0.1, zorder=2, alpha=0.4)
    # tax.ax.annotate( r'New physics:', \
    #     xy = (0.02,0.695), xycoords='data', color='k', \
    #     fontsize=10, horizontalalignment='left' )
    # tax.ax.annotate( r'Lorentz violation', \
    #     xy = (0.02,0.65), xycoords='data', color='k', \
    #     fontsize=10, horizontalalignment='left' )

    tax.plot([[0.03, 0.14], [0.07, 0.14]], color='tomato',
                ls='-', linewidth=1.5, zorder=1)
    tax.ax.annotate( r'IceCube 8 yr (68\%)', \
        xy = (0.15,0.11), xycoords='data', color='k', \
        fontsize=10, horizontalalignment='left' )
    tax.plot([[0.05, 0.095], [0.09, 0.095]], color='tan',
                ls='-', linewidth=1.5, zorder=1)
    tax.ax.annotate( r'IceCube 15 yr (68\%)', \
        xy = (0.15,0.07), xycoords='data', color='k', \
        fontsize=10, horizontalalignment='left' )
    tax.plot([[0.073, 0.05], [0.113, 0.05]], color='gold',
                ls='-', linewidth=1.5, zorder=1)
    tax.ax.annotate( r'IceCube 15 yr + Gen2 10 yr (68\%)', \
        xy = (0.15,0.03), xycoords='data', color='k', \
        fontsize=10, horizontalalignment='left' )


    # tax.ax.annotate( r'$1\sigma$', \
    #     xy = (0.58,0.56), xycoords='data', color='w', rotation=15., \
    #     fontsize=14, horizontalalignment='left',
    #     path_effects=[path_effects.withSimplePatchShadow(offset=(1,-1),
    #                 shadow_rgbFace='k', alpha=0.8)] )
    # tax.ax.annotate( r'$3\sigma$', \
    #     xy = (0.51,0.70), xycoords='data', color='w', rotation=15., \
    #     fontsize=14, horizontalalignment='left',
    #     path_effects=[path_effects.withSimplePatchShadow(offset=(1,-1),
    #                 shadow_rgbFace='k', alpha=0.8)] )


    # text = plt.text(0.092, 0.10, 'IceCube 2015', rotation=5.,
    #             fontdict={'color': 'w'},
    #             path_effects=[path_effects.withSimplePatchShadow(offset=(1,-1),
    #                 shadow_rgbFace='k', alpha=0.8)])


    ####################################################################
    # Add regions of selected compositions
    ####################################################################

    if (order.lower() == 'no'):
        if (theta23 == 'upper'):
            flavor_ratios_pion \
                = [0.2982716010150718, 0.359320837204981, 0.34240756177994724]
            flavor_ratios_muon \
                = [0.17128846521772986, 0.4533370231986065, 0.37537451158366353]
            flavor_ratios_neutron \
                = [0.5522378726097557, 0.17128846521772986, 0.2764736621725147]
        elif (theta23 == 'lower'):
            flavor_ratios_pion \
                = [0.32800656074697265, 0.34529257725052065, 0.3267008620025065]
            flavor_ratios_muon \
                = [0.2158909048155812, 0.40999341346799034, 0.3741156817164281]
            flavor_ratios_neutron \
                = [0.5522378726097557, 0.2158909048155812, 0.23187122257466325]
        elif (theta23 == 'max'):
            flavor_ratios_pion \
                = [0.31656024191908616, 0.3474455438823001, 0.3359942141986137]
            flavor_ratios_muon \
                = [0.19872142657375141, 0.4218076025365745, 0.37947097088967396]
            flavor_ratios_neutron \
                = [0.5522378726097557, 0.19872142657375141, 0.24904070081649313]
    elif (order.lower() == 'io'):
        if (theta23 == 'upper'):
            flavor_ratios_pion \
                = [0.3180471687709684, 0.3474505778780882, 0.3345022533509433]
            flavor_ratios_muon \
                = [0.2010801776404412, 0.4206357779969117, 0.378284044362647]
            flavor_ratios_neutron \
                = [0.5519811510320229, 0.2010801776404412, 0.24693867132753594]
        elif (theta23 == 'lower'):
            flavor_ratios_pion \
                = [0.34834738484344396, 0.3277980374633891, 0.32385457769316683]
            flavor_ratios_muon \
                = [0.2465305017491545, 0.36843180532050646, 0.38503769293033896]
            flavor_ratios_neutron \
                = [0.5519811510320229, 0.2465305017491545, 0.20148834721882264]
        elif (theta23 == 'max'):
            flavor_ratios_pion \
                = [0.33698743757165983, 0.3317785331561326, 0.33123402927220735]
            flavor_ratios_muon \
                = [0.22949058084147833, 0.3829225093134597, 0.3875869098450617]
            flavor_ratios_neutron \
                = [0.5519811510320229, 0.22949058084147833, 0.21852826812649873]

    lst_flavor_ratios_bf = [flavor_ratios_pion, flavor_ratios_muon, \
                            flavor_ratios_neutron]
    lst_color = ['salmon', 'orange', 'limegreen']
    lst_markers = ['o', 's', '^']

    for i in range(len(lst_flavor_ratios_bf)):
        # flavor_ratios = [lst_flavor_ratios_bf[i][0], lst_flavor_ratios_bf[i][2],
            # lst_flavor_ratios_bf[i][1]] 
        tax.scatter([lst_flavor_ratios_bf[i]], \
        # tax.scatter([flavor_ratios], \
            marker=lst_markers[i], color=lst_color[i], \
            edgecolor=lst_color[i], s=40, linewidths=0.2, zorder=4, \
            alpha=1.)

    # Annotate
    tax.scatter([[0.20, 0.98, 0.]], \
        marker='o', color='salmon', edgecolor='salmon', \
        s=40, linewidths=0.1, zorder=2, alpha=1.)
    tax.ax.annotate( r'$\pi$ decay: $\left( 1:2:0 \right)_{\rm S}$', \
        xy = (0.72,0.83), xycoords='data', color='k', \
        fontsize=10, horizontalalignment='left' )
    tax.scatter([[0.235, 0.91, 0.]], \
        marker='s', color='orange', edgecolor='orange', \
        s=40, linewidths=0.1, zorder=2, alpha=1.)
    tax.ax.annotate( r'$\mu$-damped: $\left( 0:1:0 \right)_{\rm S}$', \
        xy = (0.72,0.77), xycoords='data', color='k', \
        fontsize=10, horizontalalignment='left' )
    tax.scatter([[0.27, 0.843, 0.]], \
        marker='^', color='limegreen', edgecolor='limegreen', \
        s=40, linewidths=0.1, zorder=2, alpha=1.)
    tax.ax.annotate( r'$n$ decay: $\left( 1:0:0 \right)_{\rm S}$', \
        xy = (0.72,0.71), xycoords='data', color='k', \
        fontsize=10, horizontalalignment='left' )


    tax.clear_matplotlib_ticks()
    tax._redraw_labels()
    ternary.plt.tight_layout()

    ternary.plt.savefig(fname, dpi=200)
    print("Saved "+fname)

    plt.close()

    return


def plot_ternary_flavor_ratios_fixed_fS_combined_v0(fname_in, fname_out, 
    plot_label, order='no', theta23='upper', output_format='png'):

    # Open the plot and format it
    mpl.rcParams['xtick.labelsize']=26
    mpl.rcParams['ytick.labelsize']=26
    mpl.rcParams['legend.fontsize']=23
    mpl.rcParams['legend.borderpad']=0.4
    mpl.rcParams['axes.labelpad']=10
    mpl.rcParams['ps.fonttype']=42
    mpl.rcParams['pdf.fonttype']=42

    colors_gen = ['#2a4d69', '#4b86b4', '#adcbe3', '#D3E4F4', '#63ace5']
    colors_pi = ['darkred', 'darksalmon', 'salmon', 'lightsalmon']
    colors = [colors_gen, colors_gen, colors_gen]

    # Boundary and gridlines
    fig, tax = ternary.figure(scale=1.0)
    tax.ax.axis("off")
    fig.set_facecolor('w')

    # Draw boundary and gridlines
    tax.boundary(linewidth=1.0)
    tax.gridlines(color="gray", multiple=0.1, linewidth=0.5, ls='-', alpha=0.5)

    fontsize = 15
    label_left = r'Fraction of $\nu_\tau$, $f_{\tau, \oplus}$'
    label_right = r'Fraction of $\nu_\mu$, $f_{\mu, \oplus}$'
    label_bottom = r'Fraction of $\nu_e$, $f_{e, \oplus}$'
    tax.left_axis_label(label_left, fontsize=fontsize, offset=0.16)
    tax.right_axis_label(label_right, fontsize=fontsize, offset=0.16)
    tax.bottom_axis_label(label_bottom, fontsize=fontsize, offset=0.08)

    # Set ticks blr
    tax.ticks(axis='blr', linewidth=0.5, multiple=0.1, offset=0.022, \
        clockwise=False, tick_formats="%.1f")

    ###########################################################################
    # LIV-like region
    ###########################################################################

    # x = [0.0, 1.0, 0.0]
    # y = [0.0, 0.0, 1.0]
    # z = [1.0-x[i]-y[i] for i in range(len(x))]
    # x_rot = [ 0.5 * (2.*x[i]+y[i]) / (x[i]+y[i]+z[i]) for i in range(len(x)) ]
    # y_rot = [ sqrt(3.)/2 * (y[i]) / (x[i]+y[i]+z[i]) for i in range(len(x)) ]
    # tax.ax.fill(x_rot, y_rot, color=lst_colors[4], alpha=0.4)


    ###########################################################################
    # Decay-like region
    ###########################################################################

    # filename_in = os.getcwd()+"/flavor_ratios_earth_np_decay_3s.dat"
    # contour_sm = Read_Data_File(filename_in)
    # x = contour_sm[1]
    # y = contour_sm[0]
    # z = [1.0-x[i]-y[i] for i in range(len(x))]
    # x_rot = [ 0.5 * (2.*x[i]+y[i]) / (x[i]+y[i]+z[i]) for i in range(len(x)) ]
    # y_rot = [ sqrt(3.)/2 * (y[i]) / (x[i]+y[i]+z[i]) for i in range(len(x)) ]
    # tax.ax.fill(x_rot, y_rot, color=lst_colors[1])


    ###########################################################################
    # Standard Model region
    ###########################################################################

    """
    x = [0.0, 0.1, 0.1 ,0.0]
    y = [0.0, 0.0, 0.1, 0.1]
    z = [1.0-x[i]-y[i] for i in range(len(x))]
    tax.ax.fill(x,y)
    x_rot = [ 0.5 * (2.*x[i]+y[i]) / (x[i]+y[i]+z[i]) for i in range(len(x)) ]
    y_rot = [ sqrt(3.)/2 * (y[i]) / (x[i]+y[i]+z[i]) for i in range(len(x)) ]
    tax.ax.fill(x_rot, y_rot)
    """

    # filename_in = os.getcwd()+"/flavor_ratios_earth_full_var_std_NH_3s_outline.dat"
    # contour_sm = Read_Data_File(filename_in)
    # x = contour_sm[1]
    # y = contour_sm[0]
    # z = [1.0-x[i]-y[i] for i in range(len(x))]
    # x_rot = [ 0.5 * (2.*x[i]+y[i]) / (x[i]+y[i]+z[i]) for i in range(len(x)) ]
    # y_rot = [ sqrt(3.)/2 * (y[i]) / (x[i]+y[i]+z[i]) for i in range(len(x)) ]
    # tax.ax.fill(x_rot, y_rot, color=lst_colors[0])


    ###########################################################################
    # Standard Model regions
    ###########################################################################

    for k, fname in enumerate(fname_in):
        with open(fname+'.json') as file:
            data = json.load(file)
        for j, cl in enumerate(list(data.keys())[::-1]):
            fmu = data[cl]['f_e_Earth']
            fe = data[cl]['f_mu_Earth']
            ftau = [1.0-fe[i]-fmu[i] for i in range(len(fe))]
            fe_rot = [ 0.5 * (2.*fe[i]+fmu[i]) / (fe[i]+fmu[i]+ftau[i]) 
                for i in range(len(fe)) ]
            fmu_rot = [ sqrt(3.)/2 * (fmu[i]) / (fe[i]+fmu[i]+ftau[i]) 
                for i in range(len(fe)) ]
            ftau_rot = [1.0-fe_rot[i]-fmu_rot[i] for i in range(len(fe))]
            # tax.ax.fill(fe_rot, fmu_rot, color=colors[j])
            tax.ax.fill(fe_rot, fmu_rot, color=colors[k][j])


    ###########################################################################
    # Draw IceCube 2015 contours
    ###########################################################################

    """
    # lst_filenames_ic =  [ \
    #     os.getcwd()+"/flavor_ratios_earth_ic_2015_bf.dat", \
    #     os.getcwd()+"/flavor_ratios_earth_ic_2015_1s.dat", \
    #     # os.getcwd()+"/flavor_ratios_earth_ic_2015_2s_a.dat", \
    #     # os.getcwd()+"/flavor_ratios_earth_ic_2015_2s_b.dat", \
    #     os.getcwd()+"/flavor_ratios_earth_ic_2015_3s_a.dat", \
    #     os.getcwd()+"/flavor_ratios_earth_ic_2015_3s_b.dat"]

    lst_filenames_ic =  [ \
        os.getcwd()+"/icgen2_icrc2017_pion_68cl.dat",
        os.getcwd()+"/icgen2_icrc2017_pion_90cl.dat",
        os.getcwd()+"/icgen2_icrc2017_muon_68cl.dat",
        os.getcwd()+"/icgen2_icrc2017_muon_90cl.dat"]

    lst_color = ['gold', 'gold', 'w', 'w'] #'y'
    lst_ls = ['--', '-', '--', '-']

    # Each row: f_eE, f_mE
    for file_index in range(len(lst_filenames_ic)):

        # Read data file with IceCube flavor ratios
        lst_ratios_earth = Read_Data_File(lst_filenames_ic[file_index])
        # Recast as [f_eE, f_mE, f_tE]
        lst_ratios_earth = \
            [[lst_ratios_earth[1][i], lst_ratios_earth[2][i],\
            lst_ratios_earth[0][i]] \
            for i in range(len(lst_ratios_earth[0]))]

        tax.plot(lst_ratios_earth, color=lst_color[file_index],
                ls=lst_ls[file_index], linewidth=2.5, zorder=1)

    """

    # IceCube, 8 yr
    with open('icecube_8yr.json') as json_file:
        data = json.load(json_file)
        fraction_nue = data['data']['nue_fraction']
        fraction_numu = data['data']['numu_fraction']
        fraction_nutau = \
            [1.0-fraction_nue[i]-fraction_numu[i] \
            for i in range(len(fraction_nue))]
        confidence_level = data['data']['confidence_level']
        cs = tax.ax.contour(fraction_nue, fraction_numu, confidence_level,
                levels=[68, ], colors=None, linestyles=None, alpha=0.0)
        vertices = cs.collections[0].get_paths()[0].vertices
        x = vertices[:,0]
        y = vertices[:,1]
        z = [1.0-x[i]-y[i] for i in range(len(vertices))]
        ratios_earth = [[x[i], y[i], z[i]] for i in range(len(x))]
        tax.plot(ratios_earth, color='tomato',
                ls='-', linewidth=1.5, zorder=1)

    # IceCube, 15 yr
    with open('icecube_15yr.json') as json_file:
        data = json.load(json_file)
        fraction_nue = data['data']['nue_fraction']
        fraction_numu = data['data']['numu_fraction']
        fraction_nutau = \
            [1.0-fraction_nue[i]-fraction_numu[i] \
            for i in range(len(fraction_nue))]
        confidence_level = data['data']['confidence_level']
        cs = tax.ax.contour(fraction_nue, fraction_numu, confidence_level,
                levels=[68, ], colors=None, linestyles=None, alpha=0.0)
        vertices = cs.collections[0].get_paths()[0].vertices
        x = vertices[:,0]
        y = vertices[:,1]
        z = [1.0-x[i]-y[i] for i in range(len(vertices))]
        ratios_earth = [[x[i], y[i], z[i]] for i in range(len(x))]
        tax.plot(ratios_earth, color='tan',
                ls='-', linewidth=1.5, zorder=1)

    # IceCube-Gen2 10-yr + IceCube 15-yr
    with open('icecube+gen2-inice.json') as json_file:
        data = json.load(json_file)
        fraction_nue = data['data']['nue_fraction']
        fraction_numu = data['data']['numu_fraction']
        fraction_nutau = \
            [1.0-fraction_nue[i]-fraction_numu[i] \
            for i in range(len(fraction_nue))]
        confidence_level = data['data']['confidence_level']
        cs = tax.ax.contour(fraction_nue, fraction_numu, confidence_level,
                levels=[68, ], colors=None, linestyles=None, alpha=0.0)
        vertices = cs.collections[0].get_paths()[0].vertices
        x = vertices[:,0]
        y = vertices[:,1]
        z = [1.0-x[i]-y[i] for i in range(len(vertices))]
        ratios_earth = [[x[i], y[i], z[i]] for i in range(len(x))]
        tax.plot(ratios_earth, color='gold',
                ls='-', linewidth=1.5, zorder=1)

    ###########################################################################
    # Annotations
    ###########################################################################

    tax.ax.annotate( plot_label, \
        xy = (-0.03, 0.845), xycoords='data', color='k', \
        fontsize=10, horizontalalignment='left' )

    tax.ax.scatter([-0.01], [0.80], \
        marker='s', color=colors_gen[3], edgecolor='k', \
        s=60, linewidths=0.1, zorder=2, alpha=1.)
    tax.ax.annotate( r'68\%~C.L.', \
        xy = (0.02,0.787), xycoords='data', color='k', \
        fontsize=10, horizontalalignment='left' )

    tax.ax.scatter([-0.01], [0.75], \
        marker='s', color=colors_gen[2], edgecolor=colors_gen[2], \
        s=60, linewidths=0.1, zorder=2, alpha=1.)
    tax.ax.annotate( r'90\%~C.L.', \
        xy = (0.02,0.737), xycoords='data', color='k', \
        fontsize=10, horizontalalignment='left' )

    tax.ax.scatter([-0.01], [0.70], \
        marker='s', color=colors_gen[1], edgecolor=colors_gen[1], \
        s=60, linewidths=0.1, zorder=2, alpha=1.)
    tax.ax.annotate( r'95\%~C.L.', \
        xy = (0.02,0.687), xycoords='data', color='k', \
        fontsize=10, horizontalalignment='left' )

    tax.ax.scatter([-0.01], [0.65], \
        marker='s', color=colors_gen[0], edgecolor=colors_gen[0], \
        s=60, linewidths=0.1, zorder=2, alpha=1.)
    tax.ax.annotate( r'99.7\%~C.L.', \
        xy = (0.02,0.637), xycoords='data', color='k', \
        fontsize=10, horizontalalignment='left' )

    # tax.ax.scatter([-0.01], [0.872], \
    #     marker='s', color=colors[0], edgecolor=colors[0], \
    #     s=60, linewidths=0.1, zorder=2, alpha=1.)
    # tax.ax.annotate( r'Standard Model ($\nu$SM)', \
    #     xy = (0.02,0.86), xycoords='data', color='k', \
    #     fontsize=10, horizontalalignment='left' )

    # tax.ax.scatter([-0.01], [0.79], \
    #     marker='s', color=colors[1], edgecolor=colors[1], \
    #     s=60, linewidths=0.1, zorder=2, alpha=1.)
    # tax.ax.annotate( r'New physics:', \
    #     xy = (0.02,0.80), xycoords='data', color='k', \
    #     fontsize=10, horizontalalignment='left' )
    # tax.ax.annotate( r'$\nu$ decay-like', \
    #     xy = (0.02,0.755), xycoords='data', color='k', \
    #     fontsize=10, horizontalalignment='left' )

    # tax.ax.scatter([-0.01], [0.685], \
    #     marker='s', color=colors[4], edgecolor=colors[4], \
    #     s=60, linewidths=0.1, zorder=2, alpha=0.4)
    # tax.ax.annotate( r'New physics:', \
    #     xy = (0.02,0.695), xycoords='data', color='k', \
    #     fontsize=10, horizontalalignment='left' )
    # tax.ax.annotate( r'Lorentz violation', \
    #     xy = (0.02,0.65), xycoords='data', color='k', \
    #     fontsize=10, horizontalalignment='left' )

    tax.plot([[0.03, 0.14], [0.07, 0.14]], color='tomato',
                ls='-', linewidth=1.5, zorder=1)
    tax.ax.annotate( r'IceCube 8 yr (68\%)', \
        xy = (0.15,0.11), xycoords='data', color='k', \
        fontsize=10, horizontalalignment='left' )
    tax.plot([[0.05, 0.095], [0.09, 0.095]], color='tan',
                ls='-', linewidth=1.5, zorder=1)
    tax.ax.annotate( r'IceCube 15 yr (68\%)', \
        xy = (0.15,0.07), xycoords='data', color='k', \
        fontsize=10, horizontalalignment='left' )
    tax.plot([[0.073, 0.05], [0.113, 0.05]], color='gold',
                ls='-', linewidth=1.5, zorder=1)
    tax.ax.annotate( r'IceCube 15 yr + Gen2 10 yr (68\%)', \
        xy = (0.15,0.03), xycoords='data', color='k', \
        fontsize=10, horizontalalignment='left' )


    # tax.ax.annotate( r'$1\sigma$', \
    #     xy = (0.58,0.56), xycoords='data', color='w', rotation=15., \
    #     fontsize=14, horizontalalignment='left',
    #     path_effects=[path_effects.withSimplePatchShadow(offset=(1,-1),
    #                 shadow_rgbFace='k', alpha=0.8)] )
    # tax.ax.annotate( r'$3\sigma$', \
    #     xy = (0.51,0.70), xycoords='data', color='w', rotation=15., \
    #     fontsize=14, horizontalalignment='left',
    #     path_effects=[path_effects.withSimplePatchShadow(offset=(1,-1),
    #                 shadow_rgbFace='k', alpha=0.8)] )


    # text = plt.text(0.092, 0.10, 'IceCube 2015', rotation=5.,
    #             fontdict={'color': 'w'},
    #             path_effects=[path_effects.withSimplePatchShadow(offset=(1,-1),
    #                 shadow_rgbFace='k', alpha=0.8)])


    ####################################################################
    # Add regions of selected compositions
    ####################################################################

    if (order.lower() == 'no'):
        if (theta23 == 'upper'):
            flavor_ratios_pion \
                = [0.2982716010150718, 0.359320837204981, 0.34240756177994724]
            flavor_ratios_muon \
                = [0.17128846521772986, 0.4533370231986065, 0.37537451158366353]
            flavor_ratios_neutron \
                = [0.5522378726097557, 0.17128846521772986, 0.2764736621725147]
        elif (theta23 == 'lower'):
            flavor_ratios_pion \
                = [0.32800656074697265, 0.34529257725052065, 0.3267008620025065]
            flavor_ratios_muon \
                = [0.2158909048155812, 0.40999341346799034, 0.3741156817164281]
            flavor_ratios_neutron \
                = [0.5522378726097557, 0.2158909048155812, 0.23187122257466325]
        elif (theta23 == 'max'):
            flavor_ratios_pion \
                = [0.31656024191908616, 0.3474455438823001, 0.3359942141986137]
            flavor_ratios_muon \
                = [0.19872142657375141, 0.4218076025365745, 0.37947097088967396]
            flavor_ratios_neutron \
                = [0.5522378726097557, 0.19872142657375141, 0.24904070081649313]
    elif (order.lower() == 'io'):
        if (theta23 == 'upper'):
            flavor_ratios_pion \
                = [0.3180471687709684, 0.3474505778780882, 0.3345022533509433]
            flavor_ratios_muon \
                = [0.2010801776404412, 0.4206357779969117, 0.378284044362647]
            flavor_ratios_neutron \
                = [0.5519811510320229, 0.2010801776404412, 0.24693867132753594]
        elif (theta23 == 'lower'):
            flavor_ratios_pion \
                = [0.34834738484344396, 0.3277980374633891, 0.32385457769316683]
            flavor_ratios_muon \
                = [0.2465305017491545, 0.36843180532050646, 0.38503769293033896]
            flavor_ratios_neutron \
                = [0.5519811510320229, 0.2465305017491545, 0.20148834721882264]
        elif (theta23 == 'max'):
            flavor_ratios_pion \
                = [0.33698743757165983, 0.3317785331561326, 0.33123402927220735]
            flavor_ratios_muon \
                = [0.22949058084147833, 0.3829225093134597, 0.3875869098450617]
            flavor_ratios_neutron \
                = [0.5519811510320229, 0.22949058084147833, 0.21852826812649873]

    lst_flavor_ratios_bf = [flavor_ratios_pion, flavor_ratios_muon, \
                            flavor_ratios_neutron]
    lst_color = ['salmon', 'orange', 'limegreen']
    lst_markers = ['o', 's', '^']

    for i in range(len(lst_flavor_ratios_bf)):
        # flavor_ratios = [lst_flavor_ratios_bf[i][0], lst_flavor_ratios_bf[i][2],
            # lst_flavor_ratios_bf[i][1]] 
        tax.scatter([lst_flavor_ratios_bf[i]], \
        # tax.scatter([flavor_ratios], \
            marker=lst_markers[i], color=lst_color[i], \
            edgecolor=lst_color[i], s=40, linewidths=0.2, zorder=4, \
            alpha=1.)

    # Annotate
    tax.scatter([[0.20, 0.98, 0.]], \
        marker='o', color='salmon', edgecolor='salmon', \
        s=40, linewidths=0.1, zorder=2, alpha=1.)
    tax.ax.annotate( r'$\pi$ decay: $\left( 1:2:0 \right)_{\rm S}$', \
        xy = (0.72,0.83), xycoords='data', color='k', \
        fontsize=10, horizontalalignment='left' )
    tax.scatter([[0.235, 0.91, 0.]], \
        marker='s', color='orange', edgecolor='orange', \
        s=40, linewidths=0.1, zorder=2, alpha=1.)
    tax.ax.annotate( r'$\mu$-damped: $\left( 0:1:0 \right)_{\rm S}$', \
        xy = (0.72,0.77), xycoords='data', color='k', \
        fontsize=10, horizontalalignment='left' )
    tax.scatter([[0.27, 0.843, 0.]], \
        marker='^', color='limegreen', edgecolor='limegreen', \
        s=40, linewidths=0.1, zorder=2, alpha=1.)
    tax.ax.annotate( r'$n$ decay: $\left( 1:0:0 \right)_{\rm S}$', \
        xy = (0.72,0.71), xycoords='data', color='k', \
        fontsize=10, horizontalalignment='left' )


    tax.clear_matplotlib_ticks()
    tax._redraw_labels()
    ternary.plt.tight_layout()

    ternary.plt.savefig(fname_out, dpi=200)
    print("Saved "+fname)

    plt.close()

    return


def plot_ternary_flavor_ratios_fixed_fS_combined_v1(fname_in, fname_out, 
    plot_label, order='no', theta23='upper', year='2020', fixed_deltacp=None,
    nufit_version=None, non_unitary_bf_fS_year=None, show_benchmarks=True, 
    show_benchmark_astro_labels=True, show_benchmark_decay_labels=False, 
    output_format='png'):

    # Open the plot and format it
    mpl.rcParams['xtick.labelsize']=26
    mpl.rcParams['ytick.labelsize']=26
    mpl.rcParams['legend.fontsize']=23
    mpl.rcParams['legend.borderpad']=0.4
    mpl.rcParams['axes.labelpad']=10
    mpl.rcParams['ps.fonttype']=42
    mpl.rcParams['pdf.fonttype']=42

    colors_gen = ['#2a4d69', '#4b86b4', '#adcbe3', '#D3E4F4', '#63ace5']
    colors_010 = ['#964000', '#ff7417', '#ffbf00']
    colors_100 = ['#043927', '#00a572', '#d0f0c0']
    colors_120 = ['#420d09', '#b80f0a', '#fa8072']
    colors = [colors_010, colors_100, colors_120]

    # Boundary and gridlines
    fig, tax = ternary.figure(scale=1.0)
    tax.ax.axis("off")
    fig.set_facecolor('w')

    # Draw boundary and gridlines
    tax.boundary(linewidth=1.0)
    tax.gridlines(color="gray", multiple=0.1, linewidth=0.5, ls='-', alpha=0.5)

    fontsize = 15
    label_left = r'Fraction of $\nu_\tau$, $f_{\tau, \oplus}$'
    label_right = r'Fraction of $\nu_\mu$, $f_{\mu, \oplus}$'
    label_bottom = r'Fraction of $\nu_e$, $f_{e, \oplus}$'
    tax.left_axis_label(label_left, fontsize=fontsize, offset=0.16)
    tax.right_axis_label(label_right, fontsize=fontsize, offset=0.16)
    tax.bottom_axis_label(label_bottom, fontsize=fontsize, offset=0.08)

    # Set ticks blr
    tax.ticks(axis='blr', linewidth=0.5, multiple=0.1, offset=0.022, \
        clockwise=False, tick_formats="%.1f")

    ###########################################################################
    # Standard Model regions
    ###########################################################################

    # >>>>>
    # Plot the region for all-fS in the background
    with open(fname_in[0]+'.json') as file:
        data = json.load(file)
    fmu = data['0.997']['f_e_Earth']
    fe = data['0.997']['f_mu_Earth']
    ftau = [1.0-fe[i]-fmu[i] for i in range(len(fe))]
    fe_rot = [ 0.5 * (2.*fe[i]+fmu[i]) / (fe[i]+fmu[i]+ftau[i]) 
        for i in range(len(fe)) ]
    fmu_rot = [ sqrt(3.)/2 * (fmu[i]) / (fe[i]+fmu[i]+ftau[i]) 
        for i in range(len(fe)) ]
    ftau_rot = [1.0-fe_rot[i]-fmu_rot[i] for i in range(len(fe))]
    # tax.ax.fill(fe_rot, fmu_rot, color=colors[j])
    tax.ax.fill(fe_rot, fmu_rot, color='0.8')
    # tax.plot([[fe[i], fmu[i]] for i in range(len(fe))], 
    #     color='k', ls='-', linewidth=1.)

    # >>>>>
    if show_benchmarks:
        # Overlay the regions for the benchmarks of fS
        for k, fname in enumerate(fname_in[1:]):
            with open(fname+'.json') as file:
                data = json.load(file)
            for j, cl in enumerate(['0.997', '0.95', '0.68']): #enumerate(list(data.keys())[::-1]):
                # if not ((j == 1) or (j == 2)): continue # >>>>>
                # if not ((j == 2)): continue # >>>>>
                fmu = data[cl]['f_e_Earth']
                fe = data[cl]['f_mu_Earth']
                ftau = [1.0-fe[i]-fmu[i] for i in range(len(fe))]
                fe_rot = [ 0.5 * (2.*fe[i]+fmu[i]) / (fe[i]+fmu[i]+ftau[i]) 
                    for i in range(len(fe)) ]
                fmu_rot = [ sqrt(3.)/2 * (fmu[i]) / (fe[i]+fmu[i]+ftau[i]) 
                    for i in range(len(fe)) ]
                ftau_rot = [1.0-fe_rot[i]-fmu_rot[i] for i in range(len(fe))]
                # tax.ax.fill(fe_rot, fmu_rot, color=colors[j])
                tax.ax.fill(fe_rot, fmu_rot, color=colors[k][j])


    ###########################################################################
    # Draw IceCube 2015 contours
    ###########################################################################

    """
    # lst_filenames_ic =  [ \
    #     os.getcwd()+"/flavor_ratios_earth_ic_2015_bf.dat", \
    #     os.getcwd()+"/flavor_ratios_earth_ic_2015_1s.dat", \
    #     # os.getcwd()+"/flavor_ratios_earth_ic_2015_2s_a.dat", \
    #     # os.getcwd()+"/flavor_ratios_earth_ic_2015_2s_b.dat", \
    #     os.getcwd()+"/flavor_ratios_earth_ic_2015_3s_a.dat", \
    #     os.getcwd()+"/flavor_ratios_earth_ic_2015_3s_b.dat"]

    lst_filenames_ic =  [ \
        os.getcwd()+"/icgen2_icrc2017_pion_68cl.dat",
        os.getcwd()+"/icgen2_icrc2017_pion_90cl.dat",
        os.getcwd()+"/icgen2_icrc2017_muon_68cl.dat",
        os.getcwd()+"/icgen2_icrc2017_muon_90cl.dat"]

    lst_color = ['gold', 'gold', 'w', 'w'] #'y'
    lst_ls = ['--', '-', '--', '-']

    # Each row: f_eE, f_mE
    for file_index in range(len(lst_filenames_ic)):

        # Read data file with IceCube flavor ratios
        lst_ratios_earth = Read_Data_File(lst_filenames_ic[file_index])
        # Recast as [f_eE, f_mE, f_tE]
        lst_ratios_earth = \
            [[lst_ratios_earth[1][i], lst_ratios_earth[2][i],\
            lst_ratios_earth[0][i]] \
            for i in range(len(lst_ratios_earth[0]))]

        tax.plot(lst_ratios_earth, color=lst_color[file_index],
                ls=lst_ls[file_index], linewidth=2.5, zorder=1)

    """

    if (year == '2020' or year is None):
        # IceCube, 8 yr
        # file_ic_contour = 'icecube_8yr.json'
        file_ic_contour = 'icecube_8yr_recenter.json'
        ic_contour_label = r'IceCube 8 yr (68\%, 95\%, 99.7\% C.R.)'
    elif (year == '2030'):
        # IceCube, 15 yr
        # file_ic_contour = 'icecube_15yr.json'
        file_ic_contour = 'icecube_15yr_recenter.json'
        ic_contour_label = r'IceCube 15 yr (68\%, 95\%, 99.7\% C.R.)'
    elif (year == '2040'):
        # IceCube-Gen2 10-yr + IceCube 15-yr
        # file_ic_contour = 'icecube+gen2-inice.json'
        file_ic_contour = 'icecube+gen2-inice_recenter.json'
        ic_contour_label = r'IceCube 15 yr + Gen2 10 yr (68\%, 95\%, 99.7\% C.R.)'
    elif (year == '2040_comb'):
        # All neutrino telescopes combined
        file_ic_contour = 'combineexp_recenter.json'
        ic_contour_label = r'Combined $\nu$ telescopes (68\%, 95\%, 99.7\% C.R.)'

    # >>>>>
    tax.plot([[0.073, 0.05], [0.113, 0.05]], color='C0',
                ls='-', linewidth=1.5, zorder=1)
    tax.ax.annotate( ic_contour_label, 
        xy = (0.15,0.03), xycoords='data', color='k', 
        fontsize=10, horizontalalignment='left' )

    # >>>>>
    with open(file_ic_contour) as json_file:
        data = json.load(json_file)
        fraction_nue = data['data']['nue_fraction']
        fraction_numu = data['data']['numu_fraction']
        fraction_nutau = \
            [1.0-fraction_nue[i]-fraction_numu[i] \
            for i in range(len(fraction_nue))]
        confidence_level = data['data']['confidence_level']
        cs = tax.ax.contour(fraction_nue, fraction_numu, confidence_level,
                levels=[68., 95., 99.7], colors=None, linestyles=None, alpha=0.0)
        for j in range(len(cs.collections)):
            vertices = cs.collections[j].get_paths()[0].vertices
            x = vertices[:,0]
            y = vertices[:,1]
            z = [1.0-x[i]-y[i] for i in range(len(vertices))]
            ratios_earth = [[x[i], y[i], z[i]] for i in range(len(x))]
            tax.plot(ratios_earth, color='C0',
                    ls='-', linewidth=1., zorder=5)

    ###########################################################################
    # Annotations
    ###########################################################################

    tax.ax.annotate( plot_label, \
        xy = (-0.03, 0.845), xycoords='data', color='k', \
        fontsize=10, horizontalalignment='left' )

    # >>>>>
    tax.ax.scatter([-0.01], [0.80], \
        marker='s', color=colors_120[2], edgecolor='k', \
        s=60, linewidths=0.1, zorder=2, alpha=1.)
    tax.ax.scatter([0.03], [0.80], \
        marker='s', color=colors_010[2], edgecolor='k', \
        s=60, linewidths=0.1, zorder=2, alpha=1.)
    tax.ax.scatter([0.07], [0.80], \
        marker='s', color=colors_100[2], edgecolor='k', \
        s=60, linewidths=0.1, zorder=2, alpha=1.)
    tax.ax.annotate( r'68\%~C.R.', \
        xy = (0.10,0.787), xycoords='data', color='k', \
        fontsize=10, horizontalalignment='left' )

    # >>>>>
    tax.ax.scatter([-0.01], [0.75], \
        marker='s', color=colors_120[1], edgecolor=colors_gen[2], \
        s=60, linewidths=0.1, zorder=2, alpha=1.)
    tax.ax.scatter([0.03], [0.75], \
        marker='s', color=colors_010[1], edgecolor='k', \
        s=60, linewidths=0.1, zorder=2, alpha=1.)
    tax.ax.scatter([0.07], [0.75], \
        marker='s', color=colors_100[1], edgecolor='k', \
        s=60, linewidths=0.1, zorder=2, alpha=1.)
    tax.ax.annotate( r'95\%~C.R.', \
        xy = (0.10,0.737), xycoords='data', color='k', \
        fontsize=10, horizontalalignment='left' )

    # >>>>>
    tax.ax.scatter([-0.01], [0.70], \
        marker='s', color=colors_120[0], edgecolor=colors_gen[1], \
        s=60, linewidths=0.1, zorder=2, alpha=1.)
    tax.ax.scatter([0.03], [0.70], \
        marker='s', color=colors_010[0], edgecolor='k', \
        s=60, linewidths=0.1, zorder=2, alpha=1.)
    tax.ax.scatter([0.07], [0.70], \
        marker='s', color=colors_100[0], edgecolor='k', \
        s=60, linewidths=0.1, zorder=2, alpha=1.)
    tax.ax.annotate( r'99.7\%~C.R.', \
        xy = (0.10,0.687), xycoords='data', color='k', \
        fontsize=10, horizontalalignment='left' )

    ####################################################################
    # Add regions of selected compositions
    ####################################################################

    if (nufit_version is None):

        if (order.lower() == 'no'):
            if (theta23 == 'upper'):
                if (fixed_deltacp is None):
                    # Best-fit flavor ratios for benchmark points
                    flavor_ratios_pion \
                        = [0.2982716010150718, 0.359320837204981, 0.34240756177994724]
                    flavor_ratios_muon \
                        = [0.17128846521772986, 0.4533370231986065, 0.37537451158366353]
                    flavor_ratios_neutron \
                        = [0.5522378726097557, 0.17128846521772986, 0.2764736621725147]
                    # Best-fit flavor content of nu_1, nu_2, nu_3
                    flavor_cont_nu1 \
                        = [0.6808662622895465, 0.07368823020070894, 0.24544550750974464]
                    flavor_cont_nu2 \
                        = [0.29692740578832083, 0.3659953978327661, 0.33707719637891315]
                    flavor_cont_nu3 \
                        = [0.022206331922132744, 0.560316371966525, 0.41747729611134216]
                else:
                    if (fixed_deltacp == 0.0):
                        flavor_ratios_pion \
                            = [0.3322161727912954, 0.3480562782080484, 0.319727549000656]
                        flavor_ratios_muon \
                            = [0.22220532288206532, 0.41098175587104, 0.3668129212468944]
                        flavor_ratios_neutron \
                            = [0.5522378726097557, 0.22220532288206532, 0.22555680450817917]
                    elif (fixed_deltacp == 0.5*np.pi):
                        flavor_ratios_pion \
                            = [0.3148647971219528, 0.34795485221718947, 0.33718035066085766]
                        flavor_ratios_muon \
                            = [0.1961782593780514, 0.4238431486367586, 0.3799785919851899]
                        flavor_ratios_neutron \
                            = [0.5522378726097557, 0.1961782593780514, 0.2515838680121932]
                    elif (fixed_deltacp == np.pi):
                        flavor_ratios_pion \
                            = [0.3148647971219528, 0.34795485221718947, 0.33718035066085766]
                        flavor_ratios_muon \
                            = [0.1961782593780514, 0.4238431486367586, 0.3799785919851899]
                        flavor_ratios_neutron \
                            = [0.5522378726097557, 0.1961782593780514, 0.2515838680121932]
                    else:
                        print("Invalid value of fixed_deltacp")
                        quit()
            elif (theta23 == 'lower'):
                if (fixed_deltacp is None):
                    flavor_ratios_pion \
                        = [0.32800656074697265, 0.34529257725052065, 0.3267008620025065]
                    flavor_ratios_muon \
                        = [0.2158909048155812, 0.40999341346799034, 0.3741156817164281]
                    flavor_ratios_neutron \
                        = [0.5522378726097557, 0.2158909048155812, 0.23187122257466325]
                else:
                    print("Cannot choose fixed_deltacp != None for theta = 'lower'")
                    quit()
            elif (theta23 == 'max'):
                if (fixed_deltacp is None):
                    flavor_ratios_pion \
                        = [0.31656024191908616, 0.3474455438823001, 0.3359942141986137]
                    flavor_ratios_muon \
                        = [0.19872142657375141, 0.4218076025365745, 0.37947097088967396]
                    flavor_ratios_neutron \
                        = [0.5522378726097557, 0.19872142657375141, 0.24904070081649313]
                else:
                    print("Cannot choose fixed_deltacp != None for theta = 'max'")
                    quit()
        elif (order.lower() == 'io'):
            if (theta23 == 'upper'):
                if (fixed_deltacp is None):
                    # Best-fit flavor ratios for benchmark points
                    flavor_ratios_pion \
                        = [0.3180471687709684, 0.3474505778780882, 0.3345022533509433]
                    flavor_ratios_muon \
                        = [0.2010801776404412, 0.4206357779969117, 0.378284044362647]
                    flavor_ratios_neutron \
                        = [0.5519811510320229, 0.2010801776404412, 0.24693867132753594]
                    # Best-fit flavor content of nu_1, nu_2, nu_3
                    flavor_cont_nu1 \
                        = [0.6806215638904127, 0.15226374101164764, 0.16711469509793975]
                    flavor_cont_nu2 \
                        = [0.2970466806723904, 0.2858044688729501, 0.41714885045465944]
                    flavor_cont_nu3 \
                        = [0.02233175543719699, 0.5619317901154023, 0.41573645444740065]
                else:
                    if (fixed_deltacp == 0.0):
                        flavor_ratios_pion \
                            = [0.33181153325437795, 0.34861177147789035, 0.3195766952677316]
                        flavor_ratios_muon \
                            = [0.22172672436555554, 0.41205429503405777, 0.3662189806003866]
                        flavor_ratios_neutron \
                            = [0.5519811510320229, 0.22172672436555554, 0.22629212460242162]
                    elif (fixed_deltacp == 1.0*np.pi):
                        flavor_ratios_pion \
                            = [0.29705684749618344, 0.3605546628791725, 0.34238848962464397]
                        flavor_ratios_muon \
                            = [0.16959469572826372, 0.4560346464546269, 0.3743706578171092]
                        flavor_ratios_neutron \
                            = [0.5519811510320229, 0.16959469572826372, 0.2784241532397134]
                    elif (fixed_deltacp == 1.5*np.pi):
                        flavor_ratios_pion \
                            = [0.3144341903752807, 0.34842595462849185, 0.3371398549962274]
                        flavor_ratios_muon \
                            = [0.19566071004690963, 0.424808576919283, 0.37953071303380737]
                        flavor_ratios_neutron \
                            = [0.5519811510320229, 0.19566071004690963, 0.25235813892106745]
                    else:
                        print("Invalid value of fixed_deltacp")
                        quit()
            elif (theta23 == 'lower'):
                if (fixed_deltacp is None):
                    flavor_ratios_pion \
                        = [0.34834738484344396, 0.3277980374633891, 0.32385457769316683]
                    flavor_ratios_muon \
                        = [0.2465305017491545, 0.36843180532050646, 0.38503769293033896]
                    flavor_ratios_neutron \
                        = [0.5519811510320229, 0.2465305017491545, 0.20148834721882264]
                else:
                    print("Cannot choose fixed_deltacp != None for theta = 'lower'")
                    quit()
            elif (theta23 == 'max'):
                if (fixed_deltacp is None):
                    flavor_ratios_pion \
                        = [0.33698743757165983, 0.3317785331561326, 0.33123402927220735]
                    flavor_ratios_muon \
                        = [0.22949058084147833, 0.3829225093134597, 0.3875869098450617]
                    flavor_ratios_neutron \
                        = [0.5519811510320229, 0.22949058084147833, 0.21852826812649873]
                else:
                    print("Cannot choose fixed_deltacp != None for theta = 'max'")
                    quit()

    else: # nufit_version is not None

        # if (nufit_version == 'nufit_v1.0'):
        #     if (order.lower() == 'no'):
        #         flavor_ratios_pion \
        #             = xxx
        #         flavor_ratios_muon \
        #             = xxx
        #         flavor_ratios_neutron \
        #             = xxx
        #     elif (order.lower() == 'io'):
        #         flavor_ratios_pion \
        #             = xxx
        #         flavor_ratios_muon \
        #             = xxx
        #         flavor_ratios_neutron \
        #             = xxx
        # if (nufit_version == 'nufit_v1.1'):
        #     if (order.lower() == 'no'):
        #         flavor_ratios_pion \
        #             = xxx
        #         flavor_ratios_muon \
        #             = xxx
        #         flavor_ratios_neutron \
        #             = xxx
        #     elif (order.lower() == 'io'):
        #         flavor_ratios_pion \
        #             = xxx
        #         flavor_ratios_muon \
        #             = xxx
        #         flavor_ratios_neutron \
        #             = xxx
        # if (nufit_version == 'nufit_v1.2'):
        #     if (order.lower() == 'no'):
        #         flavor_ratios_pion \
        #             = xxx
        #         flavor_ratios_muon \
        #             = xxx
        #         flavor_ratios_neutron \
        #             = xxx
        #     elif (order.lower() == 'io'):
        #         flavor_ratios_pion \
        #             = xxx
        #         flavor_ratios_muon \
        #             = xxx
        #         flavor_ratios_neutron \
        #             = xxx
        # if (nufit_version == 'nufit_v1.3'):
        #     if (order.lower() == 'no'):
        #         flavor_ratios_pion \
        #             = xxx
        #         flavor_ratios_muon \
        #             = xxx
        #         flavor_ratios_neutron \
        #             = xxx
        #     elif (order.lower() == 'io'):
        #         flavor_ratios_pion \
        #             = xxx
        #         flavor_ratios_muon \
        #             = xxx
        #         flavor_ratios_neutron \
        #             = xxx
        if (nufit_version == 'nufit_v2.0'):
            if (order.lower() == 'no'):
                flavor_ratios_pion \
                    = [0.3556613043722851, 0.3255837945733705, 0.318754901054344]
                flavor_ratios_muon \
                    = [0.25727620733858775, 0.3597375881907619, 0.38298620447064996]
                flavor_ratios_neutron \
                    = [0.5524314984396799, 0.25727620733858775, 0.19029229422173227]
            elif (order.lower() == 'io'):
                flavor_ratios_pion \
                    = [0.3085776926666689, 0.3516915412977889, 0.3397307660355422]
                flavor_ratios_muon \
                    = [0.18670502760224328, 0.43418479814556177, 0.37911017425219495]
                flavor_ratios_neutron \
                    = [0.5523230227955201, 0.18670502760224328, 0.2609719496022368]
        if (nufit_version == 'nufit_v2.1'):
            if (order.lower() == 'no'):
                flavor_ratios_pion \
                    = [0.31499629123871653, 0.3481219497568094, 0.3368817590044737]
                flavor_ratios_muon \
                    = [0.19770980806911492, 0.42332802060065666, 0.37896217133022814]
                flavor_ratios_neutron \
                    = [0.5495692575779199, 0.19770980806911492, 0.252720934352965]
            elif (order.lower() == 'io'):
                flavor_ratios_pion \
                    = [0.3090660431370772, 0.35143623465391693, 0.33949772220900554]
                flavor_ratios_muon \
                    = [0.1890301412593759, 0.4326392813511875, 0.3783305773894363]
                flavor_ratios_neutron \
                    = [0.5491378468924799, 0.1890301412593759, 0.26183201184814403]
        if (nufit_version == 'nufit_v2.2'):
            if (order.lower() == 'no'):
                flavor_ratios_pion \
                    = [0.3541798378010116, 0.3264941108322114, 0.31932605136677694]
                flavor_ratios_muon \
                    = [0.2564473535897559, 0.3615174894534392, 0.382035156956805]
                flavor_ratios_neutron \
                    = [0.5496448062235231, 0.2564473535897559, 0.19390784018672091]
            elif (order.lower() == 'io'):
                flavor_ratios_pion \
                    = [0.31156306080426477, 0.3512832938335091, 0.3371536453622259]
                flavor_ratios_muon \
                    = [0.1925869393553972, 0.4306314710725651, 0.37678158957203756]
                flavor_ratios_neutron \
                    = [0.5495153037019999, 0.1925869393553972, 0.25789775694260275]
        if (nufit_version == 'nufit_v3.0'):
            if (order.lower() == 'no'):
                flavor_ratios_pion \
                    = [0.34569828207852565, 0.33114156111530313, 0.32316015680617094]
                flavor_ratios_muon \
                    = [0.24300229079762686, 0.37521119627414135, 0.3817865129282315]
                flavor_ratios_neutron \
                    = [0.5510902646403233, 0.24300229079762686, 0.20590744456204982]
            elif (order.lower() == 'io'):
                flavor_ratios_pion \
                    = [0.31324903482344135, 0.35147907508403603, 0.3352718900925226]
                flavor_ratios_muon \
                    = [0.19439875631305437, 0.4300192344695269, 0.3755820092174188]
                flavor_ratios_neutron \
                    = [0.5509495918442153, 0.19439875631305437, 0.2546516518427304]
        if (nufit_version == 'nufit_v3.1'):
            if (order.lower() == 'no'):
                flavor_ratios_pion \
                    = [0.3053534737546712, 0.3527746876435402, 0.34187183860178844]
                flavor_ratios_muon \
                    = [0.18301214319568435, 0.43765595986746814, 0.37933189693684743]
                flavor_ratios_neutron \
                    = [0.5500361348726449, 0.18301214319568435, 0.2669517219316705]
            elif (order.lower() == 'io'):
                flavor_ratios_pion \
                    = [0.3182553364828055, 0.34681283459069917, 0.3349318289264954]
                flavor_ratios_muon \
                    = [0.20245670395710264, 0.4189908999074975, 0.37855239613540004]
                flavor_ratios_neutron \
                    = [0.5498526015342113, 0.20245670395710264, 0.24769069450868617]
        if (nufit_version == 'nufit_v3.2'):
            if (order.lower() == 'no'):
                flavor_ratios_pion \
                    = [0.3135479127848382, 0.346028878345814, 0.34042320886934785]
                flavor_ratios_muon \
                    = [0.1953631853698009, 0.4213617248338205, 0.38327508979637864]
                flavor_ratios_neutron \
                    = [0.5499173676149128, 0.1953631853698009, 0.2547194470152863]
            elif (order.lower() == 'io'):
                flavor_ratios_pion \
                    = [0.32198676464987386, 0.3423370406991942, 0.3356761946509318]
                flavor_ratios_muon \
                    = [0.2081347790102386, 0.409438171543672, 0.38242704944608913]
                flavor_ratios_neutron \
                    = [0.5496907359291444, 0.2081347790102386, 0.24217448506061715]
        if (nufit_version == 'nufit_v4.0'):
            if (order.lower() == 'no'):
                flavor_ratios_pion \
                    = [0.2988473594059442, 0.3589894623100733, 0.3421631782839824]
                flavor_ratios_muon \
                    = [0.17459388557291625, 0.45118725067865184, 0.3742188637484317]
                flavor_ratios_neutron \
                    = [0.5473543070720001, 0.17459388557291625, 0.2780518073550838]
            elif (order.lower() == 'io'):
                flavor_ratios_pion \
                    = [0.3153279935881048, 0.34979331257813706, 0.334878693833758]
                flavor_ratios_muon \
                    = [0.19943830128706735, 0.424970818223672, 0.3755908804892607]
                flavor_ratios_neutron \
                    = [0.5471073781901797, 0.19943830128706735, 0.2534543205227526]
        if (nufit_version == 'nufit_v4.1'):
            if (order.lower() == 'no'):
                flavor_ratios_pion \
                    = [0.30438521493764115, 0.3536866096566163, 0.3419281754057424]
                flavor_ratios_muon \
                    = [0.18288455868137182, 0.4390876351442386, 0.3780278061743895]
                flavor_ratios_neutron \
                    = [0.5473865274501799, 0.18288455868137182, 0.2697289138684482]
            elif (order.lower() == 'io'):
                flavor_ratios_pion \
                    = [0.3202665854332483, 0.3448720154249937, 0.33486139914175794]
                flavor_ratios_muon \
                    = [0.20682472295246254, 0.41389566166125935, 0.3792796153862782]
                flavor_ratios_neutron \
                    = [0.5471503103948199, 0.20682472295246254, 0.24602496665271745]
        if (nufit_version == 'nufit_v5.0'):
            if (order.lower() == 'no'):
                flavor_ratios_pion \
                    = [0.29829637092300665, 0.3593068266653908, 0.34239680241160236]
                flavor_ratios_muon \
                    = [0.17144024550093245, 0.45324011724762003, 0.37531963725144746]
                flavor_ratios_neutron \
                    = [0.5520086217671551, 0.17144024550093245, 0.2765511327319123]
            elif (order.lower() == 'io'):
                flavor_ratios_pion \
                    = [0.3179873235848132, 0.34751722386393324, 0.33449545255125346]
                flavor_ratios_muon \
                    = [0.20107959603042933, 0.42073603778068525, 0.37818436618888535]
                flavor_ratios_neutron \
                    = [0.5518027786935809, 0.20107959603042933, 0.24711762527598982]

    if (non_unitary_bf_fS_year == '2020'):
        flavor_ratios_pion = [0.321788, 0.339936, 0.338276]
        flavor_ratios_muon = [0.207372, 0.403323, 0.389305]
        flavor_ratios_neutron = [0.554584, 0.210965, 0.234451]
    elif (non_unitary_bf_fS_year == '2040'):
        flavor_ratios_pion = [0.304037, 0.353711, 0.342252]
        flavor_ratios_muon = [0.180976, 0.437637, 0.381387]
        flavor_ratios_neutron = [0.553604, 0.183509, 0.262886]

    if show_benchmark_astro_labels:
    
        lst_flavor_ratios_bf = [flavor_ratios_pion, flavor_ratios_muon, \
            flavor_ratios_neutron]
        lst_color = ['salmon', 'orange', 'limegreen']
        lst_markers = ['o', 's', '^']

        for i in range(len(lst_flavor_ratios_bf)):
            tax.scatter([lst_flavor_ratios_bf[i]], \
                marker=lst_markers[i], 
                color=lst_color[i], \
                edgecolor='k', #lst_color[i], 
                s=20, #40 
                linewidths=0.5, 
                zorder=4, #4
                alpha=1.)

        # Annotate
        tax.scatter([[0.20, 0.98, 0.]], \
            marker='o', color='salmon', edgecolor='k', \
            s=40, linewidths=0.5, zorder=2, alpha=1.)
        tax.ax.annotate( r'$\pi$ decay: $\left( 1:2:0 \right)_{\rm S}$', \
            xy = (0.72,0.83), xycoords='data', color='k', \
            fontsize=10, horizontalalignment='left' )
        tax.scatter([[0.235, 0.91, 0.]], \
            marker='s', color='orange', edgecolor='k', \
            s=40, linewidths=0.5, zorder=2, alpha=1.)
        tax.ax.annotate( r'$\mu$-damped: $\left( 0:1:0 \right)_{\rm S}$', \
            xy = (0.72,0.77), xycoords='data', color='k', \
            fontsize=10, horizontalalignment='left' )
        tax.scatter([[0.267, 0.843, 0.]], \
            marker='^', color='limegreen', edgecolor='k', \
            s=40, linewidths=0.5, zorder=2, alpha=1.)
        tax.ax.annotate( r'$n$ decay: $\left( 1:0:0 \right)_{\rm S}$', \
            xy = (0.715,0.71), xycoords='data', color='k', \
            fontsize=10, horizontalalignment='left' )

    elif show_benchmark_decay_labels:
    
        lst_flavor_content_bf = [flavor_cont_nu1, flavor_cont_nu2, 
            flavor_cont_nu3]

        lst_color = ['salmon', 'limegreen', 'orange']
        lst_markers = ['o', '^', 's']

        for i in range(len(lst_flavor_content_bf)):
            tax.scatter([lst_flavor_content_bf[i]], \
                marker=lst_markers[i], 
                color=lst_color[i], \
                edgecolor='k', #lst_color[i], 
                s=20, #40 
                linewidths=0.5, 
                zorder=4, #4
                alpha=1.)

        # Annotate
        tax.scatter([[0.20, 0.98, 0.]], \
            marker='o', color='salmon', edgecolor='k', \
            s=40, linewidths=0.5, zorder=2, alpha=1.)
        tax.ax.annotate( r'$\nu_1$', \
            xy = (0.72,0.83), xycoords='data', color='k', \
            fontsize=10, horizontalalignment='left' )
        tax.scatter([[0.235, 0.91, 0.]], \
            marker='s', color='orange', edgecolor='k', \
            s=40, linewidths=0.5, zorder=2, alpha=1.)
        tax.ax.annotate( r'$\nu_2$', \
            xy = (0.72,0.77), xycoords='data', color='k', \
            fontsize=10, horizontalalignment='left' )
        tax.scatter([[0.267, 0.843, 0.]], \
            marker='^', color='limegreen', edgecolor='k', \
            s=40, linewidths=0.5, zorder=2, alpha=1.)
        tax.ax.annotate( r'$\nu_3$', \
            xy = (0.715,0.71), xycoords='data', color='k', \
            fontsize=10, horizontalalignment='left' )


    tax.clear_matplotlib_ticks()
    tax._redraw_labels()
    ternary.plt.tight_layout()

    ternary.plt.savefig(fname_out+'.'+output_format, dpi=300)
    print("Saved "+fname_out+'.'+output_format)

    plt.close()

    return


def plot_ternary_flavor_ratios_fixed_fS_combined_2020_vs_2040(fname_in, 
    fname_out, plot_label, order='no', theta23='upper', year='2020', 
    fixed_deltacp=None, nufit_version=None, non_unitary_bf_fS_year=None,
    show_benchmarks=True, show_benchmark_astro_labels=True, 
    show_benchmark_decay_labels=False, show_full_exp_labels=True, 
    output_format='png'):

    # Open the plot and format it
    mpl.rcParams['xtick.labelsize']=26
    mpl.rcParams['ytick.labelsize']=26
    mpl.rcParams['legend.fontsize']=23
    mpl.rcParams['legend.borderpad']=0.4
    mpl.rcParams['axes.labelpad']=10
    mpl.rcParams['ps.fonttype']=42
    mpl.rcParams['pdf.fonttype']=42

    if show_benchmark_astro_labels:
        colors_gen = ['#2a4d69', '#4b86b4', '#adcbe3', '#D3E4F4', '#63ace5']
        colors_010 = ['#964000', '#ff7417', '#ffbf00']
        colors_100 = ['#043927', '#00a572', '#d0f0c0']
        colors_120 = ['#420d09', '#b80f0a', '#fa8072']
        colors = [colors_010, colors_100, colors_120]
    elif show_benchmark_decay_labels:
        colors_gen = ['#2a4d69', '#4b86b4', '#adcbe3', '#D3E4F4', '#63ace5']
        colors_010 = [None, '#7E1F86', '#B9A0CF'] # nu_1
        # colors_100 = [None, '#8C6D69', '#BFACAA'] # nu_2
        # colors_100 = [None, '#DAB558', '#F4E9CD'] # nu_2
        # colors_100 = [None, '#40858C', '#7EBDC3'] # nu_2
        # colors_100 = [None, '#188B7B', '#1EA896'] # nu_2
        colors_100 = [None, '#F5D07A', '#FCF2D9'] # nu_2
        colors_120 = [None, '#05668D', '#91C4F2'] # nu_3
        colors = [colors_010, colors_100, colors_120]

    # Boundary and gridlines
    fig, tax = ternary.figure(scale=1.0)
    tax.ax.axis("off")
    fig.set_facecolor('w')

    # Draw boundary and gridlines
    tax.boundary(linewidth=1.0)
    tax.gridlines(color="gray", multiple=0.1, linewidth=0.5, ls='-', alpha=0.5)

    fontsize = 15
    label_left = r'Fraction of $\nu_\tau$, $f_{\tau, \oplus}$'
    label_right = r'Fraction of $\nu_\mu$, $f_{\mu, \oplus}$'
    label_bottom = r'Fraction of $\nu_e$, $f_{e, \oplus}$'
    tax.left_axis_label(label_left, fontsize=fontsize, offset=0.16)
    tax.right_axis_label(label_right, fontsize=fontsize, offset=0.16)
    tax.bottom_axis_label(label_bottom, fontsize=fontsize, offset=0.08)

    # Set ticks blr
    tax.ticks(axis='blr', linewidth=0.5, multiple=0.1, offset=0.022, \
        clockwise=False, tick_formats="%.1f")

    ###########################################################################
    # Standard Model regions
    ###########################################################################

    for file_index in range(len(fname_in)):

        # Plot the region for all-fS in the background
        with open(fname_in[file_index][0]+'.json') as file:
            data = json.load(file)
        fmu = data['0.997']['f_e_Earth']
        fe = data['0.997']['f_mu_Earth']
        ftau = [1.0-fe[i]-fmu[i] for i in range(len(fe))]
        fe_rot = [ 0.5 * (2.*fe[i]+fmu[i]) / (fe[i]+fmu[i]+ftau[i]) 
            for i in range(len(fe)) ]
        fmu_rot = [ sqrt(3.)/2 * (fmu[i]) / (fe[i]+fmu[i]+ftau[i]) 
            for i in range(len(fe)) ]
        ftau_rot = [1.0-fe_rot[i]-fmu_rot[i] for i in range(len(fe))]
        if file_index == 0:
            tax.ax.fill(fe_rot, fmu_rot, color='0.6')
        elif file_index == 1:
            tax.ax.fill(fe_rot, fmu_rot, color='0.8')


    for file_index in range(len(fname_in)):

        # Overlay the regions for the benchmarks of fS
        for k, fname in enumerate(fname_in[file_index][1:]):
            # if not ((k == 2) or (k == 0)): continue
            with open(fname+'.json') as file:
                data = json.load(file)
            for j, cl in enumerate(['0.997']): #enumerate(list(data.keys())[::-1]):
                fmu = data[cl]['f_e_Earth']
                fe = data[cl]['f_mu_Earth']
                ftau = [1.0-fe[i]-fmu[i] for i in range(len(fe))]
                fe_rot = [ 0.5 * (2.*fe[i]+fmu[i]) / (fe[i]+fmu[i]+ftau[i]) 
                    for i in range(len(fe)) ]
                fmu_rot = [ sqrt(3.)/2 * (fmu[i]) / (fe[i]+fmu[i]+ftau[i]) 
                    for i in range(len(fe)) ]
                ftau_rot = [1.0-fe_rot[i]-fmu_rot[i] for i in range(len(fe))]
                # tax.ax.fill(fe_rot, fmu_rot, color=colors[j])
                tax.ax.fill(fe_rot, fmu_rot, color=colors[k][file_index+1])


    ###########################################################################
    # Draw IceCube contours
    ###########################################################################

    # ls = ['--', '-']
    ls = [':', '--', '-']

    for year_index, year in enumerate(['2020', '2040_gen2', '2040_combined']):
        if (year == '2020' or year is None):
            # IceCube, 8 yr
            # file_ic_contour = 'icecube_8yr.json'
            file_ic_contour = 'icecube_8yr_recenter.json'
            # ic_contour_label = r'2020: IceCube 8 yr (99.7\% C.R.)'
            ic_contour_label = r'2020 (proj.): IC 8 yr (99.7\% C.R.)'
            xy = (0.15,0.11)
            line_legend = [[0.030, 0.14], [0.072, 0.14]]
            # xy = (0.15,0.11)
            # line_legend = [[0.050, 0.095], [0.09, 0.095]]
        elif (year == '2030'):
            # IceCube, 15 yr
            # file_ic_contour = 'icecube_15yr.json'
            file_ic_contour = 'icecube_15yr_recenter.json'
            ic_contour_label = r'IceCube 15 yr (99.7\% C.R.)'
        elif (year == '2040_gen2'):
            # IceCube-Gen2 10-yr + IceCube 15-yr
            # file_ic_contour = 'icecube+gen2-inice.json'
            file_ic_contour = 'icecube+gen2-inice_recenter.json'
            # ic_contour_label = r'2040: IceCube 15 yr + Gen2 10 yr (99.7\% C.R.)'
            ic_contour_label = r'2040 (proj.): IC 15 yr + Gen2 10 yr (99.7\% C.R.)'
            xy = (0.15,0.07)
            line_legend = [[0.050, 0.095], [0.091, 0.095]]
            # xy = (0.15,0.07)
            # line_legend = [[0.073, 0.095], [0.113, 0.095]]
        elif (year == '2040_combined'):
           # All neutrino telescopes combined
            file_ic_contour = 'combineexp_recenter.json'
            ic_contour_label = r'2040 (proj.): Combined $\nu$ telescopes (99.7\% C.R.)'
            xy = (0.15,0.03)
            # line_legend = [[0.050, 0.095], [0.09, 0.095]]
            line_legend = [[0.073, 0.05], [0.113, 0.05]]

        tax.plot(line_legend, color='C0',
            ls=ls[year_index], linewidth=1.5, zorder=2)
        tax.ax.annotate( ic_contour_label, 
            xy = xy, xycoords='data', color='k', 
            fontsize=10, horizontalalignment='left' )

        # 2015 contours
        color_2015 = '0.7' #adcbe3' #'0.5'
        # file_contour = os.getcwd()+'/ic_2015_excl_1s.dat'
        # file_contour = os.getcwd()+'/ic_2015_excl_2s.dat'
        file_contour = os.getcwd()+'/ic_2015_excl_3s.dat'
        fe, fmu = np.loadtxt(file_contour, unpack=True)
        ftau = [1.0-fe[i]-fmu[i] for i in range(len(fe))]
        fe_rot = [ 0.5 * (2.*fe[i]+fmu[i]) / (fe[i]+fmu[i]+ftau[i]) 
            for i in range(len(fe)) ]
        fmu_rot = [ sqrt(3.)/2 * (fmu[i]) / (fe[i]+fmu[i]+ftau[i]) 
            for i in range(len(fe)) ]
        ftau_rot = [1.0-fe_rot[i]-fmu_rot[i] for i in range(len(fe))]
        # tax.ax.fill(fe_rot, fmu_rot, color=color_2015,
        #     ls='-', linewidth=1., alpha=0.5, zorder=1)
        tax.ax.plot(fe_rot, fmu_rot, color=color_2015,
            ls='-.', #'-', 
            linewidth=1., 
            alpha=1.0, 
            zorder=1)
        fe_bf = 0.5
        fmu_bf = 0.5
        ftau_bf = 1.0-fe_bf-fmu_bf
        fe_bf_rot = 0.5 * (2.*fe_bf+fmu_bf) / (fe_bf+fmu_bf+ftau_bf) 
        fmu_bf_rot = sqrt(3.)/2 * (fmu_bf) / (fe_bf+fmu_bf+ftau_bf) 
        tax.ax.scatter([fe_bf_rot], [fmu_bf_rot], \
            marker='*', color=color_2015, edgecolor='k', \
            s=60, linewidths=0.5, zorder=6, alpha=1.)
        tax.ax.annotate( r'2015 (99.7\% C.L.)', \
            xy = (0.35,0.565), xycoords='data', color='0.6', \
            fontsize=10, horizontalalignment='left',
            rotation=22. )

        # > 2020 contours
        with open(file_ic_contour) as json_file:
            data = json.load(json_file)
            fraction_nue = data['data']['nue_fraction']
            fraction_numu = data['data']['numu_fraction']
            fraction_nutau = \
                [1.0-fraction_nue[i]-fraction_numu[i] \
                for i in range(len(fraction_nue))]
            confidence_level = data['data']['confidence_level']
            cs = tax.ax.contour(fraction_nue, fraction_numu, confidence_level,
                    levels=[99.7], colors=None, linestyles=None, alpha=0.0)
            for j in range(len(cs.collections)):
                vertices = cs.collections[j].get_paths()[0].vertices
                x = vertices[:,0]
                y = vertices[:,1]
                z = [1.0-x[i]-y[i] for i in range(len(vertices))]
                ratios_earth = [[x[i], y[i], z[i]] for i in range(len(x))]
                tax.plot(ratios_earth, color='C0',
                        ls=ls[year_index], linewidth=1., zorder=5)

    ###########################################################################
    # Annotations
    ###########################################################################

    tax.ax.annotate( plot_label, \
        xy = (-0.03, 0.845), xycoords='data', color='k', \
        fontsize=10, horizontalalignment='left' )

    if show_benchmark_astro_labels:
        color_2020_0 = colors_120[1]
        color_2020_1 = colors_010[1]
        color_2020_2 = colors_100[1]
        color_2040_0 = colors_120[2]
        color_2040_1 = colors_010[2]
        color_2040_2 = colors_100[2]

    if show_benchmark_decay_labels:
        color_2020_0 = colors_120[1]
        color_2020_1 = colors_100[1]
        color_2020_2 = colors_010[1]
        color_2040_0 = colors_120[2]
        color_2040_1 = colors_100[2]
        color_2040_2 = colors_010[2]

    tax.ax.scatter([-0.01], [0.80], \
        marker='s', color=color_2020_0, edgecolor='k', \
        s=60, linewidths=0.1, zorder=2, alpha=1.)
    tax.ax.scatter([0.03], [0.80], \
        marker='s', color=color_2020_1, edgecolor='k', \
        s=60, linewidths=0.1, zorder=2, alpha=1.)
    tax.ax.scatter([0.07], [0.80], \
        marker='s', color=color_2020_2, edgecolor='k', \
        s=60, linewidths=0.1, zorder=2, alpha=1.)

    tax.ax.scatter([-0.01], [0.75], \
        marker='s', color=color_2040_0, edgecolor='k', \
        s=60, linewidths=0.1, zorder=2, alpha=1.)
    tax.ax.scatter([0.03], [0.75], \
        marker='s', color=color_2040_1, edgecolor='k', \
        s=60, linewidths=0.1, zorder=2, alpha=1.)
    tax.ax.scatter([0.07], [0.75], \
        marker='s', color=color_2040_2, edgecolor='k', \
        s=60, linewidths=0.1, zorder=2, alpha=1.)


    if (show_full_exp_labels == True):
    
        tax.ax.annotate( r'2020: NuFit 5.0', \
            xy = (0.10,0.787), xycoords='data', color='k', \
            fontsize=10, horizontalalignment='left' )
        tax.ax.annotate( r'2040: JUNO', \
            xy = (0.10,0.737), xycoords='data', color='k', \
            fontsize=10, horizontalalignment='left' )
        tax.ax.annotate( r'+ DUNE',
            xy = (0.18,0.698), xycoords='data', color='k', \
            fontsize=10, horizontalalignment='left' )
        tax.ax.annotate( r'+ HK',
            xy = (0.18,0.66), xycoords='data', color='k', \
            fontsize=10, horizontalalignment='left' )
        # tax.ax.annotate( r'+ DUNE',
        #     xy = (0.12,0.695), xycoords='data', color='k', \
        #     fontsize=10, horizontalalignment='left' )
        # tax.ax.annotate( r'+ Hyper-K',
        #     xy = (0.12,0.66), xycoords='data', color='k', \
        #     fontsize=10, horizontalalignment='left' )

    else:

        tax.ax.annotate( r'2020', 
            xy = (0.10,0.787), xycoords='data', color='k', \
            fontsize=10, horizontalalignment='left' )
        tax.ax.annotate( r'2040',
            xy = (0.10,0.737), xycoords='data', color='k', \
            fontsize=10, horizontalalignment='left' )

    # tax.ax.scatter([-0.01], [0.70], \
    #     marker='s', color=colors_120[0], edgecolor=colors_gen[1], \
    #     s=60, linewidths=0.1, zorder=2, alpha=1.)
    # tax.ax.scatter([0.03], [0.70], \
    #     marker='s', color=colors_010[0], edgecolor='k', \
    #     s=60, linewidths=0.1, zorder=2, alpha=1.)
    # tax.ax.scatter([0.07], [0.70], \
    #     marker='s', color=colors_100[0], edgecolor='k', \
    #     s=60, linewidths=0.1, zorder=2, alpha=1.)
    # tax.ax.annotate( r'99.7\%~C.L.', \
    #     xy = (0.10,0.687), xycoords='data', color='k', \
    #     fontsize=10, horizontalalignment='left' )

    ####################################################################
    # Add regions of selected compositions
    ####################################################################


    if (nufit_version is None):

        if (order.lower() == 'no'):
            if (theta23 == 'upper'):
                if (fixed_deltacp is None):
                    # Best-fit flavor ratios for benchmark points
                    flavor_ratios_pion \
                        = [0.2982716010150718, 0.359320837204981, 0.34240756177994724]
                    flavor_ratios_muon \
                        = [0.17128846521772986, 0.4533370231986065, 0.37537451158366353]
                    flavor_ratios_neutron \
                        = [0.5522378726097557, 0.17128846521772986, 0.2764736621725147]
                    # Best-fit flavor content of nu_1, nu_2, nu_3
                    flavor_cont_nu1 \
                        = [0.6808662622895465, 0.07368823020070894, 0.24544550750974464]
                    flavor_cont_nu2 \
                        = [0.29692740578832083, 0.3659953978327661, 0.33707719637891315]
                    flavor_cont_nu3 \
                        = [0.022206331922132744, 0.560316371966525, 0.41747729611134216]
                else:
                    if (fixed_deltacp == 0.0):
                        flavor_ratios_pion \
                            = [0.3322161727912954, 0.3480562782080484, 0.319727549000656]
                        flavor_ratios_muon \
                            = [0.22220532288206532, 0.41098175587104, 0.3668129212468944]
                        flavor_ratios_neutron \
                            = [0.5522378726097557, 0.22220532288206532, 0.22555680450817917]
                    elif (fixed_deltacp == 0.5*np.pi):
                        flavor_ratios_pion \
                            = [0.3148647971219528, 0.34795485221718947, 0.33718035066085766]
                        flavor_ratios_muon \
                            = [0.1961782593780514, 0.4238431486367586, 0.3799785919851899]
                        flavor_ratios_neutron \
                            = [0.5522378726097557, 0.1961782593780514, 0.2515838680121932]
                    elif (fixed_deltacp == np.pi):
                        flavor_ratios_pion \
                            = [0.3148647971219528, 0.34795485221718947, 0.33718035066085766]
                        flavor_ratios_muon \
                            = [0.1961782593780514, 0.4238431486367586, 0.3799785919851899]
                        flavor_ratios_neutron \
                            = [0.5522378726097557, 0.1961782593780514, 0.2515838680121932]
                    else:
                        print("Invalid value of fixed_deltacp")
                        quit()
            elif (theta23 == 'lower'):
                if (fixed_deltacp is None):
                    flavor_ratios_pion \
                        = [0.32800656074697265, 0.34529257725052065, 0.3267008620025065]
                    flavor_ratios_muon \
                        = [0.2158909048155812, 0.40999341346799034, 0.3741156817164281]
                    flavor_ratios_neutron \
                        = [0.5522378726097557, 0.2158909048155812, 0.23187122257466325]
                else:
                    print("Cannot choose fixed_deltacp != None for theta = 'lower'")
                    quit()
            elif (theta23 == 'max'):
                if (fixed_deltacp is None):
                    flavor_ratios_pion \
                        = [0.31656024191908616, 0.3474455438823001, 0.3359942141986137]
                    flavor_ratios_muon \
                        = [0.19872142657375141, 0.4218076025365745, 0.37947097088967396]
                    flavor_ratios_neutron \
                        = [0.5522378726097557, 0.19872142657375141, 0.24904070081649313]
                else:
                    print("Cannot choose fixed_deltacp != None for theta = 'max'")
                    quit()
        elif (order.lower() == 'io'):
            if (theta23 == 'upper'):
                if (fixed_deltacp is None):
                    # Best-fit flavor ratios for benchmark points
                    flavor_ratios_pion \
                        = [0.3180471687709684, 0.3474505778780882, 0.3345022533509433]
                    flavor_ratios_muon \
                        = [0.2010801776404412, 0.4206357779969117, 0.378284044362647]
                    flavor_ratios_neutron \
                        = [0.5519811510320229, 0.2010801776404412, 0.24693867132753594]
                    # Best-fit flavor content of nu_1, nu_2, nu_3
                    flavor_cont_nu1 \
                        = [0.6806215638904127, 0.15226374101164764, 0.16711469509793975]
                    flavor_cont_nu2 \
                        = [0.2970466806723904, 0.2858044688729501, 0.41714885045465944]
                    flavor_cont_nu3 \
                        = [0.02233175543719699, 0.5619317901154023, 0.41573645444740065]
                else:
                    if (fixed_deltacp == 0.0):
                        flavor_ratios_pion \
                            = [0.33181153325437795, 0.34861177147789035, 0.3195766952677316]
                        flavor_ratios_muon \
                            = [0.22172672436555554, 0.41205429503405777, 0.3662189806003866]
                        flavor_ratios_neutron \
                            = [0.5519811510320229, 0.22172672436555554, 0.22629212460242162]
                    elif (fixed_deltacp == 1.0*np.pi):
                        flavor_ratios_pion \
                            = [0.29705684749618344, 0.3605546628791725, 0.34238848962464397]
                        flavor_ratios_muon \
                            = [0.16959469572826372, 0.4560346464546269, 0.3743706578171092]
                        flavor_ratios_neutron \
                            = [0.5519811510320229, 0.16959469572826372, 0.2784241532397134]
                    elif (fixed_deltacp == 1.5*np.pi):
                        flavor_ratios_pion \
                            = [0.3144341903752807, 0.34842595462849185, 0.3371398549962274]
                        flavor_ratios_muon \
                            = [0.19566071004690963, 0.424808576919283, 0.37953071303380737]
                        flavor_ratios_neutron \
                            = [0.5519811510320229, 0.19566071004690963, 0.25235813892106745]
                    else:
                        print("Invalid value of fixed_deltacp")
                        quit()
            elif (theta23 == 'lower'):
                if (fixed_deltacp is None):
                    flavor_ratios_pion \
                        = [0.34834738484344396, 0.3277980374633891, 0.32385457769316683]
                    flavor_ratios_muon \
                        = [0.2465305017491545, 0.36843180532050646, 0.38503769293033896]
                    flavor_ratios_neutron \
                        = [0.5519811510320229, 0.2465305017491545, 0.20148834721882264]
                else:
                    print("Cannot choose fixed_deltacp != None for theta = 'lower'")
                    quit()
            elif (theta23 == 'max'):
                if (fixed_deltacp is None):
                    flavor_ratios_pion \
                        = [0.33698743757165983, 0.3317785331561326, 0.33123402927220735]
                    flavor_ratios_muon \
                        = [0.22949058084147833, 0.3829225093134597, 0.3875869098450617]
                    flavor_ratios_neutron \
                        = [0.5519811510320229, 0.22949058084147833, 0.21852826812649873]
                else:
                    print("Cannot choose fixed_deltacp != None for theta = 'max'")
                    quit()

    else: # nufit_version is not None

        # if (nufit_version == 'nufit_v1.0'):
        #     if (order.lower() == 'no'):
        #         flavor_ratios_pion \
        #             = xxx
        #         flavor_ratios_muon \
        #             = xxx
        #         flavor_ratios_neutron \
        #             = xxx
        #     elif (order.lower() == 'io'):
        #         flavor_ratios_pion \
        #             = xxx
        #         flavor_ratios_muon \
        #             = xxx
        #         flavor_ratios_neutron \
        #             = xxx
        # if (nufit_version == 'nufit_v1.1'):
        #     if (order.lower() == 'no'):
        #         flavor_ratios_pion \
        #             = xxx
        #         flavor_ratios_muon \
        #             = xxx
        #         flavor_ratios_neutron \
        #             = xxx
        #     elif (order.lower() == 'io'):
        #         flavor_ratios_pion \
        #             = xxx
        #         flavor_ratios_muon \
        #             = xxx
        #         flavor_ratios_neutron \
        #             = xxx
        # if (nufit_version == 'nufit_v1.2'):
        #     if (order.lower() == 'no'):
        #         flavor_ratios_pion \
        #             = xxx
        #         flavor_ratios_muon \
        #             = xxx
        #         flavor_ratios_neutron \
        #             = xxx
        #     elif (order.lower() == 'io'):
        #         flavor_ratios_pion \
        #             = xxx
        #         flavor_ratios_muon \
        #             = xxx
        #         flavor_ratios_neutron \
        #             = xxx
        # if (nufit_version == 'nufit_v1.3'):
        #     if (order.lower() == 'no'):
        #         flavor_ratios_pion \
        #             = xxx
        #         flavor_ratios_muon \
        #             = xxx
        #         flavor_ratios_neutron \
        #             = xxx
        #     elif (order.lower() == 'io'):
        #         flavor_ratios_pion \
        #             = xxx
        #         flavor_ratios_muon \
        #             = xxx
        #         flavor_ratios_neutron \
        #             = xxx
        if (nufit_version == 'nufit_v2.0'):
            if (order.lower() == 'no'):
                flavor_ratios_pion \
                    = [0.3556613043722851, 0.3255837945733705, 0.318754901054344]
                flavor_ratios_muon \
                    = [0.25727620733858775, 0.3597375881907619, 0.38298620447064996]
                flavor_ratios_neutron \
                    = [0.5524314984396799, 0.25727620733858775, 0.19029229422173227]
            elif (order.lower() == 'io'):
                flavor_ratios_pion \
                    = [0.3085776926666689, 0.3516915412977889, 0.3397307660355422]
                flavor_ratios_muon \
                    = [0.18670502760224328, 0.43418479814556177, 0.37911017425219495]
                flavor_ratios_neutron \
                    = [0.5523230227955201, 0.18670502760224328, 0.2609719496022368]
        if (nufit_version == 'nufit_v2.1'):
            if (order.lower() == 'no'):
                flavor_ratios_pion \
                    = [0.31499629123871653, 0.3481219497568094, 0.3368817590044737]
                flavor_ratios_muon \
                    = [0.19770980806911492, 0.42332802060065666, 0.37896217133022814]
                flavor_ratios_neutron \
                    = [0.5495692575779199, 0.19770980806911492, 0.252720934352965]
            elif (order.lower() == 'io'):
                flavor_ratios_pion \
                    = [0.3090660431370772, 0.35143623465391693, 0.33949772220900554]
                flavor_ratios_muon \
                    = [0.1890301412593759, 0.4326392813511875, 0.3783305773894363]
                flavor_ratios_neutron \
                    = [0.5491378468924799, 0.1890301412593759, 0.26183201184814403]
        if (nufit_version == 'nufit_v2.2'):
            if (order.lower() == 'no'):
                flavor_ratios_pion \
                    = [0.3541798378010116, 0.3264941108322114, 0.31932605136677694]
                flavor_ratios_muon \
                    = [0.2564473535897559, 0.3615174894534392, 0.382035156956805]
                flavor_ratios_neutron \
                    = [0.5496448062235231, 0.2564473535897559, 0.19390784018672091]
            elif (order.lower() == 'io'):
                flavor_ratios_pion \
                    = [0.31156306080426477, 0.3512832938335091, 0.3371536453622259]
                flavor_ratios_muon \
                    = [0.1925869393553972, 0.4306314710725651, 0.37678158957203756]
                flavor_ratios_neutron \
                    = [0.5495153037019999, 0.1925869393553972, 0.25789775694260275]
        if (nufit_version == 'nufit_v3.0'):
            if (order.lower() == 'no'):
                flavor_ratios_pion \
                    = [0.34569828207852565, 0.33114156111530313, 0.32316015680617094]
                flavor_ratios_muon \
                    = [0.24300229079762686, 0.37521119627414135, 0.3817865129282315]
                flavor_ratios_neutron \
                    = [0.5510902646403233, 0.24300229079762686, 0.20590744456204982]
            elif (order.lower() == 'io'):
                flavor_ratios_pion \
                    = [0.31324903482344135, 0.35147907508403603, 0.3352718900925226]
                flavor_ratios_muon \
                    = [0.19439875631305437, 0.4300192344695269, 0.3755820092174188]
                flavor_ratios_neutron \
                    = [0.5509495918442153, 0.19439875631305437, 0.2546516518427304]
        if (nufit_version == 'nufit_v3.1'):
            if (order.lower() == 'no'):
                flavor_ratios_pion \
                    = [0.3053534737546712, 0.3527746876435402, 0.34187183860178844]
                flavor_ratios_muon \
                    = [0.18301214319568435, 0.43765595986746814, 0.37933189693684743]
                flavor_ratios_neutron \
                    = [0.5500361348726449, 0.18301214319568435, 0.2669517219316705]
            elif (order.lower() == 'io'):
                flavor_ratios_pion \
                    = [0.3182553364828055, 0.34681283459069917, 0.3349318289264954]
                flavor_ratios_muon \
                    = [0.20245670395710264, 0.4189908999074975, 0.37855239613540004]
                flavor_ratios_neutron \
                    = [0.5498526015342113, 0.20245670395710264, 0.24769069450868617]
        if (nufit_version == 'nufit_v3.2'):
            if (order.lower() == 'no'):
                flavor_ratios_pion \
                    = [0.3135479127848382, 0.346028878345814, 0.34042320886934785]
                flavor_ratios_muon \
                    = [0.1953631853698009, 0.4213617248338205, 0.38327508979637864]
                flavor_ratios_neutron \
                    = [0.5499173676149128, 0.1953631853698009, 0.2547194470152863]
            elif (order.lower() == 'io'):
                flavor_ratios_pion \
                    = [0.32198676464987386, 0.3423370406991942, 0.3356761946509318]
                flavor_ratios_muon \
                    = [0.2081347790102386, 0.409438171543672, 0.38242704944608913]
                flavor_ratios_neutron \
                    = [0.5496907359291444, 0.2081347790102386, 0.24217448506061715]
        if (nufit_version == 'nufit_v4.0'):
            if (order.lower() == 'no'):
                flavor_ratios_pion \
                    = [0.2988473594059442, 0.3589894623100733, 0.3421631782839824]
                flavor_ratios_muon \
                    = [0.17459388557291625, 0.45118725067865184, 0.3742188637484317]
                flavor_ratios_neutron \
                    = [0.5473543070720001, 0.17459388557291625, 0.2780518073550838]
            elif (order.lower() == 'io'):
                flavor_ratios_pion \
                    = [0.3153279935881048, 0.34979331257813706, 0.334878693833758]
                flavor_ratios_muon \
                    = [0.19943830128706735, 0.424970818223672, 0.3755908804892607]
                flavor_ratios_neutron \
                    = [0.5471073781901797, 0.19943830128706735, 0.2534543205227526]
        if (nufit_version == 'nufit_v4.1'):
            if (order.lower() == 'no'):
                flavor_ratios_pion \
                    = [0.30438521493764115, 0.3536866096566163, 0.3419281754057424]
                flavor_ratios_muon \
                    = [0.18288455868137182, 0.4390876351442386, 0.3780278061743895]
                flavor_ratios_neutron \
                    = [0.5473865274501799, 0.18288455868137182, 0.2697289138684482]
            elif (order.lower() == 'io'):
                flavor_ratios_pion \
                    = [0.3202665854332483, 0.3448720154249937, 0.33486139914175794]
                flavor_ratios_muon \
                    = [0.20682472295246254, 0.41389566166125935, 0.3792796153862782]
                flavor_ratios_neutron \
                    = [0.5471503103948199, 0.20682472295246254, 0.24602496665271745]
        if (nufit_version == 'nufit_v5.0'):
            if (order.lower() == 'no'):
                flavor_ratios_pion \
                    = [0.29829637092300665, 0.3593068266653908, 0.34239680241160236]
                flavor_ratios_muon \
                    = [0.17144024550093245, 0.45324011724762003, 0.37531963725144746]
                flavor_ratios_neutron \
                    = [0.5520086217671551, 0.17144024550093245, 0.2765511327319123]
            elif (order.lower() == 'io'):
                flavor_ratios_pion \
                    = [0.3179873235848132, 0.34751722386393324, 0.33449545255125346]
                flavor_ratios_muon \
                    = [0.20107959603042933, 0.42073603778068525, 0.37818436618888535]
                flavor_ratios_neutron \
                    = [0.5518027786935809, 0.20107959603042933, 0.24711762527598982]

    if (non_unitary_bf_fS_year == '2020'):
        flavor_ratios_pion = [0.321788, 0.339936, 0.338276]
        flavor_ratios_muon = [0.207372, 0.403323, 0.389305]
        flavor_ratios_neutron = [0.554584, 0.210965, 0.234451]
    elif (non_unitary_bf_fS_year == '2040'):
        flavor_ratios_pion = [0.304037, 0.353711, 0.342252]
        flavor_ratios_muon = [0.180976, 0.437637, 0.381387]
        flavor_ratios_neutron = [0.553604, 0.183509, 0.262886]

    lst_flavor_ratios_bf = [flavor_ratios_pion, flavor_ratios_muon, \
                            flavor_ratios_neutron]
    lst_flavor_content_bf = [flavor_cont_nu1, flavor_cont_nu2, flavor_cont_nu3]

    if show_benchmark_astro_labels:
    
        lst_color = ['salmon', 'orange', 'limegreen']
        lst_markers = ['o', 's', '^']

        for i in range(len(lst_flavor_ratios_bf)):
            tax.scatter([lst_flavor_ratios_bf[i]], \
                marker=lst_markers[i], 
                color=lst_color[i], \
                edgecolor='k', #lst_color[i], 
                s=20, #40 
                linewidths=0.5, 
                zorder=4, #4
                alpha=1.)

        # Annotate
        tax.scatter([[0.20, 0.98, 0.]], \
            marker='o', color='salmon', edgecolor='k', \
            s=40, linewidths=0.5, zorder=2, alpha=1.)
        tax.ax.annotate( r'$\pi$ decay: $\left( 1:2:0 \right)_{\rm S}$', \
            xy = (0.72,0.83), xycoords='data', color='k', \
            fontsize=10, horizontalalignment='left' )
        tax.scatter([[0.235, 0.91, 0.]], \
            marker='s', color='orange', edgecolor='k', \
            s=40, linewidths=0.5, zorder=2, alpha=1.)
        tax.ax.annotate( r'$\mu$-damped: $\left( 0:1:0 \right)_{\rm S}$', \
            xy = (0.72,0.77), xycoords='data', color='k', \
            fontsize=10, horizontalalignment='left' )
        tax.scatter([[0.267, 0.843, 0.]], \
            marker='^', color='limegreen', edgecolor='k', \
            s=40, linewidths=0.5, zorder=2, alpha=1.)
        tax.ax.annotate( r'$n$ decay: $\left( 1:0:0 \right)_{\rm S}$', \
            xy = (0.715,0.71), xycoords='data', color='k', \
            fontsize=10, horizontalalignment='left' )

    elif show_benchmark_decay_labels:
    

        # colors_010 = [None, '#7E1F86', '#B9A0CF'] # nu_1
        # colors_100 = [None, '#8C6D69', '#BFACAA'] # nu_2
        # # colors_100 = [None, '#DAB558', '#F4E9CD'] # nu_2
        # colors_120 = [None, '#05668D', '#91C4F2'] # nu_3

        # colors_100 = [None, '#DAB558', '#F4E9CD'] # nu_2
        
        lst_color = ['#05668D', '#DAB558', '#7E1F86']
        # lst_color = ['#05668D', '#8C6D69', '#7E1F86']
        lst_markers = ['o', 's', '^']

        for i in range(len(lst_flavor_content_bf)):
            tax.scatter([lst_flavor_content_bf[i]], \
                marker=lst_markers[i], 
                color=lst_color[i], \
                edgecolor='k', #lst_color[i], 
                s=20, #40 
                linewidths=0.5, 
                zorder=4, #4
                alpha=1.)

        # Annotate
        tax.scatter([[0.20, 0.98, 0.]], \
            marker='o', color=lst_color[0], edgecolor='k', \
            s=40, linewidths=0.5, zorder=2, alpha=1.)
        tax.ax.annotate( r'$\nu_1$', \
            xy = (0.72,0.83), xycoords='data', color='k', \
            fontsize=10, horizontalalignment='left' )
        tax.scatter([[0.235, 0.91, 0.]], \
            marker='s', color=lst_color[1], edgecolor='k', \
            s=40, linewidths=0.5, zorder=2, alpha=1.)
        tax.ax.annotate( r'$\nu_2$', \
            xy = (0.72,0.77), xycoords='data', color='k', \
            fontsize=10, horizontalalignment='left' )
        tax.scatter([[0.267, 0.843, 0.]], \
            marker='^', color=lst_color[2], edgecolor='k', \
            s=40, linewidths=0.5, zorder=2, alpha=1.)
        tax.ax.annotate( r'$\nu_3$', \
            xy = (0.715,0.71), xycoords='data', color='k', \
            fontsize=10, horizontalalignment='left' )


    tax.clear_matplotlib_ticks()
    tax._redraw_labels()
    ternary.plt.tight_layout()

    ternary.plt.savefig(fname_out+'.'+output_format, dpi=300)

    print("Saved "+fname_out+'.'+output_format)

    plt.close()

    return



"""
# Non-unitarity, submatrix assumption: 2020 vs. 2040

PATH_ALL = os.getcwd()+'/flavor_region_output_nonunitary_submatrix_post/all_fS/'
PATH_FIXED = os.getcwd()+'/flavor_region_output_nonunitary_submatrix_post/fixed_fS/'
if not os.path.exists(PATH_FIXED):
    os.mkdir(PATH_FIXED)

fnames_in = [
    [
        [
        PATH_ALL+'output_2020_contours',
        PATH_FIXED+'output_010_2020_contours', 
        PATH_FIXED+'output_100_2020_contours',
        PATH_FIXED+'output_120_2020_contours'
        ],
        [
        PATH_ALL+'output_2040_contours',
        PATH_FIXED+'output_010_2040_contours', 
        PATH_FIXED+'output_100_2040_contours',
        PATH_FIXED+'output_120_2040_contours'
        ]
    ]
]

fnames_out = [
    # PATH_FIXED+'output_nonunitarity_submatrix_2020_vs_2040_contours_gen2'
    PATH_FIXED+'output_nonunitarity_submatrix_2020_vs_2040_contours'
]

plot_labels = [
    # r'NO, all regions 99.7\% Cr.L.']
    # r'NO, upper $\theta_{23}$ octant'+'\n'+r'All regions 99.7\% Cr.L.']
    r'Non-unitary mixing'+'\n'+r'All regions 99.7\% C.R.']

theta23 = [
    'upper'
]

mass_order = [
   'NO'
]

year = [
    ['2020', '2040',]
]

for i, fname in enumerate(fnames_in):
    print(i)
    plot_ternary_flavor_ratios_fixed_fS_combined_2020_vs_2040(fnames_in[i], 
        fnames_out[i], plot_labels[i], order=mass_order[i], theta23=theta23[i],
        year=year[i], non_unitary_bf_fS_year='2040',
        show_benchmarks=True, show_benchmark_astro_labels=True,
        show_benchmark_decay_labels=False, show_full_exp_labels=False,
        output_format='pdf')
    plot_ternary_flavor_ratios_fixed_fS_combined_2020_vs_2040(fnames_in[i], 
        fnames_out[i], plot_labels[i], order=mass_order[i], theta23=theta23[i],
        year=year[i], non_unitary_bf_fS_year='2040',
        show_benchmarks=True, show_benchmark_astro_labels=True,
        show_benchmark_decay_labels=False, show_full_exp_labels=False,
        output_format='png')

quit()
"""



"""
# Non-unitarity: 2020 vs. 2040

PATH_ALL = os.getcwd()+'/flavor_region_output_nonunitary_post/all_fS/'
PATH_FIXED = os.getcwd()+'/flavor_region_output_nonunitary_post/fixed_fS/'
if not os.path.exists(PATH_FIXED):
    os.mkdir(PATH_FIXED)

fnames_in = [
    [
        [
        PATH_ALL+'output_2020_contours',
        PATH_FIXED+'output_010_2020_contours', 
        PATH_FIXED+'output_100_2020_contours',
        PATH_FIXED+'output_120_2020_contours'
        ],
        [
        PATH_ALL+'output_2040_contours',
        PATH_FIXED+'output_010_2040_contours', 
        PATH_FIXED+'output_100_2040_contours',
        PATH_FIXED+'output_120_2040_contours'
        ]
    ]
]

fnames_out = [
    PATH_FIXED+'output_nonunitarity_2020_vs_2040_contours'
]

plot_labels = [
    # r'NO, all regions 99.7\% Cr.L.']
    # r'NO, upper $\theta_{23}$ octant'+'\n'+r'All regions 99.7\% Cr.L.']
    r'Non-unitary mixing']

theta23 = [
    'upper'
]

mass_order = [
   'NO'
]

year = [
    ['2020', '2040',]
]

for i, fname in enumerate(fnames_in):
    print(i)
    plot_ternary_flavor_ratios_fixed_fS_combined_2020_vs_2040(fnames_in[i], 
        fnames_out[i], plot_labels[i], order=mass_order[i], theta23=theta23[i],
        year=year[i], non_unitary_bf_fS_year='2040',
        show_benchmarks=True, show_benchmark_astro_labels=True,
        show_benchmark_decay_labels=False, show_full_exp_labels=False,
        output_format='pdf')
    plot_ternary_flavor_ratios_fixed_fS_combined_2020_vs_2040(fnames_in[i], 
        fnames_out[i], plot_labels[i], order=mass_order[i], theta23=theta23[i],
        year=year[i], non_unitary_bf_fS_year='2040',
        show_benchmarks=True, show_benchmark_astro_labels=True,
        show_benchmark_decay_labels=False, show_full_exp_labels=False,
        output_format='png')

quit()
"""


# """
# Neutrino decay: 2020 vs. 2040

PATH_ALL = os.getcwd()+'/flavor_region_output_decay_JUNO_DUNE_HK_post/all_fS/'
PATH_FIXED = os.getcwd()+'/flavor_region_output_decay_JUNO_DUNE_HK_post/fixed_fS/'

# Old data
# PATH_ALL = os.getcwd()+'/flavor_region_output_decay_post/all_fS/'
# PATH_FIXED = os.getcwd()+'/flavor_region_output_decay_post/fixed_fS/'

if not os.path.exists(PATH_FIXED):
    os.mkdir(PATH_FIXED)

fnames_in = [
    [
        [
        PATH_ALL+'output_NO_upper_NUFIT_contours',
        PATH_FIXED+'output_001_NO_upper_NUFIT_contours', 
        PATH_FIXED+'output_010_NO_upper_NUFIT_contours',
        PATH_FIXED+'output_100_NO_upper_NUFIT_contours'
        ],
        [
        PATH_ALL+'output_NO_upper_JUNO_DUNE_HK_contours',
        PATH_FIXED+'output_001_NO_upper_JUNO_DUNE_HK_contours', 
        PATH_FIXED+'output_010_NO_upper_JUNO_DUNE_HK_contours',
        PATH_FIXED+'output_100_NO_upper_JUNO_DUNE_HK_contours'
        ]
        # Old data
        # [
        # PATH_ALL+'output_NO_upper_JUNO_DUNE_HK_contours',
        # PATH_FIXED+'output_001_NO_upper_JUNO_DUNE_HK_contours', 
        # PATH_FIXED+'output_010_NO_upper_JUNO_DUNE_HK_contours',
        # PATH_FIXED+'output_100_NO_upper_JUNO_DUNE_HK_contours'
        # ]
    ]
]

fnames_out = [
    # PATH_FIXED+'output_decay_NO_upper_2020_vs_2040_contours_gen2'
    PATH_FIXED+'output_decay_NO_upper_2020_vs_2040_contours_final'
    # PATH_FIXED+'output_decay_NO_upper_2020_vs_2040_contours_all_nui_nocont'
]

plot_labels = [
    # r'NO, all regions 99.7\% Cr.L.']
    # r'NO, upper $\theta_{23}$ octant'+'\n'+r'All regions 99.7\% Cr.L.']
    # r'$\nu$ decay, NO'+'\n'+r'All regions 99.7\% C.R.']
    r'$\nu$ decay'+'\n'+r'All regions 99.7\% C.R.']

theta23 = [
    'upper'
]

mass_order = [
   'NO'
]

year = [
    ['2020', '2040',]
]

for i, fname in enumerate(fnames_in):
    print(i)
    plot_ternary_flavor_ratios_fixed_fS_combined_2020_vs_2040(fnames_in[i], 
        fnames_out[i], plot_labels[i], order=mass_order[i], theta23=theta23[i],
        year=year[i], show_benchmarks=True, show_benchmark_astro_labels=False,
        show_benchmark_decay_labels=True, output_format='pdf')
    plot_ternary_flavor_ratios_fixed_fS_combined_2020_vs_2040(fnames_in[i], 
        fnames_out[i], plot_labels[i], order=mass_order[i], theta23=theta23[i],
        year=year[i], show_benchmarks=True, show_benchmark_astro_labels=False,
        show_benchmark_decay_labels=True, output_format='png')

quit()
# """


"""
# Fixed flavor composition at the source, combined plots, 2020 vs. 2040

PATH_ALL = os.getcwd()+'/flavor_region_output_JUNO_DUNE_HK_post/all_fS/'
PATH_FIXED = os.getcwd()+'/flavor_region_output_JUNO_DUNE_HK_post/fixed_fS/'

# Old data
# PATH_ALL = os.getcwd()+'/flavor_region_output_post/'
# PATH_FIXED = os.getcwd()+'/flavor_region_fixed_source_output_post/'

if not os.path.exists(PATH_FIXED):
    os.mkdir(PATH_FIXED)

fnames_in = [
    [
        [
        PATH_ALL+'output_NO_upper_NUFIT_contours', 
        PATH_FIXED+'output_010_NO_upper_NUFIT_contours', 
        PATH_FIXED+'output_100_NO_upper_NUFIT_contours',
        PATH_FIXED+'output_120_NO_upper_NUFIT_contours'
        ],
        [
        PATH_ALL+'output_NO_upper_JUNO_DUNE_HK_contours', 
        PATH_FIXED+'output_010_NO_upper_JUNO_DUNE_HK_contours', 
        PATH_FIXED+'output_100_NO_upper_JUNO_DUNE_HK_contours',
        PATH_FIXED+'output_120_NO_upper_JUNO_DUNE_HK_contours'
        ]
        # Old data
        # [
        # PATH_ALL+'output_NO_upper_DUNE_JUNO_contours', 
        # PATH_FIXED+'output_010_NO_upper_JUNO_DUNE_contours', 
        # PATH_FIXED+'output_100_NO_upper_JUNO_DUNE_contours',
        # PATH_FIXED+'output_120_NO_upper_JUNO_DUNE_contours'
        # ]
    ]
]

fnames_out = [
    # PATH_FIXED+'output_NO_upper_2020_vs_2040_contours_gen2'
    PATH_FIXED+'output_NO_upper_2020_vs_2040_contours_final'
    # PATH_FIXED+'output_NO_upper_2020_vs_2040_contours_all_fS_only_2020_2040_allcont'
]

plot_labels = [
    # r'NO, all regions 99.7\% Cr.L.']
    # r'NO, upper $\theta_{23}$ octant'+'\n'+r'All regions 99.7\% Cr.L.']
    r'Standard oscillations, NO'+'\n'+r'All regions 99.7\% C.R.']

theta23 = [
    'upper'
]

mass_order = [
   'NO'
]

year = [
    ['2020', '2040',]
]

for i, fname in enumerate(fnames_in):
    print(i)
    plot_ternary_flavor_ratios_fixed_fS_combined_2020_vs_2040(fnames_in[i], 
        fnames_out[i], plot_labels[i], order=mass_order[i], theta23=theta23[i],
        year=year[i], output_format='pdf')
    plot_ternary_flavor_ratios_fixed_fS_combined_2020_vs_2040(fnames_in[i], 
        fnames_out[i], plot_labels[i], order=mass_order[i], theta23=theta23[i],
        year=year[i], output_format='png')

quit()
"""


"""
# Fixed flavor composition at the source, combined plots, fixed delta_CP (v1)

PATH_ALL = os.getcwd()+'/flavor_region_output_fixed_delta_post/all_fS/'
PATH_FIXED = os.getcwd()+'/flavor_region_output_fixed_delta_post/fixed_fS/'
if not os.path.exists(PATH_FIXED):
    os.mkdir(PATH_FIXED)

fnames_in = [
    [
    PATH_ALL+'output_IO_upper_delta0p000000pi_DUNE_contours',
    PATH_FIXED+'output_010_IO_upper_delta0p000000pi_DUNE_contours',
    PATH_FIXED+'output_100_IO_upper_delta0p000000pi_DUNE_contours',
    PATH_FIXED+'output_120_IO_upper_delta0p000000pi_DUNE_contours',
    ],
    [
    PATH_ALL+'output_IO_upper_delta0p000000pi_HK_contours',
    PATH_FIXED+'output_010_IO_upper_delta0p000000pi_HK_contours',
    PATH_FIXED+'output_100_IO_upper_delta0p000000pi_HK_contours',
    PATH_FIXED+'output_120_IO_upper_delta0p000000pi_HK_contours',
    ],
    [
    PATH_ALL+'output_IO_upper_delta0p000000pi_JUNO_DUNE_contours',
    PATH_FIXED+'output_010_IO_upper_delta0p000000pi_JUNO_DUNE_contours',
    PATH_FIXED+'output_100_IO_upper_delta0p000000pi_JUNO_DUNE_contours',
    PATH_FIXED+'output_120_IO_upper_delta0p000000pi_JUNO_DUNE_contours',
    ],
    [
    PATH_ALL+'output_IO_upper_delta0p000000pi_JUNO_HK_contours',
    PATH_FIXED+'output_010_IO_upper_delta0p000000pi_JUNO_HK_contours',
    PATH_FIXED+'output_100_IO_upper_delta0p000000pi_JUNO_HK_contours',
    PATH_FIXED+'output_120_IO_upper_delta0p000000pi_JUNO_HK_contours',
    ],
    [
    PATH_ALL+'output_IO_upper_delta1p000000pi_DUNE_contours',
    PATH_FIXED+'output_010_IO_upper_delta1p000000pi_DUNE_contours',
    PATH_FIXED+'output_100_IO_upper_delta1p000000pi_DUNE_contours',
    PATH_FIXED+'output_120_IO_upper_delta1p000000pi_DUNE_contours',
    ],
    [
    PATH_ALL+'output_IO_upper_delta1p000000pi_HK_contours',
    PATH_FIXED+'output_010_IO_upper_delta1p000000pi_HK_contours',
    PATH_FIXED+'output_100_IO_upper_delta1p000000pi_HK_contours',
    PATH_FIXED+'output_120_IO_upper_delta1p000000pi_HK_contours',
    ],
    [
    PATH_ALL+'output_IO_upper_delta1p000000pi_JUNO_DUNE_contours',
    PATH_FIXED+'output_010_IO_upper_delta1p000000pi_JUNO_DUNE_contours',
    PATH_FIXED+'output_100_IO_upper_delta1p000000pi_JUNO_DUNE_contours',
    PATH_FIXED+'output_120_IO_upper_delta1p000000pi_JUNO_DUNE_contours',
    ],
    [
    PATH_ALL+'output_IO_upper_delta1p000000pi_JUNO_HK_contours',
    PATH_FIXED+'output_010_IO_upper_delta1p000000pi_JUNO_HK_contours',
    PATH_FIXED+'output_100_IO_upper_delta1p000000pi_JUNO_HK_contours',
    PATH_FIXED+'output_120_IO_upper_delta1p000000pi_JUNO_HK_contours',
    ],
    [
    PATH_ALL+'output_IO_upper_delta1p500000pi_DUNE_contours',
    PATH_FIXED+'output_010_IO_upper_delta1p500000pi_DUNE_contours',
    PATH_FIXED+'output_100_IO_upper_delta1p500000pi_DUNE_contours',
    PATH_FIXED+'output_120_IO_upper_delta1p500000pi_DUNE_contours',
    ],
    [
    PATH_ALL+'output_IO_upper_delta1p500000pi_HK_contours',
    PATH_FIXED+'output_010_IO_upper_delta1p500000pi_HK_contours',
    PATH_FIXED+'output_100_IO_upper_delta1p500000pi_HK_contours',
    PATH_FIXED+'output_120_IO_upper_delta1p500000pi_HK_contours',
    ],
    [
    PATH_ALL+'output_IO_upper_delta1p500000pi_JUNO_DUNE_contours',
    PATH_FIXED+'output_010_IO_upper_delta1p500000pi_JUNO_DUNE_contours',
    PATH_FIXED+'output_100_IO_upper_delta1p500000pi_JUNO_DUNE_contours',
    PATH_FIXED+'output_120_IO_upper_delta1p500000pi_JUNO_DUNE_contours',
    ],
    [
    PATH_ALL+'output_IO_upper_delta1p500000pi_JUNO_HK_contours',
    PATH_FIXED+'output_010_IO_upper_delta1p500000pi_JUNO_HK_contours',
    PATH_FIXED+'output_100_IO_upper_delta1p500000pi_JUNO_HK_contours',
    PATH_FIXED+'output_120_IO_upper_delta1p500000pi_JUNO_HK_contours',
    ],
    [
    PATH_ALL+'output_NO_upper_delta0p000000pi_DUNE_contours',
    PATH_FIXED+'output_010_NO_upper_delta0p000000pi_DUNE_contours',
    PATH_FIXED+'output_100_NO_upper_delta0p000000pi_DUNE_contours',
    PATH_FIXED+'output_120_NO_upper_delta0p000000pi_DUNE_contours',
    ],
    [
    PATH_ALL+'output_NO_upper_delta0p000000pi_HK_contours',
    PATH_FIXED+'output_010_NO_upper_delta0p000000pi_HK_contours',
    PATH_FIXED+'output_100_NO_upper_delta0p000000pi_HK_contours',
    PATH_FIXED+'output_120_NO_upper_delta0p000000pi_HK_contours',
    ],
    [
    PATH_ALL+'output_NO_upper_delta0p000000pi_JUNO_DUNE_contours',
    PATH_FIXED+'output_010_NO_upper_delta0p000000pi_JUNO_DUNE_contours',
    PATH_FIXED+'output_100_NO_upper_delta0p000000pi_JUNO_DUNE_contours',
    PATH_FIXED+'output_120_NO_upper_delta0p000000pi_JUNO_DUNE_contours',
    ],
    [
    PATH_ALL+'output_NO_upper_delta0p000000pi_JUNO_HK_contours',
    PATH_FIXED+'output_010_NO_upper_delta0p000000pi_JUNO_HK_contours',
    PATH_FIXED+'output_100_NO_upper_delta0p000000pi_JUNO_HK_contours',
    PATH_FIXED+'output_120_NO_upper_delta0p000000pi_JUNO_HK_contours',
    ],
    [
    PATH_ALL+'output_NO_upper_delta0p500000pi_DUNE_contours',
    PATH_FIXED+'output_010_NO_upper_delta0p500000pi_DUNE_contours',
    PATH_FIXED+'output_100_NO_upper_delta0p500000pi_DUNE_contours',
    PATH_FIXED+'output_120_NO_upper_delta0p500000pi_DUNE_contours',
    ],
    [
    PATH_ALL+'output_NO_upper_delta0p500000pi_HK_contours',
    PATH_FIXED+'output_010_NO_upper_delta0p500000pi_HK_contours',
    PATH_FIXED+'output_100_NO_upper_delta0p500000pi_HK_contours',
    PATH_FIXED+'output_120_NO_upper_delta0p500000pi_HK_contours',
    ],
    [
    PATH_ALL+'output_NO_upper_delta0p500000pi_JUNO_DUNE_contours',
    PATH_FIXED+'output_010_NO_upper_delta0p500000pi_JUNO_DUNE_contours',
    PATH_FIXED+'output_100_NO_upper_delta0p500000pi_JUNO_DUNE_contours',
    PATH_FIXED+'output_120_NO_upper_delta0p500000pi_JUNO_DUNE_contours',
    ],
    [
    PATH_ALL+'output_NO_upper_delta0p500000pi_JUNO_HK_contours',
    PATH_FIXED+'output_010_NO_upper_delta0p500000pi_JUNO_HK_contours',
    PATH_FIXED+'output_100_NO_upper_delta0p500000pi_JUNO_HK_contours',
    PATH_FIXED+'output_120_NO_upper_delta0p500000pi_JUNO_HK_contours',
    ],
    [
    PATH_ALL+'output_NO_upper_delta1p000000pi_DUNE_contours',
    PATH_FIXED+'output_010_NO_upper_delta1p000000pi_DUNE_contours',
    PATH_FIXED+'output_100_NO_upper_delta1p000000pi_DUNE_contours',
    PATH_FIXED+'output_120_NO_upper_delta1p000000pi_DUNE_contours',
    ],
    [
    PATH_ALL+'output_NO_upper_delta1p000000pi_HK_contours',
    PATH_FIXED+'output_010_NO_upper_delta1p000000pi_HK_contours',
    PATH_FIXED+'output_100_NO_upper_delta1p000000pi_HK_contours',
    PATH_FIXED+'output_120_NO_upper_delta1p000000pi_HK_contours',
    ],
    [
    PATH_ALL+'output_NO_upper_delta1p000000pi_JUNO_DUNE_contours',
    PATH_FIXED+'output_010_NO_upper_delta1p000000pi_JUNO_DUNE_contours',
    PATH_FIXED+'output_100_NO_upper_delta1p000000pi_JUNO_DUNE_contours',
    PATH_FIXED+'output_120_NO_upper_delta1p000000pi_JUNO_DUNE_contours',
    ],
    [
    PATH_ALL+'output_NO_upper_delta1p000000pi_JUNO_HK_contours',
    PATH_FIXED+'output_010_NO_upper_delta1p000000pi_JUNO_HK_contours',
    PATH_FIXED+'output_100_NO_upper_delta1p000000pi_JUNO_HK_contours',
    PATH_FIXED+'output_120_NO_upper_delta1p000000pi_JUNO_HK_contours',
    ]
]

fnames_out = [
    PATH_FIXED+'output_010_100_120_IO_upper_delta0p000000pi_DUNE_contours', 
    PATH_FIXED+'output_010_100_120_IO_upper_delta0p000000pi_HK_contours', 
    PATH_FIXED+'output_010_100_120_IO_upper_delta0p000000pi_JUNO_DUNE_contours', 
    PATH_FIXED+'output_010_100_120_IO_upper_delta0p000000pi_JUNO_HK_contours', 
    PATH_FIXED+'output_010_100_120_IO_upper_delta1p000000pi_DUNE_contours', 
    PATH_FIXED+'output_010_100_120_IO_upper_delta1p000000pi_HK_contours', 
    PATH_FIXED+'output_010_100_120_IO_upper_delta1p000000pi_JUNO_DUNE_contours', 
    PATH_FIXED+'output_010_100_120_IO_upper_delta1p000000pi_JUNO_HK_contours', 
    PATH_FIXED+'output_010_100_120_IO_upper_delta1p500000pi_DUNE_contours', 
    PATH_FIXED+'output_010_100_120_IO_upper_delta1p500000pi_HK_contours', 
    PATH_FIXED+'output_010_100_120_IO_upper_delta1p500000pi_JUNO_DUNE_contours', 
    PATH_FIXED+'output_010_100_120_IO_upper_delta1p500000pi_JUNO_HK_contours', 
    PATH_FIXED+'output_010_100_120_NO_upper_delta0p000000pi_DUNE_contours', 
    PATH_FIXED+'output_010_100_120_NO_upper_delta0p000000pi_HK_contours', 
    PATH_FIXED+'output_010_100_120_NO_upper_delta0p000000pi_JUNO_DUNE_contours', 
    PATH_FIXED+'output_010_100_120_NO_upper_delta0p000000pi_JUNO_HK_contours', 
    PATH_FIXED+'output_010_100_120_NO_upper_delta0p500000pi_DUNE_contours', 
    PATH_FIXED+'output_010_100_120_NO_upper_delta0p500000pi_HK_contours', 
    PATH_FIXED+'output_010_100_120_NO_upper_delta0p500000pi_JUNO_DUNE_contours', 
    PATH_FIXED+'output_010_100_120_NO_upper_delta0p500000pi_JUNO_HK_contours', 
    PATH_FIXED+'output_010_100_120_NO_upper_delta1p000000pi_DUNE_contours', 
    PATH_FIXED+'output_010_100_120_NO_upper_delta1p000000pi_HK_contours', 
    PATH_FIXED+'output_010_100_120_NO_upper_delta1p000000pi_JUNO_DUNE_contours', 
    PATH_FIXED+'output_010_100_120_NO_upper_delta1p000000pi_JUNO_HK_contours'
]

plot_labels = [
    r'IO, upper $\theta_{23}$ octant, $\delta_{\rm CP} = 0$,'+'\n'+r'DUNE',
    r'IO, upper $\theta_{23}$ octant, $\delta_{\rm CP} = 0$,'+'\n'+r'HK',
    r'IO, upper $\theta_{23}$ octant, $\delta_{\rm CP} = 0$,'+'\n'+r'JUNO + DUNE',
    r'IO, upper $\theta_{23}$ octant, $\delta_{\rm CP} = 0$,'+'\n'+r'JUNO + HK',
    r'IO, upper $\theta_{23}$ octant, $\delta_{\rm CP} = \pi$,'+'\n'+r'DUNE',
    r'IO, upper $\theta_{23}$ octant, $\delta_{\rm CP} = \pi$,'+'\n'+r'HK',
    r'IO, upper $\theta_{23}$ octant, $\delta_{\rm CP} = \pi$,'+'\n'+r'JUNO + DUNE',
    r'IO, upper $\theta_{23}$ octant, $\delta_{\rm CP} = \pi$,'+'\n'+r'JUNO + HK',
    r'IO, upper $\theta_{23}$ octant, $\delta_{\rm CP} = \frac{3\pi}{2}$,'+'\n'+r'DUNE',
    r'IO, upper $\theta_{23}$ octant, $\delta_{\rm CP} = \frac{3\pi}{2}$,'+'\n'+r'HK',
    r'IO, upper $\theta_{23}$ octant, $\delta_{\rm CP} = \frac{3\pi}{2}$,'+'\n'+r'JUNO + DUNE',
    r'IO, upper $\theta_{23}$ octant, $\delta_{\rm CP} = \frac{3\pi}{2}$,'+'\n'+r'JUNO + HK',
    r'NO, upper $\theta_{23}$ octant, $\delta_{\rm CP} = 0$,'+'\n'+r'DUNE',
    r'NO, upper $\theta_{23}$ octant, $\delta_{\rm CP} = 0$,'+'\n'+r'HK',
    r'NO, upper $\theta_{23}$ octant, $\delta_{\rm CP} = 0$,'+'\n'+r'JUNO + DUNE',
    r'NO, upper $\theta_{23}$ octant, $\delta_{\rm CP} = 0$,'+'\n'+r'JUNO + HK',
    r'NO, upper $\theta_{23}$ octant, $\delta_{\rm CP} = \frac{\pi}{2}$,'+'\n'+r'DUNE',
    r'NO, upper $\theta_{23}$ octant, $\delta_{\rm CP} = \frac{\pi}{2}$,'+'\n'+r'HK',
    r'NO, upper $\theta_{23}$ octant, $\delta_{\rm CP} = \frac{\pi}{2}$,'+'\n'+r'JUNO + DUNE',
    r'NO, upper $\theta_{23}$ octant, $\delta_{\rm CP} = \frac{\pi}{2}$,'+'\n'+r'JUNO + HK',
    r'NO, upper $\theta_{23}$ octant, $\delta_{\rm CP} = \pi$,'+'\n'+r'DUNE',
    r'NO, upper $\theta_{23}$ octant, $\delta_{\rm CP} = \pi$,'+'\n'+r'HK',
    r'NO, upper $\theta_{23}$ octant, $\delta_{\rm CP} = \pi$,'+'\n'+r'JUNO + DUNE',
    r'NO, upper $\theta_{23}$ octant, $\delta_{\rm CP} = \pi$,'+'\n'+r'JUNO + HK'
]

theta23 = [
    'upper', 
    'upper',
    'upper', 
    'upper',
    'upper', 
    'upper',
    'upper', 
    'upper',
    'upper', 
    'upper',
    'upper', 
    'upper',
    'upper', 
    'upper',
    'upper', 
    'upper',
    'upper', 
    'upper',
    'upper', 
    'upper',
    'upper', 
    'upper',
    'upper', 
    'upper'
]

mass_order = [
   'IO', 
   'IO', 
   'IO', 
   'IO', 
   'IO', 
   'IO', 
   'IO', 
   'IO', 
   'IO', 
   'IO', 
   'IO', 
   'IO', 
   'NO', 
   'NO', 
   'NO', 
   'NO', 
   'NO', 
   'NO', 
   'NO', 
   'NO', 
   'NO', 
   'NO', 
   'NO', 
   'NO'
]

year = [
    '2040',
    '2030',
    '2040',
    '2030',
    '2040',
    '2030',
    '2040',
    '2030',
    '2040',
    '2030',
    '2040',
    '2030',
    '2040',
    '2030',
    '2040',
    '2030',
    '2040',
    '2030',
    '2040',
    '2030',
    '2040',
    '2030',
    '2040',
    '2030'
]

fixed_deltacp = [
    0.0,
    0.0,
    0.0,
    0.0,
    np.pi,        
    np.pi,        
    np.pi,        
    np.pi,        
    1.5*np.pi,        
    1.5*np.pi,        
    1.5*np.pi,        
    1.5*np.pi,        
    0.0,
    0.0,
    0.0,
    0.0,
    0.5*np.pi,        
    0.5*np.pi,        
    0.5*np.pi,        
    0.5*np.pi,        
    np.pi,        
    np.pi,        
    np.pi,        
    np.pi,        
]

for i, fname in enumerate(fnames_in):
    print(i)
    plot_ternary_flavor_ratios_fixed_fS_combined_v1(fnames_in[i], fnames_out[i], 
        plot_labels[i], order=mass_order[i], theta23=theta23[i],
        year=year[i], fixed_deltacp=fixed_deltacp[i], output_format='png')
"""


"""
# Fixed flavor composition at the source, combined plots, different NuFit versions (v1)

PATH_ALL = os.getcwd()+'/flavor_region_output_nufit_versions_post/all_fS/'
PATH_FIXED = os.getcwd()+'/flavor_region_output_nufit_versions_post/fixed_fS/'
if not os.path.exists(PATH_FIXED):
    os.mkdir(PATH_FIXED)

fnames_in = [
    [
    PATH_ALL+'output_IO_v2p000000_contours',
    PATH_FIXED+'output_010_IO_v2p000000_contours', 
    PATH_FIXED+'output_100_IO_v2p000000_contours',
    PATH_FIXED+'output_120_IO_v2p000000_contours'
    ],
    [
    PATH_ALL+'output_IO_v2p100000_contours',
    PATH_FIXED+'output_010_IO_v2p100000_contours', 
    PATH_FIXED+'output_100_IO_v2p100000_contours',
    PATH_FIXED+'output_120_IO_v2p100000_contours'
    ],
    [
    PATH_ALL+'output_IO_v2p200000_contours',
    PATH_FIXED+'output_010_IO_v2p200000_contours', 
    PATH_FIXED+'output_100_IO_v2p200000_contours',
    PATH_FIXED+'output_120_IO_v2p200000_contours'
    ],
    [
    PATH_ALL+'output_IO_v3p000000_contours',
    PATH_FIXED+'output_010_IO_v3p000000_contours', 
    PATH_FIXED+'output_100_IO_v3p000000_contours',
    PATH_FIXED+'output_120_IO_v3p000000_contours'
    ],
    [
    PATH_ALL+'output_IO_v3p100000_contours',
    PATH_FIXED+'output_010_IO_v3p100000_contours', 
    PATH_FIXED+'output_100_IO_v3p100000_contours',
    PATH_FIXED+'output_120_IO_v3p100000_contours'
    ],
    [
    PATH_ALL+'output_IO_v3p200000_contours',
    PATH_FIXED+'output_010_IO_v3p200000_contours', 
    PATH_FIXED+'output_100_IO_v3p200000_contours',
    PATH_FIXED+'output_120_IO_v3p200000_contours'
    ],
    [
    PATH_ALL+'output_IO_v4p000000_contours',
    PATH_FIXED+'output_010_IO_v4p000000_contours', 
    PATH_FIXED+'output_100_IO_v4p000000_contours',
    PATH_FIXED+'output_120_IO_v4p000000_contours'
    ],
    [
    PATH_ALL+'output_IO_v4p100000_contours',
    PATH_FIXED+'output_010_IO_v4p100000_contours', 
    PATH_FIXED+'output_100_IO_v4p100000_contours',
    PATH_FIXED+'output_120_IO_v4p100000_contours'
    ],
    [
    PATH_ALL+'output_IO_v5p000000_contours',
    PATH_FIXED+'output_010_IO_v5p000000_contours', 
    PATH_FIXED+'output_100_IO_v5p000000_contours',
    PATH_FIXED+'output_120_IO_v5p000000_contours'
    ],
    [
    PATH_ALL+'output_NO_v2p000000_contours',
    PATH_FIXED+'output_010_NO_v2p000000_contours', 
    PATH_FIXED+'output_100_NO_v2p000000_contours',
    PATH_FIXED+'output_120_NO_v2p000000_contours'
    ],
    [
    PATH_ALL+'output_NO_v2p100000_contours',
    PATH_FIXED+'output_010_NO_v2p100000_contours', 
    PATH_FIXED+'output_100_NO_v2p100000_contours',
    PATH_FIXED+'output_120_NO_v2p100000_contours'
    ],
    [
    PATH_ALL+'output_NO_v2p200000_contours',
    PATH_FIXED+'output_010_NO_v2p200000_contours', 
    PATH_FIXED+'output_100_NO_v2p200000_contours',
    PATH_FIXED+'output_120_NO_v2p200000_contours'
    ],
    [
    PATH_ALL+'output_NO_v3p000000_contours',
    PATH_FIXED+'output_010_NO_v3p000000_contours', 
    PATH_FIXED+'output_100_NO_v3p000000_contours',
    PATH_FIXED+'output_120_NO_v3p000000_contours'
    ],
    [
    PATH_ALL+'output_NO_v3p100000_contours',
    PATH_FIXED+'output_010_NO_v3p100000_contours', 
    PATH_FIXED+'output_100_NO_v3p100000_contours',
    PATH_FIXED+'output_120_NO_v3p100000_contours'
    ],
    [
    PATH_ALL+'output_NO_v3p200000_contours',
    PATH_FIXED+'output_010_NO_v3p200000_contours', 
    PATH_FIXED+'output_100_NO_v3p200000_contours',
    PATH_FIXED+'output_120_NO_v3p200000_contours'
    ],
    [
    PATH_ALL+'output_NO_v4p000000_contours',
    PATH_FIXED+'output_010_NO_v4p000000_contours', 
    PATH_FIXED+'output_100_NO_v4p000000_contours',
    PATH_FIXED+'output_120_NO_v4p000000_contours'
    ],
    [
    PATH_ALL+'output_NO_v4p100000_contours',
    PATH_FIXED+'output_010_NO_v4p100000_contours', 
    PATH_FIXED+'output_100_NO_v4p100000_contours',
    PATH_FIXED+'output_120_NO_v4p100000_contours'
    ],
    [
    PATH_ALL+'output_NO_v5p000000_contours',
    PATH_FIXED+'output_010_NO_v5p000000_contours', 
    PATH_FIXED+'output_100_NO_v5p000000_contours',
    PATH_FIXED+'output_120_NO_v5p000000_contours'
    ]
]

fnames_out = [
    PATH_FIXED+'output_010_100_120_IO_v2p000000_contours', 
    PATH_FIXED+'output_010_100_120_IO_v2p100000_contours', 
    PATH_FIXED+'output_010_100_120_IO_v2p200000_contours', 
    PATH_FIXED+'output_010_100_120_IO_v3p000000_contours', 
    PATH_FIXED+'output_010_100_120_IO_v3p100000_contours', 
    PATH_FIXED+'output_010_100_120_IO_v3p200000_contours', 
    PATH_FIXED+'output_010_100_120_IO_v4p000000_contours', 
    PATH_FIXED+'output_010_100_120_IO_v4p100000_contours', 
    PATH_FIXED+'output_010_100_120_IO_v5p000000_contours', 
    PATH_FIXED+'output_010_100_120_NO_v2p000000_contours', 
    PATH_FIXED+'output_010_100_120_NO_v2p100000_contours', 
    PATH_FIXED+'output_010_100_120_NO_v2p200000_contours', 
    PATH_FIXED+'output_010_100_120_NO_v3p000000_contours', 
    PATH_FIXED+'output_010_100_120_NO_v3p100000_contours', 
    PATH_FIXED+'output_010_100_120_NO_v3p200000_contours', 
    PATH_FIXED+'output_010_100_120_NO_v4p000000_contours', 
    PATH_FIXED+'output_010_100_120_NO_v4p100000_contours', 
    PATH_FIXED+'output_010_100_120_NO_v5p000000_contours'
]

plot_labels = [
    # r'IO'+'\n'+r'NuFit 1.0', 
    # r'IO'+'\n'+r'NuFit 1.1',
    # r'IO'+'\n'+r'NuFit 1.2',
    # r'IO'+'\n'+r'NuFit 1.3',
    r'IO'+'\n'+r'NuFit 2.0',
    r'IO'+'\n'+r'NuFit 2.1',
    r'IO'+'\n'+r'NuFit 2.2',
    r'IO'+'\n'+r'NuFit 3.0',
    r'IO'+'\n'+r'NuFit 3.1',
    r'IO'+'\n'+r'NuFit 3.2',
    r'IO'+'\n'+r'NuFit 4.0',
    r'IO'+'\n'+r'NuFit 4.1',
    r'IO'+'\n'+r'NuFit 5.0',
    # r'NO'+'\n'+r'NuFit 1.0', 
    # r'NO'+'\n'+r'NuFit 1.1',
    # r'NO'+'\n'+r'NuFit 1.2',
    # r'NO'+'\n'+r'NuFit 1.3',
    r'NO'+'\n'+r'NuFit 2.0',
    r'NO'+'\n'+r'NuFit 2.1',
    r'NO'+'\n'+r'NuFit 2.2',
    r'NO'+'\n'+r'NuFit 3.0',
    r'NO'+'\n'+r'NuFit 3.1',
    r'NO'+'\n'+r'NuFit 3.2',
    r'NO'+'\n'+r'NuFit 4.0',
    r'NO'+'\n'+r'NuFit 4.1',
    r'NO'+'\n'+r'NuFit 5.0'
]

nufit_version = [
    # 'nufit_v1.0',
    # 'nufit_v1.1',
    # 'nufit_v1.2',
    # 'nufit_v1.3',
    'nufit_v2.0',
    'nufit_v2.1',
    'nufit_v2.2',
    'nufit_v3.0',
    'nufit_v3.1',
    'nufit_v3.2',
    'nufit_v4.0',
    'nufit_v4.1',
    'nufit_v5.0',
    # 'nufit_v1.0',
    # 'nufit_v1.1',
    # 'nufit_v1.2',
    # 'nufit_v1.3',
    'nufit_v2.0',
    'nufit_v2.1',
    'nufit_v2.2',
    'nufit_v3.0',
    'nufit_v3.1',
    'nufit_v3.2',
    'nufit_v4.0',
    'nufit_v4.1',
    'nufit_v5.0'
]

mass_order = [
    # 'IO',
    # 'IO',
    # 'IO',
    # 'IO',
    'IO',
    'IO',
    'IO',
    'IO',
    'IO',
    'IO',
    'IO',
    'IO',
    'IO',
    'NO',
    # 'NO',
    # 'NO',
    # 'NO',
    # 'NO',
    'NO',
    'NO',
    'NO',
    'NO',
    'NO',
    'NO',
    'NO',
    'NO']


for i, fname in enumerate(fnames_in):
    print(i)
    plot_ternary_flavor_ratios_fixed_fS_combined_v1(fnames_in[i], fnames_out[i], 
        plot_labels[i], order=mass_order[i], theta23=None,
        year=None, nufit_version=nufit_version[i], output_format='png')
"""


"""
# Non-unitarity: fixed flavor composition at the source, combined plots (v1)

PATH_ALL = os.getcwd()+'/flavor_region_output_nonunitary_post/all_fS/'
PATH_FIXED = os.getcwd()+'/flavor_region_output_nonunitary_post/fixed_fS/'
if not os.path.exists(PATH_FIXED):
    os.mkdir(PATH_FIXED)

fnames_in = [
    [
    PATH_ALL+'output_2020_contours',
    PATH_FIXED+'output_010_2020_contours', 
    PATH_FIXED+'output_100_2020_contours',
    PATH_FIXED+'output_120_2020_contours'
    ],
    [
    PATH_ALL+'output_2040_contours',
    PATH_FIXED+'output_010_2040_contours', 
    PATH_FIXED+'output_100_2040_contours',
    PATH_FIXED+'output_120_2040_contours'
    ],
]

fnames_out = [
    PATH_FIXED+'output_decay_010_100_120_nonunitary_2020_contours', 
    PATH_FIXED+'output_decay_010_100_120_nonunitary_2040_contours', 
]

plot_labels = [
    r'Non-unitarity, 2020',
    r'Non-unitarity, 2040'
]

theta23 = [
    'upper', 
    'upper', 
]

mass_order = [
   'NO', 
   'NO'
]

year = [
    '2020',
    '2040'
]

for i, fname in enumerate(fnames_in):
    print(i)
    plot_ternary_flavor_ratios_fixed_fS_combined_v1(fnames_in[i], fnames_out[i], 
        plot_labels[i], order=mass_order[i], theta23=theta23[i],
        year=year[i], non_unitary_bf_fS_year=year[i],
        show_benchmarks=True, show_benchmark_astro_labels=True,
        show_benchmark_decay_labels=False,
        output_format='png')
"""


"""
# Neutrino decay: fixed flavor composition at the source, combined plots (v1)

PATH_ALL = os.getcwd()+'/flavor_region_output_decay_post/all_fS/'
PATH_FIXED = os.getcwd()+'/flavor_region_output_decay_post/fixed_fS/'
if not os.path.exists(PATH_FIXED):
    os.mkdir(PATH_FIXED)

fnames_in = [
    [
    PATH_ALL+'output_IO_upper_JUNO_HK_contours',
    PATH_FIXED+'output_001_IO_upper_JUNO_HK_contours', 
    PATH_FIXED+'output_010_IO_upper_JUNO_HK_contours',
    PATH_FIXED+'output_100_IO_upper_JUNO_HK_contours'
    ],
    [
    PATH_ALL+'output_IO_upper_NUFIT_contours',
    PATH_FIXED+'output_001_IO_upper_NUFIT_contours', 
    PATH_FIXED+'output_010_IO_upper_NUFIT_contours',
    PATH_FIXED+'output_100_IO_upper_NUFIT_contours'
    ],
    [
    PATH_ALL+'output_NO_upper_JUNO_HK_contours',
    PATH_FIXED+'output_001_NO_upper_JUNO_HK_contours', 
    PATH_FIXED+'output_010_NO_upper_JUNO_HK_contours',
    PATH_FIXED+'output_100_NO_upper_JUNO_HK_contours'
    ],
    [
    PATH_ALL+'output_NO_upper_NUFIT_contours',
    PATH_FIXED+'output_001_NO_upper_NUFIT_contours', 
    PATH_FIXED+'output_010_NO_upper_NUFIT_contours',
    PATH_FIXED+'output_100_NO_upper_NUFIT_contours'
    ]
]

fnames_out = [
    PATH_FIXED+'output_decay_001_010_100_IO_upper_JUNO_HK_contours', 
    PATH_FIXED+'output_decay_001_010_100_IO_upper_NUFIT_contours', 
    PATH_FIXED+'output_decay_001_010_100_NO_upper_JUNO_HK_contours', 
    PATH_FIXED+'output_decay_001_010_100_NO_upper_NUFIT_contours'
]

plot_labels = [
    r'$\nu$ decay, IO, upper $\theta_{23}$ octant,'+'\n'+r'JUNO + HK',
    r'$\nu$ decay, IO, upper $\theta_{23}$ octant,'+'\n'+r'NuFit 5.0',
    r'$\nu$ decay, NO, upper $\theta_{23}$ octant,'+'\n'+r'JUNO + HK',
    r'$\nu$ decay, NO, upper $\theta_{23}$ octant,'+'\n'+r'NuFit 5.0',
]

theta23 = [
    'upper', 
    'upper', 
    'upper', 
    'upper'
]

mass_order = [
   'IO', 
   'IO', 
   'NO', 
   'NO'
]

year = [
    '2040',
    '2020',
    '2040',
    '2020'
]

for i, fname in enumerate(fnames_in):
    print(i)
    plot_ternary_flavor_ratios_fixed_fS_combined_v1(fnames_in[i], fnames_out[i], 
        plot_labels[i], order=mass_order[i], theta23=theta23[i],
        year=year[i], show_benchmarks=True, show_benchmark_astro_labels=False,
        show_benchmark_decay_labels=True,
        output_format='png')
"""


"""
# Fixed flavor composition at the source, combined plots (v1)

# Data for combined experiments
PATH_ALL = os.getcwd()+'/flavor_region_output_JUNO_DUNE_HK_post/all_fS/'
PATH_FIXED = os.getcwd()+'/flavor_region_output_JUNO_DUNE_HK_post/fixed_fS/'
if not os.path.exists(PATH_FIXED):
    os.mkdir(PATH_FIXED)

fnames_in = [
    [
    PATH_ALL+'output_IO_lower_JUNO_DUNE_HK_contours',
    PATH_FIXED+'output_010_IO_lower_JUNO_DUNE_HK_contours', 
    PATH_FIXED+'output_100_IO_lower_JUNO_DUNE_HK_contours',
    PATH_FIXED+'output_120_IO_lower_JUNO_DUNE_HK_contours'
    ],
    [
    PATH_ALL+'output_IO_max_JUNO_DUNE_HK_contours',
    PATH_FIXED+'output_010_IO_max_JUNO_DUNE_HK_contours', 
    PATH_FIXED+'output_100_IO_max_JUNO_DUNE_HK_contours',
    PATH_FIXED+'output_120_IO_max_JUNO_DUNE_HK_contours'
    ],
    [
    PATH_ALL+'output_IO_upper_JUNO_DUNE_HK_contours',
    PATH_FIXED+'output_010_IO_upper_JUNO_DUNE_HK_contours', 
    PATH_FIXED+'output_100_IO_upper_JUNO_DUNE_HK_contours',
    PATH_FIXED+'output_120_IO_upper_JUNO_DUNE_HK_contours'
    ],
    [
    PATH_ALL+'output_IO_upper_NUFIT_contours',
    PATH_FIXED+'output_010_IO_upper_NUFIT_contours', 
    PATH_FIXED+'output_100_IO_upper_NUFIT_contours',
    PATH_FIXED+'output_120_IO_upper_NUFIT_contours'
    ],
    [
    PATH_ALL+'output_NO_lower_JUNO_DUNE_HK_contours',
    PATH_FIXED+'output_010_NO_lower_JUNO_DUNE_HK_contours', 
    PATH_FIXED+'output_100_NO_lower_JUNO_DUNE_HK_contours',
    PATH_FIXED+'output_120_NO_lower_JUNO_DUNE_HK_contours'
    ],
    [
    PATH_ALL+'output_NO_max_JUNO_DUNE_HK_contours',
    PATH_FIXED+'output_010_NO_max_JUNO_DUNE_HK_contours', 
    PATH_FIXED+'output_100_NO_max_JUNO_DUNE_HK_contours',
    PATH_FIXED+'output_120_NO_max_JUNO_DUNE_HK_contours'
    ],
    [
    PATH_ALL+'output_NO_upper_JUNO_DUNE_HK_contours',
    PATH_FIXED+'output_010_NO_upper_JUNO_DUNE_HK_contours', 
    PATH_FIXED+'output_100_NO_upper_JUNO_DUNE_HK_contours',
    PATH_FIXED+'output_120_NO_upper_JUNO_DUNE_HK_contours'
    ],
    [
    PATH_ALL+'output_NO_upper_NUFIT_contours',
    PATH_FIXED+'output_010_NO_upper_NUFIT_contours', 
    PATH_FIXED+'output_100_NO_upper_NUFIT_contours',
    PATH_FIXED+'output_120_NO_upper_NUFIT_contours'
    ]
]

fnames_out = [
    PATH_FIXED+'output_010_100_120_IO_lower_JUNO_DUNE_HK_contours', 
    PATH_FIXED+'output_010_100_120_IO_max_JUNO_DUNE_HK_contours', 
    PATH_FIXED+'output_010_100_120_IO_upper_JUNO_DUNE_HK_contours', 
    PATH_FIXED+'output_010_100_120_IO_upper_NUFIT_contours', 
    PATH_FIXED+'output_010_100_120_NO_lower_JUNO_DUNE_HK_contours', 
    PATH_FIXED+'output_010_100_120_NO_max_JUNO_DUNE_HK_contours', 
    PATH_FIXED+'output_010_100_120_NO_upper_JUNO_DUNE_HK_contours',
    PATH_FIXED+'output_010_100_120_NO_upper_NUFIT_contours'
]

plot_labels = [
    r'IO, lower $\theta_{23}$ octant,'+'\n'+r'JUNO + DUNE + HK',
    r'IO, $\theta_{23} = 45^\circ$,'+'\n'+r'JUNO + DUNE + HK',
    r'IO, upper $\theta_{23}$ octant,'+'\n'+r'JUNO + DUNE + HK',
    r'IO, upper $\theta_{23}$ octant,'+'\n'+r'NuFit 5.0',
    r'NO, lower $\theta_{23}$ octant,'+'\n'+r'JUNO + DUNE + HK',
    r'NO, $\theta_{23} = 45^\circ$,'+'\n'+r'JUNO + DUNE + HK',
    r'NO, upper $\theta_{23}$ octant,'+'\n'+r'JUNO + DUNE + HK',
    r'NO, upper $\theta_{23}$ octant,'+'\n'+r'NuFit 5.0',
]

theta23 = [
    'lower', 
    'max', 
    'upper',
    'upper',
    'lower', 
    'max', 
    'upper',
    'upper'
]

mass_order = [
   'IO', 
   'IO', 
   'IO', 
   'IO', 
   'NO', 
   'NO', 
   'NO', 
   'NO'
]

year = [
    '2040',
    '2040',
    '2040',
    '2020',
    '2040',
    '2040',
    '2040',
    '2020'
]

for i, fname in enumerate(fnames_in):
    print(i)
    plot_ternary_flavor_ratios_fixed_fS_combined_v1(fnames_in[i], fnames_out[i], 
        plot_labels[i], order=mass_order[i], theta23=theta23[i],
        year=year[i], output_format='pdf')
    plot_ternary_flavor_ratios_fixed_fS_combined_v1(fnames_in[i], fnames_out[i], 
        plot_labels[i], order=mass_order[i], theta23=theta23[i],
        year=year[i], output_format='png')
"""


"""
# Fixed flavor composition at the source, combined plots (v1)

# Data for individual experiments
PATH_ALL = os.getcwd()+'/flavor_region_output_post/'
PATH_FIXED = os.getcwd()+'/flavor_region_fixed_source_output_post/'
if not os.path.exists(PATH_FIXED):
    os.mkdir(PATH_FIXED)

fnames_in = [
    [
    PATH_ALL+'output_IO_lower_DUNE_contours',
    PATH_FIXED+'output_010_IO_lower_DUNE_contours', 
    PATH_FIXED+'output_100_IO_lower_DUNE_contours',
    PATH_FIXED+'output_120_IO_lower_DUNE_contours'
    ],
    [
    PATH_ALL+'output_IO_lower_HK_contours',
    PATH_FIXED+'output_010_IO_lower_HK_contours', 
    PATH_FIXED+'output_100_IO_lower_HK_contours',
    PATH_FIXED+'output_120_IO_lower_HK_contours'
    ],
    [
    PATH_ALL+'output_IO_lower_DUNE_JUNO_contours',
    PATH_FIXED+'output_010_IO_lower_JUNO_DUNE_contours', 
    PATH_FIXED+'output_100_IO_lower_JUNO_DUNE_contours',
    PATH_FIXED+'output_120_IO_lower_JUNO_DUNE_contours'
    ],
    [
    PATH_ALL+'output_IO_lower_HK_JUNO_contours',
    PATH_FIXED+'output_010_IO_lower_JUNO_HK_contours', 
    PATH_FIXED+'output_100_IO_lower_JUNO_HK_contours',
    PATH_FIXED+'output_120_IO_lower_JUNO_HK_contours'
    ],
    [
    PATH_ALL+'output_IO_max_DUNE_contours',
    PATH_FIXED+'output_010_IO_max_DUNE_contours', 
    PATH_FIXED+'output_100_IO_max_DUNE_contours',
    PATH_FIXED+'output_120_IO_max_DUNE_contours'
    ],
    [
    PATH_ALL+'output_IO_max_HK_contours',
    PATH_FIXED+'output_010_IO_max_HK_contours', 
    PATH_FIXED+'output_100_IO_max_HK_contours',
    PATH_FIXED+'output_120_IO_max_HK_contours'
    ],
    [
    PATH_ALL+'output_IO_max_DUNE_JUNO_contours',
    PATH_FIXED+'output_010_IO_max_JUNO_DUNE_contours', 
    PATH_FIXED+'output_100_IO_max_JUNO_DUNE_contours',
    PATH_FIXED+'output_120_IO_max_JUNO_DUNE_contours'
    ],
    [
    PATH_ALL+'output_IO_max_HK_JUNO_contours',
    PATH_FIXED+'output_010_IO_max_JUNO_HK_contours', 
    PATH_FIXED+'output_100_IO_max_JUNO_HK_contours',
    PATH_FIXED+'output_120_IO_max_JUNO_HK_contours'
    ],
    [
    PATH_ALL+'output_IO_upper_DUNE_contours',
    PATH_FIXED+'output_010_IO_upper_DUNE_contours', 
    PATH_FIXED+'output_100_IO_upper_DUNE_contours',
    PATH_FIXED+'output_120_IO_upper_DUNE_contours'
    ],
    [
    PATH_ALL+'output_IO_upper_HK_contours',
    PATH_FIXED+'output_010_IO_upper_HK_contours', 
    PATH_FIXED+'output_100_IO_upper_HK_contours',
    PATH_FIXED+'output_120_IO_upper_HK_contours'
    ],
    [
    PATH_ALL+'output_IO_upper_DUNE_JUNO_contours',
    PATH_FIXED+'output_010_IO_upper_JUNO_DUNE_contours', 
    PATH_FIXED+'output_100_IO_upper_JUNO_DUNE_contours',
    PATH_FIXED+'output_120_IO_upper_JUNO_DUNE_contours'
    ],
    [
    PATH_ALL+'output_IO_upper_HK_JUNO_contours',
    PATH_FIXED+'output_010_IO_upper_JUNO_HK_contours', 
    PATH_FIXED+'output_100_IO_upper_JUNO_HK_contours',
    PATH_FIXED+'output_120_IO_upper_JUNO_HK_contours'
    ],
    [
    PATH_ALL+'output_IO_upper_NUFIT_JUNO_contours',
    PATH_FIXED+'output_010_IO_upper_JUNO_NUFIT_contours', 
    PATH_FIXED+'output_100_IO_upper_JUNO_NUFIT_contours',
    PATH_FIXED+'output_120_IO_upper_JUNO_NUFIT_contours'
    ],
    [
    PATH_ALL+'output_IO_upper_NUFIT_contours',
    PATH_FIXED+'output_010_IO_upper_NUFIT_contours', 
    PATH_FIXED+'output_100_IO_upper_NUFIT_contours',
    PATH_FIXED+'output_120_IO_upper_NUFIT_contours'
    ],
    [
    PATH_ALL+'output_NO_lower_DUNE_contours',
    PATH_FIXED+'output_010_NO_lower_DUNE_contours', 
    PATH_FIXED+'output_100_NO_lower_DUNE_contours',
    PATH_FIXED+'output_120_NO_lower_DUNE_contours'
    ],
    [
    PATH_ALL+'output_NO_lower_HK_contours',
    PATH_FIXED+'output_010_NO_lower_HK_contours', 
    PATH_FIXED+'output_100_NO_lower_HK_contours',
    PATH_FIXED+'output_120_NO_lower_HK_contours'
    ],
    [
    PATH_ALL+'output_NO_lower_DUNE_JUNO_contours', 
    PATH_FIXED+'output_010_NO_lower_JUNO_DUNE_contours', 
    PATH_FIXED+'output_100_NO_lower_JUNO_DUNE_contours',
    PATH_FIXED+'output_120_NO_lower_JUNO_DUNE_contours'
    ],
    [
    PATH_ALL+'output_NO_lower_HK_JUNO_contours', 
    PATH_FIXED+'output_010_NO_lower_JUNO_HK_contours', 
    PATH_FIXED+'output_100_NO_lower_JUNO_HK_contours',
    PATH_FIXED+'output_120_NO_lower_JUNO_HK_contours'
    ],
    [
    PATH_ALL+'output_NO_max_DUNE_contours', 
    PATH_FIXED+'output_010_NO_max_DUNE_contours', 
    PATH_FIXED+'output_100_NO_max_DUNE_contours',
    PATH_FIXED+'output_120_NO_max_DUNE_contours'
    ],
    [
    PATH_ALL+'output_NO_max_HK_contours', 
    PATH_FIXED+'output_010_NO_max_HK_contours', 
    PATH_FIXED+'output_100_NO_max_HK_contours',
    PATH_FIXED+'output_120_NO_max_HK_contours'
    ],
    [
    PATH_ALL+'output_NO_max_DUNE_JUNO_contours', 
    PATH_FIXED+'output_010_NO_max_JUNO_DUNE_contours', 
    PATH_FIXED+'output_100_NO_max_JUNO_DUNE_contours',
    PATH_FIXED+'output_120_NO_max_JUNO_DUNE_contours'
    ],
    [
    PATH_ALL+'output_NO_max_HK_JUNO_contours', 
    PATH_FIXED+'output_010_NO_max_JUNO_HK_contours', 
    PATH_FIXED+'output_100_NO_max_JUNO_HK_contours',
    PATH_FIXED+'output_120_NO_max_JUNO_HK_contours'
    ],
    [
    PATH_ALL+'output_NO_upper_DUNE_contours', 
    PATH_FIXED+'output_010_NO_upper_DUNE_contours', 
    PATH_FIXED+'output_100_NO_upper_DUNE_contours',
    PATH_FIXED+'output_120_NO_upper_DUNE_contours'
    ],
    [
    PATH_ALL+'output_NO_upper_HK_contours', 
    PATH_FIXED+'output_010_NO_upper_HK_contours', 
    PATH_FIXED+'output_100_NO_upper_HK_contours',
    PATH_FIXED+'output_120_NO_upper_HK_contours'
    ],
    [
    PATH_ALL+'output_NO_upper_DUNE_JUNO_contours', 
    PATH_FIXED+'output_010_NO_upper_JUNO_DUNE_contours', 
    PATH_FIXED+'output_100_NO_upper_JUNO_DUNE_contours',
    PATH_FIXED+'output_120_NO_upper_JUNO_DUNE_contours'
    ],
    [
    PATH_ALL+'output_NO_upper_HK_JUNO_contours', 
    PATH_FIXED+'output_010_NO_upper_JUNO_HK_contours', 
    PATH_FIXED+'output_100_NO_upper_JUNO_HK_contours',
    PATH_FIXED+'output_120_NO_upper_JUNO_HK_contours'
    ],
    [
    PATH_ALL+'output_NO_upper_NUFIT_JUNO_contours', 
    PATH_FIXED+'output_010_NO_upper_JUNO_NUFIT_contours', 
    PATH_FIXED+'output_100_NO_upper_JUNO_NUFIT_contours',
    PATH_FIXED+'output_120_NO_upper_JUNO_NUFIT_contours'
    ],
    [
    PATH_ALL+'output_NO_upper_NUFIT_contours', 
    PATH_FIXED+'output_010_NO_upper_NUFIT_contours', 
    PATH_FIXED+'output_100_NO_upper_NUFIT_contours',
    PATH_FIXED+'output_120_NO_upper_NUFIT_contours'
    ]
]

fnames_out = [
    PATH_FIXED+'output_010_100_120_IO_lower_DUNE_contours', 
    PATH_FIXED+'output_010_100_120_IO_lower_HK_contours', 
    PATH_FIXED+'output_010_100_120_IO_lower_JUNO_DUNE_contours', 
    PATH_FIXED+'output_010_100_120_IO_lower_JUNO_HK_contours', 
    PATH_FIXED+'output_010_100_120_IO_max_DUNE_contours', 
    PATH_FIXED+'output_010_100_120_IO_max_HK_contours', 
    PATH_FIXED+'output_010_100_120_IO_max_JUNO_DUNE_contours', 
    PATH_FIXED+'output_010_100_120_IO_max_JUNO_HK_contours', 
    PATH_FIXED+'output_010_100_120_IO_upper_DUNE_contours', 
    PATH_FIXED+'output_010_100_120_IO_upper_HK_contours', 
    PATH_FIXED+'output_010_100_120_IO_upper_JUNO_DUNE_contours', 
    PATH_FIXED+'output_010_100_120_IO_upper_JUNO_HK_contours', 
    PATH_FIXED+'output_010_100_120_IO_upper_JUNO_NUFIT_contours', 
    PATH_FIXED+'output_010_100_120_IO_upper_NUFIT_contours', 
    PATH_FIXED+'output_010_100_120_NO_lower_DUNE_contours', 
    PATH_FIXED+'output_010_100_120_NO_lower_HK_contours', 
    PATH_FIXED+'output_010_100_120_NO_lower_JUNO_DUNE_contours', 
    PATH_FIXED+'output_010_100_120_NO_lower_JUNO_HK_contours', 
    PATH_FIXED+'output_010_100_120_NO_max_DUNE_contours', 
    PATH_FIXED+'output_010_100_120_NO_max_HK_contours', 
    PATH_FIXED+'output_010_100_120_NO_max_JUNO_DUNE_contours', 
    PATH_FIXED+'output_010_100_120_NO_max_JUNO_HK_contours', 
    PATH_FIXED+'output_010_100_120_NO_upper_DUNE_contours', 
    PATH_FIXED+'output_010_100_120_NO_upper_HK_contours', 
    PATH_FIXED+'output_010_100_120_NO_upper_JUNO_DUNE_contours', 
    PATH_FIXED+'output_010_100_120_NO_upper_JUNO_HK_contours', 
    PATH_FIXED+'output_010_100_120_NO_upper_JUNO_NUFIT_contours', 
    PATH_FIXED+'output_010_100_120_NO_upper_NUFIT_contours'
]

plot_labels = [
    r'IO, lower $\theta_{23}$ octant,'+'\n'+r'DUNE',
    r'IO, lower $\theta_{23}$ octant,'+'\n'+r'HK',
    r'IO, lower $\theta_{23}$ octant,'+'\n'+r'JUNO + DUNE',
    r'IO, lower $\theta_{23}$ octant,'+'\n'+r'JUNO + HK',
    r'IO, $\theta_{23} = 45^\circ$,'+'\n'+r'DUNE',
    r'IO, $\theta_{23} = 45^\circ$,'+'\n'+r'HK',
    r'IO, $\theta_{23} = 45^\circ$,'+'\n'+r'JUNO + DUNE',
    r'IO, $\theta_{23} = 45^\circ$,'+'\n'+r'JUNO + HK',
    r'IO, upper $\theta_{23}$ octant,'+'\n'+r'DUNE',
    r'IO, upper $\theta_{23}$ octant,'+'\n'+r'HK',
    r'IO, upper $\theta_{23}$ octant,'+'\n'+r'JUNO + DUNE',
    r'IO, upper $\theta_{23}$ octant,'+'\n'+r'JUNO + HK',
    r'IO, upper $\theta_{23}$ octant,'+'\n'+r'NuFit 5.0 + JUNO',
    r'IO, upper $\theta_{23}$ octant,'+'\n'+r'NuFit 5.0',
    r'NO, lower $\theta_{23}$ octant,'+'\n'+r'DUNE',
    r'NO, lower $\theta_{23}$ octant,'+'\n'+r'HK',
    r'NO, lower $\theta_{23}$ octant,'+'\n'+r'JUNO + DUNE',
    r'NO, lower $\theta_{23}$ octant,'+'\n'+r'JUNO + HK',
    r'NO, $\theta_{23} = 45^\circ$,'+'\n'+r'DUNE',
    r'NO, $\theta_{23} = 45^\circ$,'+'\n'+r'HK',
    r'NO, $\theta_{23} = 45^\circ$,'+'\n'+r'JUNO + DUNE',
    r'NO, $\theta_{23} = 45^\circ$,'+'\n'+r'JUNO + HK',
    r'NO, upper $\theta_{23}$ octant,'+'\n'+r'DUNE',
    r'NO, upper $\theta_{23}$ octant,'+'\n'+r'HK',
    r'NO, upper $\theta_{23}$ octant,'+'\n'+r'JUNO + DUNE',
    r'NO, upper $\theta_{23}$ octant,'+'\n'+r'JUNO + HK',
    r'NO, upper $\theta_{23}$ octant,'+'\n'+r'NuFit 5.0 + JUNO',
    r'NO, upper $\theta_{23}$ octant,'+'\n'+r'NuFit 5.0']

theta23 = [
    'lower', 
    'lower', 
    'lower', 
    'lower', 
    'max', 
    'max', 
    'max', 
    'max', 
    'upper', 
    'upper', 
    'upper', 
    'upper', 
    'upper', 
    'upper',
    'lower', 
    'lower', 
    'lower', 
    'lower', 
    'max', 
    'max', 
    'max', 
    'max', 
    'upper', 
    'upper', 
    'upper', 
    'upper', 
    'upper', 
    'upper']

mass_order = [
   'IO', 
   'IO', 
   'IO', 
   'IO', 
   'IO', 
   'IO', 
   'IO', 
   'IO', 
   'IO', 
   'IO', 
   'IO', 
   'IO', 
   'IO', 
   'IO', 
   'NO', 
   'NO', 
   'NO', 
   'NO', 
   'NO', 
   'NO', 
   'NO', 
   'NO', 
   'NO', 
   'NO', 
   'NO', 
   'NO', 
   'NO', 
   'NO']

year = [
    '2030',
    '2030',
    '2040',
    '2030',
    '2040',
    '2030',
    '2040',
    '2030',
    '2040',
    '2030',
    '2040',
    '2030',
    '2030',
    '2020',
    '2040',
    '2030',
    '2040',
    '2030',
    '2040',
    '2030',
    '2040',
    '2030',
    '2040',
    '2030',
    '2040',
    '2030',
    '2030',
    '2020']

for i, fname in enumerate(fnames_in):
    print(i)
    plot_ternary_flavor_ratios_fixed_fS_combined_v1(fnames_in[i], fnames_out[i], 
        plot_labels[i], order=mass_order[i], theta23=theta23[i],
        year=year[i], output_format='png')
    plot_ternary_flavor_ratios_fixed_fS_combined_v1(fnames_in[i], fnames_out[i], 
        plot_labels[i], order=mass_order[i], theta23=theta23[i],
        year=year[i], output_format='pdf')
"""


"""
# Fixed flavor composition at the source, combined plots (v0)

PATH_IN = os.getcwd()+'/flavor_region_fixed_source_output/'
PATH_OUT = os.getcwd()+'/flavor_region_fixed_source_output_post/'
if not os.path.exists(PATH_OUT):
    os.mkdir(PATH_OUT)
# fnames = sorted(filenames_no_ext(PATH_IN)) # Keep the same file names throughout
# fnames_grouped = [fnames[n:n+3] for n in range(0, len(fnames), 3)]
# print(np.array(fnames))
# print(fnames_grouped)
# quit()

fnames_in = [
    [
    PATH_OUT+'output_010_IO_lower_DUNE_contours', 
    PATH_OUT+'output_100_IO_lower_DUNE_contours',
    PATH_OUT+'output_120_IO_lower_DUNE_contours'
    ],
    [
    PATH_OUT+'output_010_IO_lower_HK_contours', 
    PATH_OUT+'output_100_IO_lower_HK_contours',
    PATH_OUT+'output_120_IO_lower_HK_contours'
    ],
    [
    PATH_OUT+'output_010_IO_lower_JUNO_DUNE_contours', 
    PATH_OUT+'output_100_IO_lower_JUNO_DUNE_contours',
    PATH_OUT+'output_120_IO_lower_JUNO_DUNE_contours'
    ],
    [
    PATH_OUT+'output_010_IO_lower_JUNO_HK_contours', 
    PATH_OUT+'output_100_IO_lower_JUNO_HK_contours',
    PATH_OUT+'output_120_IO_lower_JUNO_HK_contours'
    ],
    [
    PATH_OUT+'output_010_IO_max_DUNE_contours', 
    PATH_OUT+'output_100_IO_max_DUNE_contours',
    PATH_OUT+'output_120_IO_max_DUNE_contours'
    ],
    [
    PATH_OUT+'output_010_IO_max_HK_contours', 
    PATH_OUT+'output_100_IO_max_HK_contours',
    PATH_OUT+'output_120_IO_max_HK_contours'
    ],
    [
    PATH_OUT+'output_010_IO_max_JUNO_DUNE_contours', 
    PATH_OUT+'output_100_IO_max_JUNO_DUNE_contours',
    PATH_OUT+'output_120_IO_max_JUNO_DUNE_contours'
    ],
    [
    PATH_OUT+'output_010_IO_max_JUNO_HK_contours', 
    PATH_OUT+'output_100_IO_max_JUNO_HK_contours',
    PATH_OUT+'output_120_IO_max_JUNO_HK_contours'
    ],
    [
    PATH_OUT+'output_010_IO_upper_DUNE_contours', 
    PATH_OUT+'output_100_IO_upper_DUNE_contours',
    PATH_OUT+'output_120_IO_upper_DUNE_contours'
    ],
    [
    PATH_OUT+'output_010_IO_upper_HK_contours', 
    PATH_OUT+'output_100_IO_upper_HK_contours',
    PATH_OUT+'output_120_IO_upper_HK_contours'
    ],
    [
    PATH_OUT+'output_010_IO_upper_JUNO_DUNE_contours', 
    PATH_OUT+'output_100_IO_upper_JUNO_DUNE_contours',
    PATH_OUT+'output_120_IO_upper_JUNO_DUNE_contours'
    ],
    [
    PATH_OUT+'output_010_IO_upper_JUNO_HK_contours', 
    PATH_OUT+'output_100_IO_upper_JUNO_HK_contours',
    PATH_OUT+'output_120_IO_upper_JUNO_HK_contours'
    ],
    [
    PATH_OUT+'output_010_IO_upper_JUNO_NUFIT_contours', 
    PATH_OUT+'output_100_IO_upper_JUNO_NUFIT_contours',
    PATH_OUT+'output_120_IO_upper_JUNO_NUFIT_contours'
    ],
    [
    PATH_OUT+'output_010_IO_upper_NUFIT_contours', 
    PATH_OUT+'output_100_IO_upper_NUFIT_contours',
    PATH_OUT+'output_120_IO_upper_NUFIT_contours'
    ],
    [
    PATH_OUT+'output_010_NO_lower_DUNE_contours', 
    PATH_OUT+'output_100_NO_lower_DUNE_contours',
    PATH_OUT+'output_120_NO_lower_DUNE_contours'
    ],
    [
    PATH_OUT+'output_010_NO_lower_HK_contours', 
    PATH_OUT+'output_100_NO_lower_HK_contours',
    PATH_OUT+'output_120_NO_lower_HK_contours'
    ],
    [
    PATH_OUT+'output_010_NO_lower_JUNO_DUNE_contours', 
    PATH_OUT+'output_100_NO_lower_JUNO_DUNE_contours',
    PATH_OUT+'output_120_NO_lower_JUNO_DUNE_contours'
    ],
    [
    PATH_OUT+'output_010_NO_lower_JUNO_HK_contours', 
    PATH_OUT+'output_100_NO_lower_JUNO_HK_contours',
    PATH_OUT+'output_120_NO_lower_JUNO_HK_contours'
    ],
    [
    PATH_OUT+'output_010_NO_max_DUNE_contours', 
    PATH_OUT+'output_100_NO_max_DUNE_contours',
    PATH_OUT+'output_120_NO_max_DUNE_contours'
    ],
    [
    PATH_OUT+'output_010_NO_max_HK_contours', 
    PATH_OUT+'output_100_NO_max_HK_contours',
    PATH_OUT+'output_120_NO_max_HK_contours'
    ],
    [
    PATH_OUT+'output_010_NO_max_JUNO_DUNE_contours', 
    PATH_OUT+'output_100_NO_max_JUNO_DUNE_contours',
    PATH_OUT+'output_120_NO_max_JUNO_DUNE_contours'
    ],
    [
    PATH_OUT+'output_010_NO_max_JUNO_HK_contours', 
    PATH_OUT+'output_100_NO_max_JUNO_HK_contours',
    PATH_OUT+'output_120_NO_max_JUNO_HK_contours'
    ],
    [
    PATH_OUT+'output_010_NO_upper_DUNE_contours', 
    PATH_OUT+'output_100_NO_upper_DUNE_contours',
    PATH_OUT+'output_120_NO_upper_DUNE_contours'
    ],
    [
    PATH_OUT+'output_010_NO_upper_HK_contours', 
    PATH_OUT+'output_100_NO_upper_HK_contours',
    PATH_OUT+'output_120_NO_upper_HK_contours'
    ],
    [
    PATH_OUT+'output_010_NO_upper_JUNO_DUNE_contours', 
    PATH_OUT+'output_100_NO_upper_JUNO_DUNE_contours',
    PATH_OUT+'output_120_NO_upper_JUNO_DUNE_contours'
    ],
    [
    PATH_OUT+'output_010_NO_upper_JUNO_HK_contours', 
    PATH_OUT+'output_100_NO_upper_JUNO_HK_contours',
    PATH_OUT+'output_120_NO_upper_JUNO_HK_contours'
    ],
    [
    PATH_OUT+'output_010_NO_upper_JUNO_NUFIT_contours', 
    PATH_OUT+'output_100_NO_upper_JUNO_NUFIT_contours',
    PATH_OUT+'output_120_NO_upper_JUNO_NUFIT_contours'
    ],
    [
    PATH_OUT+'output_010_NO_upper_NUFIT_contours', 
    PATH_OUT+'output_100_NO_upper_NUFIT_contours',
    PATH_OUT+'output_120_NO_upper_NUFIT_contours'
    ]
]

fnames_out = [
    PATH_OUT+'output_010_100_120_IO_lower_DUNE_contours', 
    PATH_OUT+'output_010_100_120_IO_lower_HK_contours', 
    PATH_OUT+'output_010_100_120_IO_lower_JUNO_DUNE_contours', 
    PATH_OUT+'output_010_100_120_IO_lower_JUNO_HK_contours', 
    PATH_OUT+'output_010_100_120_IO_max_DUNE_contours', 
    PATH_OUT+'output_010_100_120_IO_max_HK_contours', 
    PATH_OUT+'output_010_100_120_IO_max_JUNO_DUNE_contours', 
    PATH_OUT+'output_010_100_120_IO_max_JUNO_HK_contours', 
    PATH_OUT+'output_010_100_120_IO_upper_DUNE_contours', 
    PATH_OUT+'output_010_100_120_IO_upper_HK_contours', 
    PATH_OUT+'output_010_100_120_IO_upper_JUNO_DUNE_contours', 
    PATH_OUT+'output_010_100_120_IO_upper_JUNO_HK_contours', 
    PATH_OUT+'output_010_100_120_IO_upper_JUNO_NUFIT_contours', 
    PATH_OUT+'output_010_100_120_IO_upper_NUFIT_contours', 
    PATH_OUT+'output_010_100_120_NO_lower_DUNE_contours', 
    PATH_OUT+'output_010_100_120_NO_lower_HK_contours', 
    PATH_OUT+'output_010_100_120_NO_lower_JUNO_DUNE_contours', 
    PATH_OUT+'output_010_100_120_NO_lower_JUNO_HK_contours', 
    PATH_OUT+'output_010_100_120_NO_max_DUNE_contours', 
    PATH_OUT+'output_010_100_120_NO_max_HK_contours', 
    PATH_OUT+'output_010_100_120_NO_max_JUNO_DUNE_contours', 
    PATH_OUT+'output_010_100_120_NO_max_JUNO_HK_contours', 
    PATH_OUT+'output_010_100_120_NO_upper_DUNE_contours', 
    PATH_OUT+'output_010_100_120_NO_upper_HK_contours', 
    PATH_OUT+'output_010_100_120_NO_upper_JUNO_DUNE_contours', 
    PATH_OUT+'output_010_100_120_NO_upper_JUNO_HK_contours', 
    PATH_OUT+'output_010_100_120_NO_upper_JUNO_NUFIT_contours', 
    PATH_OUT+'output_010_100_120_NO_upper_NUFIT_contours'
]

plot_labels = [
    r'IO, lower $\theta_{23}$ octant,'+'\n'+r'DUNE',
    r'IO, lower $\theta_{23}$ octant,'+'\n'+r'HK',
    r'IO, lower $\theta_{23}$ octant,'+'\n'+r'JUNO + DUNE',
    r'IO, lower $\theta_{23}$ octant,'+'\n'+r'JUNO + HK',
    r'IO, $\theta_{23} = 45^\circ$,'+'\n'+r'DUNE',
    r'IO, $\theta_{23} = 45^\circ$,'+'\n'+r'HK',
    r'IO, $\theta_{23} = 45^\circ$,'+'\n'+r'JUNO + DUNE',
    r'IO, $\theta_{23} = 45^\circ$,'+'\n'+r'JUNO + HK',
    r'IO, upper $\theta_{23}$ octant,'+'\n'+r'DUNE',
    r'IO, upper $\theta_{23}$ octant,'+'\n'+r'HK',
    r'IO, upper $\theta_{23}$ octant,'+'\n'+r'JUNO + DUNE',
    r'IO, upper $\theta_{23}$ octant,'+'\n'+r'JUNO + HK',
    r'IO, upper $\theta_{23}$ octant,'+'\n'+r'NuFit 5.0 + JUNO',
    r'IO, upper $\theta_{23}$ octant,'+'\n'+r'NuFit 5.0',
    r'NO, lower $\theta_{23}$ octant,'+'\n'+r'DUNE',
    r'NO, lower $\theta_{23}$ octant,'+'\n'+r'HK',
    r'NO, lower $\theta_{23}$ octant,'+'\n'+r'JUNO + DUNE',
    r'NO, lower $\theta_{23}$ octant,'+'\n'+r'JUNO + HK',
    r'NO, $\theta_{23} = 45^\circ$,'+'\n'+r'DUNE',
    r'NO, $\theta_{23} = 45^\circ$,'+'\n'+r'HK',
    r'NO, $\theta_{23} = 45^\circ$,'+'\n'+r'JUNO + DUNE',
    r'NO, $\theta_{23} = 45^\circ$,'+'\n'+r'JUNO + HK',
    r'NO, upper $\theta_{23}$ octant,'+'\n'+r'DUNE',
    r'NO, upper $\theta_{23}$ octant,'+'\n'+r'HK',
    r'NO, upper $\theta_{23}$ octant,'+'\n'+r'JUNO + DUNE',
    r'NO, upper $\theta_{23}$ octant,'+'\n'+r'JUNO + HK',
    r'NO, upper $\theta_{23}$ octant,'+'\n'+r'NuFit 5.0 + JUNO',
    r'NO, upper $\theta_{23}$ octant,'+'\n'+r'NuFit 5.0']

theta23 = [
    'lower', 'lower', 'lower', 'lower', 'max', 'max', 'max', 'max', 
    'upper', 'upper', 'upper', 'upper', 'upper', 'upper',
    'lower', 'lower', 'lower', 'lower', 'max', 'max', 'max', 'max', 
    'upper', 'upper', 'upper', 'upper', 'upper', 'upper']

mass_order = [
   'IO', 'IO', 'IO', 'IO', 'IO', 'IO', 'IO', 'IO', 'IO', 'IO', 'IO', 'IO', 'IO', 
   'IO', 
   'NO', 'NO', 'NO', 'NO', 'NO', 'NO', 'NO', 'NO', 'NO', 'NO', 'NO', 'NO', 'NO', 
   'NO']

for i, fname in enumerate(fnames_in):
    print(i)
    plot_ternary_flavor_ratios_fixed_fS_combined_v0(fnames_in[i], fnames_out[i], 
        plot_labels[i], order=mass_order[i], theta23=theta23[i],
        output_format='png')
"""



"""
# Fixed flavor composition at the source, separate plots

PATH_IN = os.getcwd()+'/flavor_region_fixed_source_output/'
PATH_OUT = os.getcwd()+'/flavor_region_fixed_source_output_post/'
if not os.path.exists(PATH_OUT):
    os.mkdir(PATH_OUT)
fnames = sorted(filenames_no_ext(PATH_IN)) # Keep the same file names throughout
# print(np.array(fnames))
# quit()

plot_labels = [
    r'$(0:1:0)_{\rm S}$, IO, lower $\theta_{23}$ octant,'+'\n'+r'DUNE',
    r'$(0:1:0)_{\rm S}$, IO, lower $\theta_{23}$ octant,'+'\n'+r'HK',
    r'$(0:1:0)_{\rm S}$, IO, lower $\theta_{23}$ octant,'+'\n'+r'JUNO + DUNE',
    r'$(0:1:0)_{\rm S}$, IO, lower $\theta_{23}$ octant,'+'\n'+r'JUNO + HK',
    r'$(0:1:0)_{\rm S}$, IO, $\theta_{23} = 45^\circ$,'+'\n'+r'DUNE',
    r'$(0:1:0)_{\rm S}$, IO, $\theta_{23} = 45^\circ$,'+'\n'+r'HK',
    r'$(0:1:0)_{\rm S}$, IO, $\theta_{23} = 45^\circ$,'+'\n'+r'JUNO + DUNE',
    r'$(0:1:0)_{\rm S}$, IO, $\theta_{23} = 45^\circ$,'+'\n'+r'JUNO + HK',
    r'$(0:1:0)_{\rm S}$, IO, upper $\theta_{23}$ octant,'+'\n'+r'DUNE',
    r'$(0:1:0)_{\rm S}$, IO, upper $\theta_{23}$ octant,'+'\n'+r'HK',
    r'$(0:1:0)_{\rm S}$, IO, upper $\theta_{23}$ octant,'+'\n'+r'JUNO + DUNE',
    r'$(0:1:0)_{\rm S}$, IO, upper $\theta_{23}$ octant,'+'\n'+r'JUNO + HK',
    r'$(0:1:0)_{\rm S}$, IO, upper $\theta_{23}$ octant,'+'\n'+r'NuFit 5.0 + JUNO',
    r'$(0:1:0)_{\rm S}$, IO, upper $\theta_{23}$ octant,'+'\n'+r'NuFit 5.0',
    r'$(0:1:0)_{\rm S}$, NO, lower $\theta_{23}$ octant,'+'\n'+r'DUNE',
    r'$(0:1:0)_{\rm S}$, NO, lower $\theta_{23}$ octant,'+'\n'+r'HK',
    r'$(0:1:0)_{\rm S}$, NO, lower $\theta_{23}$ octant,'+'\n'+r'JUNO + DUNE',
    r'$(0:1:0)_{\rm S}$, NO, lower $\theta_{23}$ octant,'+'\n'+r'JUNO + HK',
    r'$(0:1:0)_{\rm S}$, NO, $\theta_{23} = 45^\circ$,'+'\n'+r'DUNE',
    r'$(0:1:0)_{\rm S}$, NO, $\theta_{23} = 45^\circ$,'+'\n'+r'HK',
    r'$(0:1:0)_{\rm S}$, NO, $\theta_{23} = 45^\circ$,'+'\n'+r'JUNO + DUNE',
    r'$(0:1:0)_{\rm S}$, NO, $\theta_{23} = 45^\circ$,'+'\n'+r'JUNO + HK',
    r'$(0:1:0)_{\rm S}$, NO, upper $\theta_{23}$ octant,'+'\n'+r'DUNE',
    r'$(0:1:0)_{\rm S}$, NO, upper $\theta_{23}$ octant,'+'\n'+r'HK',
    r'$(0:1:0)_{\rm S}$, NO, upper $\theta_{23}$ octant,'+'\n'+r'JUNO + DUNE',
    r'$(0:1:0)_{\rm S}$, NO, upper $\theta_{23}$ octant,'+'\n'+r'JUNO + HK',
    r'$(0:1:0)_{\rm S}$, NO, upper $\theta_{23}$ octant,'+'\n'+r'NuFit 5.0 + JUNO',
    r'$(0:1:0)_{\rm S}$, NO, upper $\theta_{23}$ octant,'+'\n'+r'NuFit 5.0',
    r'$(1:0:0)_{\rm S}$, IO, lower $\theta_{23}$ octant,'+'\n'+r'DUNE',
    r'$(1:0:0)_{\rm S}$, IO, lower $\theta_{23}$ octant,'+'\n'+r'HK',
    r'$(1:0:0)_{\rm S}$, IO, lower $\theta_{23}$ octant,'+'\n'+r'JUNO + DUNE',
    r'$(1:0:0)_{\rm S}$, IO, lower $\theta_{23}$ octant,'+'\n'+r'JUNO + HK',
    r'$(1:0:0)_{\rm S}$, IO, $\theta_{23} = 45^\circ$,'+'\n'+r'DUNE',
    r'$(1:0:0)_{\rm S}$, IO, $\theta_{23} = 45^\circ$,'+'\n'+r'HK',
    r'$(1:0:0)_{\rm S}$, IO, $\theta_{23} = 45^\circ$,'+'\n'+r'JUNO + DUNE',
    r'$(1:0:0)_{\rm S}$, IO, $\theta_{23} = 45^\circ$,'+'\n'+r'JUNO + HK',
    r'$(1:0:0)_{\rm S}$, IO, upper $\theta_{23}$ octant,'+'\n'+r'DUNE',
    r'$(1:0:0)_{\rm S}$, IO, upper $\theta_{23}$ octant,'+'\n'+r'HK',
    r'$(1:0:0)_{\rm S}$, IO, upper $\theta_{23}$ octant,'+'\n'+r'JUNO + DUNE',
    r'$(1:0:0)_{\rm S}$, IO, upper $\theta_{23}$ octant,'+'\n'+r'JUNO + HK',
    r'$(1:0:0)_{\rm S}$, IO, upper $\theta_{23}$ octant,'+'\n'+r'NuFit 5.0 + JUNO',
    r'$(1:0:0)_{\rm S}$, IO, upper $\theta_{23}$ octant,'+'\n'+r'NuFit 5.0',
    r'$(1:0:0)_{\rm S}$, NO, lower $\theta_{23}$ octant,'+'\n'+r'DUNE',
    r'$(1:0:0)_{\rm S}$, NO, lower $\theta_{23}$ octant,'+'\n'+r'HK',
    r'$(1:0:0)_{\rm S}$, NO, lower $\theta_{23}$ octant,'+'\n'+r'JUNO + DUNE',
    r'$(1:0:0)_{\rm S}$, NO, lower $\theta_{23}$ octant,'+'\n'+r'JUNO + HK',
    r'$(1:0:0)_{\rm S}$, NO, $\theta_{23} = 45^\circ$,'+'\n'+r'DUNE',
    r'$(1:0:0)_{\rm S}$, NO, $\theta_{23} = 45^\circ$,'+'\n'+r'HK',
    r'$(1:0:0)_{\rm S}$, NO, $\theta_{23} = 45^\circ$,'+'\n'+r'JUNO + DUNE',
    r'$(1:0:0)_{\rm S}$, NO, $\theta_{23} = 45^\circ$,'+'\n'+r'JUNO + HK',
    r'$(1:0:0)_{\rm S}$, NO, upper $\theta_{23}$ octant,'+'\n'+r'DUNE',
    r'$(1:0:0)_{\rm S}$, NO, upper $\theta_{23}$ octant,'+'\n'+r'HK',
    r'$(1:0:0)_{\rm S}$, NO, upper $\theta_{23}$ octant,'+'\n'+r'JUNO + DUNE',
    r'$(1:0:0)_{\rm S}$, NO, upper $\theta_{23}$ octant,'+'\n'+r'JUNO + HK',
    r'$(1:0:0)_{\rm S}$, NO, upper $\theta_{23}$ octant,'+'\n'+r'NuFit 5.0 + JUNO',
    r'$(1:0:0)_{\rm S}$, NO, upper $\theta_{23}$ octant,'+'\n'+r'NuFit 5.0',
    r'$(1:2:0)_{\rm S}$, IO, lower $\theta_{23}$ octant,'+'\n'+r'DUNE',
    r'$(1:2:0)_{\rm S}$, IO, lower $\theta_{23}$ octant,'+'\n'+r'HK',
    r'$(1:2:0)_{\rm S}$, IO, lower $\theta_{23}$ octant,'+'\n'+r'JUNO + DUNE',
    r'$(1:2:0)_{\rm S}$, IO, lower $\theta_{23}$ octant,'+'\n'+r'JUNO + HK',
    r'$(1:2:0)_{\rm S}$, IO, $\theta_{23} = 45^\circ$,'+'\n'+r'DUNE',
    r'$(1:2:0)_{\rm S}$, IO, $\theta_{23} = 45^\circ$,'+'\n'+r'HK',
    r'$(1:2:0)_{\rm S}$, IO, $\theta_{23} = 45^\circ$,'+'\n'+r'JUNO + DUNE',
    r'$(1:2:0)_{\rm S}$, IO, $\theta_{23} = 45^\circ$,'+'\n'+r'JUNO + HK',
    r'$(1:2:0)_{\rm S}$, IO, upper $\theta_{23}$ octant,'+'\n'+r'DUNE',
    r'$(1:2:0)_{\rm S}$, IO, upper $\theta_{23}$ octant,'+'\n'+r'HK',
    r'$(1:2:0)_{\rm S}$, IO, upper $\theta_{23}$ octant,'+'\n'+r'JUNO + DUNE',
    r'$(1:2:0)_{\rm S}$, IO, upper $\theta_{23}$ octant,'+'\n'+r'JUNO + HK',
    r'$(1:2:0)_{\rm S}$, IO, upper $\theta_{23}$ octant,'+'\n'+r'NuFit 5.0 + JUNO',
    r'$(1:2:0)_{\rm S}$, IO, upper $\theta_{23}$ octant,'+'\n'+r'NuFit 5.0',
    r'$(1:2:0)_{\rm S}$, NO, lower $\theta_{23}$ octant,'+'\n'+r'DUNE',
    r'$(1:2:0)_{\rm S}$, NO, lower $\theta_{23}$ octant,'+'\n'+r'HK',
    r'$(1:2:0)_{\rm S}$, NO, lower $\theta_{23}$ octant,'+'\n'+r'JUNO + DUNE',
    r'$(1:2:0)_{\rm S}$, NO, lower $\theta_{23}$ octant,'+'\n'+r'JUNO + HK',
    r'$(1:2:0)_{\rm S}$, NO, $\theta_{23} = 45^\circ$,'+'\n'+r'DUNE',
    r'$(1:2:0)_{\rm S}$, NO, $\theta_{23} = 45^\circ$,'+'\n'+r'HK',
    r'$(1:2:0)_{\rm S}$, NO, $\theta_{23} = 45^\circ$,'+'\n'+r'JUNO + DUNE',
    r'$(1:2:0)_{\rm S}$, NO, $\theta_{23} = 45^\circ$,'+'\n'+r'JUNO + HK',
    r'$(1:2:0)_{\rm S}$, NO, upper $\theta_{23}$ octant,'+'\n'+r'DUNE',
    r'$(1:2:0)_{\rm S}$, NO, upper $\theta_{23}$ octant,'+'\n'+r'HK',
    r'$(1:2:0)_{\rm S}$, NO, upper $\theta_{23}$ octant,'+'\n'+r'JUNO + DUNE',
    r'$(1:2:0)_{\rm S}$, NO, upper $\theta_{23}$ octant,'+'\n'+r'JUNO + HK',
    r'$(1:2:0)_{\rm S}$, NO, upper $\theta_{23}$ octant,'+'\n'+r'NuFit 5.0 + JUNO',
    r'$(1:2:0)_{\rm S}$, NO, upper $\theta_{23}$ octant,'+'\n'+r'NuFit 5.0']

fS = [
    '010', '010', '010', '010', '010', '010', '010', '010', '010', '010', '010',
    '010', '010', '010', 
    '010', '010', '010', '010', '010', '010', '010', '010', '010', '010', '010',
    '010', '010', '010', 
    '100', '100', '100', '100', '100', '100', '100', '100', '100', '100', '100',
    '100', '100', '100', 
    '100', '100', '100', '100', '100', '100', '100', '100', '100', '100', '100',
    '100', '100', '100', 
    '120', '120', '120', '120', '120', '120', '120', '120', '120', '120', '120',
    '120', '120', '120', 
    '120', '120', '120', '120', '120', '120', '120', '120', '120', '120', '120',
    '120', '120', '120']

theta23 = [
    'lower', 'lower', 'lower', 'lower', 'max', 'max', 'max', 'max', 
    'upper', 'upper', 'upper', 'upper', 'upper', 'upper',
    'lower', 'lower', 'lower', 'lower', 'max', 'max', 'max', 'max', 
    'upper', 'upper', 'upper', 'upper', 'upper', 'upper',
    'lower', 'lower', 'lower', 'lower', 'max', 'max', 'max', 'max', 
    'upper', 'upper', 'upper', 'upper', 'upper', 'upper',
    'lower', 'lower', 'lower', 'lower', 'max', 'max', 'max', 'max', 
    'upper', 'upper', 'upper', 'upper', 'upper', 'upper',
    'lower', 'lower', 'lower', 'lower', 'max', 'max', 'max', 'max', 
    'upper', 'upper', 'upper', 'upper', 'upper', 'upper',
    'lower', 'lower', 'lower', 'lower', 'max', 'max', 'max', 'max', 
    'upper', 'upper', 'upper', 'upper', 'upper', 'upper']

mass_order = [
   'IO', 'IO', 'IO', 'IO', 'IO', 'IO', 'IO', 'IO', 'IO', 'IO', 'IO', 'IO', 'IO', 
   'IO', 
   'NO', 'NO', 'NO', 'NO', 'NO', 'NO', 'NO', 'NO', 'NO', 'NO', 'NO', 'NO', 'NO', 
   'NO', 
   'IO', 'IO', 'IO', 'IO', 'IO', 'IO', 'IO', 'IO', 'IO', 'IO', 'IO', 'IO', 'IO', 
   'IO', 
   'NO', 'NO', 'NO', 'NO', 'NO', 'NO', 'NO', 'NO', 'NO', 'NO', 'NO', 'NO', 'NO', 
   'NO', 
   'IO', 'IO', 'IO', 'IO', 'IO', 'IO', 'IO', 'IO', 'IO', 'IO', 'IO', 'IO', 'IO', 
   'IO', 
   'NO', 'NO', 'NO', 'NO', 'NO', 'NO', 'NO', 'NO', 'NO', 'NO', 'NO', 'NO', 'NO', 
   'NO']

for i, fname in enumerate(fnames):
    print(fname+'_contours'+'... ')
    plot_ternary_flavor_ratios_fixed_fS_separate(PATH_OUT+fname+'_contours', 
        plot_labels[i], fS=fS[i], order=mass_order[i], theta23=theta23[i],
        output_format='png')
"""


"""
# Varying all flavor compositions at the source, different NuFit versions (v1)

PATH_IN = os.getcwd()+'/flavor_region_output_nufit_versions/all_fS/'
PATH_OUT = os.getcwd()+'/flavor_region_output_nufit_versions_post/all_fS/'
if not os.path.exists(PATH_OUT):
    os.mkdir(PATH_OUT)
fnames = sorted(filenames_no_ext(PATH_IN)) # Keep the same file names throughout
# print(np.array(fnames))
# quit()

plot_labels = [
    # r'IO'+'\n'+r'NuFit 1.0', 
    # r'IO'+'\n'+r'NuFit 1.1',
    # r'IO'+'\n'+r'NuFit 1.2',
    # r'IO'+'\n'+r'NuFit 1.3',
    r'IO'+'\n'+r'NuFit 2.0',
    r'IO'+'\n'+r'NuFit 2.1',
    r'IO'+'\n'+r'NuFit 2.2',
    r'IO'+'\n'+r'NuFit 3.0',
    r'IO'+'\n'+r'NuFit 3.1',
    r'IO'+'\n'+r'NuFit 3.2',
    r'IO'+'\n'+r'NuFit 4.0',
    r'IO'+'\n'+r'NuFit 4.1',
    r'IO'+'\n'+r'NuFit 5.0',
    # r'NO'+'\n'+r'NuFit 1.0', 
    # r'NO'+'\n'+r'NuFit 1.1',
    # r'NO'+'\n'+r'NuFit 1.2',
    # r'NO'+'\n'+r'NuFit 1.3',
    r'NO'+'\n'+r'NuFit 2.0',
    r'NO'+'\n'+r'NuFit 2.1',
    r'NO'+'\n'+r'NuFit 2.2',
    r'NO'+'\n'+r'NuFit 3.0',
    r'NO'+'\n'+r'NuFit 3.1',
    r'NO'+'\n'+r'NuFit 3.2',
    r'NO'+'\n'+r'NuFit 4.0',
    r'NO'+'\n'+r'NuFit 4.1',
    r'NO'+'\n'+r'NuFit 5.0'
]

nufit_version = [
    # 'nufit_v1.0',
    # 'nufit_v1.1',
    # 'nufit_v1.2',
    # 'nufit_v1.3',
    'nufit_v2.0',
    'nufit_v2.1',
    'nufit_v2.2',
    'nufit_v3.0',
    'nufit_v3.1',
    'nufit_v3.2',
    'nufit_v4.0',
    'nufit_v4.1',
    'nufit_v5.0',
    # 'nufit_v1.0',
    # 'nufit_v1.1',
    # 'nufit_v1.2',
    # 'nufit_v1.3',
    'nufit_v2.0',
    'nufit_v2.1',
    'nufit_v2.2',
    'nufit_v3.0',
    'nufit_v3.1',
    'nufit_v3.2',
    'nufit_v4.0',
    'nufit_v4.1',
    'nufit_v5.0'
]

mass_order = [
    # 'IO',
    # 'IO',
    # 'IO',
    # 'IO',
    'IO',
    'IO',
    'IO',
    'IO',
    'IO',
    'IO',
    'IO',
    'IO',
    'IO',
    'NO',
    # 'NO',
    # 'NO',
    # 'NO',
    # 'NO',
    'NO',
    'NO',
    'NO',
    'NO',
    'NO',
    'NO',
    'NO',
    'NO']

for i, fname in enumerate(fnames):
    print(fname+'_contours'+'... ')
    plot_ternary_flavor_ratios_all_fS_v1(PATH_OUT+fname+'_contours', 
        plot_labels[i], order=mass_order[i], nufit_version=nufit_version[i],
        theta23=None, year=None, fixed_deltacp=None, output_format='png')
"""


"""
# Varying all flavor compositions at the source, fixed delta_CP (v1)

PATH_IN = os.getcwd()+'/flavor_region_output_fixed_delta/all_fS/'
PATH_OUT = os.getcwd()+'/flavor_region_output_fixed_delta_post/all_fS/'
if not os.path.exists(PATH_OUT):
    os.mkdir(PATH_OUT)
fnames = sorted(filenames_no_ext(PATH_IN)) # Keep the same file names throughout
# print(np.array(fnames))
# quit()

plot_labels = [
    r'IO, upper $\theta_{23}$ octant, $\delta_{\rm CP} = 0$,'+'\n'+r'DUNE',    
    r'IO, upper $\theta_{23}$ octant, $\delta_{\rm CP} = 0$,'+'\n'+r'HK',    
    r'IO, upper $\theta_{23}$ octant, $\delta_{\rm CP} = 0$,'+'\n'+r'JUNO + DUNE',    
    r'IO, upper $\theta_{23}$ octant, $\delta_{\rm CP} = 0$,'+'\n'+r'JUNO + HK',    
    r'IO, upper $\theta_{23}$ octant, $\delta_{\rm CP} = \pi$,'+'\n'+r'DUNE',    
    r'IO, upper $\theta_{23}$ octant, $\delta_{\rm CP} = \pi$,'+'\n'+r'HK',    
    r'IO, upper $\theta_{23}$ octant, $\delta_{\rm CP} = \pi$,'+'\n'+r'JUNO + DUNE',    
    r'IO, upper $\theta_{23}$ octant, $\delta_{\rm CP} = \pi$,'+'\n'+r'JUNO + HK',    
    r'IO, upper $\theta_{23}$ octant, $\delta_{\rm CP} = \frac{3\pi}{2}$,'+'\n'+r'DUNE',    
    r'IO, upper $\theta_{23}$ octant, $\delta_{\rm CP} = \frac{3\pi}{2}$,'+'\n'+r'HK',    
    r'IO, upper $\theta_{23}$ octant, $\delta_{\rm CP} = \frac{3\pi}{2}$,'+'\n'+r'JUNO + DUNE',    
    r'IO, upper $\theta_{23}$ octant, $\delta_{\rm CP} = \frac{3\pi}{2}$,'+'\n'+r'JUNO + HK',    
    r'NO, upper $\theta_{23}$ octant, $\delta_{\rm CP} = 0$,'+'\n'+r'DUNE',    
    r'NO, upper $\theta_{23}$ octant, $\delta_{\rm CP} = 0$,'+'\n'+r'HK',    
    r'NO, upper $\theta_{23}$ octant, $\delta_{\rm CP} = 0$,'+'\n'+r'JUNO + DUNE',    
    r'NO, upper $\theta_{23}$ octant, $\delta_{\rm CP} = 0$,'+'\n'+r'JUNO + HK',    
    r'NO, upper $\theta_{23}$ octant, $\delta_{\rm CP} = \frac{\pi}{2}$,'+'\n'+r'DUNE',    
    r'NO, upper $\theta_{23}$ octant, $\delta_{\rm CP} = \frac{\pi}{2}$,'+'\n'+r'HK',    
    r'NO, upper $\theta_{23}$ octant, $\delta_{\rm CP} = \frac{\pi}{2}$,'+'\n'+r'JUNO + DUNE',    
    r'NO, upper $\theta_{23}$ octant, $\delta_{\rm CP} = \frac{\pi}{2}$,'+'\n'+r'JUNO + HK',    
    r'NO, upper $\theta_{23}$ octant, $\delta_{\rm CP} = \pi$,'+'\n'+r'DUNE',    
    r'NO, upper $\theta_{23}$ octant, $\delta_{\rm CP} = \pi$,'+'\n'+r'HK',    
    r'NO, upper $\theta_{23}$ octant, $\delta_{\rm CP} = \pi$,'+'\n'+r'JUNO + DUNE',    
    r'NO, upper $\theta_{23}$ octant, $\delta_{\rm CP} = \pi$,'+'\n'+r'JUNO + HK']

theta23 = [
    'upper',
    'upper',
    'upper',
    'upper',
    'upper',
    'upper',
    'upper',
    'upper',
    'upper',
    'upper',
    'upper',
    'upper',
    'upper',
    'upper',
    'upper',
    'upper',
    'upper',
    'upper',
    'upper',
    'upper',
    'upper',
    'upper',
    'upper',
    'upper']

mass_order = [
    'IO',
    'IO',
    'IO',
    'IO',
    'IO',
    'IO',
    'IO',
    'IO',
    'IO',
    'IO',
    'IO',
    'IO',
    'NO',
    'NO',
    'NO',
    'NO',
    'NO',
    'NO',
    'NO',
    'NO',
    'NO',
    'NO',
    'NO',
    'NO']

year = [
    '2040',
    '2030',
    '2040',
    '2030',
    '2040',
    '2030',
    '2040',
    '2030',
    '2040',
    '2030',
    '2040',
    '2030',
    '2040',
    '2030',
    '2040',
    '2030',
    '2040',
    '2030',
    '2040',
    '2030',
    '2040',
    '2030',
    '2040',
    '2030']

fixed_deltacp = [
    0.0,
    0.0,
    0.0,
    0.0,
    np.pi,        
    np.pi,        
    np.pi,        
    np.pi,        
    1.5*np.pi,        
    1.5*np.pi,        
    1.5*np.pi,        
    1.5*np.pi,        
    0.0,
    0.0,
    0.0,
    0.0,
    0.5*np.pi,        
    0.5*np.pi,        
    0.5*np.pi,        
    0.5*np.pi,        
    np.pi,        
    np.pi,        
    np.pi,        
    np.pi,        
]

for i, fname in enumerate(fnames):
    print(fname+'_contours'+'... ')
    plot_ternary_flavor_ratios_all_fS_v1(PATH_OUT+fname+'_contours', 
        plot_labels[i], order=mass_order[i], theta23=theta23[i],
        year=year[i], fixed_deltacp=fixed_deltacp[i], output_format='png')
"""



"""
# Non-unitarity: varying all flavor compositions at the source (v1)

PATH_IN = os.getcwd()+'/flavor_region_output_nonunitary/all_fS/'
PATH_OUT = os.getcwd()+'/flavor_region_output_nonunitary_post/all_fS/'
if not os.path.exists(PATH_OUT):
    os.mkdir(PATH_OUT)
fnames = sorted(filenames_no_ext(PATH_IN)) # Keep the same file names throughout
# print(np.array(fnames))
# quit()

plot_labels = [
    r'Non-unitarity, 2020',
    r'Non-unitarity, 2040'
]

theta23 = [
    'upper',
    'upper'
]

mass_order = [
    'NO',
    'NO'
]

year = [
    '2020',
    '2040',
]

for i, fname in enumerate(fnames):
    print(fname+'_contours'+'... ')
    plot_ternary_flavor_ratios_all_fS_v1(PATH_OUT+fname+'_contours', 
        plot_labels[i], order=mass_order[i], theta23=theta23[i],
        year=year[i], show_benchmarks=True, show_benchmark_astro_labels=True,
        show_benchmark_decay_labels=False, output_format='png')
"""




"""
# Neutrino decay: varying all flavor compositions at the source (v1)

PATH_IN = os.getcwd()+'/flavor_region_output_decay/all_fS/'
PATH_OUT = os.getcwd()+'/flavor_region_output_decay_post/all_fS/'
if not os.path.exists(PATH_OUT):
    os.mkdir(PATH_OUT)
fnames = sorted(filenames_no_ext(PATH_IN)) # Keep the same file names throughout
# print(np.array(fnames))
# quit()

plot_labels = [
    r'$\nu$ decay, IO, upper $\theta_{23}$ octant,'+'\n'+'HK + JUNO',
    r'$\nu$ decay, IO, upper $\theta_{23}$ octant,'+'\n'+'NuFit 5.0',
    r'$\nu$ decay, NO, upper $\theta_{23}$ octant,'+'\n'+'HK + JUNO',
    r'$\nu$ decay, NO, upper $\theta_{23}$ octant,'+'\n'+'NuFit 5.0',
]

theta23 = [
    'upper',
    'upper',
    'upper',
    'upper'
]

mass_order = [
    'IO',
    'IO',
    'NO',
    'NO'
]

year = [
    '2040',
    '2020',
    '2040',
    '2020'
]

for i, fname in enumerate(fnames):
    print(fname+'_contours'+'... ')
    plot_ternary_flavor_ratios_all_fS_v1(PATH_OUT+fname+'_contours', 
        plot_labels[i], order=mass_order[i], theta23=theta23[i],
        year=year[i], show_benchmarks=True, show_benchmark_astro_labels=False,
        show_benchmark_decay_labels=True, output_format='png')
"""


"""
# Varying all flavor compositions at the source (v1)

# Data for combined experiments
PATH_IN = os.getcwd()+'/flavor_region_output_JUNO_DUNE_HK/all_fS/'
PATH_OUT = os.getcwd()+'/flavor_region_output_JUNO_DUNE_HK_post/all_fS/'
if not os.path.exists(PATH_OUT):
    os.mkdir(PATH_OUT)
fnames = sorted(filenames_no_ext(PATH_IN)) # Keep the same file names throughout
# print(np.array(fnames))
# quit()
# fnames = ['output_NO_upper_NUFIT']

plot_labels = [
    r'IO, lower $\theta_{23}$ octant,'+'\n'+r'JUNO + DUNE + HK',    
    r'IO, $\theta_{23} = 45^\circ$,'+'\n'+r'JUNO + DUNE + HK',    
    r'IO, upper $\theta_{23}$ octant,'+'\n'+r'JUNO + DUNE + HK',    
    r'NO, lower $\theta_{23}$ octant,'+'\n'+r'JUNO + DUNE + HK',    
    r'NO, $\theta_{23} = 45^\circ$,'+'\n'+r'JUNO + DUNE + HK',    
    r'NO, upper $\theta_{23}$ octant,'+'\n'+r'JUNO + DUNE + HK',    
    r'NO, upper $\theta_{23}$ octant,'+'\n'+r'NuFit 5.0'
]

theta23 = [
    'lower',
    'max',
    'upper',
    'lower',
    'max',
    'upper',
    'upper'
]

mass_order = [
    'IO',
    'IO',
    'IO',
    'NO',
    'NO',
    'NO',
    'NO'
]

year = [
    '2040',
    '2040',
    '2040',
    '2040',
    '2040',
    '2040',
    '2020'
]

for i, fname in enumerate(fnames):
    print(fname+'_contours'+'... ')
    plot_ternary_flavor_ratios_all_fS_v1(PATH_OUT+fname+'_contours', 
        plot_labels[i], order=mass_order[i], theta23=theta23[i],
        year=year[i], output_format='png')
    plot_ternary_flavor_ratios_all_fS_v1(PATH_OUT+fname+'_contours', 
        plot_labels[i], order=mass_order[i], theta23=theta23[i],
        year=year[i], output_format='pdf')
"""


"""
# Varying all flavor compositions at the source (v1)

# Data for separate experiments
PATH_IN = os.getcwd()+'/flavor_region_output/'
PATH_OUT = os.getcwd()+'/flavor_region_output_post/'
if not os.path.exists(PATH_OUT):
    os.mkdir(PATH_OUT)
fnames = sorted(filenames_no_ext(PATH_IN)) # Keep the same file names throughout
# print(np.array(fnames))
# quit()

plot_labels = [
    r'IO, lower $\theta_{23}$ octant,'+'\n'+r'DUNE',    
    r'IO, lower $\theta_{23}$ octant,'+'\n'+r'DUNE + JUNO',
    r'IO, lower $\theta_{23}$ octant,'+'\n'+'HK',
    r'IO, lower $\theta_{23}$ octant,'+'\n'+'HK + JUNO',
    r'IO, $\theta_{23} = 45^\circ$,'+'\n'+'DUNE',
    r'IO, $\theta_{23} = 45^\circ$,'+'\n'+'DUNE + JUNO',
    r'IO, $\theta_{23} = 45^\circ$,'+'\n'+'HK',
    r'IO, $\theta_{23} = 45^\circ$,'+'\n'+'HK + JUNO',
    r'IO, upper $\theta_{23}$ octant,'+'\n'+'DUNE',
    r'IO, upper $\theta_{23}$ octant,'+'\n'+'DUNE + JUNO',
    r'IO, upper $\theta_{23}$ octant,'+'\n'+r'HK',
    r'IO, upper $\theta_{23}$ octant,'+'\n'+'HK + JUNO',
    r'IO, upper $\theta_{23}$ octant,'+'\n'+'NuFit 5.0',
    r'IO, upper $\theta_{23}$ octant,'+'\n'+'NuFit 5.0 + JUNO',
    r'IO, upper $\theta_{23}$ octant,'+'\n'+r'Gaussian distribution',
    r'NO, lower $\theta_{23}$ octant,'+'\n'+r'DUNE',    
    r'NO, lower $\theta_{23}$ octant,'+'\n'+r'DUNE + JUNO',
    r'NO, lower $\theta_{23}$ octant,'+'\n'+'HK',
    r'NO, lower $\theta_{23}$ octant,'+'\n'+'HK + JUNO',
    r'NO, $\theta_{23} = 45^\circ$,'+'\n'+'DUNE',
    r'NO, $\theta_{23} = 45^\circ$,'+'\n'+'DUNE + JUNO',
    r'NO, $\theta_{23} = 45^\circ$,'+'\n'+'HK',
    r'NO, $\theta_{23} = 45^\circ$,'+'\n'+'HK + JUNO',
    r'NO, upper $\theta_{23}$ octant,'+'\n'+'DUNE',
    r'NO, upper $\theta_{23}$ octant,'+'\n'+'DUNE + JUNO',
    r'NO, upper $\theta_{23}$ octant,'+'\n'+r'HK',
    r'NO, upper $\theta_{23}$ octant,'+'\n'+'HK + JUNO',
    r'NO, upper $\theta_{23}$ octant,'+'\n'+'NuFit 5.0',
    r'NO, upper $\theta_{23}$ octant,'+'\n'+'NuFit 5.0 + JUNO',
    r'NO, upper $\theta_{23}$ octant,'+'\n'+r'Gaussian distribution']

theta23 = [
    'lower',
    'lower',
    'lower',
    'lower',
    'max',
    'max',
    'max',
    'max',
    'upper',
    'upper',
    'upper',
    'upper',
    'upper',
    'upper',
    'upper',
    'lower',
    'lower',
    'lower',
    'lower',
    'max',
    'max',
    'max',
    'max',
    'upper',
    'upper',
    'upper',
    'upper',
    'upper',
    'upper',
    'upper']

mass_order = [
    'IO',
    'IO',
    'IO',
    'IO',
    'IO',
    'IO',
    'IO',
    'IO',
    'IO',
    'IO',
    'IO',
    'IO',
    'IO',
    'IO',
    'IO',
    'NO',
    'NO',
    'NO',
    'NO',
    'NO',
    'NO',
    'NO',
    'NO',
    'NO',
    'NO',
    'NO',
    'NO',
    'NO',
    'NO',
    'NO']

year = [
    '2040',
    '2040',
    '2030',
    '2030',
    '2040',
    '2040',
    '2030',
    '2030',
    '2040',
    '2040',
    '2030',
    '2030',
    '2020',
    '2030',
    '2020',
    '2040',
    '2040',
    '2030',
    '2030',
    '2040',
    '2040',
    '2030',
    '2030',
    '2040',
    '2040',
    '2030',
    '2030',
    '2020',
    '2030',
    '2020']

for i, fname in enumerate(fnames):
    print(fname+'_contours'+'... ')
    plot_ternary_flavor_ratios_all_fS_v1(PATH_OUT+fname+'_contours', 
        plot_labels[i], order=mass_order[i], theta23=theta23[i],
        year=year[i], output_format='png')
"""


"""
# Varying all flavor compositions at the source (v0)

PATH_IN = os.getcwd()+'/flavor_region_output/'
PATH_OUT = os.getcwd()+'/flavor_region_output_post/'
if not os.path.exists(PATH_OUT):
    os.mkdir(PATH_OUT)
fnames = filenames_no_ext(PATH_IN) # Keep the same file names throughout
# print(fnames)
# quit()

plot_labels = [
    r'IO, lower $\theta_{23}$ octant,'+'\n'+r'DUNE + JUNO',
    r'NO, lower $\theta_{23}$ octant,'+'\n'+'HK',
    r'IO, upper $\theta_{23}$ octant,'+'\n'+r'HK',
    r'NO, upper $\theta_{23}$ octant,'+'\n'+r'Gaussian distribution',
    r'IO, upper $\theta_{23}$ octant,'+'\n'+'DUNE + JUNO',
    r'NO, $\theta_{23} = 45^\circ$,'+'\n'+'HK + JUNO',
    r'IO, upper $\theta_{23}$ octant,'+'\n'+r'Gaussian distribution',
    r'NO, upper $\theta_{23}$ octant,'+'\n'+'HK',
    r'NO, upper $\theta_{23}$ octant,'+'\n'+'DUNE + JUNO',
    r'IO, lower $\theta_{23}$ octant,'+'\n'+'HK + JUNO',
    r'NO, upper $\theta_{23}$ octant,'+'\n'+'HK + JUNO',
    r'NO, $\theta_{23} = 45^\circ$,'+'\n'+'HK',
    r'IO, $\theta_{23} = 45^\circ$,'+'\n'+'HK',
    r'IO, $\theta_{23} = 45^\circ$,'+'\n'+'DUNE + JUNO',
    r'NO, lower $\theta_{23}$ octant,'+'\n'+'HK + JUNO',
    r'NO, lower $\theta_{23}$ octant,'+'\n'+'DUNE + JUNO',
    r'IO, upper $\theta_{23}$ octant,'+'\n'+'NuFit 5.0',
    r'NO, upper $\theta_{23}$ octant,'+'\n'+'NuFit 5.0 + JUNO',
    r'IO, lower $\theta_{23}$ octant,'+'\n'+'HK',
    r'IO, upper $\theta_{23}$ octant,'+'\n'+'NuFit 5.0 + JUNO',
    r'NO, $\theta_{23} = 45^\circ$,'+'\n'+'DUNE + JUNO',
    r'NO, upper $\theta_{23}$ octant,'+'\n'+'NuFit 5.0',
    r'IO, upper $\theta_{23}$ octant,'+'\n'+'HK+JUNO',
    r'IO, $\theta_{23} = 45^\circ$,'+'\n'+'HK + JUNO']

theta23 = [
    'lower',
    'lower',
    'upper',
    'upper',
    'upper',
    'max',
    'upper',
    'upper',
    'upper',
    'lower',
    'upper',
    'max',
    'max',
    'max',
    'lower',
    'lower',
    'upper',
    'upper',
    'lower',
    'upper',
    'max',
    'upper',
    'upper',
    'max']

mass_order = [
    'IO',
    'NO',
    'IO',
    'NO', 
    'IO',
    'NO',
    'IO',
    'NO',
    'NO',
    'IO',
    'NO',
    'NO',
    'IO',
    'IO',
    'NO',
    'NO',
    'IO',
    'NO',
    'IO',
    'IO',
    'NO',
    'NO',
    'IO',
    'IO']

for i, fname in enumerate(fnames):
    print(fname+'_contours'+'... ')
    plot_ternary_flavor_ratios_all_fS_v0(PATH_OUT+fname+'_contours', 
        plot_labels[i], order=mass_order[i], theta23=theta23[i],
         output_format='png')
"""


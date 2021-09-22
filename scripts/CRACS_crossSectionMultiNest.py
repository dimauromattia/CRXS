#! /usr/bin/env python

import glob
import math
import sys
import os
import argparse


import  numpy                   as      np
from    operator                import  itemgetter

import  matplotlib              as      mpl
import  matplotlib.pyplot       as      plt
import  matplotlib.patches      as      mpatches
from    matplotlib              import  rc

import  multinest_plot_functions    as mnf
import  PlotFunctions               as f

from PlotFunctions import plot2D
#from PlotFunctions import profile

#   Argument Handling
#######################################################


parser = argparse.ArgumentParser(description='CRACS - Plotting the fit results of the Cross Section Scan by MultiNest.')
parser.add_argument('--step',     help='Step .',   action='store', dest='step', type=int, default=1)

args        = parser.parse_args()
step   	    = args.step

#######################################################




print_size=10
label_size=30


def plot_1D( xlabel, ylabel, xscale='log', yscale='log', sizex=1.3, sizey=1.0):

    global print_size, label_size
    plt.close('all')
    font_props = {"size":label_size}
    rc("font", **font_props)
    mpl.rcParams['axes.linewidth'] = 2
    mpl.rcParams['mathtext.fontset']='stixsans'

    fig     = plt.figure(figsize=(print_size*sizex, print_size*sizey))
    plot    = plt.subplot2grid((1, 1), (0, 0), rowspan=1)
    plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)
    
    plot.set_xlabel ( xlabel )
    plot.set_ylabel ( ylabel )
    
    plot.set_xscale ( xscale )
    plot.set_yscale ( yscale )

    plot.tick_params('both', length=20, width=2, which='major')
    plot.tick_params('both', length=10, width=1, which='minor')

    plot.grid(b=True, which='major', alpha=0.1, linestyle='-', linewidth=2)
    plot.tick_params(axis='both', pad=10)

    return plot, fig


def profile(file, Tn, type='prod', Tn_min=1e-1, Tn_max=1e5, y_fac=1. ):
    T_prod = Tn
    data = np.genfromtxt( file )
    if not type=='prod':
        data = np.transpose(data)
    min  = 1e90;
    col  = 1
    for i in range(1,len(data[0,:])):
        diff = (data[0,i]-T_prod)*(data[0,i]-T_prod)
        if diff<min:
            min = diff
            col = i
    x_i = 1
    for i in range(1,len(data[:,0])):
        if data[i,0]<Tn_max:
            x_i = i
    x = data[1:x_i,0 ]
    y = data[1:x_i,col]*y_fac
    return x,y



################################


if step == 1:
    prefix 		= ''
    
    param          = []
    
    if prefix=='':
        l = glob.glob('MultiNest/*.txt')
        if len(l):
            prefix = (l[0].split('.txt')[0]).split('/')[1]
    print 'Use prefix: ' + prefix
    
    
    multiNest               = np.genfromtxt('MultiNest/'+prefix+'.txt')
    all_names, all_ranges   = f.readNamesAndRanges('MultiNest/'+prefix+'.ranges_plot')
    
    print all_names
    print all_ranges
    
    nCparam = 0
    
    for i in range(len(all_names)):
        if 'C' in all_names[i]:
            nCparam += 1
        param.append(i)
    
    print param

    
    mnf.prepare_triangle( nCparam + 1 )
    mnf.set_error_counter_properties(   'color',    'black'  )
    mnf.set_error_counter_properties(   'cmap',     'Greys'  )
    mnf.draw_triangle(prefix, param[:nCparam], shift_y=1)
    mnf.draw_diagonal(prefix, param[:nCparam])
    
    plt.savefig('triangle_freq_cont_Cparm.pdf')

    sys.exit(0)

    mnf.prepare_triangle( len(param) )
    mnf.draw_triangle(prefix, param, type='chiSq_scatter')
    plt.savefig('triangle_chiSq_scatter.png')

    mnf.prepare_triangle( len(param) + 1 )
    mnf.draw_triangle(prefix, param, type='error_scatter', shift_y=1 )
    mnf.set_chi_square_properties('smoothing', 0)
    mnf.draw_diagonal(prefix, param)
    plt.savefig('triangle_error_scatter.png')
    mnf.set_chi_square_properties('smoothing', 10)

    mnf.prepare_triangle( len(param) + 1 )
    mnf.set_error_counter_properties(   'color',    'black'  )
    mnf.set_error_counter_properties(   'cmap',     'Greys'  )
    mnf.draw_triangle(prefix, param, shift_y=1)
    mnf.draw_diagonal(prefix, param)

    plt.savefig('triangle_freq_cont.pdf')



if step == 2:

    print "To do: Call CRACS_crossSectionFitter.py "




#if step == 3:
#
#    plot, fig = plot_1D (r'$\mathrm{T_\bar{p}\,\,[GeV]}$', r'$\mathrm{q \, T_{\bar{p}}^{2.7}  \,  \, [GeV^{1.7}m^{-3}s^{-1}]}}$',      'log', 'log' )
#
#
#    data_q = np.genfromtxt('sourceTerm.txt', skip_header=1)
#
#    plot.fill_between   ( data_q[:,0], data_q[:,2]*np.power(data_q[:,0], 2.7), data_q[:,3]*np.power(data_q[:,0], 2.7),   lw=1,  color='black'  ,   alpha=0.25,        label=r''   )
#    plot.fill_between   ( data_q[:,0], data_q[:,4]*np.power(data_q[:,0], 2.7), data_q[:,5]*np.power(data_q[:,0], 2.7),   lw=1,  color='black'  ,   alpha=0.25,        label=r''   )
#
#    plot.errorbar       ( data_q[:,0], data_q[:,1]*np.power(data_q[:,0], 2.7),                                           lw=2,  color='black'  ,   dashes=(     ),    label=r''   )
#
#    plot.set_ylim( (1e-37, 2e-33) )
#
#    plt.savefig( 'sourceTerm.png' )

def get_min_array( list ):
    y_min = np.zeros(len(list[0]))
    new = 0
    for l in list:
        for j in range(len(l)):
            if (l[j]<y_min[j] and l[j]>1e-90) or y_min[j]<1e-90:
                y_min[j] = l[j]
    return y_min

def get_max_array( list ):
    y_max = np.zeros(len(list[0]))
    new = 0
    for l in list:
        for j in range(len(l)):
            if ( l[j]>y_max[j] ):
                y_max[j] = l[j]
    return y_max



if step == 12:


    c0='black'
#    c1='#8A0808'
#    c2='#0B610B'
#    c3='#0B0B61'

    c1='#04B404'
    c2='#0489B1'
    c3='#0404B4'
    c4='#B40404'



#    plot2D('dT_pHe_pbar_LAB_diMauro12.txt',            file2='dT_Hep_pbar_LAB_diMauro12.txt',       resfile='ratio_pHe_Hep_diMauro.pdf',  Zlabel = 'Hep/pHe', Zmin=0.8, Zmax=1.2 )
#
#
    plot2D('sourceTerm/dT_pp_pbar_LAB___best.txt',          resfile='XS_pp.pdf'   , x_label=r'$\mathrm{T_{p}\quad [GeV]}$'  )
#    plot2D('dT_Hep_pbar_LAB_diMauro12.txt',          resfile='paper_dT_Hep_pbar_LAB_diMauro12.pdf'   , x_label=r'$\mathrm{T_{He}/4\quad [GeV]}$'  )
#    plot2D('dT_pHe_pbar_LAB_diMauro12.txt',          resfile='paper_dT_pHe_pbar_LAB_diMauro12.pdf'   , x_label=r'$\mathrm{T_{p}\quad [GeV]}$'  )

    XS_parameterizations = os.environ['CRACS']+'/data/CS_lab'
    print XS_parameterizations
    #
    #   p p
    #
    plot, fig = plot_1D (r'$\mathrm{T_{p} \quad [GeV]}$',      r'$\mathrm{d\sigma/dT_{\bar{p}}\quad[m^2/GeV]}$', 'log', 'log' )

    x1_1, y1_1 = profile( XS_parameterizations+'/dT_pp_pbar_LAB_diMauro12.txt',               0.5, type='prod', y_fac=40*2.3  )
    x1_2, y1_2 = profile( XS_parameterizations+'/dT_pp_pbar_LAB_WinklerWithHypWitNbar.txt',   0.5, type='prod', y_fac=40      )
    x1_3, y1_3 = profile( XS_parameterizations+'/dT_pp_pbar_LAB_KR_PPFRAG.txt',               0.5, type='prod', y_fac=40*2.3  )
    x1_4, y1_4 = profile( XS_parameterizations+'/dT_pp_pbar_LAB_Duperray.txt',                0.5, type='prod', y_fac=40*2.3  )
    x1_5, y1_5 = profile( XS_parameterizations+'/dT_pp_pbar_LAB_TanNg.txt',                   0.5, type='prod', y_fac=40*2.3  )
    
    x2_1, y2_1 = profile( XS_parameterizations+'/dT_pp_pbar_LAB_diMauro12.txt',               20., type='prod', y_fac=2.3 )
    x2_2, y2_2 = profile( XS_parameterizations+'/dT_pp_pbar_LAB_WinklerWithHypWitNbar.txt',   20., type='prod', y_fac=1.0 )
    x2_3, y2_3 = profile( XS_parameterizations+'/dT_pp_pbar_LAB_KR_PPFRAG.txt',               20., type='prod', y_fac=2.3 )
    x2_4, y2_4 = profile( XS_parameterizations+'/dT_pp_pbar_LAB_Duperray.txt',                20., type='prod', y_fac=2.3 )
    x2_5, y2_5 = profile( XS_parameterizations+'/dT_pp_pbar_LAB_TanNg.txt',                   20., type='prod', y_fac=2.3 )
    
    x3_1, y3_1 = profile( XS_parameterizations+'/dT_pp_pbar_LAB_diMauro12.txt',               200, type='prod', y_fac=2.3 )
    x3_2, y3_2 = profile( XS_parameterizations+'/dT_pp_pbar_LAB_WinklerWithHypWitNbar.txt',   200, type='prod', y_fac=1.0 )
    x3_3, y3_3 = profile( XS_parameterizations+'/dT_pp_pbar_LAB_KR_PPFRAG.txt',               200, type='prod', y_fac=2.3 )
    x3_4, y3_4 = profile( XS_parameterizations+'/dT_pp_pbar_LAB_Duperray.txt',                200, type='prod', y_fac=2.3 )
    x3_5, y3_5 = profile( XS_parameterizations+'/dT_pp_pbar_LAB_TanNg.txt',                   200, type='prod', y_fac=2.3 )

    x4_1, y4_1 = profile( XS_parameterizations+'/dT_pp_pbar_LAB_diMauro12.txt',               500, type='prod', y_fac=2.3 )
    x4_2, y4_2 = profile( XS_parameterizations+'/dT_pp_pbar_LAB_WinklerWithHypWitNbar.txt',   500, type='prod', y_fac=1.0 )
    x4_3, y4_3 = profile( XS_parameterizations+'/dT_pp_pbar_LAB_KR_PPFRAG.txt',               500, type='prod', y_fac=2.3 )
    x4_4, y4_4 = profile( XS_parameterizations+'/dT_pp_pbar_LAB_Duperray.txt',                500, type='prod', y_fac=2.3 )
    x4_5, y4_5 = profile( XS_parameterizations+'/dT_pp_pbar_LAB_TanNg.txt',                   500, type='prod', y_fac=2.3 )


    x1,   y1   = profile('sourceTerm/dT_pp_pbar_LAB___best.txt',                            0.5, type='prod', y_fac=40*3.1  )
    x2,   y2   = profile('sourceTerm/dT_pp_pbar_LAB___best.txt',                            20., type='prod', y_fac=3.1     )
    x3,   y3   = profile('sourceTerm/dT_pp_pbar_LAB___best.txt',                            200, type='prod', y_fac=3.0     )
    x4,   y4   = profile('sourceTerm/dT_pp_pbar_LAB___best.txt',                            500, type='prod', y_fac=3.0     )

    plot.fill_between( x1_1, get_max_array([y1_1,y1_2,y1_3,y1_4,y1_5]), get_min_array([y1_1,y1_2,y1_3,y1_4,y1_5]),   lw=2,  color=c1  , alpha=0.25, label=r''   )
    plot.fill_between( x2_1, get_max_array([y2_1,y2_2,y2_3,y2_4,y2_5]), get_min_array([y2_1,y2_2,y2_3,y2_4,y2_5]),   lw=2,  color=c2  , alpha=0.25, label=r''   )
    plot.fill_between( x3_1, get_max_array([y3_1,y3_2,y3_3,y3_4,y3_5]), get_min_array([y3_1,y3_2,y3_3,y3_4,y3_5]),   lw=2,  color=c3  , alpha=0.25, label=r''   )
    plot.fill_between( x3_1, get_max_array([y4_1,y4_2,y4_3,y4_4,y4_5]), get_min_array([y4_1,y4_2,y4_3,y4_4,y4_5]),   lw=2,  color=c4  , alpha=0.25, label=r''   )
    
    plot.plot(x1,y1, color=c1, lw=3)
    plot.plot(x2,y2, color=c2, lw=3)
    plot.plot(x3,y3, color=c3, lw=3)
    plot.plot(x4,y4, color=c4, lw=3)

    plot.set_xlim( (5, 1e5)         )
    plot.set_ylim( (5e-36, 5e-31)   )

    plot.text    (1.2e1, 2e-33 , r'$\mathrm{T_{\bar{p}}= 0.5\,GeV\,\,(\times\, 40) }$',   color=c1, rotation=85 )
    plot.text    (0.5e2, 2e-34 , r'$\mathrm{T_{\bar{p}}=                    20\,GeV}$',   color=c2, rotation=80 )
    plot.text    (2.0e2, 2e-34 , r'$\mathrm{T_{\bar{p}}=                   200\,GeV}$',   color=c3, rotation=80 )
    plot.text    (  6e2, 2e-34 , r'$\mathrm{T_{\bar{p}}=                   0.5\,TeV}$',   color=c4, rotation=75 )
    
    
    plt.savefig('XS_comparison_pp_prod.png')
    
    
    
    
    plot, fig = plot_1D (r'$\mathrm{T_{p} \quad [GeV]}$',      r'$\mathrm{d\sigma/dT_{\bar{p}}\quad[m^2/GeV]}$', 'log', 'log' )
    
    x1_1, y1_1 = profile( XS_parameterizations+'/dT_pp_pbar_LAB_diMauro12.txt',                20, type='proj', y_fac=2.3  )
    x1_2, y1_2 = profile( XS_parameterizations+'/dT_pp_pbar_LAB_WinklerWithHypWitNbar.txt',    20, type='proj', y_fac=1.0  )
    x1_3, y1_3 = profile( XS_parameterizations+'/dT_pp_pbar_LAB_KR_PPFRAG.txt',                20, type='proj', y_fac=2.3  )
    x1_4, y1_4 = profile( XS_parameterizations+'/dT_pp_pbar_LAB_Duperray.txt',                 20, type='proj', y_fac=2.3  )
    x1_5, y1_5 = profile( XS_parameterizations+'/dT_pp_pbar_LAB_TanNg.txt',                    20, type='proj', y_fac=2.3  )
    
    x2_1, y2_1 = profile( XS_parameterizations+'/dT_pp_pbar_LAB_diMauro12.txt',               450, type='proj', y_fac=2.3 )
    x2_2, y2_2 = profile( XS_parameterizations+'/dT_pp_pbar_LAB_WinklerWithHypWitNbar.txt',   450, type='proj', y_fac=1.0 )
    x2_3, y2_3 = profile( XS_parameterizations+'/dT_pp_pbar_LAB_KR_PPFRAG.txt',               450, type='proj', y_fac=2.3 )
    x2_4, y2_4 = profile( XS_parameterizations+'/dT_pp_pbar_LAB_Duperray.txt',                450, type='proj', y_fac=2.3 )
    x2_5, y2_5 = profile( XS_parameterizations+'/dT_pp_pbar_LAB_TanNg.txt',                   450, type='proj', y_fac=2.3 )
    
    x3_1, y3_1 = profile( XS_parameterizations+'/dT_pp_pbar_LAB_diMauro12.txt',              6500, type='proj', y_fac=5*2.3 )
    x3_2, y3_2 = profile( XS_parameterizations+'/dT_pp_pbar_LAB_WinklerWithHypWitNbar.txt',  6500, type='proj', y_fac=5*1.0 )
    x3_3, y3_3 = profile( XS_parameterizations+'/dT_pp_pbar_LAB_KR_PPFRAG.txt',              6500, type='proj', y_fac=5*2.3 )
    x3_4, y3_4 = profile( XS_parameterizations+'/dT_pp_pbar_LAB_Duperray.txt',               6500, type='proj', y_fac=5*2.3 )
    x3_5, y3_5 = profile( XS_parameterizations+'/dT_pp_pbar_LAB_TanNg.txt',                  6500, type='proj', y_fac=5*2.3 )
    
    x1,   y1   = profile('sourceTerm/dT_pp_pbar_LAB___best.txt',                             20, type='proj', y_fac=3.1     )
    x2,   y2   = profile('sourceTerm/dT_pp_pbar_LAB___best.txt',                            450, type='proj', y_fac=3.1     )
    x3,   y3   = profile('sourceTerm/dT_pp_pbar_LAB___best.txt',                           6500, type='proj', y_fac=3.1*5   )
    
    plot.fill_between( x1_1, get_max_array([y1_1,y1_2,y1_3,y1_4,y1_5]), get_min_array([y1_1,y1_2,y1_3,y1_4,y1_5]),   lw=2,  color=c1  , alpha=0.25, label=r''   )
    plot.fill_between( x2_1, get_max_array([y2_1,y2_2,y2_3,y2_4,y2_5]), get_min_array([y2_1,y2_2,y2_3,y2_4,y2_5]),   lw=2,  color=c2  , alpha=0.25, label=r''   )
    plot.fill_between( x3_1, get_max_array([y3_1,y3_2,y3_3,y3_4,y3_5]), get_min_array([y3_1,y3_2,y3_3,y3_4,y3_5]),   lw=2,  color=c3  , alpha=0.25, label=r''   )
    
    plot.plot(x1,y1, color=c1, lw=3)
    plot.plot(x2,y2, color=c2, lw=3)
    plot.plot(x3,y3, color=c3, lw=3)
    
    plot.set_xlim( (1e-1, 1e4)         )
    plot.set_ylim( (5e-36, 5e-31)   )
    
    plot.text    ( 0.8e1, 2e-34 , r'$\mathrm{T_{p}=                   20\,GeV}$',   color=c1,     rotation=-80 )
    plot.text    (   2e2, 2e-34 , r'$\mathrm{T_{p}=                  450\,GeV}$',   color=c2,     rotation=-80 )
    plot.text    (   1e3, 2e-33 , r'$\mathrm{T_{p}=  6.5\,TeV\,\,(\times\,5) }$',   color=c3,     rotation=-70 )

    
    plt.savefig('XS_comparison_pp_proj.png')
    
    sys.exit(0)
    
    
    #
    #   p He
    #
    plot, fig = plot_1D (r'$\mathrm{T_{p} \quad [GeV]}$',      r'$\mathrm{d\sigma/dT_{\bar{p}}\quad[m^2/GeV]}$', 'log', 'log' )
    
    x1_1, y1_1 = profile( XS_parameterizations+'/dT_pHe_pbar_LAB_diMauro12.txt',               0.5, type='prod', y_fac=40*2.3  )
    x1_2, y1_2 = profile( XS_parameterizations+'/dT_pHe_pbar_LAB_WinklerWithHypWitNbar.txt',   0.5, type='prod', y_fac=40      )
    x1_3, y1_3 = profile( XS_parameterizations+'/dT_pHe_pbar_LAB_KR_PPFRAG.txt',               0.5, type='prod', y_fac=40*2.3  )
    x1_4, y1_4 = profile( XS_parameterizations+'/dT_pHe_pbar_LAB_Duperray.txt',                0.5, type='prod', y_fac=40*2.3  )
    x1_5, y1_5 = profile( XS_parameterizations+'/dT_pHe_pbar_LAB_TanNg.txt',                   0.5, type='prod', y_fac=40*2.3  )
    
    x2_1, y2_1 = profile( XS_parameterizations+'/dT_pHe_pbar_LAB_diMauro12.txt',               20., type='prod', y_fac=2.3 )
    x2_2, y2_2 = profile( XS_parameterizations+'/dT_pHe_pbar_LAB_WinklerWithHypWitNbar.txt',   20., type='prod', y_fac=1.0 )
    x2_3, y2_3 = profile( XS_parameterizations+'/dT_pHe_pbar_LAB_KR_PPFRAG.txt',               20., type='prod', y_fac=2.3 )
    x2_4, y2_4 = profile( XS_parameterizations+'/dT_pHe_pbar_LAB_Duperray.txt',                20., type='prod', y_fac=2.3 )
    x2_5, y2_5 = profile( XS_parameterizations+'/dT_pHe_pbar_LAB_TanNg.txt',                   20., type='prod', y_fac=2.3 )
    
    x3_1, y3_1 = profile( XS_parameterizations+'/dT_pHe_pbar_LAB_diMauro12.txt',               500, type='prod', y_fac=2.3 )
    x3_2, y3_2 = profile( XS_parameterizations+'/dT_pHe_pbar_LAB_WinklerWithHypWitNbar.txt',   500, type='prod', y_fac=1.0 )
    x3_3, y3_3 = profile( XS_parameterizations+'/dT_pHe_pbar_LAB_KR_PPFRAG.txt',               500, type='prod', y_fac=2.3 )
    x3_4, y3_4 = profile( XS_parameterizations+'/dT_pHe_pbar_LAB_Duperray.txt',                500, type='prod', y_fac=2.3 )
    x3_5, y3_5 = profile( XS_parameterizations+'/dT_pHe_pbar_LAB_TanNg.txt',                   500, type='prod', y_fac=2.3 )
    
    
    plot.fill_between( x1_1, get_max_array([y1_1,y1_2,y1_3,y1_4,y1_5]), get_min_array([y1_1,y1_2,y1_3,y1_4,y1_5]),   lw=2,  color=c1  , alpha=0.25, label=r''   )
    plot.fill_between( x2_1, get_max_array([y2_1,y2_2,y2_3,y2_4,y2_5]), get_min_array([y2_1,y2_2,y2_3,y2_4,y2_5]),   lw=2,  color=c2  , alpha=0.25, label=r''   )
    plot.fill_between( x3_1, get_max_array([y3_1,y3_2,y3_3,y3_4,y3_5]), get_min_array([y3_1,y3_2,y3_3,y3_4,y3_5]),   lw=2,  color=c3  , alpha=0.25, label=r''   )
    
    plot.set_xlim( (5, 1e5)         )
    plot.set_ylim( (5e-36, 5e-31)   )
    
    plot.text    (1.2e1, 2e-33 , r'$\mathrm{T_{\bar{p}}= 0.5\,GeV\,\,(\times\, 40) }$',   color=c1, rotation=85 )
    plot.text    (0.5e2, 2e-34 , r'$\mathrm{T_{\bar{p}}=                    20\,GeV}$',   color=c2, rotation=80 )
    plot.text    (  3e2, 2e-34 , r'$\mathrm{T_{\bar{p}}=                   0.5\,TeV}$',   color=c3, rotation=75 )
    
    
    plt.savefig('XS_comparison_pHe.png')

    #
    #   He p
    #
    plot, fig = plot_1D (r'$\mathrm{T_{p} \quad [GeV]}$',      r'$\mathrm{d\sigma/dT_{\bar{p}}\quad[m^2/GeV]}$', 'log', 'log' )
    
    x1_1, y1_1 = profile( XS_parameterizations+'/dT_Hep_pbar_LAB_diMauro12.txt',               0.5, type='prod', y_fac=40*2.3  )
    x1_2, y1_2 = profile( XS_parameterizations+'/dT_Hep_pbar_LAB_WinklerWithHypWitNbar.txt',   0.5, type='prod', y_fac=40      )
    x1_3, y1_3 = profile( XS_parameterizations+'/dT_Hep_pbar_LAB_KR_PPFRAG.txt',               0.5, type='prod', y_fac=40*2.3  )
    x1_4, y1_4 = profile( XS_parameterizations+'/dT_Hep_pbar_LAB_Duperray.txt',                0.5, type='prod', y_fac=40*2.3  )
    x1_5, y1_5 = profile( XS_parameterizations+'/dT_Hep_pbar_LAB_TanNg.txt',                   0.5, type='prod', y_fac=40*2.3  )
    
    x2_1, y2_1 = profile( XS_parameterizations+'/dT_Hep_pbar_LAB_diMauro12.txt',               20., type='prod', y_fac=2.3 )
    x2_2, y2_2 = profile( XS_parameterizations+'/dT_Hep_pbar_LAB_WinklerWithHypWitNbar.txt',   20., type='prod', y_fac=1.0 )
    x2_3, y2_3 = profile( XS_parameterizations+'/dT_Hep_pbar_LAB_KR_PPFRAG.txt',               20., type='prod', y_fac=2.3 )
    x2_4, y2_4 = profile( XS_parameterizations+'/dT_Hep_pbar_LAB_Duperray.txt',                20., type='prod', y_fac=2.3 )
    x2_5, y2_5 = profile( XS_parameterizations+'/dT_Hep_pbar_LAB_TanNg.txt',                   20., type='prod', y_fac=2.3 )
    
    x3_1, y3_1 = profile( XS_parameterizations+'/dT_Hep_pbar_LAB_diMauro12.txt',               500, type='prod', y_fac=2.3 )
    x3_2, y3_2 = profile( XS_parameterizations+'/dT_Hep_pbar_LAB_WinklerWithHypWitNbar.txt',   500, type='prod', y_fac=1.0 )
    x3_3, y3_3 = profile( XS_parameterizations+'/dT_Hep_pbar_LAB_KR_PPFRAG.txt',               500, type='prod', y_fac=2.3 )
    x3_4, y3_4 = profile( XS_parameterizations+'/dT_Hep_pbar_LAB_Duperray.txt',                500, type='prod', y_fac=2.3 )
    x3_5, y3_5 = profile( XS_parameterizations+'/dT_Hep_pbar_LAB_TanNg.txt',                   500, type='prod', y_fac=2.3 )
    
    
    plot.fill_between( x1_1, get_max_array([y1_1,y1_2,y1_3,y1_4,y1_5]), get_min_array([y1_1,y1_2,y1_3,y1_4,y1_5]),   lw=2,  color=c1  , alpha=0.25, label=r''   )
    plot.fill_between( x2_1, get_max_array([y2_1,y2_2,y2_3,y2_4,y2_5]), get_min_array([y2_1,y2_2,y2_3,y2_4,y2_5]),   lw=2,  color=c2  , alpha=0.25, label=r''   )
    plot.fill_between( x3_1, get_max_array([y3_1,y3_2,y3_3,y3_4,y3_5]), get_min_array([y3_1,y3_2,y3_3,y3_4,y3_5]),   lw=2,  color=c3  , alpha=0.25, label=r''   )
    
    
    
    plot.set_xlim( (5, 1e5)         )
    plot.set_ylim( (5e-36, 5e-31)   )
    
    plot.text    (1.2e1, 2e-33 , r'$\mathrm{T_{\bar{p}}= 0.5\,GeV\,\,(\times\, 40) }$',   color=c1, rotation=85 )
    plot.text    (0.5e2, 2e-34 , r'$\mathrm{T_{\bar{p}}=                    20\,GeV}$',   color=c2, rotation=80 )
    plot.text    (  3e2, 2e-34 , r'$\mathrm{T_{\bar{p}}=                   0.5\,TeV}$',   color=c3, rotation=75 )
    
    
    plt.savefig('XS_comparison_Hep.png')


    #
    #   He He
    #
    #    plot, fig = plot_1D (r'$\mathrm{T_{p} \quad [GeV]}$',      r'$\mathrm{d\sigma/dT_{\bar{p}}\quad[m^2/GeV]}$', 'log', 'log' )
    #    
    #    x1_1, y1_1 = profile( XS_parameterizations+'/dT_HeHe_pbar_LAB_diMauro12.txt',               0.5, type='prod', y_fac=40*2.3  )
    #    x1_2, y1_2 = profile( XS_parameterizations+'/dT_HeHe_pbar_LAB_WinklerWithHypWitNbar.txt',   0.5, type='prod', y_fac=40      )
    #    x1_3, y1_3 = profile( XS_parameterizations+'/dT_HeHe_pbar_LAB_KR_PPFRAG.txt',               0.5, type='prod', y_fac=40*2.3  )
    #    x1_4, y1_4 = profile( XS_parameterizations+'/dT_HeHe_pbar_LAB_Duperray.txt',                0.5, type='prod', y_fac=40*2.3  )
    #    x1_5, y1_5 = profile( XS_parameterizations+'/dT_HeHe_pbar_LAB_TanNg.txt',                   0.5, type='prod', y_fac=40*2.3  )
    #    
    #    x2_1, y2_1 = profile( XS_parameterizations+'/dT_HeHe_pbar_LAB_diMauro12.txt',               20., type='prod', y_fac=2.3 )
    #    x2_2, y2_2 = profile( XS_parameterizations+'/dT_HeHe_pbar_LAB_WinklerWithHypWitNbar.txt',   20., type='prod', y_fac=1.0 )
    #    x2_3, y2_3 = profile( XS_parameterizations+'/dT_HeHe_pbar_LAB_KR_PPFRAG.txt',               20., type='prod', y_fac=2.3 )
    #    x2_4, y2_4 = profile( XS_parameterizations+'/dT_HeHe_pbar_LAB_Duperray.txt',                20., type='prod', y_fac=2.3 )
    #    x2_5, y2_5 = profile( XS_parameterizations+'/dT_HeHe_pbar_LAB_TanNg.txt',                   20., type='prod', y_fac=2.3 )
    #    
    #    x3_1, y3_1 = profile( XS_parameterizations+'/dT_HeHe_pbar_LAB_diMauro12.txt',               500, type='prod', y_fac=2.3 )
    #    x3_2, y3_2 = profile( XS_parameterizations+'/dT_HeHe_pbar_LAB_WinklerWithHypWitNbar.txt',   500, type='prod', y_fac=1.0 )
    #    x3_3, y3_3 = profile( XS_parameterizations+'/dT_HeHe_pbar_LAB_KR_PPFRAG.txt',               500, type='prod', y_fac=2.3 )
    #    x3_4, y3_4 = profile( XS_parameterizations+'/dT_HeHe_pbar_LAB_Duperray.txt',                500, type='prod', y_fac=2.3 )
    #    x3_5, y3_5 = profile( XS_parameterizations+'/dT_HeHe_pbar_LAB_TanNg.txt',                   500, type='prod', y_fac=2.3 )
    #    
    #    
    #    plot.fill_between( x1_1, get_max_array([y1_1,y1_2,y1_3,y1_4,y1_5]), get_min_array([y1_1,y1_2,y1_3,y1_4,y1_5]),   lw=2,  color=c1  , alpha=0.25, label=r''   )
    #    plot.fill_between( x2_1, get_max_array([y2_1,y2_2,y2_3,y2_4,y2_5]), get_min_array([y2_1,y2_2,y2_3,y2_4,y2_5]),   lw=2,  color=c2  , alpha=0.25, label=r''   )
    #    plot.fill_between( x3_1, get_max_array([y3_1,y3_2,y3_3,y3_4,y3_5]), get_min_array([y3_1,y3_2,y3_3,y3_4,y3_5]),   lw=2,  color=c3  , alpha=0.25, label=r''   )
    #    
    #    plot.set_xlim( (5, 1e5)         )
    #    plot.set_ylim( (5e-36, 5e-31)   )
    #    
    #    plot.text    (1.2e1, 2e-33 , r'$\mathrm{T_{\bar{p}}= 0.5\,GeV\,\,(\times\, 40) }$',   color=c1, rotation=85 )
    #    plot.text    (0.5e2, 2e-34 , r'$\mathrm{T_{\bar{p}}=                    20\,GeV}$',   color=c2, rotation=80 )
    #    plot.text    (  7e2, 2e-34 , r'$\mathrm{T_{\bar{p}}=                   0.5\,TeV}$',   color=c3, rotation=75 )
    #    
    #    
    #    plt.savefig('XS_comparison_HeHe.png')






#! /usr/bin/env python

import glob
import math
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
from matplotlib import rc

from PlotFunctions import plot2D, plot2D_easy, profile




print_size=10
label_size=30


def plot_1D( xlabel, ylabel, xscale='log', yscale='log',  legend=''):

    global print_size, label_size
    plt.close('all')
    font_props = {"size":label_size}
    rc("font", **font_props)
    mpl.rcParams['axes.linewidth'] = 2
    mpl.rcParams['mathtext.fontset']='stixsans'

    fig     = plt.figure(figsize=(print_size*1.3, print_size))
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


import argparse

parser = argparse.ArgumentParser(description='CRACS - Plotting script for uncertainty calculator.')
parser.add_argument('--step',   help='Step .', action='store', dest='step', type=int, default=1)

args = parser.parse_args()
step   	    = args.step

# relative pbar uncertainty
if step==1:
    
    plot, fig = plot_1D (r'$\mathrm{T_{\bar{p}}\,[GeV]}$',      r'$\mathrm{\sigma_\phi/\phi}$', 'log', 'linear' )


    pbar_relative_uncertainty_parmetrization    =   np.genfromtxt(  'pbar_uncertainty/pbar_relative_uncertainty_parmetrization.txt',     skip_header=1   )
    plot.errorbar   (  pbar_relative_uncertainty_parmetrization[:,0], pbar_relative_uncertainty_parmetrization[:,1], color='black', alpha=0.5,                 lw=2 )

    plot.set_ylim   ( 0.,   0.3 )

    if len(glob.glob('pbar_uncertainty/pbar_relative_uncertainty_data.txt')):
        pbar_relative_uncertainty_data              =   np.genfromtxt(  'pbar_uncertainty/pbar_relative_uncertainty_data.txt',               skip_header=1   )
        plot.errorbar   (  pbar_relative_uncertainty_data[:,0],           pbar_relative_uncertainty_data[:,1],           color='black', fmt='o',   label='AMS-02', lw=0 )
        plot.legend     (loc='upper center', bbox_to_anchor=(0.5, 0.95), frameon=False, fontsize=0.8*label_size)

    plt.savefig     ( 'pbar_relative_uncertainty.png' )




# antiproton source term
if step==2:
    source_term_tot  =   np.genfromtxt(  'source_terms/source_term_tot.txt',   skip_header=1   )
    source_term_pp   =   np.genfromtxt(  'source_terms/source_term_pp.txt',    skip_header=1   )
    source_term_pHe  =   np.genfromtxt(  'source_terms/source_term_pHe.txt',   skip_header=1   )
    source_term_Hep  =   np.genfromtxt(  'source_terms/source_term_Hep.txt',   skip_header=1   )
    source_term_HeHe =   np.genfromtxt(  'source_terms/source_term_HeHe.txt',  skip_header=1   )

    plot, fig = plot_1D (r'$\mathrm{T_{\bar{p}}\,[GeV]}$',      r'$\mathrm{q \, T_\bar{p}^{2.7} \quad [GeV^{1.7}m^{-3}s^{-1}]}$', 'log', 'log' )
    
    plot.errorbar   (  source_term_tot[:,0],  source_term_tot [:,1]*np.power(source_term_tot [:,0], 2.7), color='black',    lw=4,  label='total'                    )
    plot.errorbar   (  source_term_pp [:,0],  source_term_pp  [:,1]*np.power(source_term_pp  [:,0], 2.7), color='black',    lw=2,  label='pp',                      )
    plot.errorbar   (  source_term_pHe[:,0],  source_term_pHe [:,1]*np.power(source_term_pHe [:,0], 2.7), color='black',    lw=2,  label='pHe',   dashes=(3,3)      )
    plot.errorbar   (  source_term_Hep[:,0],  source_term_Hep [:,1]*np.power(source_term_Hep [:,0], 2.7), color='black',    lw=2,  label='Hep',   dashes=(15,5)     )
    plot.errorbar   (  source_term_HeHe[:,0], source_term_HeHe[:,1]*np.power(source_term_HeHe[:,0], 2.7), color='black',    lw=2,  label='HeHe',  dashes=(20,5,5,5) )

    plot.set_ylim   ( 5e-38,   5e-33 )
    plot.legend     (loc='lower right', bbox_to_anchor=(0.95, 0.05), frameon=False, fontsize=0.8*label_size)
    plt.savefig     ( 'source_terms.png' )


    # relative uncertatinty

    source_term__relative_uncertainty__T_pbar__pp  =   np.genfromtxt(  'source_terms/relative_uncertainty__T_pbar__pp.txt',   skip_header=1   )
    plot, fig = plot_1D (r'$\mathrm{T_{\bar{p}}\,[GeV]}$',      r'$\mathrm{\sigma_q/q \, (pp)}$', 'log', 'linear' )
    plot.errorbar   (  source_term__relative_uncertainty__T_pbar__pp[:,0], source_term__relative_uncertainty__T_pbar__pp[:,1], color='black', alpha=0.5, lw=2, label='pp' )
    plot.set_ylim   ( 0.,   0.3 )
    plt.savefig     ( 'q_relative_uncertainty_pp.png' )





if step==3:
    # com contribution
    
    
    plot2D( '2D_contribution_LAB/pp/2D_contribution__pp.txt',  resfile='2D_contribution_LAB/pp_source_term.png', Tmax_prod=1e3, Tmax_proj=1e5, cmap_str='magma_r', Zlabel=r'$\mathrm{\frac{dq/d (log[T_p/GeV])}{q}} $', zscale='linear')
    plot2D( '2D_contribution_LAB/pp/2D_ccontribution__pp.txt', resfile='2D_contribution_LAB/pp_containence.png', Tmax_prod=1e3, Tmax_proj=1e5, cmap_str='magma',   Zlabel=r'$\mathrm{containence} $',                   zscale='linear')
    
    prof = profile('2D_contribution_LAB/pp/2D_contribution__pp.txt',  0.5,  resfile='2D_contribution_LAB/pp_contribution_profile_Tn_p.png', type='prod', dashes=(        ), Tn_min=5e0, Tn_max=1e5, y_min=1e-4, y_max=1,  label='0.5 GeV' )
    profile(       '2D_contribution_LAB/pp/2D_contribution__pp.txt',  10 ,  resfile='2D_contribution_LAB/pp_contribution_profile_Tn_p.png', type='prod', dashes=(3,3     ), Tn_min=5e0, Tn_max=1e5, y_min=1e-4, y_max=1,  label='10 GeV'  , draw_on_top=True )
    profile(       '2D_contribution_LAB/pp/2D_contribution__pp.txt', 100 ,  resfile='2D_contribution_LAB/pp_contribution_profile_Tn_p.png', type='prod', dashes=(15,5    ), Tn_min=5e0, Tn_max=1e5, y_min=1e-4, y_max=1,  label='100 GeV' , draw_on_top=True, y_label=r'$\mathrm{\frac{dq/d (log[T_p/GeV])}{q}} $', legend=True, loc='upper', x_label = r'$\mathrm{T_{p} \quad [GeV]}$')

    prof = profile('2D_contribution_LAB/pp/2D_ccontribution__pp.txt',  0.5, resfile='2D_contribution_LAB/pp_containence_profile_Tn_p.png', type='prod', dashes=(        ), Tn_min=5e0, Tn_max=1e5, y_min=0,    y_max=1, y_scale='linear',  label='0.5 GeV' )
    profile(       '2D_contribution_LAB/pp/2D_ccontribution__pp.txt',  10 , resfile='2D_contribution_LAB/pp_containence_profile_Tn_p.png', type='prod', dashes=(3,3     ), Tn_min=5e0, Tn_max=1e5, y_min=0,    y_max=1, y_scale='linear',  label='10 GeV'  , draw_on_top=True )
    profile(       '2D_contribution_LAB/pp/2D_ccontribution__pp.txt', 100 , resfile='2D_contribution_LAB/pp_containence_profile_Tn_p.png', type='prod', dashes=(15,5    ), Tn_min=5e0, Tn_max=1e5, y_min=0,    y_max=1, y_scale='linear',  label='100 GeV' , draw_on_top=True, y_label=r'$\mathrm{containence                 } $', legend=True, loc='upper', x_label = r'$\mathrm{T_{p} \quad [GeV]}$', )

    plot, fig = plot_1D (r'$\mathrm{T_{\bar{p}}\,[GeV]}$', r'$\mathrm{T_{p}\,[GeV]}$',      'log', 'log' )
    d_ccontribution = np.genfromtxt('2D_contribution_LAB/pp/2D_ccontribution__pp.txt')
    x = d_ccontribution[0,1:]
    y = d_ccontribution[1:,0]
    z = d_ccontribution[1:,1:]
    plot.contour(x, y, z, levels=[0.68,0.9, 0.95, 0.99])
    plt.savefig('2D_contribution_LAB/pp_contour.png')



# 3D Lab Proton
if step==4:

    data_Tbar=np.genfromtxt('save/save.fT_pbar.txt', skip_header=1);

    for i in range( len(data_Tbar) ):
        plot2D_easy('3D_Tp_Tpbar_eta_LAB/pp/2D_contribution_Tp_eta__'+str(i)+'__pp.txt',    x_label=r'$\mathrm{T_{p}\quad [GeV]}$', y_label=r'$\mathrm{\eta}$', x_min=1e0, x_max=1e5, y_min=-1, y_max=13, z_min=1e-7,    z_max=1e-3,    cmap_str='magma_r', resfile='3D_Tp_Tpbar_eta_LAB/pp/2D_contribution_Tp_eta__'+str(i)+'__pp.png',  z_label=r'$\mathrm{\frac{d^2q}{d (log[T_p/GeV])\,d\eta}}$',  x_scale='log', y_scale='linear', z_scale='log'   )
        plot2D_easy('3D_Tp_Tpbar_eta_LAB/pp/2D_ccontribution_Tp_eta__'+str(i)+'__pp.txt',   x_label=r'$\mathrm{T_{p}\quad [GeV]}$', y_label=r'$\mathrm{\eta}$', x_min=1e0, x_max=1e5, y_min=-1, y_max=13, z_min=0,       z_max=1,       cmap_str='magma',   resfile='3D_Tp_Tpbar_eta_LAB/pp/2D_ccontribution_Tp_eta__'+str(i)+'__pp.png', z_label='containence',                                       x_scale='log', y_scale='linear', z_scale='linear')

        plot, fig = plot_1D (r'$\mathrm{T_{p}\,[GeV]}$', r'$\mathrm{\eta}$',      'log', 'linear' )
        d_ccontribution = np.genfromtxt('3D_Tp_Tpbar_eta_LAB/pp/2D_ccontribution_Tp_eta__'+str(i)+'__pp.txt')
        y = d_ccontribution[0,1:]
        x = d_ccontribution[1:,0]
        z = np.transpose( d_ccontribution[1:,1:] )
        plot.contour(x, y, z, levels=[0.68,0.9,0.95,0.99])
        #plot.legend(loc='upper right', bbox_to_anchor=(0.95, 0.95), frameon=False, fontsize=0.8*label_size, title=r'$T_{pbar}='+str(data_Tbar[i])+'$')
        plot.text(10, 8.5,r'$\mathrm{T_{pbar}='+str("%.1f" % data_Tbar[i])+'\,GeV}$')
        plt.savefig('3D_Tp_Tpbar_eta_LAB/pp/2D_contour_Tp_eta__'+str(i)+'__pp.png')





# 3D CM Proton
if step==5:
    
    data_S=np.genfromtxt('save/save.fS.txt', skip_header=1);
    
    for i in range( len(data_S) ):
        plot2D_easy('3D_s_xR_pT/pp/2D_ccontribution_xR_pT__'+str(i)+'__pp.txt',   y_label=r'$\mathrm{p^T\quad [GeV]}$', x_label=r'$\mathrm{x_R}$', y_min=1e-1, y_max=1e2, x_min=0, x_max=1, z_min=0,       z_max=2,       cmap_str='magma',   resfile='3D_s_xR_pT/pp/2D_ccontribution_xR_pT__'+str(i)+'__pp.png', z_label='containence',   x_scale='linear', y_scale='log', z_scale='linear')

        plot2D_easy('3D_s_xR_pT/pp/2D_ccontribution_xR_pT__'+str(i)+'__pp_pos.txt',   y_label=r'$\mathrm{p^T\quad [GeV]}$', x_label=r'$\mathrm{x_R}$', y_min=1e-1, y_max=1e2, x_min=0, x_max=1, z_min=0,       z_max=2,       cmap_str='magma',   resfile='3D_s_xR_pT/pp/2D_ccontribution_xR_pT__'+str(i)+'__pp_pos.png', z_label='containence',   x_scale='linear', y_scale='log', z_scale='linear')
        plot2D_easy('3D_s_xR_pT/pp/2D_ccontribution_xR_pT__'+str(i)+'__pp_neg.txt',   y_label=r'$\mathrm{p^T\quad [GeV]}$', x_label=r'$\mathrm{x_R}$', y_min=1e-1, y_max=1e2, x_min=0, x_max=1, z_min=0,       z_max=2,       cmap_str='magma',   resfile='3D_s_xR_pT/pp/2D_ccontribution_xR_pT__'+str(i)+'__pp_neg.png', z_label='containence',   x_scale='linear', y_scale='log', z_scale='linear')

#        plot2D_easy('3D_Tp_Tpbar_eta_LAB/pp/2D_contribution_Tp_eta__'+str(i)+'__pp.txt',    x_label=r'$\mathrm{T_{p}\quad [GeV]}$', y_label=r'$\mathrm{\eta}$', x_min=1e0, x_max=1e5, y_min=-1, y_max=13, z_min=1e-7,    z_max=1e-3,    cmap_str='magma_r', resfile='3D_Tp_Tpbar_eta_LAB/pp/2D_contribution_Tp_eta__'+str(i)+'__pp.png',  z_label=r'$\mathrm{\frac{d^2q}{d (log[T_p/GeV])\,d\eta}}$',  x_scale='log', y_scale='linear', z_scale='log'   )
#        plot2D_easy('3D_Tp_Tpbar_eta_LAB/pp/2D_ccontribution_Tp_eta__'+str(i)+'__pp.txt',   x_label=r'$\mathrm{T_{p}\quad [GeV]}$', y_label=r'$\mathrm{\eta}$', x_min=1e0, x_max=1e5, y_min=-1, y_max=13, z_min=0,       z_max=1,       cmap_str='magma',   resfile='3D_Tp_Tpbar_eta_LAB/pp/2D_ccontribution_Tp_eta__'+str(i)+'__pp.png', z_label='containence',                                   x_scale='log', y_scale='linear', z_scale='linear')
#        
#        plot, fig = plot_1D (r'$\mathrm{T_{p}\,[GeV]}$', r'$\mathrm{\eta}$',      'log', 'linear' )
#        d_ccontribution = np.genfromtxt('3D_Tp_Tpbar_eta_LAB/pp/2D_ccontribution_Tp_eta__'+str(i)+'__pp.txt')
#        y = d_ccontribution[0,1:]
#        x = d_ccontribution[1:,0]
#        z = np.transpose( d_ccontribution[1:,1:] )
#        plot.contour(x, y, z, levels=[0.68,0.9,0.95,0.99])
#        plt.savefig('3D_Tp_Tpbar_eta_LAB/pp/2D_contour_Tp_eta__'+str(i)+'__pp.png')




# 3D Lab Proton
if step==6:
    
    data_Tbar=np.genfromtxt('save/save.fT_pbar.txt', skip_header=1);
    
    for i in range( len(data_Tbar) ):
        
        plot, fig = plot_1D (r'$\mathrm{T_{p}\,[GeV]}$', r'$\mathrm{\eta}$',      'log', 'linear' )
        d_ccontribution = np.genfromtxt('3D_Tp_Tpbar_eta_LAB/pp/2D_relativeUncertainty_Tp_eta__'+str(i)+'__pp.txt')
        y = d_ccontribution[0,1:]
        x = d_ccontribution[1:,0]
        z = np.transpose( d_ccontribution[1:,1:] )
        plot.contour(x, y, z, levels=[0.1])
        #plot.legend(loc='upper right', bbox_to_anchor=(0.95, 0.95), frameon=False, fontsize=0.8*label_size, title=r'$T_{pbar}='+str(data_Tbar[i])+'$')
        plot.text(10, 8.5,r'$\mathrm{T_{pbar}='+str("%.1f" % data_Tbar[i])+'\,GeV}$')
        plt.savefig('3D_Tp_Tpbar_eta_LAB/pp/2D_relativeUncertainty__Tp_eta__'+str(i)+'__pp.png')

if step==62:
    os.system('mkdir 3D_s_xR_pT/pp/uncertainty')
    data_S=np.genfromtxt('save/save.fS.txt', skip_header=1);
    for i in range( len(data_S) ):
        plot2D_easy('3D_s_xR_pT/pp/2D_relativeUncertainty_xR_pT__'+str(i)+'__pp.txt',   y_label=r'$\mathrm{p^T\quad [GeV]}$', x_label=r'$\mathrm{x_R}$', y_min=1e-1, y_max=1e2, x_min=0, x_max=1, z_min=0,       z_max=2,       cmap_str='magma',   resfile='3D_s_xR_pT/pp/2D_relativeUncertainty_xR_pT__'+str(i)+'__pp.png', z_label='uncertainty',   x_scale='linear', y_scale='log', z_scale='linear')

        plot, fig = plot_1D (r'$\mathrm{x_R}$', r'$\mathrm{p^T\quad [GeV]}$',      'linear', 'log' )
        d_ccontribution = np.genfromtxt('3D_s_xR_pT/pp/2D_relativeUncertainty_xR_pT__'+str(i)+'__pp.txt')
        y = d_ccontribution[0,1:]
        x = d_ccontribution[1:,0]
        z = np.transpose( d_ccontribution[1:,1:] )
        plot.contour(x, y, z, levels=[0.5])
        #plot.legend(loc='upper right', bbox_to_anchor=(0.95, 0.95), frameon=False, fontsize=0.8*label_size, title=r'$T_{pbar}='+str(data_Tbar[i])+'$')

        ss = np.sqrt(data_S[i])
        plot.text(0.1, 20 ,r'$\mathrm{\sqrt{s}='+str("%.1f" % ss )+'\, GeV}$')
        plt.savefig('3D_s_xR_pT/pp/uncertainty/2D_relativeUncertainty__xR_pT__'+str(i)+'__pp.png')




if step==7:
    
    CRACS = os.environ['CRACS']
    print CRACS
    
    CS_data_NA49        = np.genfromtxt( CRACS+'/data/CS_measurements/AnticicPbarNA492010conv_new.dat')
    CS_data_Allaby      = np.genfromtxt( CRACS+'/data/CS_measurements/AllabyPbarCERN1970conv_withsysnormstat.dat')
    CS_data_Antreasyan  = np.genfromtxt( CRACS+'/data/CS_measurements/AntreasyanPbarFNAL1979conv.dat')
    CS_data_Guetler     = np.genfromtxt( CRACS+'/data/CS_measurements/GuetlerPbarCERNISR1976conv.dat')
    CS_data_Capiluppi   = np.genfromtxt( CRACS+'/data/CS_measurements/CapiluppiPbarCERNISR1974conv.dat')
    CS_data_Johnson     = np.genfromtxt( CRACS+'/data/CS_measurements/JohnsonPbarFNAL1978conv.dat')
    
    marker = ['o','s','v','^']
    
    data        = [CS_data_NA49, CS_data_Allaby, CS_data_Antreasyan, CS_data_Guetler, CS_data_Capiluppi, CS_data_Johnson]
    experiment  = ['NA49',      'Allaby',       'Antreasyan',       'Guetler',       'Capiluppi',       'Johnson'      ]
    
                                 
                              
#
#                                          CapiluppiPbarCERNISR1974conv.dat           JohnsonPbarFNAL1978conv.dat
    
    data_S=np.genfromtxt('3D_s_xR_pT/pp_CS_data/sqrtS_values.txt', skip_header=1);
    for i in range( len(data_S) ):
        plot, fig = plot_1D (r'$\mathrm{x_R}$', r'$\mathrm{p^T\quad [GeV]}$',      'linear', 'log' )
        plt.subplots_adjust(left=0.15, right=0.75, top=0.9, bottom=0.15)
        d_ccontribution = np.genfromtxt('3D_s_xR_pT/pp_CS_data/2D_relativeUncertainty_xR_pT__'+str(i)+'__pp.txt')
        y = d_ccontribution[0,1:]
        x = d_ccontribution[1:,0]
        z = np.transpose( d_ccontribution[1:,1:] )
        plot.contour(x, y, z, levels=[0.5])
        #plot.legend(loc='upper right', bbox_to_anchor=(0.95, 0.95), frameon=False, fontsize=0.8*label_size, title=r'$T_{pbar}='+str(data_Tbar[i])+'$')
        
        xR=[]
        pT=[]
        err=[]
        
        sqrtS = data_S[i]
        
        cmap    = plt.get_cmap('OrRd')
        cmap.set_under('#3ADF00')
        cmap.set_over('black')
        
        norm = colors.Normalize(vmin=0.02, vmax=0.2) #LogNorm
        
        scatter = plot.scatter(xR, pT, c=err, cmap=cmap, s=100, norm=norm)
        
        
        m=0
        for k in range(len(data)):
            d = data[k]
            xR=[]
            pT=[]
            err=[]
            for j in range(len(d[:,0])):
                if np.fabs(sqrtS-d[j,0])/sqrtS<0.05:
                    pT. append(d[j,1])
                    xR. append(d[j,2])
                    err.append(d[j,4]/d[j,3])
            if len(xR):
                scatter = plot.scatter(xR, pT, c=err, cmap=cmap, s=100, norm=norm, label=experiment[k], marker=marker[m])
                m = m+1





        cbar_ax = fig.add_axes([0.8, 0.15, 0.02, 0.75])
        cbar    = fig.colorbar(scatter, cax=cbar_ax, cmap=cmap, extend='both')
        cbar.     ax.set_ylabel('relative uncertainty')

        plot.set_xlim(0,0.75)
        plot.set_ylim(2e-2,5)


        plot.legend(loc='upper right', bbox_to_anchor=(0.95, 0.95), frameon=False, fontsize=0.8*label_size)
    
        plot.text(0.1, 3 ,r'$\mathrm{\sqrt{s}='+str("%.1f" % sqrtS )+'\, GeV}$')
        plt.savefig('3D_s_xR_pT/pp_CS_data/2D_relativeUncertainty__xR_pT__'+str(i)+'__pp.png')






if step==8:
    data_Tpbar=np.genfromtxt('3D_Tp_Tpbar_eta_LAB/pp_fixed_Tp/Tp_values.txt', skip_header=1);
    for i in range( len(data_Tpbar) ):
        plot, fig = plot_1D (r'$\mathrm{T_{\bar{p}}\,[GeV]}$', r'$\mathrm{\eta}$',      'log', 'linear' )
        d_ccontribution = np.genfromtxt('3D_Tp_Tpbar_eta_LAB/pp_fixed_Tp/2D_relativeUncertainty_Tpbar_eta__'+str(i)+'__pp.txt')
        y = d_ccontribution[0,1:]
        x = d_ccontribution[1:,0]
        z = np.transpose( d_ccontribution[1:,1:] )
        plot.contour(x, y, z, levels=[0.5])
        #plot.legend(loc='upper right', bbox_to_anchor=(0.95, 0.95), frameon=False, fontsize=0.8*label_size, title=r'$T_{pbar}='+str(data_Tbar[i])+'$')
        
        ss = (data_Tpbar[i]/1e3)
        plot.text(1, 9 ,r'$\mathrm{T_{p}='+str("%.2f" % ss )+'\, TeV}$')
        plt.savefig('3D_Tp_Tpbar_eta_LAB/pp_fixed_Tp/2D_relativeUncertainty_Tpbar_eta__'+str(i)+'__pp.png')

















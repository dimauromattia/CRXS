#! /usr/bin/env python

import glob
import math
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
import matplotlib.patches as mpatches
from matplotlib import rc
from matplotlib.colors import LinearSegmentedColormap

from PlotFunctions import plot2D, plot2D_easy, profile


c0='black'
c1='#04B404'
c2='#0489B1'
c3='#0404B4'
c4='#B40404'
c5='#B40486'

cb='#EFEFFF'

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
    
    plot.set_xlabel ( xlabel, fontsize=label_size*1.4 )
    plot.set_ylabel ( ylabel, fontsize=label_size*1.4 )
    
    plot.set_xscale ( xscale )
    plot.set_yscale ( yscale )

    plot.tick_params('both', length=20, width=2, which='major', top=True, right=True, direction='in', pad=10)
    plot.tick_params('both', length=10, width=1, which='minor', top=True, right=True, direction='in', pad=10)

    plot.grid(b=True, which='major', alpha=0.1, linestyle='-', linewidth=2)
    plot.tick_params(axis='both', pad=10)

    return plot, fig

def smoothing(x, width=10, type='log'):
    if not (type=='linear' or type=='log'):
        print 'Wrong type in smoothing(): Choose either linear or log.'
        sys.exit(0)
    x_res=[]
    l    = len(x)
    for i in range(len(x)):
        s = 0
        for j in range(-width,width+1):
            if type=='log':
                s = s + np.log(x[(i+l+j)%l])
            if type=='linear':
                s = s + x[(i+l+j)%l]
        s = s/(2.*width+1)
        if type=='log':
            s = np.exp(s)
        x_res.append(s)
    return x_res

def smoothContour(plot, x, y, z, xscale='log', yscale='log', levels=[0,0.5], fc=(0.,0.,1.,0.2), ec='black', lw=2, dashes=(), alpha=1, smoothing_degree=3, label='' ):
    CS = plot.contour(x, y, z, levels=levels, linewidths=0 )
    if len(CS.collections[0].get_paths()) < 1:
        return
    
    for i in range(len(CS.collections[0].get_paths())):
    
        p = CS.collections[0].get_paths()[i]
        v = p.vertices
        xv = smoothing( v[:,0], smoothing_degree, type=xscale )
        yv = smoothing( v[:,1], smoothing_degree, type=yscale )

        #plot.plot(v[:,0],v[:,1], color='red', lw=lw, dashes=dashes)
        
        cont = mpatches.Polygon(np.transpose([xv,yv]), lw=0, fc=fc, alpha=alpha)
        plot.add_patch(cont)
        xv.append(xv[0])
        yv.append(yv[0])
        plot.plot(xv,yv, color=ec, lw=lw, dashes=dashes, label=label)

    #    CS = plot.contour(x, y, z, levels=[0.1], colors=color_list[j], linewidths=3)
    #    plot.plot([1,-1], [2,-1], color=color_list[j], dashes=dashes[j], label=label_list[j], lw=3)
    #    for c in CS.collections:
    #        c.set_dashes([(0, dashes[j])])


import argparse

parser = argparse.ArgumentParser(description='CRACS - Plotting script for uncertainty calculator.')
parser.add_argument('--step',     help='Step .',   action='store', dest='step', type=int, default=1)
parser.add_argument('--option',   help='Option .', action='store', dest='opt',  type=str, default='standard')

args = parser.parse_args()
step   	    = args.step
option      = args.opt

# relative pbar uncertainty
if step==1:
    
    plot, fig = plot_1D (r'$\mathrm{T_{\bar{p}}\,[GeV]}$',      r'$\mathrm{\sigma_\phi/\phi}$', 'log', 'linear', sizey=0.55 )

    plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.25)
    
    pbar_relative_uncertainty_parmetrization    =   np.genfromtxt(  option+'/pbar_uncertainty/pbar_source_uncertainty.txt',     skip_header=1   )
    plot.errorbar   (  pbar_relative_uncertainty_parmetrization[:,0], pbar_relative_uncertainty_parmetrization[:,1]*100, color='black', alpha=0.5,                 lw=2 )

    plot.set_xlim   ( 5e-1, 1e3 )
    plot.set_ylim   ( 0.,   30  )
    plot.set_yticks([0, 10, 20, 30])

    fmt = '%.0f%%'
    plot.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter(fmt))
    
    if len(glob.glob(option+'/pbar_uncertainty/pbar_source_uncertainty_data.txt')):
        pbar_relative_uncertainty_data              =   np.genfromtxt(  option+'/pbar_uncertainty/pbar_source_uncertainty_data.txt',               skip_header=1   )
        plot.errorbar   (  pbar_relative_uncertainty_data[:,0],           pbar_relative_uncertainty_data[:,1]*100,           color='black', fmt='o',   label='AMS-02', lw=0 )
        plot.legend     (  loc='upper center', numpoints=1, bbox_to_anchor=(0.5, 0.95), frameon=False, fontsize=0.8*label_size  )

    plt.savefig     ( option+'/pbar_flux_uncertainty.pdf' )



# 3D Lab Proton
if step==2:

    data_Tbar=np.genfromtxt('save/save.fT_pbar.txt', skip_header=1);
    
    os.system('mkdir '+option+'/contribution_LAB/contribution_plots')
    os.system('mkdir '+option+'/contribution_LAB/containence_plots')

    for i in range( len(data_Tbar) ):
        #plot2D_easy(option+'/contribution_LAB/contribution/2D_contribution_Tp_eta__'+str(i)+'.txt',    x_label=r'$\mathrm{T_{p}\quad [GeV]}$', y_label=r'$\mathrm{\eta}$', x_min=1e0, x_max=1e5, y_min=-1, y_max=13, z_min=1e-7,    z_max=1e-3,    cmap_str='magma_r', resfile=option+'/contribution_LAB/contribution_plots/2D_contribution_Tp_eta__'+str(i)+'.pdf',  z_label=r'$\mathrm{\frac{d^2q}{d (log[T_p/GeV])\,d\eta}}$',  x_scale='log', y_scale='linear', z_scale='log'   )
        
        plot, fig = plot_1D (r'$\mathrm{T_{p}\,[GeV]}$', r'$\mathrm{\eta}$',      'log', 'linear' )
        d_ccontribution = np.genfromtxt(option+'/contribution_LAB/containence//2D_containence__Tp_eta__'+str(i)+'.txt')
        y = d_ccontribution[0,1:]
        x = d_ccontribution[1:,0]
        z = np.transpose( d_ccontribution[1:,1:] )
        CS = plot.contour(x, y, z, levels=[0.9,0.99,0.999], colors='k', linewidths=2)
        fmt = {}
        strs = ['90%', '99%', '99.9%']
        for l, s in zip(CS.levels, strs):
            fmt[l] = s
        plt.clabel(CS, inline=1, fontsize=25, fmt=fmt)
        if data_Tbar[i]<1:
            plot.text(20, 10.7,r'$\mathrm{T_{pbar}='+str("%.2f" % data_Tbar[i])+'\,GeV}$')
        elif data_Tbar[i]<10:
            plot.text(20, 10.7,r'$\mathrm{T_{pbar}='+str("%.1f" % data_Tbar[i])+'\,GeV}$')
        else:
            plot.text(20, 10.7,r'$\mathrm{T_{pbar}='+str("%.0f" % data_Tbar[i])+'\,GeV}$')
        plt.savefig(option+'/contribution_LAB/containence_plots/2D_containence__Tp_eta__'+str(i)+'.pdf')

m_proton = 0.938
def xF_from_pT(pT, T_pbar_LAB, sqrtS):

    s = sqrtS*sqrtS*1.

    E_p_LAB       = (s - m_proton*m_proton)/2./m_proton
    beta          = np.sqrt(E_p_LAB - m_proton)/np.sqrt(E_p_LAB + m_proton)
    gamma         = 1./np.sqrt(1 - beta*beta)
    gammabeta     = gamma * beta



    E_pbar_LAB    = T_pbar_LAB+m_proton
    pL_pbar_LAB   = np.sqrt(  E_pbar_LAB*E_pbar_LAB - m_proton*m_proton - pT*pT  )

   #E_pbar        =  gamma     * E_pbar_LAB - gammabeta * pL_pbar_LAB
    pL_pbar       = -gammabeta * E_pbar_LAB + gamma     * pL_pbar_LAB

    x_F = 2.*pL_pbar/np.sqrt(s)

    return x_F


if step==22:
    dashes      = [(), (3,3), (6,3), (15,5), (20,5,5,5), (15,3,3,3,3,3)  ]

    #data_Tbar=np.genfromtxt('save/save.fT_pbar.txt', skip_header=1);
    data_S=np.genfromtxt('save/save.fS.txt', skip_header=1);
    
    os.system('mkdir '+option+'/contribution_CM/containence')
    
    for i in range( 0, len(data_S) ):
        #plot2D_easy(option+'/contribution_LAB/contribution/2D_contribution_Tp_eta__'+str(i)+'.txt',    x_label=r'$\mathrm{T_{p}\quad [GeV]}$', y_label=r'$\mathrm{\eta}$', x_min=1e0, x_max=1e5, y_min=-1, y_max=13, z_min=1e-7,    z_max=1e-3,    cmap_str='magma_r', resfile=option+'/contribution_LAB/contribution_plots/2D_contribution_Tp_eta__'+str(i)+'.pdf',  z_label=r'$\mathrm{\frac{d^2q}{d (log[T_p/GeV])\,d\eta}}$',  x_scale='log', y_scale='linear', z_scale='log'   )
        plot, fig = plot_1D (r'$\mathrm{x_F}$', r'$\mathrm{p_{T}\, [GeV]}$', 'linear', 'log' )
        d_ccontribution = np.genfromtxt(option+'/contribution_CM/containence/2D_contribution__xF_pT__'+str(i)+'.txt')
        y = d_ccontribution[0,1:]
        x = d_ccontribution[1:,0]
        z = np.transpose( d_ccontribution[1:,1:] )
        smoothContour(plot, x, y, z, 'linear', 'log', levels=[0.95], fc=(0.9,0.9,1.0,1.0), lw=1, smoothing_degree=3 )
        smoothContour(plot, x, y, z, 'linear', 'log', levels=[0.90], fc=(0.8,0.8,1.0,1.0), lw=1, smoothing_degree=3 )
        smoothContour(plot, x, y, z, 'linear', 'log', levels=[0.68], fc=(0.7,0.7,1.0,1.0), lw=1, smoothing_degree=3 )
        ss = np.sqrt(data_S[i])
        print ss
        if ss<1:
            plot.text(0.4, 4e-0 ,r'$\mathrm{\sqrt{s}='+str("%.2f" % ss )+'\, GeV}$', backgroundcolor=(1.,1.,1.,0.95) )
        elif ss<10:
            plot.text(0.4, 4e-0 ,r'$\mathrm{\sqrt{s}='+str("%.1f" % ss )+'\, GeV}$', backgroundcolor=(1.,1.,1.,0.95) )
        else:
            plot.text(0.4, 4e-0 ,r'$\mathrm{\sqrt{s}='+str("%.0f" % ss )+'\, GeV}$', backgroundcolor=(1.,1.,1.,0.95) )

        pT = np.arange(-2,1.01,0.1)
        pT = np.power(10,pT)
        xF = xF_from_pT(pT, 1, ss)
        for j in range(len(xF)):
            if xF[j] != xF[j]:
                xF[j] = -2.
        plot.plot(xF,pT, color=c1, dashes=dashes[1], lw=2, label=r'$\mathrm{T_{\bar{p}} = 1 \, GeV}$')

        pT = np.arange(-2,1.01,0.1)
        pT = np.power(10,pT)
        xF = xF_from_pT(pT, 10, ss)
        for j in range(len(xF)):
            if xF[j] != xF[j]:
                xF[j] = -2.
        plot.plot(xF,pT, color=c2, dashes=dashes[2], lw=2, label=r'$\mathrm{T_{\bar{p}} = 10 \, GeV}$')

        pT = np.arange(-2,1.01,0.1)
        pT = np.power(10,pT)
        xF = xF_from_pT(pT, 100, ss)
        for j in range(len(xF)):
            if xF[j] != xF[j]:
                xF[j] = -2.
        plot.plot(xF,pT, color=c3, dashes=dashes[3], lw=2, label=r'$\mathrm{T_{\bar{p}} = 100 \, GeV}$')

        plot.set_xlim( ( -1,1 ) )

        legend = plot.legend(loc='upper left', scatterpoints=1, bbox_to_anchor=(0.05, 0.95), ncol=1, fontsize=0.8*label_size)
        frame = legend.get_frame()
        frame.set_linewidth(0)
        frame.set_color( (1.,1.,1.,0.95) )

        plt.savefig(option+'/contribution_CM/containence/2D_contribution__xF_pT__'+str(i)+'.pdf')


if step==222:
    
    
    S_lhcb = 110
    S_int  = 0

    data_S=np.genfromtxt('save/save.fS.txt', skip_header=1)

    d = 1e90
    p = 0
    for i in range( len(data_S) ):
        dd = np.fabs(S_lhcb*S_lhcb - data_S[i])
        if dd < d:
            d       = dd
            S_int   = i


    print S_int

    dashes      = [(), (3,3), (6,3), (15,5), (20,5,5,5), (15,3,3,3,3,3)  ]
    os.system('mkdir '+option+'/contribution_CM/containence')

    i = S_int

    plot, fig = plot_1D (r'$\mathrm{x_F}$', r'$\mathrm{p_{T}\, [GeV]}$', 'linear', 'log' )
    d_ccontribution = np.genfromtxt(option+'/contribution_CM/containence/2D_contribution__xF_pT__'+str(i)+'.txt')
    y = d_ccontribution[0,1:]
    x = d_ccontribution[1:,0]
    z = np.transpose( d_ccontribution[1:,1:] )
    smoothContour(plot, x, y, z, 'linear', 'log', levels=[0.95], fc=(0.9,0.9,1.0,1.0), lw=1, smoothing_degree=3 )
    smoothContour(plot, x, y, z, 'linear', 'log', levels=[0.90], fc=(0.8,0.8,1.0,1.0), lw=1, smoothing_degree=3 )
    smoothContour(plot, x, y, z, 'linear', 'log', levels=[0.68], fc=(0.7,0.7,1.0,1.0), lw=1, smoothing_degree=3 )
    ss = np.sqrt(data_S[i])
    print ss
    if ss<1:
        plot.text(0.4, 4e-0 ,r'$\mathrm{\sqrt{s}='+str("%.2f" % ss )+'\, GeV}$', backgroundcolor=(1.,1.,1.,0.95) )
    elif ss<10:
        plot.text(0.4, 4e-0 ,r'$\mathrm{\sqrt{s}='+str("%.1f" % ss )+'\, GeV}$', backgroundcolor=(1.,1.,1.,0.95) )
    else:
        plot.text(-0.13, 0.33 ,r'$\mathrm{\sqrt{s}='+str("%.0f" % ss )+'\, GeV}$', backgroundcolor=(1.,1.,1.,0.95) )

    pT = np.arange(-2,1.01,0.1)
    pT = np.power(10,pT)
    xF = xF_from_pT(pT, 100, ss)
    for j in range(len(xF)):
        if xF[j] != xF[j]:
            xF[j] = -2.
    plot.plot(xF,pT, color=c1, dashes=dashes[1], lw=2, label=r'$\mathrm{T_{\bar{p}} = 100 \, GeV}$')

    pT = np.arange(-2,1.01,0.1)
    pT = np.power(10,pT)
    xF = xF_from_pT(pT, 300, ss)
    for j in range(len(xF)):
        if xF[j] != xF[j]:
            xF[j] = -2.
    plot.plot(xF,pT, color=c2, dashes=dashes[2], lw=2, label=r'$\mathrm{T_{\bar{p}} = 300 \, GeV}$')

    pT = np.arange(-2,1.01,0.1)
    pT = np.power(10,pT)
    xF = xF_from_pT(pT, 700, ss)
    for j in range(len(xF)):
        if xF[j] != xF[j]:
            xF[j] = -2.
    plot.plot(xF,pT, color=c3, dashes=dashes[3], lw=2, label=r'$\mathrm{T_{\bar{p}} = 700 \, GeV}$')

    plot.set_xlim( ( -0.15,0.15 ) )

    legend = plot.legend(loc='lower left', scatterpoints=1, bbox_to_anchor=(0.05, 0.05), ncol=1, fontsize=0.8*label_size)
    frame = legend.get_frame()
    frame.set_linewidth(0)
    frame.set_color( (1.,1.,1.,0.95) )

    plt.savefig(option+'/contribution_CM/containence/2D_contribution__xF_pT__'+str(i)+'_lhcb.pdf')



if step==23:
    
    
    data_Tbar=np.genfromtxt('save/save.fT_pbar.txt', skip_header=1);
    
    T_pbar = 50;
    T_pbar_int  = 0
    
    
    d = 1e90
    for i in range( len(data_Tbar) ):
        dd = np.fabs(T_pbar - data_Tbar[i])
        if dd < d:
            d = dd
            T_pbar_int = i


    i = T_pbar_int
    print data_Tbar[i]
    plot, fig = plot_1D (r'$\mathrm{T_{p}\,[GeV]}$', r'$\mathrm{\eta}$',      'log', 'linear' )
    d_ccontribution = np.genfromtxt(option+'/contribution_LAB/containence//2D_containence__Tp_eta__'+str(i)+'.txt')
    y = d_ccontribution[0,1:]
    x = d_ccontribution[1:,0]
    z = np.transpose( d_ccontribution[1:,1:] )
    CS = plot.contour(x, y, z, levels=[0.9,0.99,0.999], colors='k', linewidths=2)
    fmt = {}
    strs = ['90%', '99%', '99.9%']
    for l, s in zip(CS.levels, strs):
        fmt[l] = s
    plt.clabel(CS, inline=1, fontsize=25, fmt=fmt)
    if data_Tbar[i]<1:
        plot.text(20, 10.7,r'$\mathrm{T_{\bar{p}}='+str("%.2f" % data_Tbar[i])+'\,GeV}$')
    elif data_Tbar[i]<10:
        plot.text(20, 10.7,r'$\mathrm{T_{\bar{p}}='+str("%.1f" % data_Tbar[i])+'\,GeV}$')
    else:
        plot.text(20, 10.7,r'$\mathrm{T_{\bar{p}}='+str("%.0f" % data_Tbar[i])+'\,GeV}$')
    plt.savefig('2D_containence__Tp_eta__'+str(i)+'.pdf')

    option = 'Duperray'
    d_ccontribution = np.genfromtxt(option+'/contribution_LAB/containence//2D_containence__Tp_eta__'+str(i)+'.txt')
    y = d_ccontribution[0,1:]
    x = d_ccontribution[1:,0]
    z = np.transpose( d_ccontribution[1:,1:] )
    CS = plot.contour(x, y, z, levels=[0.9,0.99,0.999], colors=c1, linewidths=2)

    option = 'Winkler'
    d_ccontribution = np.genfromtxt(option+'/contribution_LAB/containence//2D_containence__Tp_eta__'+str(i)+'.txt')
    y = d_ccontribution[0,1:]
    x = d_ccontribution[1:,0]
    z = np.transpose( d_ccontribution[1:,1:] )
    CS = plot.contour(x, y, z, levels=[0.9,0.99,0.999], colors=c2, linewidths=2)

    plt.savefig('2D_containence__Tp_eta__'+str(i)+'_addXS.pdf')




# 3D Lab Proton
if step==3:
    
    data_Tbar=np.genfromtxt('save/save.fT_pbar.txt', skip_header=1)
    os.system('mkdir '+option+'/contribution_LAB/uncertainty_plots')
    
    for i in range( len(data_Tbar) ):
        
        plot, fig = plot_1D (r'$\mathrm{T_{p}\,[GeV]}$', r'$\mathrm{\eta}$',      'log', 'linear' )
        d_ccontribution = np.genfromtxt(option+'/contribution_LAB/uncertainty/2D_uncertainty__Tp_eta__'+str(i)+'.txt')
        y = d_ccontribution[0,1:]
        x = d_ccontribution[1:,0]
        z = np.transpose( d_ccontribution[1:,1:] )
        plot.contour(x, y, z, levels=[0.1], colors='k', linewidths=2)
        if data_Tbar[i]<1:
            plot.text(2.5, 10.7,r'$\mathrm{T_{pbar}='+str("%.2f" % data_Tbar[i])+'\,GeV}$')
        elif data_Tbar[i]<10:
            plot.text(2.5, 10.7,r'$\mathrm{T_{pbar}='+str("%.1f" % data_Tbar[i])+'\,GeV}$')
        else:
            plot.text(2.5, 10.7,r'$\mathrm{T_{pbar}='+str("%.0f" % data_Tbar[i])+'\,GeV}$')
        plt.savefig(option+'/contribution_LAB/uncertainty_plots/2D_uncertainty__Tp_eta__'+str(i)+'.pdf')



if step==32:

    data_S=np.genfromtxt('save/save.fS.txt', skip_header=1);
    os.system('mkdir '+option+'/contribution_CM/uncertainty_plots')

    for i in range( len(data_S) ):
        plot, fig = plot_1D (r'$\mathrm{x_R}$', r'$\mathrm{ p_{T}\quad [GeV]}$',      'log', 'log' )
        d_ccontribution = np.genfromtxt(option+'/contribution_CM/uncertainty/2D_uncertainty__xR_pT__'+str(i)+'.txt')
        y = d_ccontribution[0,1:]
        x = d_ccontribution[1:,0]
        z = np.transpose( d_ccontribution[1:,1:] )
        plot.contour(x, y, z, levels=[0.5], colors='k', linewidths=2)
        ss = np.sqrt(data_S[i])
        if ss<1:
            plot.text(0.1, 20 ,r'$\mathrm{\sqrt{s}='+str("%.2f" % ss )+'\, GeV}$')
        elif ss<10:
            plot.text(0.1, 20 ,r'$\mathrm{\sqrt{s}='+str("%.1f" % ss )+'\, GeV}$')
        else:
            plot.text(0.1, 20 ,r'$\mathrm{\sqrt{s}='+str("%.0f" % ss )+'\, GeV}$')
        plt.savefig(option+'/contribution_CM/uncertainty_plots/2D_uncertainty__xR_pT__'+str(i)+'.pdf')

if step==33:

    data_S=np.genfromtxt('save/save.fS.txt', skip_header=1);
    os.system('mkdir '+option+'/contribution_CM/uncertainty_plots')

    for i in range( len(data_S) ):
        plot, fig = plot_1D (r'$\mathrm{x_F}$', r'$\mathrm{ p_{T}\quad [GeV]}$',      'linear', 'log' )
        d_ccontribution = np.genfromtxt(option+'/contribution_CM/uncertainty/2D_uncertainty__xF_pT__'+str(i)+'.txt')
        y = d_ccontribution[0,1:]
        x = d_ccontribution[1:,0]
        z = np.transpose( d_ccontribution[1:,1:] )
        plot.contour(x, y, z, levels=[0.5], colors='k', linewidths=2)
        ss = np.sqrt(data_S[i])
        if ss<1:
            plot.text(0.1, 20 ,r'$\mathrm{\sqrt{s}='+str("%.2f" % ss )+'\, GeV}$')
        elif ss<10:
            plot.text(0.1, 20 ,r'$\mathrm{\sqrt{s}='+str("%.1f" % ss )+'\, GeV}$')
        else:
            plot.text(0.1, 20 ,r'$\mathrm{\sqrt{s}='+str("%.0f" % ss )+'\, GeV}$')
        plt.savefig(option+'/contribution_CM/uncertainty_plots/2D_uncertainty__xF_pT__'+str(i)+'.pdf')




if step==4:
    my_cdict= {'red':  ((0.0, 1.0, 1.0),
                        (0.4, 1.0, 1.0),
                        (1.0, 0.3, 0.3)),
    
              'green': ((0.0, 0.8, 0.8),
                        (0.4, 0.0, 0.0),
                        (1.0, 0.0, 0.0)),
        
              'blue':  ((0.0, 0.6, 0.6),
                        (0.4, 0.0, 0.0),
                        (1.0, 0.0, 0.0))
            }

    myMap = LinearSegmentedColormap('MyMap', my_cdict)
    plt.register_cmap(cmap=myMap)

    CRACS = os.environ['CRACS']
    print CRACS
    
    CS_data_NA49        = np.genfromtxt( CRACS+'/data/CS_measurements/pp/converted_na49.txt',        )
    CS_data_Allaby      = np.genfromtxt( CRACS+'/data/CS_measurements/pp/converted_allaby.txt',      )
    CS_data_Antreasyan  = np.genfromtxt( CRACS+'/data/CS_measurements/pp/converted_antreasyan.txt',  )
    CS_data_Guetler     = np.genfromtxt( CRACS+'/data/CS_measurements/pp/converted_guettler.txt',    )
    CS_data_Capiluppi   = np.genfromtxt( CRACS+'/data/CS_measurements/pp/converted_capiluppi.txt',   )
    CS_data_Johnson     = np.genfromtxt( CRACS+'/data/CS_measurements/pp/converted_johnson.txt',     )
    CS_data_Brahms      = np.genfromtxt( CRACS+'/data/CS_measurements/pp/converted_brahms.txt',      )
    CS_data_phenix      = np.genfromtxt( CRACS+'/data/CS_measurements/pp/converted_phenix.txt',      )
    CS_data_dekkers     = np.genfromtxt( CRACS+'/data/CS_measurements/pp/converted_dekkers.txt' ,    )


    marker = ['o','s','v','^']
    
    data        = [CS_data_NA49, CS_data_Allaby, CS_data_Antreasyan, CS_data_Guetler, CS_data_Capiluppi, CS_data_Johnson, CS_data_Brahms, CS_data_phenix, CS_data_dekkers]
    experiment  = ['NA49',      'Allaby',       'Antreasyan',       'Guettler',       'Capiluppi',       'Johnson',        'BRAHMS'       , 'Phenix', 'Dekkers'      ]
    
    
    data_S=np.genfromtxt(option+'/contribution_CM/pp_CS_data/sqrtS_values.txt', skip_header=1);
    for i in range( len(data_S) ):
        plot, fig = plot_1D (r'$\mathrm{x_R}$', r'$\mathrm{ p_{T}\quad [GeV]}$',      'log', 'log' )
        plt.subplots_adjust(left=0.15, right=0.75, top=0.9, bottom=0.15)
        d_ccontribution = np.genfromtxt(option+'/contribution_CM/pp_CS_data/2D_uncertainty__xR_pT__'+str(i)+'.txt')
        y = d_ccontribution[0,1:-2]
        x = d_ccontribution[1:-2,0]
        z = np.transpose( d_ccontribution[1:-2,1:-2] )
        
        smoothContour(plot, x, y, z, 'log', 'log', [0.,0.5], fc='#EFEFFF')
        
        #plot.contourf(x, y, z, levels=[-0.1, 0.5], colors=['#EFEFFF'], alpha=0.8,  zorder=0)
                                    
        xR=[]
        pT=[]
        err=[]
        
        sqrtS = data_S[i]
        
        cmap    = plt.get_cmap('MyMap')
        cmap.set_under('white') #3ADF00
        cmap.set_over('black')
        
        norm = colors.Normalize(vmin=3, vmax=30)
        
        scatter = plot.scatter(xR, pT, c=err, cmap=cmap, s=100, norm=norm)
        
        
        m=0
        for k in range(len(data)):
            d = data[k]
            xR=[]
            pT=[]
            err=[]
            
            xR_fine=[]
            pT_fine=[]
            for j in range(len(d[:,0])):
                if np.fabs(sqrtS-d[j,0])/sqrtS<0.05:
                    pT. append(d[j,1])
                    xR. append(d[j,2])
                    #err.append( np.sqrt( np.power(d[j,6],2) + np.power(d[j,7],2) )/d[j,5]*100)
                    err.append( np.sqrt( np.power(d[j,6],2) )/d[j,5]*100)
                    if err[-1]<3.:
                        pT_fine. append(d[j,1])
                        xR_fine. append(d[j,2])
        
            if len(xR):
                scatter = plot.scatter([1e-3], [1e-3], facecolor=cmap(0.7), s=100, norm=norm, label=experiment[k], marker=marker[m],  zorder=0)
                scatter = plot.scatter(xR, pT, c=err, cmap=cmap, facecolor='red', s=100, norm=norm, marker=marker[m],  zorder=10)
                if len(xR_fine):
                    scatter2 = plot.scatter(xR_fine, pT_fine, color='white', edgecolors='#04B404', linewidths=2, s=100, marker=marker[m],  zorder=11)
                m = m+1

        cbar_ax = fig.add_axes([0.8, 0.15, 0.02, 0.75])
        cbar    = fig.colorbar(scatter, cax=cbar_ax, cmap=cmap, extend='both', format='%.0f%%')
        cbar.     ax.set_ylabel('Statistical uncertainty')

        plot.set_xlim(2e-2,1)
        plot.set_ylim(2e-2,10)

        if sqrtS==63.:
            plot.legend(loc='lower right', scatterpoints=1, bbox_to_anchor=(0.95, 0.05), frameon=False, fontsize=0.8*label_size)
        else:
            plot.legend(loc='upper left', scatterpoints=1,  bbox_to_anchor=(0.05, 0.95), frameon=False, fontsize=0.8*label_size)
    
        plot.text(0.025, 3e-2 ,r'$\mathrm{\sqrt{s}='+str("%.1f" % sqrtS )+'\, GeV}$')
        plt.savefig(option+'/contribution_CM/pp_CS_data/2D_uncertainty__xR_pT__'+str(i)+'.pdf')


if step==5:
    data_Tp=np.genfromtxt(option+'/contribution_LAB/uncertainty_fixedTp/Tp_values.txt', skip_header=1);

    cmap    = plt.get_cmap('Blues')
    handles=[]

    plot, fig = plot_1D (r'$\mathrm{T_{\bar{p}}\,[GeV]}$', r'$\mathrm{\eta}$',      'log', 'linear' )
    for i in range( len(data_Tp) ):
        d_ccontribution = np.genfromtxt(option+'/contribution_LAB/uncertainty_fixedTp/2D_uncertainty_Tpbar_eta__'+str(i)+'.txt')
        y = d_ccontribution[0,1:-2]
        x = d_ccontribution[1:-2,0]
        z = np.transpose( d_ccontribution[1:-2,1:-2] )
        c = cmap(0.1+i*0.8/len(data_Tp))
        smoothContour(plot, x, y, z, 'log', 'linear', [0.,0.5], fc=c)
        ss = (data_Tp[i]/1e3)
        handles.append( mpatches.Patch(color=c, label=r'$\mathrm{T_{{p}}='+str("%.2f" % ss)+'\,TeV}$') )
        #plot.text(T_p_x[i], T_p_y[i], r'$\mathrm{T_{{p}}='+str("%.2f" % ss)+'\,GeV}$')

    lhcb = mpatches.Rectangle( (10, 2.5), 90,  2.5, label='LHCb', linewidth=2, linestyle='dashed', ec='black', fc=(1,0,0,0.2) )
    plot.add_patch(  lhcb  )
    handles.append(  lhcb  )
    plt.legend(loc='upper left', bbox_to_anchor=(0.05, 0.95), frameon=False, fontsize=0.8*label_size, handles=handles)
    plot.set_xlim(5e-1,5e2)
    plot.set_ylim(2,9)

#    n=20
#    mp = 0.938
#
#    lhcb_1_p  = np.ones(n+1)*np.sqrt(10*10+2*10*mp)
#    lhcb_1_pT = np.arange(np.log(0.1), np.log(5)+np.log(5/0.1)/40., np.log(5/0.1)/n)
#    lhcb_1_pT = np.exp(lhcb_1_pT)
#    lhcb_1_T   = np.sqrt( lhcb_1_p*lhcb_1_p+mp*mp)-mp
#    lhcb_1_eta = np.arccosh(lhcb_1_p/lhcb_1_pT)
#
#    lhcb_2_pT  = np.ones(n+1)*5
#    lhcb_2_p   = np.arange(np.log(10), np.log(100)+np.log(100/10)/40., np.log(100/10)/n)
#    lhcb_2_p   = np.exp(lhcb_2_p)
#    lhcb_2_T   = np.sqrt( lhcb_2_p*lhcb_2_p+mp*mp)-mp
#    lhcb_2_eta = np.arccosh(lhcb_2_p/lhcb_2_pT)
#
#
#    lhcb_3_p  = np.ones(n+1)*np.sqrt(100*100+2*100*mp)
#    lhcb_3_pT = np.arange(np.log(0.1), np.log(5)+np.log(5/0.1)/40., np.log(5/0.1)/n)
#    lhcb_3_pT = np.exp(lhcb_3_pT)
#    lhcb_3_T   = np.sqrt( lhcb_3_p*lhcb_3_p+mp*mp)-mp
#    lhcb_3_eta = np.arccosh(lhcb_3_p/lhcb_3_pT)
#
#
#    lhcb_4_pT  = np.ones(n+1)*0.1
#    lhcb_4_p   = np.arange(np.log(10), np.log(100)+np.log(100/10)/40., np.log(100/10)/n)
#    lhcb_4_p   = np.exp(lhcb_4_p)
#    lhcb_4_T   = np.sqrt( lhcb_4_p*lhcb_4_p+mp*mp)-mp
#    lhcb_4_eta = np.arccosh(lhcb_4_p/lhcb_4_pT)
#
#
#    plot.plot(lhcb_1_T, lhcb_1_eta, color='green')
#    plot.plot(lhcb_2_T, lhcb_2_eta, color='green')
#    plot.plot(lhcb_3_T, lhcb_3_eta, color='green')
#    plot.plot(lhcb_4_T, lhcb_4_eta, color='green')

    plt.savefig(option+'/contribution_LAB/uncertainty_fixedTp/2D_uncertainty__Tpbar_eta.pdf')

if step==52:
    data_Tp=np.genfromtxt(option+'/contribution_LAB/uncertainty_fixedTp/Tp_values.txt', skip_header=1);
    
    cmap    = plt.get_cmap('Blues')
    handles=[]
    
    plot, fig = plot_1D (r'$\mathrm{p_{\bar{p}}\,[GeV]}$', r'$\mathrm{ p_{T}\, [GeV]}$',      'log', 'log' )
    for i in range( len(data_Tp) ):
        d_ccontribution = np.genfromtxt(option+'/contribution_LAB/uncertainty_fixedTp/2D_uncertainty_p_pbar_p_T__'+str(i)+'.txt')
        y = d_ccontribution[0,1:-2]
        x = d_ccontribution[1:-2,0]
        z = np.transpose( d_ccontribution[1:-2,1:-2] )
        c = cmap(0.1+i*0.8/len(data_Tp))
        smoothContour(plot, x, y, z, 'log', 'log', [0.,0.5], fc=c)
        ss = (data_Tp[i]/1e3)
        if ss>=0.1:
            handles.append( mpatches.Patch(color=c, label=r'$\mathrm{T_{{p}}='+str("%.2f" % ss)+'\,TeV}$') )
        else:
            ss = ss*1000
            handles.append( mpatches.Patch(color=c, label=r'$\mathrm{T_{{p}}='+str("%.0f" % ss)+'\,GeV}$') )
    #plot.text(T_p_x[i], T_p_y[i], r'$\mathrm{T_{{p}}='+str("%.2f" % ss)+'\,GeV}$')

#    lhcb = mpatches.Rectangle( (10, 2.5), 90,  2.5, label='LHCb', linewidth=2, linestyle='dashed', ec='black', fc=(1,0,0,0.2) )
#    plot.add_patch(  lhcb  )
#    handles.append(  lhcb  )
    plt.legend(loc='upper left', bbox_to_anchor=(0.05, 0.95), frameon=False, fontsize=0.7*label_size, handles=handles, ncol=3)
    plot.set_xlim(5e-1,5e2)
    plot.set_ylim(2e-2,1e1)


    
    plt.savefig(option+'/contribution_LAB/uncertainty_fixedTp/2D_uncertainty__p_pbar_p_T.pdf')





if step==6 or step==62 or step==63:
    
    T_pbar_list = [1.1, 2, 6, 15, 40, 100, 300]
    T_pbar_int  = []

    T_pbar_x  = [10, 2e2, 7e2,  3e3,  9e3,  7e3, 4e3  ]
    T_pbar_y  = [1,  1.5, 2.6,  3.5,  4.6,  5.8, 7    ]
    
    if step==62:
        T_pbar_list = [50]
        T_pbar_int  = []
        
        T_pbar_x  = [6e3]
        T_pbar_y  = [4]

    if step==63:
        T_pbar_list =[0.1, 0.2, 0.5, 1.0]

        T_pbar_x  =[800,  100*2*2,  100*2,  100 ]
        T_pbar_y  =[1.4, 2.2,  3.0,    3.8 ]

    if option=='GAPS':
        T_pbar_list =[0.22, 0.3, 1.0, 1.4]

        T_pbar_x  =[10,  3e2,  2e2,  8e1 ]
        T_pbar_y  =[0.5, 1.4,  3,    3.8 ]


    data_Tbar=np.genfromtxt('save/save.fT_pbar.txt', skip_header=1)

    cmap    = plt.get_cmap('Blues')

    for t in T_pbar_list:
        d = 1e90
        p = 0
        for i in range( len(data_Tbar) ):
            dd = np.fabs(t - data_Tbar[i])
            if dd < d:
                d = dd
                p = i
        T_pbar_int.append(p)

    print T_pbar_int

    plot, fig = plot_1D (r'$\mathrm{T_{p}\,[GeV]}$', r'$\mathrm{\eta}$',      'log', 'linear' )
    for j in range(len(T_pbar_int)):
        i = T_pbar_int[j]
        d_ccontribution = np.genfromtxt(option+'/contribution_LAB/uncertainty/2D_uncertainty__Tp_eta__'+str(i)+'.txt')
        y = d_ccontribution[0,1:-2]
        x = d_ccontribution[1:-2,0]
        z = np.transpose( d_ccontribution[1:-2,1:-2] )
        c = cmap(0.2+j*0.7/len(T_pbar_int))
        smoothContour(plot, x, y, z, 'log', 'linear', [0.,0.5], fc=c)
        if T_pbar_list[j]<1:
            plot.text(T_pbar_x[j], T_pbar_y[j], r'$\mathrm{T_{\bar{p}}='+str("%.2f" % T_pbar_list[j])+'\,GeV}$')
        elif T_pbar_list[j]<10:
            plot.text(T_pbar_x[j], T_pbar_y[j], r'$\mathrm{T_{\bar{p}}='+str("%.1f" % T_pbar_list[j])+'\,GeV}$')
        else:
            plot.text(T_pbar_x[j], T_pbar_y[j], r'$\mathrm{T_{\bar{p}}='+str("%.0f" % T_pbar_list[j])+'\,GeV}$')

    plot.set_xlim(5e0,8e4)
    plot.set_ylim(0,10)
    if option=='GAPS':
        plot.set_xlim(5e0,5e3)
        plot.set_ylim(0,6)


    if step==62:
        plt.savefig(option+'/contribution_LAB/2D_uncertainty__Tp_eta_50.pdf')
    elif step==63 :
        plt.savefig(option+'/contribution_LAB/2D_uncertainty__Tp_eta_lowEnergy.pdf')
    else:
        plt.savefig(option+'/contribution_LAB/2D_uncertainty__Tp_eta.pdf')


    S_list =[5, 7, 12, 25, 45, 70, 87, 95, 110]
    S_int  =[]

    S_x  =[7e-1, 7e-1, 4e-1, 16e-2, 8e-2, 5e-2, 3.4e-2, 2.5e-2, 1.3e-2]
    S_y  =[0.5,  5,   10,    12,    10,   8,    7,      6,    1]

    if option=='GAPS':
        S_list =[4.6, 5.5, 8, 15, 20]

        S_x  =[8e-1, 7.5e-1, 5e-1, 3.5e-1, 2e-1]
        S_y  =[0.5,  4,      6,    4,    1,  ]


    data_S=np.genfromtxt('save/save.fS.txt', skip_header=1)

    for t in S_list:
        d = 1e90
        p = 0
        for i in range( len(data_S) ):
            dd = np.fabs(t*t - data_S[i])
            if dd < d:
                d = dd
                p = i
        S_int.append(p)

    print S_int

    plot, fig = plot_1D (r'$\mathrm{x_R}$', r'$\mathrm{ p_{T}\quad [GeV]}$',      'log', 'log' )

    for j in range(len(S_int)):
        i = S_int[j]
        d_ccontribution = np.genfromtxt(option+'/contribution_CM/uncertainty/2D_uncertainty__xR_pT__'+str(i)+'.txt')
        y = d_ccontribution[0,1:-2]
        x = d_ccontribution[1:-2,0]
        z = np.transpose( d_ccontribution[1:-2,1:-2] )
        c = cmap(0.1+j*0.8/len(S_int))
        smoothContour(plot, x, y, z, 'log', 'log', [0.,0.5], fc=c)
        ss = (S_list[j])
        if ss<1:
            plot.text(S_x[j], S_y[j], r'$\mathrm{\sqrt{s}='+str("%.2f" % ss)+'\,GeV}$', rotation=75)
        elif ss<10:
            plot.text(S_x[j], S_y[j], r'$\mathrm{\sqrt{s}='+str("%.1f" % ss)+'\,GeV}$', rotation=75)
        else:
            plot.text(S_x[j], S_y[j], r'$\mathrm{\sqrt{s}='+str("%.0f" % ss)+'\,GeV}$', rotation=75)

    plot.set_xlim(1e-2,2e0)
    plot.set_ylim(2e-2,2e1)
    if option=='GAPS':
        plot.set_xlim(1e-1,2e0)
        plot.set_ylim(2e-2,2e1)

    plt.savefig(option+'/contribution_CM/2D_uncertainty__xR_pT.pdf')

#
#
if step==66:
    
    for k in range(2):
        
        c_len = 4 + 6
        if k==0:
            S_list =[10., 7., 5.5, 4.9]
            S_x    =[0.48, 0.3 , 0.15, 0.03]
            S_y    =[0.19, 0.19, 0.19 , 0.19]
            S_r    =[45,  40,   40,   45]
            c_off  = 6
            cmap    = plt.get_cmap('Blues_r')
       
        if k==1:
            S_list =[15,   25.,  45.,  70.,   87.,  110.]
            S_x    =[0.37, 0.25, 0.08,  0.02, 0.05, -0.2]
            S_y    =[0.12, 0.17, 0.21,  0.3,  3.0,  2.00]
            S_r    =[50,   45,   55,    80,   0.,   0  ]
            c_off  = 4
            cmap    = plt.get_cmap('Blues')
        

        S_int  =[]


        data_S=np.genfromtxt('save/save.fS.txt', skip_header=1)

        for t in S_list:
            d = 1e90
            p = 0
            for i in range( len(data_S) ):
                dd = np.fabs(t*t - data_S[i])
                if dd < d:
                    d = dd
                    p = i
            S_int.append(p)

        print S_int

        plot, fig = plot_1D (r'$\mathrm{x_f}$', r'$\mathrm{ p_{T}\quad [GeV]}$',      'linear', 'log' )

        for j in range(len(S_int)):
            i = S_int[j]
            d_ccontribution = np.genfromtxt(option+'/contribution_CM/uncertainty/2D_uncertainty__xF_pT__'+str(i)+'.txt')
            y = d_ccontribution[0,1:-2]
            x = d_ccontribution[1:-2,0]
            z = np.transpose( d_ccontribution[1:-2,1:-2] )
            c = cmap(0.1+(j+c_off)*0.8/c_len)
            smoothContour(plot, x, y, z, 'linear', 'log', [0.,0.5], fc=c)
            ss = (S_list[j])
            if ss<1:
                plot.text(S_x[j], S_y[j], r'$\mathrm{\sqrt{s}='+str("%.2f" % ss)+'\,GeV}$', rotation=S_r[j], size=label_size*0.5, ha='left', va='bottom')
            elif ss<10:
                plot.text(S_x[j], S_y[j], r'$\mathrm{\sqrt{s}='+str("%.1f" % ss)+'\,GeV}$', rotation=S_r[j], size=label_size*0.5, ha='left', va='bottom')
            else:
                plot.text(S_x[j], S_y[j], r'$\mathrm{\sqrt{s}='+str("%.0f" % ss)+'\,GeV}$', rotation=S_r[j], size=label_size*0.5, ha='left', va='bottom')

        if k==1:
            plot.set_xlim(-0.4,0.8)
            plot.arrow(-0.1, 1.8,  0.1 -0.01, -1.8+0.6)
            plot.arrow( 0.1, 2.7, -0.1+0.01, -2.7+0.8)
        if k==0:
            plot.set_xlim(-0.4,0.8)
        plot.set_ylim(2e-2,2e1)

        plt.savefig(option+'/contribution_CM/2D_uncertainty__xF_pT__'+str(k)+'.pdf')



if step==7 or (step>70 and step<80) or step==8 or step==9:

    color_list  =[ c0, c1, c5, c2, c3, c4 ]
    dashes      = [(), (3,3), (6,3), (15,5), (20,5,5,5), (15,3,3,3,3,3)  ]

    if step==7:
        T_pbar      = 50
        option_list =['standard', 'Winkler', 'Duperray', '2_30', '3_100', 'pHe'] # ,'Duperray', '2_100'
        label_list  =['standard', 'Winkler', 'Duperray', r'2%$\,$/$\,$30%', r'3%$\,$/$\,$100%', r'p$\,$He'] # ,'Duperray', '2_100'
    if step==71:
        T_pbar      = 50
        option_list =['standard', 'Winkler', 'Duperray', 'pHe'] # ,'Duperray', '2_100'
        label_list  =['standard', 'Winkler', 'Duperray', r'p$\,$He'] # ,'Duperray', '2_100'
        color_list  =[ c0, c1, c5, c4 ]
        dashes      = [(), (3,3), (6,3), (15,3,3,3,3,3)  ]
    if step==72:
        T_pbar      = 50
        option_list =['standard', '2_30', '3_100'] # ,'Duperray', '2_100'
        label_list  =[r'3%$\,$/$\,$30% (standard)', r'2%$\,$/$\,$30%', r'3%$\,$/$\,$100%'] # ,'Duperray', '2_100'
        color_list  =[ c0, c2, c3 ]
        dashes      = [(), (15,5), (20,5,5,5) ]
    if step==73:
        T_pbar      = 0.1
        option_list =['0_05', '0_05__W', '0_05__D'] # ,'Duperray', '2_100'
        label_list  =['5% di Mauro', '5% Winkler', '5% Duperray'] # ,'Duperray', '2_100'
        color_list  =[ c0, c1, c5, c4 ]
        dashes      = [(), (3,3), (6,3), (15,3,3,3,3,3)  ]
    if step==74:
        T_pbar      = 0.2
        option_list =['0_05', '0_05__W', '0_05__D'] # ,'Duperray', '2_100'
        label_list  =['5% di Mauro', '5% Winkler', '5% Duperray'] # ,'Duperray', '2_100'
        color_list  =[ c0, c1, c5, c4 ]
        dashes      = [(), (3,3), (6,3), (15,3,3,3,3,3)  ]
    if step==75:
        T_pbar      = 0.5
        option_list =['0_05', '0_05__W', '0_05__D'] # ,'Duperray', '2_100'
        label_list  =['5% di Mauro', '5% Winkler', '5% Duperray'] # ,'Duperray', '2_100'
        color_list  =[ c0, c1, c5, c4 ]
        dashes      = [(), (3,3), (6,3), (15,3,3,3,3,3)  ]
    if step==76:
        T_pbar      = 1.0
        option_list =['0_05', '0_05__W', '0_05__D'] # ,'Duperray', '2_100'
        label_list  =['5% di Mauro', '5% Winkler', '5% Duperray'] # ,'Duperray', '2_100'
        color_list  =[ c0, c1, c5, c4 ]
        dashes      = [(), (3,3), (6,3), (15,3,3,3,3,3)  ]
    if step==8:
        T_pbar      = 350
        color_list  =[ c0, c1, c2, c3, c4 ]
        dashes      = [(), (3,3), (15,5), (20,5,5,5), (15,3,3,3,3,3)  ]
        option_list =['standard', '2024', '0_05'] # ,'Duperray', '2_100'
        label_list  =['standard', 'AMS-02 2024', '5%'] # ,'Duperray', '2_100'
    if step==9:
        T_pbar      = 50
        option_list =['pHe', 'pp_fixed', 'pp_fixed_6_60'] # ,'Duperray', '2_100'
        label_list  =[r'standard (pHe)', r'pp fixed', r'pp fixed, 6%$\,$/$\,$60%']


    T_pbar_int  = 0
    
    data_Tbar=np.genfromtxt('save/save.fT_pbar.txt', skip_header=1)
    
    d = 1e90
    for i in range( len(data_Tbar) ):
        dd = np.fabs(T_pbar - data_Tbar[i])
        if dd < d:
            d = dd
            T_pbar_int = i

    print T_pbar_int

    plot, fig = plot_1D (r'$\mathrm{T_{p}\,[GeV]}$', r'$\mathrm{\eta}$',      'log', 'linear' )
    for j in range(len(option_list)):
        option = option_list[j]
        i = T_pbar_int
        d_ccontribution = np.genfromtxt(option+'/contribution_LAB/uncertainty/2D_uncertainty__Tp_eta__'+str(i)+'.txt')
        y = d_ccontribution[0,1:-2]
        x = d_ccontribution[1:-2,0]
        z = np.transpose( d_ccontribution[1:-2,1:-2] )
        smoothContour(plot, x, y, z, 'log', 'linear', [0.,0.5], fc=cb, ec=color_list[j], dashes=dashes[j], label=label_list[j])

    plot.set_xlim(2e1,2e5)
    plot.set_ylim(3,10)

    if step>=73 and step<=76:
        plot.set_xlim(2e0,2e4)
        plot.set_ylim(0,6)
        j = T_pbar_int
        if data_Tbar[j]<1:
            plot.text(10, 4.5, r'$\mathrm{T_{\bar{p}}='+str("%.2f" % data_Tbar[j])+'\,GeV}$')
        elif data_Tbar[j]<10:
            plot.text(10, 4.5, r'$\mathrm{T_{\bar{p}}='+str("%.1f" % data_Tbar[j])+'\,GeV}$')
        else:
            plot.text(10, 4.5, r'$\mathrm{T_{\bar{p}}='+str("%.0f" % data_Tbar[j])+'\,GeV}$')
    
    plot.text(9e1, 9, r'$\mathrm{T_{\bar{p}}='+str("%.0f" % T_pbar)+'\,GeV}$')

    if step==7:
        plot.legend(loc='upper right', bbox_to_anchor=(0.95, 0.95), frameon=False, fontsize=0.8*label_size)
        plt.savefig('2D_uncertainty__Tp_eta__differentSetups.pdf')
    if step==71:
        plot.legend(loc='upper right', bbox_to_anchor=(0.95, 0.95), frameon=False, fontsize=0.8*label_size)
        plt.savefig('2D_uncertainty__Tp_eta__differentSetups_1.pdf')
    if step==72:
        plot.legend(loc='upper right', bbox_to_anchor=(0.95, 0.95), frameon=False, fontsize=0.8*label_size)
        plt.savefig('2D_uncertainty__Tp_eta__differentSetups_2.pdf')
    if step==73:
        plot.legend(loc='upper right', bbox_to_anchor=(0.95, 0.95), frameon=False, fontsize=0.8*label_size)
        plt.savefig('2D_uncertainty__Tp_eta__lowEnergy_0_1.pdf')
    if step==74:
        plot.legend(loc='upper right', bbox_to_anchor=(0.95, 0.95), frameon=False, fontsize=0.8*label_size)
        plt.savefig('2D_uncertainty__Tp_eta__lowEnergy_0_2.pdf')
    if step==75:
        plot.legend(loc='upper right', bbox_to_anchor=(0.95, 0.95), frameon=False, fontsize=0.8*label_size)
        plt.savefig('2D_uncertainty__Tp_eta__lowEnergy_0_5.pdf')
    if step==76:
        plot.legend(loc='upper right', bbox_to_anchor=(0.95, 0.95), frameon=False, fontsize=0.8*label_size)
        plt.savefig('2D_uncertainty__Tp_eta__lowEnergy_1_0.pdf')
    

    if step==8:
        plot.legend(loc='lower right', bbox_to_anchor=(0.95, 0.05), frameon=False, fontsize=0.8*label_size)
        plt.savefig('2D_uncertainty__Tp_eta__differentSetups_HE.pdf')
    if step==9:
        plot.legend(loc='upper right', bbox_to_anchor=(0.95, 0.95), frameon=False, fontsize=0.8*label_size)
        plt.savefig('2D_uncertainty__Tp_eta__differentSetups_ppFixed.pdf')


if step==10:
    
    
    if step==10:
        s      = 40
        option_list =['standard', 'Winkler', 'Duperray', '2_30', '3_100', 'pHe'] # ,'Duperray', '2_100'
        label_list  =['standard', 'Winkler', 'Duperray', r'2%$\,$/$\,$30%', r'3%$\,$/$\,$100%', r'p$\,$He'] # ,'Duperray', '2_100'
    if step==11:
        sys.exit(0)

    color_list  =[ c0, c1, c5, c2, c3, c4 ]
    dashes      = [(), (3,3), (6,3), (15,5), (20,5,5,5), (15,3,3,3,3,3)  ]
    
    s_int  = 0
    
    data_S=np.genfromtxt('save/save.fS.txt', skip_header=1)
    
    d = 1e90
    for i in range( len(data_S) ):
        dd = np.fabs(s*s - data_S[i])
#        print dd
#        print s*s
#        print data_S[i]
#        print ""
        if dd < d:
            d = dd
            s_int = i

    print s_int

    plot, fig = plot_1D (r'$\mathrm{x_R}$', r'$\mathrm{ p_{T}\quad [GeV]}$',      'log', 'log' )
    for j in range(len(option_list)):
        option = option_list[j]
        i = s_int
        d_ccontribution = np.genfromtxt(option+'/contribution_CM/uncertainty/2D_uncertainty__xR_pT__'+str(i)+'.txt')
        y = d_ccontribution[0,1:-2]
        x = d_ccontribution[1:-2,0]
        z = np.transpose( d_ccontribution[1:-2,1:-2] )
        smoothContour(plot, x, y, z, 'log', 'log', [0.,0.5], fc=cb, ec=color_list[j], dashes=dashes[j], label=label_list[j])

    plot.set_xlim(1e-2,2e0)
    plot.set_ylim(2e-2,2e1)

    plot.text(1.5e-2, 8, r'$\mathrm{\sqrt{s}='+str("%.0f" % s)+'\,GeV}$')
    
    if step==10:
        plot.legend(loc='upper right', bbox_to_anchor=(0.95, 0.95), frameon=False, fontsize=0.8*label_size)
        plt.savefig('2D_uncertainty__xR_eta__differentSetups.pdf')

if step==11:
    
    if step==11:
        T_p         = 4e3
        option_list =['pHe', 'pp_fixed', 'pp_fixed_6_60'] # ,'Duperray', '2_100'
        label_list  =[r'standard (pHe)', r'pp fixed', r'pp fixed, 6%$\,$/$\,$60%']
    
    color_list  =[ c0, c1, c2, c3, c4 ]
    dashes      = [(), (3,3), (15,5), (20,5,5,5), (15,3,3,3,3,3)  ]
    
    T_p_int  = 0
    
    data_Tp=np.genfromtxt(option_list[0]+'/contribution_LAB/uncertainty_fixedTp/Tp_values.txt', skip_header=1)
    
    d = 1e90
    for i in range( len(data_Tp) ):
        dd = np.fabs(T_p - data_Tp[i])
        if dd < d:
            d = dd
            T_p_int = i

    print T_p_int
    
    plot, fig = plot_1D (r'$\mathrm{T_{\bar{p}}\,[GeV]}$', r'$\mathrm{\eta}$',      'log', 'linear' )
    for j in range(len(option_list)):
        option = option_list[j]
        i = T_p_int
        d_ccontribution = np.genfromtxt(option+'/contribution_LAB/uncertainty_fixedTp/2D_uncertainty_Tpbar_eta__'+str(i)+'.txt')
        y = d_ccontribution[0,1:-2]
        x = d_ccontribution[1:-2,0]
        z = np.transpose( d_ccontribution[1:-2,1:-2] )
        smoothContour(plot, x, y, z, 'log', 'linear', [0.,0.5], fc=cb, ec=color_list[j], dashes=dashes[j], label=label_list[j])

    plot.set_xlim(5e-1,5e2)
    plot.set_ylim(2,9)

    tt = T_p/1000
    plot.text(1.5, 3.2, r'$\mathrm{T_{{p}}='+str("%.0f" % tt)+'\,TeV}$')

#    lhcb = mpatches.Rectangle( (10, 2.5), 90,  2.5, label='LHCb', linewidth=2, linestyle='dashed', ec='black', fc=(1,0,0,0.2) )
#    plot.add_patch(  lhcb  )


#    lhcb_1_p  = np.ones(20)*1
#    lhcb_1_pT = np.arrange(np.log(0.1), np.log(5)+np.log(5/0.1)/40., np.log(5/0.1)/20.)
#    lhcb_1_pT = np.exp(lhcb_1_pT)
#
#    lhcb_1_T   = np.sqrt( lhcb_1_p*lhcb_1_p+0.938^2)-0.938
#    lhcb_1_eta = np.acosh(lhcb_1_p/lhcb_1_pT)


    if step==11:
        plt.legend(loc='upper left', bbox_to_anchor=(0.05, 0.95), frameon=False, fontsize=0.8*label_size)
        plt.savefig('2D_uncertainty__Tpbar_eta__differentSetups_ppFixed.pdf')





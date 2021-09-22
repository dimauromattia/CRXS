#! /usr/bin/env python

import glob
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
from matplotlib import rc



import argparse
from PlotFunctions import plot2D


parser = argparse.ArgumentParser(description='Michael Korsmeier - Plotting script for source terms.'        )
parser.add_argument('--file',  help='Filename of table.', action='store', dest='file', type=str, default='' )
parser.add_argument('--label', help='Label.', action='store', dest='label', type=str, default=''            )
args = parser.parse_args()
file   	    = args.file
label 	    = args.label



print_size=10
label_size=30

data_diMauro = np.genfromtxt( 'pbar_SourceTerm_diMauro.txt',   skip_header=1 )
data_TanNg   = np.genfromtxt( 'pbar_SourceTerm_TanNg.txt',     skip_header=1 )
data_Winkler = np.genfromtxt( 'pbar_SourceTerm_Winkler.txt',   skip_header=1 )
data_Dup     = np.genfromtxt( 'pbar_SourceTerm_Duperray.txt',  skip_header=1 )
data_KR      = np.genfromtxt( 'pbar_SourceTerm_KR.txt',        skip_header=1 )
data_DM      = np.genfromtxt( 'pbar_SourceTerm_DM.txt',        skip_header=1 )

data_DTUNUC = np.genfromtxt( 'pbar_SourceTerm_DTUNUC.txt',        skip_header=1 )




if file!='':
    data_File    = np.genfromtxt(  file,                           skip_header=1 )


data_Dbar    = np.genfromtxt( 'Dbar_SourceTerm.txt',           skip_header=1 )
data_Dbar_DM = np.genfromtxt( 'Dbar_SourceTerm_DM.txt',        skip_header=1 )




##############################
##############################
#   Anti Deuteron
##############################
##############################



plt.close('all')
font_props = {"size":label_size}
rc("font", **font_props)
mpl.rcParams['axes.linewidth'] = 4
mpl.rcParams['mathtext.fontset']='stixsans'
fig = plt.figure(figsize=(print_size*1.5, print_size))

upperPlt = plt.subplot2grid((1, 1), (0, 0), rowspan=1)
plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)

T_power = 0;

upperPlt.errorbar( data_Dbar[:,0], data_Dbar[:, 1]*np.power(data_Dbar[:,0],T_power),   lw=4,  color='black'   , dashes=(        ), label=r'all'   )
upperPlt.errorbar( data_Dbar[:,0], data_Dbar[:, 2]*np.power(data_Dbar[:,0],T_power),   lw=2,  color='black'   , dashes=(        ), label=r'p H'   )
upperPlt.errorbar( data_Dbar[:,0], data_Dbar[:, 3]*np.power(data_Dbar[:,0],T_power),   lw=2,  color='black'   , dashes=(3, 3    ), label=r'p He'  )
upperPlt.errorbar( data_Dbar[:,0], data_Dbar[:, 4]*np.power(data_Dbar[:,0],T_power),   lw=2,  color='black'   , dashes=(14,3    ), label=r'He H'  )
upperPlt.errorbar( data_Dbar[:,0], data_Dbar[:, 5]*np.power(data_Dbar[:,0],T_power),   lw=2,  color='black'   , dashes=(20,5,5,5), label=r'He He' )

upperPlt.errorbar( data_Dbar[:,0], data_Dbar[:, 6]*np.power(data_Dbar[:,0],T_power),   lw=2,  color='black'   , dashes=(6, 6     ), label=r'$\mathrm{\bar{p}}$ H'   )
upperPlt.errorbar( data_Dbar[:,0], data_Dbar[:, 7]*np.power(data_Dbar[:,0],T_power),   lw=2,  color='black'   , dashes=(14,3,10,5), label=r'$\mathrm{\bar{p}}$ He'  )

upperPlt.set_xscale('log')
upperPlt.set_yscale('log')
upperPlt.set_xlim(2e-1, 5e2)
upperPlt.set_ylim(1e-32, 2e-28)

upperPlt.tick_params('both', length=20, width=2, which='major')
upperPlt.tick_params('both', length=10, width=1, which='minor')
upperPlt.grid(b=True, which='major', alpha=0.1, linestyle='-', linewidth=2)
upperPlt.tick_params(axis='both', pad=10)

upperPlt.set_xlabel(r'$\mathrm{T_{\bar{D}}/n [GeV/n]}$')
upperPlt.set_ylabel(r'$\mathrm{q^{(\bar{D})}  \, [(GeV/n)^{-1}m^{-3}s^{-1}]}}$')

upperPlt.legend(loc='upper right', bbox_to_anchor=(0.95, 0.95), frameon=False, fontsize=0.8*label_size)

plt.savefig('Dbar_source.pdf')

upperPlt.errorbar( data_Dbar_DM[:,0], data_Dbar_DM[:, 1]*np.power(data_Dbar[:,0],T_power),   lw=3,  color='red'   , dashes=(), label=r'$\mathrm{\bar{p}}$ He'  )

plt.savefig('Dbar_source_DM.pdf')



##############################
##############################
#   Anti Proton
##############################
##############################


#   di Mauro
##############################


plt.close('all')
font_props = {"size":label_size}
rc("font", **font_props)
mpl.rcParams['axes.linewidth'] = 4
mpl.rcParams['mathtext.fontset']='stixsans'
fig = plt.figure(figsize=(print_size*1.5, print_size))

upperPlt = plt.subplot2grid((1, 1), (0, 0), rowspan=1)
plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)

upperPlt.errorbar( data_diMauro[:,0], data_diMauro[:, 1]*np.power(data_diMauro[:,0],2.7),   lw=4,  color='black'   , dashes=(        ), label=r'all'   )
upperPlt.errorbar( data_diMauro[:,0], data_diMauro[:, 2]*np.power(data_diMauro[:,0],2.7),   lw=2,  color='black'   , dashes=(        ), label=r'p H'   )
upperPlt.errorbar( data_diMauro[:,0], data_diMauro[:, 3]*np.power(data_diMauro[:,0],2.7),   lw=2,  color='black'   , dashes=(3, 3    ), label=r'p He'  )
upperPlt.errorbar( data_diMauro[:,0], data_diMauro[:, 4]*np.power(data_diMauro[:,0],2.7),   lw=2,  color='black'   , dashes=(14,3    ), label=r'He H'  )
upperPlt.errorbar( data_diMauro[:,0], data_diMauro[:, 5]*np.power(data_diMauro[:,0],2.7),   lw=2,  color='black'   , dashes=(20,5,5,5), label=r'He He' )


upperPlt.set_xscale('log')
upperPlt.set_yscale('log')
upperPlt.set_xlim(2e-1, 1e3)
upperPlt.set_ylim(1e-26, 2e-20)

upperPlt.tick_params('both', length=20, width=2, which='major')
upperPlt.tick_params('both', length=10, width=1, which='minor')
upperPlt.grid(b=True, which='major', alpha=0.1, linestyle='-', linewidth=2)
upperPlt.tick_params(axis='both', pad=10)

upperPlt.set_xlabel(r'$\mathrm{T_{\bar{p}} [GeV]}$')
upperPlt.set_ylabel(r'$\mathrm{q^{(\bar{p})} \, T_{\bar{p}}^{2.7}  \,  \, [GeV^{1.7}m^{-3}s^{-1}]}}$')

upperPlt.legend(loc='lower right', bbox_to_anchor=(0.95, 0.05), frameon=False, fontsize=0.8*label_size)

plt.savefig('pbar_source_diMauro.pdf')

upperPlt.errorbar( data_DM[:,0], data_DM[:, 1]*np.power(data_Dbar[:,0],2.7),   lw=3,  color='red'   , dashes=(), label=r'$\mathrm{\bar{p}}$ He'  )

plt.savefig('pbar_source_diMauro_DM.pdf')



plt.close('all')
font_props = {"size":label_size}
rc("font", **font_props)
mpl.rcParams['axes.linewidth'] = 4
mpl.rcParams['mathtext.fontset']='stixsans'
fig = plt.figure(figsize=(print_size*1.5, print_size))

upperPlt = plt.subplot2grid((1, 1), (0, 0), rowspan=1)
plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)

upperPlt.errorbar( data_diMauro[:,0], data_diMauro[:, 2]/data_diMauro[:,1],   lw=2,  color='black'   , dashes=(        ), label=r'p H'   )
upperPlt.errorbar( data_diMauro[:,0], data_diMauro[:, 3]/data_diMauro[:,1],   lw=2,  color='black'   , dashes=(3, 3    ), label=r'p He'  )
upperPlt.errorbar( data_diMauro[:,0], data_diMauro[:, 4]/data_diMauro[:,1],   lw=2,  color='black'   , dashes=(14,3    ), label=r'He H'  )
upperPlt.errorbar( data_diMauro[:,0], data_diMauro[:, 5]/data_diMauro[:,1],   lw=2,  color='black'   , dashes=(20,5,5,5), label=r'He He' )

upperPlt.set_xscale('log')
upperPlt.set_yscale('linear')
upperPlt.set_xlim(2e-1, 1e3)
upperPlt.set_ylim(0, 1)

upperPlt.tick_params('both', length=20, width=2, which='major')
upperPlt.tick_params('both', length=10, width=1, which='minor')
upperPlt.grid(b=True, which='major', alpha=0.1, linestyle='-', linewidth=2)
upperPlt.tick_params(axis='both', pad=10)

upperPlt.set_xlabel(r'$\mathrm{T_{\bar{p}} [GeV]}$')
upperPlt.set_ylabel(r'relative contribution')

upperPlt.legend(loc='lower right', bbox_to_anchor=(0.95, 0.15), frameon=False, fontsize=0.8*label_size)

plt.savefig('pbar_source_contribution_diMauro.pdf')





#   Winkler
##############################

plt.close('all')
font_props = {"size":label_size}
rc("font", **font_props)
mpl.rcParams['axes.linewidth'] = 4
mpl.rcParams['mathtext.fontset']='stixsans'
fig = plt.figure(figsize=(print_size*1.5, print_size))

upperPlt = plt.subplot2grid((1, 1), (0, 0), rowspan=1)
plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)

upperPlt.errorbar( data_Winkler[:,0], data_Winkler[:, 1]*np.power(data_Winkler[:,0],2.7),   lw=4,  color='black'   , dashes=(        ), label=r'all'   )
upperPlt.errorbar( data_Winkler[:,0], data_Winkler[:, 2]*np.power(data_Winkler[:,0],2.7),   lw=2,  color='black'   , dashes=(        ), label=r'p H'   )
upperPlt.errorbar( data_Winkler[:,0], data_Winkler[:, 3]*np.power(data_Winkler[:,0],2.7),   lw=2,  color='black'   , dashes=(3, 3    ), label=r'p He'  )
upperPlt.errorbar( data_Winkler[:,0], data_Winkler[:, 4]*np.power(data_Winkler[:,0],2.7),   lw=2,  color='black'   , dashes=(14,3    ), label=r'He H'  )
upperPlt.errorbar( data_Winkler[:,0], data_Winkler[:, 5]*np.power(data_Winkler[:,0],2.7),   lw=2,  color='black'   , dashes=(20,5,5,5), label=r'He He' )


upperPlt.set_xscale('log')
upperPlt.set_yscale('log')
upperPlt.set_xlim(2e-1, 1e3)
upperPlt.set_ylim(1e-26, 2e-20)

upperPlt.tick_params('both', length=20, width=2, which='major')
upperPlt.tick_params('both', length=10, width=1, which='minor')
upperPlt.grid(b=True, which='major', alpha=0.1, linestyle='-', linewidth=2)
upperPlt.tick_params(axis='both', pad=10)

upperPlt.set_xlabel(r'$\mathrm{T_{\bar{p}} [GeV]}$')
upperPlt.set_ylabel(r'$\mathrm{q^{(\bar{p})} \, T_{\bar{p}}^{2.7}  \,  \, [GeV^{1.7}m^{-3}s^{-1}]}}$')

upperPlt.legend(loc='lower right', bbox_to_anchor=(0.95, 0.05), frameon=False, fontsize=0.8*label_size)

plt.savefig('pbar_source_Winkler.pdf')



#   Tan and Ng
##############################

plt.close('all')
font_props = {"size":label_size}
rc("font", **font_props)
mpl.rcParams['axes.linewidth'] = 4
mpl.rcParams['mathtext.fontset']='stixsans'
fig = plt.figure(figsize=(print_size*1.5, print_size))

upperPlt = plt.subplot2grid((1, 1), (0, 0), rowspan=1)
plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)

upperPlt.errorbar( data_TanNg[:,0], data_TanNg[:, 1]*np.power(data_TanNg[:,0],2.7),   lw=4,  color='black'   , dashes=(        ), label=r'all'   )
upperPlt.errorbar( data_TanNg[:,0], data_TanNg[:, 2]*np.power(data_TanNg[:,0],2.7),   lw=2,  color='black'   , dashes=(        ), label=r'p H'   )
upperPlt.errorbar( data_TanNg[:,0], data_TanNg[:, 3]*np.power(data_TanNg[:,0],2.7),   lw=2,  color='black'   , dashes=(3, 3    ), label=r'p He'  )
upperPlt.errorbar( data_TanNg[:,0], data_TanNg[:, 4]*np.power(data_TanNg[:,0],2.7),   lw=2,  color='black'   , dashes=(14,3    ), label=r'He H'  )
upperPlt.errorbar( data_TanNg[:,0], data_TanNg[:, 5]*np.power(data_TanNg[:,0],2.7),   lw=2,  color='black'   , dashes=(20,5,5,5), label=r'He He' )


upperPlt.set_xscale('log')
upperPlt.set_yscale('log')
upperPlt.set_xlim(2e-1, 1e3)
upperPlt.set_ylim(1e-26, 2e-20)

upperPlt.tick_params('both', length=20, width=2, which='major')
upperPlt.tick_params('both', length=10, width=1, which='minor')
upperPlt.grid(b=True, which='major', alpha=0.1, linestyle='-', linewidth=2)
upperPlt.tick_params(axis='both', pad=10)

upperPlt.set_xlabel(r'$\mathrm{T_{\bar{p}} [GeV]}$')
upperPlt.set_ylabel(r'$\mathrm{q^{(\bar{p})} \, T_{\bar{p}}^{2.7}  \,  \, [GeV^{1.7}m^{-3}s^{-1}]}}$')

upperPlt.legend(loc='lower right', bbox_to_anchor=(0.95, 0.05), frameon=False, fontsize=0.8*label_size)

plt.savefig('pbar_source_TanNg.pdf')


#   Duperray
##############################


plt.close('all')
font_props = {"size":label_size}
rc("font", **font_props)
mpl.rcParams['axes.linewidth'] = 4
mpl.rcParams['mathtext.fontset']='stixsans'
fig = plt.figure(figsize=(print_size*1.5, print_size))

upperPlt = plt.subplot2grid((1, 1), (0, 0), rowspan=1)
plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)

upperPlt.errorbar( data_Dup[:,0], data_Dup[:, 1]*np.power(data_Dup[:,0],2.7),   lw=4,  color='black'   , dashes=(        ), label=r'all'   )
upperPlt.errorbar( data_Dup[:,0], data_Dup[:, 2]*np.power(data_Dup[:,0],2.7),   lw=2,  color='black'   , dashes=(        ), label=r'p H'   )
upperPlt.errorbar( data_Dup[:,0], data_Dup[:, 3]*np.power(data_Dup[:,0],2.7),   lw=2,  color='black'   , dashes=(3, 3    ), label=r'p He'  )
upperPlt.errorbar( data_Dup[:,0], data_Dup[:, 4]*np.power(data_Dup[:,0],2.7),   lw=2,  color='black'   , dashes=(14,3    ), label=r'He H'  )
upperPlt.errorbar( data_Dup[:,0], data_Dup[:, 5]*np.power(data_Dup[:,0],2.7),   lw=2,  color='black'   , dashes=(20,5,5,5), label=r'He He' )


upperPlt.set_xscale('log')
upperPlt.set_yscale('log')
upperPlt.set_xlim(2e-1, 1e3)
upperPlt.set_ylim(1e-26, 2e-20)

upperPlt.tick_params('both', length=20, width=2, which='major')
upperPlt.tick_params('both', length=10, width=1, which='minor')
upperPlt.grid(b=True, which='major', alpha=0.1, linestyle='-', linewidth=2)
upperPlt.tick_params(axis='both', pad=10)

upperPlt.set_xlabel(r'$\mathrm{T_{\bar{p}} [GeV]}$')
upperPlt.set_ylabel(r'$\mathrm{q^{(\bar{p})} \, T_{\bar{p}}^{2.7}  \,  \, [GeV^{1.7}m^{-3}s^{-1}]}}$')

upperPlt.legend(loc='lower right', bbox_to_anchor=(0.95, 0.05), frameon=False, fontsize=0.8*label_size)

plt.savefig('pbar_source_Duperray.pdf')


plt.close('all')
font_props = {"size":label_size}
rc("font", **font_props)
mpl.rcParams['axes.linewidth'] = 4
mpl.rcParams['mathtext.fontset']='stixsans'
fig = plt.figure(figsize=(print_size*1.5, print_size))

upperPlt = plt.subplot2grid((1, 1), (0, 0), rowspan=1)
plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)

upperPlt.errorbar( data_Dup[:,0], data_Dup[:, 2]/data_Dup[:,1],   lw=2,  color='black'   , dashes=(        ), label=r'p H'   )
upperPlt.errorbar( data_Dup[:,0], data_Dup[:, 3]/data_Dup[:,1],   lw=2,  color='black'   , dashes=(3, 3    ), label=r'p He'  )
upperPlt.errorbar( data_Dup[:,0], data_Dup[:, 4]/data_Dup[:,1],   lw=2,  color='black'   , dashes=(14,3    ), label=r'He H'  )
upperPlt.errorbar( data_Dup[:,0], data_Dup[:, 5]/data_Dup[:,1],   lw=2,  color='black'   , dashes=(20,5,5,5), label=r'He He' )

upperPlt.set_xscale('log')
upperPlt.set_yscale('linear')
upperPlt.set_xlim(2e-1, 1e3)
upperPlt.set_ylim(0, 1)

upperPlt.tick_params('both', length=20, width=2, which='major')
upperPlt.tick_params('both', length=10, width=1, which='minor')
upperPlt.grid(b=True, which='major', alpha=0.1, linestyle='-', linewidth=2)
upperPlt.tick_params(axis='both', pad=10)

upperPlt.set_xlabel(r'$\mathrm{T_{\bar{p}} [GeV]}$')
upperPlt.set_ylabel(r'relative contribution')

upperPlt.legend(loc='lower right', bbox_to_anchor=(0.95, 0.15), frameon=False, fontsize=0.8*label_size)

plt.savefig('pbar_source_contribution_Duperray.pdf')



#   Kachelriess
##############################


plt.close('all')
font_props = {"size":label_size}
rc("font", **font_props)
mpl.rcParams['axes.linewidth'] = 4
mpl.rcParams['mathtext.fontset']='stixsans'
fig = plt.figure(figsize=(print_size*1.5, print_size))

upperPlt = plt.subplot2grid((1, 1), (0, 0), rowspan=1)
plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)

upperPlt.errorbar( data_KR    [:,0], data_KR    [:, 1]*np.power(data_KR    [:,0],2.7),   lw=4,  color='black'   , dashes=(        ), label=r'all'   )
upperPlt.errorbar( data_KR    [:,0], data_KR    [:, 2]*np.power(data_KR    [:,0],2.7),   lw=2,  color='black'   , dashes=(        ), label=r'p H'   )
upperPlt.errorbar( data_KR    [:,0], data_KR    [:, 3]*np.power(data_KR    [:,0],2.7),   lw=2,  color='black'   , dashes=(3, 3    ), label=r'p He'  )
upperPlt.errorbar( data_KR    [:,0], data_KR    [:, 4]*np.power(data_KR    [:,0],2.7),   lw=2,  color='black'   , dashes=(14,3    ), label=r'He H'  )
upperPlt.errorbar( data_KR    [:,0], data_KR    [:, 5]*np.power(data_KR    [:,0],2.7),   lw=2,  color='black'   , dashes=(20,5,5,5), label=r'He He' )

upperPlt.set_xscale('log')
upperPlt.set_yscale('log')
upperPlt.set_xlim(2e-1, 1e3)
upperPlt.set_ylim(1e-26, 2e-20)

upperPlt.tick_params('both', length=20, width=2, which='major')
upperPlt.tick_params('both', length=10, width=1, which='minor')
upperPlt.grid(b=True, which='major', alpha=0.1, linestyle='-', linewidth=2)
upperPlt.tick_params(axis='both', pad=10)

upperPlt.set_xlabel(r'$\mathrm{T_{\bar{p}} [GeV]}$')
upperPlt.set_ylabel(r'$\mathrm{q^{(\bar{p})} \, T_{\bar{p}}^{2.7}  \,  \, [GeV^{1.7}m^{-3}s^{-1}]}}$')

upperPlt.legend(loc='lower right', bbox_to_anchor=(0.95, 0.05), frameon=False, fontsize=0.8*label_size)

plt.savefig('pbar_source_KR.pdf')


plt.close('all')
font_props = {"size":label_size}
rc("font", **font_props)
mpl.rcParams['axes.linewidth'] = 4
mpl.rcParams['mathtext.fontset']='stixsans'
fig = plt.figure(figsize=(print_size*1.5, print_size))

upperPlt = plt.subplot2grid((1, 1), (0, 0), rowspan=1)
plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)

upperPlt.errorbar( data_KR[:,0], data_KR[:, 2]/data_KR[:,1],   lw=2,  color='black'   , dashes=(        ), label=r'p H'   )
upperPlt.errorbar( data_KR[:,0], data_KR[:, 3]/data_KR[:,1],   lw=2,  color='black'   , dashes=(3, 3    ), label=r'p He'  )
upperPlt.errorbar( data_KR[:,0], data_KR[:, 4]/data_KR[:,1],   lw=2,  color='black'   , dashes=(14,3    ), label=r'He H'  )
upperPlt.errorbar( data_KR[:,0], data_KR[:, 5]/data_KR[:,1],   lw=2,  color='black'   , dashes=(20,5,5,5), label=r'He He' )

upperPlt.set_xscale('log')
upperPlt.set_yscale('linear')
upperPlt.set_xlim(2e-1, 1e3)
upperPlt.set_ylim(0, 1)

upperPlt.tick_params('both', length=20, width=2, which='major')
upperPlt.tick_params('both', length=10, width=1, which='minor')
upperPlt.grid(b=True, which='major', alpha=0.1, linestyle='-', linewidth=2)
upperPlt.tick_params(axis='both', pad=10)

upperPlt.set_xlabel(r'$\mathrm{T_{\bar{p}} [GeV]}$')
upperPlt.set_ylabel(r'relative contribution')

upperPlt.legend(loc='lower right', bbox_to_anchor=(0.95, 0.15), frameon=False, fontsize=0.8*label_size)

plt.savefig('pbar_source_contribution_KR.pdf')


#   File
##############################


if file!='':
    plt.close('all')
    font_props = {"size":label_size}
    rc("font", **font_props)
    mpl.rcParams['axes.linewidth'] = 4
    mpl.rcParams['mathtext.fontset']='stixsans'
    fig = plt.figure(figsize=(print_size*1.5, print_size))

    upperPlt = plt.subplot2grid((1, 1), (0, 0), rowspan=1)
    plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)

    upperPlt.errorbar( data_File[:,0], data_File[:, 1]*np.power(data_File[:,0],2.7),   lw=4,  color='black'   , dashes=(        ), label=r'all'   )
    upperPlt.errorbar( data_File[:,0], data_File[:, 2]*np.power(data_File[:,0],2.7),   lw=2,  color='black'   , dashes=(        ), label=r'p H'   )
    upperPlt.errorbar( data_File[:,0], data_File[:, 3]*np.power(data_File[:,0],2.7),   lw=2,  color='black'   , dashes=(3, 3    ), label=r'p He'  )
    upperPlt.errorbar( data_File[:,0], data_File[:, 4]*np.power(data_File[:,0],2.7),   lw=2,  color='black'   , dashes=(14,3    ), label=r'He H'  )
    upperPlt.errorbar( data_File[:,0], data_File[:, 5]*np.power(data_File[:,0],2.7),   lw=2,  color='black'   , dashes=(20,5,5,5), label=r'He He' )

    upperPlt.set_xscale('log')
    upperPlt.set_yscale('log')
    upperPlt.set_xlim(2e-1, 1e3)
    upperPlt.set_ylim(1e-26, 2e-20)

    upperPlt.tick_params('both', length=20, width=2, which='major')
    upperPlt.tick_params('both', length=10, width=1, which='minor')
    upperPlt.grid(b=True, which='major', alpha=0.1, linestyle='-', linewidth=2)
    upperPlt.tick_params(axis='both', pad=10)

    upperPlt.set_xlabel(r'$\mathrm{T_{\bar{p}} [GeV]}$')
    upperPlt.set_ylabel(r'$\mathrm{q^{(\bar{p})} \, T_{\bar{p}}^{2.7}  \,  \, [GeV^{1.7}m^{-3}s^{-1}]}}$')

    upperPlt.legend(loc='lower right', bbox_to_anchor=(0.95, 0.05), frameon=False, fontsize=0.8*label_size)

    plt.savefig(file.split('.')[0]+'.pdf')


    plt.close('all')
    font_props = {"size":label_size}
    rc("font", **font_props)
    mpl.rcParams['axes.linewidth'] = 4
    mpl.rcParams['mathtext.fontset']='stixsans'
    fig = plt.figure(figsize=(print_size*1.5, print_size))

    upperPlt = plt.subplot2grid((1, 1), (0, 0), rowspan=1)
    plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)

    upperPlt.errorbar( data_File[:,0], data_File[:, 2]/data_File[:,1],   lw=2,  color='black'   , dashes=(        ), label=r'p H'   )
    upperPlt.errorbar( data_File[:,0], data_File[:, 3]/data_File[:,1],   lw=2,  color='black'   , dashes=(3, 3    ), label=r'p He'  )
    upperPlt.errorbar( data_File[:,0], data_File[:, 4]/data_File[:,1],   lw=2,  color='black'   , dashes=(14,3    ), label=r'He H'  )
    upperPlt.errorbar( data_File[:,0], data_File[:, 5]/data_File[:,1],   lw=2,  color='black'   , dashes=(20,5,5,5), label=r'He He' )

    upperPlt.set_xscale('log')
    upperPlt.set_yscale('linear')
    upperPlt.set_xlim(2e-1, 1e3)
    upperPlt.set_ylim(0, 1)

    upperPlt.tick_params('both', length=20, width=2, which='major')
    upperPlt.tick_params('both', length=10, width=1, which='minor')
    upperPlt.grid(b=True, which='major', alpha=0.1, linestyle='-', linewidth=2)
    upperPlt.tick_params(axis='both', pad=10)

    upperPlt.set_xlabel(r'$\mathrm{T_{\bar{p}} [GeV]}$')
    upperPlt.set_ylabel(r'relative contribution')

    upperPlt.legend(loc='lower right', bbox_to_anchor=(0.95, 0.15), frameon=False, fontsize=0.8*label_size)

    plt.savefig('pbar_source_contribution_'+file.split('.')[0]+'.pdf')








##############################
##############################
#   Anti Proton Comparison
##############################
##############################


#   pp
##############################


plt.close('all')
font_props = {"size":label_size}
rc("font", **font_props)
mpl.rcParams['axes.linewidth'] = 4
mpl.rcParams['mathtext.fontset']='stixsans'
fig = plt.figure(figsize=(print_size*1.5, print_size))

upperPlt = plt.subplot2grid((1, 1), (0, 0), rowspan=1)
plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)

upperPlt.errorbar( data_diMauro[:,0], data_diMauro[:, 2]*np.power(data_diMauro[:,0],2.7),   lw=2,  color='black'   , dashes=(        ), label=r'di Mauro'   )
upperPlt.errorbar( data_Winkler[:,0], data_Winkler[:, 2]*np.power(data_Winkler[:,0],2.7),   lw=2,  color='black'   , dashes=(        ), label=r'Winkler'    )
upperPlt.errorbar( data_TanNg  [:,0], data_TanNg  [:, 2]*np.power(data_TanNg  [:,0],2.7),   lw=2,  color='black'   , dashes=(14,3    ), label=r'Tan&Ng'     )
upperPlt.errorbar( data_KR     [:,0], data_KR     [:, 2]*np.power(data_KR     [:,0],2.7),   lw=2,  color='black'   , dashes=(3,3     ), label=r'KMO')
upperPlt.errorbar( data_Dup    [:,0], data_Dup    [:, 2]*np.power(data_Dup    [:,0],2.7),   lw=2,  color='black'   , dashes=(20,5,5,5), label=r'Duperray'   )
if file!='':
    upperPlt.errorbar( data_File[:,0], data_File  [:, 2]*np.power(data_File   [:,0],2.7),   lw=2,  color='red'     , dashes=(        ), label=r''+label   )


upperPlt.set_xscale('log')
upperPlt.set_yscale('log')
upperPlt.set_xlim(2e-1, 1e3)
upperPlt.set_ylim(1e-26, 2e-20)

upperPlt.tick_params('both', length=20, width=2, which='major')
upperPlt.tick_params('both', length=10, width=1, which='minor')
upperPlt.grid(b=True, which='major', alpha=0.1, linestyle='-', linewidth=2)
upperPlt.tick_params(axis='both', pad=10)

upperPlt.set_xlabel(r'$\mathrm{T_{\bar{p}} [GeV]}$')
upperPlt.set_ylabel(r'$\mathrm{q^{(\bar{p})} \, T_{\bar{p}}^{2.7}  \,  \, [GeV^{1.7}m^{-3}s^{-1}]}}$')

upperPlt.legend(loc='lower right', bbox_to_anchor=(0.95, 0.05), frameon=False, fontsize=0.8*label_size)

plt.savefig('pbar_source_pp.pdf')
upperPlt.set_xlim(2e0, 1e3)
upperPlt.set_ylim(1e-26, 2e-20)
plt.savefig('pbar_source_pp_soom.pdf')

#   pHe
##############################


plt.close('all')
font_props = {"size":label_size}
rc("font", **font_props)
mpl.rcParams['axes.linewidth'] = 4
mpl.rcParams['mathtext.fontset']='stixsans'
fig = plt.figure(figsize=(print_size*1.5, print_size))

upperPlt = plt.subplot2grid((1, 1), (0, 0), rowspan=1)
plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)

upperPlt.errorbar( data_diMauro[:,0], data_diMauro[:, 3]*np.power(data_diMauro[:,0],2.7),   lw=2,  color='black'   , dashes=(        ), label=r'di Mauro'   )
upperPlt.errorbar( data_TanNg  [:,0], data_TanNg  [:, 3]*np.power(data_TanNg  [:,0],2.7),   lw=2,  color='black'   , dashes=(14,3    ), label=r'Tan&Ng'     )
upperPlt.errorbar( data_KR     [:,0], data_KR     [:, 3]*np.power(data_KR     [:,0],2.7),   lw=2,  color='black'   , dashes=(3,3     ), label=r'KMO')
upperPlt.errorbar( data_Dup    [:,0], data_Dup    [:, 3]*np.power(data_Dup    [:,0],2.7),   lw=2,  color='black'   , dashes=(20,5,5,5), label=r'Duperray'   )
if file!='':
    upperPlt.errorbar( data_File[:,0], data_File  [:, 3]*np.power(data_File   [:,0],2.7),   lw=2,  color='red'     , dashes=(        ), label=r''+label   )


upperPlt.set_xscale('log')
upperPlt.set_yscale('log')
upperPlt.set_xlim(2e-1, 1e3)
upperPlt.set_ylim(1e-26, 2e-20)

upperPlt.tick_params('both', length=20, width=2, which='major')
upperPlt.tick_params('both', length=10, width=1, which='minor')
upperPlt.grid(b=True, which='major', alpha=0.1, linestyle='-', linewidth=2)
upperPlt.tick_params(axis='both', pad=10)

upperPlt.set_xlabel(r'$\mathrm{T_{\bar{p}} [GeV]}$')
upperPlt.set_ylabel(r'$\mathrm{q^{(\bar{p})} \, T_{\bar{p}}^{2.7}  \,  \, [GeV^{1.7}m^{-3}s^{-1}]}}$')

upperPlt.legend(loc='lower right', bbox_to_anchor=(0.95, 0.05), frameon=False, fontsize=0.8*label_size)

plt.savefig('pbar_source_pHe.pdf')
upperPlt.set_xlim(2e0, 1e3)
upperPlt.set_ylim(1e-26, 2e-20)
plt.savefig('pbar_source_pHe_soom.pdf')


#   Hep
##############################



plt.close('all')
font_props = {"size":label_size}
rc("font", **font_props)
mpl.rcParams['axes.linewidth'] = 4
mpl.rcParams['mathtext.fontset']='stixsans'
fig = plt.figure(figsize=(print_size*1.5, print_size))

upperPlt = plt.subplot2grid((1, 1), (0, 0), rowspan=1)
plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)

upperPlt.errorbar( data_diMauro[:,0], data_diMauro[:, 4]*np.power(data_diMauro[:,0],2.7),   lw=2,  color='black'   , dashes=(        ), label=r'di Mauro'   )
upperPlt.errorbar( data_TanNg  [:,0], data_TanNg  [:, 4]*np.power(data_TanNg  [:,0],2.7),   lw=2,  color='black'   , dashes=(14,3    ), label=r'Tan&Ng'     )
upperPlt.errorbar( data_KR     [:,0], data_KR     [:, 4]*np.power(data_KR     [:,0],2.7),   lw=2,  color='black'   , dashes=(3,3     ), label=r'KMO')
upperPlt.errorbar( data_Dup    [:,0], data_Dup    [:, 4]*np.power(data_Dup    [:,0],2.7),   lw=2,  color='black'   , dashes=(20,5,5,5), label=r'Duperray'   )
if file!='':
    upperPlt.errorbar( data_File[:,0], data_File  [:, 4]*np.power(data_File   [:,0],2.7),   lw=2,  color='red'     , dashes=(        ), label=r''+label   )


upperPlt.set_xscale('log')
upperPlt.set_yscale('log')
upperPlt.set_xlim(2e-1, 1e3)
upperPlt.set_ylim(1e-26, 2e-20)

upperPlt.tick_params('both', length=20, width=2, which='major')
upperPlt.tick_params('both', length=10, width=1, which='minor')
upperPlt.grid(b=True, which='major', alpha=0.1, linestyle='-', linewidth=2)
upperPlt.tick_params(axis='both', pad=10)

upperPlt.set_xlabel(r'$\mathrm{T_{\bar{p}} [GeV]}$')
upperPlt.set_ylabel(r'$\mathrm{q^{(\bar{p})} \, T_{\bar{p}}^{2.7}  \,  \, [GeV^{1.7}m^{-3}s^{-1}]}}$')

upperPlt.legend(loc='lower right', bbox_to_anchor=(0.95, 0.05), frameon=False, fontsize=0.8*label_size)

plt.savefig('pbar_source_Hep.pdf')
upperPlt.set_xlim(2e0, 1e3)
upperPlt.set_ylim(1e-26, 2e-20)
plt.savefig('pbar_source_Hep_soom.pdf')



# All

c0='black'
c1='#8A0808'
c2='#0B610B'
c3='#0B0B61'


plt.close('all')
font_props = {"size":label_size}
rc("font", **font_props)
mpl.rcParams['axes.linewidth'] = 2
mpl.rcParams['mathtext.fontset']='stixsans'
fig = plt.figure(figsize=(print_size*1.5, print_size))

upperPlt = plt.subplot2grid((1, 1), (0, 0), rowspan=1)
plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)



#pp
upperPlt.errorbar( data_diMauro[:,0], data_diMauro[:, 2]*np.power(data_diMauro[:,0],2.7),   lw=3,  color=c0   , dashes=(        ), label=r'di Mauro'   )
upperPlt.errorbar( data_Winkler[:,0], data_Winkler[:, 2]*np.power(data_Winkler[:,0],2.7),   lw=3,  color=c0   , dashes=(16,4,3,3,3,4), label=r'Winkler'    )
upperPlt.errorbar( data_KR     [:,0], data_KR     [:, 2]*np.power(data_KR     [:,0],2.7),   lw=3,  color=c0   , dashes=(3,3     ), label=r'KMO')
upperPlt.errorbar( data_Dup    [:,0], data_Dup    [:, 2]*np.power(data_Dup    [:,0],2.7),   lw=3,  color=c0   , dashes=(20,5,5,5), label=r'Duperray'   )
upperPlt.errorbar( data_TanNg  [:,0], data_TanNg  [:, 2]*np.power(data_TanNg  [:,0],2.7),   lw=3,  color=c0   , dashes=(14,3    ), label=r'Tan&Ng'     )
if file!='':
    upperPlt.errorbar( data_File[:,0], data_File  [:, 2]*np.power(data_File   [:,0],2.7),   lw=2,  color='red'     , dashes=(        ), label=r''+label   )
#pp
upperPlt.errorbar( data_diMauro[:,0], data_diMauro[:, 2]*np.power(data_diMauro[:,0],2.7),   lw=3,  color=c1   , dashes=(        ) )
upperPlt.errorbar( data_Winkler[:,0], data_Winkler[:, 2]*np.power(data_Winkler[:,0],2.7),   lw=3,  color=c1   , dashes=(16,4,3,3,3,4) )
upperPlt.errorbar( data_TanNg  [:,0], data_TanNg  [:, 2]*np.power(data_TanNg  [:,0],2.7),   lw=3,  color=c1   , dashes=(14,3    ) )
upperPlt.errorbar( data_KR     [:,0], data_KR     [:, 2]*np.power(data_KR     [:,0],2.7),   lw=3,  color=c1   , dashes=(3,3     ) )
upperPlt.errorbar( data_Dup    [:,0], data_Dup    [:, 2]*np.power(data_Dup    [:,0],2.7),   lw=3,  color=c1   , dashes=(20,5,5,5) )
if file!='':
    upperPlt.errorbar( data_File[:,0], data_File  [:, 2]*np.power(data_File   [:,0],2.7),   lw=2,  color='red'     , dashes=(        ), label=r''+label   )


#pHe
upperPlt.errorbar( data_diMauro[:,0], data_diMauro[:, 3]*np.power(data_diMauro[:,0],2.7)*2e-1,   lw=3,  color=c2   , dashes=(        )  )
upperPlt.errorbar( data_Winkler[:,0], data_Winkler[:, 3]*np.power(data_Winkler[:,0],2.7)*2e-1,   lw=3,  color=c2   , dashes=(16,4,3,3,3,4)  )
#upperPlt.errorbar( data_DTUNUC [:,0], data_DTUNUC [:, 3]*np.power(data_DTUNUC [:,0],2.7)*2e-1,   lw=3,  color=c2   , dashes=(14,3    )  )
upperPlt.errorbar( data_KR     [:,0], data_KR     [:, 3]*np.power(data_KR     [:,0],2.7)*2e-1,   lw=3,  color=c2   , dashes=(3,3     )  )
upperPlt.errorbar( data_Dup    [:,0], data_Dup    [:, 3]*np.power(data_Dup    [:,0],2.7)*2e-1,   lw=3,  color=c2   , dashes=(20,5,5,5)  )
if file!='':
    upperPlt.errorbar( data_File[:,0], data_File  [:, 3]*np.power(data_File   [:,0],2.7),   lw=2,  color='red'     , dashes=(        )   )


#HeP
upperPlt.errorbar( data_diMauro[:,0], data_diMauro[:, 4]*np.power(data_diMauro[:,0],2.7)*1e-2,   lw=3,  color=c3   , dashes=(        )  )
upperPlt.errorbar( data_Winkler[:,0], data_Winkler[:, 4]*np.power(data_Winkler[:,0],2.7)*1e-2,   lw=3,  color=c3   , dashes=(16,4,3,3,3,4)  )
#upperPlt.errorbar( data_DTUNUC [:,0], data_DTUNUC [:, 4]*np.power(data_DTUNUC [:,0],2.7)*1e-2,   lw=3,  color=c3   , dashes=(14,3    )  )
upperPlt.errorbar( data_KR     [:,0], data_KR     [:, 4]*np.power(data_KR     [:,0],2.7)*1e-2,   lw=3,  color=c3   , dashes=(3,3     )  )
upperPlt.errorbar( data_Dup    [:,0], data_Dup    [:, 4]*np.power(data_Dup    [:,0],2.7)*1e-2,   lw=3,  color=c3   , dashes=(20,5,5,5)  )
if file!='':
    upperPlt.errorbar( data_File[:,0], data_File  [:, 4]*np.power(data_File   [:,0],2.7)*1e-1,   lw=2,  color='red'     , dashes=(        )  )


upperPlt.text    (   1.2e2,  2.5e-33 , r'$p\,p$',                                 color=c1, rotation=5 )
upperPlt.text    (   1.2e2, 19e-35 , r'$p\,\mathrm{He} \quad (\times 0.2)$',    color=c2, rotation=5  )
upperPlt.text    (   1.2e2, 10e-36 , r'$\mathrm{He}\,p \quad (\times 0.01)$',   color=c3, rotation=5  )


upperPlt.set_xscale('log')
upperPlt.set_yscale('log')
upperPlt.set_xlim(2e-1, 1e3)
upperPlt.set_ylim(1e-26, 2e-20)

upperPlt.tick_params('both', length=20, width=2, which='major')
upperPlt.tick_params('both', length=10, width=1, which='minor')
upperPlt.grid(b=True, which='major', alpha=0.1, linestyle='-', linewidth=2)
upperPlt.tick_params(axis='both', pad=10)

upperPlt.set_xlabel(r'$\mathrm{T_{\bar{p}} [GeV]}$', fontsize=label_size*1.4 )
upperPlt.set_ylabel(r'$\mathrm{q \, T_{\bar{p}}^{2.7}  \,  \, [GeV^{1.7}m^{-3}s^{-1}]}}$', fontsize=label_size*1.4 )

upperPlt.legend(loc='lower right', bbox_to_anchor=(0.95, 0.05), frameon=False, fontsize=0.7*label_size)

plt.savefig('pbar_source_all.pdf')






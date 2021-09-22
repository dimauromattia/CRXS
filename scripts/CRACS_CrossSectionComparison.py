#! /usr/bin/env python

import sys
import numpy as np

sys.path.append("/usr/local/software/antiDeuteronM/scripts")
from PlotFunctions import plot2D
from PlotFunctions import profile

import argparse
from PlotFunctions import plot2D


parser = argparse.ArgumentParser(description='Michael Korsmeier - Plotting script for source terms.')
parser.add_argument('--file_end', help='Filename of table, e.g. KR for dT_Hep_pbar_LAB_KR.txt', action='store', dest='file', type=str, default='')
parser.add_argument('--label',    help='Label.',   action='store',   dest='label',  type=str,   default='')
parser.add_argument('--factor',   help='Factor .', action='store',   dest='factor', type=float, default=2.3)
args = parser.parse_args()
file   	    = args.file
label 	    = args.label
fac         = args.factor

pp_file  = 'dT_pp_pbar_LAB_'  + file + '.txt'
Hep_file = 'dT_Hep_pbar_LAB_' + file + '.txt'
pHe_file = 'dT_pHe_pbar_LAB_' + file + '.txt'

HeHe_file= 'dT_HeHe_pbar_LAB_' + file + '.txt'


plot2D('dT_pp_Hebar_LAB.txt',            resfile='dT_pp_Hebar_LAB.png'     )

plot2D('dT_HeHe_pbar_LAB_diMauro12.txt',            resfile='dT_HeHe_pbar_LAB_diMauro12.png'     )

plot2D('dT_pp_pbar_LAB_WinklerWithHypWitNbar.txt',             resfile='dT_pp_pbar_LAB_WinklerWithHypWitNbar.png'       )
plot2D('dT_HeHe_pbar_LAB_WinklerWithHypWitNbar.txt',           resfile='dT_HeHe_pbar_LAB_WinklerWithHypWitNbar.png'     )
#plot2D('dT_pHe_pbar_LAB_WinklerWithHypWitNbar.txt',            resfile='dT_pHe_pbar_LAB_WinklerWithHypWitNbar.png'      )
#plot2D('dT_Hep_pbar_LAB_WinklerWithHypWitNbar.txt',            resfile='dT_Hep_pbar_LAB_WinklerWithHypWitNbar.png'      )



#
#   He He
#

p_proj = profile('dT_HeHe_pbar_LAB_DTUNUC_USIN.txt',  20, resfile='proj_dT_HeHe_pbar.png', type='proj', dashes=(3, 3    ), Tn_min=1e-1, Tn_max=1e4, y_fac = 2.3,  label='DTUNUC'   )
profile(         'dT_HeHe_pbar_LAB_Duperray_USIN.txt',20, resfile='proj_dT_HeHe_pbar.png', type='proj', dashes=(20,5,5,5), Tn_min=1e-1, Tn_max=9e1, y_fac = 2.3,   label='Duperray', draw_on_top=True )
profile(         'dT_HeHe_pbar_LAB_diMauro12.txt',    20, resfile='proj_dT_HeHe_pbar.png', type='proj', dashes=(        ), Tn_min=1e-1, Tn_max=1e4, y_fac = 2.3 ,   label='di Mauro', draw_on_top=True )
profile(         'dT_HeHe_pbar_LAB_WinklerWithHypWitNbar.txt',    20, resfile='proj_dT_HeHe_pbar.png', type='proj', dashes=(10,10), Tn_min=1e-1, Tn_max=1e4, y_fac = 2.3 ,   label='Winkler', draw_on_top=True )
if file!='':
    profile(     HeHe_file,                           20, resfile='proj_dT_HeHe_pbar.png', type='proj', dashes=(        ), Tn_min=1e-1, Tn_max=1e4, y_fac = fac,   label=label,      draw_on_top=True, color='red', alpha=0.5 )


profile(        'dT_HeHe_pbar_LAB_DTUNUC_USIN.txt',  500  , resfile='proj_dT_HeHe_pbar.png', type='proj', dashes=(3, 3    ), Tn_min=1e-1, Tn_max=9e1, y_fac = 2.3,      draw_on_top=True )
profile(        'dT_HeHe_pbar_LAB_Duperray_USIN.txt',500  , resfile='proj_dT_HeHe_pbar.png', type='proj', dashes=(20,5,5,5), Tn_min=1e-1, Tn_max=9e1, y_fac = 2.3,      draw_on_top=True )
profile(        'dT_HeHe_pbar_LAB_diMauro12.txt',    500  , resfile='proj_dT_HeHe_pbar.png', type='proj', dashes=(        ), Tn_min=1e-1, Tn_max=1e4, y_fac = 2.3,      draw_on_top=True )
profile(        'dT_HeHe_pbar_LAB_WinklerWithHypWitNbar.txt',    500, resfile='proj_dT_HeHe_pbar.png', type='proj', dashes=(10,10), Tn_min=1e-1, Tn_max=1e4, y_fac = 2.3 ,   draw_on_top=True )
if file!='':
    profile(    HeHe_file,                           500  , resfile='proj_dT_HeHe_pbar.png', type='proj', dashes=(        ), Tn_min=1e-1, Tn_max=1e4, y_fac = fac,      draw_on_top=True, color='red', alpha=0.5,  )


p_proj.text    (   8e0, 1e-32 , r'$\mathrm{T_{He}/4=                   20\,GeV}$',   color='black',     rotation=-70 )
p_proj.text    (   1e2, 2e-32 , r'$\mathrm{T_{He}/4=                  500\,GeV}$',   color='black',     rotation=-70 )
p_proj.text    (   8e2, 5e-32 , r'$\mathrm{T_{He}/4=  8\,TeV\,\,(\times\,10)  }$',   color='black',     rotation=-60 )



if file!='':
    profile( HeHe_file,                           8000  , resfile='proj_dT_HeHe_pbar.png', type='proj', dashes=(        ), Tn_min=1e-1, Tn_max=1e4, y_fac = 10*fac, draw_on_top=True, color='red', alpha=0.5 )
profile(     'dT_HeHe_pbar_LAB_DTUNUC_USIN.txt',  8000  , resfile='proj_dT_HeHe_pbar.png', type='proj', dashes=(3, 3    ), Tn_min=1e-1, Tn_max=9e1, y_fac = 10*2.3, draw_on_top=True)
profile(     'dT_HeHe_pbar_LAB_Duperray_USIN.txt',8000  , resfile='proj_dT_HeHe_pbar.png', type='proj', dashes=(20,5,5,5), Tn_min=1e-1, Tn_max=9e1, y_fac = 10*2.3, draw_on_top=True)
profile(     'dT_HeHe_pbar_LAB_WinklerWithHypWitNbar.txt',8000  , resfile='proj_dT_HeHe_pbar.png', type='proj', dashes=(10,10), Tn_min=1e-1, Tn_max=9e1, y_fac = 10*2.3, draw_on_top=True)
profile(     'dT_HeHe_pbar_LAB_diMauro12.txt',    8000  , resfile='proj_dT_HeHe_pbar.png', type='proj', dashes=(        ), Tn_min=1e-1, Tn_max=1e4, y_fac = 10*2.3, draw_on_top=True, legend=True, x_label = r'$\mathrm{T_{\bar{p}} \quad [GeV]}$', y_min=5e-34, y_max=1e-29 )




p_prod = profile('dT_HeHe_pbar_LAB_DTUNUC_USIN.txt',  0.5, resfile='prod_dT_HeHe_pbar.png', type='prod', dashes=(3, 3    ), Tn_min=5, Tn_max=1e4, y_fac = 2.3, label='DTUNUC'  )
profile(         'dT_HeHe_pbar_LAB_Duperray_USIN.txt',0.5, resfile='prod_dT_HeHe_pbar.png', type='prod', dashes=(20,5,5,5), Tn_min=5, Tn_max=1e4, y_fac = 2.3, label='Duperray' , draw_on_top=True)
profile(         'dT_HeHe_pbar_LAB_WinklerWithHypWitNbar.txt',    0.5, resfile='prod_dT_HeHe_pbar.png', type='prod', dashes=(10,10), Tn_min=5, Tn_max=1e4, y_fac = 2.3, label='Winkler' , draw_on_top=True)
profile(         'dT_HeHe_pbar_LAB_diMauro12.txt',    0.5, resfile='prod_dT_HeHe_pbar.png', type='prod', dashes=(        ), Tn_min=5, Tn_max=1e4, y_fac = 2.3, label='di Mauro' , draw_on_top=True)
if file!='':
    profile(     HeHe_file,                           0.5, resfile='prod_dT_HeHe_pbar.png', type='prod', dashes=(        ), Tn_min=5, Tn_max=1e4, y_fac = fac, label=label,       draw_on_top=True,  color='red', alpha=0.5 )


if file!='':
    profile(HeHe_file,                            10  , resfile='prod_dT_HeHe_pbar.png', type='prod', dashes=(        ), Tn_min=5, Tn_max=1e4, y_fac = fac, draw_on_top=True, color='red', alpha=0.5  )
profile(    'dT_HeHe_pbar_LAB_DTUNUC_USIN.txt',   10  , resfile='prod_dT_HeHe_pbar.png', type='prod', dashes=(3, 3    ), Tn_min=5, Tn_max=1e4, y_fac = 2.3, draw_on_top=True)
profile(    'dT_HeHe_pbar_LAB_Duperray_USIN.txt', 10  , resfile='prod_dT_HeHe_pbar.png', type='prod', dashes=(20,5,5,5), Tn_min=5, Tn_max=1e4, y_fac = 2.3, draw_on_top=True)
profile(    'dT_HeHe_pbar_LAB_diMauro12.txt',     10  , resfile='prod_dT_HeHe_pbar.png', type='prod', dashes=(        ), Tn_min=5, Tn_max=1e4, y_fac = 2.3, draw_on_top=True)
profile(    'dT_HeHe_pbar_LAB_WinklerWithHypWitNbar.txt',     10  , resfile='prod_dT_HeHe_pbar.png', type='prod', dashes=(10,10   ), Tn_min=5, Tn_max=1e4, y_fac = 2.3, draw_on_top=True)


p_prod.text    (  9e0, 2e-34 , r'$\mathrm{T_{\bar{p}}=                   0.5\,GeV}$',   color='black', rotation=80 )
p_prod.text    (  2e1, 2e-34 , r'$\mathrm{T_{\bar{p}}=                    10\,GeV}$',   color='black', rotation=80 )
p_prod.text    (  9e1, 2e-34 , r'$\mathrm{T_{\bar{p}}=                    90\,GeV}$',   color='black', rotation=80 )

if file!='':
    profile(HeHe_file,                           90  , resfile='prod_dT_HeHe_pbar.png', type='prod', dashes=(        ), Tn_min=5, Tn_max=1e4, y_fac = fac, draw_on_top=True,  color='red', alpha=0.5)
profile(    'dT_HeHe_pbar_LAB_DTUNUC_USIN.txt',  90  , resfile='prod_dT_HeHe_pbar.png', type='prod', dashes=(3, 3    ), Tn_min=5, Tn_max=1e4, y_fac = 2.3, draw_on_top=True)
profile(    'dT_HeHe_pbar_LAB_diMauro12.txt',    90  , resfile='prod_dT_HeHe_pbar.png', type='prod', dashes=(        ), Tn_min=5, Tn_max=1e4, y_fac = 2.3, draw_on_top=True)
profile(    'dT_HeHe_pbar_LAB_WinklerWithHypWitNbar.txt',    90  , resfile='prod_dT_HeHe_pbar.png', type='prod', dashes=(10,10), Tn_min=5, Tn_max=1e4, y_fac = 2.3, draw_on_top=True)
profile(    'dT_HeHe_pbar_LAB_Duperray_USIN.txt',90  , resfile='prod_dT_HeHe_pbar.png', type='prod', dashes=(20,5,5,5), Tn_min=5, Tn_max=1e4, y_fac = 2.3, draw_on_top=True, legend=True, x_label = r'$\mathrm{T_{He}/4 \quad [GeV/4]}$' )


#
#   He p
#

p_proj = profile('dT_Hep_pbar_LAB_DTUNUC_USIN.txt',  20, resfile='proj_dT_Hep_pbar.png', type='proj', dashes=(3, 3    ), Tn_min=1e-1, Tn_max=1e4, y_fac = 2.3,  label='DTUNUC'   )
profile(         'dT_Hep_pbar_LAB_Duperray.txt',     20, resfile='proj_dT_Hep_pbar.png', type='proj', dashes=(20,5,5,5), Tn_min=1e-1, Tn_max=1e4, y_fac = 2.3,   label='Duperray', draw_on_top=True )
profile(         'dT_Hep_pbar_LAB_diMauro.txt',      20, resfile='proj_dT_Hep_pbar.png', type='proj', dashes=(        ), Tn_min=1e-1, Tn_max=1e4, y_fac = 2.3,   label='di Mauro', draw_on_top=True )
if file!='':
    profile(     Hep_file,                           20, resfile='proj_dT_Hep_pbar.png', type='proj', dashes=(        ), Tn_min=1e-1, Tn_max=1e4, y_fac = fac,   label=label,      draw_on_top=True, color='red', alpha=0.5 )


profile(        'dT_Hep_pbar_LAB_DTUNUC_USIN.txt', 500  , resfile='proj_dT_Hep_pbar.png', type='proj', dashes=(3, 3    ), Tn_min=1e-1, Tn_max=9e1, y_fac = 2.3,      draw_on_top=True )
profile(        'dT_Hep_pbar_LAB_Duperray.txt',    500  , resfile='proj_dT_Hep_pbar.png', type='proj', dashes=(20,5,5,5), Tn_min=1e-1, Tn_max=1e4, y_fac = 2.3,      draw_on_top=True )
profile(        'dT_Hep_pbar_LAB_diMauro.txt',     500  , resfile='proj_dT_Hep_pbar.png', type='proj', dashes=(        ), Tn_min=1e-1, Tn_max=1e4, y_fac = 2.3,      draw_on_top=True )
if file!='':
    profile(    Hep_file,                          500  , resfile='proj_dT_Hep_pbar.png', type='proj', dashes=(        ), Tn_min=1e-1, Tn_max=1e4, y_fac = fac,      draw_on_top=True, color='red', alpha=0.5,  )


p_proj.text    (   4e0, 1e-32 , r'$\mathrm{T_{He}/4=                   20\,GeV}$',   color='black',     rotation=-70 )
p_proj.text    (   5e1, 2e-32 , r'$\mathrm{T_{He}/4=                  500\,GeV}$',   color='black',     rotation=-70 )
p_proj.text    (   4e2, 5e-32 , r'$\mathrm{T_{He}/4=  8\,TeV\,\,(\times\,10)  }$',   color='black',     rotation=-60 )



if file!='':
    profile( Hep_file,                          8000  , resfile='proj_dT_Hep_pbar.png', type='proj', dashes=(        ), Tn_min=1e-1, Tn_max=1e4, y_fac = 10*fac, draw_on_top=True, color='red', alpha=0.5 )
profile(     'dT_Hep_pbar_LAB_DTUNUC_USIN.txt', 8000  , resfile='proj_dT_Hep_pbar.png', type='proj', dashes=(3, 3    ), Tn_min=1e-1, Tn_max=9e1, y_fac = 10*2.3, draw_on_top=True)
profile(     'dT_Hep_pbar_LAB_Duperray.txt',    8000  , resfile='proj_dT_Hep_pbar.png', type='proj', dashes=(20,5,5,5), Tn_min=1e-1, Tn_max=1e4, y_fac = 10*2.3, draw_on_top=True)
profile(     'dT_Hep_pbar_LAB_diMauro.txt',     8000  , resfile='proj_dT_Hep_pbar.png', type='proj', dashes=(        ), Tn_min=1e-1, Tn_max=1e4, y_fac = 10*2.3, draw_on_top=True, legend=True, x_label = r'$\mathrm{T_{\bar{p}} \quad [GeV]}$', y_min=5e-34, y_max=1e-29 )




p_prod = profile('dT_Hep_pbar_LAB_DTUNUC_USIN.txt',  0.5, resfile='prod_dT_Hep_pbar.png', type='prod', dashes=(3, 3    ), Tn_min=5, Tn_max=1e4, y_fac = 2.3, label='DTUNUC'  )
profile(         'dT_Hep_pbar_LAB_Duperray.txt',     0.5, resfile='prod_dT_Hep_pbar.png', type='prod', dashes=(20,5,5,5), Tn_min=5, Tn_max=1e4, y_fac = 2.3, label='Duperray' , draw_on_top=True)
profile(         'dT_Hep_pbar_LAB_diMauro.txt',      0.5, resfile='prod_dT_Hep_pbar.png', type='prod', dashes=(        ), Tn_min=5, Tn_max=1e4, y_fac = 2.3, label='di Mauro' , draw_on_top=True)
if file!='':
    profile(     Hep_file,                           0.5, resfile='prod_dT_Hep_pbar.png', type='prod', dashes=(        ), Tn_min=5, Tn_max=1e4, y_fac = fac, label=label,       draw_on_top=True,  color='red', alpha=0.5 )


if file!='':
    profile(Hep_file,                          10  , resfile='prod_dT_Hep_pbar.png', type='prod', dashes=(        ), Tn_min=5, Tn_max=1e4, y_fac = fac, draw_on_top=True, color='red', alpha=0.5  )
profile(    'dT_Hep_pbar_LAB_DTUNUC_USIN.txt', 10  , resfile='prod_dT_Hep_pbar.png', type='prod', dashes=(3, 3    ), Tn_min=5, Tn_max=1e4, y_fac = 2.3, draw_on_top=True)
profile(    'dT_Hep_pbar_LAB_Duperray.txt',    10  , resfile='prod_dT_Hep_pbar.png', type='prod', dashes=(20,5,5,5), Tn_min=5, Tn_max=1e4, y_fac = 2.3, draw_on_top=True)
profile(    'dT_Hep_pbar_LAB_diMauro.txt',     10  , resfile='prod_dT_Hep_pbar.png', type='prod', dashes=(        ), Tn_min=5, Tn_max=1e4, y_fac = 2.3, draw_on_top=True)


p_prod.text    (  9e0, 2e-34 , r'$\mathrm{T_{\bar{p}}=                   0.5\,GeV}$',   color='black', rotation=80 )
p_prod.text    (  2e1, 2e-34 , r'$\mathrm{T_{\bar{p}}=                    10\,GeV}$',   color='black', rotation=80 )
p_prod.text    (  9e1, 2e-34 , r'$\mathrm{T_{\bar{p}}=                    90\,GeV}$',   color='black', rotation=80 )

if file!='':
    profile(Hep_file,                          90  , resfile='prod_dT_Hep_pbar.png', type='prod', dashes=(        ), Tn_min=5, Tn_max=1e4, y_fac = fac, draw_on_top=True,  color='red', alpha=0.5)
profile(    'dT_Hep_pbar_LAB_DTUNUC_USIN.txt', 90  , resfile='prod_dT_Hep_pbar.png', type='prod', dashes=(3, 3    ), Tn_min=5, Tn_max=1e4, y_fac = 2.3, draw_on_top=True)
profile(    'dT_Hep_pbar_LAB_diMauro.txt',     90  , resfile='prod_dT_Hep_pbar.png', type='prod', dashes=(        ), Tn_min=5, Tn_max=1e4, y_fac = 2.3, draw_on_top=True)
profile(    'dT_Hep_pbar_LAB_Duperray.txt',    90  , resfile='prod_dT_Hep_pbar.png', type='prod', dashes=(20,5,5,5), Tn_min=5, Tn_max=1e4, y_fac = 2.3, draw_on_top=True, legend=True, x_label = r'$\mathrm{T_{He}/4 \quad [GeV/4]}$' )

#
#   p He
#


p_proj = profile('dT_pHe_pbar_LAB_DTUNUC_USIN.txt',  20, resfile='proj_dT_pHe_pbar.png', type='proj', dashes=(3, 3    ), Tn_min=1e-1, Tn_max=1e4, y_fac = 2.3, label='DTUNUC'  )
profile(         'dT_pHe_pbar_LAB_Duperray.txt',     20, resfile='proj_dT_pHe_pbar.png', type='proj', dashes=(20,5,5,5), Tn_min=1e-1, Tn_max=1e4, y_fac = 2.3, label='Duperray',  draw_on_top=True)
profile(         'dT_pHe_pbar_LAB_diMauro.txt',      20, resfile='proj_dT_pHe_pbar.png', type='proj', dashes=(        ), Tn_min=1e-1, Tn_max=1e4, y_fac = 2.3, label='di Mauro', draw_on_top=True)
if file!='':
    profile(     pHe_file,                           20, resfile='proj_dT_pHe_pbar.png', type='proj', dashes=(        ), Tn_min=1e-1, Tn_max=1e4, y_fac = fac, label=label,       draw_on_top=True, color='red', alpha=0.5 )

if file!='':
    profile(     pHe_file,                 500  , resfile='proj_dT_pHe_pbar.png', type='proj', dashes=(        ), Tn_min=1e-1, Tn_max=1e4, y_fac = fac,      draw_on_top=True, color='red', alpha=0.5  )
profile('dT_pHe_pbar_LAB_DTUNUC_USIN.txt', 500  , resfile='proj_dT_pHe_pbar.png', type='proj', dashes=(3, 3    ), Tn_min=1e-1, Tn_max=9e1, y_fac = 2.3,      draw_on_top=True )
profile('dT_pHe_pbar_LAB_Duperray.txt',    500  , resfile='proj_dT_pHe_pbar.png', type='proj', dashes=(20,5,5,5), Tn_min=1e-1, Tn_max=1e4, y_fac = 2.3,      draw_on_top=True )
profile('dT_pHe_pbar_LAB_diMauro.txt',     500  , resfile='proj_dT_pHe_pbar.png', type='proj', dashes=(        ), Tn_min=1e-1, Tn_max=1e4, y_fac = 2.3,      draw_on_top=True )

p_proj.text    (   4e0, 1e-32 , r'$\mathrm{T_{p}=                   20\,GeV}$',   color='black',     rotation=-70 )
p_proj.text    (   5e1, 2e-32 , r'$\mathrm{T_{p}=                  500\,GeV}$',   color='black',     rotation=-70 )
p_proj.text    (   4e2, 5e-32 , r'$\mathrm{T_{p}=  8\,TeV\,\,(\times\,10)  }$',   color='black',     rotation=-60 )


if file!='':
    profile(     pHe_file,                 8000  , resfile='proj_dT_pHe_pbar.png', type='proj', dashes=(        ), Tn_min=1e-1, Tn_max=1e4, y_fac = 10*fac, draw_on_top=True, color='red', alpha=0.5 )
profile('dT_pHe_pbar_LAB_DTUNUC_USIN.txt', 8000  , resfile='proj_dT_pHe_pbar.png', type='proj', dashes=(3, 3    ), Tn_min=1e-1, Tn_max=9e1, y_fac = 10*2.3, draw_on_top=True)
profile('dT_pHe_pbar_LAB_Duperray.txt',    8000  , resfile='proj_dT_pHe_pbar.png', type='proj', dashes=(20,5,5,5), Tn_min=1e-1, Tn_max=1e4, y_fac = 10*2.3, draw_on_top=True)
profile('dT_pHe_pbar_LAB_diMauro.txt',     8000  , resfile='proj_dT_pHe_pbar.png', type='proj', dashes=(        ), Tn_min=1e-1, Tn_max=1e4, y_fac = 10*2.3, draw_on_top=True, legend=True, x_label = r'$\mathrm{T_{\bar{p}} \quad [GeV]}$', y_min=5e-34, y_max=1e-29 )


p_prod = profile('dT_pHe_pbar_LAB_DTUNUC_USIN.txt',  0.5, resfile='prod_dT_pHe_pbar.png', type='prod', dashes=(3, 3    ), Tn_min=5, Tn_max=1e4, y_fac = 2.3, label='DTUNUC'  )
profile(         'dT_pHe_pbar_LAB_Duperray.txt',     0.5, resfile='prod_dT_pHe_pbar.png', type='prod', dashes=(20,5,5,5), Tn_min=5, Tn_max=1e4, y_fac = 2.3, label='Duperray' , draw_on_top=True)
profile(         'dT_pHe_pbar_LAB_diMauro.txt',      0.5, resfile='prod_dT_pHe_pbar.png', type='prod', dashes=(        ), Tn_min=5, Tn_max=1e4, y_fac = 2.3, label='di Mauro' , draw_on_top=True)
if file!='':
    profile(      pHe_file,                          0.5, resfile='prod_dT_pHe_pbar.png', type='prod', dashes=(        ), Tn_min=5, Tn_max=1e4, y_fac = fac, label=label      , draw_on_top=True, color='red', alpha=0.5 )


if file!='':
    profile(     pHe_file,                 10  , resfile='prod_dT_pHe_pbar.png', type='prod', dashes=(        ), Tn_min=5, Tn_max=1e4, y_fac = fac, draw_on_top=True, color='red', alpha=0.5 )
profile('dT_pHe_pbar_LAB_DTUNUC_USIN.txt', 10  , resfile='prod_dT_pHe_pbar.png', type='prod', dashes=(3, 3    ), Tn_min=5, Tn_max=1e4, y_fac = 2.3, draw_on_top=True)
profile('dT_pHe_pbar_LAB_Duperray.txt',    10  , resfile='prod_dT_pHe_pbar.png', type='prod', dashes=(20,5,5,5), Tn_min=5, Tn_max=1e4, y_fac = 2.3, draw_on_top=True)
profile('dT_pHe_pbar_LAB_diMauro.txt',     10  , resfile='prod_dT_pHe_pbar.png', type='prod', dashes=(        ), Tn_min=5, Tn_max=1e4, y_fac = 2.3, draw_on_top=True)



p_prod.text    (  9e0, 2e-34 , r'$\mathrm{T_{\bar{p}}=                   0.5\,GeV}$',   color='black', rotation=80 )
p_prod.text    (  2e1, 2e-34 , r'$\mathrm{T_{\bar{p}}=                    10\,GeV}$',   color='black', rotation=80 )
p_prod.text    (  9e1, 2e-34 , r'$\mathrm{T_{\bar{p}}=                    90\,GeV}$',   color='black', rotation=80 )

if file!='':
    profile(     pHe_file,                 90  , resfile='prod_dT_pHe_pbar.png', type='prod', dashes=(        ), Tn_min=5, Tn_max=1e4, y_fac = fac, draw_on_top=True, color='red', alpha=0.5 )
profile('dT_pHe_pbar_LAB_DTUNUC_USIN.txt', 90  , resfile='prod_dT_pHe_pbar.png', type='prod', dashes=(3, 3    ), Tn_min=5, Tn_max=1e4, y_fac = 2.3, draw_on_top=True)
profile('dT_pHe_pbar_LAB_diMauro.txt',     90  , resfile='prod_dT_pHe_pbar.png', type='prod', dashes=(        ), Tn_min=5, Tn_max=1e4, y_fac = 2.3, draw_on_top=True)
profile('dT_pHe_pbar_LAB_Duperray.txt',    90  , resfile='prod_dT_pHe_pbar.png', type='prod', dashes=(20,5,5,5), Tn_min=5, Tn_max=1e4, y_fac = 2.3, draw_on_top=True, legend=True, x_label = r'$\mathrm{T_{p} \quad [GeV]}$' )

#
#   p p
#

p_prod = profile('dT_pp_pbar_LAB_diMauro.txt',      0.5, resfile='prod_dT_pp_pbar.png', type='prod', dashes=(        ), Tn_min=5, Tn_max=1e5, y_fac = 50*2.3, label='di Mauro' )
profile(         'dT_pp_pbar_LAB_TanNg.txt',        0.5, resfile='prod_dT_pp_pbar.png', type='prod', dashes=(15,5    ), Tn_min=5, Tn_max=1e5, y_fac = 50*2.3, label='Tan & Ng' , draw_on_top=True)
profile(         'dT_pp_pbar_LAB_DTUNUC_USIN.txt',  0.5, resfile='prod_dT_pp_pbar.png', type='prod', dashes=(3, 3    ), Tn_min=5, Tn_max=1e4, y_fac = 50*2.3, label='DTUNUC'   , draw_on_top=True)
profile(         'dT_pp_pbar_LAB_Duperray.txt',     0.5, resfile='prod_dT_pp_pbar.png', type='prod', dashes=(20,5,5,5), Tn_min=5, Tn_max=1e4, y_fac = 50*2.3, label='Duperray' , draw_on_top=True)
if file!='':
    profile(     pp_file,                           0.5, resfile='prod_dT_pp_pbar.png', type='prod', dashes=(        ), Tn_min=5, Tn_max=1e4, y_fac = 50*fac, label=label ,      draw_on_top=True, color='red', alpha=0.5)

if file!='':
    profile(     pp_file,                 20  , resfile='prod_dT_pp_pbar.png', type='prod', dashes=(        ), Tn_min=5, Tn_max=1e5, y_fac = fac, draw_on_top=True, color='red', alpha=0.5 )
profile('dT_pp_pbar_LAB_diMauro.txt',     20  , resfile='prod_dT_pp_pbar.png', type='prod', dashes=(        ), Tn_min=5, Tn_max=1e5, y_fac = 2.3, draw_on_top=True)
profile('dT_pp_pbar_LAB_TanNg.txt',       20  , resfile='prod_dT_pp_pbar.png', type='prod', dashes=(15,5    ), Tn_min=5, Tn_max=1e5, y_fac = 2.3, draw_on_top=True)
profile('dT_pp_pbar_LAB_DTUNUC_USIN.txt', 20  , resfile='prod_dT_pp_pbar.png', type='prod', dashes=(3, 3    ), Tn_min=5, Tn_max=1e4, y_fac = 2.3, draw_on_top=True)
profile('dT_pp_pbar_LAB_Duperray.txt',    20  , resfile='prod_dT_pp_pbar.png', type='prod', dashes=(20,5,5,5), Tn_min=5, Tn_max=1e4, y_fac = 2.3, draw_on_top=True)


p_prod.text    (1.2e1, 2e-33 , r'$\mathrm{T_{\bar{p}}= 0.5\,GeV\,\,(\times\, 50) }$',   color='black', rotation=85 )
p_prod.text    (0.5e2, 2e-34 , r'$\mathrm{T_{\bar{p}}=                    20\,GeV}$',   color='black', rotation=80 )
p_prod.text    (  7e2, 2e-34 , r'$\mathrm{T_{\bar{p}}=                   0.5\,TeV}$',   color='black', rotation=75 )

if file!='':
    profile(     pp_file,                500  , resfile='prod_dT_pp_pbar.png', type='prod', dashes=(        ), Tn_min=5, Tn_max=1e5, y_fac = fac, draw_on_top=True, color='red', alpha=0.5 )
profile('dT_pp_pbar_LAB_TanNg.txt',      500  , resfile='prod_dT_pp_pbar.png', type='prod', dashes=(15,5    ), Tn_min=5, Tn_max=1e5, y_fac = 2.3, draw_on_top=True)
profile('dT_pp_pbar_LAB_DTUNUC_USIN.txt',500  , resfile='prod_dT_pp_pbar.png', type='prod', dashes=(3, 3    ), Tn_min=5, Tn_max=1e4, y_fac = 2.3, draw_on_top=True)
profile('dT_pp_pbar_LAB_Duperray.txt',   500  , resfile='prod_dT_pp_pbar.png', type='prod', dashes=(20,5,5,5), Tn_min=5, Tn_max=1e4, y_fac = 2.3, draw_on_top=True)
profile('dT_pp_pbar_LAB_diMauro.txt',    500  , resfile='prod_dT_pp_pbar.png', type='prod', dashes=(        ), Tn_min=5, Tn_max=1e5, y_fac = 2.3, draw_on_top=True, legend=True, x_label = r'$\mathrm{T_{p} \quad [GeV]}$' )


p_proj = profile('dT_pp_pbar_LAB_diMauro.txt',       20, resfile='proj_dT_pp_pbar.png', type='proj', dashes=(        ), y_fac = 2.3, label='di Mauro' )
profile(         'dT_pp_pbar_LAB_TanNg.txt',         20, resfile='proj_dT_pp_pbar.png', type='proj', dashes=(15,5    ), y_fac = 2.3, label='Tan & Ng' , draw_on_top=True)
profile(         'dT_pp_pbar_LAB_DTUNUC_USIN.txt',   20, resfile='proj_dT_pp_pbar.png', type='proj', dashes=(3, 3    ), y_fac = 2.3, label='DTUNUC'   , draw_on_top=True)
profile(         'dT_pp_pbar_LAB_Duperray.txt',      20, resfile='proj_dT_pp_pbar.png', type='proj', dashes=(20,5,5,5), y_fac = 2.3, label='Duperray' , draw_on_top=True)
if file!='':
    profile(     pp_file,                            20, resfile='proj_dT_pp_pbar.png', type='proj', dashes=(        ), y_fac = fac, label=label ,      draw_on_top=True, color='red', alpha=0.5)


if file!='':
    profile(     pp_file,                  500, resfile='proj_dT_pp_pbar.png', type='proj', dashes=(        ), y_fac = fac,                    draw_on_top=True, color='red', alpha=0.5)
profile('dT_pp_pbar_LAB_TanNg.txt',        500, resfile='proj_dT_pp_pbar.png', type='proj', dashes=(15,5    ), y_fac = 2.3,                    draw_on_top=True)
profile('dT_pp_pbar_LAB_DTUNUC_USIN.txt',  500, resfile='proj_dT_pp_pbar.png', type='proj', dashes=(3, 3    ), y_fac = 2.3,                    draw_on_top=True, Tn_max=0.9e2 )
profile('dT_pp_pbar_LAB_Duperray.txt',     500, resfile='proj_dT_pp_pbar.png', type='proj', dashes=(20,5,5,5), y_fac = 2.3,                    draw_on_top=True )
profile('dT_pp_pbar_LAB_diMauro.txt',      500, resfile='proj_dT_pp_pbar.png', type='proj', dashes=(        ), y_fac = 2.3,                    draw_on_top=True )

p_proj.text    ( 0.8e1, 2e-34 , r'$\mathrm{T_{p}=                   20\,GeV}$',   color='black',     rotation=-80 )
p_proj.text    (   2e2, 2e-34 , r'$\mathrm{T_{p}=                  500\,GeV}$',   color='black',     rotation=-80 )
p_proj.text    (   1e3, 2e-33 , r'$\mathrm{T_{p}=  8\,TeV\,\,(\times\,5)   }$',   color='black',     rotation=-70 )

if file!='':
    profile(     pp_file,                 8000, resfile='proj_dT_pp_pbar.png', type='proj', dashes=(        ), y_fac = 5*fac,                    draw_on_top=True, color='red', alpha=0.5 )
profile('dT_pp_pbar_LAB_TanNg.txt',       8000, resfile='proj_dT_pp_pbar.png', type='proj', dashes=(15,5    ), y_fac = 5*2.3,                    draw_on_top=True)
profile('dT_pp_pbar_LAB_DTUNUC_USIN.txt', 8000, resfile='proj_dT_pp_pbar.png', type='proj', dashes=(3, 3    ), y_fac = 5*2.3,                    draw_on_top=True, Tn_max=0.9e2 )
profile('dT_pp_pbar_LAB_Duperray.txt',    8000, resfile='proj_dT_pp_pbar.png', type='proj', dashes=(20,5,5,5), y_fac = 5*2.3,                    draw_on_top=True  )
profile('dT_pp_pbar_LAB_diMauro.txt',     8000, resfile='proj_dT_pp_pbar.png', type='proj', dashes=(        ), y_fac = 5*2.3,                    draw_on_top=True, legend=True, x_label = r'$\mathrm{T_{\bar{p}} \quad [GeV]}$' )


if file!='':
    plot2D(file,            resfile='dT_pp_pbar_LAB_'  + file + '.png'     )
    plot2D('dT_pp_pbar_LAB_diMauro12.txt',            file2=file,       resfile='ratio_pp_diMauro12_'+label+'.png',  Zlabel = label+'/diMauro12' )
    sys.exit(0)



plot2D('dT_pp_pbar_LAB_diMauro.txt',            resfile='dT_pp_pbar_LAB_diMauro.png'     )
plot2D('dT_pp_pbar_LAB_diMauro12.txt',          resfile='dT_pp_pbar_LAB_diMauro12.png'   )
plot2D('dT_pp_pbar_LAB_TanNg.txt',              resfile='dT_pp_pbar_LAB_TanNg.png'       )

plot2D('dT_pp_pbar_LAB_KR_PPFRAG.txt',          resfile='dT_pp_pbar_LAB_KR.png'          )
plot2D('dT_Hep_pbar_LAB_KR_PPFRAG.txt',         resfile='dT_Hep_pbar_LAB_KR.png'         )
plot2D('dT_pHe_pbar_LAB_KR_PPFRAG.txt',         resfile='dT_pHe_pbar_LAB_KR.png'         )

plot2D('dT_pp_pbar_LAB_DTUNUC_USIN.txt',        resfile='dT_pp_pbar_LAB_DTUNUC.png',        Tmax_proj=1e4, Tmax_prod=1e2    )
plot2D('dT_Hep_pbar_LAB_DTUNUC_USIN.txt',       resfile='dT_Hep_pbar_LAB_DTUNUC.png',       Tmax_proj=1e4, Tmax_prod=1e2    )
plot2D('dT_pHe_pbar_LAB_DTUNUC_USIN.txt',       resfile='dT_pHe_pbar_LAB_DTUNUC.png',       Tmax_proj=1e4, Tmax_prod=1e2    )
plot2D('dT_HeHe_pbar_LAB_DTUNUC_USIN.txt',      resfile='dT_HeHe_pbar_LAB_DTUNUC.png',      Tmax_proj=1e4, Tmax_prod=1e2    )


plot2D('dT_pp_pbar_LAB_Duperray.txt',           resfile='dT_pp_pbar_LAB_Duperray.png',      Tmax_proj=1e4, Tmax_prod=1e2    )
plot2D('dT_Hep_pbar_LAB_Duperray.txt',          resfile='dT_Hep_pbar_LAB_Duperray.png',     Tmax_proj=1e4, Tmax_prod=1e2    )
plot2D('dT_pHe_pbar_LAB_Duperray.txt',          resfile='dT_pHe_pbar_LAB_Duperray.png',     Tmax_proj=1e4, Tmax_prod=1e2    )
#plot2D('dT_HeHe_pbar_LAB_Duperray_USIN.txt',     resfile='dT_HeHe_pbar_LAB_Duperray.png',    Tmax_proj=1e4, Tmax_prod=1e2    )


plot2D('dT_pp_pbar_LAB_diMauro.txt',            file2='dT_pp_pbar_LAB_diMauro12.txt',       resfile='ratio_pp_diMauro12_diMauro.png',  Zlabel = 'diMauro12/diMauro13' )
plot2D('dT_pp_pbar_LAB_diMauro.txt',            file2='dT_pp_pbar_LAB_TanNg.txt',           resfile='ratio_pp_TanNg_diMauro.png',      Zlabel = 'Tan&Ng/diMauro'      )
plot2D('dT_pp_pbar_LAB_diMauro.txt',            file2='dT_pp_pbar_LAB_Duperray.txt',        resfile='ratio_pp_Duperray_diMauro.png',   Zlabel = 'Duperray/diMauro',   )
plot2D('dT_pp_pbar_LAB_diMauro.txt',            file2='dT_pp_pbar_LAB_DTUNUC_USIN.txt',     resfile='ratio_pp_DTUNUC_diMauro.png',     Zlabel = 'DTUNUC/diMauro',       Tmax_proj=1e4,  Tmax_prod=1e2         )
plot2D('dT_pp_pbar_LAB_diMauro.txt',            file2='dT_pp_pbar_LAB_KR_PPFRAG.txt',       resfile='ratio_pp_KR_diMauro.png',         Zlabel = 'KR/diMauro',           zscale='log',   Zmin=1e-1, Zmax=1e1   )


plot2D('dT_pp_pbar_LAB_DTUNUC_USIN.txt',        file2='dT_pp_pbar_LAB_Duperray.txt',        resfile='ratio_pp_DTUNUC_Duperray.png',     Zlabel='Duperray/DTUNUC (pp)',  Tmax_proj=1e4,  Tmax_prod=1e2   )
plot2D('dT_Hep_pbar_LAB_DTUNUC_USIN.txt',       file2='dT_Hep_pbar_LAB_Duperray.txt',       resfile='ratio_Hep_DTUNUC_Duperray.png',    Zlabel='Duperray/DTUNUC (Hep)', Tmax_proj=1e4,  Tmax_prod=1e2   )
plot2D('dT_pHe_pbar_LAB_DTUNUC_USIN.txt',       file2='dT_pHe_pbar_LAB_Duperray.txt',       resfile='ratio_pHe_DTUNUC_Duperray.png',    Zlabel='Duperray/DTUNUC (pHe)', Tmax_proj=1e4,  Tmax_prod=1e2   )
#plot2D('dT_HeHe_pbar_LAB_DTUNUC_USIN.txt',      file2='dT_HeHe_pbar_LAB_Duperray_USIN.txt',  resfile='ratio_HeHe_DTUNUC_Duperray.png',   Zlabel='Duperray/DTUNUC (HeHe)',Tmax_proj=1e4,  Tmax_prod=1e2 )


plot2D('dT_Hep_pbar_LAB_DTUNUC_USIN.txt',       file2='dT_pHe_pbar_LAB_DTUNUC_USIN.txt',    resfile='ratio_pHe_to_Hep_DTUNUC.png',      Zlabel='pHe/Hep (DTUNUC)',      Tmax_proj=1e4,  Tmax_prod=1e2 )






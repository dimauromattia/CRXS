#! /usr/bin/env python

import sys
import glob
import numpy as np
import matplotlib.pyplot as plt

sys.path.append("/usr/local/software/antiDeuteronM/scripts")
from PlotFunctions import plot2D
from PlotFunctions import profile

import argparse
from PlotFunctions import plot2D


parser = argparse.ArgumentParser(description='Michael Korsmeier - Plotting script for source terms.')
parser.add_argument('--file_end', help='Filename of table, e.g. KR for dT_Hep_pbar_LAB_KR.txt', action='store', dest='file', type=str, default='')
parser.add_argument('--label',    help='Label.', action='store',   dest='label',  type=str,   default='')
parser.add_argument('--factor',   help='Factor .', action='store', dest='factor', type=float, default=2.3)
args = parser.parse_args()
file   	    = args.file
label 	    = args.label
fac         = args.factor

pp_file  = 'dT_pp_pbar_LAB_'   + file + '.txt'
Hep_file = 'dT_Hep_pbar_LAB_'  + file + '.txt'
pHe_file = 'dT_pHe_pbar_LAB_'  + file + '.txt'

HeHe_file= 'dT_HeHe_pbar_LAB_' + file + '.txt'

c0='black'
c1='#8A0808'
c2='#0B610B'
c3='#0B0B61'


plot2D('dT_Hep_pbar_LAB_Winkler.txt',         resfile='dT_Hep_pbar_LAB_Winkler.pdf'         )
plot2D('dT_pHe_pbar_LAB_Winkler.txt',         resfile='dT_pHe_pbar_LAB_Winkler.pdf'         )



plot2D('dT_pHe_pbar_LAB_diMauro12.txt',            file2='dT_Hep_pbar_LAB_diMauro12.txt',       resfile='ratio_pHe_Hep_diMauro.pdf',  Zlabel = 'Hep/pHe', Zmin=0.8, Zmax=1.2 )


plot2D('dT_pp_pbar_LAB_diMauro12.txt',          resfile='paper_dT_pp_pbar_LAB_diMauro12.pdf'   , x_label=r'$\mathrm{T_{p}\quad [GeV]}$'  )
plot2D('dT_Hep_pbar_LAB_diMauro12.txt',          resfile='paper_dT_Hep_pbar_LAB_diMauro12.pdf'   , x_label=r'$\mathrm{T_{He}/4\quad [GeV]}$'  )
plot2D('dT_pHe_pbar_LAB_diMauro12.txt',          resfile='paper_dT_pHe_pbar_LAB_diMauro12.pdf'   , x_label=r'$\mathrm{T_{p}\quad [GeV]}$'  )




#
#   p p
#

p_prod = profile('dT_pp_pbar_LAB_diMauro12.txt',    0.5, resfile='paper_prod_dT_pp_pbar.pdf', type='prod', color=c0, dashes=(        ), Tn_min=5, Tn_max=1e5, y_fac = 40*2.3, label='di Mauro' )
profile('dT_pp_pbar_LAB_WinklerWithHypWitNbar.txt', 0.5, resfile='paper_prod_dT_pp_pbar.pdf', type='prod', color=c0, dashes=(15,5,3,3,3,5), Tn_min=5, Tn_max=1e5, y_fac = 40*1.0, label='Winkler', draw_on_top=True)

profile(         'dT_pp_pbar_LAB_KR_PPFRAG.txt',    0.5, resfile='paper_prod_dT_pp_pbar.pdf', type='prod', color=c0, dashes=(3, 3    ), Tn_min=5, Tn_max=1e5, y_fac = 40*2.3, label='KMO'   , draw_on_top=True)
profile(         'dT_pp_pbar_LAB_Duperray.txt',     0.5, resfile='paper_prod_dT_pp_pbar.pdf', type='prod', color=c0, dashes=(20,5,5,5), Tn_min=5, Tn_max=1e5, y_fac = 40*2.3, label='Duperray' , draw_on_top=True)
profile(         'dT_pp_pbar_LAB_TanNg.txt',        0.5, resfile='paper_prod_dT_pp_pbar.pdf', type='prod', color=c0, dashes=(12,3    ), Tn_min=5, Tn_max=1e5, y_fac = 40*2.3, label='Tan & Ng' , draw_on_top=True)

if file!='' and len(pp):
    profile(     pp_file,                           0.5, resfile='paper_prod_dT_pp_pbar.pdf', type='prod', dashes=(        ), Tn_min=5, Tn_max=1e5, y_fac = 40*fac, label=label ,      draw_on_top=True, color='red', alpha=0.5)
profile(         'dT_pp_pbar_LAB_diMauro12.txt',    0.5, resfile='paper_prod_dT_pp_pbar.pdf', type='prod', color=c1, dashes=(        ), Tn_min=5, Tn_max=1e5, y_fac = 40*2.3, draw_on_top=True )
profile('dT_pp_pbar_LAB_WinklerWithHypWitNbar.txt', 0.5, resfile='paper_prod_dT_pp_pbar.pdf', type='prod', color=c1, dashes=(15,5,3,3,3,5), Tn_min=5, Tn_max=1e5, y_fac =40*1.0, draw_on_top=True)
profile(         'dT_pp_pbar_LAB_TanNg.txt',        0.5, resfile='paper_prod_dT_pp_pbar.pdf', type='prod', color=c1, dashes=(12,3    ), Tn_min=5, Tn_max=1e5, y_fac = 40*2.3, draw_on_top=True)
profile(         'dT_pp_pbar_LAB_KR_PPFRAG.txt',    0.5, resfile='paper_prod_dT_pp_pbar.pdf', type='prod', color=c1, dashes=(3, 3    ), Tn_min=5, Tn_max=1e5, y_fac = 40*2.3, draw_on_top=True)
profile(         'dT_pp_pbar_LAB_Duperray.txt',     0.5, resfile='paper_prod_dT_pp_pbar.pdf', type='prod', color=c1, dashes=(20,5,5,5), Tn_min=5, Tn_max=1e5, y_fac = 40*2.3, draw_on_top=True)



if file!='' and len(pp):
    profile(     pp_file,                 20  , resfile='paper_prod_dT_pp_pbar.pdf', type='prod', dashes=(        ), Tn_min=5, Tn_max=1e5, y_fac = fac, draw_on_top=True, color='red', alpha=0.5 )
profile('dT_pp_pbar_LAB_diMauro12.txt',   20  , resfile='paper_prod_dT_pp_pbar.pdf', type='prod', color=c2, dashes=(        ), Tn_min=5, Tn_max=1e5, y_fac = 2.3, draw_on_top=True)
profile('dT_pp_pbar_LAB_WinklerWithHypWitNbar.txt',20, resfile='paper_prod_dT_pp_pbar.pdf', type='prod', color=c2, dashes=(15,5,3,3,3,5), Tn_min=5, Tn_max=1e5, y_fac = 1.0, draw_on_top=True)
profile('dT_pp_pbar_LAB_TanNg.txt',       20  , resfile='paper_prod_dT_pp_pbar.pdf', type='prod', color=c2, dashes=(12,3    ), Tn_min=5, Tn_max=1e5, y_fac = 2.3, draw_on_top=True)
profile('dT_pp_pbar_LAB_KR_PPFRAG.txt',   20  , resfile='paper_prod_dT_pp_pbar.pdf', type='prod', color=c2, dashes=(3, 3    ), Tn_min=5, Tn_max=1e5, y_fac = 2.3, draw_on_top=True)
profile('dT_pp_pbar_LAB_Duperray.txt',    20  , resfile='paper_prod_dT_pp_pbar.pdf', type='prod', color=c2, dashes=(20,5,5,5), Tn_min=5, Tn_max=1e5, y_fac = 2.3, draw_on_top=True)


p_prod.text    (1.2e1, 2e-33 , r'$\mathrm{T_{\bar{p}}= 0.5\,GeV\,\,(\times\, 40) }$',   color=c1, rotation=85 )
p_prod.text    (0.5e2, 2e-34 , r'$\mathrm{T_{\bar{p}}=                    20\,GeV}$',   color=c2, rotation=80 )
p_prod.text    (  7e2, 2e-34 , r'$\mathrm{T_{\bar{p}}=                   0.5\,TeV}$',   color=c3, rotation=75 )

if file!='' and len(pp):
    profile(     pp_file,                500  , resfile='paper_prod_dT_pp_pbar.pdf', type='prod', dashes=(        ), Tn_min=5, Tn_max=1e5, y_fac = fac, draw_on_top=True, color='red', alpha=0.5 )
profile('dT_pp_pbar_LAB_TanNg.txt',      500  , resfile='paper_prod_dT_pp_pbar.pdf', type='prod', color=c3, dashes=(12,3    ), Tn_min=5, Tn_max=1e5, y_fac = 2.3, draw_on_top=True)
profile('dT_pp_pbar_LAB_WinklerWithHypWitNbar.txt',500, resfile='paper_prod_dT_pp_pbar.pdf', type='prod', color=c3, dashes=(15,5,3,3,3,5), Tn_min=5, Tn_max=1e5, y_fac = 1.0, draw_on_top=True)
profile('dT_pp_pbar_LAB_KR_PPFRAG.txt',  500  , resfile='paper_prod_dT_pp_pbar.pdf', type='prod', color=c3, dashes=(3, 3    ), Tn_min=5, Tn_max=1e5, y_fac = 2.3, draw_on_top=True)
profile('dT_pp_pbar_LAB_Duperray.txt',   500  , resfile='paper_prod_dT_pp_pbar.pdf', type='prod', color=c3, dashes=(20,5,5,5), Tn_min=5, Tn_max=1e5, y_fac = 2.3, draw_on_top=True)
profile('dT_pp_pbar_LAB_diMauro12.txt',  500  , resfile='paper_prod_dT_pp_pbar.pdf', type='prod', color=c3, dashes=(        ), Tn_min=5, Tn_max=1e5, y_fac = 2.3, draw_on_top=True, legend=True, x_label = r'$\mathrm{T_{p} \quad [GeV]}$' )


p_proj = profile('dT_pp_pbar_LAB_diMauro12.txt',     20, resfile='paper_proj_dT_pp_pbar.pdf', type='proj', color=c0, dashes=(        ), y_fac = 2.3, label='di Mauro' )
profile('dT_pp_pbar_LAB_WinklerWithHypWitNbar.txt',  20, resfile='paper_proj_dT_pp_pbar.pdf', type='proj', color=c0, dashes=(15,5,3,3,3,5), y_fac = 1.0, label='Winkler' , draw_on_top=True)

profile(         'dT_pp_pbar_LAB_KR_PPFRAG.txt',     20, resfile='paper_proj_dT_pp_pbar.pdf', type='proj', color=c0, dashes=(3, 3    ), y_fac = 2.3, label='KMO'   , draw_on_top=True)
profile(         'dT_pp_pbar_LAB_Duperray.txt',      20, resfile='paper_proj_dT_pp_pbar.pdf', type='proj', color=c0, dashes=(20,5,5,5), y_fac = 2.3, label='Duperray' , draw_on_top=True)
profile(         'dT_pp_pbar_LAB_TanNg.txt',         20, resfile='paper_proj_dT_pp_pbar.pdf', type='proj', color=c0, dashes=(12,3    ), y_fac = 2.3, label='Tan & Ng' , draw_on_top=True)

if file!='' and len(pp):
    profile(     pp_file,                            20, resfile='paper_proj_dT_pp_pbar.pdf', type='proj', dashes=(        ), y_fac = fac, label=label ,      draw_on_top=True, color='red', alpha=0.5)
profile(         'dT_pp_pbar_LAB_diMauro12.txt',     20, resfile='paper_proj_dT_pp_pbar.pdf', type='proj', color=c1, dashes=(        ), y_fac = 2.3, draw_on_top=True )
profile('dT_pp_pbar_LAB_WinklerWithHypWitNbar.txt',  20, resfile='paper_proj_dT_pp_pbar.pdf', type='proj', color=c1, dashes=(15,5,3,3,3,5), y_fac = 1.0, draw_on_top=True)
profile(         'dT_pp_pbar_LAB_TanNg.txt',         20, resfile='paper_proj_dT_pp_pbar.pdf', type='proj', color=c1, dashes=(12,3    ), y_fac = 2.3, draw_on_top=True )
profile(         'dT_pp_pbar_LAB_KR_PPFRAG.txt',     20, resfile='paper_proj_dT_pp_pbar.pdf', type='proj', color=c1, dashes=(3, 3    ), y_fac = 2.3, draw_on_top=True )
profile(         'dT_pp_pbar_LAB_Duperray.txt',      20, resfile='paper_proj_dT_pp_pbar.pdf', type='proj', color=c1, dashes=(20,5,5,5), y_fac = 2.3, draw_on_top=True )


if file!='' and len(pp):
    profile(     pp_file,                  450, resfile='paper_proj_dT_pp_pbar.pdf', type='proj', dashes=(        ), y_fac = fac,                    draw_on_top=True, color='red', alpha=0.5)
profile('dT_pp_pbar_LAB_TanNg.txt',        450, resfile='paper_proj_dT_pp_pbar.pdf', type='proj', color=c2, dashes=(12,3    ), y_fac = 2.3,                    draw_on_top=True)
profile('dT_pp_pbar_LAB_WinklerWithHypWitNbar.txt',  450, resfile='paper_proj_dT_pp_pbar.pdf', type='proj', color=c2, dashes=(15,5,3,3,3,5), y_fac = 1.0, draw_on_top=True)
profile('dT_pp_pbar_LAB_KR_PPFRAG.txt',    450, resfile='paper_proj_dT_pp_pbar.pdf', type='proj', color=c2, dashes=(3, 3    ), y_fac = 2.3,                    draw_on_top=True, Tn_max=0.9e2 )
profile('dT_pp_pbar_LAB_Duperray.txt',     450, resfile='paper_proj_dT_pp_pbar.pdf', type='proj', color=c2, dashes=(20,5,5,5), y_fac = 2.3,                    draw_on_top=True )
profile('dT_pp_pbar_LAB_diMauro12.txt',    450, resfile='paper_proj_dT_pp_pbar.pdf', type='proj', color=c2, dashes=(        ), y_fac = 2.3,                    draw_on_top=True )

p_proj.text    ( 0.8e1, 2e-34 , r'$\mathrm{T_{p}=                   20\,GeV}$',   color=c1,     rotation=-80 )
p_proj.text    (   2e2, 2e-34 , r'$\mathrm{T_{p}=                  450\,GeV}$',   color=c2,     rotation=-80 )
p_proj.text    (   1e3, 2e-33 , r'$\mathrm{T_{p}=  6.5\,TeV\,\,(\times\,5)   }$',   color=c3,     rotation=-70 )

if file!='' and len(pp):
    profile(     pp_file,                 6500, resfile='paper_proj_dT_pp_pbar.pdf', type='proj', dashes=(        ), y_fac = 5*fac,                    draw_on_top=True, color='red', alpha=0.5 )
profile('dT_pp_pbar_LAB_WinklerWithHypWitNbar.txt',  6500, resfile='paper_proj_dT_pp_pbar.pdf', type='proj', color=c3, dashes=(15,5,3,3,3,5), y_fac = 5.0, draw_on_top=True)
profile('dT_pp_pbar_LAB_TanNg.txt',       6500, resfile='paper_proj_dT_pp_pbar.pdf', type='proj', color=c3, dashes=(12,3    ), y_fac = 5*2.3,                    draw_on_top=True)
profile('dT_pp_pbar_LAB_KR_PPFRAG.txt',   6500, resfile='paper_proj_dT_pp_pbar.pdf', type='proj', color=c3, dashes=(3, 3    ), y_fac = 5*2.3,                    draw_on_top=True, Tn_max=0.9e2 )
profile('dT_pp_pbar_LAB_Duperray.txt',    6500, resfile='paper_proj_dT_pp_pbar.pdf', type='proj', color=c3, dashes=(20,5,5,5), y_fac = 5*2.3,                    draw_on_top=True  )
profile('dT_pp_pbar_LAB_diMauro12.txt',   6500, resfile='paper_proj_dT_pp_pbar.pdf', type='proj', color=c3, dashes=(        ), y_fac = 5*2.3,                    draw_on_top=True, legend=True, x_label = r'$\mathrm{T_{\bar{p}} \quad [GeV]}$' )






#
#   p He
#


p_proj =profile(         'dT_pHe_pbar_LAB_diMauro.txt',      20, resfile='paper_proj_dT_pHe_pbar.pdf', type='proj', color=c0, dashes=(        ), Tn_min=1e-1, Tn_max=1e4, y_fac = 2.3, label='di Mauro')
profile(         'dT_pHe_pbar_LAB_WinklerWithHypWitNbar.txt',      20, resfile='proj_dT_pHe_pbar.pdf', type='proj', color=c0, dashes=(15,5,3,3,3,5), Tn_min=1e-1, Tn_max=1e4, y_fac = 1.0, label='Winkler', draw_on_top=True)
profile(         'dT_pHe_pbar_LAB_KR_PPFRAG.txt',    20, resfile='paper_proj_dT_pHe_pbar.pdf', type='proj', color=c0, dashes=(3, 3    ), Tn_min=1e-1, Tn_max=1e4, y_fac = 2.3, label='KMO', draw_on_top=True  )
profile(         'dT_pHe_pbar_LAB_Duperray.txt',     20, resfile='paper_proj_dT_pHe_pbar.pdf', type='proj', color=c0, dashes=(20,5,5,5), Tn_min=1e-1, Tn_max=1e4, y_fac = 2.3, label='Duperray',  draw_on_top=True)
profile(         'dT_pHe_pbar_LAB_DTUNUC_USIN.txt',  20, resfile='paper_proj_dT_pHe_pbar.pdf', type='proj', color=c0, dashes=(12,3    ), Tn_min=1e-1, Tn_max=1e4, y_fac = 2.3, label='DTUNUC',  draw_on_top=True)


if file!='' and len(pHe):
    profile(     pHe_file,                           20, resfile='paper_proj_dT_pHe_pbar.pdf', type='proj', dashes=(        ), Tn_min=1e-1, Tn_max=1e4, y_fac = fac, label=label,       draw_on_top=True, color='red', alpha=0.5 )
profile(         'dT_pHe_pbar_LAB_KR_PPFRAG.txt',    20, resfile='paper_proj_dT_pHe_pbar.pdf', type='proj', color=c1, dashes=(3, 3    ), Tn_min=1e-1, Tn_max=1e4, y_fac = 2.3, draw_on_top=True)
profile(         'dT_pHe_pbar_LAB_Duperray.txt',     20, resfile='paper_proj_dT_pHe_pbar.pdf', type='proj', color=c1, dashes=(20,5,5,5), Tn_min=1e-1, Tn_max=1e4, y_fac = 2.3, draw_on_top=True)
profile(         'dT_pHe_pbar_LAB_diMauro.txt',      20, resfile='paper_proj_dT_pHe_pbar.pdf', type='proj', color=c1, dashes=(        ), Tn_min=1e-1, Tn_max=1e4, y_fac = 2.3, draw_on_top=True)
profile(         'dT_pHe_pbar_LAB_WinklerWithHypWitNbar.txt',      20, resfile='proj_dT_pHe_pbar.pdf', type='proj', color=c1, dashes=(15,5,3,3,3,5), Tn_min=1e-1, Tn_max=1e4, y_fac = 1.0, draw_on_top=True)
profile(         'dT_pHe_pbar_LAB_DTUNUC_USIN.txt',  20, resfile='paper_proj_dT_pHe_pbar.pdf', type='proj', color=c1, dashes=(12,3    ), Tn_min=1e-1, Tn_max=1e4, y_fac = 2.3,  draw_on_top=True)


if file!='' and len(pHe):
    profile(     pHe_file,                 450  , resfile='paper_proj_dT_pHe_pbar.pdf', type='proj', dashes=(        ), Tn_min=1e-1, Tn_max=1e4, y_fac = fac,      draw_on_top=True, color='red', alpha=0.5  )
profile('dT_pHe_pbar_LAB_KR_PPFRAG.txt',   450  , resfile='paper_proj_dT_pHe_pbar.pdf', type='proj', color=c2, dashes=(3, 3    ), Tn_min=1e-1, Tn_max=1e4, y_fac = 2.3,      draw_on_top=True )
profile('dT_pHe_pbar_LAB_Duperray.txt',    450  , resfile='paper_proj_dT_pHe_pbar.pdf', type='proj', color=c2, dashes=(20,5,5,5), Tn_min=1e-1, Tn_max=1e4, y_fac = 2.3,      draw_on_top=True )
profile('dT_pHe_pbar_LAB_diMauro.txt',     450  , resfile='paper_proj_dT_pHe_pbar.pdf', type='proj', color=c2, dashes=(        ), Tn_min=1e-1, Tn_max=1e4, y_fac = 2.3,      draw_on_top=True )
profile('dT_pHe_pbar_LAB_WinklerWithHypWitNbar.txt',     450, resfile='proj_dT_pHe_pbar.pdf', type='proj', color=c2, dashes=(15,5,3,3,3,5), Tn_min=1e-1, Tn_max=1e4, y_fac = 1.0, draw_on_top=True)
profile('dT_pHe_pbar_LAB_DTUNUC_USIN.txt', 450  , resfile='paper_proj_dT_pHe_pbar.pdf', type='proj', color=c2, dashes=(12,3    ), Tn_min=1e-1, Tn_max=1e4, y_fac = 2.3,  draw_on_top=True)


p_proj.text    (   4e0, 1e-32 , r'$\mathrm{T_{p}=                   20\,GeV}$',   color=c1,     rotation=-70 )
p_proj.text    (   5e1, 2e-32 , r'$\mathrm{T_{p}=                  450\,GeV}$',   color=c2,     rotation=-70 )
p_proj.text    (   4e2, 5e-32 , r'$\mathrm{T_{p}=  6.5\,TeV\,\,(\times\,10)  }$', color=c3,     rotation=-60 )


if file!='' and len(pHe):
    profile(     pHe_file,                 6500  , resfile='paper_proj_dT_pHe_pbar.pdf', type='proj', dashes=(        ), Tn_min=1e-1, Tn_max=1e4, y_fac = 10*fac, draw_on_top=True, color='red', alpha=0.5 )
profile('dT_pHe_pbar_LAB_KR_PPFRAG.txt',   6500  , resfile='paper_proj_dT_pHe_pbar.pdf', type='proj', color=c3, dashes=(3, 3    ), Tn_min=1e-1, Tn_max=1e4, y_fac = 10*2.3, draw_on_top=True)
profile('dT_pHe_pbar_LAB_Duperray.txt',    6500  , resfile='paper_proj_dT_pHe_pbar.pdf', type='proj', color=c3, dashes=(20,5,5,5), Tn_min=1e-1, Tn_max=1e4, y_fac = 10*2.3, draw_on_top=True)
profile('dT_pHe_pbar_LAB_WinklerWithHypWitNbar.txt',      6500, resfile='proj_dT_pHe_pbar.pdf', type='proj', color=c3, dashes=(15,5,3,3,3,5), Tn_min=1e-1, Tn_max=1e4, y_fac = 10*1.0, draw_on_top=True)
profile('dT_pHe_pbar_LAB_diMauro.txt',     6500  , resfile='paper_proj_dT_pHe_pbar.pdf', type='proj', color=c3, dashes=(        ), Tn_min=1e-1, Tn_max=1e4, y_fac = 10*2.3, draw_on_top=True, legend=True, x_label = r'$\mathrm{T_{\bar{p}} \quad [GeV]}$', y_min=5e-34, y_max=1e-29 )



#8033
#8283

sys.exit(0)

#
#   He He
#

HeHe = glob.glob( HeHe_file )
pHe  = glob.glob( pHe_file  )
Hep  = glob.glob( Hep_file  )
pp   = glob.glob( pp_file   )


p_proj = profile('dT_HeHe_pbar_LAB_KR_PPFRAG.txt',  20, resfile='proj_dT_HeHe_pbar.pdf', type='proj', dashes=(3, 3    ), Tn_min=1e-1, Tn_max=1e4, y_fac = 2.3,  label='KMO'   )
#profile(         'dT_HeHe_pbar_LAB_Duperray_USIN.txt',20, resfile='proj_dT_HeHe_pbar.pdf', type='proj', dashes=(20,5,5,5), Tn_min=1e-1, Tn_max=1e4, y_fac = 2.3,   label='Duperray', draw_on_top=True )
profile(         'dT_pp_pbar_LAB_diMauro.txt',        20, resfile='proj_dT_HeHe_pbar.pdf', type='proj', dashes=(        ), Tn_min=1e-1, Tn_max=1e4, y_fac = 16 ,   label='di Mauro', draw_on_top=True )
if file!='' and len(HeHe):
    profile(     HeHe_file,                           20, resfile='proj_dT_HeHe_pbar.pdf', type='proj', dashes=(        ), Tn_min=1e-1, Tn_max=1e4, y_fac = fac,   label=label,      draw_on_top=True, color='red', alpha=0.5 )


profile(        'dT_HeHe_pbar_LAB_KR_PPFRAG.txt',  500  , resfile='proj_dT_HeHe_pbar.pdf', type='proj', dashes=(3, 3    ), Tn_min=1e-1, Tn_max=1e4, y_fac = 2.3,      draw_on_top=True )
#profile(        'dT_HeHe_pbar_LAB_Duperray_USIN.txt',500  , resfile='proj_dT_HeHe_pbar.pdf', type='proj', dashes=(20,5,5,5), Tn_min=1e-1, Tn_max=1e4, y_fac = 2.3,      draw_on_top=True )
profile(        'dT_pp_pbar_LAB_diMauro.txt',        500  , resfile='proj_dT_HeHe_pbar.pdf', type='proj', dashes=(        ), Tn_min=1e-1, Tn_max=1e4, y_fac = 16 ,      draw_on_top=True )
if file!='' and len(HeHe):
    profile(    HeHe_file,                           500  , resfile='proj_dT_HeHe_pbar.pdf', type='proj', dashes=(        ), Tn_min=1e-1, Tn_max=1e4, y_fac = fac,      draw_on_top=True, color='red', alpha=0.5,  )


p_proj.text    (   8e0, 1e-32 , r'$\mathrm{T_{He}/4=                   20\,GeV}$',   color='black',     rotation=-70 )
p_proj.text    (   1e2, 2e-32 , r'$\mathrm{T_{He}/4=                  500\,GeV}$',   color='black',     rotation=-70 )
p_proj.text    (   8e2, 5e-32 , r'$\mathrm{T_{He}/4=  8\,TeV\,\,(\times\,10)  }$',   color='black',     rotation=-60 )



if file!='' and len(HeHe):
    profile( HeHe_file,                           8000  , resfile='proj_dT_HeHe_pbar.pdf', type='proj', dashes=(        ), Tn_min=1e-1, Tn_max=1e4, y_fac = 10*fac, draw_on_top=True, color='red', alpha=0.5 )
profile(     'dT_HeHe_pbar_LAB_KR_PPFRAG.txt',  8000  , resfile='proj_dT_HeHe_pbar.pdf', type='proj', dashes=(3, 3    ), Tn_min=1e-1, Tn_max=1e4, y_fac = 10*2.3, draw_on_top=True)
#profile(     'dT_HeHe_pbar_LAB_Duperray_USIN.txt',8000  , resfile='proj_dT_HeHe_pbar.pdf', type='proj', dashes=(20,5,5,5), Tn_min=1e-1, Tn_max=1e4, y_fac = 10*2.3, draw_on_top=True)
profile(     'dT_pp_pbar_LAB_diMauro.txt',        8000  , resfile='proj_dT_HeHe_pbar.pdf', type='proj', dashes=(        ), Tn_min=1e-1, Tn_max=1e4, y_fac = 10*16 , draw_on_top=True, legend=True, x_label = r'$\mathrm{T_{\bar{p}} \quad [GeV]}$', y_min=5e-34, y_max=1e-29 )




p_prod = profile('dT_HeHe_pbar_LAB_KR_PPFRAG.txt',  0.5, resfile='prod_dT_HeHe_pbar.pdf', type='prod', dashes=(3, 3    ), Tn_min=5, Tn_max=1e4, y_fac = 2.3, label='KMO'  )
#profile(         'dT_HeHe_pbar_LAB_Duperray_USIN.txt', 0.5, resfile='prod_dT_HeHe_pbar.pdf', type='prod', dashes=(20,5,5,5), Tn_min=5, Tn_max=1e4, y_fac = 2.3, label='Duperray' , draw_on_top=True)
profile(         'dT_pp_pbar_LAB_diMauro.txt',        0.5, resfile='prod_dT_HeHe_pbar.pdf', type='prod', dashes=(        ), Tn_min=5, Tn_max=1e4, y_fac = 16 , label='di Mauro' , draw_on_top=True)
if file!='' and len(HeHe):
    profile(     HeHe_file,                           0.5, resfile='prod_dT_HeHe_pbar.pdf', type='prod', dashes=(        ), Tn_min=5, Tn_max=1e4, y_fac = fac, label=label,       draw_on_top=True,  color='red', alpha=0.5 )


if file!='' and len(HeHe):
    profile(HeHe_file,                            10  , resfile='prod_dT_HeHe_pbar.pdf', type='prod', dashes=(        ), Tn_min=5, Tn_max=1e4, y_fac = fac, draw_on_top=True, color='red', alpha=0.5  )
profile(    'dT_HeHe_pbar_LAB_KR_PPFRAG.txt',   10  , resfile='prod_dT_HeHe_pbar.pdf', type='prod', dashes=(3, 3    ), Tn_min=5, Tn_max=1e4, y_fac = 2.3, draw_on_top=True)
#profile(    'dT_HeHe_pbar_LAB_Duperray_USIN.txt', 10  , resfile='prod_dT_HeHe_pbar.pdf', type='prod', dashes=(20,5,5,5), Tn_min=5, Tn_max=1e4, y_fac = 2.3, draw_on_top=True)
profile(    'dT_pp_pbar_LAB_diMauro.txt',         10  , resfile='prod_dT_HeHe_pbar.pdf', type='prod', dashes=(        ), Tn_min=5, Tn_max=1e4, y_fac = 16 , draw_on_top=True)


p_prod.text    (  9e0, 2e-34 , r'$\mathrm{T_{\bar{p}}=                   0.5\,GeV}$',   color='black', rotation=80 )
p_prod.text    (  2e1, 2e-34 , r'$\mathrm{T_{\bar{p}}=                    10\,GeV}$',   color='black', rotation=80 )
p_prod.text    (  9e1, 2e-34 , r'$\mathrm{T_{\bar{p}}=                    90\,GeV}$',   color='black', rotation=80 )

if file!='' and len(HeHe):
    profile(HeHe_file,                           90  , resfile='prod_dT_HeHe_pbar.pdf', type='prod', dashes=(        ), Tn_min=5, Tn_max=1e4, y_fac = fac, draw_on_top=True,  color='red', alpha=0.5)
profile(    'dT_HeHe_pbar_LAB_KR_PPFRAG.txt',  90  , resfile='prod_dT_HeHe_pbar.pdf', type='prod', dashes=(3, 3    ), Tn_min=5, Tn_max=1e4, y_fac = 2.3, draw_on_top=True)
#profile(    'dT_HeHe_pbar_LAB_Duperray_USIN.txt',90  , resfile='prod_dT_HeHe_pbar.pdf', type='prod', dashes=(20,5,5,5), Tn_min=5, Tn_max=1e4, y_fac = 2.3, draw_on_top=True )
profile(    'dT_pp_pbar_LAB_diMauro.txt',        90  , resfile='prod_dT_HeHe_pbar.pdf', type='prod', dashes=(        ), Tn_min=5, Tn_max=1e4, y_fac = 16 , draw_on_top=True, legend=True, x_label = r'$\mathrm{T_{He}/4 \quad [GeV/4]}$')


#
#   He p
#

p_proj = profile('dT_Hep_pbar_LAB_KR_PPFRAG.txt',  20, resfile='proj_dT_Hep_pbar.pdf', type='proj', dashes=(3, 3    ), Tn_min=1e-1, Tn_max=1e4, y_fac = 2.3,  label='KMO'   )
profile(         'dT_Hep_pbar_LAB_Duperray.txt',     20, resfile='proj_dT_Hep_pbar.pdf', type='proj', dashes=(20,5,5,5), Tn_min=1e-1, Tn_max=1e4, y_fac = 2.3,   label='Duperray', draw_on_top=True )
profile(         'dT_Hep_pbar_LAB_diMauro.txt',      20, resfile='proj_dT_Hep_pbar.pdf', type='proj', dashes=(        ), Tn_min=1e-1, Tn_max=1e4, y_fac = 2.3,   label='di Mauro', draw_on_top=True )
if file!='' and len(Hep):
    profile(     Hep_file,                           20, resfile='proj_dT_Hep_pbar.pdf', type='proj', dashes=(        ), Tn_min=1e-1, Tn_max=1e4, y_fac = fac,   label=label,      draw_on_top=True, color='red', alpha=0.5 )


profile(        'dT_Hep_pbar_LAB_KR_PPFRAG.txt', 500  , resfile='proj_dT_Hep_pbar.pdf', type='proj', dashes=(3, 3    ), Tn_min=1e-1, Tn_max=1e4, y_fac = 2.3,      draw_on_top=True )
profile(        'dT_Hep_pbar_LAB_Duperray.txt',    500  , resfile='proj_dT_Hep_pbar.pdf', type='proj', dashes=(20,5,5,5), Tn_min=1e-1, Tn_max=1e4, y_fac = 2.3,      draw_on_top=True )
profile(        'dT_Hep_pbar_LAB_diMauro.txt',     500  , resfile='proj_dT_Hep_pbar.pdf', type='proj', dashes=(        ), Tn_min=1e-1, Tn_max=1e4, y_fac = 2.3,      draw_on_top=True )
if file!='' and len(Hep):
    profile(    Hep_file,                          500  , resfile='proj_dT_Hep_pbar.pdf', type='proj', dashes=(        ), Tn_min=1e-1, Tn_max=1e4, y_fac = fac,      draw_on_top=True, color='red', alpha=0.5,  )


p_proj.text    (   4e0, 1e-32 , r'$\mathrm{T_{He}/4=                   20\,GeV}$',   color='black',     rotation=-70 )
p_proj.text    (   5e1, 2e-32 , r'$\mathrm{T_{He}/4=                  500\,GeV}$',   color='black',     rotation=-70 )
p_proj.text    (   4e2, 5e-32 , r'$\mathrm{T_{He}/4=  8\,TeV\,\,(\times\,10)  }$',   color='black',     rotation=-60 )



if file!='' and len(Hep):
    profile( Hep_file,                          8000  , resfile='proj_dT_Hep_pbar.pdf', type='proj', dashes=(        ), Tn_min=1e-1, Tn_max=1e4, y_fac = 10*fac, draw_on_top=True, color='red', alpha=0.5 )
profile(     'dT_Hep_pbar_LAB_KR_PPFRAG.txt',   8000  , resfile='proj_dT_Hep_pbar.pdf', type='proj', dashes=(3, 3    ), Tn_min=1e-1, Tn_max=1e4, y_fac = 10*2.3, draw_on_top=True)
profile(     'dT_Hep_pbar_LAB_Duperray.txt',    8000  , resfile='proj_dT_Hep_pbar.pdf', type='proj', dashes=(20,5,5,5), Tn_min=1e-1, Tn_max=1e4, y_fac = 10*2.3, draw_on_top=True)
profile(     'dT_Hep_pbar_LAB_diMauro.txt',     8000  , resfile='proj_dT_Hep_pbar.pdf', type='proj', dashes=(        ), Tn_min=1e-1, Tn_max=1e4, y_fac = 10*2.3, draw_on_top=True, legend=True, x_label = r'$\mathrm{T_{\bar{p}} \quad [GeV]}$', y_min=5e-34, y_max=1e-29 )




p_prod = profile('dT_Hep_pbar_LAB_KR_PPFRAG.txt',  0.5, resfile='prod_dT_Hep_pbar.pdf', type='prod', dashes=(3, 3    ), Tn_min=5, Tn_max=1e4, y_fac = 2.3, label='KMO'  )
profile(         'dT_Hep_pbar_LAB_Duperray.txt',     0.5, resfile='prod_dT_Hep_pbar.pdf', type='prod', dashes=(20,5,5,5), Tn_min=5, Tn_max=1e4, y_fac = 2.3, label='Duperray' , draw_on_top=True)
profile(         'dT_Hep_pbar_LAB_diMauro.txt',      0.5, resfile='prod_dT_Hep_pbar.pdf', type='prod', dashes=(        ), Tn_min=5, Tn_max=1e4, y_fac = 2.3, label='di Mauro' , draw_on_top=True)
if file!='' and len(Hep):
    profile(     Hep_file,                           0.5, resfile='prod_dT_Hep_pbar.pdf', type='prod', dashes=(        ), Tn_min=5, Tn_max=1e4, y_fac = fac, label=label,       draw_on_top=True,  color='red', alpha=0.5 )


if file!='' and len(Hep):
    profile(Hep_file,                          10  , resfile='prod_dT_Hep_pbar.pdf', type='prod', dashes=(        ), Tn_min=5, Tn_max=1e4, y_fac = fac, draw_on_top=True, color='red', alpha=0.5  )
profile(    'dT_Hep_pbar_LAB_KR_PPFRAG.txt',   10  , resfile='prod_dT_Hep_pbar.pdf', type='prod', dashes=(3, 3    ), Tn_min=5, Tn_max=1e4, y_fac = 2.3, draw_on_top=True)
profile(    'dT_Hep_pbar_LAB_Duperray.txt',    10  , resfile='prod_dT_Hep_pbar.pdf', type='prod', dashes=(20,5,5,5), Tn_min=5, Tn_max=1e4, y_fac = 2.3, draw_on_top=True)
profile(    'dT_Hep_pbar_LAB_diMauro.txt',     10  , resfile='prod_dT_Hep_pbar.pdf', type='prod', dashes=(        ), Tn_min=5, Tn_max=1e4, y_fac = 2.3, draw_on_top=True)


p_prod.text    (  9e0, 2e-34 , r'$\mathrm{T_{\bar{p}}=                   0.5\,GeV}$',   color='black', rotation=80 )
p_prod.text    (  2e1, 2e-34 , r'$\mathrm{T_{\bar{p}}=                    10\,GeV}$',   color='black', rotation=80 )
p_prod.text    (  9e1, 2e-34 , r'$\mathrm{T_{\bar{p}}=                    90\,GeV}$',   color='black', rotation=80 )

if file!='' and len(Hep):
    profile(Hep_file,                          90  , resfile='prod_dT_Hep_pbar.pdf', type='prod', dashes=(        ), Tn_min=5, Tn_max=1e4, y_fac = fac, draw_on_top=True,  color='red', alpha=0.5)
profile(    'dT_Hep_pbar_LAB_KR_PPFRAG.txt',   90  , resfile='prod_dT_Hep_pbar.pdf', type='prod', dashes=(3, 3    ), Tn_min=5, Tn_max=1e4, y_fac = 2.3, draw_on_top=True)
profile(    'dT_Hep_pbar_LAB_diMauro.txt',     90  , resfile='prod_dT_Hep_pbar.pdf', type='prod', dashes=(        ), Tn_min=5, Tn_max=1e4, y_fac = 2.3, draw_on_top=True)
profile(    'dT_Hep_pbar_LAB_Duperray.txt',    90  , resfile='prod_dT_Hep_pbar.pdf', type='prod', dashes=(20,5,5,5), Tn_min=5, Tn_max=1e4, y_fac = 2.3, draw_on_top=True, legend=True, x_label = r'$\mathrm{T_{He}/4 \quad [GeV/4]}$' )

#
#   p He
#




p_proj = profile('dT_pHe_pbar_LAB_KR_PPFRAG.txt',    20, resfile='proj_dT_pHe_pbar.pdf', type='proj', color=c0, dashes=(3, 3    ), Tn_min=1e-1, Tn_max=1e4, y_fac = 2.3, label='KMO'  )
profile(         'dT_pHe_pbar_LAB_Duperray.txt',     20, resfile='proj_dT_pHe_pbar.pdf', type='proj', color=c0, dashes=(20,5,5,5), Tn_min=1e-1, Tn_max=1e4, y_fac = 2.3, label='Duperray',  draw_on_top=True)
profile(         'dT_pHe_pbar_LAB_diMauro.txt',      20, resfile='proj_dT_pHe_pbar.pdf', type='proj', color=c0, dashes=(        ), Tn_min=1e-1, Tn_max=1e4, y_fac = 2.3, label='di Mauro', draw_on_top=True)
if file!='' and len(pHe):
    profile(     pHe_file,                           20, resfile='proj_dT_pHe_pbar.pdf', type='proj', dashes=(        ), Tn_min=1e-1, Tn_max=1e4, y_fac = fac, label=label,       draw_on_top=True, color='red', alpha=0.5 )
profile(         'dT_pHe_pbar_LAB_KR_PPFRAG.txt',    20, resfile='proj_dT_pHe_pbar.pdf', type='proj', color=c1, dashes=(3, 3    ), Tn_min=1e-1, Tn_max=1e4, y_fac = 2.3, draw_on_top=True)
profile(         'dT_pHe_pbar_LAB_Duperray.txt',     20, resfile='proj_dT_pHe_pbar.pdf', type='proj', color=c1, dashes=(20,5,5,5), Tn_min=1e-1, Tn_max=1e4, y_fac = 2.3, draw_on_top=True)
profile(         'dT_pHe_pbar_LAB_diMauro.txt',      20, resfile='proj_dT_pHe_pbar.pdf', type='proj', color=c1, dashes=(        ), Tn_min=1e-1, Tn_max=1e4, y_fac = 2.3, draw_on_top=True)

if file!='' and len(pHe):
    profile(     pHe_file,                 500  , resfile='proj_dT_pHe_pbar.pdf', type='proj', dashes=(        ), Tn_min=1e-1, Tn_max=1e4, y_fac = fac,      draw_on_top=True, color='red', alpha=0.5  )
profile('dT_pHe_pbar_LAB_KR_PPFRAG.txt',   500  , resfile='proj_dT_pHe_pbar.pdf', type='proj', color=c2, dashes=(3, 3    ), Tn_min=1e-1, Tn_max=1e4, y_fac = 2.3,      draw_on_top=True )
profile('dT_pHe_pbar_LAB_Duperray.txt',    500  , resfile='proj_dT_pHe_pbar.pdf', type='proj', color=c2, dashes=(20,5,5,5), Tn_min=1e-1, Tn_max=1e4, y_fac = 2.3,      draw_on_top=True )
profile('dT_pHe_pbar_LAB_diMauro.txt',     500  , resfile='proj_dT_pHe_pbar.pdf', type='proj', color=c2, dashes=(        ), Tn_min=1e-1, Tn_max=1e4, y_fac = 2.3,      draw_on_top=True )

p_proj.text    (   4e0, 1e-32 , r'$\mathrm{T_{p}=                   20\,GeV}$',   color=c1,     rotation=-70 )
p_proj.text    (   5e1, 2e-32 , r'$\mathrm{T_{p}=                  500\,GeV}$',   color=c2,     rotation=-70 )
p_proj.text    (   4e2, 5e-32 , r'$\mathrm{T_{p}=  8\,TeV\,\,(\times\,10)  }$',   color=c3,     rotation=-60 )


if file!='' and len(pHe):
    profile(     pHe_file,                 8000  , resfile='proj_dT_pHe_pbar.pdf', type='proj', dashes=(        ), Tn_min=1e-1, Tn_max=1e4, y_fac = 10*fac, draw_on_top=True, color='red', alpha=0.5 )
profile('dT_pHe_pbar_LAB_KR_PPFRAG.txt',   8000  , resfile='proj_dT_pHe_pbar.pdf', type='proj', color=c3, dashes=(3, 3    ), Tn_min=1e-1, Tn_max=1e4, y_fac = 10*2.3, draw_on_top=True)
profile('dT_pHe_pbar_LAB_Duperray.txt',    8000  , resfile='proj_dT_pHe_pbar.pdf', type='proj', color=c3, dashes=(20,5,5,5), Tn_min=1e-1, Tn_max=1e4, y_fac = 10*2.3, draw_on_top=True)
profile('dT_pHe_pbar_LAB_diMauro.txt',     8000  , resfile='proj_dT_pHe_pbar.pdf', type='proj', color=c3, dashes=(        ), Tn_min=1e-1, Tn_max=1e4, y_fac = 10*2.3, draw_on_top=True, legend=True, x_label = r'$\mathrm{T_{\bar{p}} \quad [GeV]}$', y_min=5e-34, y_max=1e-29 )


p_prod = profile('dT_pHe_pbar_LAB_KR_PPFRAG.txt',    0.5, resfile='prod_dT_pHe_pbar.pdf', type='prod', color=c0, dashes=(3, 3    ), Tn_min=5, Tn_max=1e4, y_fac = 2.3, label='KMO'  )
profile(         'dT_pHe_pbar_LAB_Duperray.txt',     0.5, resfile='prod_dT_pHe_pbar.pdf', type='prod', color=c0, dashes=(20,5,5,5), Tn_min=5, Tn_max=1e4, y_fac = 2.3, label='Duperray' , draw_on_top=True)
profile(         'dT_pHe_pbar_LAB_diMauro.txt',      0.5, resfile='prod_dT_pHe_pbar.pdf', type='prod', color=c0, dashes=(        ), Tn_min=5, Tn_max=1e4, y_fac = 2.3, label='di Mauro' , draw_on_top=True)
if file!='' and len(pHe):
    profile(      pHe_file,                          0.5, resfile='prod_dT_pHe_pbar.pdf', type='prod', dashes=(        ), Tn_min=5, Tn_max=1e4, y_fac = fac, label=label      , draw_on_top=True, color='red', alpha=0.5 )
profile(         'dT_pHe_pbar_LAB_KR_PPFRAG.txt',    0.5, resfile='prod_dT_pHe_pbar.pdf', type='prod', color=c1, dashes=(3, 3    ), Tn_min=5, Tn_max=1e4, y_fac = 2.3, draw_on_top=True)
profile(         'dT_pHe_pbar_LAB_Duperray.txt',     0.5, resfile='prod_dT_pHe_pbar.pdf', type='prod', color=c1, dashes=(20,5,5,5), Tn_min=5, Tn_max=1e4, y_fac = 2.3, draw_on_top=True)
profile(         'dT_pHe_pbar_LAB_diMauro.txt',      0.5, resfile='prod_dT_pHe_pbar.pdf', type='prod', color=c1, dashes=(        ), Tn_min=5, Tn_max=1e4, y_fac = 2.3, draw_on_top=True)


if file!='' and len(pHe):
    profile(     pHe_file,                 10  , resfile='prod_dT_pHe_pbar.pdf', type='prod', dashes=(        ), Tn_min=5, Tn_max=1e4, y_fac = fac, draw_on_top=True, color='red', alpha=0.5 )
profile('dT_pHe_pbar_LAB_KR_PPFRAG.txt',   10  , resfile='prod_dT_pHe_pbar.pdf', type='prod', color=c2, dashes=(3, 3    ), Tn_min=5, Tn_max=1e4, y_fac = 2.3, draw_on_top=True)
profile('dT_pHe_pbar_LAB_Duperray.txt',    10  , resfile='prod_dT_pHe_pbar.pdf', type='prod', color=c2, dashes=(20,5,5,5), Tn_min=5, Tn_max=1e4, y_fac = 2.3, draw_on_top=True)
profile('dT_pHe_pbar_LAB_diMauro.txt',     10  , resfile='prod_dT_pHe_pbar.pdf', type='prod', color=c2, dashes=(        ), Tn_min=5, Tn_max=1e4, y_fac = 2.3, draw_on_top=True)



p_prod.text    (  9e0, 2e-34 , r'$\mathrm{T_{\bar{p}}=                   0.5\,GeV}$',   color=c1, rotation=80 )
p_prod.text    (  2e1, 2e-34 , r'$\mathrm{T_{\bar{p}}=                    10\,GeV}$',   color=c2, rotation=80 )
p_prod.text    (  9e1, 2e-34 , r'$\mathrm{T_{\bar{p}}=                    90\,GeV}$',   color=c3, rotation=80 )

if file!='' and len(pHe):
    profile(     pHe_file,                 90  , resfile='prod_dT_pHe_pbar.pdf', type='prod', dashes=(        ), Tn_min=5, Tn_max=1e4, y_fac = fac, draw_on_top=True, color='red', alpha=0.5 )
profile('dT_pHe_pbar_LAB_KR_PPFRAG.txt',   90  , resfile='prod_dT_pHe_pbar.pdf', type='prod', color=c3, dashes=(3, 3    ), Tn_min=5, Tn_max=1e4, y_fac = 2.3, draw_on_top=True)
profile('dT_pHe_pbar_LAB_diMauro.txt',     90  , resfile='prod_dT_pHe_pbar.pdf', type='prod', color=c3, dashes=(        ), Tn_min=5, Tn_max=1e4, y_fac = 2.3, draw_on_top=True)
profile('dT_pHe_pbar_LAB_Duperray.txt',    90  , resfile='prod_dT_pHe_pbar.pdf', type='prod', color=c3, dashes=(20,5,5,5), Tn_min=5, Tn_max=1e4, y_fac = 2.3, draw_on_top=True, legend=True, x_label = r'$\mathrm{T_{p} \quad [GeV]}$' )

#
#   p p
#

p_prod = profile('dT_pp_pbar_LAB_diMauro12.txt',    0.5, resfile='prod_dT_pp_pbar.pdf', type='prod', color=c0, dashes=(        ), Tn_min=5, Tn_max=1e5, y_fac = 40*2.3, label='di Mauro' )
profile(         'dT_pp_pbar_LAB_TanNg.txt',        0.5, resfile='prod_dT_pp_pbar.pdf', type='prod', color=c0, dashes=(12,3    ), Tn_min=5, Tn_max=1e5, y_fac = 40*2.3, label='Tan & Ng' , draw_on_top=True)
profile(         'dT_pp_pbar_LAB_KR_PPFRAG.txt',    0.5, resfile='prod_dT_pp_pbar.pdf', type='prod', color=c0, dashes=(3, 3    ), Tn_min=5, Tn_max=1e5, y_fac = 40*2.3, label='KMO'   , draw_on_top=True)
profile(         'dT_pp_pbar_LAB_Duperray.txt',     0.5, resfile='prod_dT_pp_pbar.pdf', type='prod', color=c0, dashes=(20,5,5,5), Tn_min=5, Tn_max=1e5, y_fac = 40*2.3, label='Duperray' , draw_on_top=True)
if file!='' and len(pp):
    profile(     pp_file,                           0.5, resfile='prod_dT_pp_pbar.pdf', type='prod', dashes=(        ), Tn_min=5, Tn_max=1e5, y_fac = 40*fac, label=label ,      draw_on_top=True, color='red', alpha=0.5)
profile(         'dT_pp_pbar_LAB_diMauro12.txt',    0.5, resfile='prod_dT_pp_pbar.pdf', type='prod', color=c1, dashes=(        ), Tn_min=5, Tn_max=1e5, y_fac = 40*2.3, draw_on_top=True )
profile(         'dT_pp_pbar_LAB_TanNg.txt',        0.5, resfile='prod_dT_pp_pbar.pdf', type='prod', color=c1, dashes=(12,3    ), Tn_min=5, Tn_max=1e5, y_fac = 40*2.3, draw_on_top=True)
profile(         'dT_pp_pbar_LAB_KR_PPFRAG.txt',    0.5, resfile='prod_dT_pp_pbar.pdf', type='prod', color=c1, dashes=(3, 3    ), Tn_min=5, Tn_max=1e5, y_fac = 40*2.3, draw_on_top=True)
profile(         'dT_pp_pbar_LAB_Duperray.txt',     0.5, resfile='prod_dT_pp_pbar.pdf', type='prod', color=c1, dashes=(20,5,5,5), Tn_min=5, Tn_max=1e5, y_fac = 40*2.3, draw_on_top=True)



if file!='' and len(pp):
    profile(     pp_file,                 20  , resfile='prod_dT_pp_pbar.pdf', type='prod', dashes=(        ), Tn_min=5, Tn_max=1e5, y_fac = fac, draw_on_top=True, color='red', alpha=0.5 )
profile('dT_pp_pbar_LAB_diMauro12.txt',   20  , resfile='prod_dT_pp_pbar.pdf', type='prod', color=c2, dashes=(        ), Tn_min=5, Tn_max=1e5, y_fac = 2.3, draw_on_top=True)
profile('dT_pp_pbar_LAB_TanNg.txt',       20  , resfile='prod_dT_pp_pbar.pdf', type='prod', color=c2, dashes=(12,3    ), Tn_min=5, Tn_max=1e5, y_fac = 2.3, draw_on_top=True)
profile('dT_pp_pbar_LAB_KR_PPFRAG.txt',   20  , resfile='prod_dT_pp_pbar.pdf', type='prod', color=c2, dashes=(3, 3    ), Tn_min=5, Tn_max=1e5, y_fac = 2.3, draw_on_top=True)
profile('dT_pp_pbar_LAB_Duperray.txt',    20  , resfile='prod_dT_pp_pbar.pdf', type='prod', color=c2, dashes=(20,5,5,5), Tn_min=5, Tn_max=1e5, y_fac = 2.3, draw_on_top=True)


p_prod.text    (1.2e1, 2e-33 , r'$\mathrm{T_{\bar{p}}= 0.5\,GeV\,\,(\times\, 40) }$',   color=c1, rotation=85 )
p_prod.text    (0.5e2, 2e-34 , r'$\mathrm{T_{\bar{p}}=                    20\,GeV}$',   color=c2, rotation=80 )
p_prod.text    (  7e2, 2e-34 , r'$\mathrm{T_{\bar{p}}=                   0.5\,TeV}$',   color=c3, rotation=75 )

if file!='' and len(pp):
    profile(     pp_file,                500  , resfile='prod_dT_pp_pbar.pdf', type='prod', dashes=(        ), Tn_min=5, Tn_max=1e5, y_fac = fac, draw_on_top=True, color='red', alpha=0.5 )
profile('dT_pp_pbar_LAB_TanNg.txt',      500  , resfile='prod_dT_pp_pbar.pdf', type='prod', color=c3, dashes=(12,3    ), Tn_min=5, Tn_max=1e5, y_fac = 2.3, draw_on_top=True)
profile('dT_pp_pbar_LAB_KR_PPFRAG.txt',  500  , resfile='prod_dT_pp_pbar.pdf', type='prod', color=c3, dashes=(3, 3    ), Tn_min=5, Tn_max=1e5, y_fac = 2.3, draw_on_top=True)
profile('dT_pp_pbar_LAB_Duperray.txt',   500  , resfile='prod_dT_pp_pbar.pdf', type='prod', color=c3, dashes=(20,5,5,5), Tn_min=5, Tn_max=1e5, y_fac = 2.3, draw_on_top=True)
profile('dT_pp_pbar_LAB_diMauro12.txt',  500  , resfile='prod_dT_pp_pbar.pdf', type='prod', color=c3, dashes=(        ), Tn_min=5, Tn_max=1e5, y_fac = 2.3, draw_on_top=True, legend=True, x_label = r'$\mathrm{T_{p} \quad [GeV]}$' )


p_proj = profile('dT_pp_pbar_LAB_diMauro12.txt',     20, resfile='proj_dT_pp_pbar.pdf', type='proj', color=c0, dashes=(        ), y_fac = 2.3, label='di Mauro' )
profile(         'dT_pp_pbar_LAB_TanNg.txt',         20, resfile='proj_dT_pp_pbar.pdf', type='proj', color=c0, dashes=(12,3    ), y_fac = 2.3, label='Tan & Ng' , draw_on_top=True)
profile(         'dT_pp_pbar_LAB_KR_PPFRAG.txt',     20, resfile='proj_dT_pp_pbar.pdf', type='proj', color=c0, dashes=(3, 3    ), y_fac = 2.3, label='KMO'   , draw_on_top=True)
profile(         'dT_pp_pbar_LAB_Duperray.txt',      20, resfile='proj_dT_pp_pbar.pdf', type='proj', color=c0, dashes=(20,5,5,5), y_fac = 2.3, label='Duperray' , draw_on_top=True)
if file!='' and len(pp):
    profile(     pp_file,                            20, resfile='proj_dT_pp_pbar.pdf', type='proj', dashes=(        ), y_fac = fac, label=label ,      draw_on_top=True, color='red', alpha=0.5)
profile(         'dT_pp_pbar_LAB_diMauro12.txt',     20, resfile='proj_dT_pp_pbar.pdf', type='proj', color=c1, dashes=(        ), y_fac = 2.3, draw_on_top=True )
profile(         'dT_pp_pbar_LAB_TanNg.txt',         20, resfile='proj_dT_pp_pbar.pdf', type='proj', color=c1, dashes=(12,3    ), y_fac = 2.3, draw_on_top=True )
profile(         'dT_pp_pbar_LAB_KR_PPFRAG.txt',     20, resfile='proj_dT_pp_pbar.pdf', type='proj', color=c1, dashes=(3, 3    ), y_fac = 2.3, draw_on_top=True )
profile(         'dT_pp_pbar_LAB_Duperray.txt',      20, resfile='proj_dT_pp_pbar.pdf', type='proj', color=c1, dashes=(20,5,5,5), y_fac = 2.3, draw_on_top=True )


if file!='' and len(pp):
    profile(     pp_file,                  500, resfile='proj_dT_pp_pbar.pdf', type='proj', dashes=(        ), y_fac = fac,                    draw_on_top=True, color='red', alpha=0.5)
profile('dT_pp_pbar_LAB_TanNg.txt',        500, resfile='proj_dT_pp_pbar.pdf', type='proj', color=c2, dashes=(12,3    ), y_fac = 2.3,                    draw_on_top=True)
profile('dT_pp_pbar_LAB_KR_PPFRAG.txt',    500, resfile='proj_dT_pp_pbar.pdf', type='proj', color=c2, dashes=(3, 3    ), y_fac = 2.3,                    draw_on_top=True, Tn_max=0.9e2 )
profile('dT_pp_pbar_LAB_Duperray.txt',     500, resfile='proj_dT_pp_pbar.pdf', type='proj', color=c2, dashes=(20,5,5,5), y_fac = 2.3,                    draw_on_top=True )
profile('dT_pp_pbar_LAB_diMauro12.txt',    500, resfile='proj_dT_pp_pbar.pdf', type='proj', color=c2, dashes=(        ), y_fac = 2.3,                    draw_on_top=True )

p_proj.text    ( 0.8e1, 2e-34 , r'$\mathrm{T_{p}=                   20\,GeV}$',   color=c1,     rotation=-80 )
p_proj.text    (   2e2, 2e-34 , r'$\mathrm{T_{p}=                  500\,GeV}$',   color=c2,     rotation=-80 )
p_proj.text    (   1e3, 2e-33 , r'$\mathrm{T_{p}=  8\,TeV\,\,(\times\,5)   }$',   color=c3,     rotation=-70 )

if file!='' and len(pp):
    profile(     pp_file,                 8000, resfile='proj_dT_pp_pbar.pdf', type='proj', dashes=(        ), y_fac = 5*fac,                    draw_on_top=True, color='red', alpha=0.5 )
profile('dT_pp_pbar_LAB_TanNg.txt',       8000, resfile='proj_dT_pp_pbar.pdf', type='proj', color=c3, dashes=(12,3    ), y_fac = 5*2.3,                    draw_on_top=True)
profile('dT_pp_pbar_LAB_KR_PPFRAG.txt',   8000, resfile='proj_dT_pp_pbar.pdf', type='proj', color=c3, dashes=(3, 3    ), y_fac = 5*2.3,                    draw_on_top=True, Tn_max=0.9e2 )
profile('dT_pp_pbar_LAB_Duperray.txt',    8000, resfile='proj_dT_pp_pbar.pdf', type='proj', color=c3, dashes=(20,5,5,5), y_fac = 5*2.3,                    draw_on_top=True  )
profile('dT_pp_pbar_LAB_diMauro12.txt',   8000, resfile='proj_dT_pp_pbar.pdf', type='proj', color=c3, dashes=(        ), y_fac = 5*2.3,                    draw_on_top=True, legend=True, x_label = r'$\mathrm{T_{\bar{p}} \quad [GeV]}$' )




if file!='' and len(pp):
    if label!='':
        plot2D(pp_file,            resfile='dT_pp_pbar_LAB_'  + file + '.pdf'     )
        plot2D('dT_pp_pbar_LAB_diMauro12.txt',            file2=pp_file,       resfile='ratio_pp_diMauro12_'+label+'.pdf',  Zlabel = label+'/diMauro12' )
    sys.exit(0)



plot2D('dT_pp_pbar_LAB_diMauro.txt',            resfile='dT_pp_pbar_LAB_diMauro.pdf'     , x_label=r'$\mathrm{T_{p}\quad [GeV]}$'  )
plot2D('dT_pp_pbar_LAB_diMauro12.txt',          resfile='dT_pp_pbar_LAB_diMauro12.pdf'   , x_label=r'$\mathrm{T_{p}\quad [GeV]}$'  )
plot2D('dT_pp_pbar_LAB_TanNg.txt',              resfile='dT_pp_pbar_LAB_TanNg.pdf'       , x_label=r'$\mathrm{T_{p}\quad [GeV]}$'  )

plot2D('dT_pp_pbar_LAB_KR_PPFRAG.txt',          resfile='dT_pp_pbar_LAB_KR.pdf'          , x_label=r'$\mathrm{T_{p}\quad [GeV]}$'  )
plot2D('dT_Hep_pbar_LAB_KR_PPFRAG.txt',         resfile='dT_Hep_pbar_LAB_KR.pdf'         )
plot2D('dT_pHe_pbar_LAB_KR_PPFRAG.txt',         resfile='dT_pHe_pbar_LAB_KR.pdf'         )



plot2D('dT_pp_pbar_LAB_Duperray.txt',           resfile='dT_pp_pbar_LAB_Duperray.pdf'      , x_label=r'$\mathrm{T_{p}\quad [GeV]}$'  )
plot2D('dT_Hep_pbar_LAB_Duperray.txt',          resfile='dT_Hep_pbar_LAB_Duperray.pdf'     )
plot2D('dT_pHe_pbar_LAB_Duperray.txt',          resfile='dT_pHe_pbar_LAB_Duperray.pdf'     )


plot2D('dT_pp_pbar_LAB_diMauro.txt',            file2='dT_pp_pbar_LAB_diMauro12.txt',       resfile='ratio_pp_diMauro12_diMauro.pdf',  Zlabel = 'diMauro12/diMauro13' )
plot2D('dT_pp_pbar_LAB_diMauro.txt',            file2='dT_pp_pbar_LAB_TanNg.txt',           resfile='ratio_pp_TanNg_diMauro.pdf',      Zlabel = 'Tan&Ng/diMauro'      )
plot2D('dT_pp_pbar_LAB_diMauro.txt',            file2='dT_pp_pbar_LAB_Duperray.txt',        resfile='ratio_pp_Duperray_diMauro.pdf',   Zlabel = 'Duperray/diMauro',   )
plot2D('dT_pp_pbar_LAB_diMauro.txt',            file2='dT_pp_pbar_LAB_KR_PPFRAG.txt',       resfile='ratio_pp_KR_diMauro.pdf',         Zlabel = 'KR/diMauro',           zscale='log',   Zmin=1e-1, Zmax=1e1   )








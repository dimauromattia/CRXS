#! /usr/bin/env python

import glob
import math
import numpy as np

import XS_transformations as t


p_col = [200, 300, 400]

def readTab( file, col, err_scale, tab=1 ):
    global f, fScale
    global p_col
    p = p_col[col]
    
    s = t.s__from_pLAB(p)
    
    line = '\n#\n#\n# Tab '+str(tab)+'\n'
    f.write( line )


    line  = '#*  '
    line += str('sqrt(s)'   ).ljust(20)+' '
    line += str('pT'        ).ljust(20)+' '
    line += str('xR'        ).ljust(20)+' '
    line += str('theta_LAB' ).ljust(20)+' '
    line += str('pT'        ).ljust(20)+' '
    line += str('inv CS'    ).ljust(20)+' '
    line += str('CS_stat'   ).ljust(20)+' '
    line += str('CS_sys'    ).ljust(20)+' '
    line += str('err_scale' ).ljust(20)+'\n'
    f.write( line )



    tab1 = np.genfromtxt( file, skip_header=9 )
    for i in range( len(tab1[:,0]) ):
        
        #factor = np.sqrt(t.fMass_proton*t.fMass_proton+p_pbar*p_pbar)/p_pbar/p_pbar*1e3
        
        pT   =  tab1[i,0]
        CS   =  tab1[i,3+3*col]*1e27
        stat =  tab1[i,4+3*col]*1e27
        
        if CS!=CS:
            continue
        
        sys  =  0
        scale = err_scale
        
        theta  = 77*1e-3
        
        eta = t.theta_to_eta(theta)
        
        p_pbar = pT*np.cosh(eta)

        s, pT, xR = t.convert_LAB_to_CM__from__s_ppbar_eta(s, p_pbar, eta)
        
        line = '   '
        line += str(np.sqrt(s)  ).ljust(20)+' '
        line += str(pT          ).ljust(20)+' '
        line += str(xR          ).ljust(20)+' '
        line += str(theta       ).ljust(20)+' '
        line += str(pT          ).ljust(20)+' '
        line += str(CS          ).ljust(20)+' '
        line += str(stat        ).ljust(20)+' '
        line += str(sys         ).ljust(20)+' '
        line += str(scale       ).ljust(20)+'\n'
        f.write( line )

    f.write( '#\n#End Tab '+str(tab)+'\n#' )


fScale=0.20

def write(data_dir, res_dir):
    global f, fScale
    f = open(res_dir+'/converted_antreasyan.txt', 'w')

    line = '# Antreasyan et al, DOI:10.1103/PhysRevD.19.764\n'
    f.write( line )
    line = '# p + p -> pbar + X \n#'
    f.write( line )

    readTab( data_dir+'/antreasyan_8.txt', 0, fScale, 8   )
    readTab( data_dir+'/antreasyan_8.txt', 1, fScale, 8   )
    readTab( data_dir+'/antreasyan_8.txt', 2, fScale, 8   )


    f.close()















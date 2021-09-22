#! /usr/bin/env python

import glob
import math
import numpy as np

import CRACS_CS_transformations as t

def readTab( file, err_scale, tab=1 ):
    global f, fScale

    line = '\n#\n#\n# Tab '+str(tab)+'\n'
    f.write( line )


    line  = '#*  '
    line += str('sqrt(s)'   ).ljust(20)+' '
    line += str('pT'        ).ljust(20)+' '
    line += str('xR'        ).ljust(20)+' '
    line += str('theta_LAB' ).ljust(20)+' '
    line += str('xR'        ).ljust(20)+' '
    line += str('inv CS'    ).ljust(20)+' '
    line += str('CS_stat'   ).ljust(20)+' '
    line += str('CS_sys'    ).ljust(20)+' '
    line += str('err_scale' ).ljust(20)+'\n'
    f.write( line )



    tab1 = np.genfromtxt( file, skip_header=9 )
    for i in range( len(tab1[:,0]) ):
        
        p_p    =    tab1[i,0]
        s      =    t.s__from_pLAB(p_p)
        p_pbar =    tab1[i,1]
        
        factor = np.sqrt(t.fMass_proton*t.fMass_proton+p_pbar*p_pbar)/p_pbar/p_pbar
        
        CS   =  tab1[i,4]*factor
        stat =  tab1[i,5]*factor
        
        if CS!=CS:
            continue
        
        sys  =  0
        scale = err_scale
        
        theta  = tab1[i,3]*1e-3
        

        s, pT, xR = t.convert_LAB_to_CM__from__s_ppbar_theta(s, p_pbar, theta)
        
        line = '   '
        line += str(np.sqrt(s)  ).ljust(20)+' '
        line += str(pT          ).ljust(20)+' '
        line += str(xR          ).ljust(20)+' '
        line += str(theta       ).ljust(20)+' '
        line += str(xR          ).ljust(20)+' '
        line += str(CS          ).ljust(20)+' '
        line += str(stat        ).ljust(20)+' '
        line += str(sys         ).ljust(20)+' '
        line += str(scale       ).ljust(20)+'\n'
        f.write( line )

    f.write( '#\n#End Tab '+str(tab)+'\n#' )


fScale=0.10

def write(data_dir, res_dir):
    global f, fScale
    f = open(res_dir+'/converted_dekkers.txt', 'w')

    line = '# Dekkers et al, doi: doi.org/10.1103/PhysRev.137.B962\n'
    f.write( line )
    line = '# p + p -> pbar + X \n#'
    f.write( line )

    readTab( data_dir+'/dekkers.txt', fScale, 1   )


    f.close()















#! /usr/bin/env python

import glob
import math
import numpy as np
import os
CRACS = os.environ['CRACS']+'/scripts/CS_transformation_pp'

print CRACS

import sys
sys.path.append(CRACS)

import CRACS_CS_transformations as t


def readTab( file, p_p, theta_LAB, err_scale, tab=1 ):
    global f, fScale

    line = '\n#\n#\n# Tab '+str(tab)+'\n'
    f.write( line )


    line  = '#*  '
    line += str('sqrt(s)'   ).ljust(20)+' '
    line += str('pT'        ).ljust(20)+' '
    line += str('xR'        ).ljust(20)+' '
    line += str('theta_LAB' ).ljust(20)+' '
    line += str('p_LAB'     ).ljust(20)+' '
    line += str('inv CS'    ).ljust(20)+' '
    line += str('CS_stat'   ).ljust(20)+' '
    line += str('CS_sys'    ).ljust(20)+' '
    line += str('err_scale' ).ljust(20)+'\n'
    f.write( line )



    tab1 = np.genfromtxt( file, skip_header=9 )
    for i in range( len(tab1[:,0]) ):
       
        E       =  np.sqrt(p_p*p_p+t.fMass_proton*t.fMass_proton)
       
       
        p_pbar  =  tab1[i,0]
        E_pbar  =  np.sqrt(p_pbar*p_pbar+t.fMass_proton*t.fMass_proton)
        
        factor = E_pbar/p_pbar/p_pbar
       
        
        CS   =  tab1[i,3+5*(tab-1)]*factor
        stat =  tab1[i,4+5*(tab-1)]*factor
        sys  =  tab1[i,6+5*(tab-1)]*factor


        if CS!=CS:
            continue
        
        scale = err_scale
        
        s = 2*(E+t.fMass_proton)*t.fMass_proton
        
        theta  = 0
        
        
        
        
        s_2, pT, xR = t.convert_LAB_to_CM__from__s_ppbar_theta(s, p_pbar, theta)
        
        print np.sqrt(s)
        print pT
        print xR
        
        line = '   '
        line += str(np.sqrt(s)  ).ljust(20)+' '
        line += str(pT          ).ljust(20)+' '
        line += str(xR          ).ljust(20)+' '
        line += str(theta       ).ljust(20)+' '
        line += str(p_pbar      ).ljust(20)+' '
        line += str(CS          ).ljust(20)+' '
        line += str(stat        ).ljust(20)+' '
        line += str(sys         ).ljust(20)+' '
        line += str(scale       ).ljust(20)+'\n'
        f.write( line )

    f.write( '#\n#End Tab '+str(tab)+'\n#' )


fScale=0.10  # assumption, not given

def write(data_dir, res_dir):
    global f, fScale
    f = open(res_dir+'/converted_amann.txt', 'w')

    line = '# Amann et al, \n'
    f.write( line )
    line = '# p + C -> pbar + X \n#'
    f.write( line )

    readTab( data_dir + '/amann.txt',  10, 0, fScale, 1   )
    readTab( data_dir + '/amann.txt',  18, 0, fScale, 2   )
    readTab( data_dir + '/amann.txt',  24, 0, fScale, 3   )
    
    f.close()















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


def readTab( file, E, theta_LAB, err_scale, tab=1 ):
    global f, fScale

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



    tab1 = np.genfromtxt( file, skip_header=10 )
    for i in range( len(tab1[:,0]) ):
        
        
        CS   =  tab1[i,3]
        stat =  tab1[i,4]
        sys  =  0

        pT   =  tab1[i,0]
        
        if CS!=CS:
            continue
        
        sys  =  0
        scale = err_scale
        
        s = 2*(E+t.fMass_proton)*t.fMass_proton
        
        theta  = theta_LAB*6.2831853072/360.
        
        eta = t.theta_to_eta(theta)
        
        p_pbar = pT/np.sin(theta)

        s_2, pT_2, xR = t.convert_LAB_to_CM__from__s_ppbar_eta(s, p_pbar, eta)
        
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


fScale=0.10 # assumption, not given

def write(data_dir, res_dir):
    global f, fScale
    f = open(res_dir+'/converted_abramov.txt', 'w')

    line = '# Abramov et al, \n'
    f.write( line )
    line = '# p + C -> pbar + X \n#'
    f.write( line )

    readTab( data_dir + '/abramov.txt',  70, 12.01, fScale, 1   )
    
    f.close()















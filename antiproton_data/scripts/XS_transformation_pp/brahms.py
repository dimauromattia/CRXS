#! /usr/bin/env python

import glob
import math
import numpy as np

import XS_transformations as t



def readTab( file, err_scale, tab=1 ):
    global f, fScale
    global p_col
    
    s = 200*200
    
    line = '\n#\n#\n# Tab '+str(tab)+'\n'
    f.write( line )


    line  = '#*  '
    line += str('sqrt(s)'   ).ljust(20)+' '
    line += str('pT'        ).ljust(20)+' '
    line += str('xR'        ).ljust(20)+' '
    line += str('y_LAB'     ).ljust(20)+' '
    line += str('pT'        ).ljust(20)+' '
    line += str('inv CS'    ).ljust(20)+' '
    line += str('CS_stat'   ).ljust(20)+' '
    line += str('CS_sys'    ).ljust(20)+' '
    line += str('err_scale' ).ljust(20)+'\n'
    f.write( line )



    tab1 = np.genfromtxt( file, skip_header=10 )
    for i in range( len(tab1[:,0]) ):
        
        #factor = np.sqrt(t.fMass_proton*t.fMass_proton+p_pbar*p_pbar)/p_pbar/p_pbar*1e3
        
        pT   =  tab1[i,0]
        CS   =  tab1[i,3]
        stat =  tab1[i,4]
        sys  =  tab1[i,6]
        
        print( sys/CS )
        
        if CS!=CS:
            continue
        
        scale = err_scale
        
        y   = 2.95
        
        eta = t.y_to_eta(y, pT)
        
        #print( 'eta:')
        #print( eta)
        
        #print( 'pT:')
        #print( pT)
        
        p_pbar = pT*np.cosh(eta)        
        xR = t.x_R__from_s_ppbar(s, p_pbar)

        line = '   '
        line += str(np.sqrt(s)  ).ljust(20)+' '
        line += str(pT          ).ljust(20)+' '
        line += str(xR          ).ljust(20)+' '
        line += str(y           ).ljust(20)+' '
        line += str(pT          ).ljust(20)+' '
        line += str(CS          ).ljust(20)+' '
        line += str(stat        ).ljust(20)+' '
        line += str(sys         ).ljust(20)+' '
        line += str(scale       ).ljust(20)+'\n'
        f.write( line )

    f.write( '#\n#End Tab '+str(tab)+'\n#' )

fScale=0.10

def write(data_dir, res_dir):
    global f, fScale
    f = open(res_dir+'/converted_brahms.txt', 'w')

    line = '# BRAHMS, DOI:10.1103/PhysRevLett.98.252001\n'
    f.write( line )
    line = '# p + p -> pbar + X \n#'
    f.write( line )

    readTab( data_dir+'/brahms6.txt', fScale, 6   )



    f.close()















#! /usr/bin/env python

import glob
import math
import numpy as np

import CRACS_CS_transformations as t



def readTab( file, p, pT, err_scale, tab=1 ):
    global f, fScale
    global p_col
    
    s = t.s__from_pLAB(p)
    
    line = '\n#\n#\n# Tab '+str(tab)+'\n'
    f.write( line )


    line  = '#*  '
    line += str('sqrt(s)'   ).ljust(20)+' '
    line += str('pT'        ).ljust(20)+' '
    line += str('xR'        ).ljust(20)+' '
    line += str('pT'        ).ljust(20)+' '
    line += str('xR'        ).ljust(20)+' '
    line += str('inv CS'    ).ljust(20)+' '
    line += str('CS_stat'   ).ljust(20)+' '
    line += str('CS_sys'    ).ljust(20)+' '
    line += str('err_scale' ).ljust(20)+'\n'
    f.write( line )



    tab1 = np.genfromtxt( file, skip_header=10 )
    for i in range( len(tab1[:,0]) ):
        
        #factor = np.sqrt(t.fMass_proton*t.fMass_proton+p_pbar*p_pbar)/p_pbar/p_pbar*1e3
        
        xR   =  tab1[i,0]
        CS   =  tab1[i,3]*1e27
        stat =  tab1[i,4]*1e27
        sys  =  0
        
        if CS!=CS:
            continue
        
        scale = err_scale
        
        
        line = '   '
        line += str(np.sqrt(s)  ).ljust(20)+' '
        line += str(pT          ).ljust(20)+' '
        line += str(xR          ).ljust(20)+' '
        line += str(pT          ).ljust(20)+' '
        line += str(xR          ).ljust(20)+' '
        line += str(CS          ).ljust(20)+' '
        line += str(stat        ).ljust(20)+' '
        line += str(sys         ).ljust(20)+' '
        line += str(scale       ).ljust(20)+'\n'
        f.write( line )

    f.write( '#\n#End Tab '+str(tab)+'\n#' )


fScale=0.07

def write(data_dir, res_dir):
    global f, fScale
    f = open(res_dir+'/converted_johnson.txt', 'w')

    line = '# Johnson et al, DOI:10.1103/PhysRevLett.39.1173\n'
    f.write( line )
    line = '# p + p -> pbar + X \n#'
    f.write( line )

    readTab(  data_dir+'/johnson90.txt', 100, 0.25, fScale, 90   )
    readTab(  data_dir+'/johnson91.txt', 100, 0.5 , fScale, 91   )
    readTab(  data_dir+'/johnson92.txt', 100, 0.75, fScale, 92   )
    readTab(  data_dir+'/johnson93.txt', 200, 0.25, fScale, 93   )
    readTab(  data_dir+'/johnson94.txt', 200, 0.5 , fScale, 94   )
    readTab(  data_dir+'/johnson95.txt', 200, 0.75, fScale, 95   )
    readTab(  data_dir+'/johnson96.txt', 400, 0.25, fScale, 96   )
    readTab(  data_dir+'/johnson97.txt', 400, 0.5 , fScale, 97   )
    readTab(  data_dir+'/johnson98.txt', 400, 0.75, fScale, 98   )




    f.close()















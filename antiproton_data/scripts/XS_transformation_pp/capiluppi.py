#! /usr/bin/env python

import glob
import math
import numpy as np

import XS_transformations as t


def readTab( file, err_scale, tab=1 ):
    global f, fScale

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



    tab1 = np.genfromtxt( file, skip_header=14 )
    for i in range( len(tab1[:,0]) ):
        CS   =  tab1[i,4]
        stat =  tab1[i,5]
        sys  =  0
        scale = err_scale
        
        s   = tab1[i,0]*tab1[i,0]
        pT  = tab1[i,1]
        xF  = tab1[i,2]


        xR  = t.x_R__from_s_pT_xF(s, pT, xF)
        
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



fScale=0.10

def write(data_dir, res_dir):
    global f, fScale
    f = open(res_dir+'/converted_capiluppi.txt', 'w')

    line = '# Capiluppi et al, doi 10.1016/0550-3213(74)90484-2\n'
    f.write( line )
    line = '# p + p -> pbar + X \n #'
    f.write( line )

    readTab( data_dir+'/capiluppi.txt', fScale   )

    f.close()















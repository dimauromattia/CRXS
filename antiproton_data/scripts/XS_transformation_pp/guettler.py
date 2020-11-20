#! /usr/bin/env python

import glob
import math
import numpy as np

import XS_transformations as t


def readTab( file, s, err_scale, tab=1 ):
    global f, fScale

    line = '\n#\n#\n# Tab '+str(tab)+'\n'
    f.write( line )


    line  = '#*  '
    line += str('sqrt(s)'   ).ljust(20)+' '
    line += str('pT'        ).ljust(20)+' '
    line += str('xR'        ).ljust(20)+' '
    line += str('y'         ).ljust(20)+' '
    line += str('pT'        ).ljust(20)+' '
    line += str('inv CS'    ).ljust(20)+' '
    line += str('CS_stat'   ).ljust(20)+' '
    line += str('CS_sys'    ).ljust(20)+' '
    line += str('err_scale' ).ljust(20)+'\n'
    f.write( line )



    tab1 = np.genfromtxt( file, skip_header=9 )
    for i in range( len(tab1[:,0]) ):
        CS   =  tab1[i,3]
        stat =  tab1[i,4]
        sys  =  0
        scale = err_scale
        
        pT  = tab1[i,0]
        
        xR = t.x_R__from_s_ppbar ( s, pT )
        
        line = '   '
        line += str(np.sqrt(s)  ).ljust(20)+' '
        line += str(pT          ).ljust(20)+' '
        line += str(xR          ).ljust(20)+' '
        line += str(0           ).ljust(20)+' '
        line += str(pT          ).ljust(20)+' '
        line += str(CS          ).ljust(20)+' '
        line += str(stat        ).ljust(20)+' '
        line += str(sys         ).ljust(20)+' '
        line += str(scale       ).ljust(20)+'\n'
        f.write( line )

    f.write( '#\n#End Tab '+str(tab)+'\n#' )




fScale=0.04

def write(data_dir, res_dir):
    global f, fScale
    f = open(res_dir+'/converted_guettler.txt', 'w')

    line = '# guettler et al, DOI:10.1016/0550-3213(76)90313-8\n'
    f.write( line )
    line = '# p + p -> pbar + X \n'
    f.write( line )

    readTab( data_dir+'/guettler_23.txt', 23.*23., fScale   )
    readTab( data_dir+'/guettler_31.txt', 31.*31., fScale   )
    readTab( data_dir+'/guettler_45.txt', 45.*45., fScale   )
    readTab( data_dir+'/guettler_53.txt', 53.*53., fScale   )
    readTab( data_dir+'/guettler_63.txt', 63.*63., fScale   )


    f.close()















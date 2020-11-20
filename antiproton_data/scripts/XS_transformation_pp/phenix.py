#! /usr/bin/env python

import glob
import math
import numpy as np

import XS_transformations as t


fScale=0.10

def write(data_dir, res_dir):
    global f, fScale
    f = open(res_dir+'/converted_phenix.txt', 'w')

    line = '# Phenix collaboration, doi: 10.1103/PhysRevC.83.064903, Tab 14\n'
    f.write( line )
    line = '# p + p -> pbar + X, Feed-down weak decay corrections are applied\n'
    f.write( line )
    line = '# Additional scale uncertainty: 9.7% \n'
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



    tab1 = np.genfromtxt(  data_dir+'/phenix_200.txt', skip_header=8 )
    for i in range( len(tab1[:,0]) ):
        CS   =  tab1[i,8]
        stat =  tab1[i,9]
        sys  =  tab1[i,11]
        scale = 0.097
        
        s   = 200.*200
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

    f.write( '#\n#End Tab 14\n#' )


    line = '\n#\n#\n#\n#\n# Tab 18\n'
    f.write( line )
    line = '# p + p -> pbar + X, Feed-down weak decay corrections are applied\n'
    f.write( line )
    line = '# Additional scale uncertainty: 11% \n'
    f.write( line )

    line  = '#*  '
    line += str('sqrt(s)' ).ljust(20)+' '
    line += str('pT'        ).ljust(20)+' '
    line += str('xR'        ).ljust(20)+' '
    line += str('y'         ).ljust(20)+' '
    line += str('pT'        ).ljust(20)+' '
    line += str('inv CS'    ).ljust(20)+' '
    line += str('CS_stat'   ).ljust(20)+' '
    line += str('CS_sys'    ).ljust(20)+' '
    line += str('err_scale' ).ljust(20)+'\n'
    f.write( line )


    tab2 = np.genfromtxt( data_dir+'/phenix_62.txt', skip_header=8 )
    for i in range( len(tab2[:,0]) ):
        CS   =  tab2[i,8]
        stat =  tab2[i,9]
        sys  =  tab2[i,11]
        scale = 0.11
        
        s   = 62.4*62.4
        pT  = tab2[i,0]
        
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

    f.write( '#\n#End Tab 18\n#' )
    f.close()















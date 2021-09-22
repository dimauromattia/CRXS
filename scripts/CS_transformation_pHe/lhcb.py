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

fMass_proton = 0.9382720813

def readTab( file, err_scale, tab=1 ):
    global f, fScale

    line = '\n#\n#\n# Tab '+str(tab)+'\n'
    f.write( line )


    line  = '#*  '
    line += str('sqrt(s)'   ).ljust(20)+' '
    line += str('pT'        ).ljust(20)+' '
    line += str('xR'        ).ljust(20)+' '
    line += str('pT'        ).ljust(20)+' '
    line += str('xf'        ).ljust(20)+' '
    line += str('inv CS'    ).ljust(20)+' '
    line += str('CS_stat'   ).ljust(20)+' '
    line += str('CS_sys'    ).ljust(20)+' '
    line += str('err_scale' ).ljust(20)+'\n'
    f.write( line )



    tab1 = np.genfromtxt( file  )
    for i in range( len(tab1[:,0]) ):
        CS   =  tab1[i,7]
        stat =  tab1[i,8]
        sys  =  tab1[i,9]
        scale = err_scale
        
        p   = tab1[i,4]
        pT  = tab1[i,5]
        
        p_p = 6.5e3
        
        s   = t.s__from_pLAB(p_p)
        
        ss, pT_pbar, x_R, pL_pbar  = t.convert_LAB_to_CM__from__s_ppbar_theta(  s, p, np.arcsin(pT/p), with_pL=True  )
        
        if pL_pbar<0:
            x_R = x_R*(-1.)
        
        x_F = 2*pL_pbar/np.sqrt(s)
        
#        print s
#        print ss
#        print pT
#        print pT_pbar
#        print pL_pbar
#        
#        print x_R

        pL = np.sqrt(p*p-pT*pT)
        
        E = np.sqrt(p*p+fMass_proton*fMass_proton)
        
        factor = E * pL / pT / p / 2 / np.pi * 1e-3
        
        CS   *= factor
        stat *= factor
        sys  *= factor

        line = '   '
        line += str(np.sqrt(s)  ).ljust(20)+' '
        line += str(pT          ).ljust(20)+' '
        line += str(x_R         ).ljust(20)+' '
        line += str(pT          ).ljust(20)+' '
        line += str(x_F         ).ljust(20)+' '
        line += str(CS          ).ljust(20)+' '
        line += str(stat        ).ljust(20)+' '
        line += str(sys         ).ljust(20)+' '
        line += str(scale       ).ljust(20)+'\n'
        f.write( line )

    f.write( '#\n#End Tab '+str(tab)+'\n#' )




fScale=0.06

def write(data_dir, res_dir):
    global f, fScale
    f = open(res_dir+'/converted_lhcb.txt', 'w')

    line = '# LHCb \n'
    f.write( line )
    line = '# p + He -> pbar + X \n #'
    f.write( line )

    readTab( data_dir+'/lhcb.txt', fScale   )


    f.close()















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



# TO BE IMPLEMENTED!!
#
#def readTab( file, p_p, theta_LAB, err_scale, tab=1 ):
#    global f, fScale
#
#    line = '\n#\n#\n# Tab '+str(tab)+'\n'
#    f.write( line )
#
#
#    line  = '#*  '
#    line += str('sqrt(s)'   ).ljust(20)+' '
#    line += str('pT'        ).ljust(20)+' '
#    line += str('xR'        ).ljust(20)+' '
#    line += str('theta_LAB' ).ljust(20)+' '
#    line += str('p'         ).ljust(20)+' '
#    line += str('inv CS'    ).ljust(20)+' '
#    line += str('CS_stat'   ).ljust(20)+' '
#    line += str('CS_sys'    ).ljust(20)+' '
#    line += str('err_scale' ).ljust(20)+'\n'
#    f.write( line )
#
#
#
#    tab1 = np.genfromtxt( file, skip_header=8 )
#    for i in range( len(tab1[:,0]) ):
#        
#        E = np.sqrt(p_p*p_p + t.fMass_proton*t.fMass_proton)
#        
#        CS   =  tab1[i,3]*1000
#        stat =  tab1[i,4]*1000
#        sys  =  0
#
#        p_pbar   =  tab1[i,0]
#        
#        if CS!=CS:
#            continue
#        
#        scale = err_scale
#        
#        s = 2*(E+t.fMass_proton)*t.fMass_proton
#        
#        theta  = theta_LAB.
#        
#        eta = t.theta_to_eta(theta)
#        
#        s_2, pT, xR = t.convert_LAB_to_CM__from__s_ppbar_eta(s, p_pbar, eta)
#        
#        line = '   '
#        line += str(np.sqrt(s)  ).ljust(20)+' '
#        line += str(pT          ).ljust(20)+' '
#        line += str(xR          ).ljust(20)+' '
#        line += str(theta       ).ljust(20)+' '
#        line += str(p_pbar      ).ljust(20)+' '
#        line += str(CS          ).ljust(20)+' '
#        line += str(stat        ).ljust(20)+' '
#        line += str(sys         ).ljust(20)+' '
#        line += str(scale       ).ljust(20)+'\n'
#        f.write( line )
#
#    f.write( '#\n#End Tab '+str(tab)+'\n#' )
#
#
#fScale=0.20
#
#def write(data_dir, res_dir):
#    global f, fScale
#    f = open(res_dir+'/converted_lepikhin.txt', 'w')
#
#    line = '# lepikhin et al, \n'
#    f.write( line )
#    line = '# http://hepdata.cedar.ac.uk/view/ins256006 \n'
#    f.write( line )
#    line = '# p + Be -> pbar + X \n#'
#    f.write( line )
#
#    readTab( data_dir + '/lepikhin.txt',  10.1, 3.5, fScale, 1   )
#    
#    f.close()
#














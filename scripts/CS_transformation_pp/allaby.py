#! /usr/bin/env python

import glob
import math
import numpy as np

import CRACS_CS_transformations as t




def readTab( file, s, p_pbar, err_scale, tab=1 ):
    global f, fScale

    line = '\n#\n#\n# Tab '+str(tab)+'\n'
    f.write( line )


    line  = '#*  '
    line += str('sqrt(s)'   ).ljust(20)+' '
    line += str('pT'        ).ljust(20)+' '
    line += str('xR'        ).ljust(20)+' '
    line += str('p_pbar_LAB').ljust(20)+' '
    line += str('theta_LAB' ).ljust(20)+' '
    line += str('inv CS'    ).ljust(20)+' '
    line += str('CS_stat'   ).ljust(20)+' '
    line += str('CS_sys'    ).ljust(20)+' '
    line += str('err_scale' ).ljust(20)+'\n'
    f.write( line )



    tab1 = np.genfromtxt( file, skip_header=9 )
    for i in range( len(tab1[:,0]) ):
        
        factor = np.sqrt(t.fMass_proton*t.fMass_proton+p_pbar*p_pbar)/p_pbar/p_pbar*1e3
        
        CS   =  tab1[i,15]*factor
        stat =  tab1[i,16]*factor
        
        if CS!=CS:
            continue
        
        sys  =  0
        scale = err_scale
        
        theta  = tab1[i,0]*1e-3
        
        eta = t.theta_to_eta(theta)

        s, pT, xR = t.convert_LAB_to_CM__from__s_ppbar_eta(s, p_pbar, eta)
        
        line = '   '
        line += str(np.sqrt(s)  ).ljust(20)+' '
        line += str(pT          ).ljust(20)+' '
        line += str(xR          ).ljust(20)+' '
        line += str(p_pbar      ).ljust(20)+' '
        line += str(theta       ).ljust(20)+' '
        line += str(CS          ).ljust(20)+' '
        line += str(stat        ).ljust(20)+' '
        line += str(sys         ).ljust(20)+' '
        line += str(scale       ).ljust(20)+'\n'
        f.write( line )

    f.write( '#\n#End Tab '+str(tab)+'\n#' )


def readTab_t( file, s, theta, err_scale, tab=1 ):
    global f, fScale
    
    line = '\n#\n#\n# Tab '+str(tab)+'\n'
    f.write( line )
    
    
    line  = '#*  '
    line += str('sqrt(s)'   ).ljust(20)+' '
    line += str('pT'        ).ljust(20)+' '
    line += str('xR'        ).ljust(20)+' '
    line += str('p_pbar_LAB').ljust(20)+' '
    line += str('theta_LAB' ).ljust(20)+' '
    line += str('inv CS'    ).ljust(20)+' '
    line += str('CS_stat'   ).ljust(20)+' '
    line += str('CS_sys'    ).ljust(20)+' '
    line += str('err_scale' ).ljust(20)+'\n'
    f.write( line )
    
    
    
    tab1 = np.genfromtxt( file, skip_header=9 )
    for i in range( len(tab1[:,0]) ):
        
        p_pbar  = tab1[i,0]*1e-3

        factor = np.sqrt(t.fMass_proton*t.fMass_proton+p_pbar*p_pbar)/p_pbar/p_pbar*1e3
        
        CS   =  tab1[i,15]*factor
        stat =  tab1[i,16]*factor
        
        if CS!=CS:
            continue
        
        sys  =  0
        scale = err_scale
        
        
        eta = t.theta_to_eta(theta)
        
        s, pT, xR = t.convert_LAB_to_CM__from__s_ppbar_eta(s, p_pbar, eta)
        
        line = '   '
        line += str(np.sqrt(s)  ).ljust(20)+' '
        line += str(pT          ).ljust(20)+' '
        line += str(xR          ).ljust(20)+' '
        line += str(p_pbar      ).ljust(20)+' '
        line += str(theta       ).ljust(20)+' '
        line += str(CS          ).ljust(20)+' '
        line += str(stat        ).ljust(20)+' '
        line += str(sys         ).ljust(20)+' '
        line += str(scale       ).ljust(20)+'\n'
        f.write( line )
    
    f.write( '#\n#End Tab '+str(tab)+'\n#' )


fScale=0.15

def write(data_dir, res_dir):
    global f, fScale
    f = open(res_dir+'/converted_allaby.txt', 'w')

    line = '# Allaby et al, DOI:10.17182/hepdata.1345\n'
    f.write( line )
    line = '# p + p -> pbar + X \n#'
    f.write( line )

    readTab( data_dir + '/allaby_p4.txt',  6.15*6.15, 4.5, fScale, 1   )
    readTab( data_dir + '/allaby_p6.txt',  6.15*6.15, 6.0, fScale, 2   )
    readTab( data_dir + '/allaby_p8.txt',  6.15*6.15, 8.0, fScale, 3   )
    readTab( data_dir + '/allaby_p10.txt', 6.15*6.15, 10., fScale, 4   )
    readTab( data_dir + '/allaby_p11.txt', 6.15*6.15, 11., fScale, 5   )
    readTab( data_dir + '/allaby_p12.txt', 6.15*6.15, 12., fScale, 6   )
    readTab( data_dir + '/allaby_p13.txt', 6.15*6.15, 13., fScale, 7   )
    readTab( data_dir + '/allaby_p14.txt', 6.15*6.15, 14., fScale, 8   )


    #readTab_t( data_dir + '/allaby_t12.txt', 6.15*6.15, 0.0125, fScale, 11   )

    f.close()















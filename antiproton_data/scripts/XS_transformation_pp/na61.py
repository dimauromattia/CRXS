#! /usr/bin/env python

import glob
import math
import numpy as np

import XS_transformations as t


def readTab( file, pLab, err_scale, tab=1 ):
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



    tab1 = np.genfromtxt( file )
    
    for i in range( len(tab1[:,0]) ):
        
        y = tab1[i,0]
        
        #print( y )
        
        for j in range(10):
            s   = t.s__from_pLAB(pLab)
            pT = 0.05 + j*0.1                   # FIXME: Check that this corresponds to the plot in the paper!
            factor = 1./6.2831853072/pT*t.tot_pp__diMauro( s )
            CS = tab1[i,1+j*3] * factor
            e1 = tab1[i,2+j*3] * factor
            e2 = tab1[i,3+j*3] * factor
            
            if CS!=CS:
                continue
            
            stat =  np.sqrt(e1*e1+e2*e2)
            sys  =  0.1*CS
            scale = err_scale
        
        
            eta = t.y_to_eta(y, pT)
            
            sig = np.sign(eta)
            if sig == 0:
                sig = 1
            
            ppbar = pT*np.cosh(eta)
            
            xR  = t.x_R__from_s_ppbar(s, ppbar)*sig
            
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



def readTab_scanfromplot( file, err_scale, tab=1 ):
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



    tab1 = np.genfromtxt( file )

    for i in range( len(tab1[:,0]) ):

        s   = t.s__from_pLAB(tab1[i,0])
        pT = tab1[i,2]
        y = tab1[i,1]

        factor = 1./6.2831853072/pT*t.tot_pp__diMauro( s )
        CS = tab1[i,3] * factor
        e1 = tab1[i,4] * factor
        e2 = tab1[i,5] * factor

        if CS!=CS:
            continue

        stat =  e1  # np.sqrt(e1*e1+e2*e2)
        sys  =  e2  # 0.1*CS
        scale = err_scale


        eta = t.y_to_eta(y, pT)

        ppbar = pT*np.cosh(eta)

        sig = np.sign(eta)
        if sig == 0:
            sig = 1

        xR  = t.x_R__from_s_ppbar(s, ppbar)*sig

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






fScale=0.05

def write(data_dir, res_dir):
    global f, fScale
    f = open(res_dir+'/converted_na61.txt', 'w')

    line = '# NA61 et al,  arXiv:1705.02467 \n'
    f.write( line )
    line = '# p + p -> pbar + X \n #'
    f.write( line )
    
    readTab_scanfromplot  ( data_dir+'/NA61_31GeV.txt',                     fScale  )
    readTab_scanfromplot  ( data_dir+'/NA61_40GeV.txt',                     fScale  )
    readTab_scanfromplot  ( data_dir+'/NA61_80GeV.txt',                     fScale  )
    readTab               ( data_dir+'/na61__pp_to_pbar__158GeV.txt', 158,  fScale  )

    f.close()















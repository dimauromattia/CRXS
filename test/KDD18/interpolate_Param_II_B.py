#
#       Import numpy and interpolation functions from scipy
#
import  numpy   as      np
import  os
from    scipy   import  interpolate

#***************************************************************#
#       Read table and define function for interpolation        #
#                                                               #
#                    ** Param II B **                           #
#                                                               #
#***************************************************************#

os.system('head -n 20 KDD18/XS_table_Param_II_B.dat')

#
#   Read table
#
param_II_B    = np.genfromtxt('KDD18/XS_table_Param_II_B.dat')

#
#   Define number of columns for the project kinetic energy per nucleon (Tn)
#   and the kinetic energy of the produced antiproton (Tpbar).
#
n_Tn    = 30*7+1
n_Tpbar = 30*5+1

#
#   Extract the list values of Tn and Tbar (in log)
#
log_Tn    = np.log10(param_II_B[::n_Tpbar,0])
log_Tpbar = np.log10(param_II_B[:n_Tpbar, 1])

#
#   Transform the 1D table of cross sections (XS) into a 2D table.
#   For example the cross section at Tn[i], Tpbar[j] is now
#   given by XS_pp[i,j]. (Add +1e-90 to avoid log(0).)
#
log_XS_pp   = np.log( np.array(np.split(param_II_B[:,2],n_Tn))+1e-90 )
log_XS_pHe  = np.log( np.array(np.split(param_II_B[:,3],n_Tn))+1e-90 )
log_XS_Hep  = np.log( np.array(np.split(param_II_B[:,8],n_Tn))+1e-90 )
log_XS_HeHe = np.log( np.array(np.split(param_II_B[:,9],n_Tn))+1e-90 )

#
#   Create interpolation fuction on the tables. Interpolatoin is done in log of all
#   quantities to improve interpolation quality of the strongly variing XS.
#
f_log_XS_pp   = interpolate.interp2d(log_Tn, log_Tpbar, np.transpose(log_XS_pp  ), kind='cubic')
f_log_XS_pHe  = interpolate.interp2d(log_Tn, log_Tpbar, np.transpose(log_XS_pHe ), kind='cubic')
f_log_XS_Hep  = interpolate.interp2d(log_Tn, log_Tpbar, np.transpose(log_XS_Hep ), kind='cubic')
f_log_XS_HeHe = interpolate.interp2d(log_Tn, log_Tpbar, np.transpose(log_XS_HeHe), kind='cubic')

##  Energy differential production cross section of antiprotons in the Galaxy by proton-proton collisions.
#
#   @param  Tp      Kinetic energy of the proton.
#   @param  Tpbar   Kinetic energy of the produced antiproton.
#
#   \f$ d\sigma/dT_{\bar{p}} (p+p\rightarrow\bar{p}+X \f$. Antineutrons and hyperon-induced production is included.
#
def dT_pp_pbar(Tp, Tpbar):
    global f_log_XS_pp
    if np.amin(Tp)<0 or np.amax(Tp)>1e7 or np.amin(Tpbar)<0.1 or np.amax(Tpbar)>1e4:
        print( 'Warning: In function dT_pp_pbar(Tp, Tpbar) on of the input energies is out of range.\nT_p='+str(THe_n)+'\n T_pbar='+str(Tpbar)+'\nReturn 0.')
        return 0
    xs = f_log_XS_pp(np.log10(Tp), np.log10(Tpbar))
    return np.exp(xs)

##  Energy differential production cross section of antiprotons in the Galaxy by proton-helium collisions.
#
#   @param  Tp      Kinetic energy of the proton.
#   @param  Tpbar   Kinetic energy of the produced antiproton.
#
#   \f$ d\sigma/dT_{\bar{p}} (p+He\rightarrow\bar{p}+X \f$. Antineutrons and hyperon-induced production is included.
#
def dT_pHe_pbar(Tp, Tpbar):
    global f_log_XS_pHe
    if np.amin(Tp)<0 or np.amax(Tp)>1e7 or np.amin(Tpbar)<0.1 or np.amax(Tpbar)>1e4:
        print( 'Warning: In function dT_pHe_pbar(Tp, Tpbar) on of the input energies is out of range.\nT_p='+str(THe_n)+'\n T_pbar='+str(Tpbar)+'\nReturn 0.' )
        return 0
    xs = f_log_XS_pHe(np.log10(Tp), np.log10(Tpbar))
    return np.exp(xs)

##  Energy differential production cross section of antiprotons in the Galaxy by helium-helium collisions.
#
#   @param  T_He/n  Kinetic energy per nucleon of the incident helium.
#   @param  Tpbar   Kinetic energy of the produced antiproton.
#
#   \f$ d\sigma/dT_{\bar{p}} (He+p\rightarrow\bar{p}+X \f$. Antineutrons and hyperon-induced production is included.
#
def dT_Hep_pbar(THe_n, Tpbar):
    global f_log_XS_Hep
    if np.amin(THe_n)<0 or np.amax(THe_n)>1e7 or np.amin(Tpbar)<0.1 or np.amax(Tpbar)>1e4:
        print( 'Warning: In function dT_Hep_pbar(Tp, Tpbar) on of the input energies is out of range.\nT_He/n='+str(THe_n)+'\n T_pbar='+str(Tpbar)+'\nReturn 0.' )
        return 0
    xs = f_log_XS_Hep(np.log10(THe_n), np.log10(Tpbar))
    return np.exp(xs)

##  Energy differential production cross section of antiprotons in the Galaxy by proton-proton collisions.
#
#   @param  T_He/n  Kinetic energy per nucleon of the incident helium.
#   @param  Tpbar   Kinetic energy of the produced antiproton.
#
#   \f$ d\sigma/dT_{\bar{p}} (He+He\rightarrow\bar{p}+X \f$. Antineutrons and hyperon-induced production is included.
#
def dT_HeHe_pbar(THe_n, Tpbar):
    global f_log_XS_HeHe
    if np.amin(THe_n)<0 or np.amax(THe_n)>1e7 or np.amin(Tpbar)<0.1 or np.amax(Tpbar)>1e4:
        print( 'Warning: In function dT_HeHe_pbar(Tp, Tpbar) on of the input energies is out of range.\nT_He/n='+str(THe_n)+'\n T_pbar='+str(Tpbar)+'\nReturn 0.' )
        return 0
    xs = f_log_XS_HeHe(np.log10(THe_n), np.log10(Tpbar))
    return np.exp(xs)

################################################################################################
#
#   Now you can simply use the the functions dT_.._pbar. They accept array-like arguments
#
################################################################################################

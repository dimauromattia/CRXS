#! /usr/bin/env python3

import numpy                         as np

from scipy import integrate

#
#   Use matplotlib to produce plots
#
import matplotlib.pyplot    as       plt
import matplotlib           as       mpl
import matplotlib.colors    as       colors
from   matplotlib           import   rc

#
#   Function to set some nice layout
#
print_size=10
label_size=35
def plot_1D( xlabel, ylabel, xscale='log', yscale='log', sizex=1.3, sizey=1.0):
    global print_size, label_size
    plt.close('all')
    font_props = {"size":label_size}
    rc("font", **font_props)
    mpl.rcParams['axes.linewidth'] = 2
    fig     = plt.figure(figsize=(print_size*sizex, print_size*sizey))
    plot    = plt.subplot2grid((1, 1), (0, 0), rowspan=1)
    plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)
    plot.set_xlabel ( xlabel )
    plot.set_ylabel ( ylabel )
    plot.set_xscale ( xscale )
    plot.set_yscale ( yscale )
    plot.tick_params('both', length=20, width=2, which='major', top=True, right=True, direction='in', pad=10)
    plot.tick_params('both', length=10, width=1, which='minor', top=True, right=True, direction='in', pad=10)
    plot.grid(b=True, which='major', alpha=0.1, linestyle='-', linewidth=2)
    return plot, fig


from scipy import interpolate

dat_pp_winkler    = np.genfromtxt('pp_cross.dat.txt')
winkler_pp_Tp     =  dat_pp_winkler[::251,0]
winkler_pp_Tpbar  =  dat_pp_winkler[:251, 1]
winkler_pp_xs     =  dat_pp_winkler[:,2].reshape( len(winkler_pp_Tp), 251 )
_f_winkler_pp      = interpolate.interp2d( np.log(winkler_pp_Tpbar), np.log(winkler_pp_Tp), np.log(np.maximum(winkler_pp_xs, 1e-99)) )

def f_winkler_pp(  Tpbar, Tp ):
    return np.exp(_f_winkler_pp( np.log(Tpbar), np.log(Tp) ) )


dat_pHe_winkler    = np.genfromtxt('pHe_cross.dat.txt')
winkler_pHe_Tp     =  dat_pHe_winkler[::251,0]
winkler_pHe_Tpbar  =  dat_pHe_winkler[:251, 1]
winkler_pHe_xs     =  dat_pHe_winkler[:,2].reshape( len(winkler_pHe_Tp), 251 )
_f_winkler_pHe     = interpolate.interp2d( np.log(winkler_pHe_Tpbar), np.log(winkler_pHe_Tp), np.log(np.maximum(winkler_pHe_xs, 1e-99)) )

def f_winkler_pHe( Tpbar, Tp ):
    return np.exp(_f_winkler_pHe( np.log(Tpbar), np.log(Tp) ) )


dat_Hep_winkler    = np.genfromtxt('Hep_cross.dat.txt')
winkler_Hep_Tp     =  dat_Hep_winkler[::251,0]
winkler_Hep_Tpbar  =  dat_Hep_winkler[:251, 1]
winkler_Hep_xs     =  dat_Hep_winkler[:,2].reshape( len(winkler_Hep_Tp), 251 )
_f_winkler_Hep     = interpolate.interp2d( np.log(winkler_Hep_Tpbar), np.log(winkler_Hep_Tp), np.log(np.maximum(winkler_Hep_xs, 1e-99)) )

def f_winkler_Hep( Tpbar, Tp ):
    return np.exp(_f_winkler_Hep( np.log(Tpbar), np.log(Tp) ) )


dat_HeHe_winkler    = np.genfromtxt('HeHe_cross.dat.txt')
winkler_HeHe_Tp     =  dat_HeHe_winkler[::251,0]
winkler_HeHe_Tpbar  =  dat_HeHe_winkler[:251, 1]
winkler_HeHe_xs     =  dat_HeHe_winkler[:,2].reshape( len(winkler_HeHe_Tp), 251 )
_f_winkler_HeHe     = interpolate.interp2d( np.log(winkler_HeHe_Tpbar), np.log(winkler_HeHe_Tp), np.log(np.maximum(winkler_HeHe_xs, 1e-99)) )

def f_winkler_HeHe( Tpbar, Tp ):
    return np.exp(_f_winkler_HeHe( np.log(Tpbar), np.log(Tp) ) )




import CRXS.XS_wrapper               as XS
import KDD18.interpolate_Param_II_B  as kdd18





##########################################################################
#   Plot a profile of the cross sections d sigma/d T (T=200 GeV, Tpbar)
##########################################################################

#
#   p + p -> pbar + X
#

dE_AA_pbar_LAB_incNbarAndHyperon = np.vectorize(XS.dE_AA_pbar_LAB_incNbarAndHyperon)
dE_AA_pbar_LAB                   = np.vectorize(XS.dE_AA_pbar_LAB                  )


list_T      = [20, 100, 1000, 5000]

list_Tpbar  = [1, 10, 100, 500]


#T      = 200
name = [ 'pp',      'pHe',     'Hep',     'HeHe'     ]
nuc  = [ [1,0,1,0], [1,0,4,2], [4,2,1,0] , [4,2,4,2] ]


m_proton = 0.938

def deltaHyperon(s):
    
    C1 = 0.31
    C2 = 0.30
    C3 = 21316.
    C4 = 0.9

    factor   = 0.81

    hyperon = C1 + C2/(1+pow(C3/s,C4))
    hyperon       *= factor
    return hyperon

def deltaIsospin(s):

    C14 = 0.114
    C15 = 20736.
    C16 = 0.51

    return C14/(1+pow(s/C15, C16));

for i, n in enumerate(name):

    A_p = nuc[i][0]
    N_p = nuc[i][1]
    
    A_t = nuc[i][2]
    N_t = nuc[i][3]
    
    
    print(n)
    print(A_t)
    print(N_t)

    for i, T in enumerate(list_T):
        
        s = 4*m_proton**2 + 2 * T * m_proton;
        
        fac =  2 + 2 * deltaHyperon(s) + deltaIsospin(s)
        
        print( 'T:     %f \n  s:   %f \n  fac: %f ' % (T, s, fac) )
        
        
        Tpbar  = np.power(10, np.arange(-1,np.log10(T),0.05) )

        if A_t == 1 and A_p == 1:
            xs_kdd18    = kdd18.dT_pp_pbar                    (T, Tpbar)[:,0]
            xs_winker_tab = f_winkler_pp(Tpbar, T) *1e-31
        if A_t == 4:
            xs_kdd18    = kdd18.dT_pHe_pbar                   (T, Tpbar)[:,0]
            xs_winker_tab = f_winkler_pHe(Tpbar, T) *1e-31

        if A_p == 4:
            xs_kdd18      = kdd18.dT_Hep_pbar                   (T, Tpbar)[:,0]
            xs_winker_tab = f_winkler_Hep(Tpbar, T) *1e-31

        if A_p == 4 and A_t == 4:
            xs_kdd18      = kdd18.dT_HeHe_pbar                   (T, Tpbar)[:,0]
            xs_winker_tab = f_winkler_HeHe(Tpbar, T) *1e-31

        
        xs_wrapper           = dE_AA_pbar_LAB_incNbarAndHyperon   (T, Tpbar, A_projectile=A_p, N_projectile=N_p, A_target=A_t, N_target=N_t, parametrization='KORSMEIER_II') * 1e-31   # factor 1e-31, conversion from mbarn to m^2
        xs_wrapper_III       = dE_AA_pbar_LAB_incNbarAndHyperon   (T, Tpbar, A_projectile=A_p, N_projectile=N_p, A_target=A_t, N_target=N_t, parametrization='KORSMEIER_III')* 1e-31   # factor 1e-31, conversion from mbarn to m^2
        xs_wrapper_winkler   = dE_AA_pbar_LAB_incNbarAndHyperon   (T, Tpbar, A_projectile=A_p, N_projectile=N_p, A_target=A_t, N_target=N_t, parametrization='WINKLER'     ) * 1e-31   # factor 1e-31, conversion from mbarn to m^2
        xs_wrapper_winkler_II= dE_AA_pbar_LAB_incNbarAndHyperon   (T, Tpbar, A_projectile=A_p, N_projectile=N_p, A_target=A_t, N_target=N_t, parametrization='WINKLER_II'  ) * 1e-31   # factor 1e-31, conversion from mbarn to m^2
        xs_wrapper_noHyp     = dE_AA_pbar_LAB                     (T, Tpbar, A_projectile=A_p, N_projectile=N_p, A_target=A_t, N_target=N_t, parametrization='KORSMEIER_II') * 1e-31   # factor 1e-31, conversion from mbarn to m^2

        
        
        plot, fig = plot_1D ( r'$T_{\bar{p}}\;[\mathrm{GeV}]$', r'$d\sigma/dT_{\bar{p}} \;[\mathrm{m^2GeV^{-1}}]$', 'log', 'log' )
        plt.subplots_adjust(left=0.25, right=0.9, top=0.9, bottom=0.15)
        plot.plot(Tpbar, xs_kdd18,             c='#aa0000', lw=3,                         label='KDD18 (table)' )
        plot.plot(Tpbar, xs_wrapper,           c='#0000aa', lw=3, dashes=(5,5),           label='CRXS         ' )
        plot.plot(Tpbar, xs_wrapper_III,       c='#000000', lw=3, dashes=(1,1),           label='CRXS (with square)' )
        #plot.plot(Tpbar, xs_wrapper_noHyp*fac, c='#00aa00', lw=3, dashes=(2,2),           label=r'CRXS, no Iso, no Hyp $\times$ %.3f    ' % fac )

        plot.plot(Tpbar, xs_wrapper_winkler_II,c='#666666', lw=3, dashes=(5,2,2,2,2,2),   label='CRXS, Winkler II      ' )
        plot.plot(Tpbar, xs_wrapper_winkler,   c='#00aaaa', lw=3, dashes=(10,2),          label='CRXS, Winkler         ' )
        plot.plot(Tpbar, xs_winker_tab,        c='#aaaa00', lw=3, dashes=(5,2,2,2),       label='Winkler (table)       ' )

        plot.legend( loc='upper left', bbox_to_anchor=(0.03, 0.97), frameon=False, ncol=1, fontsize=0.5*label_size )
        plot.text(0.7, 2e-34, r'$T_p=%.2f\;\mathrm{GeV}$'%T)
        plot.set_xlim( ( 1e-1,  1e2  ) )
        plot.set_ylim( ( 5e-35, 8e-32) )
        plt.savefig('profile_'+n+'_T_%i.png'%T)

        
        plot, fig = plot_1D ( r'$T_{\bar{p}}\;[\mathrm{GeV}]$', r'$d\sigma/dT_{\bar{p}} \;[\mathrm{m^2GeV^{-1}}]$', 'log', 'linear' )
        plt.subplots_adjust(left=0.25, right=0.9, top=0.9, bottom=0.15)

        plot.plot(Tpbar, xs_kdd18/xs_wrapper,                c='#aa0000', lw=3,                       label='KDD18 (table/CRXS)'    )
        plot.plot(Tpbar, xs_wrapper_III/xs_wrapper,          c='#000000', lw=3, dashes=(1,1),         label='CRXS  (square/no square)'    )
        plot.plot(Tpbar, xs_kdd18/(xs_wrapper_noHyp*fac),    c='#00aa00', lw=3, dashes=(2,2),         label='KDD18 (table/ (CRXS no hyp+iso x %.3f )'  % fac   )
        plot.plot(Tpbar, xs_winker_tab/xs_wrapper_winkler,   c='#aaaa00', lw=3, dashes=(5,2,2,2),     label='Winkler (table/CRXS) ' )
        plot.plot(Tpbar, xs_winker_tab/xs_wrapper_winkler_II,c='#666666', lw=3, dashes=(5,2,2,2,2,2), label='Winkler II (table/CRXS)'  )

        
        plot.legend( loc='upper left', bbox_to_anchor=(0.03, 0.97), frameon=False, ncol=1, fontsize=0.5*label_size )
        plot.text(0.7, 0.85, r'$T_p=%.2f\;\mathrm{GeV}$'%T)
        plot.set_xlim( ( 1e-1,  T  ) )
        plot.set_ylim( ( 0.8, 1.2) )
        plt.savefig('ratio_'+n+'_T_%i.png'%T)


    for i, Tpbar in enumerate(list_Tpbar):
        
        
        T      = np.power(10, np.arange(0,5,0.01) )
        
        if A_t == 1 and A_p==1:
            #print(kdd18.dT_pp_pbar(T, Tpbar))
            xs_kdd18    = kdd18.dT_pp_pbar                    (T, Tpbar)
            xs_winker_tab = f_winkler_pp(Tpbar, T)  * 1e-31
        if A_t == 4:
            xs_kdd18    = kdd18.dT_pHe_pbar                   (T, Tpbar)
            xs_winker_tab = f_winkler_pHe(Tpbar, T) * 1e-31
        if A_p == 4:
            xs_kdd18      = kdd18.dT_Hep_pbar                 (T, Tpbar)
            xs_winker_tab = f_winkler_Hep(Tpbar, T) *1e-31

        if A_p == 4 and A_t == 4:
            xs_kdd18      = kdd18.dT_HeHe_pbar                 (T, Tpbar)
            xs_winker_tab = f_winkler_HeHe(Tpbar, T) *1e-31


        xs_winker_tab = xs_winker_tab[:,0]
        
        
        xs_wrapper           = dE_AA_pbar_LAB_incNbarAndHyperon   (T, Tpbar, A_projectile=A_p, N_projectile=N_p, A_target=A_t, N_target=N_t, parametrization='KORSMEIER_II') * 1e-31   # factor 1e-31, conversion from mbarn to m^2
        xs_wrapper_winkler   = dE_AA_pbar_LAB_incNbarAndHyperon   (T, Tpbar, A_projectile=A_p, N_projectile=N_p, A_target=A_t, N_target=N_t, parametrization='WINKLER'     ) * 1e-31   # factor 1e-31, conversion from mbarn to m^2
        xs_wrapper_winkler_II= dE_AA_pbar_LAB_incNbarAndHyperon   (T, Tpbar, A_projectile=A_p, N_projectile=N_p, A_target=A_t, N_target=N_t, parametrization='WINKLER_II'  ) * 1e-31   # factor 1e-31, conversion from mbarn to m^2
        xs_wrapper_noHyp     = dE_AA_pbar_LAB                     (T, Tpbar, A_projectile=A_p, N_projectile=N_p, A_target=A_t, N_target=N_t, parametrization='KORSMEIER_II') * 1e-31   # factor 1e-31, conversion from mbarn to m^2
        

        plot, fig = plot_1D ( r'$T\;[\mathrm{GeV}]$', r'$T^{-1.7} \;\; d\sigma/dT_{\bar{p}} \;[\mathrm{m^2GeV^{-1+1.7}}]$', 'log', 'log' )
        plt.subplots_adjust(left=0.25, right=0.9, top=0.9, bottom=0.15)
        plot.plot(T, T**-1.7 * xs_kdd18,             c='#aa0000', lw=3,                          label='KDD18 (table)' )
        plot.plot(T, T**-1.7 * xs_wrapper,           c='#0000aa', lw=3, dashes=(5,5),            label='CRXS         ' )
        plot.plot(T, T**-1.7 * xs_wrapper_noHyp,     c='#00aa00', lw=3, dashes=(2,2),            label='CRXS, no Iso, no Hyp  ' )
        plot.plot(T, T**-1.7 * xs_wrapper_winkler,   c='#00aaaa', lw=3, dashes=(10,2),           label='CRXS, Winkler         ' )
        plot.plot(T, T**-1.7 * xs_wrapper_winkler_II,c='#666666', lw=3, dashes=(5,2,2,2,2,2),    label='CRXS, Winkler II      ' )
        plot.plot(T, T**-1.7 * xs_winker_tab,        c='#aaaa00', lw=3, dashes=(5,2,2,2),        label='Winkler (table)       ' )
        
        plot.legend( loc='upper left', bbox_to_anchor=(0.03, 0.97), frameon=False, ncol=1, fontsize=0.5*label_size )

        max = np.amax(T**-1.7*xs_wrapper)

        plot.text(1e2, max/1e2*5, r'$T_{\bar{p}}=%.2f\;\mathrm{GeV}$'%Tpbar)
        plot.set_xlim( ( 1e+0,  1e5  ) )
        plot.set_ylim( ( max/1e2, max*5) )
        plt.savefig('profile_'+n+'_Tpbar_%i.png'%Tpbar)
        
        
        plot, fig = plot_1D ( r'$T\;[\mathrm{GeV}]$', r'$d\sigma/dT_{\bar{p}} \;[\mathrm{m^2GeV^{-1}}]$', 'log', 'linear' )
        plt.subplots_adjust(left=0.25, right=0.9, top=0.9, bottom=0.15)

        plot.plot(T, xs_kdd18/xs_wrapper,                c='#aa0000', lw=3,                       label='KDD18 (table/CRXS)'  )
        plot.plot(T, xs_winker_tab/xs_wrapper_winkler,   c='#aaaa00', lw=3, dashes=(5,2,2,2),     label='Winkler (table/CRXS)'  )
        plot.plot(T, xs_winker_tab/xs_wrapper_winkler_II,c='#666666', lw=3, dashes=(5,2,2,2,2,2), label='Winkler II (table/CRXS)'  )

        plot.legend( loc='upper left', bbox_to_anchor=(0.03, 0.97), frameon=False, ncol=1, fontsize=0.5*label_size )
        plot.text(100, 0.85, r'$T_\bar{p}=%.2f\;\mathrm{GeV}$'%Tpbar)
        plot.set_xlim( ( 1e+0,  1e5  ) )
        plot.set_ylim( ( 0.8, 1.2) )
        plt.savefig('ratio_'+n+'_Tpbar_%i.png'%Tpbar)



# fake sourceterm

    if A_t == 1 and A_p==1:
        _xs_kdd18    = kdd18.dT_pp_pbar
        xs_winker    = f_winkler_pp
    if A_t == 4:
        _xs_kdd18    = kdd18.dT_pHe_pbar
        xs_winker    = f_winkler_pHe
    if A_p == 4:
        _xs_kdd18    = kdd18.dT_Hep_pbar
        xs_winker    = f_winkler_Hep
    if A_p == 4 and A_t == 4:
        print('HelloWorld!!!!!')
        _xs_kdd18    = kdd18.dT_HeHe_pbar
        xs_winker    = f_winkler_HeHe
    


    def f_xs_kdd18(T, Tpbar):
        return _xs_kdd18(T, Tpbar)[0]


    T_pbar_vec = np.power(10, np.arange(-0.5, 3, 0.2))

    st_wrapper            = []
    st_wrapper_III        = []
    st_wrapper_winkler    = []
    st_wrapper_winkler_II = []
    st_winker_tab         = []
    st_kdd18              = []

    for Tpbar in T_pbar_vec:
        #        print(Tpbar)
        st = integrate.quad(   lambda x: np.exp(x)**-1.7 * dE_AA_pbar_LAB_incNbarAndHyperon   (np.exp(x), Tpbar, A_projectile=A_p, N_projectile=N_p, A_target=A_t, N_target=N_t, parametrization='KORSMEIER_II') * 1e-31, np.log(Tpbar), np.log(1e5)  )
        st_wrapper.append           (st[0])
        st = integrate.quad(   lambda x: np.exp(x)**-1.7 * dE_AA_pbar_LAB_incNbarAndHyperon   (np.exp(x), Tpbar, A_projectile=A_p, N_projectile=N_p, A_target=A_t, N_target=N_t, parametrization='KORSMEIER_III')* 1e-31, np.log(Tpbar), np.log(1e5)  )
        st_wrapper_III.append       (st[0])
        st = integrate.quad(   lambda x: np.exp(x)**-1.7 * dE_AA_pbar_LAB_incNbarAndHyperon   (np.exp(x), Tpbar, A_projectile=A_p, N_projectile=N_p, A_target=A_t, N_target=N_t, parametrization='WINKLER'     ) * 1e-31, np.log(Tpbar), np.log(1e5)  )
        st_wrapper_winkler.append   (st[0])
        st = integrate.quad(   lambda x: np.exp(x)**-1.7 * dE_AA_pbar_LAB_incNbarAndHyperon   (np.exp(x), Tpbar, A_projectile=A_p, N_projectile=N_p, A_target=A_t, N_target=N_t, parametrization='WINKLER_II'  ) * 1e-31, np.log(Tpbar), np.log(1e5)  )
        st_wrapper_winkler_II.append(st[0])
        st = integrate.quad(   lambda x: np.exp(x)**-1.7 * f_xs_kdd18                         (np.exp(x), Tpbar                                                                                            ),         np.log(Tpbar), np.log(1e5)  )
        st_kdd18.append             (st[0])
        st = integrate.quad(   lambda x: np.exp(x)**-1.7 * xs_winker                          (Tpbar    , np.exp(x)                                                                                        ) * 1e-31, np.log(Tpbar), np.log(1e5)  )
        st_winker_tab.append        (st[0])

    print(st_winker_tab)

    st_wrapper            = np.array(st_wrapper            )
    st_wrapper_III        = np.array(st_wrapper_III        )
    st_wrapper_winkler    = np.array(st_wrapper_winkler    )
    st_wrapper_winkler_II = np.array(st_wrapper_winkler_II )
    st_winker_tab         = np.array(st_winker_tab         )
    st_kdd18              = np.array(st_kdd18              )

    if A_t == 1 and A_p==1:
        st_wrapper_pp                = st_wrapper.copy()
        st_wrapper_winkler_pp        = st_wrapper_winkler.copy()
        st_wrapper_winkler_II_pp     = st_wrapper_winkler_II.copy()
        st_winker_tab_pp             = st_winker_tab.copy()
        st_kdd18_pp                  = st_kdd18.copy()

    if A_t == 1 and A_p==4:
        st_wrapper_Hep                = st_wrapper.copy()
        st_wrapper_winkler_Hep        = st_wrapper_winkler.copy()
        st_wrapper_winkler_II_Hep     = st_wrapper_winkler_II.copy()
        st_winker_tab_Hep             = st_winker_tab.copy()
        st_kdd18_Hep                  = st_kdd18.copy()


    if A_t == 4 and A_p==4:
        print('HelloWorld')
        st_wrapper_HeHe                = st_wrapper.copy()
        st_wrapper_winkler_HeHe        = st_wrapper_winkler.copy()
        st_wrapper_winkler_II_HeHe     = st_wrapper_winkler_II.copy()
        st_winker_tab_HeHe             = st_winker_tab.copy()
        st_kdd18_HeHe                  = st_kdd18.copy()
        print(st_winker_tab_HeHe)
    


    plot, fig = plot_1D ( r'$T_{\bar{p}}\;[\mathrm{GeV}]$', r'$ T^{2.7} \times $(fake source term)', 'log', 'log' )
    plt.subplots_adjust(left=0.25, right=0.9, top=0.9, bottom=0.15)

    plot.plot(T_pbar_vec, T_pbar_vec**2.7 * st_kdd18,             c='#aa0000', lw=3,                          label='KDD18 (table)' )
    plot.plot(T_pbar_vec, T_pbar_vec**2.7 * st_wrapper_III,       c='#000000', lw=3, dashes=(1,1),            label='CRXS (square)' )
    plot.plot(T_pbar_vec, T_pbar_vec**2.7 * st_wrapper,           c='#0000aa', lw=3, dashes=(5,5),            label='CRXS         ' )
    plot.plot(T_pbar_vec, T_pbar_vec**2.7 * st_wrapper_winkler,   c='#00aaaa', lw=3, dashes=(10,2),           label='CRXS, Winkler         ' )
    plot.plot(T_pbar_vec, T_pbar_vec**2.7 * st_wrapper_winkler_II,c='#666666', lw=3, dashes=(5,2,2,2,2,2),    label='CRXS, Winkler II      ' )
    plot.plot(T_pbar_vec, T_pbar_vec**2.7 * st_winker_tab,        c='#aaaa00', lw=3, dashes=(5,2,2,2),        label='Winkler (table)       ' )


    plot.legend( loc='upper left', bbox_to_anchor=(0.03, 0.97), frameon=False, ncol=1, fontsize=0.5*label_size )
    plt.savefig('fake_sourceterm_'+n+'.png')


    plot, fig = plot_1D ( r'$T_{\bar{p}}\;[\mathrm{GeV}]$', r'ratio', 'log', 'linear' )
    plt.subplots_adjust(left=0.25, right=0.9, top=0.9, bottom=0.15)

    plot.plot(T_pbar_vec, st_kdd18/st_wrapper,                     c='#aa0000', lw=3,                   label='KDD18  (table/CRXS)   '        )
    plot.plot(T_pbar_vec, st_winker_tab/st_wrapper_winkler,        c='#aaaa00', lw=3, dashes=(5,2,2,2), label='Winkler (table/CRXS)  '        )
    plot.plot(T_pbar_vec, st_winker_tab/st_wrapper_winkler_II,     c='#0000aa', lw=3, dashes=(5,5),     label='Winkler (table/CRXS with square)  ' )
    plot.plot(T_pbar_vec, st_wrapper_III/st_wrapper,               c='#000000', lw=3, dashes=(1,1),     label='CRXS (square/no square)  ' )

    plot.plot(T_pbar_vec, st_kdd18/st_winker_tab,                  c='#666666', lw=3, dashes=(10,2),    label='Table  (KDD18/Winkler)   '        )


    plot.legend( loc='upper left', bbox_to_anchor=(0.03, 0.97), frameon=False, ncol=1, fontsize=0.5*label_size )
    plot.set_ylim( ( 0.8, 1.2) )
    plt.savefig('fake_sourceterm_ratio_'+n+'.png')

    plot, fig = plot_1D ( r'$T_{\bar{p}}\;[\mathrm{GeV}]$', r'ratio', 'log', 'linear' )
    plt.subplots_adjust(left=0.25, right=0.9, top=0.9, bottom=0.15)
    
    plot.plot(T_pbar_vec, st_wrapper_winkler/st_wrapper_winkler_II, c='#0000aa', lw=3, dashes=(5,5), label='Winkler I/Winkler II  ' )
    
    plot.legend( loc='upper left', bbox_to_anchor=(0.03, 0.97), frameon=False, ncol=1, fontsize=0.5*label_size )
    plot.set_ylim( ( 0.8, 1.2) )
    plt.savefig('fake_sourceterm_ratio_Winkler_'+n+'.png')



plot, fig = plot_1D ( r'$T_{\bar{p}}\;[\mathrm{GeV}]$', r'ratio', 'log', 'linear' )
plt.subplots_adjust(left=0.25, right=0.9, top=0.9, bottom=0.15)

plot.plot(T_pbar_vec, st_wrapper_Hep/st_wrapper_pp,                        c='#aa0000', lw=3, dashes=(), label='Korsmeier II Hep / Korsmeier II pp ' )
plot.plot(T_pbar_vec, st_wrapper_winkler_II_Hep/st_wrapper_winkler_II_pp,  c='#0000aa', lw=3, dashes=(5,5), label='Winkler II Hep / Winkler II pp  ' )

plot.legend( loc='upper left', bbox_to_anchor=(0.03, 0.97), frameon=False, ncol=1, fontsize=0.5*label_size )
#plot.set_ylim( ( 0., 1.0) )
plt.savefig('fake_sourceterm_ratio_tests.png')


print(st_winker_tab_HeHe)
print(st_winker_tab_pp)
plot, fig = plot_1D ( r'$T_{\bar{p}}\;[\mathrm{GeV}]$', r'ratio', 'log', 'linear' )
plt.subplots_adjust(left=0.25, right=0.9, top=0.9, bottom=0.15)

plot.plot(T_pbar_vec, st_winker_tab_HeHe/st_winker_tab_pp,  c='#0000aa', lw=3,               label='Winkler Tab. HeHe / Winkler Tab. pp  ' )
plot.plot(T_pbar_vec, st_wrapper_HeHe/st_wrapper_pp,        c='#aa0000', lw=3, dashes=(5,5), label='Korsmeier II HeHe / Korsmeier II pp  ' )

plot.legend( loc='upper left', bbox_to_anchor=(0.03, 0.97), frameon=False, ncol=1, fontsize=0.5*label_size )
plot.set_ylim( ( 9, 16) )
plt.savefig('fake_sourceterm_ratio_testsHeHe.png')





#        plot.text(0.7, 0.85, r'$T_p=%.2f\;\mathrm{GeV}$'%T)
#        plot.set_xlim( ( 1e-1,  T  ) )
#        plot.set_ylim( ( 0.8, 1.2) )





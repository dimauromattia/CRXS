#! /usr/bin/env python3

import numpy                         as np
import CRXS.XS_wrapper               as XS


dE_AA_pbar_LAB_incNbarAndHyperon = np.vectorize(XS.dE_AA_pbar_LAB_incNbarAndHyperon)

s=''
s += '#\n'
s += '#        *************************    \n'
s += '#        * Sublementary material *    \n'
s += '#        *************************    \n'
s += '#    \n'
s += '#        * -------------------------------------------------------------------- *    \n'
s += '#        * Production cross sections of cosmic antiprotons                      *    \n'
s += '#        *            in the light of new data from NA49 and NA61 measurements  *    \n'
s += '#        * -------------------------------------------------------------------- *    \n'
s += '#    \n'
s += '#          Michael Korsmeier, Fiorenza Donato, and Mattia Di Mauro    \n'
s += '#         *-------------------------------------------------------*    \n'
s += '#    \n'
s += '#                    Published in: PRD    \n'
s += '#                    arXiv:        1802.03030 [astro-ph.HE]    \n'
s += '#    \n'
s += '#    \n'
s += '#          IF YOU USE THIS TABLE, PLEASE CITE THe PAPER.    \n'
s += '#    \n'
s += '#    \n'
s += '#          We provide the energy differential cross section $d\sigma_{ij}/dT$    \n'
s += '#          for cosmic-ray (CR) component i and intersellar medium (ISM)    \n'
s += '#          component j. Here j are protons (Z=1,A=1) and Helium (Z=2,A=4).    \n'
s += '#          The first column contains the kinetic energy per nucleon of the    \n'
s += '#          incident CR $T_{projectile}$ and the second column the kinetic    \n'
s += '#          energy of the antiproton $T_{ar{p}}, both in GeV.    \n'
s += '#          The following columns contain the energy differential cross    \n'
s += '#          secions of various possible CR isotopes in units of m^2/GeV.    \n'
s += '#          These are total cross section, i.e. including antineutrons and    \n'
s += '#          antihyperons.    \n'
s += '#    \n'
s += '#    \n'
s += '#               *-----------------------------------------------*    \n'
s += '#               | This corresponds to Param. II-B in the paper. |    \n'
s += '#               *-----------------------------------------------*    \n'
s += '#    \n'
s += '#          We recommend to use Param. II-B because of the better fit    \n'
s += '#          behaviour at high energies.    \n'
s += '#    \n'
s += '#    \n'
s += '#    \n'
s += '#**********************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************      \n'
s += '#  T_{proj}/n [GeV]      T_{pbar} [GeV]           Z= 1, A= 1 (ISM p )      Z= 1, A= 1 (ISM He)      Z= 1, A= 2 (ISM p )      Z= 1, A= 2 (ISM He)      Z= 2, A= 3 (ISM p )      Z= 2, A= 3 (ISM He)      Z= 2, A= 4 (ISM p )      Z= 2, A= 4 (ISM He)      Z= 6, A=12 (ISM p )      Z= 6, A=12 (ISM He)      Z= 6, A=13 (ISM p )      Z= 6, A=13 (ISM He)      Z= 7, A=14 (ISM p )      Z= 7, A=14 (ISM He)      Z= 7, A=15 (ISM p )      Z= 7, A=15 (ISM He)      Z= 8, A=16 (ISM p )      Z= 8, A=16 (ISM He)      Z= 8, A=17 (ISM p )      Z= 8, A=17 (ISM He)      Z= 8, A=18 (ISM p )      Z= 8, A=18 (ISM He)           \n'
s += '#**********************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************      \n'


list_Z1 =                                             [    1 ,                      1 ,                      1 ,                      1 ,                      2 ,                      2 ,                      2 ,                      2 ,                      6 ,                      6 ,                      6 ,                      6 ,                      7 ,                      7 ,                      7 ,                      7 ,                      8 ,                      8 ,                      8 ,                      8 ,                      8 ,                      8         ]
list_A1 =                                             [    1 ,                      1 ,                      2 ,                      2 ,                      3 ,                      3 ,                      4 ,                      4 ,                      12,                      12,                      13,                      13,                      14,                      14,                      15,                      15,                      16,                      16,                      17,                      17,                      18,                      18        ]
list_Z2 =                                             [    1 ,                      2 ,                      1 ,                      2 ,                      1 ,                      2 ,                      1 ,                      2 ,                      1 ,                      2 ,                      1 ,                      2 ,                      1 ,                      2 ,                      1 ,                      2 ,                      1 ,                      2 ,                      1 ,                      2 ,                      1 ,                      2         ]
list_A2 =                                             [    1 ,                      4 ,                      1 ,                      4 ,                      1 ,                      4 ,                      1 ,                      4 ,                      1 ,                      4 ,                      1 ,                      4 ,                      1 ,                      4 ,                      1 ,                      4 ,                      1 ,                      4 ,                      1 ,                      4 ,                      1 ,                      4         ]


vT_pbar = np.power( 10, np.arange( -1, 4+1./60, 1./30) )
vTn     = np.power( 10, np.arange(  0, 7+1./60, 1./30) )

#dE_AA_pbar_LAB_incNbarAndHyperon   (T, Tpbar, A_projectile=A_p, N_projectile=N_p, A_target=A_t, N_target=N_t, parametrization='KORSMEIER_II') * 1e-31   # factor 1e-31, conversion from mbarn to m^2

f = open('XS_table_Param_II_B.dat','w')
f.write(s)
s  = ''
for Tn in vTn:
    print(Tn)
    for T_pbar in vT_pbar:
        s += ' %-23.6e ' % Tn
        s += ' %-23.6e ' % T_pbar
        for i, Z1 in enumerate(list_Z1):
            Z2 = list_Z2[i]
            A1 = list_A1[i]
            A2 = list_A2[i]
            xs = dE_AA_pbar_LAB_incNbarAndHyperon   (Tn, T_pbar, A_projectile=A1, N_projectile=A1-Z1, A_target=A2, N_target=A2-Z2, parametrization='KORSMEIER_II') * 1e-31   # factor 1e-31, conversion from mbarn to m^2
            s += ' %-23.6e ' % xs
        s += '\n'
f.write(s)
f.close()




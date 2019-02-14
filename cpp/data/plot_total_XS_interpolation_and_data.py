#! /usr/bin/env python

import numpy as np

import MKastro.basic.PlotFunctions as pf
import matplotlib.pyplot           as plt

m_p = 0.938270
m_D = 1.875600

XS_pbarD_tot    = np.genfromtxt( 'pbardeut_total.dat.txt' , skip_header=11  )
XS_pbarp_ela    = np.genfromtxt( 'pbarp_elastic.dat.txt'  , skip_header=11  )
XS_pbarp_tot    = np.genfromtxt( 'pbarp_total.dat.txt'    , skip_header=11  )

XS_pbarD_nar    = np.genfromtxt( 'pbard_nar.dat.txt'  )

plot, fig = pf.new_plot(r'$T/n       \quad  \mathrm{[GeV/n]}$', r'$\sigma  \quad  \mathrm{[mbarn]}$', 'log', 'log')

T_D      = m_D/m_p * ( np.sqrt( XS_pbarD_tot[:,1]*XS_pbarD_tot[:,1] + m_p*m_p ) - m_p  )
T_p      =           ( np.sqrt( XS_pbarp_tot[:,1]*XS_pbarp_tot[:,1] + m_p*m_p ) - m_p  )
T_p_ela  =           ( np.sqrt( XS_pbarp_ela[:,1]*XS_pbarp_ela[:,1] + m_p*m_p ) - m_p  )
T_D_nar  = m_D/m_p * ( np.sqrt( XS_pbarD_nar[:,0]*XS_pbarD_nar[:,0] + m_p*m_p ) - m_p  )

plot.errorbar(T_D/2,   XS_pbarD_tot[:,4], xerr=0., yerr=np.sqrt( np.power(XS_pbarD_tot[:,5],2)+np.power(0.01*XS_pbarD_tot[:,4]*XS_pbarD_tot[:,6],2) ), fmt='o', markersize=3, label=r'$D\bar{p}$'      )
plot.errorbar(T_p,     XS_pbarp_tot[:,4], xerr=0., yerr=np.sqrt( np.power(XS_pbarp_tot[:,5],2)+np.power(0.01*XS_pbarp_tot[:,4]*XS_pbarp_tot[:,6],2) ), fmt='o', markersize=3, label=r'$p\bar{p}$'      )
plot.errorbar(T_p_ela, XS_pbarp_ela[:,4], xerr=0., yerr=np.sqrt( np.power(XS_pbarp_ela[:,5],2)+np.power(0.01*XS_pbarp_ela[:,4]*XS_pbarp_ela[:,6],2) ), fmt='o', markersize=3, label=r'$p\bar{p}$ (el)' )

plot.errorbar(T_D_nar[ 0: 1]/2,   XS_pbarD_nar[ 0: 1,1], xerr=0., yerr=np.sqrt( np.power(XS_pbarD_nar[ 0: 1,2],2) ), fmt='o', label=r'$D\bar{p}\rightarrow D\bar{p} \pi^+\pi^- \pi^0$'      )
plot.errorbar(T_D_nar[ 1:12]/2,   XS_pbarD_nar[ 1:12,1], xerr=0., yerr=np.sqrt( np.power(XS_pbarD_nar[ 1:12,2],2) ), fmt='o', label=r'$D\bar{p}\rightarrow D\bar{p} \pi^+\pi^-$'            )
plot.errorbar(T_D_nar[12:13]/2,   XS_pbarD_nar[12:13,1], xerr=0., yerr=np.sqrt( np.power(XS_pbarD_nar[12:13,2],2) ), fmt='o', label=r'$D\bar{p}\rightarrow D\bar{n} \pi^-\pi^- \pi^0$'      )
plot.errorbar(T_D_nar[13:14]/2,   XS_pbarD_nar[13:14,1], xerr=0., yerr=np.sqrt( np.power(XS_pbarD_nar[13:14,2],2) ), fmt='o', label=r'$D\bar{p}\rightarrow D\bar{\Lambda}^{--} \pi^+$'      )
plot.errorbar(T_D_nar[14:18]/2,   XS_pbarD_nar[14:18,1], xerr=0., yerr=np.sqrt( np.power(XS_pbarD_nar[14:18,2],2) ), fmt='o', label=r'$D\bar{p}\rightarrow D\bar{p} \pi^0$'                 )
plot.errorbar(T_D_nar[18:23]/2,   XS_pbarD_nar[18:23,1], xerr=0., yerr=np.sqrt( np.power(XS_pbarD_nar[18:23,2],2) ), fmt='o', label=r'$D\bar{p}\rightarrow D\bar{n} \pi^-$'                 )


p_XS_pbarD_tot    = np.genfromtxt( 'table_dpbar_tot.txt'   )
p_XS_pbarp_ela    = np.genfromtxt( 'table_ppbar_el.txt'    )
p_XS_pbarp_tot    = np.genfromtxt( 'table_ppbar_tot.txt'   )

plot.plot(p_XS_pbarD_tot[:,0], p_XS_pbarD_tot[:,1], marker='', linestyle='-' )
plot.plot(p_XS_pbarp_ela[:,0], p_XS_pbarp_ela[:,1], marker='', linestyle='-' )
plot.plot(p_XS_pbarp_tot[:,0], p_XS_pbarp_tot[:,1], marker='', linestyle='-' )

print(p_XS_pbarD_tot[:,0], p_XS_pbarD_tot[:,1])
print( p_XS_pbarp_tot[:,0] )
print( p_XS_pbarp_tot[:,1] )

dpbar_nar    = np.genfromtxt( 'table_dpbar_nar.txt'   )

plot.plot(dpbar_nar[:,0], dpbar_nar[:,1], marker='', linestyle='-' )
plot.fill_between(dpbar_nar[:,0], dpbar_nar[:,2], dpbar_nar[:,3], alpha=0.2)

handles, labels = plot.get_legend_handles_labels()

handles1 = handles[0:3]
labels1  = labels [0:3]

handles2 = handles[3:]
labels2  = labels [3:]


leg1 = plot.legend( handles1, labels1,  loc='upper center', bbox_to_anchor=(0.5, 0.95 ),  numpoints=1, ncol=3, frameon=False, fontsize=10)
leg2 = plot.legend( handles2, labels2,  loc='lower right',  bbox_to_anchor=(0.95, 0.05 ), numpoints=1, ncol=2, frameon=False, fontsize=10)


plot.add_artist(leg1)
plot.add_artist(leg2)

plot.set_ylim(1e-2, 1e3)
plt.savefig('Total_xs.png')


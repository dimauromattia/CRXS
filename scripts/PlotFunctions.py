import glob
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
from matplotlib import rc


#powerLawParam = np.array(  [[8631.08,   -1.47971,   0.123597,   0.0300829,   3.86131,   352.54,   0.613356],   [300.37,   -1.82257,   0.0956891,   -0.017557,   1.80785,   86.9094,   0.9792],   [8.27335,   -2.14854,   0.0721835,   0.0114956,   1.54107,   90.9645,   0.994855],   [12.0917,   -1.65672,   0.0740929,   -0.056748,   2.55923,   103.934,   1.01341],   [0.817645,   -3.82582,   0.287464,   0.187906,   4.04722,   352.54,   0.637847]]  )
#
#def brokenPowerLaw(E, p):
#    if E!=E:
#        return 0
#    if E<=0:
#        return 0
#    if E<=p[5]:
#        delta21 = p[2]-p[1]
#        ret =  p[0]  *  np.power(  E, -p[1] )   /   np.power( p[4], -p[1] ) * np.power(  np.power(E, 1./p[6])+np.power(p[4], 1./p[6]), -p[6]*delta21 )   / np.power( pow(2, p[6])*p[4], -delta21 )
#        #print ret
#        return  ret
#    return     brokenPowerLaw(p[5], p) * pow( E, -p[3] ) / pow( p[5], -p[3] )


def plotTest():
    print_size=10
    label_size=25
    

    plt.close('all')
    font_props = {"size":label_size}
    rc("font", **font_props)
    mpl.rcParams['axes.linewidth'] = 2
    mpl.rcParams['mathtext.fontset']='stixsans'
    fig = plt.figure(figsize=(print_size*1.3, print_size))
    
    upperPlt = plt.subplot2grid((1, 1), (0, 0), rowspan=1)
    plt.subplots_adjust(left=0.15, right=0.8, top=0.9, bottom=0.15)
    
    T = np.arange(-1,3.3+0.01,0.01)
    T = np.power(10,T)
    
    E = []
    for t in T:
        E.append(brokenPowerLaw(t, powerLawParam[0,:]))
    upperPlt.plot(T,E)
    upperPlt.set_xscale('log')
    upperPlt.set_yscale('log')
    
    plt.savefig("test.png")




def plot2D(file, file2='', Tmin_proj=1e0, Tmax_proj=1e6, Tmin_prod=1e-1, Tmax_prod=1e4, Zmin=-1, Zmax=-1, cmap_str='d', resfile='resFile', Zlabel=r'$\mathrm{d\sigma/dT_{\bar{p}}\quad[m^2/GeV]}$', Z1scale=1, Z2scale=1, zscale='', x_label=r'$\mathrm{T_{proj}/n\quad [GeV/n]}$',  y_label=r'$\mathrm{T_{\bar{p}}\quad [GeV]}$'): #, fluxfolding=''

    print_size=10
    label_size=30


    plt.close('all')
    font_props = {"size":label_size}
    rc("font", **font_props)
    mpl.rcParams['axes.linewidth'] = 2
    mpl.rcParams['mathtext.fontset']='stixsans'
    fig = plt.figure(figsize=(print_size*1.3, print_size))

    upperPlt = plt.subplot2grid((1, 1), (0, 0), rowspan=1)

    plt.subplots_adjust(left=0.15, right=0.75, top=0.9, bottom=0.15)


    data = np.genfromtxt( file )
    x = np.ones(len(data[:,0]))
    x[:-1] = data[1:,0 ]*np.sqrt(data[1,0]/data[2,0])
    x[ -1] = data[-1,0 ]/np.sqrt(data[1,0]/data[2,0])

    y = np.ones(len(data[0,:]))
    y[:-1] = data[0, 1:]*np.sqrt(data[0,1]/data[0,2])
    y[ -1] = data[0 ,-1]/np.sqrt(data[0,1]/data[0,2])


    z = np.transpose( data[1:,1:] )
    if file2!='':
        data2 = np.genfromtxt( file2 )
        z1 = np.transpose( data [1:,1:] )*Z1scale
        z2 = np.transpose( data2[1:,1:] )*Z2scale
        z = z2/z1
        for i in range(len(z[:,0])):
            for j in range(len(z[0,:])):
                if z[i,j] != z[i,j]:
                    z[i,j] = 1

#    if fluxfolding=='proton':
#        for i in range(len(z[:,0])):
#            for j in range(len(z[0,:])):
#                z[i,j] = z[i,j] * brokenPowerLaw(x[j], powerLawParam[0,:]) * np.power(x[j], -1.7) * np.power(y[i], 2.7)
#    if fluxfolding=='helium':
#        for i in range(len(z[:,0])):
#            for j in range(len(z[0,:])):
#                z[i,j] = z[i,j] * brokenPowerLaw(x[j], powerLawParam[1,:]) * np.power(x[j], -1.7) * np.power(y[i], 2.7)



    if cmap_str=='d':
        if file2=='':
            cmap_str = 'magma_r'
        else:
            cmap_str = 'seismic'
    cmap    = plt.get_cmap(cmap_str)

    if file2=='':
        if Zmin==-1 and Zmax==-1:
            max = np.amax(z)
            norm = colors.LogNorm(vmin=max*1e-4, vmax=max)
        else:
            norm = colors.LogNorm(vmin=Zmin, vmax=Zmax)
    else:
        if Zmin==-1 and Zmax==-1:
            norm = colors.Normalize(vmin=0, vmax=2)
        else:
            norm = colors.Normalize(vmin=Zmin, vmax=Zmax)
    if zscale=='log':
        if Zmin==-1 and Zmax==-1:
            norm = colors.LogNorm()
        else:
            norm = colors.LogNorm(vmin=Zmin, vmax=Zmax)
    if zscale=='linear':
        if Zmin==-1 and Zmax==-1:
            norm = colors.Normalize()
        else:
            norm = colors.Normalize(vmin=Zmin, vmax=Zmax)


#    print x
#    print y
    myPlt0 = upperPlt.pcolormesh(  x, y, z, norm=norm, cmap=cmap  )


    upperPlt.set_xscale('log')
    upperPlt.set_yscale('log')
    upperPlt.set_xlim(Tmin_proj, Tmax_proj)
    upperPlt.set_ylim(Tmin_prod, Tmax_prod)


    cbar_ax = fig.add_axes([0.8, 0.15, 0.02, 0.75])
    cbar = fig.colorbar(myPlt0, cax=cbar_ax, cmap=cmap)
    cbar.ax.set_ylabel(Zlabel, fontsize=label_size*1.4 )


    upperPlt.tick_params('both', length=20, width=2, which='major')
    upperPlt.tick_params('both', length=10, width=1, which='minor')

    #upperPlt.grid(b=True, which='major', alpha=0.5, linestyle='-', linewidth=2)

    upperPlt.tick_params(axis='both', pad=10)

    upperPlt.set_xlabel( x_label, fontsize=label_size*1.4  )
    upperPlt.set_ylabel( y_label, fontsize=label_size*1.4  )

    #upperPlt.legend(loc='upper left', fontsize=20)


    plt.savefig(resfile)


def profile(file, Tn, type='prod', file2='', Tn_min=1e-1, Tn_max=1e4, y_min=-1, y_max=-1, color='black', alpha=None, resfile='resFile', y_label=r'$\mathrm{d\sigma/dT_{\bar{p}}\quad[m^2/GeV]}$', x_label='', y_fac=1, y2_fac=1, y_scale='', label='', dashes=(), draw_on_top=False, legend=False, loc='', cs_power=0):
    
    T_prod = Tn
    
    print_size=10
    label_size=30
    if not draw_on_top:
        plt.close('all')
    font_props = {"size":label_size}
    rc("font", **font_props)
    mpl.rcParams['axes.linewidth'] = 2
    mpl.rcParams['mathtext.fontset']='stixsans'

    if not draw_on_top:
        fig = plt.figure(figsize=(print_size*1.3, print_size))
        upperPlt = plt.subplot2grid((1, 1), (0, 0), rowspan=1)
    else:
        global fig, upperPlt
    plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.15)
    data = np.genfromtxt( file )

    if not type=='prod':
        data = np.transpose(data)

    min  = 1e90;
    col  = 1
    for i in range(1,len(data[0,:])):
        diff = (data[0,i]-T_prod)*(data[0,i]-T_prod)
        if diff<min:
            min = diff
            col = i

    x_i = 1
    for i in range(1,len(data[:,0])):
        if data[i,0]<Tn_max:
            x_i = i

    if file2=='':
        if y_min==-1 and y_max==-1:
            y_min=5e-36
            y_max=5e-31
        if y_scale == '':
            y_scale = 'log'
    else:
        if y_min==-1 and y_max==-1:
            y_min=0
            y_max=2
        if y_scale == '':
            y_scale = 'linear'

    x = data[1:x_i,0 ]
    if file2=='':
        y = data[1:x_i,col]*y_fac
    else:
        data2 = np.genfromtxt( file2 )
        y1 = data [1:x_i,col]*y_fac
        y2 = data2[1:x_i,col]*y2_fac
        y = y2/y1
        for i in range(len(y)):
            if y[i] != y[i]:
                y[i] = 1

    
    myPlt0 = upperPlt.plot(  x, y*np.power(x, cs_power), color=color, label=label, lw=3, dashes=dashes, alpha=alpha  )

    upperPlt.set_xscale('log')
    upperPlt.set_yscale(y_scale)

    upperPlt.set_xlim(Tn_min, Tn_max)
    if y_min!='free':
        upperPlt.set_ylim(y_min, y_max)

    upperPlt.tick_params('both', length=20, width=2, which='major')
    upperPlt.tick_params('both', length=10, width=1, which='minor')

    upperPlt.grid(b=True, which='major', alpha=0.1, linestyle='-', linewidth=2)
    upperPlt.tick_params(axis='both', pad=10)
    if x_label=='':
        x_label = r'$\mathrm{T_{proj}/n \quad [GeV/n]}$'
        if not type=='prod':
            x_label = r'$\mathrm{T_{prod}/n \quad [GeV/n]}$'
    upperPlt.set_xlabel(x_label, fontsize=label_size*1.4 )
    upperPlt.set_ylabel(y_label, fontsize=label_size*1.4 )
    if legend:
        if type=='proj' or loc=='upper':
            upperPlt.legend(loc='upper right', bbox_to_anchor=(0.95, 0.95), frameon=False, fontsize=0.7*label_size)
        else:
            upperPlt.legend(loc='lower right', bbox_to_anchor=(0.95, 0.05), frameon=False, fontsize=0.7*label_size)
    plt.savefig(resfile)
    return upperPlt




def plot2D_easy(file, x_label=r'$\mathrm{T_{proj}/n\quad [GeV/n]}$', y_label=r'$\mathrm{T_{\bar{p}}\quad [GeV]}$', x_min=1e0, x_max=1e6, y_min=1e-1, y_max=1e4, z_min=-1, z_max=-1, cmap_str='magma', resfile='resFile', z_label='', x_scale='linear', y_scale='linear', z_scale='linear'): #, fluxfolding=''
    
    print_size=10
    label_size=30
    
    
    plt.close('all')
    font_props = {"size":label_size}
    rc("font", **font_props)
    mpl.rcParams['axes.linewidth'] = 2
    mpl.rcParams['mathtext.fontset']='stixsans'
    fig = plt.figure(figsize=(print_size*1.3, print_size))
    
    upperPlt = plt.subplot2grid((1, 1), (0, 0), rowspan=1)
    
    plt.subplots_adjust(left=0.15, right=0.75, top=0.9, bottom=0.15)
    
    
    data = np.genfromtxt( file )
    x = np.ones(len(data[:,0]))
    
    if x_scale=='log':
        x[:-1] = data[1:,0 ]*np.sqrt(data[1,0]/data[2,0])
        x[ -1] = data[-1,0 ]/np.sqrt(data[1,0]/data[2,0])
    else:
        x[:-1] = data[1:,0 ]-(data[2,0]-data[1,0])/2.
        x[ -1] = data[-1,0 ]+(data[2,0]-data[1,0])/2.

    y = np.ones(len(data[0,:]))
    if y_scale=='log':
        y[:-1] = data[0, 1:]*np.sqrt(data[0,1]/data[0,2])
        y[ -1] = data[0 ,-1]/np.sqrt(data[0,1]/data[0,2])
    else:
        y[:-1] = data[0, 1:]-(data[0,2]-data[0,1])/2.
        y[ -1] = data[0 ,-1]+(data[0,2]-data[0,1])/2.

    z = np.transpose( data[1:,1:] )

    if z_scale=='log':
        if z_min==-1 and z_max==-1:
            norm = colors.LogNorm()
        else:
            norm = colors.LogNorm(vmin=z_min, vmax=z_max)
    if z_scale=='linear':
        if z_min==-1 and z_max==-1:
            norm = colors.Normalize()
        else:
            norm = colors.Normalize(vmin=z_min, vmax=z_max)

    cmap    = plt.get_cmap(cmap_str)
    myPlt0 = upperPlt.pcolormesh(  x, y, z, norm=norm, cmap=cmap  )
    
    
    upperPlt.set_xscale(x_scale)
    upperPlt.set_yscale(y_scale)
    upperPlt.set_xlim(x_min, x_max)
    upperPlt.set_ylim(y_min, y_max)
    
    
    cbar_ax = fig.add_axes([0.8, 0.15, 0.02, 0.75])
    cbar = fig.colorbar(myPlt0, cax=cbar_ax, cmap=cmap)
    cbar.ax.set_ylabel(z_label)
    
    
    upperPlt.tick_params('both', length=20, width=2, which='major')
    upperPlt.tick_params('both', length=10, width=1, which='minor')
    
    #upperPlt.grid(b=True, which='major', alpha=0.5, linestyle='-', linewidth=2)
    
    upperPlt.tick_params(axis='both', pad=10)
    
    upperPlt.set_xlabel(x_label)
    upperPlt.set_ylabel(y_label)
    
    #upperPlt.legend(loc='upper left', fontsize=20)
    
    
    plt.savefig(resfile)




############################################################################################
#
#       Function used for plotting MultiNest results
#
############################################################################################


def paramTitle(name):
    if 'example'                          in name:    return r'Example!'
    return name


def paramScale(name):
    if 'example'                          in name:    return 'log'
    return 'linear'


def scale(sc, val):
    if sc=='log':
        return np.log10(val)
    else:
        return val


def readNamesAndRanges(file):
    all_names  = []
    all_ranges = []

    lines = [line.rstrip('\n') for line in open(file)]
    for l in lines:
        ll = l.split(' ')
        while '' in ll:
            ll.remove('')
        all_names .append(                ll[0]                     )
        if ll[1]=='N' or ll[2]=='N':
            all_ranges.append(  [0,0]                               )
            continue
        all_ranges.append(  [float       (ll[1]), float(ll[2])]     )

    return all_names, all_ranges













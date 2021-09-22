import  math
import  os
import  numpy                   as      np
from    operator                import  itemgetter

import  matplotlib              as      mpl
import  matplotlib.pyplot       as      plt
import  matplotlib.patches      as      mpatches
from    matplotlib.colors       import  BoundaryNorm
from    matplotlib.colors       import  colorConverter
from    matplotlib.ticker       import  MaxNLocator
from    matplotlib.path         import  Path
from    matplotlib              import  rc


from    PlotFunctions           import  scale
from    PlotFunctions           import  paramTitle
from    PlotFunctions           import  paramScale
from    PlotFunctions           import  readNamesAndRanges


leg_handles =   []
leg_labels  =   []

fig         =   -1
plotArray   =   -1
global_data =   -1
global_ChiSq=   -1
global_Scatter = -1

def get_legend_handels():
    global leg_handles, leg_labels
    m_leg_handles = leg_handles
    m_leg_labels  = leg_labels
    leg_handles =   []
    leg_labels  =   []
    return m_leg_handles, m_leg_labels


def prepare_triangle(npar, label_size=10, print_size=15):

    global leg_handles, leg_labels
    global fig, plotArray
    global global_data
    print 'prepare'
    #
    # Read multinest files and parameter information
    #
    plt.close('all')
    leg_handles=[]
    leg_labels=[]
    #
    # General settings: fontsize and plotgrid
    #
    font_props = {"size":label_size}
    rc("font", **font_props)
    mpl.rcParams['axes.linewidth'] = 1
    mpl.rcParams['mathtext.fontset']='stixsans'

    dim = npar - 1
    fig, plot = plt.subplots(dim, dim)
    plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)
    fig.set_size_inches(print_size, print_size)
    
    
    if dim==1:
        plotArray=[]
        plotArray.append([plot])
    else:
        plotArray=plot

    return fig, plotArray


chiSq_scatter_properties = { 'cmap':'magma', 's':3, 'chiSqRange':20}

def set_chiSq_scatter_properties(key, value):
    global chiSq_scatter_properties
    chiSq_scatter_properties[key] = value


error_scatter_properties = {  's':3, 'colors':['#641E16','#2980B9','#A2D9CE']}#A2D9CE

def set_error_scatter_properties(key, value):
    global error_scatter_properties
    error_scatter_properties[key] = value



def draw_triangle(  multinest_file_prefix, param, column_chi2=[1], column_fact='', type='error_contour', color='black', label='', shift_x=0, shift_y=0, max_n_locator=3, dir='.' ):
    print 'draw'
    global fig, plotArray, leg_labels, leg_handles
    global global_data, global_ChiSq
    #cmap_styles = 'blues'


    dim = len(param) - 1
    all_scales = []
    all_labels = []


    all_names, all_ranges = readNamesAndRanges(dir+'/MultiNest/'+multinest_file_prefix+'.ranges_plot')
    
    for n in all_names:
        all_labels.append(  paramTitle   (n)   )
        all_scales.append(  paramScale   (n)   )
    

    names  = []
    scales = []
    labels = []
    ranges = []

    for i in param:
        names .append( all_names [i] )
        labels.append( all_labels[i] )
        scales.append( all_scales[i] )
        ranges.append( all_ranges[i] )

    nLines = sum(1 for line in open(dir+'/MultiNest/'+multinest_file_prefix+'.txt'))
    global_data   = np.genfromtxt( dir+'/MultiNest/'+multinest_file_prefix+'.txt', skip_header=1 )

    global_ChiSq = np.zeros(len(global_data[:,0]))


    for i in range(len(column_chi2)):
        col = column_chi2[i]
        fac = 1.0
        if column_fact!='':
            fac = column_fact[i]
        d = global_data[:,col]*fac
        global_ChiSq = global_ChiSq + d

    if 'scatter' in type:
        ind = np.argsort( global_ChiSq )[::-1]
#        print global_ChiSq
#        print ind

        global_ChiSq = global_ChiSq[ind]
        for i in range(len(global_data[0,:])):
            d = global_data[:,i]
            d = d[ind]



    if label!='':
        leg_handles.append( mpatches.Patch(color=color) )
        leg_labels. append( label )

    # * * * * * * * * * * * * * * * * * * * * * #
    # Loop over upper triangel:  nothing        #
    # * * * * * * * * * * * * * * * * * * * * * #
    for i in range( 0, len( plotArray )  ):
        for j in range( 0, len( plotArray[i] )  ):
            plotArray[i][j].axis('off')


    # * * * * * * * * * * * * * * * * * * * * * #
    # Loop both triangles:                      #
    #    - plot in lower triangle               #
    #    - remove plots in upper triangle       #
    # * * * * * * * * * * * * * * * * * * * * * #
    for i in range(0, len(param)):
    # * * * * * * * * * * * * * * * * * * * * * #
    # Loop upper triangel:                      #
    #    - frequentist: error contours          #
    #                   (Delta X2 = 2.3 )       #
    #                   (Delta X2 = 6.2 )       #
    # * * * * * * * * * * * * * * * * * * * * * #
        #print i
        for j in range(0, i):
            
            #print str(i) + "   " + str(j)
            iP = i -1 + shift_y
            jP = j    + shift_x
            
            plotArray[iP][jP].axis('on')
           
            if type=='error_contour':
                draw_error_contour( plotArray[iP][jP], param[j], param[i], ranges[j], ranges[i], chiSq='global')
            
            if type=='chiSq_scatter':
                chiSq_min = global_ChiSq[-1]
                
                global chiSq_scatter_properties
                cmap        =   chiSq_scatter_properties    [  'cmap'  ]
                s           =   chiSq_scatter_properties    [  's'  ]
                chiSqRange  =   chiSq_scatter_properties    [  'chiSqRange'  ]

                cmap = mpl.cm.get_cmap(cmap)
                
                im = plotArray[iP][jP].scatter(  scale(scales[j], global_data[:,2+param[j]]) , scale(scales[i], global_data[:,2+param[i]]), c=global_ChiSq[:], s=s, lw=0, vmin=chiSq_min, vmax=chiSq_min+chiSqRange, cmap=cmap)
            if type=='error_scatter':
                chiSq_min = global_ChiSq[-1]
                
                global error_scatter_properties
                colors      =   error_scatter_properties    [  'colors'  ]
                s           =   error_scatter_properties    [  's'  ]
                
                n=[]
                last_chiSq = global_ChiSq[0]
                not3, not2, not1  = True, True, True
                
                for k in range(0, len(global_ChiSq)):
                    chiSq = global_ChiSq[k]
                    if chiSq-chiSq_min < 11.83 and not3:
                        not3 = False
                        n.append(k)
                    if chiSq-chiSq_min < 6.18 and not2:
                        not2 = False
                        n.append(k)
                    if chiSq-chiSq_min < 2.3 and not1:
                        not1 = False
                        n.append(k)
                        break
                n.append(len(global_ChiSq))
                
                im = plotArray[iP][jP].scatter(  scale(scales[j], global_data[n[0]:n[1],2+param[j]]) , scale(scales[i], global_data[n[0]:n[1],2+param[i]]), c=colors[2], s=s, lw=0 )
                im = plotArray[iP][jP].scatter(  scale(scales[j], global_data[n[1]:n[2],2+param[j]]) , scale(scales[i], global_data[n[1]:n[2],2+param[i]]), c=colors[1], s=s, lw=0 )
                im = plotArray[iP][jP].scatter(  scale(scales[j], global_data[n[2]:n[3],2+param[j]]) , scale(scales[i], global_data[n[2]:n[3],2+param[i]]), c=colors[0], s=s, lw=0 )


            #
            # General axes settings needed for both error plots
            #
            plotArray[iP][jP].set_xlim( [ scale( scales[j], ranges[j][0]), scale( scales[j], ranges[j][1]) ]  )
            plotArray[iP][jP].set_ylim( [ scale( scales[i], ranges[i][0]), scale( scales[i], ranges[i][1]) ]  )

            plotArray[iP][jP].xaxis.get_major_formatter().set_powerlimits((-3, 4))
            plotArray[iP][jP].yaxis.get_major_formatter().set_powerlimits((-3, 4))

            plotArray[iP][jP].xaxis.set_major_locator(MaxNLocator(max_n_locator))
            plotArray[iP][jP].yaxis.set_major_locator(MaxNLocator(max_n_locator))
            for tick in plotArray[iP][jP].get_xticklabels():
                tick.set_rotation(45)
                tick.set_horizontalalignment('right')
            plotArray[iP][jP].xaxis.set_tick_params(length=6, width=1)
            plotArray[iP][jP].yaxis.set_tick_params(length=6, width=1)
            if j==0:
                plotArray[iP][jP].set_ylabel( labels[i] )
                if scales[i]=='log':
                    plotArray[iP][jP].set_ylabel( 'log(' + labels[i] + ')' )
            if i==len(param)-1:
                plotArray[iP][jP].set_xlabel( labels[j] )
                if scales[j]=='log':
                    plotArray[iP][jP].set_xlabel( 'log(' + labels[j] + ')' )
            if i<len(param)-1:
                plotArray[iP][jP].xaxis.set_tick_params( labelsize=0  )
                plotArray[iP][jP].xaxis.get_offset_text().set_fontsize(0)
            if j>0:
                plotArray[iP][jP].yaxis.set_tick_params( labelsize=0  )
                plotArray[iP][jP].yaxis.get_offset_text().set_fontsize(0)

        # * * * * * * * * * * * * * * * * * * * * * #
        # Loop over upper triangel:  nothing        #
        # * * * * * * * * * * * * * * * * * * * * * #
        for j in range(i+1, dim):
            plotArray[i][j].axis('off')


def draw_diagonal(  multinest_file_prefix, param, column_chi2=[1], column_fact='', color='black', label='', max_n_locator=3, dir='.'):
    print 'draw'
    global fig, plotArray, leg_labels, leg_handles
    global global_data, global_ChiSq
    #cmap_styles = 'blues'


    dim = len(param) - 1

    all_scales = []
    all_labels = []
    all_names, all_ranges = readNamesAndRanges(dir+'/MultiNest/'+multinest_file_prefix+'.ranges_plot')
    
    for n in all_names:
        all_labels.append(  paramTitle   (n)   )
        all_scales.append(  paramScale   (n)   )
    

    names  = []
    scales = []
    labels = []
    ranges = []

    for i in param:
        names .append( all_names [i] )
        labels.append( all_labels[i] )
        scales.append( all_scales[i] )
        ranges.append( all_ranges[i] )

    nLines = sum(1 for line in open(dir+'/MultiNest/'+multinest_file_prefix+'.txt'))
    global_data   = np.genfromtxt( dir+'/MultiNest/'+multinest_file_prefix+'.txt', skip_header=1 )

    global_ChiSq = np.ones(len(global_data[:,0]))

    for i in range(len(column_chi2)):
        col = column_chi2[i]
        fac = 1.0
        if column_fact!='':
            fac = column_fact[i]
        d = global_data[:,col]*fac
        global_ChiSq = global_ChiSq + d

    if label!='':
        leg_handles.append( mpatches.Patch(color=color) )
        leg_labels. append( label )

    # * * * * * * * * * * * * * * * * * * * * * #
    # Loop over diagonal:                       #
    # * * * * * * * * * * * * * * * * * * * * * #
    for i in range(0, len(param) ):
        #
        # General input needed for both error types
        #
        # Get parameter range
        x_min = scale( scales[i], ranges[i][0] )
        x_max = scale( scales[i], ranges[i][1] )
        # Calculate bin left edges
        
        plotArray[i][i].axis('on')
        draw_chi_square( plotArray[i][i], param[i], ranges[i], chiSq='global')

        plotArray[i][i].set_ylim( [0, 10]  )
        if i == 0:
            plotArray[i][i].set_ylabel(r'$\Delta \chi^2$')

        plotArray[i][i].xaxis.set_tick_params(length=6, width=1)
        plotArray[i][i].yaxis.set_tick_params(length=6, width=1)
        plotArray[i][i].xaxis.get_major_formatter().set_powerlimits((-3, 4))
        plotArray[i][i].set_xlim( [ scale( scales[i], ranges[i][0]), scale( scales[i], ranges[i][1]) ]  )
        plotArray[i][i].xaxis.set_major_locator(MaxNLocator(max_n_locator))
        plotArray[i][i].yaxis.set_major_locator(MaxNLocator(max_n_locator))
        for tick in plotArray[i][i].get_xticklabels():
            tick.set_rotation(45)
            tick.set_horizontalalignment('right')
        

        if i == len(param)-1:
            plotArray[i][i].set_xlabel( labels[i] )
            if scales[i]=='log':
                plotArray[i][i].set_xlabel(r''+'log('+labels[i]+')'+'')
        if i < len(param)-1:
            plotArray[i][i].xaxis.set_tick_params( labelsize = 0  )
            plotArray[i][i].xaxis.get_offset_text().set_fontsize(0)
        if i >0:
            plotArray[i][i].yaxis.set_tick_params( labelsize = 0  )
            plotArray[i][i].yaxis.get_offset_text().set_fontsize(0)




error_counter_properties = { 'scale1':'linear', 'scale2':'linear', 'grid_size':40, 'color':'black', 'fc':'', 'cmap':'', 'alpha':1.0, 'label':'', 'fill':True, 'smoothing':10, 'column_fact':'', 'zorder_fill':1, 'zorder_line':10, 'draw_power10':False}
def reset_error_counter_properties():
    global error_counter_properties
    error_counter_properties = { 'scale1':'linear', 'scale2':'linear', 'grid_size':40, 'color':'black', 'fc':'', 'cmap':'', 'alpha':1.0, 'label':'', 'fill':True, 'smoothing':10, 'column_fact':'', 'zorder_fill':1, 'zorder_line':10, 'draw_power10':False}
def set_error_counter_properties(key, value):
    global error_counter_properties
    error_counter_properties[key] = value

def draw_error_contour( ax, param1, param2, ranges1, ranges2, data='global', chiSq=1, label='' ):
    global global_data, global_ChiSq
    if data=='global':
        data    = global_data
    if chiSq=='global':
        chiSq   = global_ChiSq
    else:
        chiSq   = data[:,chiSq]

    global error_counter_properties
    scale1      =   error_counter_properties[  'scale1'       ]
    scale2      =   error_counter_properties[  'scale2'       ]
    grid_size   =   error_counter_properties[  'grid_size'    ]
    color       =   error_counter_properties[  'color'        ]
    fc          =   error_counter_properties[  'fc'           ]
    cmap        =   error_counter_properties[  'cmap'         ]
    alpha       =   error_counter_properties[  'alpha'        ]
    label       =   error_counter_properties[  'label'        ]
    fill        =   error_counter_properties[  'fill'         ]
    smoothing   =   error_counter_properties[  'smoothing'    ]
    zorder_fill =   error_counter_properties[  'zorder_fill'  ]
    zorder_line =   error_counter_properties[  'zorder_line'  ]
    draw_power10=   error_counter_properties[  'draw_power10' ]


    #'linear', scale2='linear', grid_size = 40, color='black', fc='', cmap='', label='', fill=True, smoothing=10, column_fact='', zorder_fill=1, zorder_line=100

    chiSq_min = np.amin(chiSq)

    x_min = scale( scale1, ranges1[0] )
    x_max = scale( scale1, ranges1[1] )
    y_min = scale( scale2, ranges2[0] )
    y_max = scale( scale2, ranges2[1] )
    # Get gird of 2D parameter space (bin edges)
    dx = 1.*(  x_max - x_min  )/( grid_size )
    dy = 1.*(  y_max - y_min  )/( grid_size )
    x, y = np.mgrid[  slice(x_min, x_max + dx/2, dx), slice(y_min, y_max + dy/2, dy)  ]

    # Get contour values (delta chi square = 1 and 2)
    v  = [0, chiSq_min+2.3, chiSq_min+6.2, chiSq_min+11.8]
    # initialize the bin content with -1
    c = np.ones( (grid_size+4,grid_size+4) )
    c = -1*c
    c_max = 0


    # marginalize chi suqare in each bin
    for k in range(0, len(data) ):
        chi = chiSq[k]
        x_ = scale( scale1, data[k,2+param1] )
        y_ = scale( scale2, data[k,2+param2] )
        n1 = int( (x_-x_min)/dx )
        n2 = int( (y_-y_min)/dy )
        if n1<0 or n2<0 or n1>=grid_size or n2>=grid_size:
            continue
        c_ = c[n1+2][n2+2]
        if chi<c_ or c_<0:
            c[n1+2][n2+2] = chi
        if chi>c_max:
            c_max = chi
    for i1 in range(0, grid_size+4):
        for i2 in range(0, grid_size+4):
            if c[i1][i2]<0:
                c[i1][i2] = chiSq_min+30
    for it in range(0,smoothing):
        for i1 in range(0, grid_size):
            for i2 in range(0, grid_size):
                j1 = i1+2
                j2 = i2+2
                vec = []
                vec.append(c[j1,j2])
                vec.append( 0.5*c[j1-1,j2  ]+0.5*c[j1+1,j2  ] )
                vec.append( 0.5*c[j1  ,j2-1]+0.5*c[j1  ,j2+1] )

                c[j1,j2] = min(vec)

    if label!='':
        myfc = np.asarray(colorConverter.to_rgba(fc))
        myfc[3] = alpha
        leg_handles.append( mpatches.Patch(color=color, fc=myfc, lw=2) )
        leg_labels. append( label )

    # contours are *point* based plots, so convert our bound into point centers
    # plot contours
    CS = ax.contour(x[:-1, :-1] + dx/2., y[:-1, :-1] + dy/2., c[2:-2, 2:-2], v, colors=color)
    plt.setp(CS.collections , linewidth=1, zorder=zorder_line)
    if fill:
        colors=[color,color,color]
        if fc!='':
            colors=[fc, fc, fc]
        if cmap!='':
            cmap = mpl.cm.get_cmap(cmap)
            colors=[cmap(0.9), cmap(0.6), cmap(0.3)]

        CSF = ax.contourf(x[:-1, :-1] + dx/2., y[:-1, :-1] + dy/2., c[2:-2, 2:-2], v, lw=2, colors=colors, alpha=alpha, zorder=zorder_fill,  )
        plt.setp(CSF.collections , zorder=zorder_fill)
        if draw_power10:
            k = 3
            for j in range(len(CSF.collections)):
                i = len(CSF.collections)-1-j
                for p in CSF.collections[i].get_paths():
                    k = k-1
                    limit = p.vertices
                    limit = np.power(10, limit)
                    limit[:,0] *= 1e-3
                    path = Path(limit, p.codes)
                    fill_patch = mpatches.PathPatch(path, fc=colors[k%3], lw=0, zorder=zorder_fill, alpha=alpha)
                    ax.add_patch(fill_patch)
    if draw_power10:
        k = 0
        for j in range(len(CS.collections)):
            i = len(CS.collections)-1-j
            for p in CS.collections[i].get_paths():
                k = k+1
                limit = p.vertices
                limit = np.power(10, limit)
                limit[:,0] *= 1e-3
                path = Path(limit, p.codes)
                e_patch = mpatches.PathPatch(path, lw=2, ec=color, zorder=zorder_line, fill=False)
                ax.add_patch(e_patch)


    return CS


chi_square_contour_properties = { 'scale1':'linear', 'scale2':'linear', 'grid_size':100, 'cmap':'magma', 'alpha':1.0, 'label':'', 'column_fact':'', 'zorder':1, 'range_delta_chi2':50}
def reset_chi_square_contour_properties():
    global chi_square_contour_properties
    chi_square_contour_properties = { 'scale1':'linear', 'scale2':'linear', 'grid_size':100, 'cmap':'magma', 'alpha':1.0, 'label':'', 'column_fact':'', 'zorder':1, 'range_delta_chi2':50}
def set_chi_square_contour_properties(key, value):
    global chi_square_contour_properties
    chi_square_contour_properties[key] = value

def draw_chi_square_contour( ax, param1, param2, ranges1, ranges2, data='global', chiSq=1, label='' ):
    global global_data, global_ChiSq
    if data=='global':
        data    = global_data
    if chiSq=='global':
        chiSq   = global_ChiSq
    else:
        chiSq   = data[:,chiSq]

    global chi_square_contour_properties
    scale1      =   chi_square_contour_properties[  'scale1'       ]
    scale2      =   chi_square_contour_properties[  'scale2'       ]
    grid_size   =   chi_square_contour_properties[  'grid_size'    ]
    cmap        =   chi_square_contour_properties[  'cmap'         ]
    alpha       =   chi_square_contour_properties[  'alpha'        ]
    label       =   chi_square_contour_properties[  'label'        ]
    zorder      =   chi_square_contour_properties[  'zorder'       ]
    range_delta_chi2      =   chi_square_contour_properties[  'range_delta_chi2'       ]


    chiSq_min = np.amin(chiSq)
    
    x_min = scale( scale1, ranges1[0] )
    x_max = scale( scale1, ranges1[1] )
    y_min = scale( scale2, ranges2[0] )
    y_max = scale( scale2, ranges2[1] )
    # Get gird of 2D parameter space (bin edges)
    dx = 1.*(  x_max - x_min  )/( grid_size )
    dy = 1.*(  y_max - y_min  )/( grid_size )
    x, y = np.mgrid[  slice(x_min, x_max + dx/2, dx), slice(y_min, y_max + dy/2, dy)  ]

    # Get contour values
    v = np.arange(chiSq_min, chiSq_min+range_delta_chi2, 1.*range_delta_chi2/200)

    # initialize the bin content with -1
    c = np.ones( (grid_size,grid_size) )
    c = -1*c
    c_max = 0
    # marginalize chi suqare in each bin
    for k in range(0, len(chiSq) ):
        chi = chiSq[k]
        x_ = scale( scale1, data[k,2+param1] )
        y_ = scale( scale2, data[k,2+param2] )
        n1 = int( (x_-x_min)/dx )
        n2 = int( (y_-y_min)/dy )
        if n1<0 or n2<0 or n1>=grid_size or n2>=grid_size:
            continue
        c_ = c[n1][n2]
        if chi<c_ or c_<0:
            c[n1][n2] = chi
        if chi>c_max:
            c_max = chi
    for i1 in range(0, grid_size):
        for i2 in range(0, grid_size):
            if c[i1][i2]<0:
                c[i1][i2] = c_max
    if 2==2:
        toWrite = str('# log10(m_DM/MeV)').ljust(20,' ')+str('log10(sv/(cm3/s))').ljust(20,' ')+str('chiSquare').ljust(20,' ')
        for i1 in range(0, grid_size):
            for i2 in range(0, grid_size):
                toWrite += '\n'+str(x[i1,i2]).ljust(20,' ')+str(y[i1,i2]).ljust(20,' ')+str(c[i1,i2]).ljust(20,' ')
        target = open('likelihood.txt', 'w')
        target.write(toWrite)
        target.close()
    # contours are *point* based plots, so convert our bound into point centers
    # plot contours
    ax.contourf(x[:-1, :-1] + dx/2., y[:-1, :-1] + dy/2., c, v, extend='max', cmap=cmap)




chi_square_properties = { 'scale':'linear', 'grid_size':40, 'color':'black', 'alpha':1.0, 'label':'', 'smoothing':10, 'zorder':1 }
def reset_chi_square_properties():
    global chi_square_properties
    chi_square_properties = { 'scale':'linear', 'grid_size':40, 'color':'black', 'alpha':1.0, 'label':'', 'smoothing':10, 'zorder':1 }
def set_chi_square_properties(key, value):
    global chi_square_properties
    chi_square_properties[key] = value

def draw_chi_square( ax, param, ranges, data='global', chiSq=1, iteration=0, max_n_iterations=5     ):
    dataO   = data
    chiSqO  = chiSq
    global global_data, global_ChiSq
    if data=='global':
        data    = global_data
    if chiSq=='global':
        chiSq   = global_ChiSq
    else:
        chiSq   = data[:,chiSq]

    global chi_square_properties
    scaleV      =   chi_square_properties[  'scale'        ]
    grid_size   =   chi_square_properties[  'grid_size'    ]
    color       =   chi_square_properties[  'color'        ]
    alpha       =   chi_square_properties[  'alpha'        ]
    label       =   chi_square_properties[  'label'        ]
    smoothing   =   chi_square_properties[  'smoothing'    ]
    zorder      =   chi_square_properties[  'zorder'       ]
    

    chiSq_min = np.amin(chiSq)

    # Get parameter range

    index_best  = np.argmin(chiSq)

    chiSq_min   = chiSq[index_best]
    x_best      = data[index_best,param+2]
    if ranges=='sm':
        ranges = [x_best-200, x_best+200]
    if ranges=='norm':
        ranges = [x_best-0.15, x_best+0.15]
    x_min = scale( scaleV, ranges[0] )
    x_max = scale( scaleV, ranges[1] )
    
    # Calculate bin left edges
    dx = (  x_max - x_min  )/( grid_size )
    if dx==0:
        print "Warning. Parameter: '" + str(label) + "' has range zero. Can't produce plot."
        return
    bins = np.arange( x_min, x_max-dx/2, dx )
    # Initialize bin content with -1
    cont = np.ones ( grid_size ) # 0: total chiSquare, i: chiSquare for model i
    cont = cont*(-1)
    chiMax = 0
    chiMin = -1
    pBest  = -1
    nBest  = -1
    # Marginalize chi square in eacht bin
    for k in range(0, len(chiSq) ):
        chi = chiSq[k]
        x_  = scale(scaleV, data[k, param+2])
        n_  = (x_-x_min)/dx
        if n_!=n_:
            print "Warning. Parameter: '" + str(label) + "' has range zero. Can't produce plot."
            return
        n   = int( n_ )
        if n < 0 or n>=grid_size:
            continue
        c   = cont[n]
        if chi>chiMax:
            chiMax = chi
        if chi<chiMin or chiMin<0:
            chiMin = chi
            pBest  = x_
            nBest  = n
        if chi<c or c<0:
            cont[n] = chi

    for i in range(0, len(cont)):
        if cont[i]<0:
            cont[i] = chiMax

    for r in range(0,smoothing):
        for i in range(0,grid_size):
            s_chiL = 30 + chiMin
            s_chi  = cont[i]
            s_chiU = 30 + chiMin
            if i>0:
                s_chiL = cont[i-1]
            if i<len(cont)-1:
                s_chiU = cont[i+1]
            s_chi_int = s_chiL*0.5+s_chiU*0.5
            if ( s_chi_int < s_chi):
                cont[i] = s_chi_int


    x=[]
    y=[]

    for i in range(0, len(cont)):
        x.append(x_min+i*dx+dx/2)
        #x.append(x_min+i*dx+dx)
        y.append(cont[i]-chiMin)
        #y.append(cont[i]-chiMin)



    sigma1_l = 0
    sigma1_u = 0
    for i in range(nBest, grid_size):
        if (cont[i]-chiMin)>9:
            break
        sigma1_u = x_min+dx/2+i*dx -pBest

    for i in range(0, nBest):
        j = nBest - i-1
        if (cont[j]-chiMin)>9:
            break
        sigma1_l = pBest -( x_min+dx/2+j*dx )

    bin_width   = (ranges[1]-ranges[0])/grid_size

    sigma1_l = max(sigma1_l, bin_width)
    sigma1_u = max(sigma1_u, bin_width)

    sigma_width = (sigma1_u + sigma1_l)

    if iteration<max_n_iterations:
        if (bin_width*10>sigma_width):
            #print "Go to "+str(iteration)+". iteration"
            x,y = draw_chi_square( ax, param, [x_best-2*sigma1_l, x_best+2*sigma1_u], dataO, chiSqO, iteration+1, max_n_iterations )


    if iteration==0:
        ax.plot(x, y, color=color, lw=2, dashes=[])

    return x,y


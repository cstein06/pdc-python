
from numpy import *
import matplotlib.pyplot as pp

from globals import *
from matplotlib.pyplot import imshow
import matplotlib.colors as mc_

#def plot_all(mes, th, ic1, ic2, ss = None, sample_f = 1.0, 
#             logss = False, sqrtmes = False, plotf = None):

def plot_all(**args):
    '''Plots nxn graphics, with confidence intervals and threshold. 
       If ss == True, plots ss in the diagonal.
       Already expects data in power form: abs(x)^2'''
    
    read_args(args)
    
    print '\nPlotting...'
    
    mes = res_.mes.copy()
    
    th = None
    if res_.th is not None:
        th = res_.th.copy()
        
    ic1 = None
    ic2 = None
    if res_.ic1 is not None:
        ic1 = res_.ic1.copy()
        ic2 = res_.ic2.copy()
    #else:
    #    ic1 = zeros(mes.shape)
    #    ic2 = zeros(mes.shape)
    
    ss = None
    if res_.ss is not None:
        ss = res_.ss.copy()
    
    if pr_.logss and ss is not None:
        ss = log(ss)
       
    if pr_.sqrtmes: 
        mes = sqrt(mes)
        th = sqrt(th)
        ic1 = sqrt(ic1)
        ic2 = sqrt(ic2)
        
        
    n,n,nf = mes.shape
    
    #print 'mes', mes.min(), mes.max()
    
    #pp.ion()
    
    if not pr_.power:
        print 'Taking squared power for plotting!'
        if mes.dtype != 'complex':
            print 'But it was already non-complex data...'
        mes = (mes*mes.conj()).real
        
    if mes.dtype == 'complex':
        print 'Plotting complex data, something seems to be wrong.'

    if pr_.plot_title is not None:
        pp.suptitle(pr_.plot_title)
    else:
        pp.suptitle(mnames_[pr_.alg])

    x = pr_.sample_f*arange(nf)/(2.0*nf)
    for i in range(n):
        for j in range(n):
            pp.subplot(n,n,i*n+j+1)
            
            if (i < n-1):
                pp.xticks([])
            if (j > 0):
                pp.yticks([])
                
            if (i == n-1):
                try:
                    if pr_.plot_labels != None:
                        pp.xlabel(pr_.plot_labels[j])
                    else:
                        pp.xlabel(str(j+1))
                except:
                    print '\nProblem with plot labels.'
                    pp.xlabel(str(j+1))
                    
            if (j == 0):
                try:
                    if pr_.plot_labels != None:
                        pp.ylabel(pr_.plot_labels[i])
                    else:
                        pp.ylabel(str(i+1))
                except:
                    print '\nProblem with plot labels.'
                    pp.ylabel(str(i+1))
                
            if i == j and not pr_.plot_diag:
                #pp.xticks([])
                if pr_.plotf != None:
                    pp.xlim([0, pr_.plotf])
                else:
                    pp.xlim([0, pr_.sample_f/2.0])
                continue

            pp.subplot(n,n,i*n+j+1)
            #over = mes[i,j][mes[i,j]>th[i,j]]
            #overx = x[mes[i,j]>th[i,j]]
            #over = mes[i,j]
            #overx = x
            
            #Old code
            #under = mes[i,j][mes[i,j]<=th[i,j]]
            #underx = x[mes[i,j]<=th[i,j]]
            #pp.plot(x, th[i,j], 'r:', x, ic1[i,j], 'k:', x, ic2[i,j], 'k:', 
            #        overx, over, 'b-', underx, under, 'r-')
            
            if pr_.plot_th and th is not None:
                pp.plot(x, th[i,j], 'r:')
            if pr_.plot_ic and ic1 is not None:
                pp.plot(x, ic1[i,j], 'k:')
                pp.plot(x, ic2[i,j], 'k:')
            
            if pr_.plot_color != None:
                pp.plot(x, mes[i,j,:], color=pr_.plot_color)
            else:
                pp.plot(x, mes[i,j,:])
            #pp.plot(x, th[i,j], 'r:', 
            #        overx, over, 'b-')
            
            if pr_.plot_th and th is not None:
                #Complicated code for underthreshold painting
                k = 0
                while(k < nf):
                    while(mes[i,j,k] >= th[i,j,k]):
                        k = k+1
                        if (k == nf): break
                    if (k == nf): break
                    kold = k
                    while(mes[i,j,k] < th[i,j,k]):
                        k = k+1
                        if (k == nf): break
                    pp.plot(x[kold:k], mes[i,j,kold:k], 'r-')
            
            pp.ylim(-0.05,1.05)
            
            
            if pr_.plotf is not None:
                pp.xlim([0, pr_.plotf])
                
        if (pr_.ss and ss is not None):
            ax = pp.subplot(n,n,i*n+i+1).twinx()
            ax.plot(pr_.sample_f*arange(nf)/(2.0*nf), ss[i,i,:], color='g')
            if pr_.logss:
                ax.set_ylim(ymin = ss[i,i,:].min(), ymax = ss[i,i,:].max())
            else:
                ax.set_ylim(ymin = 0, ymax = ss[i,i,:].max())
                
            if (i < n-1):
                ax.set_xticks([])
            
            if pr_.plotf != None:
                ax.set_xlim(xmin = 0, xmax = pr_.plotf)
            else:
                ax.set_xlim(xmin = 0, xmax = pr_.sample_f/2.0)
        
        pp.draw()
    #pp.show()
    
def plot_coherogram(res, states = None):
    
    nwin,n,n,nf = res.shape
    
    if pr_.plot_title is not None:
        pp.suptitle(pr_.plot_title)
    else:
        pp.suptitle(mnames_[pr_.alg])

    #pp.ion()
    #print res.shape
    #states = states[:400]
    #res = res[:400,:,:,:120]
    
    if size(states) != res.shape[0]:
        print 'states doesn\'t match  res size.'
    
    for i in range(n):
        pp.subplot(n+1,n,i+1)
        if states is not None:
            #saux = zeros([5, size(states)])
            #saux[0] = states
            
            saux = states.reshape(1,-1,1)
            
            srgb = repeat(saux, 3, axis = 2)
            
            for i in arange(saux.shape[1]):
                if srgb[0,i,0] > 0:
                    col = pr_.state_colors[int(srgb[0,i,0]-1)]
                else:
                    col = mc_.colorConverter.to_rgb('k')
                srgb[0,i] = mc_.colorConverter.to_rgb(col)
            
            imshow(srgb, origin='lower', extent=(0,5,0,3))
            pp.xticks([])
            pp.yticks([])
        
        
    for i in range(n):
        for j in range(n):
            
            pp.subplot(n+1,n,n+i*n+j+1)
            
            if (i == n-1):
                try:
                    if pr_.plot_labels != None:
                        pp.xlabel(pr_.plot_labels[j])
                    else:
                        pp.xlabel(str(j+1))
                except:
                    print '\nProblem with plot labels.'
                    pp.xlabel(str(j+1))
                    
            if (j == 0):
                try:
                    if pr_.plot_labels != None:
                        pp.ylabel(pr_.plot_labels[i])
                    else:
                        pp.ylabel(str(i+1))
                except:
                    print '\nProblem with plot labels.'
                    pp.ylabel(str(i+1))
                    
            if (i < n-1):
                pp.xticks([])
            if (j > 0):
                pp.yticks([])
            
#            if i == j:
#                if states is not None:
#                    imshow(states.reshape(1,-1), extent=(0,1,0,1))
#                    pp.xticks([])
#                    pp.yticks([])
#                continue
            
            
            if pr_.plotf is not None:
                auxres = res[:,:,:,:pr_.plotf]
            else:
                auxres = res
                
            if i == j and pr_.ss and pr_.logss:
                imshow(log(auxres[:,i,j,:]).T, origin='lower', extent=(0,5,0,3))
            else:
                imshow(auxres[:,i,j,:].T, origin='lower', extent=(0,5,0,3))
            
    
    pp.draw()
    
    
#    
##pdc, ss = None, sample_f = 1.0, power = True, logss = False
#def pdc_plot(mes = None, ss = None, **args):
#    '''Plots nxn graphics. 
#       mes(n,n,nf) -> data
#       
#       If ss == True, plots ss in the diagonal.
#       Expects data in complex form. Does: abs(x)^2 before plotting.'''
#    
#    read_args(args)
#    
#    print '\nPlotting...'
#    
#    if mes == None:
#        mes = res_.mes.copy()
#        
#    if res_.ss == None and pr_.ss:
#        ss = res_.ss.copy()
#       
#    n,n,nf = mes.shape
#    
#    if not pr_.power:
#        print 'Taking squared power for plotting!'
#        if mes.dtype != 'complex':
#            print 'But it was already non-complex data...'
#        mes = (mes*mes.conj()).real
#        
#    if mes.dtype == 'complex':
#        print 'Plotting complex data, something seems to be wrong.'
#    
#    if pr_.logss and pr_.ss:
#        ss = log(ss)
#        
#    
#        #print ss.shape
#        #print ss[0,0,:10]
#    
#    #pp.ion()
#    
#    for i in range(n):
#        for j in range(n):
#            pp.subplot(n,n,i*n+j+1)
#            if pr_.plot_color != None:
#                pp.plot(pr_.sample_f*arange(nf)/(2.0*nf), mes[i,j,:], color=pr_.plot_color)
#            else:
#                pp.plot(pr_.sample_f*arange(nf)/(2.0*nf), mes[i,j,:])
#            pp.ylim(-0.05,1.05)
#            if (i < n-1):
#                pp.xticks([])
#            else:
#                try:
#                    if pr_.plot_labels != None:
#                        pp.xlabel(pr_.plot_labels[j])
#                    else:
#                        pp.xlabel(str(j+1))
#                except:
#                    print '\nProblem with plot labels.'
#                    pp.xlabel(str(j+1))
#                
#            if (j > 0):
#                pp.yticks([])
#            else:
#                try:
#                    if pr_.plot_labels != None:
#                        pp.ylabel(pr_.plot_labels[i])
#                    else:
#                        pp.ylabel(str(i+1))
#                except:
#                    print '\nProblem with plot labels.'
#                    pp.ylabel(str(i+1))
#            
#            if pr_.plotf != None:
#                pp.xlim([0, pr_.plotf])
#        if (pr_.ss):
#            ax = pp.subplot(n,n,i*n+i+1).twinx()
#            ax.plot(pr_.sample_f*arange(nf)/(2.0*nf), ss[i,i,:], color='g')
#            if pr_.logss:
#                ax.set_ylim(ymin = ss[i,i,:].min(), ymax = ss[i,i,:].max())
#            else:
#                ax.set_ylim(ymin = 0, ymax = ss[i,i,:].max())
#                
#            if (i < n-1):
#                ax.set_xticks([])
#            if pr_.plotf != None:
#                ax.set_xlim(xmin = 0, xmax = pr_.plotf)
#    
#        pp.draw()
#    #pp.show()

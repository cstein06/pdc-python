
from numpy import *
import matplotlib.pyplot as pp

from globals import *

#def plot_all(mes, th, ic1, ic2, ss = None, sample_f = 1.0, 
#             logss = False, sqrtmes = False, plotf = None):

def plot_all():
    '''Plots nxn graphics, with confidence intervals and threshold. 
       If ss == True, plots ss in the diagonal.
       Already expects data in power form: abs(x)^2'''
    
    print '\nPlotting...'
    
    mes = res_.mes.copy()
    th = res_.th.copy()
    ic1 = res_.ic1.copy()
    ic2 = res_.ic2.copy()
    
    if res_.ss != None:
        ss = res_.ss.copy()
    
    if pr_.logss:
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

    x = pr_.sample_f*arange(pr_.nf)/(2.0*pr_.nf)
    for i in range(n):
        for j in range(n):
            pp.subplot(n,n,i*n+j+1)
            
            if (i < n-1):
                pp.xticks([])
            if (j > 0):
                pp.yticks([])
                
            if i == j and not pr_.plot_diag:
                pp.xticks([])
                if pr_.plotf != None:
                    pp.xlim([0, pr_.plotf])
                else:
                    pp.xlim([0, pr_.sample_f/2.0])
                continue

            pp.subplot(n,n,i*n+j+1)
            #over = mes[i,j][mes[i,j]>th[i,j]]
            #overx = x[mes[i,j]>th[i,j]]
            over = mes[i,j]
            overx = x
            
            #Old code
            #under = mes[i,j][mes[i,j]<=th[i,j]]
            #underx = x[mes[i,j]<=th[i,j]]
            #pp.plot(x, th[i,j], 'r:', x, ic1[i,j], 'k:', x, ic2[i,j], 'k:', 
            #        overx, over, 'b-', underx, under, 'r-')
            
            #pp.plot(x, th[i,j], 'r:', x, ic1[i,j], 'k:', x, ic2[i,j], 'k:', 
            #        overx, over, 'b-')
            pp.plot(x, th[i,j], 'r:', 
                    overx, over, 'b-')
            
            #Complicated code for underthreshold painting
            k = 0
            while(k < pr_.nf):
                while(mes[i,j,k] >= th[i,j,k]):
                    k = k+1
                    if (k == pr_.nf): break
                if (k == pr_.nf): break
                kold = k
                while(mes[i,j,k] < th[i,j,k]):
                    k = k+1
                    if (k == pr_.nf): break
                pp.plot(x[kold:k], mes[i,j,kold:k], 'r-')
            
            pp.ylim(-0.05,1.05)
            
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
            
            if pr_.plotf != None:
                pp.xlim([0, pr_.plotf])
                
        if (pr_.ss):
            ax = pp.subplot(n,n,i*n+i+1).twinx()
            ax.plot(pr_.sample_f*arange(pr_.nf)/(2.0*pr_.nf), ss[i,i,:], color='g')
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
    
    
#pdc, ss = None, sample_f = 1.0, power = True, logss = False
def pdc_plot(mes = None, ss = None, **args):
    '''Plots nxn graphics. 
       mes(n,n,nf) -> data
       
       If ss == True, plots ss in the diagonal.
       Expects data in complex form. Does: abs(x)^2 before plotting.'''
    
    read_args(args)
    
    print '\nPlotting...'
    
    if mes == None:
        mes = res_.mes.copy()
        
    if res_.ss == None and pr_.ss:
        ss = res_.ss.copy()
       
    n,n,nf = mes.shape
    
    if not pr_.power:
        print 'Taking squared power for plotting!'
        if mes.dtype != 'complex':
            print 'But it was already non-complex data...'
        mes = (mes*mes.conj()).real
        
    if mes.dtype == 'complex':
        print 'Plotting complex data, something seems to be wrong.'
    
    if pr_.logss and pr_.ss:
        ss = log(ss)
        
    
        #print ss.shape
        #print ss[0,0,:10]
    
    #pp.ion()
    
    for i in range(n):
        for j in range(n):
            pp.subplot(n,n,i*n+j+1)
            if pr_.plot_color != None:
                pp.plot(pr_.sample_f*arange(nf)/(2.0*nf), mes[i,j,:], color=pr_.plot_color)
            else:
                pp.plot(pr_.sample_f*arange(nf)/(2.0*nf), mes[i,j,:])
            pp.ylim(-0.05,1.05)
            if (i < n-1):
                pp.xticks([])
            else:
                try:
                    if pr_.plot_labels != None:
                        pp.xlabel(pr_.plot_labels[j])
                    else:
                        pp.xlabel(str(j+1))
                except:
                    print '\nProblem with plot labels.'
                    pp.xlabel(str(j+1))
                
            if (j > 0):
                pp.yticks([])
            else:
                try:
                    if pr_.plot_labels != None:
                        pp.ylabel(pr_.plot_labels[i])
                    else:
                        pp.ylabel(str(i+1))
                except:
                    print '\nProblem with plot labels.'
                    pp.ylabel(str(i+1))
            
            if pr_.plotf != None:
                pp.xlim([0, pr_.plotf])
        if (pr_.ss):
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
    
        pp.draw()
    #pp.show()

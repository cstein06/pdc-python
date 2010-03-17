
from numpy import *
import matplotlib.pyplot as pp

#def plot_all(mes, th, ic1, ic2, ss = None, sample_f = 1.0, 
#             logss = False, sqrtmes = False, plotf = None):

def plot_all(res_, pr_):
    '''Plots nxn graphics, with confidence intervals and threshold. 
       If ss == True, plots ss in the diagonal.
       Already expects data in power form: abs(x)^2'''
    
    r = res_.copy()
    
    if pr_.logss:
        r.ss = log(r.ss)
       
    if pr_.sqrtmes: 
        r.mes = sqrt(r.mes)
        r.th = sqrt(r.th)
        r.ic1 = sqrt(r.ic1)
        r.ic2 = sqrt(r.ic2)
        
    n,n,nf = r.mes.shape
    
    print 'mes', r.mes.min(), r.mes.max()
    
    #pp.ion()

    x = pr_.sample_f*arange(pr_.nf)/(2.0*pr_.nf)
    for i in range(n):
        for j in range(n):
            pp.subplot(n,n,i*n+j+1)
            #over = mes[i,j][mes[i,j]>th[i,j]]
            #overx = x[mes[i,j]>th[i,j]]
            over = r.mes[i,j]
            overx = x
            
            #Old code
            #under = mes[i,j][mes[i,j]<=th[i,j]]
            #underx = x[mes[i,j]<=th[i,j]]
            #pp.plot(x, th[i,j], 'r:', x, ic1[i,j], 'k:', x, ic2[i,j], 'k:', 
            #        overx, over, 'b-', underx, under, 'r-')
            
            pp.plot(x, r.th[i,j], 'r:', x, r.ic1[i,j], 'k:', x, r.ic2[i,j], 'k:', 
                    overx, over, 'b-')
            
            #Complicated code for underthreshold painting
            k = 0
            while(k < pr_.nf):
                while(r.mes[i,j,k] >= r.th[i,j,k]):
                    k = k+1
                    if (k == pr_.nf): break
                if (k == pr_.nf): break
                kold = k
                while(r.mes[i,j,k] < r.th[i,j,k]):
                    k = k+1
                    if (k == pr_.nf): break
                pp.plot(x[kold:k], r.mes[i,j,kold:k], 'r-')
            
            pp.ylim(-0.05,1.05)
            if (i < n-1):
                pp.xticks([])
            if (j > 0):
                pp.yticks([])
            
            if pr_.plotf != None:
                pp.xlim([0, pr_.plotf])
                
        if (pr_.ss):
            ax = pp.subplot(n,n,i*n+i+1).twinx()
            ax.plot(pr_.sample_f*arange(pr_.nf)/(2.0*pr_.nf), r.ss[i,i,:], color='g')
            if pr_.logss:
                ax.set_ylim(ymin = r.ss[i,i,:].min(), ymax = r.ss[i,i,:].max())
            else:
                ax.set_ylim(ymin = 0, ymax = r.ss[i,i,:].max())
                
            if (i < n-1):
                ax.set_xticks([])
            
            if pr_.plotf != None:
                ax.set_xlim(xmin = 0, xmax = pr_.plotf)
        
        pp.draw()
    #pp.show()
    
#pdc, ss = None, sample_f = 1.0, power = True, logss = False
def pdc_plot(res_, pr_):
    '''Plots nxn graphics. 
       If ss == True, plots ss in the diagonal.
       Expects data in complex form. Does: abs(x)^2 before plotting.'''
       
    r = res_.copy()
       
    n,n,nf = r.mes.shape
    if pr_.power:
        r.mes = r.mes*r.mes.conj()
    
    if pr_.logss and pr_.ss:
        r.ss = log(r.ss)
    
        #print r.ss.shape
        #print r.ss[0,0,:10]
    
    #pp.ion()
    
    for i in range(n):
        for j in range(n):
            pp.subplot(n,n,i*n+j+1)
            pp.plot(pr_.sample_f*arange(pr_.nf)/(2.0*nf), r.mes[i,j,:])
            pp.ylim(-0.05,1.05)
            if (i < n-1):
                pp.xticks([])
            if (j > 0):
                pp.yticks([])
        if (pr_.ss):
            ax = pp.subplot(n,n,i*n+i+1).twinx()
            ax.plot(pr_.sample_f*arange(pr_.nf)/(2.0*nf), r.ss[i,i,:], color='g')
            if pr_.logss:
                ax.set_ylim(ymin = r.ss[i,i,:].min(), ymax = r.ss[i,i,:].max())
            else:
                ax.set_ylim(ymin = 0, ymax = r.ss[i,i,:].max())
                
            if (i < n-1):
                ax.set_xticks([])
    
        pp.draw()
    #pp.show()

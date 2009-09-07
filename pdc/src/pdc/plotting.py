
from numpy import *
import matplotlib.pyplot as pp


def plot_all(mes, th, ic1, ic2, ss = None, nf = 64, sample_f = 1.0):
    '''Plots nxn graphics, with confidence intervals and threshold. 
       If ss == True, plots ss in the diagonal.'''
    x = sample_f*arange(nf)/(2.0*nf)
    n = mes.shape[0]
    for i in range(n):
        for j in range(n):
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
            
            pp.plot(x, th[i,j], 'r:', x, ic1[i,j], 'k:', x, ic2[i,j], 'k:', 
                    overx, over, 'b-')
            
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
            if (i < n-1):
                pp.xticks([])
            if (j > 0):
                pp.yticks([])
        if (ss != None):
            ax = pp.subplot(n,n,i*n+i+1).twinx()
            ax.plot(sample_f*arange(nf)/(2.0*nf), ss[i,i,:], color='g')
            ax.set_ylim(ymin = 0, ymax = ss[i,i,:].max())
            if (i < n-1):
                ax.set_xticks([])
    pp.show()
    
def pdc_plot(pdc, ss = None, nf = 64, sample_f = 1.0):
    '''Plots nxn graphics. 
       If ss == True, plots ss in the diagonal.'''
    n = pdc.shape[0]
    pdc = pdc*pdc.conj()
    for i in range(n):
        for j in range(n):
            pp.subplot(n,n,i*n+j+1)
            pp.plot(sample_f*arange(nf)/(2.0*nf), pdc[i,j,:])
            pp.ylim(-0.05,1.05)
            if (i < n-1):
                pp.xticks([])
            if (j > 0):
                pp.yticks([])
        if (ss != None):
            ax = pp.subplot(n,n,i*n+i+1).twinx()
            ax.plot(sample_f*arange(nf)/(2.0*nf), ss[i,i,:], color='g')
            ax.set_ylim(ymin = 0, ymax = ss[i,i,:].max())
            if (i < n-1):
                ax.set_xticks([])
    pp.show()

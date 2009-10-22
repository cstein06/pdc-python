
from numpy import *
from matplotlib import pyplot as pp
    
def plot_hbm09(mes, th, ic1, ic2, ss = None, nf = 64, sample_f = 1.0):
    
    x = sample_f*arange(nf)/(2.0*nf)
    n = mes.shape[0]
    for i in range(n):
        for j in range(n):
            if (i == j): continue
            pp.subplot(1,2,j+1)
            #over = mes[i,j][mes[i,j]>th[i,j]]
            #overx = x[mes[i,j]>th[i,j]]
            over = mes[i,j]
            overx = x
            under = mes[i,j][mes[i,j]<=th[i,j]]
            underx = x[mes[i,j]<=th[i,j]]
            pp.plot(x, th[i,j], 'r:', 
                    overx, over, 'b-', underx, under, 'r-')
            pp.ylim(-0.05,1.05)
            pp.ylabel ('PDC')
            pp.xlabel('Frequency (Hz)')
            if (j == 0):
                pp.title('Parietal -> Occipital')
            else:
                pp.title('Occipital -> Parietal')
    pp.show()

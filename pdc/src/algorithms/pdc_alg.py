
from numpy import *
import matplotlib.pyplot as pp
from scipy.linalg import inv

import cProfile

from data_simulation.ar_data import ar_data
import algorithms.ar_fit as ar_fit
import algorithms.asymp as as_

def list_to_array(data):
    d = data[0].reshape(1,-1)
    for i in range(1,len(data)):
        d = concatenate([d, data[i].reshape(1,-1)], axis = 0)
    return d

def pdc_ss_coh(data, maxp = 30, nf = 64, detrend = True):

    if(type(data) == 'list'):
        data = list_to_array(data)
        
    if (detrend):
        data = data - mean(data, axis = 1).reshape(-1,1) #TODO: usar signal.detrend?
    
    A, er = ar_fit.nstrand(data, maxp)
    return abs(pdc_alg(A, er, nf))**2, abs(ss_alg(A, er, nf))**2, abs(coh_alg(A, er, nf))**2


def pdc(data, maxp = 30, nf = 64, detrend = True, SS = True, metric = 'gen'):
    '''Generates spectral PDC matrix from data array
    
      Input: 
        data(n, m) - data matrix (n - number of signals, m - data length)
        maxp - maximum order for estimated AR model
        nf - frequency resolution
        detrend - Shall the data be detrended
        
      Output:
        PDC(n, n, nf) - PDC matrix
        ss(n, n, nf) - Parametric cross spectral matrix
    '''
    if(type(data) == 'list'):
        d = data[0].reshape(1,-1)
        for i in range(size(data)):
            d = concatenate([d, data[i].reshape(1,-1)], axis = 0)
        data = d
        
    if (detrend):
        data = data - mean(data, axis = 1).reshape(-1,1) #TODO: usar signal.detrend?
    A, er = ar_fit.nstrand(data, maxp)
    print A
    if (SS):
        return pdc_alg(A, er, nf, metric = metric), ss_alg(A, er, nf)
    else:
        return pdc_alg(A, er, nf, metric = metric)



def A_to_f(A, nf = 64):
    '''Calculates A(f), in the frequency domain
    
    Input:
        A(n, n, r) - recurrence matrix (n - number of signals, r - model order)
        nf - frequency resolution
        
    Output:
        AL(nf, n, n)
    '''
    
    n, n, r = A.shape
    
    # expoentes eh o array de expoentes da fft, com todas frequencias por todos lags
    expoentes = (-1j*pi*kron(arange(nf),(arange(r)+1.0))/nf).reshape(nf,r)
    # Af faz a multiplicacoes da exp(ar) pela matrix A, para todas frequencias
    # as funcoes repeat e transpose sao truques para possibilitar o calculo vetorial
    Af = (A.reshape(n,n,1,r).repeat(nf, axis=2)*exp(expoentes)).transpose([2,0,1,3])
    # o fft soma o valor para todos os lags
    AL = eye(n) - sum(Af, axis = 3)
    
    return AL
    
def pc_alg(A, e_cov, nf = 64):
    n, n, r = A.shape
    
    e_cov = mat(e_cov)
    AL = A_to_f(A, nf)
    pc = empty(AL.shape, dtype = 'complex')
    for i in range(nf):
        ALi = mat(AL[i])
        ps = ALi.T*e_cov.I*ALi.conj()
        d = ps.diagonal()
        m = kron(d,d).reshape(n,n)
        pc[i] = ps/sqrt(m)
    return pc.transpose(1,2,0)
    
def ss_alg(A, e_cov, nf = 64):
    n, n, r = A.shape
    
    AL = A_to_f(A, nf)
    ss = empty(AL.shape, dtype = 'complex')
    for i in range(nf):
        H = mat(AL[i]).I
        ss[i] = H*e_cov*H.T.conj()
    return ss.transpose(1,2,0)

def ss_coh():
    n, n, r = A.shape
    
    AL = A_to_f(A, nf)
    coh = empty(AL.shape, dtype = 'complex')
    ss = empty(AL.shape, dtype = 'complex')
    for i in range(nf):
        H = mat(AL[i]).I
        ss[i] = H*e_cov*H.T.conj()
        d = ss[i].diagonal()
        m = kron(d,d).reshape(n,n)
        coh[i] = ss[i]/sqrt(m)
    return ss.transpose(1,2,0), coh.transpose(1,2,0)

def coh_alg(A, e_cov, nf = 64):
    n, n, r = A.shape
    
    AL = A_to_f(A, nf)
    coh = empty(AL.shape, dtype = 'complex')
    for i in range(nf):
        H = mat(AL[i]).I
        ss = H*e_cov*H.T.conj()
        d = ss.diagonal()
        m = kron(d,d).reshape(n,n)
        coh[i] = ss/sqrt(m)
    return coh.transpose(1,2,0)

def coh(data, maxp = 30, nf = 64, detrend = True, SS = True):
    if(type(data) == 'list'):
        d = data[0].reshape(1,-1)
        for i in range(size(data)):
            d = concatenate(d, data[i].reshape(1,-1, axis = 0))
        data = d
        
    if (detrend):
        data = data - mean(data, axis = 1).reshape(-1,1) #TODO: usar signal.detrend?
    A, er = ar_fit.nstrand(data, maxp)
    return coh_alg(A,er,nf)


def pdc_alg(A, e_cov, nf = 64, metric = 'gen'):
    '''Generates spectral general (estatis. norm) PDC matrix from AR matrix
    
      Input: 
        A(n, n, r) - recurrence matrix (n - number of signals, r - model order)
        e_cov(n, n) - error covariance matrix
        nf - frequency resolution
    
      Output:
        PDC(n, n, nf) - PDC matrix
    '''
    
    n, n, r = A.shape
    if metric == 'euc':
        nornum = ones(n)
        norden = identity(n)
    elif metric == 'diag':
        nornum = 1/diag(e_cov)
        norden = diag(1/diag(e_cov))
    else: #metric == 'gen'
        nornum = 1/diag(e_cov)
        norden = inv(e_cov)
    
    AL = A_to_f(A, nf)
    
    ALT = AL.transpose([0,2,1])
    #dPDC = sum(dot(ALT,norden)*ALT.conj(), axis = -1).reshape(nf,-1)
    dPDC = sum(dot(ALT,norden)*ALT.conj(), axis = -1).reshape(nf,-1)
    nPDC = AL*sqrt(nornum).reshape(-1,1)
    PDC = nPDC/sqrt(abs(dPDC)).reshape(nf,1,n).repeat(n, axis = 1)
    return PDC.transpose(1,2,0)

def dtf_one_alg(A, er, nf = 64):
    '''Generates spectral non-norm. DTF matrix from AR matrix
    
      Input: 
        A(n, n, r) - recurrence matrix (n - number of signals, r - model order)
        e_cov(n, n) - error covariance matrix
        nf - frequency resolution
    
      Output:
        DTF(n, n, nf) - PDC matrix
    '''
    
    n, n, r = A.shape
    #nor = ones(n) # dtf_one nao tem normalizacao
    
    AL = A_to_f(A, nf)
    HL = empty(AL.shape, dtype=complex)
    for i in range(nf):
        HL[i] = inv(AL[i])
        
    # normalizacao por sum(ai ai* sig)
    dDTF = sum(HL*HL.conj(), axis = 2)
    nDTF = HL
    DTF = nDTF/sqrt(abs(dDTF)).reshape(nf,n,1).repeat(n, axis = 2)
    
    return DTF.transpose(1,2,0)

def plot_all(mes, th, ic1, ic2, nf = 64, sample_f = 1.0):
    
    x = sample_f*arange(nf)/(2.0*nf)
    n = mes.shape[0]
    for i in range(n):
        for j in range(n):
            pp.subplot(n,n,i*n+j+1)
            #over = mes[i,j][mes[i,j]>th[i,j]]
            #overx = x[mes[i,j]>th[i,j]]
            over = mes[i,j]
            overx = x
            under = mes[i,j][mes[i,j]<=th[i,j]]
            underx = x[mes[i,j]<=th[i,j]]
            pp.plot(x, th[i,j], 'r:', x, ic1[i,j], 'k:', x, ic2[i,j], 'k:', 
                    overx, over, 'b-', underx, under, 'r-')
            pp.ylim(-0.05,1.05)
            if (i < n-1):
                pp.xticks([])
            if (j > 0):
                pp.yticks([])
        #if (ss != None):
        #    ax = pp.subplot(n,n,i*n+i+1).twinx()
        #    ax.plot(sample_f*arange(nf)/(2.0*nf), ss[i,i,:], color='g')
        #    ax.set_ylim(ymin = 0, ymax = ss[i,i,:].max())
        #    if (i < n-1):
        #        ax.set_xticks([])
    pp.show()
    
def pdc_plot(pdc, ss = None, nf = 64, sample_f = 1.0):
    
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

def pdc_ass_and_plot(data, maxp = 5, nf = 64, sample_f = 1, 
                     ss = True, alpha = 0.05, metric = 'gen'):
 
    if(type(data) == type([])):
        data = list_to_array(data)
        
    #Estimate AR parameters with Nuttall-Strand
    Aest, erest = ar_fit.nstrand(data, maxp = maxp)
    print  'A:', Aest
    print 'evar:', erest
    #Calculate the connectivity and statistics
    mes, th, ic1, ic2 = as_.asymp_pdc(data, Aest, nf, erest, 
                                   maxp, alpha = alpha, metric = metric)
    
    plot_all(mes, th, ic1, ic2, nf = nf)

def pdc_and_plot(data, maxp = 30, nf = 64, sample_f = 1, ss = True, metric = 'gen'):
    
    if(type(data) == type([])):
        data = list_to_array(data)
    
    pdc_, ss_ = pdc(data, maxp, nf, metric = metric)
    if(not ss):
        ss_ = None
    pdc_plot(pdc_, ss_, nf, sample_f)
    


from numpy import *
import matplotlib.pyplot as pp
from scipy.linalg import inv
import scipy.signal as sig

import cProfile

from data_simulation.ar_data import ar_data
import algorithms.ar_fit as ar_fit
import algorithms.asymp as as_
from algorithms.plotting import *


def list_to_array(data):
    '''Converts a list to an array'''
    d = data[0].reshape(1,-1)
    for i in range(1,len(data)):
        d = concatenate([d, data[i].reshape(1,-1)], axis = 0)
    return d

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
    '''Calculates the Partial Coherence
        A -> autoregressive matrix
        e_cov -> residues
        nf -> number of frequencies
        '''
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
    '''Calculates the Spectral density (SS)
        A -> autoregressive matrix
        e_cov -> residues
        nf -> number of frequencies
        '''
    n, n, r = A.shape
    
    AL = A_to_f(A, nf)
    ss = empty(AL.shape, dtype = 'complex')
    for i in range(nf):
        H = mat(AL[i]).I
        ss[i] = H*e_cov*H.T.conj()
    return ss.transpose(1,2,0)

def ss_coh_alg(A, e_cov, nf = 64):
    '''Calculates the Spectral density (SS) and Coherence (coh)
        A -> autoregressive matrix
        e_cov -> residues
        nf -> number of frequencies
        '''
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
    '''Calculates the Coherence (coh)
        A -> autoregressive matrix
        e_cov -> residues
        nf -> number of frequencies
        '''
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
    '''Generates spectral not normalized DTF matrix from AR matrix
    
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


def pdc_ss_coh(data, maxp = 30, nf = 64, detrend = True):
    '''Interface that returns the PDC, SS and coh'''

    if(type(data) == 'list'):
        data = list_to_array(data)
        
    if (detrend):
        data = sig.detrend(data)
    
    A, er = ar_fit.ar_fit(data, maxp)
    return abs(pdc_alg(A, er, nf))**2, abs(ss_alg(A, er, nf))**2, abs(coh_alg(A, er, nf))**2


def pdc(data, maxp = 30, nf = 64, detrend = True, ss = True, metric = 'gen'):
    '''Generates spectral PDC matrix from data array
    
      Input: 
        data(n, m) - data matrix (n - number of signals, m - data length)
        maxp - maximum order for estimated AR model
        nf - frequency resolution
        detrend - Shall the data be detrended
        SS - Shall calculate the SS also
        metric - which PDC to use ('euc', 'diag' or 'gen')
        
      Output:
        PDC(n, n, nf) - PDC matrix
        ss(n, n, nf) - Parametric cross spectral matrix
    '''
    
    if(type(data) == type([])):
        data = list_to_array(data)
        
    if (detrend):
        data = sig.detrend(data) 
        
    A, er = ar_fit.nstrand(data, maxp)
    
    if (ss):
        return pdc_alg(A, er, nf, metric = metric), ss_alg(A, er, nf)
    else:
        return pdc_alg(A, er, nf, metric = metric)



def coh(data, maxp = 30, nf = 64, detrend = True, ss = True):
    '''Interface that calculate the Coherence from data'''
    
    if(type(data) == type([])):
        data = list_to_array(data)
        
    if (detrend):
        data = sig.detrend(data)
    A, er = ar_fit.ar_fit(data, maxp)
    
    if (ss):
        return coh_alg(A, er, nf), ss_alg(A, er, nf)
    else:
        return coh_alg(A, er, nf)


def dtf(data, maxp = 30, nf = 64, detrend = True, ss = True):
    '''Interface that calculate the Coherence from data'''
    
    if(type(data) == type([])):
        data = list_to_array(data)
        
    if (detrend):
        data = sig.detrend(data)
    A, er = ar_fit.ar_fit(data, maxp)
    
    
    if (ss):
        return dtf_one_alg(A, er, nf), ss_alg(A, er, nf)
    else:
        return dtf_one_alg(A, er, nf)

def ss(data, maxp = 30, nf = 64, detrend = True, ss = True):
    '''Interface that calculate the Coherence from data'''
    
    if(type(data) == type([])):
        data = list_to_array(data)
        
    if (detrend):
        data = sig.detrend(data)
    A, er = ar_fit.ar_fit(data, maxp)
    return ss_alg(A,er,nf)

def pc(data, maxp = 30, nf = 64, detrend = True, ss = True):
    '''Interface that calculate the Coherence from data'''
    
    if(type(data) == type([])):
        data = list_to_array(data)
        
    if (detrend):
        data = sig.detrend(data)
    A, er = ar_fit.ar_fit(data, maxp)
    
    if (ss):
        return pc_alg(A, er, nf), ss_alg(A, er, nf)
    else:
        return pc_alg(A, er, nf)

def pdc_ass_and_plot(data, maxp = 5, nf = 64, sample_f = 1, 
                     ss = True, alpha = 0.05, metric = 'gen', detrend = True):
    '''Interface that calculates PDC from data, calculates asymptotics statistics and plots everything.'''
    
    if(type(data) == type([])):
        data = list_to_array(data)
        
    if (detrend):
        data = sig.detrend(data)
        
    #Estimate AR parameters with Nuttall-Strand
    Aest, erest = ar_fit.ar_fit(data, maxp)
    print  'A:', Aest
    erest = (erest+erest.T)/2   #TODO: conferir isso.
    print 'evar:', erest
    #Calculate the connectivity and statistics
    mes, th, ic1, ic2 = as_.asymp_pdc(data, Aest, nf, erest, 
                                   maxp, alpha = alpha, metric = metric)
    if (ss == True):
        ssm = ss_alg(Aest, erest, nf)
    else:
        ssm = None
    
    plot_all(mes, th, ic1, ic2, nf = nf, ss = ssm, sample_f = sample_f)
    

def measure_ass_and_plot(data, measure, maxp = 5, nf = 64, sample_f = 1, 
                 ss = True, alpha = 0.05, detrend = True):
    '''Interface that calculates measure from data, calculates asymptotics statistics and plots everything.
       measure: 'dtf', 'coh', 'ss', 'pc'
       '''
    
    if(type(data) == type([])):
        data = list_to_array(data)
        
    if (detrend):
        data = sig.detrend(data)
        
    #Estimate AR parameters with Nuttall-Strand
    Aest, erest = ar_fit.ar_fit(data, maxp)
    print  'A:', Aest
    erest = (erest+erest.T)/2   #TODO: conferir isso.
    print 'evar:', erest
    #Calculate the connectivity and statistics
    if (measure == 'dtf'):
        mes, th, ic1, ic2 = as_.asymp_dtf_one(data, Aest, nf, erest, 
                                          maxp, alpha = alpha)
    if (measure == 'coh'):
        mes, th, ic1, ic2 = as_.asymp_coh(data, Aest, nf, erest, 
                                          maxp, alpha = alpha)
    if (measure == 'ss'):
        mes, th, ic1, ic2 = as_.asymp_ss(data, Aest, nf, erest, 
                                          maxp, alpha = alpha)
    if (measure == 'pc'):
        mes, th, ic1, ic2 = as_.asymp_pc(data, Aest, nf, erest, 
                                          maxp, alpha = alpha)
        
    if (ss == True):
        ssm = ss_alg(Aest, erest, nf)
    else:
        ssm = None
    
    plot_all(mes, th, ic1, ic2, nf = nf, ss = ssm, sample_f = sample_f)


def measure_and_plot(data, measure, maxp = 30, nf = 64, sample_f = 1, ss = True):
    '''Interface that calculates PDC from data and plots it'''
    if(type(data) == type([])):
        data = list_to_array(data)
        
    
    if (measure == 'dtf'):
        alg = dtf
    if (measure == 'coh'):
        alg = coh
    if (measure == 'ss'):
        alg = ss
    if (measure == 'pc'):
        alg = pc
    
    if (ss):
        mea, ss_ = alg(data, maxp, nf, ss = True)
    else:
        mea = alg(data, maxp, nf, ss = False)
        ss_ = None
        
    pdc_plot(mea, ss_, nf, sample_f)

def pdc_and_plot(data, maxp = 30, nf = 64, sample_f = 1, ss = True, metric = 'gen'):
    '''Interface that calculates PDC from data and plots it'''
    if(type(data) == type([])):
        data = list_to_array(data)
    
    pdc_, ss_ = pdc(data, maxp, nf, metric = metric)
    if(not ss):
        ss_ = None
    pdc_plot(pdc_, ss_, nf, sample_f)
    

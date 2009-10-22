
from numpy import *
import matplotlib.pyplot as pp
from scipy.linalg import inv
import scipy.signal as sig
from scipy.stats import f

import cProfile

from pdc.ar_data import ar_data
import pdc.ar_fit as ar_fit
import pdc.asymp as as_
from pdc.plotting import *
import pdc.bootstrap as bt_

def list_to_array(data):
    '''Converts a list to an array'''
    d = data[0].reshape(1,-1)
    for i in range(1,len(data)):
        d = concatenate([d, data[i].reshape(1,-1)], axis = 0)
    return d

def pre_data(data, normalize = True, detrend = True):
    if (detrend):
        data = sig.detrend(data)
        
    if (normalize):
        data = data/std(data, axis = 1).reshape(-1,1)
        
    return data

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

def dtf_alg(A, er, nf = 64):
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
        return dtf_alg(A, er, nf), ss_alg(A, er, nf)
    else:
        return dtf_alg(A, er, nf)

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
    

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%
#%  Computes granger causality index
#%
#%  Input: 
#%    D(n, N) - data (n channels)
#%    MaxIP - externaly defined maximum IP
#%
#%  Output:
#%    Gr(n, n) - Granger causalit index
#%    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#function [Gr] = alg_ganger(u, maxIP)
#
#[n N] = size(u);
#
#[IP,pf,A,pb,B,ef,eb,vaic,Vaicv] = mvar(u,maxIP,[0 0]);
#
#va = diag(pf);
#
#va_n = zeros(n, n);
#
#for iu = 1:n
#  aux_u = u;
#  aux_u(iu,:) = [];
#  [IP,pf,A,pb,B,ef,eb,vaic,Vaicv] = mvar(aux_u,maxIP,[0 0]);
#  aux = diag(pf)';
#  va_n(iu,:) = cat(2, aux(1:iu-1), 0, aux(iu:n-1));
#end
#
#Gr = zeros(n, n);
#for iu = 1:n
#  for ju = 1:n
#    if (iu == ju) continue; end
#    Gr(iu,ju) = log(va_n(ju,iu)/va(iu));
#  end
#end

def gci(data, maxp = 30, detrend = True):
    
    n = data.shape[0]
    
    if (detrend):
        data = sig.detrend(data)
        
    A0, er0 = ar_fit.ar_fit(data, maxp)
    va0 = diag(er0)
    
    gci = zeros([n,n])
    for i in arange(n): 
        aux_data = delete(data, i, 0)
        A1, er1 = ar_fit.ar_fit(aux_data, maxp)
        va1 = diag(er1) 
        va1 = insert(va1, i, 0)
        gci[:,i] = log(float64(va1)/va0)
        
    return gci

#def gct(data, maxp = 30, detrend = True):
#    #TODO: esta errado.
#    n,T = data.shape
#    
#    if (detrend):
#        data = sig.detrend(data)
#        
#    A0, er0 = ar_fit.ar_fit(data, maxp)
#    va0 = diag(er0)
#    
#    p = A0.shape[2] #TODO: p pode variar depois. fixar para A1?
#    print p 
#    
#    gci = zeros([n,n])
#    for i in arange(n): 
#        aux_data = delete(data, i, 0)
#        A1, er1 = ar_fit.ar_fit(aux_data, maxp)
#        va1 = float64(diag(er1)) 
#        va1 = insert(va1, i, 0)
#        gci[:,i] = ((va1-va0)/(n*p))/(va0/(T-n*p-1))
#    
#    gct = f.cdf(gci, n*p, T-n*p-1)
#        
#    return gct
        

def pdc_full(data, maxp = 5, nf = 64, sample_f = 1, 
                   ss = True, alpha = 0.05, metric = 'gen', 
                   detrend = True, normalize = False, stat = 'asymp', n_boot = 1000):
    '''Interface that calculates PDC from data, calculates asymptotics statistics and plots everything.'''
    
    if(type(data) == type([])):
        data = list_to_array(data)
    
    n,nd = data.shape
        
    data = pre_data(data, normalize, detrend)
        
    #Estimate AR parameters with Nuttall-Strand
    Aest, erest = ar_fit.ar_fit(data, maxp)
    print  'A:', Aest
    #erest = (erest+erest.T)/2   #TODO: conferir isso. porque nao eh sempre simetrico?
    print 'evar:', erest
    #Calculate the connectivity and statistics
    if stat == 'asymp':
        mes, th, ic1, ic2 = as_.asymp_pdc(data, Aest, nf, erest, 
                                   maxp, alpha = alpha, metric = metric)
    elif stat == 'boot':
        mes, th, ic1, ic2 = bt_.bootstrap(pdc_alg, nd, n_boot, Aest, erest, 
                                          nf, alpha = alpha, metric = metric)
    else:
        mes = pdc_alg(Aest, erest, nf, metric)
        th = zeros(mes.shape)
        ic1 = zeros(mes.shape)
        ic2 = zeros(mes.shape)
         
    if (ss == True):
        ssm = ss_alg(Aest, erest, nf)
    else:
        ssm = None
        
    plot_all(mes, th, ic1, ic2, nf = nf, ss = ssm, sample_f = sample_f)
    
def coh_full(data, maxp = 5, nf = 64, sample_f = 1, 
             ss = True, alpha = 0.05, detrend = True, normalize = False, stat = 'asymp', n_boot = 1000):
    measure_full(data, 'coh', maxp, nf, sample_f, ss, alpha, detrend, normalize, stat = stat, n_boot = n_boot)

def dtf_full(data, maxp = 5, nf = 64, sample_f = 1, 
             ss = True, alpha = 0.05, detrend = True, normalize = False, stat = 'asymp', n_boot = 1000):
    measure_full(data, 'dtf', maxp, nf, sample_f, ss, alpha, detrend, normalize, stat = stat, n_boot = n_boot)

def ss_full(data, maxp = 5, nf = 64, sample_f = 1, 
             ss = True, alpha = 0.05, detrend = True, normalize = False, stat = 'asymp', n_boot = 1000):
    measure_full(data, 'ss', maxp, nf, sample_f, ss, alpha, detrend, normalize, stat = stat, n_boot = n_boot)

def pc_full(data, maxp = 5, nf = 64, sample_f = 1, 
             ss = True, alpha = 0.05, detrend = True, normalize = False, stat = 'asymp', n_boot = 1000):
    measure_full(data, 'pc', maxp, nf, sample_f, ss, alpha, detrend, normalize, stat = stat, n_boot = n_boot)

def measure_full(data, measure, maxp = 5, nf = 64, sample_f = 1, 
                 ss = True, alpha = 0.05, detrend = True, 
                 normalize = False, stat = 'asymp', n_boot = 1000):
    '''Interface that calculates measure from data, calculates asymptotics statistics and plots everything.
       measure: 'dtf', 'coh', 'ss', 'pc'
       '''
    
    if(type(data) == type([])):
        data = list_to_array(data)
    
    n,nd = data.shape
        
    data = pre_data(data, normalize, detrend)
        
    #Estimate AR parameters with Nuttall-Strand
    Aest, erest = ar_fit.ar_fit(data, maxp)
    print  'A:', Aest
    #erest = (erest+erest.T)/2   #TODO: conferir isso.
    print 'evar:', erest
    #Calculate the connectivity and statistics
    
    
    if stat == 'asymp': 
        if (measure == 'dtf'):
            mes, th, ic1, ic2 = as_.asymp_dtf(data, Aest, nf, erest, 
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
    elif stat == 'boot':
        methcall = globals()[measure + '_alg']
        mes, th, ic1, ic2 = bt_.bootstrap(methcall, nd, n_boot, Aest, erest, 
                                          nf, alpha = alpha)
    else:
        methcall = var()[measure + '_alg']
        mes = methcall(Aest, erest, nf)
        th = zeros(mes.shape)
        ic1 = zeros(mes.shape)
        ic2 = zeros(mes.shape)
        
    if (ss == True):
        ssm = ss_alg(Aest, erest, nf)
    else:
        ssm = None
    
    plot_all(mes, th, ic1, ic2, nf = nf, ss = ssm, sample_f = sample_f)


def pdc_and_plot(data, maxp = 30, nf = 64, sample_f = 1, ss = True, metric = 'gen'):
    '''Interface that calculates PDC from data and plots it'''
    if(type(data) == type([])):
        data = list_to_array(data)
    
    pdc_, ss_ = pdc(data, maxp, nf, metric = metric)
    if(not ss):
        ss_ = None
    pdc_plot(pdc_, ss_, nf, sample_f)

def coh_and_plot(data, maxp = 30, nf = 64, sample_f = 1, ss = True):
    measure_and_plot(data, 'coh', maxp, nf, sample_f, ss)
    
def dtf_and_plot(data, maxp = 30, nf = 64, sample_f = 1, ss = True):
    measure_and_plot(data, 'dtf', maxp, nf, sample_f, ss)
    
def ss_and_plot(data, maxp = 30, nf = 64, sample_f = 1, ss = True):
    measure_and_plot(data, 'ss', maxp, nf, sample_f, ss)
    
def pc_and_plot(data, maxp = 30, nf = 64, sample_f = 1, ss = True):
    measure_and_plot(data, 'pc', maxp, nf, sample_f, ss)
    
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
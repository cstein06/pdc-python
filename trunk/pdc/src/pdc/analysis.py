# -*- coding:utf-8 -*-

__all__ = ['pre_data', 'A_to_f', 'measure_full', 'measure_and_plot', 'measure',
           'pdc', 'coh', 'dtf', 'pc', 'ss', 
           'pdc_full', 'coh_full', 'dtf_full', 'pc_full', 'ss_full',
           'pdc_and_plot', 'coh_and_plot', 'dtf_and_plot', 'pc_and_plot', 'ss_and_plot', 
           'gci', 'gct', 'igct', 'white_test']

from numpy import *
import matplotlib.pyplot as pp
from scipy.linalg import inv
import scipy.signal as sig
from scipy.stats import f

import cProfile
import time

from pdc.ar_data import ar_data
import pdc.ar_fit as ar_fit
import pdc.asymp as as_
import pdc.plotting as pl_
import pdc.bootstrap as bt_
from pdc.globals import *
from pdc.globals import mnames_

def list_to_array(data):
    '''Converts a list to an array'''
    d = data[0].reshape(1,-1)
    for i in range(1,len(data)):
        d = concatenate([d, data[i].reshape(1,-1)], axis = 0)
    return d

def pre_data(data, normalize = False, detrend = True):
    
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

def pc_alg(A, e_cov, nf = 64, metric = 'dummy'):
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
    
def ss_alg(A, e_cov, nf = 64, metric = 'dummy'):
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
        
    #print ss[5]
    return ss.transpose(1,2,0)

def ss_coh_alg(A, e_cov, nf = 64, metric = 'dummy'):
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

def coh_alg(A, e_cov, nf = 64, metric = 'dummy'):
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
    
    #print metric
    
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

def dtf_alg(A, er, nf = 64, metric = 'dummy'):
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


#def pdc_ss_coh(data, maxp = 30, nf = 64, detrend = True):
#    '''Interface that returns the PDC, SS and coh'''
#
#    if(type(data) == 'list'):
#        data = list_to_array(data)
#        
#    data = pre_data(data, pr_.normalize, pr_.detrend)
#    
#    A, er = ar_fit.ar_fit(data, maxp)
#    return abs(pdc_alg(A, er, nf))**2, abs(ss_alg(A, er, nf))**2, abs(coh_alg(A, er, nf))**2


def arfit(data, **args):
    '''Interface that calculate the Coherence from data'''
    read_args(args)
    
    return measure(data, alg = 'ar')

#        
#def read_args(args):
#    
#    for a,b in args.items():
#        globals()[a+'_'] = b


#maxp = 30, nf = 64, detrend = True, normalize = False, 
#        fixp = False, ss = True, metric = 'diag', power = False
def pdc(data, **args):
    '''Generates spectral PDC matrix from data array
    
      Input: 
        data(n, m) - data matrix (n - number of signals, m - data length)
        maxp - maximum order for estimated AR model
        nf - frequency resolution
        detrend - Shall the data be detrended
        SS - Shall calculate the SS also
        metric - which PDC to use ('euc', 'diag' or 'gen')
        power - returns abs(PDC)^2
        
      Output:
        PDC(n, n, nf) - PDC matrix
        ss(n, n, nf) - Parametric cross spectral matrix
    '''
    read_args(args)
    
    return measure(data, alg = 'pdc')

def coh(data, **args):
    '''Interface that calculate the Coherence from data'''
    read_args(args)
    
    return measure(data, alg = 'coh')
    
def dtf(data, **args):
    '''Interface that calculate the Dtf from data'''
    read_args(args)
    
    return measure(data, alg = 'dtf')
    

def ss(data, **args):
    '''Interface that calculate the SS from data'''
    read_args(args)
    
    return measure(data, alg = 'ss')
    

def pc(data, **args):
    '''Interface that calculate the PC from data'''
    read_args(args)
    
    return measure(data, alg = 'pc')
    
#maxp = 30, nf = 64, detrend = True, normalize = False, 
#        fixp = False, ss = True, metric = 'diag', power = False
def measure(data, **args):
    '''Generates spectral measure from data array
    
      Input: 
        data(n, m) - data matrix (n - number of signals, m - data length)
        maxp - maximum order for estimated AR model
        nf - frequency resolution
        detrend - Shall the data be detrended
        SS - Shall calculate the SS also
        metric - which PDC to use ('euc', 'diag' or 'gen')
        power - returns abs(PDC)^2
        
      Output:
        PDC(n, n, nf) - PDC matrix
        ss(n, n, nf) - Parametric cross spectral matrix
    ''' 
    
    read_args(args)
    
    auxp = pr_.do_plot
    pr_.do_plot = False
    
    auxs = pr_.stat
    pr_.stat = 'no'
    
    ret = measure_full(data)[0]
    
    pr_.stat = auxs
    pr_.do_plot = auxp
    
    return ret
#    
#    if(type(data) == type([])):
#        data = list_to_array(data)
#    
#    data = pre_data(data, pr_.normalize, pr_.detrend)
#        
#    crit = 0 #AIC
#    if pr_.fixp:
#        crit = 1
#        
#        
#    if pr_.v:
#        print 'Will calculate the', mnames_[pr_.alg], 'of the data'
#        print 'Dimensions of the data:', data.shape
#        if pr_.fixp:
#            print 'Using fixed VAR order'
#        else:
#            print 'Using AIC for VAR order' 
#        
#        print '\nEstimating VAR'
#    
#    
#    res_.A, res_.er = ar_fit.ar_fit(data, pr_.maxp, criterion=crit)
#    
#    
#    if pr_.v:
#        print '\nVAR estimaded. Order:', res_.A.shape
#    #print 'data:', data.shape
#    #print 'A:', res_.A.shape
#    #print er
#    
#    #print A
#    
#    method = globals()[pr_.alg + '_alg']
#    
#    res_.mes = method(res_.A, res_.er, nf = pr_.nf, metric = pr_.metric)
#    
#    
#    if pr_.v:
#        print '\n', mnames_[pr_.alg], 'estimaded.'
#        if pr_.alg == 'pdc':
#            print 'Used metric:', pr_.metric
#        print '\nNumber of frequencies:', pr_.nf
#        print 'From 0 to ', pr_.sample_f/2.0, 'Hz\n'
#    
#    
#    if pr_.power:
#        if pr_.v:
#            print 'Calculates squared power of the measure\n'
#        res_.mes = (res_.mes*res_.mes.conj()).real #todo: check power consistency of everything
#    
#    if pr_.ss:
#        if pr_.v:
#            print 'Calculates also the spectrum\n'
#        res_.ss = ss_alg(res_.A, res_.er, pr_.nf)
#    
#    res_.data_shape = data.shape
#    res_.alg = pr_.alg
#    
#    
#    if pr_.do_log:
#        if pr_.v:
#            print 'Logging the results in file:', res_.log_file   
#        pr_.time = time.ctime()
#        log_results()
#    
#    if pr_.ss:
#        return res_.mes, res_.ss
#    else:
#        return res_.mes
#    

        

#maxp = 5, nf = 64, sample_f = 1, 
#ss = True, alpha = 0.05, detrend = True, normalize = False, 
#stat = 'asymp', n_boot = 1000, fixp = False, metric = None   
def pdc_full(data, **args):
    read_args(args)
    return measure_full(data, alg = 'pdc')
        
def coh_full(data, **args):
    read_args(args)
    return measure_full(data, alg = 'coh')

def dtf_full(data, **args):
    read_args(args)
    return measure_full(data, alg = 'dtf')
    
def ss_full(data, **args):
    read_args(args)
    return measure_full(data, alg = 'ss')
    
def pc_full(data, **args):
    read_args(args)
    return measure_full(data, alg = 'pc')

def measure_full(data, **args):
    '''Interface that calculates some measure from data, calculates asymptotics statistics and plots everything.
       
       data 
       alg: 'pdc', 'dtf', 'coh', 'ss', 'pc'
       
    Possible parameters:
       maxp = 30, nf = 64, sample_f = 1, 
       ss = True, alpha = 0.05, metric = 'diag', 
       detrend = True, normalize = False, 
       stat = 'asymp', n_boot = 1000, fixp = False,
       plotf = None, *
    '''
    
    read_args(args)
    
    pr_.power == True
    
    if(type(data) == type([])):
        data = list_to_array(data)
    
    n,nd = data.shape
        
    data = pre_data(data, pr_.normalize, pr_.detrend)
        
    #Estimate AR parameters with Nuttall-Strand
    crit = 0 #AIC
    if pr_.fixp:
        crit = 1
    
    if pr_.v:
        print 'Will calculate the', mnames_[pr_.alg], 'of the data'
        print 'Dimensions of the data:', data.shape
        if pr_.fixp:
            print 'Using fixed VAR order'
        else:
            print 'Using AIC for VAR order' 
        print 'Error type I used:', pr_.alpha*100, '%'
        
        print '\nEstimating VAR'
    
    #Estimate AR parameters with Nuttall-Strand
    res_.A, res_.er = ar_fit.ar_fit(data)
    
    res_.p = res_.A.shape[2]
    
    #print res_.A, res_.er
         
    if pr_.v:
        print '\nVAR estimaded. Order:', res_.A.shape
        
    if pr_.alg == 'ar':
        return [(res_.A, res_.er)]
   
    #print  'A:', Aest
    #erest = (erest+erest.T)/2   #TODO: conferir isso.
    #print 'evar:', erest
    
    #Calculate the connectivity and statistics
    
    if pr_.stat == 'asymp': 
        if pr_.v:
            print 'Calculating asymptotic statistics'
        as_method = vars(as_)['asymp_' + pr_.alg]
        res_.mes, res_.th, res_.ic1, res_.ic2 = \
            as_method(data, res_.A, pr_.nf, res_.er, 
                      res_.p, alpha = pr_.alpha, metric = pr_.metric)
    elif pr_.stat == 'boot':
        if pr_.v:
            print 'Calculating bootstrap statistics'
        alg_method = globals()[pr_.alg + '_alg']
        res_.mes, res_.th, res_.ic1, res_.ic2 = \
            bt_.bootstrap(alg_method, nd, pr_.n_boot, res_.A, res_.er, 
            pr_.nf, alpha = pr_.alpha, metric = pr_.metric)
    else:
        if pr_.v:
            print 'Choosing no statistics'
        alg_method = globals()[pr_.alg + '_alg']
        res_.mes = alg_method(res_.A, res_.er, pr_.nf, metric = pr_.metric)
        #th = zeros(mes.shape)
        #ic1 = zeros(mes.shape)
        #ic2 = zeros(mes.shape)
    
        if pr_.power:
            if pr_.v:
                print 'Calculates squared power of the measure\n'
            res_.mes = (res_.mes*res_.mes.conj()).real #todo: check power consistency of everything
    
    
    if pr_.v:
        print '\n', mnames_[pr_.alg], 'estimaded'
        if pr_.alg == 'pdc':
            print 'Used metric:', pr_.metric
        print '\nNumber of frequencies:', pr_.nf
        print 'From 0 to ', (pr_.sample_f/2.0)*(pr_.nf-1.0)/pr_.nf, 'Hz\n'
        
        
    if (pr_.ss == True):
        if pr_.v:
            print 'Calculates also the spectrum\n'
        ssm = ss_alg(res_.A, res_.er, nf = pr_.nf)
    else:
        ssm = None
        
    res_.data_shape = data.shape
    res_.alg = pr_.alg
    res_.ss = ssm
     
    if pr_.do_log:   
        log_results()
    
    if pr_.do_plot:
        pl_.plot_all()
        
        
    if pr_.v:
        print '\nDone!\n'
    
    return res_.mes, res_.th, res_.ic1, res_.ic2

#data, maxp = 30, nf = 64, sample_f = 1, ss = True, metric = 'gen',
#                 detrend = True, normalize = False, fixp = False, logss = False
def pdc_and_plot(data, **args):
    read_args(args)
    return measure_and_plot(data, alg = 'pdc')
    
#maxp = 30, nf = 64, sample_f = 1, ss = True,
#detrend = True, normalize = False, fixp = False, power = True):
def coh_and_plot(data, **args):
    read_args(args)
    return measure_and_plot(data, alg = 'coh')
    
def dtf_and_plot(data, **args):
    read_args(args)
    return measure_and_plot(data, alg = 'coh')
    
def ss_and_plot(data, **args):
    read_args(args)
    return measure_and_plot(data, alg = 'coh')
        
def pc_and_plot(data, **args):
    read_args(args)
    return measure_and_plot(data, alg = 'coh')
    
def measure_and_plot(data, **args):
    '''Interface that calculates alg from data and plots it'''
    
    read_args(args)
    
    aux = pr_.stat
    pr_.stat = 'no'
    
    measure_full(data)
    
    pr_.stat = aux
    
#    if(type(data) == type([])):
#        data = list_to_array(data)
#    
#    method = globals()[pr_.alg]
#        
#    method(data)
#    
#    pl_.plot_all(plot_ic = False, plot_th = False)
#    #pl_.pdc_plot()
#    
    
    

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

def gci(data, **args):
    
    read_args(args)
    
    if(type(data) == type([])):
        data = list_to_array(data)
        
    n = data.shape[0]
    
    data = pre_data(data, pr_.normalize, pr_.detrend)
        
    A0, er0 = ar_fit.ar_fit(data)
    va0 = diag(er0)
    
    gci = zeros([n,n])
    for i in arange(n): 
        aux_data = delete(data, i, 0)
        A1, er1 = ar_fit.ar_fit(aux_data)
        va1 = diag(er1) 
        va1 = insert(va1, i, 0)
        gci[:,i] = log(float64(va1)/va0)
        
    return gci


def gct(data,**args):
    '''Asymptotic statistics for Wald statistic of the GC in time
    '''
    
    read_args(args)
    
    if(type(data) == type([])):
        data = list_to_array(data)
    
    data = pre_data(data, pr_.normalize, pr_.detrend)
    
    A, e_var = ar_fit.ar_fit(data)
        
    return as_.asymp_gct(data, A, e_var)

def igct(data, **args):
    '''Asymptotic statistics for Wald statistic of instantaneous GC '''

    read_args(args)
    
    if(type(data) == type([])):
        data = list_to_array(data)
    
    data = pre_data(data, pr_.normalize, pr_.detrend)
    
    A, e_var = ar_fit.ar_fit(data)
    
    n, nd = data.shape
        
    return as_.asymp_igct(e_var, nd)

def white_test(data, h = 20, **args):
    
    read_args(args)
    
    A, res = ar_fit.ar_fit(data, return_ef=True)
    
    p = A.shape[2]
    
    return as_.asymp_white(data, res, p, h)



#def gct(data, maxp = 30, detrend = True):
#    #TODO: esta errado, apagar.
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

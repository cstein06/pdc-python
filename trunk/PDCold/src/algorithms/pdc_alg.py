
from numpy import *
import matplotlib.pyplot as pp

import cProfile

from rpy2.robjects import r as r_
import rpy2.rinterface as ri_
import rpy2.robjects as ro_

from data_simulation.ar_data import ar_data


def pdc(data, maxp = 30, nf = 64, detrend = True, SS = True):
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
            d = concatenate(d, data[i].reshape(1,-1, axis = 0))
        data = d
        
    if (detrend):
        data = data - mean(data, axis = 1).reshape(-1,1) #TODO: usar signal.detrend?
    A, er = ar_fit(data, maxp)
    if (SS):
        return pdc_alg(A, er, nf), ss(A, er, nf)
    else:
        return pdc_alg(A, er, nf)

def ar_fit(data, maxp = 30):
    '''Estimates multivariate AR fit for data
    
      Input: 
        data(n, m) - data matrix (n - number of signals, m - data length)
        maxp - maximum order for estimated AR model
    
      Output:
        A(n, n, p) - estimated AR model
        er(n, n) - covariance of residuals
    '''
    
    ri_.initr()
    
    ri_.globalEnv["data"] = ri_.FloatSexpVector(data.ravel())
    ri_.globalEnv["dim"] = ri_.IntSexpVector(data.shape[::-1])
    ro_.globalEnv["maxp"] = maxp
    r_('data <- ar(array(data, dim), order.max = maxp)')
    
    A = array(r_('data$ar')).transpose(1,2,0) #TODO: conferir A e A.T em todos.
    er = array(r_('cov(data$resid[-seq(data$order),])'))
    #print 'Model order: ', array(r_('data$order'))[0]
    return A, er

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
    
def ss(A, e_cov, nf = 64):
    n, n, r = A.shape
    
    AL = A_to_f(A, nf)
    ss = empty(AL.shape, dtype = 'complex')
    for i in range(nf):
        H = mat(AL[i]).I
        ss[i] = H*e_cov*H.T
    return ss.transpose(1,2,0)

def coh_alg(A, e_cov, nf = 64):
    n, n, r = A.shape
    
    AL = A_to_f(A, nf)
    coh = empty(AL.shape, dtype = 'complex')
    for i in range(nf):
        H = mat(AL[i]).I
        ss = H*e_cov*H.T
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
    A, er = ar_fit(data, maxp)
    return coh_alg(A,er,nf)
    
def pdc_alg(A, e_cov, nf = 64):
    '''Generates spectral PDC matrix from AR matrix
    
      Input: 
        A(n, n, r) - recurrence matrix (n - number of signals, r - model order)
        e_cov(n, n) - error covariance matrix
        nf - frequency resolution
    
      Output:
        PDC(n, n, nf) - PDC matrix
    '''
    # TODO pdc esta errado. normalizacao nao esta OK.
    
    n, n, r = A.shape
    nor = 1/e_cov.diagonal() # normaliza pelo inverso da diagonal
    
    AL = A_to_f(A, nf)
    
    # normalizacao por sum(ai ai* sig)
    dPDC = dot(nor,AL*AL.conj())
    nPDC = AL*sqrt(nor).reshape(-1,1)
    PDC = nPDC/sqrt(abs(dPDC)).reshape(nf,1,n).repeat(n, axis = 1)
    
    return PDC.transpose(1,2,0)

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

def pdc_and_plot(data, maxp = 30, nf = 64, sample_f = 1, ss = True):
    pdc_, ss_ = pdc(data, maxp, nf)
    if(not ss):
        ss_ = None
    pdc_plot(pdc_, ss_, nf, sample_f)

def teste():
    A = array([[4,3],[0,3]], dtype=float).reshape(2,2,1)/10
    er = array([2,5], dtype = float)
    print 'A:', A
    print 'er:', er
    data = ar_data(A, er, 100000)
    Aest, er_est = ar_fit(data, maxp = 10)
    print 'A estimado:', Aest # ver dimensoes do A
    print 'er estimado:', er_est
    pdc = pdc_alg(Aest,er_est)
    pdc_plot(pdc)
    
def teste2():
    A = array([[4,3],[0,3]], dtype=float).reshape(2,2,1)/10
    er = array([2,5], dtype = float)
    print 'A:', A
    print 'er:', er
    data = ar_data(A, er, 1000)
    p, ss = pdc(data)
    pdc_plot(p, ss)
    
def teste_vel():
    a = random.randn(100)
    b = random.randn(100,100)
    for i in range(100):
        pdc([a,b[i]])

if __name__ == "__main__":
    pass
    # TODO: testar qdo n == r, pode haver conversao automatica falha
    
    # Testa ar_data
#    A = array([[4,0],[0,6]], dtype=float).reshape(2,2,1)/10
#    er = array([2,5], dtype = float)
#    data = ar_data(A, er, 10000)
    #print var(data, axis=1), mean(data, axis=1), 'var teorico:', er[0]/(1-0.4**2), er[1]/(1-0.6**2), 
    
    # Testa Ar_fit
#    A = ar_fit(random.randn(2,1000), maxp = 2)[1:]
#    print A.shape
    
    # Testa PDC
#    A = array([[4,0],[4,4]], dtype=float).reshape(2,2,1)/10
#    e = array([[1,0],[0,1]], dtype=float)
#    pdc = pdc_alg(A,e)
#    pdc_plot(pdc)

    # Teste completo
    A = array([[4,3],[0,3]], dtype=float).reshape(2,2,1)/10
    er = array([2,5], dtype = float)
#    print 'A:', A
#    print 'er:', er
    data = ar_data(A, er, 100)
    Aest, er_est = ar_fit(data, maxp = 1)
#    print 'A estimado:', Aest # ver dimensoes do A
#    print 'er estimado:', er_est
#    pdc = pdc_alg(Aest,er_est)
#    pdc_plot(pdc)
    
    #cProfile.run('teste_vel()', sort = 'time')
    #teste2()
    
    print coh(Aest, er_est)
    
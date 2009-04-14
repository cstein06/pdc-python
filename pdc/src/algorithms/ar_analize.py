from numpy import *
import matplotlib.pyplot as pp

import cProfile

from rpy2.robjects import r as r_
import rpy2.rinterface as ri_
import rpy2.robjects as ro_

from data_simulation.ar_data import ar_data

class Ar_analyze:
    
    def __init__(self, data, maxp = 30, nf = 64):
        '''Ar analysis
            maxp - maximum order for estimated AR model
            nf - frequency resolution'''
            
        self.data_ = data
        self.maxp_ = maxp
        self.nf_ = nf
        self.A_ = None
        self.er_ = None
        self.pdc_ = None
        self.dtf_ = None
        self.ss_ = None
    
    def pdc(self):
        '''Generates spectral PDC matrix from data array
        
          Input: 
            data(n, m) - data matrix (n - number of signals, m - data length)
            
          Output:
            PDC(n, n, nf) - PDC matrix
        '''
        
        if(self.pdc_ != None):
            return self.pdc_
        if(self.A_ == None):
            ar_fit()
        pdc_alg()
        return self.pdc_
    
    def dtf(data, ):
        pass
    
    def A_to_f(A, nf = 64):
        '''Calculates A(f), in the frequency domain
        
        Input:
            A(n, n, r) - recurrence matrix (n - number of signals, r - model order)
            nf - frequency resolution
            
        Output:
            AL(nf, n, n)
        '''
        
        n, n, r = A.shape
        
        # expoentes eh o array de expoentes do fft, com todas frequencias por todos lags
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
        ss = empty(AL.shape)
        for i in range(nf):
            H = mat(AL[i]).I
            ss[i] = H*e_cov*H.T
        return ss.transpose(1,2,0)
        
    def pdc_alg(A, e_cov, nf = 64):
        '''Generates spectral PDC matrix from AR matrix
        
          Input: 
            A(n, n, r) - recurrence matrix (n - number of signals, r - model order)
            e_cov(n, n) - error covariance matrix
            nf - frequency resolution
        
          Output:
            PDC(n, n, nf) - PDC matrix
        '''
        
        n, n, r = A.shape
        nor = 1/e_cov.diagonal() # normaliza pelo inverso da diagonal
        
        if(self.AL_ == None):
            self.AL_ = A_to_f(A, nf)
        AL = self.AL_
        
        # normalizacao por sum(ai ai* sig)
        dPDC = dot(nor,AL*AL.conj())
        nPDC = AL*sqrt(nor).reshape(-1,1)
        PDC = nPDC/sqrt(abs(dPDC)).reshape(nf,1,n).repeat(n, axis = 1)
        
        return PDC.transpose(1,2,0)
    
    def ar_fit(self, maxp = 30):
        '''Estimates multivariate AR fit for data
        
          Input: 
            data(n, m) - data matrix (n - number of signals, m - data length)
            maxp - maximum order for estimated AR model
        
          Output:
            A(n, n, p) - estimated AR model
        '''
        
        ri_.initr()
        
        ri_.globalEnv["data"] = ri_.FloatSexpVector(data.ravel())
        ri_.globalEnv["dim"] = ri_.IntSexpVector(data.shape[::-1])
        ro_.globalEnv["maxp"] = maxp
        r_('data <- ar(array(data, dim), order.max = maxp)')
        
        self.A = array(r_('data$ar')).transpose(2,1,0) #conferir A e A.T em todos.
        self.er = array(r_('cov(data$resid[-1,])'))
        return self.A, self.er
    
        


if __name__ == "__main__":
    A = array([[4,3],[0,3]], dtype=float).reshape(2,2,1)/10
    er = array([2,5], dtype = float)
    print 'A:', A
    print 'er:', er
    data = ar_data(A, er, 100000)
    p, ss = pdc(data)
    pdc_plot(p, ss)
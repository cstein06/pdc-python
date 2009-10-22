from numpy import *

from rpy2.robjects import r as r_
import rpy2.rinterface as ri_
import rpy2.robjects as ro_

def ar_data(A, er = None, m = 1000):
    '''Simulate ar-model from A matrix
    
      Input: 
        A(n, n, p) - AR model (n - number of signals, p - model order)
        er(n) - variance of innovations
        m - length of simulated time-series
    
      Output:
        data(n, m) - simulated time-series
    '''
    if er == None:
        er = ones(A.shape[0])
    
    ri_.initr()
    r_('library(dse1)')
    
    n = A.shape[0]
    A = concatenate((eye(n).reshape(n,n,1), -A), axis = 2)
    ri_.globalEnv["A"] = ri_.FloatSexpVector(A.ravel())
    ri_.globalEnv["dim"] = ri_.IntSexpVector(A.shape[::-1])
    ri_.globalEnv["er"] = ri_.FloatSexpVector(er)
    ro_.globalEnv["m"] = m
    ro_.globalEnv["n"] = n
    return array(r_('simulate(ARMA(A = array(A, dim), B = diag(1,n)), sd = sqrt(er), sampleT = m)$output')).T

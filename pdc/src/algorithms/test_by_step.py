from numpy import *
import matplotlib.pyplot as pp
from scipy.stats import chi2
import time

import algorithms.pdc_alg as pdc_
import algorithms.assym as ass_
from data_simulation.ar_data import ar_data
from algorithms.pdc_alg import ar_fit

from scipy.stats import cov as cov
from scipy.linalg import cholesky
from scipy.linalg import eigh
from scipy.linalg import inv as inv

def ass_alpha(x, e_var, p, n):
    
    gammai = inv(ass_.bigautocorr(x, p))
    omega = kron(gammai,e_var)
    return omega/n
    
def test_alpha(nm = 100, nd = 100, A = None, er = None, maxp = 2):
    
    if A == None:
        #A = array([[[4,-4],[4,-2]],[[0,0],[0,3]]], dtype=float).reshape(2,2,2)/10
        a21 = 0
        A = array([[[0.2, 0],[-0.4, -0.2],[0.3,0]], 
                   [[a21, 0],[0.8,-0.1],[0.4,0]],
                   [[0,0.5],[-0.1,0.2],[0.4,0.1]]], dtype = float) #Ex artigo daniel
        maxp = 2
    if er == None:
        #er = array([[0.7,0.3],[0.3,2]], dtype = float)
        er = identity(3)
        
    n = A.shape[0]
    varass = empty([nm, n**2*maxp, n**2*maxp])
    al = empty([nm, n**2*maxp])
    time.clock()
    for i in range(nm):
        data = ar_data(A, er, nd)
        Aest, erest = ar_fit(data, maxp = maxp)  
        varass[i] = ass_alpha(data, erest, maxp, nd)
        al[i] = Aest.transpose([2,1,0]).ravel()
      
        if (i%10 == 0):
            print 'nm:', i, 'time:', time.clock()
    varal = cov(al)
    #pp.hist(th10, bins = 50)
    #pp.plot(sorted(pdc10), chi2.ppf(linspace(0.5/montei,1-0.5/montei,montei), 2))
    #pp.show()
    #print h0size/float(montei)
    return varass, varal



def test_dinv(nm = 100):
    #e_var = array([[1, 0.1, 0.2, 0.1], [0, 2, 0.3, 0.2], 
    #         [0.1, 0.3, 1, 0.3], [0.2, 0, 0, 3]])
    #e_var = e_var+e_var.T
    e_var = diag(array([1, 0.1, 0.2, 0.1, 0.3, 2, 0.3, 0.2]))/100
    n = 2
    #e_var = diag(array([0.001,0.0005]))
    #n = 1
    print e_var
    hv = empty([nm, 2*n**2])
    ha = empty([nm, 2*n**2, 2*n**2])
    #mean = ones([2*n**2])
    mean = array([1, 0.3, 0.2, 2, 0.8, -0.3, 0.1, 1])
    #mean = [1,1]
    for i in range(nm):
        
        a = random.multivariate_normal(mean, e_var)
        al = a[0:n**2].reshape(n,n).T + 1j*a[n**2:2*n**2].reshape(n,n).T
        h = inv(al)
        hvt = ass_.vec(h)
        hv[i] = ass_.cat(hvt.real, hvt.imag, 0).reshape(-1)
        #print ass_.fdh_da(mat(al), 2).shape
        dh_da = ass_.fdh_da(mat(al), n)
        #return dh_da
        ha[i] = dot(dot(dh_da, e_var), dh_da.T)
        #print a
        #print al
        #print h
        #print hvt
        #print hv
        #print dh_da
        #print ha[i]
    vara = cov(hv)
    return vara, ha
        
        
        
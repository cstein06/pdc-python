from numpy import *
import matplotlib.pyplot as pp
from scipy.stats import chi2
import time
from scipy import randn

import algorithms.pdc_alg as pdc_
import algorithms.assym as ass_
from data_simulation.ar_data import ar_data
import algorithms.ar_fit as ar_fit

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
        #A = array([[[4],[4]],[[0],[3]]], dtype=float).reshape(2,2,2)/10
        #maxp=1
        A = array([[[4,-4],[4,-2]],[[0,0],[0,3]]], dtype=float).reshape(2,2,2)/10
        maxp=2
        #a21 = 0
        #A = array([[[0.2, 0],[-0.4, -0.2],[0.3,0]], 
        #           [[a21, 0],[0.8,-0.1],[0.4,0]],
        #           [[0,0.5],[-0.1,0.2],[0.4,0.1]]], dtype = float) #Ex artigo daniel
        #maxp = 2
    if er == None:
        er = array([[0.7,0.3],[0.3,2]], dtype = float)
        #er = identity(3)
        
    n = A.shape[0]
    varass = empty([nm, n**2*maxp, n**2*maxp])
    al = empty([nm, n**2*maxp])
    time.clock()
    for i in range(nm):
        data = ar_data(A, er, nd)
        Aest, erest = ar_fit.nstrand(data, maxp = maxp)  
        varass[i] = ass_alpha(data, erest, maxp, nd)
        al[i] = Aest.transpose([2,1,0]).ravel()
      
        if (i%10 == 0):
            print 'nm:', i, 'time:', time.clock()
    varal = cov(al)
    return varass, varal

def debig_de(e_var):
    n = e_var.shape[0]
    
    return dot(kron(ass_.TT(2*n, n), ass_.I(n*2*n)),
           kron(ass_.I(n), kron(ass_.vec(ass_.I(2*n)), ass_.I(n))))
    #p = 2
    #return dot(kron(ass_.TT(n), ass_.I(p*p)),
    #       kron(ass_.I(n), kron(ass_.vec(ass_.I(p)), ass_.I(n))))
    #return dot(kron(ass_.I(p), kron(ass_.TT(n), ass_.I(n))), \
    #           kron(ass_.vec(ass_.I(p)), ass_.I(n*n))) 
    
def ass_evar(e_var, nd):
    n = e_var.shape[0]
    di = ass_.Dup(n).I
    return 2*di*kron(e_var, e_var)*di.T/nd
    
def test_evar(nm = 100, nd = 100, A = None, er = None, maxp = 2):
    
    if A == None:
        A = array([[[4,-4],[4,-2]],[[0,0],[0,3]]], dtype=float).reshape(2,2,2)/10
        #a21 = 0
        #A = array([[[0.2, 0],[-0.4, -0.2],[0.3,0]], 
        #           [[a21, 0],[0.8,-0.1],[0.4,0]],
        #           [[0,0.5],[-0.1,0.2],[0.4,0.1]]], dtype = float) #Ex artigo daniel
        maxp = 2
    if er == None:
        er = array([[0.7,0.3],[0.3,2]], dtype = float)
        #er = identity(3)
        
    n = A.shape[0]
    vares = empty([nm, (n*(n+1))/2, (n*(n+1))/2])
    al = empty([nm, (n*(n+1))/2])
    erm = empty([nm, n, n])
    time.clock()
    for i in range(nm):
        data = ar_data(A, er, nd)
        Aest, erest = ar_fit.nstrand(data, maxp = maxp)  
        vares[i] = ass_evar(erest, nd)
        al[i] = ass_.vech(erest)
        erm[i] = erest
      
        if (i%10 == 0):
            print 'nm:', i, 'time:', time.clock()
    varel = cov(al)
    return vares, varel, erm

def test_eaind(nm = 100, nd = 100, A = None, er = None, maxp = 2):
    
    if A == None:
        A = array([[[4,-4],[4,-2]],[[0,0],[0,3]]], dtype=float).reshape(2,2,2)/10
        #a21 = 0
        #A = array([[[0.2, 0],[-0.4, -0.2],[0.3,0]], 
        #           [[a21, 0],[0.8,-0.1],[0.4,0]],
        #           [[0,0.5],[-0.1,0.2],[0.4,0.1]]], dtype = float) #Ex artigo daniel
        maxp = 2
    if er == None:
        er = array([[0.7,0.3],[0.3,2]], dtype = float)
        #er = identity(3)
        
    n = A.shape[0]
    vares = empty([nm, (n*(n+1))/2, (n*(n+1))/2])
    ale = empty([nm, (n*(n+1))/2])
    ala = empty([nm, n**2*maxp])
    erm = empty([nm, n, n])
    time.clock()
    for i in range(nm):
        data = ar_data(A, er, nd)
        Aest, erest = ar_fit.nstrand(data, maxp = maxp)  
        vares[i] = ass_evar(erest, nd)
        ala[i] = Aest.transpose([2,1,0]).ravel()
        #print ala[i]
        ale[i] = ass_.vech(erest)
        #print ale[i]
        erm[i] = erest
        if (i%10 == 0):
            print 'nm:', i, 'time:', time.clock()
    varel = cov(ale)
    big = ass_.cat(ala, ale, 1)
    varbig = cov(big)
    return vares, varel, erm, varbig 

def test_dinv(nm = 100):
    #e_var = diag(array([1, 0.1, 0.2, 0.1, 0.3, 2, 0.3, 0.2]))/100 + randn(8,8)/10000
   # e_var = diag(randn(2*3**2))/100 + randn(2*3**2,2*3**2)/10000
    #e_var = array([[1, 0.1, 0.2, 0.1], [0, 2, 0.3, 0.2], 
    #         [0.1, 0.3, 1, 0.3], [0.2, 0, 0, 3]])
    #e_var = e_var+e_var.T
    #e_var = diag(array([1, 0.1, 0.2, 0.1, 0.3, 2, 0.3, 0.2]))/100
    e_var = diag(array([1, 0.3, 0.2, 2, 0.8, -0.3, 0.1, 1, 2, 0.8, -0.3, 0.1, 1, 2, 0.8, -0.3, 0.1, 1]))/1000
    #e_var = identity(18)/1000
    n = 3
    #e_var = diag(array([0.001,0.0005]))
    #n = 1
    #print e_var
    hv = empty([nm, 2*n**2])
    ha = empty([nm, 2*n**2, 2*n**2])
    #mea = ones([2*n**2])
    #mea = array([1, 0.3, 0.2, 2, 0.8, -0.3, 0.1, 1])
    #mea = array([1, 0.3, 0.2, 2, 0.8, -0.3, 0.1, 1, 2, 0.8, -0.3, 0.1, 1, 2, 0.8, -0.3, 0.1, 1])
    mea = array(ass_.cat(ass_.vec(identity(3)), ass_.vec(zeros([3,3])), 0)).reshape(-1)
    #print mea
    #mea = [1,1]
    for i in range(nm):
        
        a = random.multivariate_normal(mea, e_var)
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
        if (i%1000 == 0):
            print 'nm:', i, 'time:', time.clock()
    vara = cov(hv)
    return vara, mean(ha, axis = 0)
        
        
        
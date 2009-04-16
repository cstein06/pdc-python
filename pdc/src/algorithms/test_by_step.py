from numpy import *
import matplotlib.pyplot as pp
from scipy.stats import chi2
import time
from scipy import randn

import algorithms.pdc_alg as pdc_
import algorithms.asymp as ass_
from data_simulation.ar_data import ar_data
import algorithms.ar_fit as ar_fit

from scipy.stats import cov as cov
from scipy.linalg import cholesky
from scipy.linalg import eigh
from scipy.linalg import inv as inv
from scipy.stats import chi2

def compare_matlab_pdc_one(nd = 100, nf = 5, metric = 'euc'):
    ''' Compara resultado do pdc e assym pdc com o matlab '''
    
    A = array([[4,3],[0,3]], dtype=float).reshape(2,2,1)/10
    #A = array([[[4,-4],[3,3]],[[0,0],[0,3]]], dtype=float).reshape(2,2,2)/20
    er = array([[0.7,0.1],[0.1,2]], dtype = float)
    IP = A.shape[2]
    u1 = sin(linspace(0,10,nd)).reshape(1,-1)
    u2 = sin(linspace(0,13,nd)).reshape(1,-1)
    u = concatenate((u1, u2), axis = 0)
    
    pdc = pdc_.pdc_alg(A, er, nf, metric = metric)
    #print abs(pdc)**2
    
    Af = pdc_.A_to_f(A, nf = nf)
    th, ic1, ic2 = ass_.assym_pdc(u, Af, er, IP, alpha = 0.05, metric = metric)
    print 'pdc', abs(pdc)**2
    print 'th', th
    print 'ic1', ic1
    print 'ic2', ic2

def test_pdc():
    A = array([[[4,-4],[3,3]],[[0,0],[0,3]]], dtype=float).reshape(2,2,2)/20
    er = array([[0.7,0],[0,2]], dtype = float)
    pdc = pdc_.pdc_one_alg(A, er, nf = 5)
    #pdc_plot(pdc)
    #print abs(pdc)**2
    return pdc
    
def test_alpha_var(montei = 100, nd = 100, A = None, er = None, maxp = 2):
    # not ready.
    if A == None:
        A = array([[[4,-4],[3,3]],[[0,0],[0,3]]], dtype=float).reshape(2,2,2)/20
    if er == None:
        er = array([[0.7,0],[0,2]], dtype = float)
    n = A.shape[0]
    p = A.shape[2]
    alpha = empty(montei, n, n, p)
    for i in arange(montei):
        data = ar_data(A, er, nd)
        alpha[i], erest = nstrand(data, maxp = maxp)
        
def test_assym_pdc_semr():
    Aest = array([[4,3],[0,3]], dtype=float).reshape(2,2,1)/10
    erest = array([[0.7,0],[0,2]], dtype = float)
    data = ar_data(Aest, erest, 1000)
    maxp = 1
    nf = 5
    Af = pdc_.A_to_f(Aest, nf = nf)
    pdc = pdc_.pdc_alg(Aest, erest)
    th, ic1, ic2 = ass_.assym_pdc(data, Af, erest, maxp)
    print 'pdc', pdc
    print 'th', th
    print 'ic1', ic1
    print 'ic2', ic2


def testR():
    n = 2
    tR1 = kron(I(n),kron(array([[1],[0]]), I(n)))
    tR2 = kron(I(n),kron(array([[0],[1]]), I(n)))
    tR = mat(cat(tR1, tR2, 1))
    tS = mat(tR.T)
    
    a = arange(8*8).reshape(8,8)
    print a
    print tR
    print dot(a, tR)
    print dot(dot(tR, a), tR)

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
        

def test_pc():
    Aest = array([[4,2],[0,3]], dtype=float).reshape(2,2,1)/10
    #erest = array([[0.7,0],[0,2]], dtype = float)
    erest = identity(2)
    data = ar_data(Aest, erest, 1000)
    maxp = 1
    nf = 2
    Af = pdc_.A_to_f(Aest, nf = nf)
    pc = pdc_.pc(Aest, erest, nf = nf)
    th, ic1, ic2, v1, v2 = ass_.assym_pc(data, Af, erest, maxp)
    print 'pc', abs(pc)**2
    print 'th', th
    print 'ic1', ic1
    print 'ic2', ic2
    print 'pca', (ic1+ic2)/2
    

def test_coh():
    Aest = array([[4,2],[0,3]], dtype=float).reshape(2,2,1)/10
    #erest = array([[0.7,0],[0,2]], dtype = float)
    erest = identity(2)
    data = ar_data(Aest, erest, 1000)
    maxp = 1
    nf = 2
    Af = pdc_.A_to_f(Aest, nf = nf)
    coh = pdc_.coh_alg(Aest, erest, nf = nf)
    th, ic1, ic2, v1, v2 = ass_.assym_coh(data, Af, erest, maxp)
    print 'pc', abs(coh)**2
    print 'th', th
    print 'ic1', ic1
    print 'ic2', ic2
    print 'pca', (ic1+ic2)/2
    
def test_ss():
    #Aest = array([[4,2],[0,3]], dtype=float).reshape(2,2,1)/10
    Aest = array([[[0.2, 0.4],[0.3, 0.2],[0.3,-0.2]], 
                   [[0, 0],[0.8,-0.1],[0,0]],
                   [[0, 0],[0,0],[0.4,0.1]]], dtype = float) 
    #erest = array([[0.7,0],[0,2]], dtype = float)
    erest = identity(3)
    data = ar_data(Aest, erest, 1000)
    maxp = 1
    nf = 2
    Af = pdc_.A_to_f(Aest, nf = nf)
    ss = pdc_.ss(Aest, erest, nf = nf)
    th, ic1, ic2, v1, v2 = ass_.assym_ss(data, Af, erest, maxp)
    print 'ss', abs(ss)**2
    print 'th', th
    print 'ic1', ic1
    print 'ic2', ic2
    print 'ssa', (ic1+ic2)/2
        
def test_pdc_normalizations():
    A = array([[4,3],[0,3]], dtype=float).reshape(2,2,1)/10
    er = array([[0.7,0.1],[0.1,2]], dtype = float)
    #er = identity(2)
    nf = 5
    #print pdc_.pdc_one_alg(A, er, nf = nf)
    #print pdc_.pdc_gen_alg(A, er, nf = nf)
    #print pdc_.pdc_diag_alg(A, er, nf = nf)
    print pdc_.pdc_alg(A, er, metric = 'euc', nf = nf)
    print pdc_.pdc_alg(A, er, metric = 'diag', nf = nf)
    print pdc_.pdc_alg(A, er, metric = 'gen', nf = nf)

def test_patnaik(d, mci = 1000):
    ''' Gera varias amostras sob a distribuicao patnaik e sobre a 
    distribuicao exata, a fim de comparar as distribuicoes '''
    
    patdf = sum(d)**2/sum(d**2)
    patden = sum(d)/sum(d**2)
    if (d.size != 2):
        print 'expecting d.size = 2'
        print d
    ds = d.size
    pat = empty(mci)
    rpat = zeros(mci)
    for i in arange(mci):
        pat[i] = chi2.rvs(patdf)/patden
        for j in arange(ds):
            rpat[i] = rpat[i] + chi2.rvs(1)*d[j]
    return pat, rpat
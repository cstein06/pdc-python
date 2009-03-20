from numpy import *
import matplotlib.pyplot as pp
from scipy.stats import chi2
import time

import algorithms.pdc_alg as pdc_
import algorithms.assym as ass_
from data_simulation.ar_data import ar_data
from algorithms.pdc_alg import ar_fit


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
        alpha[i], erest = ar_fit(data, maxp = maxp)
        
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
    
def test_assym_pdc_old(nm = 100, nd = 100, A = None, er = None, 
                      maxp = 2, nf = 20, alpha = 0.05, metric = 'euc'):
    ''' Testa se test H0 e H1 do pdc esta de acordo com o alpha.
    
    o th retorna threshold de h0, ic1 e ic2 o intervalo de confianca, 
    h0size o tamanho do teste sob h0; e h11 e h12 o tamanho sob h1. '''
    
    if A == None:
        #A = array([[[4,-4],[4,-2]],[[0,0],[0,3]]], dtype=float).reshape(2,2,2)/10
        a21 = 0
        A = array([[[0.2, 0],[-0.4, -0.2],[0.3,0]], 
                   [[a21, 0],[0.8,-0.1],[0.4,0]],
                   [[0,0.5],[-0.1,0.2],[0.4,0.1]]], dtype = float) #Ex artigo daniel
    if er == None:
        #er = array([[0.7,0.3],[0.3,2]], dtype = float)
        #er = identity(3)
        er = array([[1,0.3,0.2],[0.3,2,0.6], [0.2, 0.6, 2]], dtype = float)
    n = A.shape[0]
    pdc = empty([nm, n, n, nf])
    th  = empty([nm, n, n, nf])
    ic1 = empty([nm, n, n, nf])
    ic2 = empty([nm, n, n, nf])
    varass = empty([nm, n, n, nf])
    varass2 = empty([nm, n, n, nf])
    time.clock()
    for i in range(nm):
        data = ar_data(A, er, nd)
        Aest, erest = ar_fit(data, maxp = maxp)
        pdc[i] = abs(pdc_.pdc_alg(Aest, erest, nf = nf, metric = metric))**2
        th[i], ic1[i], ic2[i], varass[i], varass2[i] = ass_.assym_pdc(data, pdc_.A_to_f(Aest, nf = nf), erest, 
                                               maxp, alpha = alpha, metric = metric)
        if (i%10 == 0):
            print 'nm:', i, 'time:', time.clock()

    h0size = sum((pdc < th), axis = 0)/float(nm)
    pdcr = abs(pdc_.pdc_alg(A, er, nf = nf, metric = metric))**2
    h11 = sum((pdcr < ic1), axis = 0)/float(nm)
    h12 = sum((pdcr > ic2), axis = 0)/float(nm)
    #pp.hist(th10, bins = 50)
    #pp.plot(sorted(pdc10), chi2.ppf(linspace(0.5/montei,1-0.5/montei,montei), 2))
    #pp.show()
    #print h0size/float(montei)
    return pdc, th, ic1, ic2, h0size, h11, h12, pdcr, varass, varass2

def test_assym_pdc(nm = 100, nth = 50, nd = 100, A = None, er = None, 
                      maxp = 2, nf = 20, alpha = 0.05, metric = 'euc'):
    ''' Testa se test H0 e H1 do pdc esta de acordo com o alpha.
    
    o th retorna threshold de h0, ic1 e ic2 o intervalo de confianca, 
    h0size o tamanho do teste sob h0; e h11 e h12 o tamanho sob h1. '''
    
    if A == None:
        #A = array([[[4,-4],[4,-2]],[[0,0],[0,3]]], dtype=float).reshape(2,2,2)/10
        a21 = 0
        A = array([[[0.2, 0],[-0.4, -0.2],[0.3,0]], 
                   [[a21, 0],[0.8,-0.1],[0.4,0]],
                   [[0,0.5],[-0.1,0.2],[0.4,0.1]]], dtype = float) #Ex artigo daniel
    if er == None:
        #er = array([[0.7,0.3],[0.3,2]], dtype = float)
        #er = identity(3)
        er = array([[1,0.3,0.2],[0.3,2,0.6], [0.2, 0.6, 2]], dtype = float)
    n = A.shape[0]
    pdc = empty([nm, n, n, nf])
    th  = empty([nm, n, n, nf])
    ic1 = empty([nm, n, n, nf])
    ic2 = empty([nm, n, n, nf])
    varass = empty([nm, n, n, nf])
    varass2 = empty([nm, n, n, nf])
    time.clock()
    for i in range(nm):
        data = ar_data(A, er, nd)
        Aest, erest = ar_fit(data, maxp = maxp)
        pdc[i] = abs(pdc_.pdc_alg(Aest, erest, nf = nf, metric = metric))**2
        th[i], ic1[i], ic2[i], varass[i], varass2[i] = ass_.assym_pdc(data, pdc_.A_to_f(Aest, nf = nf), erest, 
                                               maxp, alpha = alpha, metric = metric)
        if (i%10 == 0):
            print 'nm:', i, 'time:', time.clock()

    h0size = sum((pdc < th), axis = 0)/float(nm)
    pdcr = abs(pdc_.pdc_alg(A, er, nf = nf, metric = metric))**2
    h11 = sum((pdcr < ic1), axis = 0)/float(nm)
    h12 = sum((pdcr > ic2), axis = 0)/float(nm)
    #pp.hist(th10, bins = 50)
    #pp.plot(sorted(pdc10), chi2.ppf(linspace(0.5/montei,1-0.5/montei,montei), 2))
    #pp.show()
    #print h0size/float(montei)
    return pdc, th, ic1, ic2, h0size, h11, h12, pdcr, varass, varass2


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

def test_assym_dtf(nm = 100, nd = 100, A = None, er = None, 
                  maxp = 2, nf = 20, alpha = 0.05, metric = 'euc'):
    ''' Testa se test H0 e H1 do pdc esta de acordo com o alpha.
    
    o th retorna threshold de h0, ic1 e ic2 o intervalo de confianca, 
    h0size o tamanho do teste sob h0; e h11 e h12 o tamanho sob h1. '''
    
    if A == None:
        #A = array([[[4,-4],[6,-3]],[[0,0],[0,3]]], dtype=float).reshape(2,2,2)/10
        #a21 = 0
        #A = array([[[0.2, 0],[-0.4, -0.2],[0.3,0]], 
        #           [[a21, 0],[0.8,-0.1],[0.4,0]],
         #          [[0,0.5],[-0.1,0.2],[0.4,0.1]]], dtype = float) #Ex artigo daniel
        
        A = array([[[0.4, -0.6],[-0.3, -0.2],[0.4,-0.2]], 
                   [[-0.2,0.4],[0.8,-0.1],[0.4,0.3]],
                   [[0,0],[0,0],[0.4,0.1]]], dtype = float) #Ex dtf = 0
    if er == None:
        #er = array([[1,0],[0,1]], dtype = float)
        er = identity(3)
    n = A.shape[0]
    dtf = empty([nm, n, n, nf])
    th  = empty([nm, n, n, nf])
    ic1 = empty([nm, n, n, nf])
    ic2 = empty([nm, n, n, nf])
    varass = empty([nm, n, n, nf])
    varass2 = empty([nm, n, n, nf])
    time.clock()
    for i in range(nm):
        data = ar_data(A, er, nd)
        Aest, erest = ar_fit(data, maxp = maxp)
        dtf[i] = abs(pdc_.dtf_one_alg(Aest, erest, nf = nf))**2
        th[i], ic1[i], ic2[i], varass[i], varass2[i] = ass_.assym_dtf_one(data, pdc_.A_to_f(Aest, nf = nf), erest, 
                                               maxp, alpha = alpha)
        if (i%10 == 0):
            print 'nm:', i, 'time:', time.clock()

    h0size = sum((dtf < th), axis = 0)/float(nm)
    dtfr = abs(pdc_.dtf_one_alg(A, er, nf = nf))**2
    h11 = sum((dtfr < ic1), axis = 0)/float(nm)
    h12 = sum((dtfr > ic2), axis = 0)/float(nm)
    #pp.hist(th10, bins = 50)
    #pp.plot(sorted(pdc10), chi2.ppf(linspace(0.5/montei,1-0.5/montei,montei), 2))
    #pp.show()
    #print h0size/float(montei)
    return dtf, th, ic1, ic2, h0size, h11, h12, dtfr, varass, varass2
    
        
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
    pat = empty(mci)
    rpat = empty(mci)
    for i in arange(mci):
        pat[i] = chi2.rvs(patdf)/patden
        rpat[i] = chi2.rvs(1)*d[0] + chi2.rvs(1)*d[1]
    return pat, rpat

if __name__ == "__main__":
    #test_assym_pdc_semr();
    #test_assym_pdc_th_bootstrap()
    test_pdc_normalizations()
    

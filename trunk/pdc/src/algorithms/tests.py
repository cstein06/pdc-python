from numpy import *
import matplotlib.pyplot as pp
from scipy.stats import chi2

import algorithms.pdc_alg as pdc_
import algorithms.assym as ass_
from data_simulation.ar_data import ar_data
from algorithms.pdc_alg import ar_fit

from numpy import *

def compare_matlab_pdc_one():
    
    #A = array([[4,3],[0,3]], dtype=float).reshape(2,2,1)/10
    A = array([[[4,-4],[3,3]],[[0,0],[0,3]]], dtype=float).reshape(2,2,2)/20
    er = array([[0.7,0],[0,2]], dtype = float)
    IP = A.shape[2]
    nd = 100
    u1 = sin(linspace(0,10,nd)).reshape(1,-1)
    u2 = sin(linspace(0,13,nd)).reshape(1,-1)
    u = concatenate((u1, u2), axis = 0)
    nf = 5
    
    pdc = pdc_.pdc_one_alg(A, er, nf)
    #print abs(pdc)**2
    
    Af = pdc_.A_to_f(A, nf = nf)
    th, ic1, ic2 = ass_.assym_pdc_one(u, Af, er, IP, alpha = 0.05)
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
    
def test_assym_pdc_th(montei = 100, nd = 100, A = None, er = None, maxp = 2, nf = 20, alpha = 0.05):
    if A == None:
        A = array([[[4,-4],[3,3]],[[0,0],[0,3]]], dtype=float).reshape(2,2,2)/20
    if er == None:
        er = array([[0.7,0],[0,2]], dtype = float)
    n = A.shape[0]
    pdc = empty([montei, n, n, nf])
    th  = empty([montei, n, n, nf])
    ic1 = empty([montei, n, n, nf])
    ic2 = empty([montei, n, n, nf])
    for i in range(montei):
        data = ar_data(A, er, nd)
        Aest, erest = ar_fit(data, maxp = maxp)
        pdc[i] = abs(pdc_.pdc_one_alg(Aest, erest, nf = nf))**2
        th[i], ic1[i], ic2[i] = ass_.assym_pdc_one(data, pdc_.A_to_f(Aest, nf = nf), erest, maxp, alpha = alpha)

    h0size = sum((pdc < th), axis = 0)/float(montei)
    pdcr = abs(pdc_.pdc_one_alg(A, er, nf = nf))**2
    #pp.hist(th10, bins = 50)
    #pp.plot(sorted(pdc10), chi2.ppf(linspace(0.5/montei,1-0.5/montei,montei), 2))
    #pp.show()
    #print h0size/float(montei)
    return pdc, th, ic1, ic2, h0size, pdcr
    
        
def test_assym_pdc_ic():
    booti = 1000
    A = array([[4,3],[0,3]], dtype=float).reshape(2,2,1)/10
    er = array([[0.7,0],[0,2]], dtype = float)
    pdcr = pdc_.pdc_one_alg(A, nf = nf)
    maxp = 1
    nf = 8
    icsize = empty(64)
    for i in range(booti):
        data = ar_data(A, er, 1000)
        Aest, erest = ar_fit(A, maxp = maxp)
        pdc = pdc_.pdc_one_alg(Aest, nf = nf)
        th, ic1, ic2 = ass_.assym_pdc(data, A_to_f(Aest, nf = nf), erest, maxp)
        for j in range(nf):
            if abs(pdcr[0,1,j])**2 < ic2[j] and abs(pdcr[0,1,j])**2 > ic1[j] :
                icsize[j] = icsize[j] + 1
    
    print icsize/float(booti)
        
def test_pdc_normalizations():
    A = array([[4,3],[0,3]], dtype=float).reshape(2,2,1)/10
    #er = array([[0.7,0],[0,2]], dtype = float)
    er = identity(2)
    nf = 5
    print pdc_.pdc_one_alg(A, er, nf = nf)
    print pdc_.pdc_gen_alg(A, er, nf = nf)
    print pdc_.pdc_diag_alg(A, er, nf = nf)

def test_patnaik(d, mci = 1000):
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
    test_assym_pdc_th_bootstrap()
    #test_pdc_normalizations()
    

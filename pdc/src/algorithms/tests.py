from numpy import *
import matplotlib.pyplot as pp
from scipy.stats import chi2
import time

import algorithms.pdc_alg as pdc_
import algorithms.assym as ass_
from data_simulation.ar_data import ar_data
from algorithms.ar_fit import nstrand


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
    
def test_assym_pdc_old(nm = 100, nd = 100, A = None, er = None, 
                      maxp = 2, nf = 20, alpha = 0.05, metric = 'euc'):
    ''' Testa se test H0 e H1 do pdc esta de acordo com o alpha.
    
    o th retorna threshold de h0, ic1 e ic2 o intervalo de confianca, 
    h0size o tamanho do teste sob h0; e h11 e h12 o tamanho sob h1. '''
    
    if A == None:
        A = array([[[4,-4],[4,-2]],[[0,0],[0,3]]], dtype=float).reshape(2,2,2)/10
        #a21 = 0
        ##A = array([[[0.2, 0],[-0.4, -0.2],[0.3,0]], 
        #           [[a21, 0],[0.8,-0.1],[0.4,0]],
        #           [[0,0.5],[-0.1,0.2],[0.4,0.1]]], dtype = float) #Ex artigo daniel
    if er == None:
        er = array([[0.7,0.3],[0.3,2]], dtype = float)
        #er = identity(3)
        #er = array([[1,0.3,0.2],[0.3,2,0.6], [0.2, 0.6, 2]], dtype = float)
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
        Aest, erest = nstrand(data, maxp = maxp)
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

def test_assym_pdc(np = 100, nth = 50, nd = 100, A = None, er = None, 
                      maxp = 2, nf = 20, alpha = 0.05, metric = 'euc'):
    ''' Testa se test H0 e H1 do pdc esta de acordo com o alpha.
    
    o th retorna threshold de h0, ic1 e ic2 o intervalo de confianca, 
    h0size o tamanho do teste sob h0; e h11 e h12 o tamanho sob h1. '''
    
    if A == None:
        A = array([[[4,-4],[4,-2]],[[0,0],[0,3]]], dtype=float).reshape(2,2,2)/10
        #a21 = 0
        #A = array([[[0.2, 0],[-0.4, -0.2],[0.3,0]], 
        #           [[a21, 0],[0.8,-0.1],[0.4,0]],
        #           [[0,0.5],[-0.1,0.2],[0.4,0.1]]], dtype = float) #Ex artigo daniel
    if er == None:
        er = array([[0.7,0.3],[0.3,2]], dtype = float)
        #er = identity(3)
        #er = array([[1,0.3,0.2],[0.3,2,0.6], [0.2, 0.6, 2]], dtype = float)
    n = A.shape[0]
    pdc = empty([np+nth, n, n, nf])
    th  = empty([nth, n, n, nf])
    ic1 = empty([nth, n, n, nf])
    ic2 = empty([nth, n, n, nf])
    varass = empty([nth, n, n, nf])
    varass2 = empty([nth, n, n, nf])
    time.clock()
    for i in range(nth):
        data = ar_data(A, er, nd)
        Aest, erest = nstrand(data, maxp = maxp)
        pdc[i] = abs(pdc_.pdc_alg(Aest, erest, nf = nf, metric = metric))**2
        th[i], ic1[i], ic2[i], varass[i], varass2[i] = ass_.assym_pdc(data, pdc_.A_to_f(Aest, nf = nf), erest, 
                                               maxp, alpha = alpha, metric = metric)
        if (i%10 == 0):
            print 'nth:', i, 'time:', time.clock()

    for i in range(np):        
        if (i%10 == 0):
            print 'nm1:', i, 'time:', time.clock()
        data = ar_data(A, er, nd)        
        if (i%10 == 0):
            print 'nm2:', i, 'time:', time.clock()
        Aest, erest = nstrand(data, maxp = maxp)
        if (i%10 == 0):
            print 'nm3:', i, 'time:', time.clock()
        pdc[i+nth] = abs(pdc_.pdc_alg(Aest, erest, nf = nf, metric = metric))**2
        if (i%10 == 0):
            print 'np4:', i, 'time:', time.clock()
        

    h0size = sum((pdc < mean(th, axis=0)), axis = 0)/float(np+nth)
    pdcr = abs(pdc_.pdc_alg(A, er, nf = nf, metric = metric))**2
    h11 = sum((pdcr < ic1), axis = 0)/float(nth)
    h12 = sum((pdcr > ic2), axis = 0)/float(nth)
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
        #er = array([[0.7,0.3],[0.3,2]], dtype = float)
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
        #print 'nm:', i, 'time:', time.clock()
        data = ar_data(A, er, nd)
        #print 'data', 'time:', time.clock()
        Aest, erest = nstrand(data, maxp = maxp)
        #print 'nstrand', 'time:', time.clock()
        dtf[i] = abs(pdc_.dtf_one_alg(Aest, erest, nf = nf))**2
        
        th[i], ic1[i], ic2[i], varass[i], varass2[i] = ass_.assym_dtf_one(data, pdc_.A_to_f(Aest, nf = nf), erest, 
                                               maxp, alpha = alpha)
        #print 'ass', 'time:', time.clock()
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

def test_assym_pc(nm = 100, nd = 100, A = None, er = None, 
                  maxp = 2, nf = 20, alpha = 0.05):
    ''' Testa se test H0 e H1 do pdc esta de acordo com o alpha.
    
    o th retorna threshold de h0, ic1 e ic2 o intervalo de confianca, 
    h0size o tamanho do teste sob h0; e h11 e h12 o tamanho sob h1. '''
    
    if A == None:
        #A = array([[[4,-4],[6,-3]],[[0,0],[0,3]]], dtype=float).reshape(2,2,2)/10
        A = array([[[0.2, 0],[0, 0],[0.3,-0.2]], 
                   [[0, 0],[0.8,-0.1],[0.4,-0.1]],
                   [[0, 0],[-0.1,0.2],[0.4,0.1]]], dtype = float) 
        
        #A = array([[[0.4, -0.6],[-0.3, -0.2],[0.4,-0.2]], 
        #           [[-0.2,0.4],[0.8,-0.1],[0.4,0.3]],
        #           [[0,0],[0,0],[0.4,0.1]]], dtype = float) #Ex dtf = 0
    if er == None:
        #er = array([[1,0],[0,1]], dtype = float)
        er = identity(3)
        #er = array([[0.7,0.3, 0], [0.3, 1.2, 0.4], [0, 0.4, 2]], dtype = float)
        
    n = A.shape[0]
    pc = empty([nm, n, n, nf])
    th  = empty([nm, n, n, nf])
    ic1 = empty([nm, n, n, nf])
    ic2 = empty([nm, n, n, nf])
    varass = empty([nm, n, n, nf])
    varass2 = empty([nm, n, n, nf])
    time.clock()
    for i in range(nm):
        #print time.clock()
        data = ar_data(A, er, nd)
        Aest, erest = nstrand(data, maxp = maxp)
        pc[i] = abs(pdc_.pc(Aest, erest, nf = nf))**2
        #print time.clock()
        th[i], ic1[i], ic2[i], varass[i], varass2[i] = ass_.assym_pc(data, pdc_.A_to_f(Aest, nf = nf), erest, 
                                               maxp, alpha = alpha)
        #print time.clock()
        if (i%10 == 0):
            print 'nm:', i, 'time:', time.clock()

    h0size = sum((pc < th), axis = 0)/float(nm)
    pcr = abs(pdc_.pc(A, er, nf = nf))**2
    h11 = sum((pcr < ic1), axis = 0)/float(nm)
    h12 = sum((pcr > ic2), axis = 0)/float(nm)
    #pp.hist(th10, bins = 50)
    #pp.plot(sorted(pdc10), chi2.ppf(linspace(0.5/montei,1-0.5/montei,montei), 2))
    #pp.show()
    #print h0size/float(montei)
    return pc, th, ic1, ic2, h0size, h11, h12, pcr, varass, varass2

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

def test_assym_coh(nm = 100, nd = 100, A = None, er = None, 
                  maxp = 2, nf = 20, alpha = 0.05):
    ''' Testa se test H0 e H1 do pdc esta de acordo com o alpha.
    
    o th retorna threshold de h0, ic1 e ic2 o intervalo de confianca, 
    h0size o tamanho do teste sob h0; e h11 e h12 o tamanho sob h1. '''
    
    if A == None:
        #A = array([[[4,-4],[0,0]],[[0,0],[0,3]]], dtype=float).reshape(2,2,2)/10
        A = array([[[0.2, 0.4],[0.3, 0.2],[0.3,-0.2]], 
                   [[0, 0],[0.8,-0.1],[0,0]],
                   [[0, 0],[0,0],[0.4,0.1]]], dtype = float) 
        
        #A = array([[[0.4, -0.6],[-0.3, -0.2],[0.4,-0.2]], 
        #           [[-0.2,0.4],[0.8,-0.1],[0.4,0.3]],
        #           [[0,0],[0,0],[0.4,0.1]]], dtype = float) #Ex dtf = 0
    if er == None:
        #er = array([[1,0],[0,1]], dtype = float)
        er = identity(3)
        #er = array([[0.7,0.3, 0], [0.3, 1.2, 0.4], [0, 0.4, 2]], dtype = float)
        
    n = A.shape[0]
    coh = empty([nm, n, n, nf])
    th  = empty([nm, n, n, nf])
    ic1 = empty([nm, n, n, nf])
    ic2 = empty([nm, n, n, nf])
    varass = empty([nm, n, n, nf])
    varass2 = empty([nm, n, n, nf])
    time.clock()
    for i in range(nm):
        data = ar_data(A, er, nd)
        Aest, erest = nstrand(data, maxp = maxp)
        coh[i] = abs(pdc_.coh_alg(Aest, erest, nf = nf))**2
        th[i], ic1[i], ic2[i], varass[i], varass2[i] = ass_.assym_coh(data, pdc_.A_to_f(Aest, nf = nf), erest, 
                                               maxp, alpha = alpha)
        if (i%10 == 0):
            print 'nm:', i, 'time:', time.clock()

    h0size = sum((coh < th), axis = 0)/float(nm)
    cohr = abs(pdc_.coh_alg(A, er, nf = nf))**2
    h11 = sum((cohr < ic1), axis = 0)/float(nm)
    h12 = sum((cohr > ic2), axis = 0)/float(nm)
    #pp.hist(th10, bins = 50)
    #pp.plot(sorted(pdc10), chi2.ppf(linspace(0.5/montei,1-0.5/montei,montei), 2))
    #pp.show()
    #print h0size/float(montei)
    return coh, th, ic1, ic2, h0size, h11, h12, cohr, varass, varass2

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

def test_assym_ss(nm = 100, nd = 100, A = None, er = None, 
                  maxp = 2, nf = 20, alpha = 0.05):
    ''' Testa se test H0 e H1 do pdc esta de acordo com o alpha.
    
    o th retorna threshold de h0, ic1 e ic2 o intervalo de confianca, 
    h0size o tamanho do teste sob h0; e h11 e h12 o tamanho sob h1. '''
    
    if A == None:
        #A = array([[[4,-4],[6,-3]],[[0,0],[0,3]]], dtype=float).reshape(2,2,2)/10
        A = array([[[0.2, 0.4],[0.3, 0.2],[0.3,-0.2]], 
                   [[0, 0],[0.8,-0.1],[0,0]],
                   [[0, 0],[0,0],[0.4,0.1]]], dtype = float) 
        
        #A = array([[[0.4, -0.2],[-0.3, -0.2],[0.4,-0.2]], 
        #           [[-0.2,0.2],[0.6,-0.1],[0.2,0.3]],
        #           [[0,0],[0.2,0],[0.3,0.1]]], dtype = float) #Ex dtf = 0
    if er == None:
        #er = array([[1,0],[0,1]], dtype = float)
        er = identity(3)
        #er = array([[0.7,0.3, 0], [0.3, 1.2, 0.4], [0, 0.4, 2]], dtype = float)
        
    n = A.shape[0]
    ss = empty([nm, n, n, nf])
    th  = empty([nm, n, n, nf])
    ic1 = empty([nm, n, n, nf])
    ic2 = empty([nm, n, n, nf])
    varass = empty([nm, n, n, nf])
    varass2 = empty([nm, n, n, nf])
    time.clock()
    for i in range(nm):
        data = ar_data(A, er, nd)
        Aest, erest = nstrand(data, maxp = maxp)
        ss[i] = abs(pdc_.ss(Aest, erest, nf = nf))**2
        th[i], ic1[i], ic2[i], varass[i], varass2[i] = ass_.assym_ss(data, pdc_.A_to_f(Aest, nf = nf), erest, 
                                               maxp, alpha = alpha)
        if (i%10 == 0):
            print 'nm:', i, 'time:', time.clock()

    h0size = sum((ss < th), axis = 0)/float(nm)
    ssr = abs(pdc_.ss(A, er, nf = nf))**2
    h11 = sum((ssr < ic1), axis = 0)/float(nm)
    h12 = sum((ssr > ic2), axis = 0)/float(nm)
    #pp.hist(th10, bins = 50)
    #pp.plot(sorted(pdc10), chi2.ppf(linspace(0.5/montei,1-0.5/montei,montei), 2))
    #pp.show()
    #print h0size/float(montei)
    return ss, th, ic1, ic2, h0size, h11, h12, ssr, varass, varass2


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


def test_assym_coh_test(nm = 100, nd = 100, A = None, er = None, 
                  maxp = 2, nf = 20, alpha = 0.05):
    ''' Testa se test H0 e H1 do pdc esta de acordo com o alpha.
    
    o th retorna threshold de h0, ic1 e ic2 o intervalo de confianca, 
    h0size o tamanho do teste sob h0; e h11 e h12 o tamanho sob h1. '''
    
    if A == None:
        #A = array([[[4,-4],[0,0]],[[0,0],[0,3]]], dtype=float).reshape(2,2,2)/10
        A = array([[[0.4, -0.2],[-0.3, -0.2],[0.4,-0.2]], 
                   [[-0.2,0.2],[0.6,-0.1],[0.2,0.3]],
                   [[0,0],[0.2,0],[0.3,0.1]]], dtype = float) #Ex dtf = 0
        
        #A = array([[[0.4, -0.6],[-0.3, -0.2],[0.4,-0.2]], 
        #           [[-0.2,0.4],[0.8,-0.1],[0.4,0.3]],
        #           [[0,0],[0,0],[0.4,0.1]]], dtype = float) #Ex dtf = 0
    if er == None:
        #er = array([[1,0],[0,1]], dtype = float)
        er = identity(3)
        #er = array([[0.7,0.3, 0], [0.3, 1.2, 0.4], [0, 0.4, 2]], dtype = float)
        
    n = A.shape[0]
    coh = empty([nm, n, n, nf])
    th  = empty([nm, n, n, nf])
    ic1 = empty([nm, n, n, nf])
    ic2 = empty([nm, n, n, nf])
    varass = empty([nm, n, n, nf])
    varass2 = empty([nm, n, n, nf])
    time.clock()
    for i in range(nm):
        data = ar_data(A, er, nd)
        Aest, erest = nstrand(data, maxp = maxp)
        coh[i] = ass_.coh_test(pdc_.A_to_f(Aest, nf = nf), er)
        th[i], ic1[i], ic2[i], varass[i], varass2[i] = ass_.assym_coh_test(data, pdc_.A_to_f(Aest, nf = nf), er, 
                                               maxp, alpha = alpha, At = pdc_.A_to_f(A, nf = nf))
        if (i%10 == 0):
            print 'nm:', i, 'time:', time.clock()

    h0size = sum((coh < th), axis = 0)/float(nm)
    cohr = ass_.coh_test(pdc_.A_to_f(A, nf = nf), er)
    h11 = sum((cohr < ic1), axis = 0)/float(nm)
    h12 = sum((cohr > ic2), axis = 0)/float(nm)
    #pp.hist(th10, bins = 50)
    #pp.plot(sorted(pdc10), chi2.ppf(linspace(0.5/montei,1-0.5/montei,montei), 2))
    #pp.show()
    #print h0size/float(montei)
    return coh, th, ic1, ic2, h0size, h11, h12, cohr, varass, varass2


if __name__ == "__main__":
    #test_assym_pdc_semr();
    #test_assym_pdc_th_bootstrap()
    test_pdc_normalizations()
    

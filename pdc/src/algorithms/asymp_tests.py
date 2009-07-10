#encoding:UTF-8

from numpy import *
import matplotlib.pyplot as pp
from scipy.stats import chi2
import scipy.stats as st
import time

import algorithms.asymp as ass_
import algorithms.pdc_alg as pdc_
from data_simulation.ar_data import ar_data
from data_simulation.ar_data import ar_models
from algorithms.ar_fit import nstrand

def test_asymp(asymp_func, method_func, nm = 100, nd = 100, A = None, er = None, 
               nf = 20, alpha = 0.05, metric = None):
    ''' Testa se test H0 e H1 esta de acordo com o alpha.
    
    o th retorna threshold de h0, ic1 e ic2 o intervalo de confianca, 
    h0size o tamanho do teste sob h0; e h11 e h12 o tamanho sob h1. 
    Input:
        asymp_func -> função para estatistica assintotica. (exs: ass_.asymp_pc, ass_.asymp_pdc, ...)
        method_func -> função correspondente a estatística. (exs: pdc_.pc_alg, pdc_.pdc_alg, ...)
        nm -> número de iterações do Monte Carlo
        nd -> tamanho dos dados gerados
        A -> matriz do VAR
        er -> covariância da inovação
        nf -> quantas frequências calcular
        alpha -> tamanho do intervalo de confiança
        metrica -> no caso do pdc, dizer qual metrica
    '''
    if A == None:
        A = array([[[0.2, 0],[0, 0],[0.3,-0.2]], 
                   [[0, 0],[0.8,-0.1],[0.4,-0.1]],
                   [[0, 0],[0,0],[0.4,0.1]]], dtype = float) 
    if er == None:
        er = identity(3)
    
    n = A.shape[0]
    mes = empty([nm, n, n, nf])
    th  = empty([nm, n, n, nf])
    ic1 = empty([nm, n, n, nf])
    ic2 = empty([nm, n, n, nf])
    maxp = A.shape[2]
    tbegin = time.clock()
    for i in range(nm):
        if (i%10 == 0):
            print 'nm:', i, 'time:', time.clock()-tbegin
        #Generate data from AR
        data = ar_data(A, er, nd)
        #Estimate AR parameters with Nuttall-Strand
        Aest, erest = nstrand(data, maxp = maxp)
        #Calculate the connectivity and statistics
        if (metric == None):
            mes[i], th[i], ic1[i], ic2[i] = asymp_func(data, Aest, nf, erest, 
                                                   maxp, alpha = alpha)
        else:
            mes[i], th[i], ic1[i], ic2[i] = asymp_func(data, Aest, nf, erest, 
                                                   maxp, alpha = alpha, metric = metric)

    h0size = sum((mes < th), axis = 0)/float(nm) 
    mesr = abs(method_func(A, er, nf = nf))**2 #Valor real da medida
    h11 = sum((mesr < ic1), axis = 0)/float(nm)
    h12 = sum((mesr > ic2), axis = 0)/float(nm)
    return mes, th, ic1, ic2, h0size, h11, h12, mesr

def test_coh():
    '''Para usar o test_asymp, deve-se fornecer uma matriz A e er.
       A deve ser uma matriz (nxnxp), onde p é a ordem do var.
       er deve ser (nxn), sendo a covariância da inovação.
       
       Deve dizer qual estatistica e medida vai testar. Existem:
       pdc -> ass_.asymp_pdc, pdc_.pdc_alg
       dtf -> ass_.asymp_dtf_one, pdc_.dtf_one_alg
       pc -> ass_.asymp_pc, pdc_.pc_alg
       coh -> ass_.asymp_coh, pdc_.coh_alg
       ss -> ass_.asymp_ss, pdc_.ss_alg
       
       Para o PDC deve informar a metrica: metric = 'euc', 'diag' ou 'gen'. 
       Para os outros não deve ser esepecificado.
       
       nm eh numero de iterações do Monte Carlo
       nd eh o tamanho dos dados gerados
       nf eh quantas frequências calcular
       alpha eh tamanho do intervalo de confiança
    '''
       
    A = array([[[0.2, 0],[0, 0],[0.3,-0.2]], 
               [[0, 0],[0.8,-0.1],[0.4,-0.1]],
               [[0, 0],[-0.1,0.2],[0.4,0.1]]], dtype = float) 
    
    er = identity(3)
    #er = array([[0.7,0.3, 0], [0.3, 1.2, 0.4], [0, 0.4, 2]], dtype = float)
    coh, th, ic1, ic2, h0, h11, h12, cohr = test_asymp(ass_.asymp_coh, pdc_.coh_alg, nm = 100, nd = 10000, A = A, er = er, 
                                                       nf = 10, alpha = 0.05)

    #pdc, th, ic1, ic2, h0, h11, h12, pdcr = test_asymp(ass_.asymp_pdc, pdc_.pdc_alg, nm = 100, nd = 10000, A = A, er = er, 
    #                                                   nf = 10, alpha = 0.05, metric = 'diag')
    print h0
    print h11+h12
    

    
def teste_simples():
    A = array([[[0.2, 0],[0.3,-0.2],[0.3,-0.2]], 
                   [[0, 0],[0.8,-0.1],[0.4,-0.1]],
                   [[0, 0],[0.3,0.2],[0.4,0.1]]], dtype = float) 
    er = identity(3)
    nd = 100
    nf = 20
    alpha = 0.05
    n = A.shape[0]
    maxp = A.shape[2]
    metric = 'gen'
    
    #Generate data from AR
    data = ar_data(A, er, nd)
    #Estimate AR parameters with Nuttall-Strand
    Aest, erest = nstrand(data, maxp = maxp)
    #Calculate the connectivity and statistics
    mes, th, ic1, ic2 = ass_.asymp_pdc(data, Aest, nf, erest, 
                                   maxp, alpha = alpha, metric = metric)
    
    pdc_.plot_all(mes, th, ic1, ic2, nf = nf)

def bootstrap(method_func, nd, nm, A, er, 
              nf, alpha = 0.05, metric = None):
    '''Faz histograma e estatistica bootstrap'''
    
    n = A.shape[0]
    mes = empty([nm, n, n, nf])
    #th  = empty([nm, n, n, nf])
    #ic1 = empty([nm, n, n, nf])
    #ic2 = empty([nm, n, n, nf])
    maxp = A.shape[2]
    tbegin = time.clock()
    for i in range(nm):
        if (i%100 == 0):
            print 'nm:', i, 'time:', time.clock()-tbegin
        #Generate data from AR
        data = ar_data(A, er, nd)
        #Estimate AR parameters with Nuttall-Strand
        Aest, erest = nstrand(data, maxp = maxp)
        #Calculate the connectivity and statistics
        if (metric == None):
            mes[i] = abs(method_func(Aest, erest, nf = nf))**2
        else:
            mes[i] = abs(method_func(Aest, erest, nf = nf, metric = metric))**2

    
    so = sort(mes, axis = 0)
    ic1 = so[(alpha/2)*nm]
    ic2 = so[(1-alpha/2)*nm]
    bvar = var(mes, axis = 0)
    
    return mes, bvar, ic1, ic2

def compare_bootstrap_asymp():
    
    A, er = ar_models(1)
    maxp = A.shape[2]    
    nd = 200
    nm = 500
    nf = 5
    alpha = 0.05
    meth = pdc_.pdc_alg
    asymp_func = ass_.asymp_pdc
    metric = 'gen'
    
    #Generate data from AR
    data = ar_data(A, er, nd)
    #Estimate AR parameters with Nuttall-Strand
    Aest, erest = nstrand(data, maxp = maxp)
    
    #Generate data from AR
    data2 = ar_data(A, er, nd)
    #Estimate AR parameters with Nuttall-Strand
    Aest2, erest2 = nstrand(data2, maxp = maxp)
    
    mes = abs(meth(A, er, nf))**2
    mesest = abs(meth(Aest, erest, nf))**2
    
    #er = array([[0.7,0.3, 0], [0.3, 1.2, 0.4], [0, 0.4, 2]], dtype = float)
    mesb, bvar, ic1, ic2 = bootstrap(meth, nd = nd, nm = nm, A = A, er = er, 
                                     nf = nf, alpha = alpha, metric = None)
    mesa, tha, ic1a, ic2a = asymp_func(data, A, nf, er, 
                                       maxp, alpha = alpha)
    
    mesa2, tha2, ic1a2, ic2a2 = asymp_func(data2, A, nf, er, 
                                           maxp, alpha = alpha)
    
    
    
    print mes
    print mesest
    print 'ic1b', ic1
    print 'ic1bv', mesest - sqrt(bvar)*st.norm.ppf(1-alpha/2.0)
    print 'ic1a', ic1a
    print 'ic2b', ic2
    print 'ic2bv', mesest + sqrt(bvar)*st.norm.ppf(1-alpha/2.0)
    print 'ic2a', ic2a
    print bvar.shape
    
    

def test_bootstrap():
    #A = array([[[0.2, 0],[0, 0],[0.3,-0.2]], 
    #           [[0, 0],[0.8,-0.1],[0.4,-0.1]],
    #           [[0, 0],[-0.1,0.2],[0.4,0.1]]], dtype = float) 
    
    #er = identity(3)
    
    A = array([[[4,-4],[3,3]],[[0,0],[0,3]]], dtype=float).reshape(2,2,2)/20
    er = array([[0.7,0],[0,2]], dtype = float)
    
    nd = 500
    nm = 200
    maxp = 2
    nf = 4
    alpha = 0.05
    meth = pdc_.coh_alg
    asymp_func = ass_.asymp_coh
    
    #Generate data from AR
    data = ar_data(A, er, nd)
    #Estimate AR parameters with Nuttall-Strand
    Aest, erest = nstrand(data, maxp = maxp)
    
    mes = abs(meth(A, er, nf))**2
    mesest = abs(meth(Aest, erest, nf))**2
    
    #er = array([[0.7,0.3, 0], [0.3, 1.2, 0.4], [0, 0.4, 2]], dtype = float)
    coh, bvar, ic1, ic2 = bootstrap(meth, nd = nd, nm = nm, A = A, er = er, 
                                     nf = nf, alpha = alpha, metric = None)
    mesa, tha, ic1a, ic2a = asymp_func(data, Aest, nf, erest, 
                                       maxp, alpha = alpha)
    print mes
    print mesest
    print 'ic1b', ic1
    print 'ic1bv', mesest - sqrt(bvar)*st.norm.ppf(1-alpha/2.0)
    print 'ic1a', ic1a
    print 'ic2b', ic2
    print 'ic2bv', mesest + sqrt(bvar)*st.norm.ppf(1-alpha/2.0)
    print 'ic2a', ic2a
    print bvar.shape

def test_bootstrap_MxN_loop(method_func, ns = None, ms = None, A = None, er = None, 
                       nf = 10, alpha = 0.05, metric = None):
    '''Testa qual M é necessário para bootstrap razoável'''
    if A == None:
        A = array([[[0.2, 0],[0, 0],[0.3,-0.2]], 
                   [[0, 0],[0.8,-0.1],[0.4,-0.1]],
                   [[0, 0],[0,0],[0.4,0.1]]], dtype = float) 
    if er == None:
        er = identity(3)
        
    if ns == None:
        ns = array([10, 30, 50, 100, 200, 500, 1000])
    if ms == None:
        ms = array([100, 300, 1000, 3000])
    
    n = A.shape[0]
        
    bvar = empty([size(ns), size(ms), n, n, nf])
    ic1 = empty([size(ns), size(ms), n, n, nf])
    ic2 = empty([size(ns), size(ms), n, n, nf])
    for i in range(size(ns)):
        for j in range(size(ms)):
            dummy, bvar[i,j], ic1[i,j], ic2[i,j] = bootstrap(method_func, ns[i], ms[j], 
                                                             A, er, nf, alpha, metric)
    
    return bvar, ic1, ic2            
            
            
def test_bootstrap_MxN():
    '''Testa qual M é necessário para bootstrap razoável'''
    
    #ns = array([10, 30, 50, 100, 200, 500, 1000])
    #ms = array([100, 300, 1000, 3000])
    
    ns = array([200])
    ms = array([100, 300, 1000, 3000])
    
    A = array([[[4,-4],[3,3]],[[0,0],[0,3]]], dtype=float).reshape(2,2,2)/20
    er = array([[0.7,0],[0,2]], dtype = float)
    
    nf = 4
    
    bvar, ic1, ic2 = test_bootstrap_MxN_loop(pdc_.pdc_alg, ns = ns, ms = ms, nf = nf, A = A, er = er, metric = 'gen')

    print 'bvar', bvar.transpose([0,2,3,4,1])
    #print 'ic1, ic1, ic2
    

if __name__ == "__main__":
    #test_coh()
    #teste_simples()
    #test_bootstrap()
    #test_bootstrap_MxN()
    compare_bootstrap_asymp()
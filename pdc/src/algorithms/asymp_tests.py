#encoding:UTF-8

from numpy import *
import matplotlib.pyplot as pp
from scipy.stats import chi2
import time

import algorithms.asymp as ass_
import algorithms.pdc_alg as pdc_
from data_simulation.ar_data import ar_data
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
    pdc, th, ic1, ic2, h0, h11, h12, pdcr = test_asymp(ass_.asymp_pdc, pdc_.pdc_alg, nm = 100, nd = 10000, A = A, er = er, 
                                                       nf = 10, alpha = 0.05, metric = 'diag')
    #coh, th, ic1, ic2, h0, h11, h12, cohr = test_asymp(ass_.asymp_coh, pdc_.coh_alg, nm = 100, nd = 10000, A = A, er = er, 
    #                                                   nf = 10, alpha = 0.05)

    print h0
    print h11+h12
    
def plot_all(mes, th, ic1, ic2, nf = 64, sample_f = 1.0):
    
    x = sample_f*arange(nf)/(2.0*nf)
    n = mes.shape[0]
    for i in range(n):
        for j in range(n):
            pp.subplot(n,n,i*n+j+1)
            #over = mes[i,j][mes[i,j]>th[i,j]]
            #overx = x[mes[i,j]>th[i,j]]
            over = mes[i,j]
            overx = x
            under = mes[i,j][mes[i,j]<=th[i,j]]
            underx = x[mes[i,j]<=th[i,j]]
            pp.plot(x, th[i,j], 'r:', x, ic1[i,j], 'k:', x, ic2[i,j], 'k:', 
                    overx, over, 'b-', underx, under, 'r-')
            pp.ylim(-0.05,1.05)
            if (i < n-1):
                pp.xticks([])
            if (j > 0):
                pp.yticks([])
        #if (ss != None):
        #    ax = pp.subplot(n,n,i*n+i+1).twinx()
        #    ax.plot(sample_f*arange(nf)/(2.0*nf), ss[i,i,:], color='g')
        #    ax.set_ylim(ymin = 0, ymax = ss[i,i,:].max())
        #    if (i < n-1):
        #        ax.set_xticks([])
    pp.show()
    
def teste_simples():
    A = array([[[0.2, 0],[0.3,-0.2],[0.3,-0.2]], 
                   [[0, 0],[0.8,-0.1],[0.4,-0.1]],
                   [[0, 0],[0.3,0.2],[0.4,0.1]]], dtype = float) 
    er = identity(3)
    nd = 500
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
    x = arange(nf)/(2.0*nf)
    
    plot_all(mes, th, ic1, ic2, nf = nf)
    #pp.plot(x, mes[0,2], 'k-', x, th[0,2], 'r-', x, ic1[0,2], 'b-', x, ic2[0,2], 'b-')
    #pp.show()

    

if __name__ == "__main__":
    #test_coh()
    teste_simples()
    

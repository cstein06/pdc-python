# -*- coding:utf-8 -*-
"""
Created on 28/04/2010

Simulações para qualificação

Testa ditribuicao e tamanho da estat asymp

@author: Carlos Stein
"""

from pdc import *
from pdc.tests.asymp_tests import test_asymp

from numpy import *
from numpy.random import multivariate_normal as mnorm
import matplotlib.pyplot as pp
import scipy.stats as st
import cPickle

pp.rcParams['ps.usedistiller'] = 'xpdf' #para ter texto no .ps salvo da figura
pp.rcParams['text.usetex'] = True


def gen_data1(nd, dummy = 100):
    '''Adaptado de Baccala e Sameshima 2001, alterando conexao 0->1'''
    p = 3
    x = zeros([5,nd+dummy])
    e = mnorm(zeros(7), identity(7), nd+dummy).T
        
    for i in arange(p, nd+dummy):
        x[0,i] = 0.95*sqrt(2)*x[0,i-1] - 0.9025*x[0,i-2] + e[0,i]
        x[1,i] = 0.5*x[1,i-2] + e[1,i]
        x[2,i] = -0.4*x[0,i-3] + e[2,i]
        x[3,i] = -0.5*x[0,i-2] + 0.25*sqrt(2)*x[3,i-1] + 0.25*sqrt(2)*x[4,i-1] + e[3,i]
        x[4,i] = -0.25*sqrt(2)*x[3,i-1] + 0.25*sqrt(2)*x[4,i-1] + e[4,i]

    return x[:,dummy:]

def model1():
    '''Adaptado de Baccala e Sameshima 2001, alterando conexao 0->1'''
    e = identity(5)
    A = zeros([5,5,3])
    A[0,0,0] = 0.95*sqrt(2)
    A[0,0,1] = -0.9025
    A[1,1,1] = 0.5
    A[2,0,2] = -0.4
    A[3,0,1] = -0.5
    A[3,3,0] = 0.25*sqrt(2)
    A[3,4,0] = 0.25*sqrt(2)
    A[4,3,0] = -0.25*sqrt(2)
    A[4,4,0] = 0.25*sqrt(2)
    return A, e


def monte_carlo1(m = 10000, nd = 20000, bignd = 200000):
    '''Roda todas medidas m vezes para model1 e guarda em res. 
       Roda asymp e guarda em resa.
    '''
    
    mes = ['coh', 'pdc', 'dtf', 'pc']
    
    pr_.metric = 'diag'
    pr_.fixp = True
    pr_.maxp = 3
    pr_.ss = False
    pr_.nf = 5
    pr_.v = False
    
    n = 5
    
    Ao, eo = model1()
    
    sumA = 0
    sume = 0
    
    res = zeros([4,m,n,n,pr_.nf])
    for i in arange(m):
        
        if (i+1)%5 == 0:
            print 'iter', i+1
        
        #data = gen_data1(nd)
        data = ar_data(Ao, eo, nd = nd)
        
        A, e = nstrand(data, pr_.maxp)
    
        res[0,i] = abs(coh_alg(A, e, pr_.nf))**2
        res[1,i] = abs(pdc_alg(A, e, pr_.nf))**2
        res[2,i] = abs(dtf_alg(A, e, pr_.nf))**2
        res[3,i] = abs(pc_alg(A, e, pr_.nf))**2
    
        sumA += A
        sume += e
    
        pr_.v = False
        
    Am = sumA/m
    em = sume/m
        
    #med = res.mean(axis = 1)
    #std = res.std(axis = 1)
    
    resa = zeros([4,4,n,n,pr_.nf])
    
    #data = gen_data1(nd)
    data = ar_data(Ao, eo, nd = nd)
    for i in arange(4):
        f = globals()['asymp_' + mes[i]]
        r = f(data, Am, pr_.nf, e_var = em, p = pr_.maxp, metric = pr_.metric)
        resa[i,0] = r[0]
        resa[i,1] = sqrt(res_.varass1)
        resa[i,2] = res_.patden*nd
        resa[i,3] = res_.patdf
    
    return res, resa

def monte_carlo_geral(Ao, eo, m = 10000, nd = 20000, bignd = 200000):
    '''Roda m vezes e guarda em res. 
       Roda asymp e guarda em resa.
       setar pr_ antes de rodas
    '''
    
    pr_.ss = False
    
    n = Ao.shape[0]
    
    sumA = 0
    sume = 0
    
    mes = globals()[pr_.alg + '_alg']
    asy = globals()['asymp_' + pr_.alg]
    
    res = zeros([m,n,n,pr_.nf])
    for i in arange(m):
        
        if (i+1)%5 == 0:
            print 'iter', i+1
        
        #data = gen_data1(nd)
        data = ar_data(Ao, eo, nd = nd)
        
        A, e = nstrand(data, pr_.maxp)
    
        res[i] = abs(mes(A, e, pr_.nf))**2
    
        sumA += A
        sume += e
    
        pr_.v = False
        
    Am = sumA/m
    em = sume/m
        
    #med = res.mean(axis = 1)
    #std = res.std(axis = 1)
    
    resa = zeros([4,n,n,pr_.nf])
    
    #data = gen_data1(nd)
    data = ar_data(Ao, eo, nd = nd)
    r = asy(data, Am, pr_.nf, e_var = em, p = pr_.maxp, metric = pr_.metric)
    resa[0] = r[0]
    resa[1] = sqrt(res_.varass1)
    resa[2] = res_.patden*nd
    resa[3] = res_.patdf
    
    return res, resa
            
def qqnorm(res, resa, ax):
    '''qqplot para H1, com a normal
    ja espera apenas uma frequencia e conexao.
    res tem amostral, resa tem asymp (mes, var1, patden, patdf)'''
    
    m = res.shape[0]
    
    x = linspace(1.0/m, 1-1.0/m, m)
    y = st.norm.ppf(x, loc = resa[0], scale = resa[1])
    
    tic = array(['0.01', '0.05', '0.25', '0.75', '0.95', '0.99'])
    ytic = y[int32(tic.astype(float)*m)]
    
#    xmi = min(y.min(), res[:,3,0,2].min())
#    xma = max(y.max(), res[:,3,0,2].max())
    xmi = res.min()
    xma = res.max()

    ax.plot(sorted(res), y, 'k+')
    pp.plot([xmi,xma],[xmi,xma], 'b')
    
    ax.yaxis.set_ticklabels(tic)
    ax.yaxis.set_ticks(ytic)
    
    a,b = pp.xlim()
    a1 = a + 0.2*(b-a)
    a2 = a + 0.8*(b-a)
    ax.xaxis.set_ticks([a1,a2])
    
    pp.ylim(y[0], y[-1])
    
def qqchi2(res, resa, ax):
    '''qqplot para H0, com a chi2 ponderada
    ja espera apenas uma frequencia e conexao.
    res tem amostral, resa tem asymp (mes, var1, patden, patdf)'''
    
    m = res.shape[0]
    
    x = linspace(1.0/m, 1-1.0/m, m)
    y = st.chi2.ppf(x, resa[3], scale = 1/(resa[2]*2.0))
    
    tic = array(['0.01', '0.50', '0.75', '0.95', '0.99'])
    ytic = y[int32(tic.astype(float)*m)]
    
    #xmi = min(y.min(), res[:,4,0,2].min())
    #xma = max(y.max(), res[:,4,0,2].max())
    xmi = res.min()
    xma = res.max()
    
    pp.plot(sorted(res), y, 'k+')
    pp.plot([xmi,xma],[xmi,xma], 'r')
    
    ax.yaxis.set_ticklabels(tic)
    ax.yaxis.set_ticks(ytic)
    
    a,b = pp.xlim()
    a1 = 0.0
    a2 = 0.8*b
    ax.xaxis.set_ticks([a1,a2])
    
    pp.ylim(y[0], y[-1])
    
def qqplots1(res, resa):
    '''plota figura para o texto com quantis de cada medida 
    para 3 pares, em 0.2hz.
    separa em dois graficos'''
    
    mes = ['coh', 'pdc', 'dtf', 'pc']
    
    f = 2
    
    pares = array([[1,0], [3,0], [4,0]])
    
    fig = pp.figure()
    
    for i in arange(4):
        for j in arange(3):
            ax = fig.add_subplot(3,4,i+j*4+1)
            qqnorm(res[i,:,pares[j,0],pares[j,1],f], 
                   resa[i,:,pares[j,0],pares[j,1],f], 
                   ax)
    
            if j == 2:
                pp.xlabel(mes[i].upper())
                #todo colocar titulo sampled, mini titulos com alg
            if i == 0:
                pp.ylabel('Quantile for Normal')
                ax.yaxis.set_ticks([])
                #todo colocar titulo quantile, mini titulos com pares
    
    fig = pp.figure()
    
    for i in arange(4):
        for j in arange(3):
            ax = fig.add_subplot(3,4,i+j*4+1)
            qqchi2(res[i,:,pares[j,0],pares[j,1],f], 
                   resa[i,:,pares[j,0],pares[j,1],f], 
                   ax)
            
            if j == 2:
                pp.xlabel(mes[i].upper())
            if i == 0:
                pp.ylabel('Quantile for Chi^2')
            else:
                ax.yaxis.set_ticks([])
    
#    
    pp.show()
    
    
def qqplots1_oneplot(res, resa):
    '''plota figura para o texto com quantis de cada medida 
    para 3 pares, em 0.2hz.
    plota normal para significativos, chi2 para h0 (linha vermelha)'''
    
    #seria melhor usar A e E originais, faz mais sentido
    #usar original testa teoria, media testa estimacao da estatistica
#    resa = zeros([4,4,n,n,pr_.nf])
#    
#    for i in arange(4):
#        f = globals()['asymp_' + mes[i]]
#        r = f(data, Am, pr_.nf, e_var = em, p = pr_.maxp, metric = pr_.metric)
#        resa[i,0] = r[0]
#        resa[i,1] = sqrt(res_.varass1)
#        resa[i,2] = res_.patden*nd
#        resa[i,3] = res_.patdf
        
    
    
    mes = ['coh', 'pdc', 'dtf', 'pc']
    
    f = 2
    
    pares = array([[1,0], [3,0], [4,0]])
    
    fig = pp.figure()
    
    for i in arange(4):
        for j in arange(3):
            ax = fig.add_subplot(3,4,i+j*4+1)
            
            r1 = res[i,:,pares[j,0],pares[j,1],f]
            ra1= resa[i,:,pares[j,0],pares[j,1],f]
            if ra1[0] > 1E-2:
                qqnorm(r1, ra1, ax)
            else:
                qqchi2(r1, ra1, ax)
    
            if j == 2:
                pp.xlabel(mes[i].upper())
                #todo colocar titulo sampled, mini titulos com alg
            if i == 0:
                pp.ylabel(("$%d \\to %d$") % (pares[j,1]+1,pares[j,0]+1))
            else:
                ax.yaxis.set_ticks([])
                #todo colocar titulo quantile, mini titulos com pares
    
    #pp.show()
    
    
def hists1(res, resa):
    '''compara histograma sampled com asintotica'''
    
    res1 = res[0,:,3,0,2]
    resa1 = resa[0,:,3,0,2]
    
    pp.figure()
    
    pp.suptitle(r'Sampled Coherence histogram and asymptotic statistics')
    
    pp.subplot(1,2,1)
    
    bins = 60
    
    pp.hist(res1, bins = bins, normed = True)
    
    x = linspace(resa1[0]-3*resa1[1], resa1[0]+3*resa1[1], 300)
    
    
    pp.plot(x, st.norm.pdf(x, loc = resa1[0], scale = resa1[1]))
    
    pp.title(r'$1 \to 4 (Normal)$')
    
    pp.subplot(1,2,2)
    #pp.figure()
    
    res1 = res[0,:,1,0,2]
    resa1 = resa[0,:,1,0,2]
    
    pp.hist(res1, bins = bins, normed = True)
    
    x = linspace(0, 3*resa1[3]/resa1[2], 300)
    
    pp.plot(x, st.chi2.pdf(x, resa1[3], scale = 1.0/(resa1[2]*2)))
    
    pp.title(r'$1 \to 5 (weigthed \chi^2$)')
        
    #pp.draw()

root = '/home/stein/producao/quali/simulados/'

def tab1(res, resa):
    '''acha tamanho do erro I sob h0 e h1'''
    
    m = res.shape[1]
    
    patdf = resa[:,3,1,0,2]
    patden = resa[:,2,1,0,2]
    
    alpha = 0.01
    th1 = st.chi2.ppf(1-alpha, patdf)/(patden*2)
    
    alpha = 0.05
    th5 = st.chi2.ppf(1-alpha, patdf)/(patden*2)
    
    mes = resa[:,0,3,0,2]
    pst = resa[:,1,3,0,2]
    
    alpha = 0.01
    ic111 = mes - pst*st.norm.ppf(1-alpha/2.0)
    ic211 = mes + pst*st.norm.ppf(1-alpha/2.0)
    
    alpha = 0.05
    ic115 = mes - pst*st.norm.ppf(1-alpha/2.0)
    ic215 = mes + pst*st.norm.ppf(1-alpha/2.0)
    
    mes = resa[:,0,4,0,2]
    pst = resa[:,1,4,0,2]
    
    alpha = 0.01
    ic121 = mes - pst*st.norm.ppf(1-alpha/2.0)
    ic221 = mes + pst*st.norm.ppf(1-alpha/2.0)
    
    alpha = 0.05
    ic125 = mes - pst*st.norm.ppf(1-alpha/2.0)
    ic225 = mes + pst*st.norm.ppf(1-alpha/2.0)
    
    #li = [th1, th5, ic11, ic21, ic15, ic25]
    
    #print li
    eth1 = zeros(4)
    eth5 = zeros(4)
    eic11 = zeros(4)
    eic15 = zeros(4)
    eic21 = zeros(4)
    eic25 = zeros(4)
    for i in arange(4):
        eth1[i] = sum(res[i,:,1,0,2] > th1[i])/float(m)
        eth5[i] = sum(res[i,:,1,0,2] > th5[i])/float(m)
        eic11[i] = sum((res[i,:,3,0,2] < ic111[i]) + (res[i,:,3,0,2] > ic211[i]))/float(m)
        eic15[i] = sum((res[i,:,3,0,2] < ic115[i]) + (res[i,:,3,0,2] > ic215[i]))/float(m)
        eic21[i] = sum((res[i,:,4,0,2] < ic121[i]) + (res[i,:,4,0,2] > ic221[i]))/float(m)
        eic25[i] = sum((res[i,:,4,0,2] < ic125[i]) + (res[i,:,4,0,2] > ic225[i]))/float(m)
    
    r = [eth1, eth5, eic11, eic15, eic21, eic25]
    
    print array(r).T
    
    return r

def size_and_power(res, resa, alpha = 0.05):
    '''acha tamanho do erro I sob h0 (ou poder) e h1, variando o nd'''
    
    m = res.shape[0]
    
    patdf = resa[3]
    patden = resa[2]
    
    th = st.chi2.ppf(1-alpha, patdf)/(patden*2)
    
    mes = resa[0]
    pst = resa[1]
    
    ic1 = mes - pst*st.norm.ppf(1-alpha/2.0)
    ic2 = mes + pst*st.norm.ppf(1-alpha/2.0)
    
    eth = sum(res[:] > th)/float(m)
    eic = sum((res[:] < ic1) + (res[:] > ic2))/float(m)
        
    r = [eth, eic]
    
    print 'size/power', array(r)
    
    return r

     
           
def tab1_hard():
    #testa estimacao da estatistica por dados, nao a estatistica exata
    asymp = asymp_coh
    method_func = coh_alg
    
    A,e = model1()
    
    res = test_asymp(asymp, method_func, nm = 200, nd = 100, 
                     A = A, er = e, 
                     nf = 5, alpha = 0.05, metric = None)
    
    print 1-res[4][1,0,2], res[5][3,0,2], res[6][3,0,2]
    
            
def figuras1():      
    #res, resa = monte_carlo1(m = 200, nd = 2000, bignd = 2000)
    [res, resa, du, du, du] = load_data1()
    
    #qqplots1(res,resa)
    qqplots1_oneplot(res, resa)
    hists1(res,resa)
    tab1(res,resa)
    table_varios_nd()
    #tab1_hard()
    
    pp.show()
    
def load_data1(file = 'data1.pic'):
    f = open(root + file, 'rb')
    res = cPickle.load(f)
    resa = cPickle.load(f)
    m = cPickle.load(f)
    nd = cPickle.load(f)
    bignd = cPickle.load(f)
    return res, resa, m, nd, bignd
    
def save_data1(m = 10000, nd = 20000, bignd = 200000, file = 'data1.pic'):
    res, resa = monte_carlo1(m = m, nd = nd, bignd = bignd)
    
    f = open(root + 'mini' + file, 'wb')
    res[:,:300].dump(f)
    resa.dump(f)
    cPickle.dump(m, f)
    cPickle.dump(nd, f)
    cPickle.dump(bignd, f)
    f.close()
    
    f = open(root + file, 'wb')
    res.dump(f)
    resa.dump(f)
    cPickle.dump(m, f)
    cPickle.dump(nd, f)
    cPickle.dump(bignd, f)
    f.close()

def save_data_geral(res, resa, m, nd, bignd, file = 'data1.pic'):
    
    f = open(root + 'mini' + file, 'wb')
    res[:,:300].dump(f)
    resa.dump(f)
    cPickle.dump(m, f)
    cPickle.dump(nd, f)
    cPickle.dump(bignd, f)
    f.close()
    
    f = open(root + file, 'wb')
    res.dump(f)
    resa.dump(f)
    cPickle.dump(m, f)
    cPickle.dump(nd, f)
    cPickle.dump(bignd, f)
    f.close()

nds = [50, 100, 200, 500, 2000, 5000]

def table_varios_nd():
        
    for i in arange(len(nds)):
        [res, resa, du, du, du] = load_data1('data1_nd' + str(nds[i]) + '.pic')
        
        for j in arange(4):
            print 'alg', j, 'nd', nds[i]
            size_and_power(res[j,:,1,0,2], resa[j,:,1,0,2], alpha = 0.01)
            size_and_power(res[j,:,3,0,2], resa[j,:,3,0,2], alpha = 0.01)
            size_and_power(res[j,:,4,0,2], resa[j,:,4,0,2], alpha = 0.01)
        
def save_varios_nd(m = 10000, bignd = 200000):
    
    for i in arange(len(nds)):
        save_data1(m = m, bignd = bignd, nd = nds[i], file = 'data1_nd' + str(nds[i]) + '.pic')

def simple_plot1():
    
    A,e = model1()
    data = ar_data(A,e)
    pr_.alg = 'pdc'
    measure_full(data)
    
    pp.show()
    
def all_plots1():
    
    A,e = model1()
    pr_.reuse_A = True
    res_.A = A
    res_.er = e
    
    pp.figure()
    data = ar_data(A,e)

    pr_.alg = 'coh'
    measure_and_plot(data)
    
    pp.figure()
#    data = ar_data(A,e)

    pr_.alg = 'pdc'
    measure_and_plot(data)
    
    pp.figure()
#    data = ar_data(A,e)

    pr_.alg = 'dtf'
    measure_and_plot(data)
    
    pp.figure()
#    data = ar_data(A,e)

    pr_.alg = 'pc'
    measure_and_plot(data)
    
    
    pp.show()
    
def gen_big_Ar(p = 100, n = 3):
    Aa = ((rand([n,n])<0.25)*rand([n,n])*0.2 + 
          identity(n)*rand([n,n])*0.4)
         
    ra = 0.3
    
    A = zeros(n,n,p)
    
    for i in arange(p):
        A[:,:] = Aa * ra**i * randn([n,n])/3
    
    return Aa

ps = [5, 10, 20, 50]
ns = [3]

a = b43.tr()

def test_varios_p():
    n = 3
    A = ps
    for i in ps:
        3
        

if __name__ == '__main__':
    pass#
    
    figuras1()
    #all_plots1()
    #save_data1()
    #print 'loading'
    #r = load_data1()
    #print r[3]
    #save_varios_nd(m = 10000, bignd = 200000)
    #table_varios_nd()
    #simple_plot1()
    
    

            
            
    
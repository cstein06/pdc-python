'''
Created on Jun 30, 2010

@author: stein
'''

# -*- coding:utf-8 -*-

# In this file we calculate the asymptotic statistics for all measures, including 
# first and second order asymptotic aproximations.

__all__ = ['asymp_pdc', 'asymp_dtf', 'asymp_coh', 'asymp_pc', 'asymp_ss']

from numpy import *
import scipy.stats as st
from scipy.stats import cov as cov
from scipy.linalg import cholesky
from scipy.linalg import eigh
from scipy.linalg import inv as inv
from matplotlib.pyplot import xcorr
from numpy.dual import eig
from numpy.linalg.linalg import LinAlgError
import time
import sys

from pdc.params import pr_, res_
from pdc.analysis import ss_alg

# These functions are used to make code more readable.
vec = lambda x: mat(x.ravel('F')).T
O = lambda n: mat(zeros([n, n], dtype=float))
I = lambda n: mat(identity(n, dtype=float))
cat = lambda a, b, ax: concatenate((a, b), axis = ax)
mdiag = lambda a: mat(diag(diag(a)))
diagtom = lambda a: mat(diag(array(a).reshape(-1)))
#corrlag = lambda a,b,lag: cov(a[k:],b[:-k])/a.size
#GInv = lambda a: dot(inv(dot(a.T, a)),a.T)


def A_to_f(A, nf = 64):
    '''Calculates A(f), in the frequency domain.
    
    Input:
        A(n, n, r) - recurrence matrix (n - number of signals, r - model order)
        nf - frequency resolution
        
    Output:
        AL(nf, n, n)
    '''
    
    n, n, r = A.shape
    
    # expoentes eh o array de expoentes da fft, com todas frequencias por todos lags
    expoentes = (-1j*pi*kron(arange(nf),(arange(r)+1.0))/nf).reshape(nf,r)
    # Af faz a multiplicacoes da exp(ar) pela matrix A, para todas frequencias
    # as funcoes repeat e transpose sao truques para possibilitar o calculo vetorial
    Af = (A.reshape(n,n,1,r).repeat(nf, axis=2)*exp(expoentes)).transpose([2,0,1,3])
    # o fft soma o valor para todos os lags
    AL = eye(n) - sum(Af, axis = 3)
    
    return AL

def vech(a):
    ''' Returns vech(A)'''
    n = a.shape[0]
    v = empty((n*(n+1))/2)
    cont = 0
    for j in range(n):
        for i in range(n):
            if i >= j:
                v[cont] = a[i,j]
                cont = cont + 1
    return v
            
def Dup(n):
    '''D*vech(A) = vec(A), with symmetric A'''
    d = matrix(zeros([n*n, (n*(n+1))/2]))
    count = 0
    for j in range(n):
        for i in range(n):
            if i >= j:
                d[j*n+i, count]=1
                count = count+1
            else:
                d[j*n+i, :]=d[i*n+j, :]
    return d

def TT(a,b):
    ''' TT(a,b)*vec(B) = vec(B.T), where B is (a x b).'''
    t = O(a*b)
    for i in range(a):
        for j in range(b):
            t[i*b+j, j*a+i] = 1 
    return t

def fTe(n):
    '''Derivative of kron(I(2n), A) by A'''
    return dot(kron(TT(2*n, n), I(n*2*n)),
           kron(I(n), kron(vec(I(2*n)), I(n))))
    

def fTe2(n):
    '''Derivative of kron(I(2), A) by A'''
    return dot(kron(TT(2, n), I(n*2)),
           kron(I(n), kron(vec(I(2)), I(n))))

def xlag(x, lag):
    if(lag == 0):
        return x.copy();
    xl = matrix(zeros(x.shape))
    xl[:, lag:] = x[:, :-lag]
    return xl

def xlags(x, lags):
    
    xl = zeros([lags] + list(x.shape))
    
    xl[0] = x
    for i in arange(1,lags):
        xl[i, :, i:] = x[:, :-i]
        
    return xl

def bigautocorr(x, p):
    '''Autocorrelation. Data in rows. From order 0 to p-1.
    Output: nxn blocks of autocorr of lags i. (Nuttall Strand matrix)'''
    n, nd = x.shape
    gamma = zeros([n*p, n*p])
    xl = xlags(x, p)
    for i in arange(p):
        for j in arange(p):
            gamma[i*n:i*n+n, j*n:j*n+n] = dot(xl[i], xl[j].T)/nd

    return gamma
    
def fdh_da(Af, n):
    '''Derivative of vec(H) by vec(A), with H = A^-1 and complex A.'''
    ha = Af.I
    h = -kron(ha.T, ha)
    
    h1 = cat(h.real, -h.imag, 1)
    h2 = cat(h.imag, h.real, 1)
    hh = cat(h1, h2, 0)
    return hh
    
def fIij(i, j, n):
    '''Returns Iij of the formula'''
    Iij = zeros(n**2)
    Iij[n*j+i] = 1
    Iij = diag(Iij)
    return mat(kron(I(2), Iij))

def fIj(j, n):
    '''Returns Ij of the formula'''
    Ij = zeros(n)
    Ij[j] = 1
    Ij = diag(Ij)
    Ij = kron(Ij, I(n))
    return mat(kron(I(2), Ij))

def fIi(i, n):
    '''Returns Ii of the formula'''
    Ii = zeros(n)
    Ii[i] = 1
    Ii = diag(Ii)
    Ii = kron(I(n), Ii)
    return mat(kron(I(2), Ii))

def fCa(f, p, n):
    '''Returns C* of the formula'''
    C1 = cos(-2*pi*f*mat(arange(1, p+1)))
    S1 = sin(-2*pi*f*mat(arange(1, p+1)))    
    C2 = cat(C1.reshape(1, -1), S1.reshape(1, -1), 0)
    return kron(C2, identity(n**2))    

def fChol(omega):
    # Try Cholesky factorization
    try:
        L = mat(cholesky(omega, lower=1))
    # If there's a small negative eigenvalue, diagonalize
    except LinAlgError:
        val, vec = eigh(omega)
        if sum(val<-1E-10) > 0:
            print 'negative eig. in omega:', val[val<0]
        #print omega
        L = zeros(vec.shape)
        for i in range(len(val)):
            if val[i]<0.:
                val[i]=0.
            L[:,i] = vec[:,i]*sqrt(val[i])
        #print 'L', L
        
    return L

def fEig(L, G2):
    '''Returns the eigenvalues'''

    #L = mat(cholesky(omega, lower=1))
    D = L.T*G2*L
    d = eigh(D, eigvals_only=True)
    d = d[abs(d) > 1E-8]
    if (d.size > 2):
        print 'more than two chi-square in the sum:'
        print d
    return d

def asymp_pdc(x):
    '''Asymptotic statistics for the three PDC formulations
        x -> data
        res_.A -> autoregressive matrix
        res_.er -> residues
    '''
    nf = pr_.nf
    n,n,p = res_.A.shape
    n, nd = x.shape
    
    x = mat(x)
    er = mat(res_.er)
    Af = A_to_f(res_.A, nf)
    
    
    res = empty([8, n, n, nf])
    
    gamma = mat(bigautocorr(x, p))
    omega = kron(gamma.I, er)

    if pr_.metric == 'orig':
        Sn = kron(I(2*n), I(n))
        Sd = kron(I(2*n), I(n))
        En = 0
        Ed = 0
    elif pr_.metric == 'gen':
        Sn = kron(I(2*n), mdiag(er))
        Sd = kron(I(2*n), mdiag(er))
        En = -diagtom(vec(Sn.I*Sn.I))
        Ed = -diagtom(vec(Sn.I*Sn.I))
    elif pr_.metric == 'info':
        Sn = kron(I(2*n), mdiag(er))
        Sd = kron(I(2*n), er)
        En = -diagtom(vec(Sn.I*Sn.I))
        Ed = -kron(Sd.I.T, Sd.I)
    else:
        print 'metric doesn\'t exist. using gen.'
        Sn = kron(I(2*n), mdiag(er))
        Sd = kron(I(2*n), mdiag(er))
        En = 0
        Ed = 0
    
    Teta = fTe(n)
    EnT = En*Teta
    EdT = Ed*Teta
    D = Dup(n)
    omega_e = 2*D*D.I*kron(er, er)*D.I.T*D.T
    
    for ff in range(nf):
        #print 'ff', ff
        
        f = ff/(2.0*nf)
        Ca = fCa(f, p, n)
        omega_a = Ca*omega*Ca.T
        
        L = fChol(omega_a)
            
        a = vec(Af[ff, :, :])
        a = cat(a.real, a.imag, 0)
        
        for i in range(n):
            for j in range(n):
                Iij = fIij(i, j, n)
                Ij = fIj(j, n)           
                
                num = a.T*Iij*Sn.I*Iij*a
                den = a.T*Ij*Sd.I*Ij*a
                
                pdc = num/den
                Ge = pdc*(kron((Iij*a).T, a.T*Iij)*EnT/num - 
                          kron((Ij*a).T, a.T*Ij)*EdT/den)
                vare = Ge*omega_e*Ge.T
               
                Ga = 2*pdc*(a.T*Iij*Sn.I*Iij/num - 
                            a.T*Ij*Sd.I*Ij/den)
                vara = Ga*omega_a*Ga.T
                
                var1 = (vara + vare)/nd
                ic1 = pdc - sqrt(var1)*st.norm.ppf(1-pr_.alpha/2.0)
                ic2 = pdc + sqrt(var1)*st.norm.ppf(1-pr_.alpha/2.0)
                
                G2a = 2*Iij*Sn.I*Iij/den                
                #print omega, eigh(omega, eigvals_only=True)
                d = fEig(L, G2a)
                patdf = sum(d)**2/sum(d**2)
                patden = sum(d)/sum(d**2)
                
                th = st.chi2.ppf(1-pr_.alpha, patdf)/(patden*2*nd)
                var2 = 2*patdf/(patden*2*nd)**2
                
                res[:,i,j,ff] = [pdc[0,0], th, ic1, ic2, patdf, patden, var1, var2]
        
    res_.mes = res[0]
    res_.th = res[1]
    res_.ic1 = res[2]
    res_.ic2 = res[3]
    res_.asy = res[4:]

    return res



def asymp_pdt(x):
    '''Asymptotic statistics for the three PDC formulations
        x -> data
        res_.A -> autoregressive matrix
        res_.er -> residues
    '''
    nf = pr_.nf
    n,n,p = res_.A.shape
    n, nd = x.shape
    
    x = mat(x)
    er = mat(res_.er)
    Af = A_to_f(res_.A, nf)
    
    res = zeros([8, n, n, nf])
    
    gamma = mat(bigautocorr(x, p))
    omega = kron(gamma.I, er)
    

    for ff in range(nf):
        #print 'ff', ff
        
        f = ff/(2.0*nf)
        Ca = fCa(f, p, n)
        omega_a = Ca*omega*Ca.T
        
        L = fChol(omega_a)
            
        a = vec(Af[ff, :, :])
        a = cat(a.real, a.imag, 0)
        
        H = mat(Af[ff]).I
        f = H*er*H.T.conj()
        
        for i in range(n):
            su = 0.0
            for j in arange(n):
                su += (abs(Af[ff,i,j])**2)*f[j,j]
            su += er[i,i]
            
            for j in range(n):
                
                Iij = fIij(i, j, n)
                
                pdt = (abs(Af[ff,i,j])**2)*f[j,j]/su
                
                G2a = 2*Iij*f[j,j]/f[i,i]              
                #print omega, eigh(omega, eigvals_only=True)
                d = fEig(L, G2a)
                patdf = sum(d)**2/sum(d**2)
                patden = sum(d)/sum(d**2)
                
                th = st.chi2.ppf(1-pr_.alpha, patdf)/(patden*2*nd)
                var2 = 2*patdf/(patden*2*nd)**2
                
                res[:,i,j,ff] = [pdt, th, 0, 0, patdf, patden, 0, var2]
                
    res_.mes = res[0]
    res_.th = res[1]
    res_.asy = res[4:]

    return res


def asymp_dtf(x):
    '''Asymptotic statistics for the DTF
    '''
    nf = pr_.nf
    n,n,p = res_.A.shape
    n,nd = x.shape
    
    x = mat(x)
    er = mat(res_.er)
    Af = A_to_f(res_.A, nf)    
    
    res = empty([8, n, n, nf])
    
    gamma = mat(bigautocorr(x, p))
    omega = kron(gamma.I, er)

    if pr_.metric == 'dtf':
        Sn = kron(I(2*n), I(n))
        Sd = kron(I(2*n), I(n))
        En = 0
        Ed = 0
    elif pr_.metric == 'gen':
        Sn = kron(I(2*n), mdiag(er))
        Sd = kron(I(2*n), mdiag(er))
        En = -diagtom(vec(Sn.I*Sn.I))
        Ed = -diagtom(vec(Sn.I*Sn.I))
    elif pr_.metric == 'dc':
        Sn = kron(I(2*n), mdiag(er))
        Sd = kron(I(2*n), er)
        En = -diagtom(vec(Sn.I*Sn.I))
        Ed = -kron(Sd.I.T, Sd.I)
    else:
        print 'metric doesn\'t exist. using gen.'
        Sn = kron(I(2*n), mdiag(er))
        Sd = kron(I(2*n), mdiag(er))
        En = 0
        Ed = 0
    
    Teta = fTe(n)
    EnT = En*Teta
    EdT = Ed*Teta
    D = Dup(n)
    omega_e = 2*D*D.I*kron(er, er)*D.I.T*D.T
    
    for ff in range(nf):
        f = ff/(2.0*nf)
        Ca = fCa(f, p, n)
        
        Hf = mat(Af[ff, :, :]).I
        h = vec(Hf)
        h = cat(h.real, h.imag, 0)
        
        dhda = fdh_da(mat(Af[ff, :, :]), n)
        omega_h = dhda*Ca*omega*Ca.T*dhda.T
        L = fChol(omega_h)
       
        for i in range(n):
            for j in range(n):
                Iij = fIij(i, j, n)
                Ii = fIi(i, n)
                
                num = h.T*Iij*Sn*Iij*h
                den = h.T*Ii*Sd*Ii*h
                dtf = num/den
    
                Ge = dtf*(kron((Iij*h).T, h.T*Iij)*EnT/num - 
                          kron((Ii*h).T, h.T*Ii)*EdT/den)
                vare = Ge*omega_e*Ge.T
                
                Gh = 2*dtf*(h.T*Iij*Sn*Iij/num - 
                            h.T*Ii*Sd*Ii/den)
                vara = Gh*omega_h*Gh.T
                
                var1 = (vara + vare)/nd
                ic1 = dtf - sqrt(var1)*st.norm.ppf(1-pr_.alpha/2)
                ic2 = dtf + sqrt(var1)*st.norm.ppf(1-pr_.alpha/2)
                
                G2 = 2*Iij*Sn*Iij/den   
                
                d = fEig(L, G2)
                patdf = sum(d)**2/sum(d**2)
                patden = sum(d)/sum(d**2)
                
                th = st.chi2.ppf(1-pr_.alpha, patdf)/(patden*2*nd)
                var2 = 2*patdf/(patden*2*nd)**2
                
                res[:,i,j,ff] = [dtf, th, ic1, ic2, patdf, patden, var1, var2]
                
    res_.mes = res[0]
    res_.th = res[1]
    res_.ic1 = res[2]
    res_.ic2 = res[3]
    res_.asy = res[4:]
    
    return res

def fc(i, n):
    vi = zeros(n)
    vi[i] = 1
    vi.resize(n,1)
    return mat(kron(vi, I(n)))

def fl(i, n):
    vi = zeros(n)
    vi[i] = 1
    vi.resize(n,1)
    return mat(kron(I(n), vi))

def asymp_pc(x):
    '''Asymptotic statistics for the PC
        x -> data
        res_.A -> autoregressive matrix
        res_.er -> residues
    '''
    nf = pr_.nf
    n,n,p = res_.A.shape
    n,nd = x.shape
    
    x = mat(x)
    er = mat(res_.er)
    Af = A_to_f(res_.A, nf)    
    
    res = empty([8, n, n, nf])
    
    gamma = mat(bigautocorr(x, p))
    omega_a = kron(gamma.I, er)
    
    D = Dup(n)
    omega_e = 2*D*D.I*kron(er, er)*D.I.T*D.T
    
    S = kron(I(2*n), er)
    Teta = fTe(n)
    E = -kron(S.I.T, S.I)
    ET = E*Teta
    
    for ff in range(nf):
        #print 'ff', ff
        f = ff/(2.0*nf)
        Ca = fCa(f, p, n)
        
        a = vec(Af)
        a = cat(a.real, a.imag, 0)
        
        omega_a = Ca*omega_a*Ca.T
        
        Rf = empty([n,n])
        If = empty([n,n])
        RGa = empty([n,n,size(a)])
        RGe = empty([n,n,size(a)])
        IGa = empty([n,n,size(a)])
        IGe = empty([n,n,size(a)])
        for i in range(n):
            for j in range(n):
                Ii = fIj(i, n)
                Ij = fIj(j, n)
                Ija = fIjast(j, n)
                
                Rf[i,j] = a.T*Ii*S*Ij*a
                If[i,j] = a.T*Ii*S*Ija*a
                
                RGa[i,j] = a.T*(Ii*S*Ij + Ij*S*Ii)
                RGe[i,j] = kron((Ij*a).T, a.T*Ii)*ET
                IGa[i,j] = a.T*(Ii*S*Ija + Ija*S*Ii)
                IGe[i,j] = kron((Ija*a).T, a.T*Ii)*ET
                
        for i in range(n):
            for j in range(n):
                pc = (Rf[i,j]**2 + If[i,j]**2)/(Rf[i,i]*Rf[j,j])
                
                omega_fh = vstack((RGa[i,j], IGa[i,j], 
                                   RGa[i,i], RGa[j,j]))
                omega_fe = vstack((RGe[i,j], IGe[i,j], 
                                   RGe[i,i], RGe[j,j]))
                omega_f = (omega_fh*omega_a*omega_fh.T +  
                          omega_fe*omega_e*omega_fe.T) 
            
                pcG = array([2*Rf[i,j]/(Rf[i,i]*Rf[j,j]), 2*If[i,j]/(Rf[i,i]*Rf[j,j]),
                              -pc/Rf[i,i], -pc/Rf[j,j]])
                
                var1 = pcG*omega_f*pcG.T   
                ic1 = pc - sqrt(var1)*st.norm.ppf(1-pr_.alpha/2.0)
                ic2 = pc + sqrt(var1)*st.norm.ppf(1-pr_.alpha/2.0)
                
                omega_f = omega_f[:2,:2]
                
                pcG2 = array([[2/(Rf[i,i]*Rf[j,j]),0], 
                               [0,2/(Rf[i,i]*Rf[j,j])]])
                
                L = fChol(omega_f)        
                d = fEig(L, pcG2)
                #d1 = fEig(omega, G2)
                #d2 = fEig(omega_evar, d2coh_dev2)
                #d = concatenate((d1,d2)) #TODO: conferir
                patdf = sum(d)**2/sum(d**2)
                patden = sum(d)/sum(d**2)
                th = st.chi2.ppf(1-pr_.alpha, patdf)/(patden*2*nd)
                var2 = 2*patdf/(patden*2*nd)**2
                
                res[:,i,j,ff] = [pc[0,0], th, ic1, ic2, patdf, patden, var1, var2]
                
    res_.mes = res[0]
    res_.th = res[1]
    res_.ic1 = res[2]
    res_.ic2 = res[3]
    res_.asy = res[4:]
    
    return res



def asymp_coh(x):
    '''Asymptotic statistics for the Coherence (coh)
        x -> data
        res_.A -> autoregressive matrix
        res_.er -> residues
    '''
    nf = pr_.nf
    n,n,p = res_.A.shape
    n,nd = x.shape
    
    x = mat(x)
    er = mat(res_.er)
    Af = A_to_f(res_.A, nf)    
    
    res = empty([8, n, n, nf])
    
    gamma = mat(bigautocorr(x, p))
    omega_a = kron(gamma.I, er)
    
    D = Dup(n)
    omega_e = 2*D*D.I*kron(er, er)*D.I.T*D.T
    
    S = kron(I(2*n), er)
    Teta = fTe(n)
    #E = ?
    ET = E*Teta
    
    for ff in range(nf):
        #print 'ff', ff
        f = ff/(2.0*nf)
        Ca = fCa(f, p, n)
        
        Hf = mat(Af[ff, :, :]).I
        h = vec(Hf)
        h = cat(h.real, h.imag, 0)
        
        dhda = fdh_da(mat(Af[ff, :, :]), n)
        
        omega_h = dhda*Ca*omega_a*Ca.T*dhda.T
        
        Rf = empty([n,n])
        If = empty([n,n])
        RGh = empty([n,n,size(h)])
        RGe = empty([n,n,size(h)])
        IGh = empty([n,n,size(h)])
        IGe = empty([n,n,size(h)])
        for i in range(n):
            for j in range(n):
                Ii = fIi(i, n)
                Ij = fIi(j, n)
                Ija = fIiast(j, n)
                
                Rf[i,j] = h.T*Ii*S*Ij*h
                If[i,j] = h.T*Ii*S*Ija*h
                
                RGh[i,j] = h.T*(Ii*S*Ij + Ij*S*Ii)
                RGe[i,j] = kron((Ij*h).T, h.T*Ii)*ET
                IGh[i,j] = h.T*(Ii*S*Ija + Ija*S*Ii)
                IGe[i,j] = kron((Ija*h).T, h.T*Ii)*ET
                
        for i in range(n):
            for j in range(n):
                coh = (Rf[i,j]**2 + If[i,j]**2)/(Rf[i,i]*Rf[j,j])
                
                omega_fh = vstack((RGh[i,j], IGh[i,j], 
                                   RGh[i,i], RGh[j,j]))
                omega_fe = vstack((RGe[i,j], IGe[i,j], 
                                   RGe[i,i], RGe[j,j]))
                omega_f = (omega_fh*omega_h*omega_fh.T +  
                          omega_fe*omega_e*omega_fe.T) 
            
                cohG = array([2*Rf[i,j]/(Rf[i,i]*Rf[j,j]), 2*If[i,j]/(Rf[i,i]*Rf[j,j]),
                              -coh/Rf[i,i], -coh/Rf[j,j]])
                
                var1 = cohG*omega_f*cohG.T   
                ic1 = coh - sqrt(var1)*st.norm.ppf(1-pr_.alpha/2.0)
                ic2 = coh + sqrt(var1)*st.norm.ppf(1-pr_.alpha/2.0)
                
                omega_f = omega_f[:2,:2]
                
                cohG2 = array([[2/(Rf[i,i]*Rf[j,j]),0], 
                               [0,2/(Rf[i,i]*Rf[j,j])]])
                
                L = fChol(omega_f)        
                d = fEig(L, cohG2)
                #d1 = fEig(omega, G2)
                #d2 = fEig(omega_evar, d2coh_dev2)
                #d = concatenate((d1,d2)) #TODO: conferir
                patdf = sum(d)**2/sum(d**2)
                patden = sum(d)/sum(d**2)
                th = st.chi2.ppf(1-pr_.alpha, patdf)/(patden*2*nd)
                var2 = 2*patdf/(patden*2*nd)**2
                
                res[:,i,j,ff] = [coh, th, ic1, ic2, patdf, patden, var1, var2]
                
    res_.mes = res[0]
    res_.th = res[1]
    res_.ic1 = res[2]
    res_.ic2 = res[3]
    res_.asy = res[4:]
    
    return res




def asymp_ss(x):
    '''Asymptotic statistics for the Coherence (coh)
        x -> data
        res_.A -> autoregressive matrix
        res_.er -> residues
    '''
    nf = pr_.nf
    n,n,p = res_.A.shape
    n,nd = x.shape
    
    x = mat(x)
    er = mat(res_.er)
    Af = A_to_f(res_.A, nf)    
    
    res = empty([8, n, n, nf])
    
    gamma = mat(bigautocorr(x, p))
    omega_a = kron(gamma.I, er)
    
    D = Dup(n)
    omega_e = 2*D*D.I*kron(er, er)*D.I.T*D.T
    
    S = kron(I(2*n), er)
    Teta = fTe(n)
    #E = ?
    
    for ff in range(nf):
        #print 'ff', ff
        f = ff/(2.0*nf)
        Ca = fCa(f, p, n)
        
        Hf = mat(Af[ff, :, :]).I
        h = vec(Hf)
        h = cat(h.real, h.imag, 0)
        
        dhda = fdh_da(mat(Af[ff, :, :]), n)
        
        omega_h = dhda*Ca*omega_a*Ca.T*dhda.T
        
        Rf = empty([n,n])
        If = empty([n,n])
        RGh = empty([n,n,size(h)])
        RGe = empty([n,n,size(h)])
        IGh = empty([n,n,size(h)])
        IGe = empty([n,n,size(h)])
        for i in range(n):
            for j in range(n):
                Ii = fIi(i, n)
                Ij = fIi(j, n)
                Ija = fIiast(j, n)
                
                Rf[i,j] = h.T*Ii*S*Ij*h
                If[i,j] = h.T*Ii*S*Ija*h
                
                RGh[i,j] = h.T*(Ii*S*Ij + Ij*S*Ii)
                RGe[i,j] = kron((Ij*h).T, h.T*Ii)*E*Teta
                IGh[i,j] = h.T*(Ii*S*Ija + Ija*S*Ii)
                IGe[i,j] = kron((Ija*h).T, h.T*Ii)*E*Teta
                
        for i in range(n):
            for j in range(n):
                ss = Rf[i,j]**2 + If[i,j]**2
                
                omega_fh = vstack((RGh[i,j], IGh[i,j]))
                omega_fe = vstack((RGe[i,j], IGe[i,j]))
                omega_f = (omega_fh*omega_h*omega_fh.T +  
                          omega_fe*omega_e*omega_fe.T) 
            
                ssG = array([2*Rf[i,j], 2*If[i,j]])
                
                var1 = ssG*omega_f*ssG.T   
                ic1 = ss - sqrt(var1)*st.norm.ppf(1-pr_.alpha/2.0)
                ic2 = ss + sqrt(var1)*st.norm.ppf(1-pr_.alpha/2.0)
                
                ssG2 = array([[2,0], 
                               [0,2]])
                
                L = fChol(omega_f)        
                d = fEig(L, ssG2)
                #d1 = fEig(omega, G2)
                #d2 = fEig(omega_evar, d2coh_dev2)
                #d = concatenate((d1,d2)) #TODO: conferir
                patdf = sum(d)**2/sum(d**2)
                patden = sum(d)/sum(d**2)
                th = st.chi2.ppf(1-pr_.alpha, patdf)/(patden*2*nd)
                var2 = 2*patdf/(patden*2*nd)**2
                
                res[:,i,j,ff] = [ss, th, ic1, ic2, patdf, patden, var1, var2]
                
    res_.mes = res[0]
    res_.th = res[1]
    res_.ic1 = res[2]
    res_.ic2 = res[3]
    res_.asy = res[4:]
    
    return res




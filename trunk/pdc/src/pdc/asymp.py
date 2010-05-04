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
import pdc.params as gl
import time
import sys

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
    ''' Returns vech(v) '''
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

def fdebig_de(n):
    '''Derivative of kron(I(2n), A) by A'''
    return dot(kron(TT(2*n, n), I(n*2*n)),
           kron(I(n), kron(vec(I(2*n)), I(n))))
    

def fdebig_de_small(n):
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

def bigautocorr_old(x, p):
    '''Autocorrelation. Data in rows. From order 0 to p-1.
    Output: nxn blocks of autocorr of lags i. (Nuttall Strand matrix)'''
    y = x[:]
    nd = x.shape[1]
    for i in arange(1, p):
        y = concatenate((y, xlag(x, i)), axis=0)
    return dot(y, y.T)/nd
    #return cov(y.T)

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

def bigautocorr_old_wrong(x, p):
    '''Autocorrelation. Data in rows. From order 0 to p-1.
    Output: nxn blocks of autocorr of lags i. (Nuttall Strand matrix)
    
    its not positive definite... it's wrong.
    '''
    n, nd = x.shape
    gamma = zeros([n*p, n*p])
    for i in arange(p):
        g = dot(x, xlag(x, i).T)/(nd-i)
        for j in arange(p-i):
            gamma[(i+j)*n:(i+j)*n+n, j*n:j*n+n] = g.T
            gamma[j*n:j*n+n, (i+j)*n:(i+j)*n+n] = g
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

def asymp_pdc(x, A, nf, e_var, p, metric = 'gen', alpha = 0.05):
    '''Asymptotic statistics for the three PDC formulations
        x -> data
        A -> autoregressive matrix
        nf -> number of frequencies
        e_var -> residues
        p -> model order
        metric -> witch PDC (gPDC = 'gen', dPDC = 'diag', PDC = 'euc')
        alpha -> confidence margin
    '''
    
    x = mat(x)
    e_var = mat(e_var)
    Af = A_to_f(A, nf)
    
    n, nd = x.shape
    
    th = empty([n, n, nf])
    ic1 = empty([n, n, nf])
    ic2 = empty([n, n, nf])
    pdc = empty([n, n, nf])
    varass1 = empty([n, n, nf])
    varass2 = empty([n, n, nf])
    patdfr = empty([n, n, nf])
    patdenr = empty([n, n, nf])
    
    gamma = mat(bigautocorr(x, p))
    gammai = inv(gamma)
    omega = kron(gammai, e_var)
    
    omega_evar = 2*Dup(n).I*kron(e_var, e_var)*Dup(n).I.T
    
   

    #print 'alpha', alpha
    
    for ff in range(nf):
        
        #print 'ff', ff
        
        f = ff/(2.0*nf)
        
        Ca = fCa(f, p, n)
        
        omega2 = Ca*omega*Ca.T
        
        L = fChol(omega2)
            
        a = vec(Af[ff, :, :])
        a = cat(a.real, a.imag, 0)
        #a = vec(cat(I(n), O(n), 1)) - dot(Ca, al)
        
        for i in range(n):
            for j in range(n):
                
                Iij = fIij(i, j, n)
                Ij = fIj(j, n)
                
                'No caso diag ou gen, deve-se acrescentar o evar na formula'
                if metric == 'euc':
                    Iije = Iij
                    Ije = Ij
                
                elif metric == 'diag':
                    evar_d = mdiag(e_var)
                    evar_d_big = kron(I(2*n), evar_d)
                    Iije = Iij*evar_d_big.I
                    Ije = Ij*evar_d_big.I
                
                else: #metric == 'gen' 
                    evar_d = mdiag(e_var)
                    evar_d_big = kron(I(2*n), evar_d)
                    Iije = Iij*evar_d_big.I
                    
                    evar_big = kron(I(2*n), e_var)
                    Ije = Ij*evar_big.I*Ij                
                
                num = a.T*Iije*a
                den = a.T*Ije*a
                pdc[i, j, ff] = num/den
                
                'Acrescenta derivada em relacao a evar'
                if metric == 'euc':
                    dpdc_dev = mat(zeros((n*(n+1))/2))
                
                elif metric == 'diag':
                    
                    #todo: tirar partes que nao dependem de f do loop.
                    
                    if i == 0 and j == 0 and ff == 0:
                        evar_d = mdiag(e_var)
                        evar_d_big = kron(I(2*n), evar_d)
                        inv_ed = evar_d_big.I
                        
                        'derivada de vec(Ed-1) por vecE'
                        de_deh = Dup(n)
                        debig_de = fdebig_de(n)
                        dedinv_dev = diagtom(vec(-inv_ed*inv_ed))
                        dedinv_deh = dedinv_dev*debig_de*de_deh
                    
                    'derivada do num por vecE'
                    dnum_dev = kron((Iij*a).T, a.T)*dedinv_deh
                    'derivada do den por vecE'
                    dden_dev = kron((Ij*a).T, a.T)*dedinv_deh
                    dpdc_dev = (den*dnum_dev - num*dden_dev)/(den**2)
                
                else: # metric == 'gen'
                    
                    if i == 0 and j == 0 and ff == 0:
                        evar_d = mdiag(e_var)
                        evar_d_big = kron(I(2*n), evar_d)
                        inv_ed = evar_d_big.I
                        
                        evar_big = kron(I(2*n), e_var)
                        inv_e = evar_big.I
    
                        'derivada de vec(Ed-1) por vecE'
                        de_deh = Dup(n)
                        debig_de = fdebig_de(n)
                        
                        dedinv_devd = diagtom(vec(-inv_ed*inv_ed)) 
                        dedinv_dehd = dedinv_devd*debig_de*de_deh
                        
                        dedinv_dev = -kron(inv_e.T, inv_e)
                        dedinv_deh = dedinv_dev*debig_de*de_deh
                    
                    'derivada do num por vecE'
                    dnum_dev = kron((Iij*a).T, a.T)*dedinv_dehd
                    'derivada do den por vecE'
                    dden_dev = kron((Ij*a).T, a.T*Ij)*dedinv_deh
                    dpdc_dev = (den*dnum_dev - num*dden_dev)/(den**2)
                    
                G1a = 2*a.T*Iije/den - 2*num*a.T*Ije/(den**2)
                G1 = -G1a*Ca

                varalpha = G1*omega*G1.T
                varevar = dpdc_dev*omega_evar*dpdc_dev.T
                varass1[i, j, ff] = (varalpha + varevar)/nd
                
                ic1[i, j, ff] = pdc[i, j, ff] - sqrt(varass1[i, j, ff])*st.norm.ppf(1-alpha/2.0)
                ic2[i, j, ff] = pdc[i, j, ff] + sqrt(varass1[i, j, ff])*st.norm.ppf(1-alpha/2.0)
                
                G2a = 2*Iije/den
                #G2 = Ca.T*G2a*Ca
                
                #print omega, eigh(omega, eigvals_only=True)
                d = fEig(L, G2a)
                
                patdf = sum(d)**2/sum(d**2)
                patden = sum(d)/sum(d**2)
                th[i, j, ff] = st.chi2.ppf(1-alpha, patdf)/(patden*2*nd)
                varass2[i, j, ff] = 2*patdf/(patden*2*nd)**2
                patdfr[i, j, ff] = patdf
                patdenr[i, j, ff] = patden
                
                
                
                #if (i == 1 and j == 0 and ff == 3):
                #    patdfr = patdf
                #    patdenr = patden TODO retirar

    gl.res_.varass1 = varass1
    gl.res_.patden = patdenr
    gl.res_.patdf = patdfr

    return pdc, th, ic1, ic2


def asymp_dtf(x, A, nf, e_var, p, alpha = 0.05, metric = 'dummy'):
    '''Asymptotic statistics for the DTF
        x -> data
        A -> autoregressive matrix
        nf -> number of frequencies
        e_var -> residues
        p -> model order
        alpha -> confidence margin
    '''
    
    x = mat(x)
    e_var = mat(e_var)
    Af = A_to_f(A, nf)
    
    n, nd = x.shape
    
    th = empty([n, n, nf])
    ic1 = empty([n, n, nf])
    ic2 = empty([n, n, nf])
    dtf = empty([n, n, nf])
    varass1 = empty([n, n, nf])
    varass2 = empty([n, n, nf])
    patdfr = empty([n, n, nf])
    patdenr = empty([n, n, nf])
    
    gammai = inv(bigautocorr(x, p))
    omega = kron(gammai, e_var)
    
    for ff in range(nf):
        f = ff/(2.0*nf)
        
        Ca = fCa(f, p, n)
        
        Hf = mat(Af[ff, :, :]).I
        h = vec(Hf)
        h = cat(h.real, h.imag, 0)
        
        dhda = fdh_da(mat(Af[ff, :, :]), n)
        
        #L = fChol(omega)
        
        omega2 = dhda*Ca*omega*Ca.T*dhda.T
        L = fChol(omega2)
       
        for i in range(n):
            for j in range(n):
                
                Iij = fIij(i, j, n)
                
                Ii = fIi(i, n)
                
                num = h.T*Iij*h
                den = h.T*Ii*h
                dtf[i, j, ff] = num/den
    
                G1a = 2*h.T*Iij/den - 2*num*h.T*Ii/(den**2)
                G1 = -G1a*dhda*Ca
                
                varass1[i, j, ff] = G1*omega*G1.T/nd
                ic1[i, j, ff] = dtf[i, j, ff] - sqrt(varass1[i, j, ff])*st.norm.ppf(1-alpha/2)
                ic2[i, j, ff] = dtf[i, j, ff] + sqrt(varass1[i, j, ff])*st.norm.ppf(1-alpha/2)
                
                G2a = 2*Iij/den
                #G2 = Ca.T*dhda.T*G2a*dhda*Ca 
                G2 = G2a
                
                d = fEig(L, G2)
                patdf = sum(d)**2/sum(d**2)
                patden = sum(d)/sum(d**2)
                
                th[i, j, ff] = st.chi2.ppf(1-alpha, patdf)/(patden*2*nd)
                varass2[i, j, ff] = 2*patdf/(patden*2*nd)**2
                patdfr[i, j, ff] = patdf
                patdenr[i, j, ff] = patden
                
    gl.res_.varass1 = varass1
    gl.res_.patden = patdenr
    gl.res_.patdf = patdfr
    
    return dtf, th, ic1, ic2

def fc(i, n):
    vi = zeros(n)
    vi[i] = 1
    vi.resize(n,1)
    return mat(kron(vi, I(n)))

def fk1(e_var_inv, i,j,n):
    ci = fc(i,n)
    cj = fc(j,n)
    
    return kron(I(2), ci)*e_var_inv*kron(I(2),cj).T
    
def fk2(e_var_inv, i,j,n):
    ci = fc(i,n)
    cj = fc(j,n)
    
    return kron(I(2), ci)*e_var_inv*kron(array([[0,1],[-1,0]]),cj).T
    
    
def asymp_pc(x, A, nf, e_var, p, alpha = 0.05, metric = 'dummy'):
    '''Asymptotic statistics for the Partial Coherence (PC)
        x -> data
        A -> autoregressive matrix
        nf -> number of frequencies
        e_var -> residues
        p -> model order
        alpha -> confidence margin
    '''
       
    x = mat(x)
    e_var = mat(e_var)
    Af = A_to_f(A, nf)
    
    n, nd = x.shape
    
    th = empty([n, n, nf])
    ic1 = empty([n, n, nf])
    ic2 = empty([n, n, nf])
    pc = empty([n, n, nf])
    varass1 = empty([n, n, nf]) #TODO: retirar varass's
    varass2 = empty([n, n, nf])
    patdfr = empty([n, n, nf])
    patdenr = empty([n, n, nf])
    
    gammai = inv(bigautocorr(x, p))
    omega = kron(gammai, e_var)
    #print bigautocorr(x, p)
    
    omega_evar = 2*Dup(n).I*kron(e_var, e_var)*Dup(n).I.T
    
    for ff in range(nf):
        f = ff/(2.0*nf)
        
        Ca = fCa(f, p, n)
        
        
        omega2 = Ca*omega*Ca.T
        
        ehs = omega_evar.shape[0]
        oas = omega2.shape[0]
        omegabig = cat(cat(omega2, zeros([oas, ehs]), 1),
                      cat(zeros([ehs, oas]), omega_evar, 1), 0) 

        L = fChol(omegabig)
        
        a = vec(Af[ff, :, :])
        a = cat(a.real, a.imag, 0)
        
        for i in range(n):
            for j in range(n):
                
                evar_big = kron(I(2), e_var)
                evar_big_inv = evar_big.I
                ei = evar_big_inv
                
                k1ij = fk1(evar_big_inv, i, j, n)       
                k2ij = fk2(evar_big_inv, i, j, n)     
                k1ii = fk1(evar_big_inv, i, i, n)      
                k1jj = fk1(evar_big_inv, j, j, n)  
                
                num = (a.T*k1ij*a)**2+(a.T*k2ij*a)**2
                den = (a.T*k1ii*a)*(a.T*k1jj*a)
                pc[i, j, ff] = num/den
                
                #Acrescenta derivada em relacao a evar
                ci = fc(i,n)
                cj = fc(j,n)
                
                a11i = a.T*kron(I(2), ci)
                a12j = kron(I(2),cj).T*a
                a21i = a.T*kron(I(2), ci)
                a22j = kron(array([[0,1],[-1,0]]),cj).T*a
                
                a11j = a.T*kron(I(2), cj)
                a12i = kron(I(2),ci).T*a
                
                num1 = a11i*ei*a12j
                num2 = a21i*ei*a22j
                den1 = a11i*ei*a12i
                den2 = a11j*ei*a12j
                
                dnum1_dei = kron(a12j.T, a11i)
                dnum2_dei = kron(a22j.T, a21i)
                dden1_dei = kron(a12i.T, a11i)
                dden2_dei = kron(a12j.T, a11j)
                
                #check den1*den2 == den e num1**2+num2**2 = num. checked!.
                dpc_dei = (2*num1*dnum1_dei + 2*num2*dnum2_dei)/den + \
                          - num*(dden1_dei/den1 + dden2_dei/den2)/den
                          
                'derivada de vec(Ed-1) por vecE'
                de_deh = Dup(n)
                debig_de = fdebig_de_small(n)
                dedinv_dev = -kron(ei.T, ei)
                dedinv_deh = dedinv_dev*debig_de*de_deh
                
                dpc_dev = dpc_dei*dedinv_deh
                
                d2pc_dei2 = 2*(dnum1_dei.T*dnum1_dei +  \
                               dnum2_dei.T*dnum2_dei)/den
                               
                d2pc_dev2 = dedinv_deh.T*d2pc_dei2*dedinv_deh
                
                #novas derivadas:
                dnum1_deh = dnum1_dei*dedinv_deh
                dnum2_deh = dnum2_dei*dedinv_deh
                #dnum1_da = a.T*(k1ij+k1ij.T)*Ca
                #dnum2_da = a.T*(k2ij+k2ij.T)*Ca
                dnum1_da = a.T*(k1ij+k1ij.T)
                dnum2_da = a.T*(k2ij+k2ij.T)
                d2pc_dade = 2*(dnum1_da.T*dnum1_deh + dnum2_da.T*dnum2_deh)/den
                d2pc_deda = d2pc_dade.T
                
                
                #fim da derivada por e_var
                
                dnum_da = 2*((a.T*k1ij*a)*(a.T*(k1ij+k1ij.T)) + \
                             (a.T*k2ij*a)*(a.T*(k2ij+k2ij.T)))
                dden_da = 2*((a.T*k1jj*a)*(a.T*k1ii) + \
                             (a.T*k1ii*a)*(a.T*k1jj))
                G1a = dnum_da/den - num*dden_da/den**2
                G1 = -G1a*Ca
                
            
                varalpha = G1*omega*G1.T 
                varevar = dpc_dev*omega_evar*dpc_dev.T
                
                varass1[i, j, ff] = (varalpha + varevar)/nd
                
                ic1[i, j, ff] = pc[i, j, ff] - sqrt(varass1[i, j, ff])*st.norm.ppf(1-alpha/2.0)
                ic2[i, j, ff] = pc[i, j, ff] + sqrt(varass1[i, j, ff])*st.norm.ppf(1-alpha/2.0)
                
                d1 = (k1ij+k1ij.T)*a
                d2 = (k2ij+k2ij.T)*a
                G2a = 2*(d1*d1.T + d2*d2.T)/den
                #G2 = Ca.T*G2a*Ca
                G2 = G2a 
                
                Gbig = cat(cat(G2, d2pc_dade, 1),
                           cat(d2pc_deda, d2pc_dev2, 1), 0) 
                
                d = fEig(L, Gbig)
                
                #d1 = fEig(omega, G2)
                #d2 = fEig(omega_evar, d2pc_dev2)
                #d = concatenate((d1,d2)) #TODO: conferir
                #print 'd', ff, i, j, d
                patdf = sum(d)**2/sum(d**2)
                patden = sum(d)/sum(d**2)
                th[i, j, ff] = st.chi2.ppf(1-alpha, patdf)/(patden*2*nd)
                varass2[i, j, ff] = 2*patdf/(patden*2*nd)**2
                patdfr[i, j, ff] = patdf
                patdenr[i, j, ff] = patden
                
    gl.res_.varass1 = varass1
    gl.res_.patden = patdenr
    gl.res_.patdf = patdfr
                
    return pc, th, ic1, ic2

def fl(i, n):
    vi = zeros(n)
    vi[i] = 1
    vi.resize(n,1)
    return mat(kron(I(n), vi))

def fkl1(evar_big, i,j,n):
    li = fl(i,n)
    lj = fl(j,n)
    
    return kron(I(2), li)*evar_big*kron(I(2),lj).T

def fkl2(evar_big, i,j,n):
    li = fl(i,n)
    lj = fl(j,n)
    
    return kron(I(2), li)*evar_big*kron(array([[0,1],[-1,0]]),lj).T

def asymp_coh(x, A, nf, e_var, p, alpha = 0.05, metric = 'dummy'):
    '''Asymptotic statistics for the Coherence (coh)
        x -> data
        A -> autoregressive matrix
        nf -> number of frequencies
        e_var -> residues
        p -> model order
        alpha -> confidence margin
    '''
       
    x = mat(x)
    e_var = mat(e_var)
    Af = A_to_f(A, nf)
    
    n, nd = x.shape
    
    th = empty([n, n, nf])
    ic1 = empty([n, n, nf])
    ic2 = empty([n, n, nf])
    coh = empty([n, n, nf])
    varass1 = empty([n, n, nf])
    varass2 = empty([n, n, nf])
    patdfr = empty([n, n, nf])
    patdenr = empty([n, n, nf])
    
    gammai = inv(bigautocorr(x, p))
    omega = kron(gammai, e_var)
    #print bigautocorr(x, p)
    
    omega_evar = 2*Dup(n).I*kron(e_var, e_var)*Dup(n).I.T
    
    
    for ff in range(nf):
        #print 'ff', ff
        
        f = ff/(2.0*nf)
        
        Ca = fCa(f, p, n)
        
        Hf = mat(Af[ff, :, :]).I
        h = vec(Hf)
        h = cat(h.real, h.imag, 0)
        
        dhda = fdh_da(mat(Af[ff, :, :]), n)
        
        omega2 = dhda*Ca*omega*Ca.T*dhda.T
        
        ehs = omega_evar.shape[0]
        oas = omega2.shape[0]
        omegabig = cat(cat(omega2, zeros([oas, ehs]), 1),
                   cat(zeros([ehs, oas]), omega_evar, 1), 0)
         
        L = fChol(omegabig)
        
        for i in range(n):
            for j in range(n):
                
                evar_big = kron(I(2), e_var)
                #evar_inv = e_var.I
                #ei = evar_inv
                
                k1ij = fkl1(evar_big, i, j, n)       
                k2ij = fkl2(evar_big, i, j, n)     
                k1ii = fkl1(evar_big, i, i, n)      
                k1jj = fkl1(evar_big, j, j, n)  
                
                num = (h.T*k1ij*h)**2+(h.T*k2ij*h)**2
                den = (h.T*k1ii*h)*(h.T*k1jj*h)
                coh[i, j, ff] = num/den
                
                #Acrescenta derivada em relacao a evar
                li = fl(i,n)
                lj = fl(j,n)
                
                h11i = h.T*kron(I(2), li)
                h12j = kron(I(2),lj).T*h
                h21i = h.T*kron(I(2), li)
                h22j = kron(array([[0,1],[-1,0]]),lj).T*h
                
                h11j = h.T*kron(I(2), lj)
                h12i = kron(I(2),li).T*h
                
                num1 = h11i*evar_big*h12j
                num2 = h21i*evar_big*h22j
                den1 = h11i*evar_big*h12i
                den2 = h11j*evar_big*h12j
                
                dnum1_de = kron(h12j.T, h11i)
                dnum2_de = kron(h22j.T, h21i)
                dden1_de = kron(h12i.T, h11i)
                dden2_de = kron(h12j.T, h11j)
                
                #check den1*den2 == den e num1**2+num2**2 = num. checked!.
                dcoh_de = (2*num1*dnum1_de + 2*num2*dnum2_de)/den + \
                          - num*(dden1_de/den1 + dden2_de/den2)/den
                          
                'derivada de Ebig por vecE'
                de_deh = Dup(n)
                debig_de = fdebig_de_small(n)
                ded_deh = debig_de*de_deh
                
                dcoh_dev = dcoh_de*ded_deh
                
                d2coh_de2 = 2*(dnum1_de.T*dnum1_de +  \
                               dnum2_de.T*dnum2_de)/den
                               
                d2coh_dev2 = ded_deh.T*d2coh_de2*ded_deh
                
                #novas derivadas:
                dnum1_deh = dnum1_de*ded_deh
                dnum2_deh = dnum2_de*ded_deh
                dnum1_dh = h.T*(k1ij+k1ij.T)
                dnum2_dh = h.T*(k2ij+k2ij.T)
                #dnum1_dh = h.T*(k1ij+k1ij.T)*dhda*Ca
                #dnum2_dh = h.T*(k2ij+k2ij.T)*dhda*Ca
                d2coh_dhde = 2*(dnum1_dh.T*dnum1_deh + dnum2_dh.T*dnum2_deh)/den
                d2coh_dedh = d2coh_dhde.T
                
                #fim da derivada por e_var
                
                dnum_dh = 2*((h.T*k1ij*h)*(h.T*(k1ij+k1ij.T)) + \
                             (h.T*k2ij*h)*(h.T*(k2ij+k2ij.T)))
                dden_dh = 2*((h.T*k1jj*h)*(h.T*k1ii) + \
                             (h.T*k1ii*h)*(h.T*k1jj))
                G1a = dnum_dh/den - num*dden_dh/den**2
                G1 = -G1a*dhda*Ca
                
            
                varalpha = G1*omega*G1.T 
                varevar = dcoh_dev*omega_evar*dcoh_dev.T
                
                #print varalpha, varevar
                varass1[i, j, ff] = (varalpha + varevar)/nd
                
                
                ic1[i, j, ff] = coh[i, j, ff] - sqrt(varass1[i, j, ff])*st.norm.ppf(1-alpha/2.0)
                ic2[i, j, ff] = coh[i, j, ff] + sqrt(varass1[i, j, ff])*st.norm.ppf(1-alpha/2.0)
                
                d1 = (k1ij+k1ij.T)*h
                d2 = (k2ij+k2ij.T)*h
                G2a = 2*(d1*d1.T + d2*d2.T)/den
                #G2 = Ca.T*dhda.T*G2a*dhda*Ca 
                G2 = G2a
                
                Gbig = cat(cat(G2, d2coh_dhde, 1),
                           cat(d2coh_dedh, d2coh_dev2, 1), 0) 
               
                
                d = fEig(L, Gbig)
                
                #d1 = fEig(omega, G2)
                #d2 = fEig(omega_evar, d2coh_dev2)
                #d = concatenate((d1,d2)) #TODO: conferir
                patdf = sum(d)**2/sum(d**2)
                patden = sum(d)/sum(d**2)
                th[i, j, ff] = st.chi2.ppf(1-alpha, patdf)/(patden*2*nd)
                varass2[i, j, ff] = 2*patdf/(patden*2*nd)**2
                patdfr[i, j, ff] = patdf
                patdenr[i, j, ff] = patden
                
    gl.res_.varass1 = varass1
    gl.res_.patden = patdenr
    gl.res_.patdf = patdfr
                
    return coh, th, ic1, ic2

def asymp_ss(x, A, nf, e_var, p, alpha = 0.05, metric = 'dummy'):
    '''Asymptotic statistics for the Spectral density (SS)
        x -> data
        A -> autoregressive matrix
        nf -> number of frequencies
        e_var -> residues
        p -> model order
        alpha -> confidence margin
    '''
       
    x = mat(x)
    e_var = mat(e_var)
    Af = A_to_f(A, nf)
    
    n, nd = x.shape
    
    th = empty([n, n, nf])
    ic1 = empty([n, n, nf])
    ic2 = empty([n, n, nf])
    ss = empty([n, n, nf])
    varass1 = empty([n, n, nf])
    varass2 = empty([n, n, nf])
    
    gammai = inv(bigautocorr(x, p))
    omega = kron(gammai, e_var)
    #print bigautocorr(x, p)
    omega_evar = 2*Dup(n).I*kron(e_var, e_var)*Dup(n).I.T
    
    for ff in range(nf):
        f = ff/(2.0*nf)
        
        Ca = fCa(f, p, n)
        
        Hf = mat(Af[ff, :, :]).I
        h = vec(Hf)
        h = cat(h.real, h.imag, 0)
        
        dhda = fdh_da(mat(Af[ff, :, :]), n)
        
        omega2 = dhda*Ca*omega*Ca.T*dhda.T
        
        ehs = omega_evar.shape[0]
        oas = omega2.shape[0]
        omegabig = cat(cat(omega2, zeros([oas, ehs]), 1),
                   cat(zeros([ehs, oas]), omega_evar, 1), 0)
    
        L = fChol(omegabig) 
    
        
        for i in range(n):
            for j in range(n):
                
                evar_big = kron(I(2), e_var)
                #evar_big_inv = evar_big.I
                #ei = evar_big_inv
                #evar_inv = e_var.I
                #ei = evar_inv
                
                k1ij = fkl1(evar_big, i, j, n)       
                k2ij = fkl2(evar_big, i, j, n)
                
                num = (h.T*k1ij*h)**2+(h.T*k2ij*h)**2
                ss[i, j, ff] = num
                
                #Acrescenta derivada em relacao a evar
                li = fl(i,n)
                lj = fl(j,n)
                
                h11i = h.T*kron(I(2), li)
                h12j = kron(I(2),lj).T*h
                h21i = h.T*kron(I(2), li)
                h22j = kron(array([[0,1],[-1,0]]),lj).T*h
                
                num1 = h11i*evar_big*h12j
                num2 = h21i*evar_big*h22j
                
                dnum1_de = kron(h12j.T, h11i)
                dnum2_de = kron(h22j.T, h21i)
                
                #check den1*den2 == den e num1**2+num2**2 = num. checked!.
                dss_de = 2*num1*dnum1_de + 2*num2*dnum2_de
                          
                'derivada de vec(Ed-1) por vecE'
                de_deh = Dup(n)
                debig_de = fdebig_de_small(n)
                #dedinv_dev = -kron(ei.T, ei)
                ded_deh = debig_de*de_deh
                #dedinv_deh = dedinv_dev*de_deh
                
                dss_dev = dss_de*ded_deh
                
                d2ss_de2 = 2*(dnum1_de.T*dnum1_de +  \
                               dnum2_de.T*dnum2_de)
                               
                d2ss_dev2 = ded_deh.T*d2ss_de2*ded_deh
                
                dnum1_deh = dnum1_de*ded_deh
                dnum2_deh = dnum2_de*ded_deh
                dnum1_dh = h.T*(k1ij+k1ij.T)
                dnum2_dh = h.T*(k2ij+k2ij.T)
                #dnum1_dh = h.T*(k1ij+k1ij.T)*dhda*Ca
                #dnum2_dh = h.T*(k2ij+k2ij.T)*dhda*Ca
                d2ss_dhde = 2*(dnum1_dh.T*dnum1_deh + dnum2_dh.T*dnum2_deh)
                d2ss_dedh = d2ss_dhde.T
                
                #fim da derivada por e_var
                
                dnum_dh = 2*((h.T*k1ij*h)*(h.T*(k1ij+k1ij.T)) + \
                             (h.T*k2ij*h)*(h.T*(k2ij+k2ij.T)))
                G1a = dnum_dh
                G1 = -G1a*dhda*Ca
                
            
                varalpha = G1*omega*G1.T 
                varevar = dss_dev*omega_evar*dss_dev.T
                varass1[i, j, ff] = (varalpha + varevar)/nd
                
                
                ic1[i, j, ff] = ss[i, j, ff] - sqrt(varass1[i, j, ff])*st.norm.ppf(1-alpha/2.0)
                ic2[i, j, ff] = ss[i, j, ff] + sqrt(varass1[i, j, ff])*st.norm.ppf(1-alpha/2.0)
                
                d1 = (k1ij+k1ij.T)*h
                d2 = (k2ij+k2ij.T)*h
                G2a = 2*(d1*d1.T + d2*d2.T)
                #G2 = Ca.T*dhda.T*G2a*dhda*Ca
                G2 = G2a 
                
                Gbig = cat(cat(G2, d2ss_dhde, 1),
                           cat(d2ss_dedh, d2ss_dev2, 1), 0) 
                d = fEig(L, Gbig)
                
                #d1 = fEig(omega, G2)
                #d2 = fEig(omega_evar, d2ss_dev2)
                #d = concatenate((d1,d2)) #TODO: conferir

                patdf = sum(d)**2/sum(d**2)
                patden = sum(d)/sum(d**2)
                th[i, j, ff] = st.chi2.ppf(1-alpha, patdf)/(patden*2*nd)
                varass2[i, j, ff] = 2*patdf/(patden*2*nd)**2
                patdfr[i, j, ff] = patdf
                patdenr[i, j, ff] = patden
                
    gl.res_.varass1 = varass1
    gl.res_.patden = patdenr
    gl.res_.patdf = patdfr
                
    return ss, th, ic1, ic2

def fCij(i, j, n, p):
    '''Returns Cij of the formula'''
    Cij = mat(zeros([1,n*n]))
    Cij[0,(j*n + i)] = 1
    Cij = kron(I(p), Cij)
    return mat(Cij)

def asymp_gct(x, A, e_var):
    '''Asymptotic statistics for Wald statistic of the GC in time
        x -> data
        A
        e_var
    '''
    
    x = mat(x)
    e_var = mat(e_var)
    
    n, nd = x.shape
    n, n, p = A.shape
    
    wt = empty([n, n])
    
    gammai = inv(bigautocorr(x, p))
    omega = kron(gammai, e_var)
    
    a = vec(A)
    
    for i in range(n):
        for j in range(n):
            Cij = fCij(i,j,n,p)
            wt[i,j] = (nd*1.0/p)*(Cij*a).T*((Cij*omega*Cij.T).I)*(Cij*a)
            
    pv = 1-st.f.cdf(wt, p, nd-n*p-1)
    
    return pv, wt

def asymp_igct(e_var, nd):
    '''Asymptotic statistics for Wald statistic of instantaneous GC
        e_var
        nd
    '''
    
    e_var = mat(e_var)
    
    n, n = e_var.shape
    
    wt = zeros([n, n])
    
    for i in arange(n):
        for j in arange(n):
            Cij = zeros([n,n])
            Cij[i,j] = 1
            Cij[j,i] = 1
            Cij = mat(vech(Cij))
            eij = e_var[i,j]
            Di = Dup(n).I
            wt[i,j] = nd*eij*(2*Cij*Di*kron(e_var,e_var)*Di.T*Cij.T).I*eij
            
            
    pv = 1-st.chi2.cdf(wt, 1)
    
    return pv, wt


def asymp_white(x, res, p, h = 20):

    n,nd = res.shape
    
    x = empty([n,n,h+1])
    
    for i in arange(n):
        for j in arange(n):
            dum, aux, dum, dum = xcorr(res[i], res[j], maxlags = h, normed = False)
            x[i,j] = aux[h:]
    
    s0 = mat(inv(x[:,:,0]))
    s = 0
    for i in arange(1,h+1):
        s = s + (1.0/(nd-i))*trace(x[:,:,i].T*s0*x[:,:,i]*s0)
    s = s*nd**2
    print s
    print (h-p)*n**2
    pv = 1-st.chi2.cdf(s,(h-p)*n**2)
    
    return pv, s


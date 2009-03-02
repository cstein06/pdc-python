from numpy import *
import scipy.stats as st
from scipy.stats import cov
from scipy.linalg import cholesky
from scipy.linalg import eigh
from scipy.linalg import inv

from algorithms.assym import bigautocorr

#corrlag = lambda a,b,lag: cov(a[k:],b[:-k])/a.size

vec = lambda x: x.ravel('F')
O = lambda n: zeros([n,n],  dtype=float)
I = lambda n: identity(n, dtype=float)
cat = lambda a,b, ax: concatenate((a,b), axis = ax)


def dtf_one_alg(A, er, nf = 64):
    '''Generates spectral non-norm. DTF matrix from AR matrix
    
      Input: 
        A(n, n, r) - recurrence matrix (n - number of signals, r - model order)
        e_cov(n, n) - error covariance matrix
        nf - frequency resolution
    
      Output:
        DTF(n, n, nf) - PDC matrix
    '''
    
    n, n, r = A.shape
    nor = ones(n) # dtf_one nao tem normalizacao
    
    AL = A_to_f(A, nf)
    HL = empty(AL.shape)
    for i in range(nf):
        HL[i] = inv(HL[i])
    
    # normalizacao por sum(ai ai* sig)
    dDTF = dot(nor,HL*HL.conj())
    nDTF = HL*sqrt(nor).reshape(-1,1)
    DTF = nDTF/sqrt(abs(dDTF)).reshape(nf,1,n).repeat(n, axis = 1)
    
    return DTF.transpose(1,2,0)

def dh_da(a, n):
    
    aa1 = cat(Af.real, -Af.imag, 1)
    aa2 = cat(Af.imag, Af.real, 1)
    aa = cat(aa1, aa2, 0)
    ha = inv(aa)
    dvechda = -kron(ha.T,ha)
    dhda = dvechda[:2*n**2,:2*n**2]
    
    tR1 = kron(I(n),kron(array([[1],[0]]), I(n)))
    tR2 = kron(I(n),kron(array([[0],[1]]), I(n)))
    tR = cat(tR1, tR2, 1)
    tS = tR.T
    
    dhda = dot(dot(dhda, tR), tS)
    
    return dhda

def assym_dtf_one(x, Af, e_var, p, alpha = 0.05):
    
    n, nd = x.shape
    nf = Af.shape[0]
    
    th = empty([n,n,nf])
    ic1 = empty([n,n,nf])
    ic2 = empty([n,n,nf])
    
    gammai = inv(bigautocorr(x, p))
    omega = kron(gammai,e_var)
    
    for ff in range(nf):
        f = ff/(2.0*nf)
        C1 = cos(-2*pi*f*arange(1,p+1))
        S1 = sin(-2*pi*f*arange(1,p+1))
        
        C2 = cat(C1.reshape(1,-1),S1.reshape(1,-1),0)
        
        Ca = kron(C2,identity(n**2))
        
        a = vec(Af[ff,:,:])
        a = cat(a.real, a.imag, 0)
        #a = vec(cat(I(n), O(n), 1)) - dot(Ca, al)
        
        dhda = dh_da(a, n)
        
        for i in range(n):
            for j in range(n):
                
                Iij = zeros(n**2)
                Iij[n*j+i] = 1
                Iij = diag(Iij)
                Iij = kron(I(2), Iij)
                
                Ii = zeros(n)
                Ii[i] = 1
                Ii = diag(Ii)
                Ii = kron(I(n), Ii)
                Ii = kron(I(2), Ii)
                
                num = dot(dot(h, Iij), h)
                den = dot(dot(h, Ii), h)
                pdc = num/den
    
                G1a = 2*dot(Iij,h)/den - 2*dot(Ii,h)*num/(den**2)
                G2a = 2*Iij/den
                
                G1 = dot(dot(G1a,dhda),-Ca)
                G2 = dot(dot(Ca.T,dot(dot(dhda.T,G2a),dhda)),Ca)
                
                varass1 = dot(dot(G1,omega),G1)/nd
                ic1[i,j,ff] = pdc - sqrt(varass1)*st.norm.ppf(1-alpha)
                ic2[i,j,ff] = pdc + sqrt(varass1)*st.norm.ppf(1-alpha)
                
                L = cholesky(omega, lower=1)
                D = dot(dot(L.T, G2), L)
                d = eigh(D, eigvals_only=True)
                d = d[abs(d) > 1E-8]
                
                #patdf = 0.5*sum(d)**2/(2*sum(d**2))
                #patden = sum(d)/(2*sum(d**2))
                patdf = sum(d)**2/sum(d**2)
                patden = sum(d)/sum(d**2)
                #print 1-alpha, patdf, patden*n, d
                th[i,j,ff] = st.chi2.ppf(1-alpha, patdf)/(patden*nd)
                
    return th, ic1, ic2
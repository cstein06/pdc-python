from numpy import *
import scipy.stats as st
from scipy.stats import cov as cov
from scipy.linalg import cholesky
from scipy.linalg import eigh
from scipy.linalg import inv as inv

#corrlag = lambda a,b,lag: cov(a[k:],b[:-k])/a.size

vec = lambda x: x.ravel('F')
O = lambda n: zeros([n,n],  dtype=float)
I = lambda n: identity(n, dtype=float)
cat = lambda a,b, ax: concatenate((a,b), axis = ax)

def xlag(x, lag):
    if(lag == 0):
        return x.copy;
    xl = zeros(x.shape)
    xl[:,:-lag] = x[:,lag:]
    return xl

# Autocorrelation. Data in rows. From order 0 to p-1.
# Output: nxn blocks of autocorr of lags i. (Nuttall Strand matrix)
def bigautocorr(x, p):
    y = x[:]
    for i in arange(1,p):
        y = concatenate((y, xlag(y,i)))
    return cov(y.T)

def assym_pdc_one(x, Af, e_var, p, alpha = 0.05):
    
    n, nd = x.shape
    nf = Af.shape[0]
    
    th = empty([n,n,nf])
    ic1 = empty([n,n,nf])
    ic2 = empty([n,n,nf])
    
    #print bigautocorr(x, p)
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
        
        for i in range(n):
            for j in range(n):
                
                Iij = zeros(n**2)
                Iij[n*j+i] = 1
                Iij = diag(Iij)
                Iij = kron(I(2), Iij)
                
                Ij = zeros(n)
                Ij[j] = 1
                Ij = diag(Ij)
                Ij = kron(Ij, I(n))
                Ij = kron(I(2), Ij)
                
                num = dot(dot(a, Iij), a)
                den = dot(dot(a, Ij), a)
                pdc = num/den
    
                G1a = 2*dot(Iij,a)/den - 2*dot(Ij,a)*num/(den**2)
                G2a = 2*Iij/den
                
                G1 = -dot(G1a,Ca)
                G2 = dot(dot(Ca.T,G2a),Ca)
                
                varass1 = dot(dot(G1,omega),G1)/nd
                ic1[i,j,ff] = pdc - sqrt(varass1)*st.norm.ppf(1-alpha)
                ic2[i,j,ff] = pdc + sqrt(varass1)*st.norm.ppf(1-alpha)
                
                L = cholesky(omega, lower=1)
                D = dot(dot(L.T, G2), L)
                d = eigh(D, eigvals_only=True)
                d = d[abs(d) > 1E-8]
                if (d.size > 2):
                    print 'assym tem mais de 2 chi2'
                
                #patdf = 0.5*sum(d)**2/(2*sum(d**2))
                #patden = sum(d)/(2*sum(d**2))
                patdf = sum(d)**2/sum(d**2)
                patden = sum(d)/sum(d**2)
                #print 1-alpha, patdf, patden*n, d
                th[i,j,ff] = st.chi2.ppf(1-alpha, patdf)/(patden*nd)
                
    return th, ic1, ic2
            
    

from numpy import *
import scipy.stats as st
from scipy.stats import cov as cov
from scipy.linalg import cholesky
from scipy.linalg import eigh
from scipy.linalg import inv as inv

#corrlag = lambda a,b,lag: cov(a[k:],b[:-k])/a.size

vec = lambda x: mat(x.ravel('F')).T
O = lambda n: mat(zeros([n, n], dtype=float))
I = lambda n: mat(identity(n, dtype=float))
cat = lambda a, b, ax: concatenate((a, b), axis = ax)
mdiag = lambda a: mat(diag(diag(a)))
diagtom = lambda a: mat(diag(array(a).reshape(-1)))

#GInv = lambda a: dot(inv(dot(a.T, a)),a.T)

def Dup(n):
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

def TT(n):
    ''' TT(n)*vec(B) = T(n,n)*vecB = vec(B.T)'''
    t = O(n**2)
    for i in range(n):
        for j in range(n):
            t[i*n+j, j*n+i] = 1 
    return t

def xlag(x, lag):
    if(lag == 0):
        return x.copy();
    xl = matrix(zeros(x.shape))
    xl[:, lag:] = x[:, :-lag]
    return xl

# Autocorrelation. Data in rows. From order 0 to p-1.
# Output: nxn blocks of autocorr of lags i. (Nuttall Strand matrix)
def bigautocorr(x, p):
    y = x[:]
    nd = x.shape[1]
    for i in arange(1, p):
        y = concatenate((y, xlag(y, i)), axis=0)
    return dot(y, y.T)/nd
    #return cov(y.T)
    
    
def fdh_da_old(Af, n):
    
    aa1 = cat(Af.real, -Af.imag, 1)
    aa2 = cat(Af.imag, Af.real, 1)
    aa = cat(aa1, aa2, 0)
    print 'aa', aa
    ha = aa.I
    print hae
    dvechda = -kron(ha.T, ha)
    dhdaa = dvechda[:2*n**2, :2*n**2]
    print dvechda
    print dhda
    
    tR1 = kron(I(n), kron(array([[1], [0]]), I(n)))
    tR2 = kron(I(n), kron(array([[0], [1]]), I(n)))
    tR = mat(cat(tR1, tR2, 1))
    
    dhda = tR*dhda*tR
    
    return dhda

    
def fdh_da(Af, n):
    
    aa1 = cat(Af.real, -Af.imag, 1)
    aa2 = cat(Af.imag, Af.real, 1)
    aa = cat(aa1, aa2, 0)
    #print 'aa', aa
    ha = aa.I
    #print ha
    dvechda = -kron(ha.T, ha)
    dhdaa = dvechda[:2*n**2, :]
    #print dvechda
    #print dhdaa
    
    tR1 = kron(I(n), kron(array([[1], [0]]), I(n)))
    tR2 = kron(I(n), kron(array([[0], [1]]), I(n)))
    tR = mat(cat(tR1, tR2, 1))
    
    tS = mat(kron(array([[0, -1], [1, 0]]), I(n**2)))
    
    daada = cat(tR, dot(tR, tS), 0)
    
    dhda = tR*dhdaa*daada
    
    return dhda
    
def fIij(i, j, n):
    Iij = zeros(n**2)
    Iij[n*j+i] = 1
    Iij = diag(Iij)
    return mat(kron(I(2), Iij))

def fIj(j, n):
    Ij = zeros(n)
    Ij[j] = 1
    Ij = diag(Ij)
    Ij = kron(Ij, I(n))
    return mat(kron(I(2), Ij))

def fIi(i, n):
    Ii = zeros(n)
    Ii[i] = 1
    Ii = diag(Ii)
    Ii = kron(I(n), Ii)
    return mat(kron(I(2), Ii))

def fCa(f, p, n):
    C1 = cos(-2*pi*f*mat(arange(1, p+1)))
    S1 = sin(-2*pi*f*mat(arange(1, p+1)))    
    C2 = cat(C1.reshape(1, -1), S1.reshape(1, -1), 0)
    return kron(C2, identity(n**2))    

def fEig(omega, G2):
    L = mat(cholesky(omega, lower=1))
    D = L.T*G2*L
    d = eigh(D, eigvals_only=True)
    d = d[abs(d) > 1E-8]
    if (d.size > 2):
        print 'assym tem mais de 2 chi2:'
        print d
    return d

def assym_pdc(x, Af, e_var, p, metric = 'gen', alpha = 0.05):
    
    x = mat(x)
    e_var = mat(e_var)
    
    n, nd = x.shape
    nf = Af.shape[0]
    
    th = empty([n, n, nf])
    ic1 = empty([n, n, nf])
    ic2 = empty([n, n, nf])
    varass1 = empty([n, n, nf])
    varass2 = empty([n, n, nf])
    
    gammai = inv(bigautocorr(x, p))
    omega = kron(gammai, e_var)
    #print bigautocorr(x, p)
    
    for ff in range(nf):
        f = ff/(2.0*nf)
        
        Ca = fCa(f, p, n)
        
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
                pdc = num/den
                
                
                'Acrescenta derivada em relacao a evar'
                if metric == 'euc':
                    dpdc_dev = mat(zeros((n*(n+1))/2))
                
                elif metric == 'diag':
                    evar_d = mdiag(e_var)
                    evar_d_big = kron(I(2*n), evar_d)
                    inv_ed = evar_d_big.I
                    
                    'derivada de vec(Ed-1) por vecE'
                    de_deh = Dup(n)
                    debig_de = kron(TT(n), I(4*n**2)) * \
                               kron(I(n), kron(vec(I(2*n)), I(n)))
                    dedinv_dev = diagtom(vec(-inv_ed*inv_ed))
                    dedinv_deh = dedinv_dev*debig_de*de_deh
                    
                    'derivada do num por vecE'
                    dnum_dev = kron((Iij*a).T, a.T)*dedinv_deh
                    'derivada do den por vecE'
                    dden_dev = kron((Ij*a).T, a.T)*dedinv_deh
                    dpdc_dev = (den*dnum_dev - num*dden_dev)/(den**2)
                
                else: # metric == 'gen'
                    evar_d = mdiag(e_var)
                    evar_d_big = kron(I(2*n), evar_d)
                    inv_ed = evar_d_big.I
                    
                    evar_big = kron(I(2*n), e_var)
                    inv_e = evar_big.I

                    'derivada de vec(Ed-1) por vecE'
                    de_deh = Dup(n)
                    debig_de = kron(TT(n), I(2*n*2*n)) * \
                               kron(I(n), kron(vec(I(2*n)), I(n)))
                    
                    dedinv_devd = diagtom(vec(-inv_ed*inv_ed)) #TODO: conferir
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
                
                omega_evar = 2*Dup(n).I*kron(e_var, e_var)*Dup(n).I.T
            
                varalpha = G1*omega*G1.T #TODO: conferir
                varevar = dpdc_dev*omega_evar*dpdc_dev.T
                varass1[i, j, ff] = (varalpha + varevar)/nd
                
                
                ic1[i, j, ff] = pdc - sqrt(varass1[i, j, ff])*st.norm.ppf(1-alpha/2.0)
                ic2[i, j, ff] = pdc + sqrt(varass1[i, j, ff])*st.norm.ppf(1-alpha/2.0)
                
                G2a = 2*Iije/den
                G2 = Ca.T*G2a*Ca #TODO: conferir!
                
                d = fEig(omega, G2)
                patdf = sum(d)**2/sum(d**2)
                patden = sum(d)/sum(d**2)
                th[i, j, ff] = st.chi2.ppf(1-alpha, patdf)/(patden*2*nd)
                varass2[i, j, ff] = 2*patdf/(patden*2*nd)**2
                
    return th, ic1, ic2, varass1, varass2


def assym_dtf_one(x, Af, e_var, p, alpha = 0.05):
    
    x = mat(x)
    e_var = mat(e_var)
    
    n, nd = x.shape
    nf = Af.shape[0]
    
    th = empty([n, n, nf])
    ic1 = empty([n, n, nf])
    ic2 = empty([n, n, nf])
    varass1 = empty([n, n, nf])
    varass2 = empty([n, n, nf])
    
    gammai = inv(bigautocorr(x, p))
    omega = kron(gammai, e_var)
    
    for ff in range(nf):
        f = ff/(2.0*nf)
        
        Ca = fCa(f, p, n)
        
        Hf = mat(Af[ff, :, :]).I
        h = vec(Hf)
        h = cat(h.real, h.imag, 0)
        #a = vec(cat(I(n), O(n), 1)) - dot(Ca, al)
        
        dhda = fdh_da(mat(Af[ff, :, :]), n)
       
        for i in range(n):
            for j in range(n):
                
                Iij = fIij(i, j, n)
                
                Ii = fIi(i, n)
                
                num = h.T*Iij*h
                den = h.T*Ii*h
                dtf = num/den
    
                G1a = 2*h.T*Iij/den - 2*num*h.T*Ii/(den**2)
                G1 = -G1a*dhda*Ca
                
                varass1[i, j, ff] = G1*omega*G1.T/nd
                ic1[i, j, ff] = dtf - sqrt(varass1[i, j, ff])*st.norm.ppf(1-alpha/2)
                ic2[i, j, ff] = dtf + sqrt(varass1[i, j, ff])*st.norm.ppf(1-alpha/2)
                
                G2a = 2*Iij/den
                G2 = Ca.T*dhda.T*G2a*dhda*Ca #TODO: conferir o .T
                
                d = fEig(omega, G2)
                patdf = sum(d)**2/sum(d**2)
                patden = sum(d)/sum(d**2)
                
                th[i, j, ff] = st.chi2.ppf(1-alpha, patdf)/(patden*2*nd)
                varass2[i, j, ff] = 2*patdf/(patden*2*nd)**2
                
    return th, ic1, ic2, varass1, varass2

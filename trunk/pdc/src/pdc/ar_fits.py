# -*- coding:utf-8 -*-

'''Feito pelo Gilson, precisa revisar. parece bater c matlab.'''

__all__ = ['nstrand', 'yule_walker', 'ar_estim', 'adapt_ar']

from numpy import *
from numpy.linalg import pinv, inv, eig
from scipy.linalg.basic import det
from scipy.linalg.decomp import schur

from pdc.params import *
from matplotlib import pyplot

eps = finfo(float).eps.item()

def lyap(A, B, C=[]):
    '''
    % LYAP Lyapunov equation.
    %  X = LYAP(A,C) solves the special form of the Lyapunov matrix equation:
    %     A*X + X*A' = -C
    %
    %  X = LYAP(A,B,C) solves the general form of the Lyapunov matrix equation:
    %     A*X + X*B = -C
    '''
    if C == []:
        C = B
        B = A.transpose()
    [ma,na] = A.shape
    [mb,nb] = B.shape
    [mc,nc] = C.shape
    if not(ma == na) or not(mb == nb) or not(mc == ma) or not(nc == mb): 
        print 'Dimensions do not agree.'
        return
    #check if problem has any complex inputs
    real_flg = 1
    if A.imag.any() or B.imag.any() or C.imag.any():
        real_flg = 0
    # perform schur decomposition on A and B (note: complex schur form forced by
    # adding small complex part so ua and ub are complex upper triangular)
    #i = 0+1j
    #ta,ua = schur(A+eps*eps*i*A)
    #tb,ub = schur(B+eps*eps*i*B)
    
    ta,ua = schur(A, 'complex')
    tb,ub = schur(B, 'complex')
    
    
    # check all combinations of ua(i,i)+ub(j,j) for zero
    p1 = ta.diagonal().transpose()
    p1.resize(1,p1.size)
    p1 = dot(ones((mb,1)),p1)
    p2 = tb.diagonal()
    p2.resize(p2.size,1)
    p2 = dot(p2,ones((1,ma)))
    summ = abs((p1 + p2)/(abs(p1) + abs(p2)))
    #summ = abs(p1) + abs(p2)
    if (summ < 1000*eps).any() or (isnan(summ)).any():
        print 'Solution is not unique.'
        return
    # transform C
    ucu = dot(dot(-ua.transpose().conj(),C),ub)
    # solve for first column of transformed solution
    y = zeros((ma,mb), complex)
    ema = eye(ma)
    y[:,0] = dot(inv(ta+dot(ema,tb[0,0])), ucu[:,0])
    
    # solve for remaining columns of transformed solution
    for k in range(1,mb):
        km1 = range(0, k)
        y[:,k] = dot(inv(ta+dot(ema,tb[k,k])), (ucu[:,k]-dot(y[:,km1],tb[km1,k])))
    
    # find untransformed solution 
    X = dot(dot(ua,y),ub.transpose().conj())
    # ignore complex part if real inputs (better be small)
    
    if real_flg:
        X = X.real
    return X

def calc_ef(data, A):
    
    n, nd = data.shape
    
    n,n,p = A.shape
    
    x = zeros(data.shape)
    s = zeros(data.shape)
    for i in arange(p,nd):
        #x = zeros(n)
        for j in arange(p):
            x[:,i] += dot(A[:,:,j], data[:,i-j-1])
        s[:,i] = data[:,i]-x[:,i]
        
    return sum(s, axis = 1)

def nstrand(u, p = None, return_ef = False):
    '''
    %   Calculate the coeficients of multi-channel auto-regressive matrix using
    %   Nuttall-Strand algorithm (a generalization of single channel harmonic
    %                             method)
    %
    %   Input parameters:
    %     IP     - Ordem of autoregressive model (integer)
    %     u      - Complex matrix with NUMCHS channels of sample data
    %
    %   Output parameters:
    %     PF     - Covariance matrix of NUMCHS x NUMCHS of linear forward
    %              prediction error
    %     A      - Complex array of forward linear prediction matrix
    %              coefficients
    %     PB     - Complex backward linear prediction error covariance array
    %     B      - Complex array of backward linear prediction matrix
    %              coefficients
    '''
    [lx,cx]=u.shape
    if lx > cx:
        print ('Input matrix is probably transposed.')
        return
    
    if p is None:
        p = pr_.maxp
    
    NUMCHS=lx      #% Number of channels.
    MAXORDER=200   #% Maximum order of AR model allowed for calculation.
    N=max(u.shape) #% N - Number of samples per channel.
    IP = p

    #    Initialization
    ISTAT=0
    if (IP > MAXORDER):
        ISTAT=3
        print('IP > 200')
        return
    
    ef=u.copy()                    #% Eq. (15.91)
    eb=u.copy()                    #% Eq. (15.91)
    pf=dot(u, u.transpose())       #% Eq. (15.90)
    pb=array(pf)                   #% Eq. (15.90)
    M=0
    #    Main Loop
    while 1:
        #%  Update estimated covariance errors  Eq. (15.89)
        pfhat=dot(ef[:,M+1:N],ef[:,M+1:N].transpose())
        pbhat=dot(eb[:,M:N-1],eb[:,M:N-1].transpose())
        pfbhat=dot(ef[:,M+1:N],eb[:,M:N-1].transpose())
        M=M+1
        #%  Calculate estimated partial correlation matrix - Eq. (15.98)
        #%             (Nuttall-Strand algorithm only)
        RHO=lyap(dot(pfhat,inv(pf)),dot(inv(pb),pbhat),-2*pfbhat)
        #%  Update forward and backward reflection coeficients
        #%  Eqs. (15.73),(15.74),(15.78) (algoritmo de Nuttall-Strand)
        AM=dot(-RHO,inv(pb))
        BM=dot(-RHO.transpose(),inv(pf))
        dimA=AM.shape[0]
        dimB=BM.shape[0]
        if M == 1:
            A=zeros((1,dimA,dimA),float)
            B=zeros((1,dimB,dimB),float)
        else:
            #print 'M', M
            #print 'A', A
            #print 'AM', AM
            #print 'AA', A.resize((M,dimA,dimA))
            #print 'dimA', dimA
            A=resize(A,(M,dimA,dimA))
            B=resize(B,(M,dimB,dimB))
        A[M-1,:,:]=AM
        B[M-1,:,:]=BM
        #%  Update forward and backward covariance error  - Eqs. (15.75),(15.76)
        pf=pf-dot(dot(AM,BM),pf)
        pb=pb-dot(dot(BM,AM),pb)
        
        #%  Update forward and backward predictor coeficients - Eqs.(15.84),(15.85)
        if not (M == 1):
            for K in range(1,M):
                temp1=A[K-1,:,:].copy()
                A[K-1,:,:]=A[K-1,:,:]+dot(AM,B[M-K-1,:,:])
                B[M-K-1,:,:]=B[M-K-1,:,:]+dot(BM,temp1)
        #%  Update residues
        #%  existe erro no calculo dos residuos
        Tef=array(ef)
        ef[0:NUMCHS,range(N-1,M-1,-1)]=ef[:,range(N-1,M-1,-1)]+dot(AM,eb[:,range(N-2,M-2,-1)])
        eb[0:NUMCHS,range(N-1,M-1,-1)]=eb[:,range(N-2,M-2,-1)]+dot(BM,Tef[:,range(N-1,M-1,-1)])
        
        #%  Verify if model order is adequate
        if M == IP:
            A=-A
            B=-B
            break
    
    #print 'pf:', pf
    #pf = abs(pf) #TODO: conferir o pf, pq as vezes da matriz toda negativa?
    
    if (return_ef):
        return A.transpose(1,2,0), ef #TODO: porque o abs?
    else:
        #return pf,A,pb,B,ef,eb,ISTAT 
        return A.transpose(1,2,0), pf/N

#def xlag(x, lags):
#    
#    if(lag == 0):
#        return x.copy();
#    xl = matrix(zeros(x.shape))
#    xl[:, lag:] = 
    
def xlags(x, lags):
    
    xl = zeros([lags] + list(x.shape))
    
    xl[0] = x
    for i in arange(1,lags):
        xl[i, :, i:] = x[:, :-i]
        
    return xl

def yule_walker(x, p = None, return_ef = False):
    
    if p is None:
        p = pr_.maxp
    
    xlag = xlags(x, p+1)
    
    n, nd = x.shape
    gamma = zeros([n*(p+1), n*(p+1)])
    for i in arange(p+1):
        for j in arange(p+1):
            gamma[i*n:i*n+n, j*n:j*n+n] = dot(xlag[i], xlag[j].T)/nd
            
    yz = gamma[:n*1,n*1:n*(p+1)]
    igamma = inv(gamma[:n*p,:n*p])
    A = dot(yz, igamma)
    er = gamma[:n*1,:n*1] - dot(dot(yz, igamma), yz.T)

    if return_ef:
        print '\nreturn_ef for YW not implemented yet.'

    return A.reshape(n,p,n).transpose(0,2,1), er
    

#def ar_fit_old(u, MaxIP = 0, alg='ns', criterion=0, return_ef = False):
#    '''
#    %
#    %[IP,pf,A,pb,B,ef,eb,vaic,Vaicv] = mvar(u,MaxIP,alg,criterion)
#    %
#    % input: u     - data rows
#    %        MaxIP - externaly defined maximum IP (default = 30)
#    %        alg   - for algorithm (0: Nutall-Strand)
#    %        criterion for order choice - 0: AIC; 1: fixed order in MaxIP
#    %                                     2(not yet): estimate up to max order
#    %                                     Negative(not yet) - keep criterion changes
#    %
#    '''
#    StopFlag=0
#    [nSegLength,nChannels] = u.transpose().shape
#
#    if criterion<0:
#        stopFlag=1
#        criterion=abs(criterion)
#    
#    if criterion==1:
#        [npf, na, npb, nb, nef, neb, ISTAT]=nstrand(u,MaxIP,False)
#        if (not return_ef):
#            return na.transpose(1,2,0), npf/nSegLength
#        else:
#            return na.transpose(1,2,0), nef
#        
#    if alg == 'yw':
#        fitalg = yule_walker
#    else:
#        fitalg = nstrand
#    
#    vaicv=0
#    if MaxIP == 0:
#        MaxOrder = 30;
#        UpperboundOrder = round(3*sqrt(nSegLength)/nChannels)
#        #% Marple Jr. page 409
#        #% Suggested by Nuttall, 1976.
#        UpperboundOrder = min([MaxOrder, UpperboundOrder])
#    else:
#        MaxOrder=MaxIP
#        UpperboundOrder=MaxIP
#    
#    #print 'MaxOrder limited to ', MaxOrder
#       
#    IP=1
#    Vaicv=zeros((MaxOrder+1,1), float)
#    while IP <= UpperboundOrder:
#        if (IP > 10):
#            print 'Testando ar_estim com P =', IP
#            
#        [npf, na, npb, nb, nef, neb, ISTAT]=fitalg(u,IP,False)
#        
#        vaic=max(u.shape)*log(det(npf))+2*nChannels*nChannels*IP;
#        
#        Vaicv[IP,0]=vaic
#
#        if (vaic>vaicv) and not (IP==1):
#            vaic=vaicv
#            if not StopFlag:
#                break
#        #% Akaike condition
#        vaicv=vaic
#        pf = array(npf)
#        A  = array(na)
#        ef = array(nef)
##        if alg==0:
##            B  = array(nb)
##            eb = array(neb)
##            pb = array(npb)
#        IP=IP+1
#
#    IP=IP-1
#    vaic=vaicv
#    Vaicv=Vaicv[range(1,IP+1),0]
#    Vaicv.shape = (Vaicv.size,1)
#
#    #return IP,pf,A,pb,B,ef,eb,vaic,Vaicv
#    if (not return_ef):
#        return A.transpose(1,2,0), pf/nSegLength
#    else:
#        return A.transpose(1,2,0), ef
    
def ar_estim(data, return_ef = False, **args):
    '''
    %
    %[IP,pf,A,pb,B,ef,eb,vaic,Vaicv] = mvar(u,MaxIP,alg,criterion)
    %
    '''
    
    read_args(args)
    
    [n,nd] = data.shape

    if pr_.ar_estim == 'yw':
        if pr_.v:
            print '\nUsing Yule-Walker VAR estimator.'
        fitalg = yule_walker
    elif pr_.ar_estim == 'ns':
        if pr_.v:
            print '\nUsing Nutall-Strand VAR estimator.'
        fitalg = nstrand
    else:
        if pr_.v:
            print '\nVar estimator invalid, using Nutall-Strand.'
        fitalg = nstrand
    
    if pr_.fixp:
        return fitalg(data,pr_.maxp,return_ef)
    
    #if round(3*sqrt(nd)/n) < pr_.maxp:
    #    print 'Using lower maxp, seems high for this data'
    #    maxp = round(3*sqrt(nd)/n)
       
    aic = zeros(pr_.maxp+1, dtype=float)
    for ip in arange(1,pr_.maxp+1):
        if (ip > 10):
            if pr_.v:
                print 'Testando ar_estim com P =', ip
            
        #[na, npf, npb, nb, nef, neb, ISTAT]=fitalg(data,ip,return_ef)
        na, npf = fitalg(data,ip,return_ef)
        
        aic[ip] = nd*log(det(npf))+2*n*n*ip;
        
        if aic[ip] > aic[ip-1] and ip > 1:
            if not pr_.test_allp:
                break
            
        pf = npf
        A  = na
#       ef = nef
#       B  = array(nb)
#       eb = array(neb)
#       pb = array(npb)
    
    if pr_.test_allp:
        if pr_.v:
            print '\nUsing global AIC minimum:', argmin(aic[1:])+1
        return fitalg(data,argmin(aic[1:])+1)
    
    #return ip-1,pf,A,pb,B,ef,eb,vaic,Vaicv
    return A, pf

#def R_YW(data, maxp = 30):
#    '''Estimates multivariate AR fit for data, using Yule-walker of R package.
#    
#      Input: 
#        data(n, m) - data matrix (n - number of signals, m - data length)
#        maxp - maximum order for estimated AR model
#    
#      Output:
#        A(n, n, p) - estimated AR model
#        er(n, n) - covariance of residuals
#    '''
#    
#    from rpy2.robjects import r as r_
#    import rpy2.rinterface as ri_
#    import rpy2.robjects as ro_
#
#    if (data.ndim == 1):
#        data.resize(1, data.shape[0])
#    
#    ri_.initr()
#    
#    ri_.globalEnv["data"] = ri_.FloatSexpVector(data.ravel())
#    ri_.globalEnv["dim"] = ri_.IntSexpVector(data.shape[::-1])
#    ro_.globalEnv["maxp"] = maxp
#    r_('data <- ar(array(data, dim), order.max = maxp, aic=FALSE)')
#    
#    A = array(r_('data$ar')).transpose(1,2,0) #TODO: conferir A e A.T em todos.
#    er = array(r_('cov(data$resid[-seq(data$order),])'))
#    #print 'Model order: ', array(r_('data$order'))[0]
#    
#    print A
#    
#    return A, er


def adapt_ar (data, p, se = 100):
    '''data(m,n,nd) -> data, m = #trials, n = #channels, nd = #samples'''
    
    m,n,nd = data.shape
    
    A = zeros([nd,n,n*p])
    er = zeros([nd,n,n])
    era = zeros([nd,n,n])
    C = mat(identity(n*p))
        
    #Stein modification: (usar cf em torno de 0.02 para um AR(2,2,2))
    cf = 2.0*m/float64(se+m)
    #print cf
    
    for i in arange(p,nd):
        C = (1.0/(1.0-cf))*C
        if i == p:
            Wt = mat(data[:,:,i-1::-1].transpose(0,2,1).reshape(m,-1))
        else:
            Wt = mat(data[:,:,i-1:i-p-1:-1].transpose(0,2,1).reshape(m,-1))
        
        for j in arange(m):
            C = C*(identity(n*p) - Wt[j].T*Wt[j]*C/(Wt[j]*C*Wt[j].T + 1.0))
        
        K = Wt*C
        #Yn = data(:,:,i)
        Z = data[:,:,i] - Wt*A[i-1].T
        A[i] = A[i-1] + Z.T*K
        
        if cf > 0:
            er[i] = (1-cf)*er[i-1] + (cf/float64(m))*Z.T*Z
            era[i] = er[i]/(1.0-(1.0-cf)**i)
        else:
            er[i] = (n/(n+1.0))*er[i-1] + (1.0/(m*(n+1.0)))*Z.T*Z
            era[i] = er[i]
        
        #if (i%(nd/10) == 0 and i*m > 200):
        #    print 'nd', i
        #    print 'C', C
        #    #print 'Z', Z
        #    print 'dA', Z.T*K 
        

    er = era
    A = A.reshape(nd,n,p,n).transpose(0,1,3,2)

    return A, er


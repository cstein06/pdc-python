'''Feito pelo Gilson, precisa revisar. parece bater c matlab.'''

from numpy import array,zeros,eye,kron,dot,exp,pi,cos,sin,diag,ones,tril,resize,finfo
from numpy.core.fromnumeric import reshape
from numpy.core.umath import isnan
from numpy.lib.scimath import log, sqrt
from numpy.linalg import pinv, inv, eig
from scipy.linalg.basic import det
from scipy.linalg.decomp import schur

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
    i = 0+1j
    ta,ua = schur(A+eps*eps*i*A)
    tb,ub = schur(B+eps*eps*i*B)
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

def nstrand(u, maxp = 30):
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
    NUMCHS=lx      #% Number of channels.
    MAXORDER=200   #% Maximum order of AR model allowed for calculation.
    N=max(u.shape) #% N - Number of samples per channel.
    IP = maxp

    #    Initialization
    ISTAT=0
    if (IP > MAXORDER):
       ISTAT=3
       print('IP > 200')
       return
    
    ef=array(u)                    #% Eq. (15.91)
    eb=array(u)                    #% Eq. (15.91)
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
             temp1=A[K-1,:,:]
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
    pf = abs(pf) #TODO: conferir o pf
    return A.transpose(1,2,0), abs(pf)/N #pf,A,pb,B,ef,eb,ISTAT 

def ar_fit(u,MaxIP=0,alg=1,criterion=1):
    '''
    %
    %[IP,pf,A,pb,B,ef,eb,vaic,Vaicv] = mvar(u,MaxIP,alg,criterion)
    %
    % input: u     - data rows
    %        MaxIP - externaly defined maximum IP
    %        alg   - for algorithm (=1 Nutall-Strand), (=2 - mlsm) ,
    %                              (=3 - Vieira Morf), (=4 - QR artfit)
    %        criterion for order choice - 1: AIC; 2: Hanna-Quinn, 3 Schwartz,
    %                                     4 FPE, 5 - fixed order in MaxIP
    %                                     6 - xdelay,  7,8 - Broersen
    %                                     9 - derivative (future) (envelope)
    %                                     10 - estimate up to max order
    %                                     Negative - keep criterion changes
    % output: See Data Structure Data base Base
    %
    '''
    StopFlag=0
    [nSegLength,nChannels] = u.transpose().shape

    if criterion<0:
        stopFlag=1
        criterion=abs(criterion)
    
    vaicv=0
    if MaxIP == 0:
       MaxOrder=30;
       print 'MaxOrder limited to ', MaxOrder
       UpperboundOrder = round(3*sqrt(nSegLength)/nChannels)
       #% Marple Jr. page 409
       #% Suggested by Nuttall, 1976.
       UpperboundOrder = min([MaxOrder, UpperboundOrder])
    else:
       MaxOrder=MaxIP
       UpperboundOrder=MaxIP
       print 'MaxOrder limited to ', MaxOrder
       
    IP=1
    Vaicv=zeros((MaxOrder+1,1), float)
    while IP <= UpperboundOrder:
        [npf, na, npb, nb, nef, neb, ISTAT]=nstrand(u,IP)
        
        vaic=max(u.shape)*log(det(npf))+2*nChannels*nChannels*IP;
        
        Vaicv[IP,0]=vaic

        if (vaic>vaicv) and not (IP==1):
            vaic=vaicv
            if not StopFlag:
                break
        #% Akaike condition
        vaicv=vaic
        pf = array(npf)
        A  = array(na)
        ef = array(nef)
        if alg==1:
           B  = array(nb)
           eb = array(neb)
           pb = array(npb)
        else:
           #% review status for backward prediction in clmsm
           B  = []
           eb = []
           pb = []
        IP=IP+1

    IP=IP-1
    vaic=vaicv
    Vaicv=Vaicv[range(1,IP+1),0]
    Vaicv.shape = (Vaicv.size,1)

    return IP,pf,A,pb,B,ef,eb,vaic,Vaicv

def R_YW(data, maxp = 30):
    '''Estimates multivariate AR fit for data, using Yule-walker of R package.
    
      Input: 
        data(n, m) - data matrix (n - number of signals, m - data length)
        maxp - maximum order for estimated AR model
    
      Output:
        A(n, n, p) - estimated AR model
        er(n, n) - covariance of residuals
    '''
    
    from rpy2.robjects import r as r_
    import rpy2.rinterface as ri_
    import rpy2.robjects as ro_

    if (data.ndim == 1):
        data.resize(1, data.shape[0])
    
    ri_.initr()
    
    ri_.globalEnv["data"] = ri_.FloatSexpVector(data.ravel())
    ri_.globalEnv["dim"] = ri_.IntSexpVector(data.shape[::-1])
    ro_.globalEnv["maxp"] = maxp
    r_('data <- ar(array(data, dim), order.max = maxp)')
    
    A = array(r_('data$ar')).transpose(1,2,0) #TODO: conferir A e A.T em todos.
    er = array(r_('cov(data$resid[-seq(data$order),])'))
    #print 'Model order: ', array(r_('data$order'))[0]
    return A, er
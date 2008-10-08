from numpy import array,zeros,eye,kron,dot,exp,pi,cos,sin,diag,ones,tril,fliplr
from numpy.linalg import pinv, inv, eig
from numpy.lib.scimath import sqrt
from scipy.stats.distributions import gamma
from matplotlib.mlab import find
import numarray.ieeespecial as ieee

def pdcn(A,pf,nFreqs,nSegLength,gamma,option,metric,aValue):
    '''
    [SS,Lpdc,LTra,Lpdcvinf,Lpdcvsup,Lpatnaik]
    %Compute Partial^Directed^Coherence from MAR estimate
    %
    %[SS Lpdc LTra Lpdcvinf Lpdcvsup Lpatnaik] = ...
    %                        pdcn(A,pf,nFreqs,nSegLength,gamma,option,metric,aValue)
    % input: A          - MAR model matrix
    %        pf         - covariance matrix
    %        nSegLength - signal sample points
    %        nFreqs     - 
    %        gamma      -
    %        option     -
    %        metric     -
    %        aValue     -
    %
    % output: SS       - Power spectrum
    %         Lpdc     - normalized PDC
    %         LTra     - PDC above threshold value
    %         Lpdcvinf - pdc confidence interval
    %         Lpdcvsup - superior and inferior
    %         Lpatnaik - threshold calculated by Patnaik approximation
    '''
    [SS,Lpdc,Lpatnaik,fvar] = prpardd(A,pf,nFreqs,nSegLength,gamma,option,metric,aValue)

    #% Significant PDC values on frequency scale
    Ltemp = (abs(Lpdc)**2-Lpatnaik**2>0)*Lpdc+(abs(Lpdc)**2-Lpatnaik**2<=0)*(-1)
    #Ltemp(ind2sub(size(Ltemp),find(Ltemp==-1))) = NaN;
    for i in range(0,Ltemp.shape[0]):
	for j in range(0,Ltemp.shape[1]):
		for k in range(0,Ltemp.shape[2]):
			if Ltemp[i,j,k] == -1.0:
				Ltemp[i,j,k] = ieee.nan
    LTra = array(Ltemp)
    #Ltemp(ind2sub(size(Ltemp),find(Ltemp~=-1))) = 1;
    for i in range(0,Ltemp.shape[0]):
	for j in range(0,Ltemp.shape[1]):
		for k in range(0,Ltemp.shape[2]):
			if not (Ltemp[i,j,k] == 1.0):
				Ltemp[i,j,k] = 1
    #% Upper and lower variances are only defined for significant PDC range.
    #% Variance upper
    Lpdcvsup = sqrt(abs(Lpdc)**2 + fvar)*Ltemp
    #% Variance lower
    Lpdcvinf = sqrt(abs(Lpdc)**2 - fvar)*Ltemp
    return [SS,Lpdc,LTra,Lpdcvinf,Lpdcvsup,Lpatnaik]
    
def prpardd(A,pf,nFreqs,np,G,option=2,metric=1,avalue=0.95):
    '''
    [SS,L,th,fvar]
    %Compute Partial^Directed^Coherence from series j to i.
    %
    %function [SS,L]=prpardd(A,pf,nFreqs,option,metric)
    %
    % input: A - power spectral density estimate of the i-th time-series
    %        pf - covariance matrix
    %        nFreqs - number of point in [0,fs/2]
    %        np - number of sample points
    %        G
    %        option - 0 - L=[classic coherence], 1 (L=dircor), 2  L=pardir,
    %                 3   L=transfer.   4- partial coherence
    %        metric   1 - Euclidean,
    %                 2 - diagonal,
    %                 3 - statistical,
    %                 4 - unnormalized,
    %                 5 - weighted by Si(f) receiver
    %                 6 - weighted by Sj(f) sender
    %                 7 - joint weighted
    %        avalue = .95 default for distribution
    %        Note: Pdcn - metric=2 - option=2
    %
    % output: SS
    %         SS - power spectrum for each sequence.
    %         According to option
    %         L = G -(dircor) or
    %         L= PG (pardir
    %         L= H - transfer function gain - modified 11:55 AM 3/22/2000 (LAB)
    %         th - Threshold value with (1-avalue) significance level.
    %         fvar - confidence interval variance for significant PDC|DTF|DC|Coh
    %                
    '''
    mflag=0
    [nIP,nChannels,n0]=A.shape
    SS=zeros((nFreqs,nChannels,nChannels),complex)
    L=zeros((nFreqs,nChannels,nChannels),complex)
    nL=zeros((nFreqs,nChannels,nChannels),complex)
    dL=zeros((nFreqs,nChannels,nChannels),complex)
    th=zeros((nFreqs,nChannels,nChannels),float)
    fvar=zeros((nFreqs,nChannels,nChannels),float)

    optionr=0
    if (metric==5)or(metric==6)or(metric==7):
       mflag=metric
       metric=1

    if option==4: #%Partial coherence
       metric=3
       optionr=1
       option=2

    if metric == 1:
        #% Euclidean
        pu=eye(nChannels)
    elif metric == 2:
        #% Diagonal
        pu=diag(diag(pf))
        #print '\npu', pu
    elif metric == 3:
        #% Statistical
        pu=pf #% needs revision
    elif metric == 4:
        #% Unnormalized
        pu=eye(nChannels)

    if option==2:
        pu=pinv(pu) #% For PDC
        #print '\npu', pu

    Omega=kron(inv(G),dot(dot(sqrt(pu),pf),dot(sqrt(pu),np)))
    #print '\nOmega', Omega
    Omega_a=kron(array([[1, 1],[1, 1]],float),Omega)
    #print '\nOmega_a', Omega_a

    for iFreqs in range(0,nFreqs):
        f=iFreqs/(2.0*nFreqs)
        AL=eye(nChannels,nChannels)
        for j in range(1,nIP+1):
          AL=AL-dot(A[j-1,:,:],exp(-sqrt(-1.0)*2.0*j*pi*f))
        AT=pinv(AL)
        #print '\nAT', AT

        #% Matrix C
        Cmat=zeros((0,0),float)
        Smat=zeros((0,0),float)

        for r in range(1,nIP+1):
            divect=2*pi*f*r
            cvector=cos(divect)
            svector=-sin(divect)
            C1 = eye(nChannels**2)*cvector
            S1 = eye(nChannels**2)*svector
            C2 = array(Cmat)
            S2 = array(Smat)
            Cmat=zeros((nChannels**2,C2.shape[1]+C1.shape[1]),float)
            Smat=zeros((nChannels**2,S2.shape[1]+S1.shape[1]),float)
            if C2.shape[1] > 0:
                Cmat[:,0:C2.shape[1]]=C2
            if S2.shape[1] > 0:
                Smat[:,0:S2.shape[1]]=S2
            Cmat[:,C2.shape[1]:C2.shape[1]+C1.shape[1]]=C1
            Smat[:,S2.shape[1]:S2.shape[1]+S1.shape[1]]=S1
        Zmat=zeros(Cmat.shape, float)
        Cjoint=zeros((2*Cmat.shape[0],2*Cmat.shape[1]), float)
        Cjoint[0:Cmat.shape[0],0:Cmat.shape[1]]=Cmat
        Cjoint[0:Cmat.shape[0],Cmat.shape[1]:2*Cmat.shape[1]]=Zmat
        Cjoint[Cmat.shape[0]:2*Cmat.shape[0],0:Cmat.shape[1]]=Zmat
        Cjoint[Cmat.shape[0]:2*Cmat.shape[0],Cmat.shape[1]:2*Cmat.shape[1]]=-Smat
        Ct=dot(dot(Cjoint,Omega_a),Cjoint.transpose())
        if option==1:
            TH=-kron(AT.transpose(),AT)
            THr=TH.real
            THi=TH.imag
            TH =zeros((THr.shape[0]+THi.shape[0],THr.shape[1]+THi.shape[1]),complex)
            TH[0:THr.shape[0],0:THr.shape[1]]=THr
            TH[0:THi.shape[0],THr.shape[1]:THr.shape[1]+THi.shape[1]]=-THi
            TH[THr.shape[0]:THr.shape[0]+THi.shape[0],0:THi.shape[1]]=THi
            TH[THi.shape[0]:THr.shape[0]+THi.shape[0],THi.shape[1]:THr.shape[1]+THi.shape[1]]=THr
            Ct=dot(dot(TH.transpose(),Ct),TH)

        #print '\nCt', Ct
        SU=dot(dot(AT,pf),AT.transpose())
        #print '\nSU', SU
        SS[iFreqs,:,:]=SU
        
        if option==0:  #% Classical Coherence computation
            #%
            #%--vector version--
            #%
            #% SV=inv(sqrt(diag(diag(SU))));
            #% L(:,:,iFreqs+1)=SV*SU*SV;
            #%
            for iu in range(1,nChannels+1):
                for ju in range(1,nChannels+1):
                    L[iFreqs,iu-1,ju-1]=SU[iu-1,ju-1]/sqrt(SU[iu-1,iu-1]*SU[ju-1,ju-1])
        
        else: #% option~=0
            for iu in range(1,nChannels+1):
                for ju in range(1,nChannels+1):
                    #% eigenvalue computation
                    Co1 = Ct[(ju-1)*nChannels+iu-1,(ju-1)*nChannels+iu-1]
                    Co2 = Ct[(ju-1)*nChannels+iu-1,(nChannels+ju-1)*nChannels+iu-1]
                    Co3 = Ct[(nChannels+ju-1)*nChannels+iu-1,(ju-1)*nChannels+iu-1]
                    Co4 = Ct[(nChannels+ju-1)*nChannels+iu-1,(nChannels+ju-1)*nChannels+iu-1]
                    Co=array([[Co1, Co2],[Co3, Co4]])
                    v=eig(Co)[0]

                    Pat=gamma.ppf(avalue,.5*sum(v)**2/(2*dot(v.transpose(),v)),scale=2)
                    #print '\nPat', Pat

                    if option==1: #% DTF/DC
                        
                        if metric==1:
                            SP=dot(AT,AT.transpose())
                            L[iFreqs,iu-1,ju-1]=AT[iu-1,ju-1]*sqrt(pu[ju-1,ju-1])/sqrt(SP[iu-1,iu-1]) #% Coherences
                            nL[iFreqs,iu-1,ju-1]=AT(iu-1,ju-1)*sqrt(pu[ju-1,ju-1])
                            dL[iFreqs,iu-1,ju-1]=SP[iu-1,iu-1]
                            th[iFreqs,iu-1,ju-1]=Pat/np*abs(dL[iFreqs,iu-1,ju-1])**2

                            dLu=abs(dL[iFreqs,iu-1,ju-1])
                            nLu=abs(nL[iFreqs,iu-1,ju-1])**2
                            th[iFreqs,iu-1,ju-1]=sqrt(Pat/((sum(v)/(2*dot(v.transpose(),v)))*(np*abs(dL[iFreqs,iu-1,ju-1]))))
                            #%
                            I_ij=zeros((nChannels**2,nChannels**2), float)
                            I_j=zeros((nChannels**2,nChannels**2), float)
                            I_ij[(ju-1)*nChannels+iu-1,(ju-1)*nChannels+iu-1]=1

                            I_t=find(range(1,nChannels**2+1)%nChannels+iu-1==1)

                            for i in I_t:
                              for j in I_t:
                                  I_j[i,j]=1

                            I_j=diag(diag(I_j))
                            Icj=zeros(array(I_j.shape)*2,float)
                            Ic=zeros(array(I_ij.shape)*2,float)
                            Icj[0:I_j.shape[0],0:I_j.shape[1]]=I_j
                            Icj[I_j.shape[0]:2*I_j.shape[0],I_j.shape[1]:2*I_j.shape[1]]=I_j
                            Ic[0:I_ij.shape[0],0:I_ij.shape[1]]=I_ij
                            Ic[I_ij.shape[0]:2*I_ij.shape[0],I_ij.shape[1]:2*I_ij.shape[1]]=I_ij

                            #%
                            rr=AT.reshape(nChannels**2,1)
                            a=zeros((rr.shape[0],2*rr.shape[1]),float)
                            a[:,0:rr.shape[1]]=rr.real
                            a[:,rr.shape[1]:2*rr.shape[1]]=rr.imag

                            G=2*(dot(Ic,a)/dLu - (dot(Icj,a)*nLu)/(dLu**2))

                            fvar[iFreqs,iu-1,ju-1]=dot(dot(G.transpose(),Ct),G)

                        elif metric==2:
                            SP=dot(dot(AT,pu),AT.transpose())
                            L[iFreqs,iu-1,ju-1]=AT[iu-1,ju-1]*sqrt(pu[ju-1,ju-1])/sqrt(SP[iu-1,iu-1])
                            nL[iFreqs,iu-1,ju-1]=AT[iu-1,ju-1]*sqrt(pu[ju-1,ju-1])
                            dL[iFreqs,iu-1,ju-1]=SP[iu-1,iu-1]

                            dLu=abs(dL[iFreqs,iu-1,ju-1])
                            nLu=abs(nL[iFreqs,iu-1,ju-1])**2
                            th[iFreqs,iu-1,ju-1]=sqrt(Pat/((sum(v)/(2*dot(v.transpose(),v)))*(np*abs(dL[iFreqs,iu,ju]))))
                            #%
                            I_ij=zeros((nChannels**2,nChannels**2), float)
                            I_j=zeros((nChannels**2,nChannels**2), float)
                            I_ij[(ju-1)*nChannels+iu-1,(ju-1)*nChannels+iu-1]=1
                            #%1
                            I_t=find(range(1,nChannels^2+1)%nChannels+iu-1==1)
                            for i in I_t:
                              for j in I_t:
                                  I_j[i,j]=1
                            #% 1 review
                            I_j=diag(diag(I_j))
                            Icj=zeros(array(I_j.shape)*2,float)
                            Ic=zeros(array(I_ij.shape)*2,float)
                            Icj[0:I_j.shape[0],0:I_j.shape[1]]=I_j
                            Icj[I_j.shape[0]:2*I_j.shape[0],I_j.shape[1]:2*I_j.shape[1]]=I_j
                            Ic[0:I_ij.shape[0],0:I_ij.shape[1]]=I_ij
                            Ic[I_ij.shape[0]:2*I_ij.shape[0],I_ij.shape[1]:2*I_ij.shape[1]]=I_ij
                            #%
                            rr=dot(pu,AT).reshape(nChannels**2,1)
                            a=zeros((rr.shape[0],2*rr.shape[1]),float)
                            a[:,0:rr.shape[1]]=rr.real
                            a[:,rr.shape[1]:2*rr.shape[1]]=rr.imag
                            G=2*(dot(Ic,a)/dLu - (dot(Icj,a)*nLu)/(dLu**2))
                            fvar[iFreqs,iu-1,ju-1]=dot(dot(G.transpose(),Ct),G)
                            #%========================================================
                        elif metric==3:
                            #% review
                            L[iFreqs,iu-1,ju-1]=AT[iu-1,ju-1]*sqrt(pf[ju-1,ju-1])/sqrt(SU[iu-1,iu-1]) 
                            nL[iFreqs,iu-1,ju-1]=AT[iu-1,ju-1]*sqrt(pu[ju-1,ju-1])
                            dL[iFreqs,iu-1,ju-1]=sqrt(SU[iu-1,iu-1])
                            th[iFreqs,iu-1,ju-1]=Pat/((sum(v)**2/(2*dot(v.transpose(),v)))*(np*abs(dL[iFreqs,iu-1,ju-1])**2))
                            #%iu,ju
                            #%pause ==================================================
                        elif metric==4:
                            L[iFreqs,iu-1,ju-1]=AT[iu-1,ju-1]
                    #% ==========PDC/PDCn===========================================
                    if option==2:
                        if metric==4:
                            L[iFreqs,iu-1,ju-1]=AL[iu-1,ju-1]
                        elif metric==2: #% PDCn

                            L[iFreqs,iu-1,ju-1]=sqrt(pu[iu-1,iu-1])*AL[iu-1,ju-1]/sqrt(dot(dot(AL[:,ju-1].transpose(),pu),AL[:,ju-1]))
                            nL[iFreqs,iu-1,ju-1]=sqrt(pu[iu-1,iu-1])*AL[iu-1,ju-1]
                            dL[iFreqs,iu-1,ju-1]=dot(dot(AL[:,ju-1].transpose(),pu),AL[:,ju-1])
                            #%
                            th[iFreqs,iu-1,ju-1]=Pat/pu[iu-1,ju-1]*np*abs(dL[iFreqs,iu-1,ju-1])
                            #%
                            dLu=abs(dL[iFreqs,iu-1,ju-1])
                            nLu=abs(nL[iFreqs,iu-1,ju-1])**2
                            th[iFreqs,iu-1,ju-1]=sqrt(Pat/((sum(v)/(2*dot(v.transpose(),v)))*(np*abs(dL[iFreqs,iu-1,ju-1]))))
                            #%
                            I_ij=zeros((nChannels**2,nChannels**2), float)
                            I_j=zeros((nChannels**2,nChannels**2), float)
                            I_ij[(ju-1)*nChannels+iu-1,(ju-1)*nChannels+iu-1]=1
                            I_j[(ju-1)*nChannels:ju*nChannels,(ju-1)*nChannels:ju*nChannels]=1
                            I_j=diag(diag(I_j))
                            Icj=zeros(array(I_j.shape)*2,float)
                            Ic=zeros(array(I_ij.shape)*2,float)
                            Icj[0:I_j.shape[0],0:I_j.shape[1]]=I_j
                            Icj[I_j.shape[0]:2*I_j.shape[0],I_j.shape[1]:2*I_j.shape[1]]=I_j
                            Ic[0:I_ij.shape[0],0:I_ij.shape[1]]=I_ij
                            Ic[I_ij.shape[0]:2*I_ij.shape[0],I_ij.shape[1]:2*I_ij.shape[1]]=I_ij
                            #%
                            rr=dot(sqrt(pu),AL).transpose().reshape(nChannels**2,1)
                            #a=zeros((rr.shape[0],2*rr.shape[1]),float)
                            a=zeros((2*rr.shape[0],1),float)
                            a[0:rr.shape[0],:]=rr.real
                            a[rr.shape[0]:2*rr.shape[0]]=rr.imag
                            #print '\nCoi1', dot(Icj,a)*nLu/(dLu**2)
                            #print '\nCoi2', dot(Ic,a)/dLu
                            #print '\nCoi3', 2*(dot(Ic,a)/dLu - (dot(Icj,a)*nLu)/(dLu**2))
                            G=2*(dot(Ic,a)/dLu - (dot(Icj,a)*nLu)/(dLu**2))

                            #%% Confidence interval computation for generalized PDC
                            #%% by DY Takahashi 21July2007
                            foq = diag(sqrt(pu))
                            foq.shape=(foq.size,1)
                            foq = kron(ones((2*nChannels,1)),foq)
                            weight=diag(foq.reshape(foq.size,))
                            #print '\nG', G
                            #print '\nweight', dot(dot(dot(dot(G.transpose(),weight.transpose()),Ct),weight),G)
                            fvar[iFreqs,iu-1,ju-1]=dot(dot(dot(dot(G.transpose(),weight.transpose()),Ct),weight),G)
                            #%norminv(1-aValue/2,0,1)/sqrt(np);

                        elif metric==1:  #% PDC original
                            L[iFreqs,iu-1,ju-1]=AL[iu-1,ju-1]/sqrt(dot(dot(AL[:,ju-1].transpose(),pu),AL[:,ju-1]))
                            nL[iFreqs,iu-1,ju-1]=AL[iu-1,ju-1]
                            dL[iFreqs,iu-1,ju-1]=dot(dot(AL[:,ju-1].transpose(),pu),AL[:,ju-1])
                            #%
                            dLu=abs(dL[iFreqs,iu-1,ju-1])
                            nLu=abs(nL[iFreqs,iu-1,ju-1])**2
                            th[iFreqs,iu-1,ju-1]=sqrt(Pat/((sum(v)/(2*dot(v.transpose(),v)))*(np*abs(dL[iFreqs,iu-1,ju-1]))))
                            #%
                            I_ij=zeros((nChannels**2,nChannels**2), float)
                            I_j=zeros((nChannels**2,nChannels**2), float)
                            I_ij[(ju-1)*nChannels+iu-1,(ju-1)*nChannels+iu-1]=1
                            I_j[(ju-1)*nChannels:ju*nChannels,(ju-1)*nChannels:ju*nChannels]=1
                            I_j=diag(diag(I_j))
                            Icj=zeros(array(I_j.shape)*2,float)
                            Ic=zeros(array(I_ij.shape)*2,float)
                            Icj[0:I_j.shape[0],0:I_j.shape[1]]=I_j
                            Icj[I_j.shape[0]:2*I_j.shape[0],I_j.shape[1]:2*I_j.shape[1]]=I_j
                            Ic[0:I_ij.shape[0],0:I_ij.shape[1]]=I_ij
                            Ic[I_ij.shape[0]:2*I_ij.shape[0],I_ij.shape[1]:2*I_ij.shape[1]]=I_ij
                            #%
                            rr=AL.reshape(nChannels**2,1)
                            a=zeros((rr.shape[0],2*rr.shape[1]),float)
                            a[:,0:rr.shape[1]]=rr.real
                            a[:,rr.shape[1]:2*rr.shape[1]]=rr.imag
                            G=2*(dot(Ic,a)/dLu - (dot(Icj,a)*nLu)/(dLu**2))
                            fvar[iFreqs,iu-1,ju-1]=dot(dot(G.transpose(),Ct),G)
                            #%*norminv(1-avalue/2,0,1)/sqrt(np)
                            #% Partial directed Coherences
                    #%======= Transfer function =====================================
                    if option==3:
                        L[iFreqs,iu-1,ju-1]=-AL[iu-1,ju-1]/AL[iu-1,iu-1] #% Baccala et 1998
                    #%======= Partial coherence =====================================
                    if optionr==1: #% Partial coherence
                        L[iFreqs,iu-1,ju-1]=dot(dot(AL[:,iu-1].transpose(),pu),AL[:,ju-1])/sqrt(dot((dot(dot(AL[:,iu-1].transpose(),pu,AL[:,iu-1]))),(dot(dot(AL[:,ju-1].transpose(),pu),AL[:,ju-1]))))

    if mflag>0:
        for iu in range(1,nChannels+1):
            for ju in range(1,nChannels+1):
                if mflag==5:
                    L[:,iu-1,ju-1]=L[:,iu-1,ju-1]*SS[:,ju-1,ju-1]/max(SS[:,ju-1,ju-1])
                elif mflag==6:
                    L[:,iu-1,ju-1]=L[:,iu-1,ju-1]*SS[:,iu-1,iu-1]/max(SS[:,iu-1,iu-1])
                elif mflag==7:
                    L[:,iu-1,ju-1]=L[:,iu-1,ju-1]*SS[:,iu-1,iu-1]*SS[:,ju-1,ju-1]/(max(SS[:,ju-1,ju-1])*max(SS[:,iu-1,iu-1]))

    if option == 0:
        print 'option = 0: Classical Coherence'
    elif option == 1:
        print 'option = 1: DTF/DC'
    elif option ==  2:
        if metric == 2:
            print 'option = 2 metric = 2: Generalized Partial Directed Coherence'
        elif metric == 1:
            print 'option = 2 metric = 1: Partial Directed Coherence original'
        else:
            print 'option = 2 metric <> 0 or 1: Unexpected'
    elif option == 3:
        print 'option = 3: Transfer function'
    elif option == 4:
        print 'option = 4: Partial Coherence'
    else:
        print 'option = 0: Partial Coherence'
          
    return [SS,L,th,fvar]

def granmaty(pf,N,significance):
    '''
    % Test Granger causality structure
    %
    %[Tr,Va,v,th,pValue]=granmaty(SU,N,significance);
    % Program to test granger causality structure
    %
    % input: N (number of points)
    %        pf - covariance of modelling errors
    %        significance - test significance level
    %
    % output: Tr -test result matrix (i,j) entry=1 j->i causality cannot
    %             be rejected
    %         Va - test value matrix
    %         v  - degrees of freedom
    %         th - threshold value
    %         pValue - test p-value
    %
    % % 01/30/1998 - L.A.B.
    %
    % disp('Instantaneous Granger causality test: ');
    % significance
    '''
    [n, m]=pf.shape
    Va=zeros((n,m),float)
    Tr=zeros((n,m),float)
    CO=zeros((n,m),float)
    pValue=zeros((n,m),float)
    for i in range(1,n+1):
        for j in range(1,n+1):
            if i>j:
                CO[i-1,j-1]=1
                [Tr[i-1,j-1],Va[i-1,j-1],v,th,pValue[i-1,j-1]]=instata(CO,pf,N,significance)
                Tr[j-1,i-1]=Tr[i-1,j-1]
                Va[j-1,i-1]=Va[i-1,j-1]
                CO[i-1,j-1]=0
                pValue[j-1,i-1]=pValue[i-1,j-1]
    return [Tr,Va,v,th,pValue]

#%==========================================================================
def instata(CO,pf,N,significance):
    '''
    % Test for instataneous causality
    % input: CO - matrix describing the structure for testing - 1 position to test.
    %        pf - residual covariance
    %        N - number of poinst
    %
    % output: y - test result - 0 instantaneous causality rejected - 1 not rejected
    %         value - test value
    %         v - degrees of freedom # constraints.
    %         th -threschold
    '''
    si=vech(pf)
    CO=tril(CO)
    [m, n]=CO.shape
    lb=max(si.shape)
    Ct=vech(CO)
    #print Ct.transpose()
    Ct1=zeros(Ct.transpose().shape,float)
    Ctf=zeros((0,0),float)
    l=sum(Ct.transpose())
    for i in range(1,max(Ct.shape)+1):
        if Ct[i-1]==1:
            Ct1[0,i-1]=1
            C2 = array(Ctf)
            Ctf=zeros((C2.shape[0]+Ct1.shape[0],Ct1.shape[1]),float)
            if C2.shape[0]>0:
                Ctf[0:C2.shape[0],:]=C2
            Ctf[C2.shape[0]:C2.shape[0]+Ct1.shape[0],:]=Ct1
            Ct1=zeros(Ct.transpose().shape,float)
    C=array(Ctf)
    ln=max(pf.shape)
    D=pinv(dmatrix(ln))
    value=N*dot(dot(dot((dot(C,si)).transpose(),inv(2*dot(dot(dot(dot(C,D),kron(pf,pf)),D.transpose()),C.transpose()))),C),si)
    v=2 #% Chi-square distribution degree of freedom.
    #th=chi2inv(significance,v)
    th=gamma.ppf(significance,v/2,scale=2)
    y=value>=th
    #pValue=1-chi2cdf(value,v) #% p-value of instantaneous Granger causality test
    pValue=1-gamma.cdf(value,v/2,scale=2)
    return [y,value,v,th,pValue]

def vech(Y):
    y=zeros((0,0),float)
    [m, n]=Y.shape
    for i in range(1,m+1):
        y1 = array(y)
        y = zeros((y1.shape[0]+Y[i-1:n,i-1].shape[0],1),float)
        if y1.shape[0]>0:
                y[0:y1.shape[0],:]=y1
        y[y1.shape[0]:y1.shape[0]+Y[i-1:n,i-1].shape[0],:]=Y[i-1:n,i-1].reshape(n-i+1,1)
        #y=[y ;Y(i:n,i)]
    return y
#%==========================================================================
#%  vec=D vech
#%
#%  01/30/1998 - L.A.B.
#%
def dmatrix(m):
    '''
    return D
    '''
    D=zeros((m*m,m*(m+1)/2),float)
    u=zeros((0,0),float)
    v=zeros((0,0),float)
    for j in range(1,m+1):
        for i in range(1,m+1):
            u1 = array(u)
            u = zeros((u1.shape[0]+1,2),float)
            if u1.shape[0]>0:
                u[0:u1.shape[0],:]=u1
            u[u1.shape[0]:u1.shape[0]+1,:]=array([i,j],float)
            #u=[u ;[i j]];
            if j<=i:
                v1 = array(v)
                v = zeros((v1.shape[0]+1,2),float)
                if v1.shape[0]>0:
                    v[0:v1.shape[0],:]=v1
                v[v1.shape[0]:v1.shape[0]+1,:]=array([i,j],float)
                #v=[v ;[i j]]
                
             
    w=fliplr(v)
    for i in range(1,m*m+1):
        for j in range(1,m*(m+1)/2+1):
            if sum(u[i-1,:]==v[j-1,:])==2:
                D[i-1,j-1]=1
        for j in range(1,m*(m+1)/2+1):
            if sum(u[i-1,:]==w[j-1,:])==2:
                D[i-1,j-1]=1
    return D
                                            
A = array([[[0.8280,-12.1097],[0.0009, -0.1258]],[[-0.0724,13.4798],[ 0.0036,-0.1391]],[[-0.3561,-36.4805],[0.0013,-0.0735]]], float)
pf = array([[568.0873,-1.5815],[-1.5815,0.0474]],float)
nFreqs = 64
nSegLength = 37
G =array([[9.9104,0.0149,7.2786,-0.0097,2.9629,-0.0243],
[0.0149,0.0004,0.0397,0.0001,0.0449,0.0000],
[7.2786,0.0397,9.9080,0.0149,7.2635,-0.0098],
[-0.0097,0.0001,0.0149,0.0004,0.0393,0.0001],
[2.9629,0.0449,7.2635,0.0393,9.8163,0.0142],
[-0.0243,0.0000,-0.0098,0.0001,0.0142,0.0004]],float)
G=1.0e+004*G
option = 2
metric = 2
aValue = 0.9500
#[SS,Lpdc,LTra,Lpdcvinf,Lpdcvsup,Lpatnaik] = pdcn(A,pf,nFreqs,nSegLength,G,option,metric,aValue)
#[SS,Lpdc,Lpatnaik,fvar] = prpardd(A,pf,nFreqs,nSegLength,G,option,metric,aValue)
pf = array([[0.2413,-0.0556],[-0.0556,0.7496]],float)
nSegLength = 20
igt_signif = 0.9500
[Tr,Va,v,th,pValue]=granmaty(pf,nSegLength,igt_signif)


p<-500   #dimension
n<-300    #sample size
ncycle<-200   #circle times

#True Factor Number
r=5; rmax<-10

## for Calculating ACT
under_m <- function(n, j, z, hatEigValues, p){
  rho = (p-j)/(n-1)
  return( -(1-rho)/z + rho*m(n,j,z, hatEigValues, p) )
}
m <- function(n, j, z, hatEigValues, p){
  return( (p - j)^(-1) * ( sum( (hatEigValues[(j+1):p] - z)^(-1) ) + ( (3*hatEigValues[j] + hatEigValues[j+1])/4 - z )^(-1) ) )
}


##  initinalize
estFN_by_ACT = rep(0,ncycle) #estimated Factor Number by ACT
estFN_by_ACTT = rep(0,ncycle)
estFN_by_all=matrix(0,ncol=13,nrow=ncycle) #estimated Factor Numbers using all methods



#### For other Methods
{
    T = n; N = p
    PC1f=PC2f=PC3f=IC1f=IC2f=IC3f=zf=ON1f=ON2f=ON3f=NONf=ERf=GRf<-rep(0,ncycle)
    v<-rep(0, rmax); bpn<-(p+n)/(p*n); cpn<-min(p,n)
    kfactor<-1:rmax
    
    rho=0; beta0=0; J=0; nc=sqrt((1-rho)^2/(1+2*J*beta0^2)); 
    cm11<-matrix(rep(0,N^2),N,N)
    for(i in 1:N-1){
      for(j in (i+1):min((i+J),N)){
        cm11[i,j]<-beta0
      }
    }
    cm1<-cm11+t(cm11)+diag(rep(1,N))
    cm2<-matrix(rep(0,T^2),T,T)
    for(i in 1:T){
      for(j in i:T){
        cm2[i,j]<-rho^(j-i)
      }
    }
    hm1=hm2=hm3<-matrix(rep(0,N*rmax),N,rmax)
    for(i in 1:rmax){
      for(j in i:N){
        hm1[j,i]<-1
      }
    }		
    for(i in 1:rmax){
      for(j in (i+1):N){
        hm2[j,i]<-1
      }
    }
    for(i in 1:rmax){
      for(j in (i+2):N){
        hm3[j,i]<-1
      }
    }
}


for (KKK in 1:ncycle) {
  
  ###  Define Factor Model
  ## Model 2
  #beta<-matrix(rnorm(p*r,0,1),p,r); c<-runif(p,0,20); diagc<-diag(sqrt(c))*3
  
  ## Model 3
  #beta<-matrix(rnorm(p*r,0,1),p,r); diagc=diag(rep(sqrt(1*r),p))*6
  
  ## Model 4
  #beta<-matrix(rnorm(p*r,0,0.2),p,r); for (j in 1:r){beta[j,j]=1}
  #diagc<-diag(sqrt( runif(p,0,5.5) ))
  
  ## Custom Model1  # When loading matrix is diagonal
  #beta<-diag(1, nrow = p, ncol = r); diagc<-diag(sqrt( runif(p,0,5.5) ))
  
  ## Custom Model2  # no error term
  beta<-matrix(rnorm(p*r,0,1),p,r); diagc<-diag(rep(0, p))
  
  ## Custom Model3.1  # Multicollinearity, one column
  #A=matrix(NA, nrow = p, ncol = r)
  #A[,1:(r-1)] = matrix( rnorm(p*(r-1),0,1), nrow = p, ncol = r-1)
  #A[,r]=runif(r-1) %*% t( A[,1:r-1] )  #make last column be a linear combination of other columns.
  #beta=A
  #diagc<-diag(sqrt( runif(p,0,20) ))*3
  
  ## Custom Model3.2  # Multicollinearity, two columns
  #A=matrix(NA, nrow = p, ncol = r)
  #A[,1:(r-2)] = matrix( rnorm(p*(r-2),0,1), nrow = p, ncol = r-2)
  #A[,(r-1)]=runif(r-2) %*% t( A[,1:(r-2)] )  #make column be a linear combination of other columns.
  #A[,r]=runif(r-2) %*% t( A[,1:(r-2)] )  #make column be a linear combination of other columns.
  #beta=A
  #diagc<-diag(sqrt( runif(p,0,20) ))*3
  
  ## Custom Model4.1  # different scale for last column
  #A=matrix(NA, nrow = p, ncol = r)
  #A[,1:(r-1)] = matrix( rnorm(p*(r-1),0,1), nrow = p, ncol = r-1)
  #A[,r] = runif(p, 0, 0.01)  #different scale
  #beta=A
  #diagc<-diag(sqrt( runif(p,0,20) ))*3
  
  ## Custom Model4.2  # different scale for last two columns
  A=matrix(NA, nrow = p, ncol = r)
  A[,1:(r-2)] = matrix( rnorm(p*(r-2),0,1), nrow = p, ncol = r-2)
  A[,(r-1)] = runif(p, 0, 0.01)  #different scale
  A[,r] = runif(p, 0, 0.01)  #different scale
  beta=A
  diagc<-diag(sqrt( runif(p,0,20) ))*3
  
  
  #Population Correlation Matrix
  RR=cov2cor(beta%*%t(beta)+(diagc)^2) 
  lambdaRR=eigen(RR)$values #eigenvalues of Population Correlation Matrix
  
  
  ###  Data Generating
  #####   Gaussian
  factors<-matrix(rnorm(r*n,0,1),r,n); epsilon<-matrix(rnorm(n*p,0,1),p,n) 
  #####   Uniform 
  #factors<-matrix(runif(r*n,0,2),r,n)*sqrt(3); epsilon<-matrix(runif(n*p,0,2),p,n)*sqrt(3)
  
  Data = beta%*%factors + diagc%*%epsilon
  X = Data - matrix( rep(t(apply(Data,1,mean)),n), ncol=n ) #make Data centered
  
  
  #### Estimate the number of factors
  ## Method 1: the method of zheng: estFN_by_all[,13]
  sampleCovMat = cov(t(X)); hatRR=cov2cor(sampleCovMat) #Sample Correlation Matrix
  lambdaHatRR = eigen(hatRR)$values #eigenvalues of Sample Correlation Matrix
  
  lambdaCorrected = c() #Corrected eigenvalues by ACT
  for(j in 1:rmax){
    lambdaCorrected[j] = -1/under_m(n, j, lambdaHatRR[j], lambdaHatRR, p)
  }
  
  if( all((lambdaCorrected > (1 + sqrt(p/(n-1)))) == F) ){
    estFN_by_ACT[KKK] = 0
    estFN_by_all[KKK,13] = estFN_by_ACT[KKK]
  } else{
    estFN_by_ACT[KKK] = tail( which( lambdaCorrected > (1 + sqrt(p/(n-1))) ), n=1 )
    estFN_by_all[KKK,13] = estFN_by_ACT[KKK]
  } #estimated Factor Number is max j which lambdaCorrected > (1 + sqrt(p/(n-1))) satisfy.
  
  
  ##corpcor
  #hatRR_byCorpcor = cor.shrink(hatRR)
  #lambdaHatRR_byCorpcor = eigen(hatRR_byCorpcor)$values
  #cbind(lambdaHatRR[1:50], lambdaHatRR_byCorpcor[1:50])
  
  
  #### estimating using other methods
  {
      ## Method 2: PC1, PC2, PC3, IC1, IC2, IC3  
      ## ’‚∂Œcode «∏˘æ›Onasti 2010µƒ¬€Œƒ÷––¥µƒ°£
      v<-rep(0,rmax)
      kfactor<-1:rmax
      bNT<-(N+T)/(N*T)
      cNT<-min(N,T)
      bev<-eigen((X%*%t(X))/T)$values
      for(k in 1:rmax){
        v[k]<-sum(bev[(k+1):N])
      }
      
      PC1<-v-v[rmax]*bNT*log(bNT)*kfactor
      PC2<-v+v[rmax]*bNT*log(cNT)*kfactor
      PC3<-v+v[rmax]*log(cNT)/cNT*kfactor
      
      IC1<-log(v)-bNT*log(bNT)*kfactor
      IC2<-log(v)+bNT*log(cNT)*kfactor
      IC3<-log(v)+log(cNT)/cNT*kfactor
      
      PC1f[KKK]<-which.min(PC1); estFN_by_all[KKK,1]=PC1f[KKK]
      PC2f[KKK]<-which.min(PC2); estFN_by_all[KKK,2]=PC2f[KKK]
      PC3f[KKK]<-which.min(PC3); estFN_by_all[KKK,3]=PC3f[KKK]
      
      IC1f[KKK]<-which.min(IC1); estFN_by_all[KKK,4]=IC1f[KKK]
      IC2f[KKK]<-which.min(IC2); estFN_by_all[KKK,5]=IC2f[KKK]
      IC3f[KKK]<-which.min(IC3); estFN_by_all[KKK,6]=IC3f[KKK]
      
      
      ## Method 3: the method of onatski
      oev<-eigen((X%*%t(X))/T)$values
      ow<-2^(2/3)/(2^(2/3)-1)
      #	ormax<-1.55*min(N^(2/5),T^(2/5))	
      delte1<-max(N^(-1/2),T^(-1/2))
      delte2<-0
      delte3<-max(N^(-2/3),T^(-2/3))
      ou<-ow*oev[rmax+1]+(1-ow)*oev[2*rmax+1]
      ON1f[KKK]<-sum(ifelse(oev>(1+delte1)*ou,1,0)); estFN_by_all[KKK,7]=ON1f[KKK]
      ON2f[KKK]<-sum(ifelse(oev>(1+delte2)*ou,1,0)); estFN_by_all[KKK,8]=ON2f[KKK]
      ON3f[KKK]<-sum(ifelse(oev>(1+delte3)*ou,1,0)); estFN_by_all[KKK,9]=ON3f[KKK]
      
      ## Method 4: the method of onatski2
      nols<-4
      noev<-eigen((X%*%t(X))/T)$values
      oJ<-rep(1,(nols+1))
      ed<-noev[1:rmax]-noev[2:(rmax+1)]
      noj<-rmax+1
      for(j in 1:4){
        noy<-noev[noj:(noj+nols)]
        nox<-(noj+seq(-1,(nols-1),1))^(2/3)
        nobeta<-sum((nox-mean(nox))*(noy-mean(noy)))/sum((nox-mean(nox))^2)
        nodelte<-2*abs(nobeta)
        noer<-ifelse(max(ed-nodelte)<0,0,max(which(ed>=nodelte)))
        noj<-noer+1
      }
      NONf[KKK]<-noer; estFN_by_all[KKK,10]=NONf[KKK]
      
      ## Method 5: the method of horenstein
      hev<-eigen((X%*%t(X))/(T*N))$values
      er1<-hev[1:rmax]
      er2<-hev[2:(rmax+1)]
      gr1<-hev%*%hm1
      gr2<-hev%*%hm2
      gr3<-hev%*%hm3
      er<-er1/er2
      gr<-log(gr1/gr2)/log(gr2/gr3)
      ERf[KKK]<-which.max(er); estFN_by_all[KKK,11]=ERf[KKK]
      GRf[KKK]<-which.max(gr); estFN_by_all[KKK,12]=GRf[KKK]
  }
  
  ##: mean of K 
  brate1N<-mean(PC1f[1:KKK])
  brate2N<-mean(PC2f[1:KKK])
  brate3N<-mean(PC3f[1:KKK])
  brate4N<-mean(IC1f[1:KKK])
  brate5N<-mean(IC2f[1:KKK])
  brate6N<-mean(IC3f[1:KKK])
  
  orate1N<-mean(ON1f[1:KKK])
  orate2N<-mean(ON2f[1:KKK])
  orate3N<-mean(ON3f[1:KKK])
  
  norateN<-mean(NONf[1:KKK])
  hrate1N<-mean(ERf[1:KKK])
  hrate2N<-mean(GRf[1:KKK])
  zrateN<-mean(estFN_by_ACT[1:KKK])
  
  ##: accuracy of selection oder K 
  brate1<-sum(ifelse(PC1f==r,1,0))/KKK
  brate2<-sum(ifelse(PC2f==r,1,0))/KKK
  brate3<-sum(ifelse(PC3f==r,1,0))/KKK
  brate4<-sum(ifelse(IC1f==r,1,0))/KKK
  brate5<-sum(ifelse(IC2f==r,1,0))/KKK
  brate6<-sum(ifelse(IC3f==r,1,0))/KKK
  
  orate1<-sum(ifelse(ON1f==r,1,0))/KKK
  orate2<-sum(ifelse(ON2f==r,1,0))/KKK
  orate3<-sum(ifelse(ON3f==r,1,0))/KKK
  
  norate<-sum(ifelse(NONf==r,1,0))/KKK
  hrate1<-sum(ifelse(ERf==r,1,0))/KKK
  hrate2<-sum(ifelse(GRf==r,1,0))/KKK
  zrate<-sum(ifelse(estFN_by_ACT==r,1,0))/KKK
  
  ##  
  brate1U<-sum(ifelse(PC1f>r,1,0))/KKK
  brate2U<-sum(ifelse(PC2f>r,1,0))/KKK
  brate3U<-sum(ifelse(PC3f>r,1,0))/KKK
  brate4U<-sum(ifelse(IC1f>r,1,0))/KKK
  brate5U<-sum(ifelse(IC2f>r,1,0))/KKK
  brate6U<-sum(ifelse(IC3f>r,1,0))/KKK
  
  orate1U<-sum(ifelse(ON1f>r,1,0))/KKK
  orate2U<-sum(ifelse(ON2f>r,1,0))/KKK
  orate3U<-sum(ifelse(ON3f>r,1,0))/KKK
  
  norateU<-sum(ifelse(NONf>r,1,0))/KKK
  hrate1U<-sum(ifelse(ERf>r,1,0))/KKK
  hrate2U<-sum(ifelse(GRf>r,1,0))/KKK
  zrateU<-sum(ifelse(estFN_by_ACT>r,1,0))/KKK
  
  ##  
  brate1L<-sum(ifelse(PC1f[1:KKK]<r,1,0))/KKK
  brate2L<-sum(ifelse(PC2f[1:KKK]<r,1,0))/KKK
  brate3L<-sum(ifelse(PC3f[1:KKK]<r,1,0))/KKK
  brate4L<-sum(ifelse(IC1f[1:KKK]<r,1,0))/KKK
  brate5L<-sum(ifelse(IC2f[1:KKK]<r,1,0))/KKK
  brate6L<-sum(ifelse(IC3f[1:KKK]<r,1,0))/KKK
  
  orate1L<-sum(ifelse(ON1f[1:KKK]<r,1,0))/KKK
  orate2L<-sum(ifelse(ON2f[1:KKK]<r,1,0))/KKK
  orate3L<-sum(ifelse(ON3f[1:KKK]<r,1,0))/KKK
  
  norateL<-sum(ifelse(NONf[1:KKK]<r,1,0))/KKK
  hrate1L<-sum(ifelse(ERf[1:KKK]<r,1,0))/KKK
  hrate2L<-sum(ifelse(GRf[1:KKK]<r,1,0))/KKK
  zrateL<-sum(ifelse(estFN_by_ACT[1:KKK]<r,1,0))/KKK
  
  rate<-c(brate1,brate2,brate3,brate4,brate5,brate6,orate1,orate2,orate3,norate,hrate1,hrate2,zrate)
  #          PC1    PC2     PC3      IC1    IC2    IC3       ON1    ON2     ON3      NON    ER      GR           ours  
  rateU<-c(brate1U,brate2U,brate3U,brate4U,brate5U,brate6U,orate1U,orate2U,orate3U,norateU,hrate1U,hrate2U,zrateU)
  rateL<-c(brate1L,brate2L,brate3L,brate4L,brate5L,brate6L,orate1L,orate2L,orate3L,norateL,hrate1L,hrate2L,zrateL)
  rateN<-c(brate1N,brate2N,brate3N,brate4N,brate5N,brate6N,orate1N,orate2N,orate3N,norateN,hrate1N,hrate2N,zrateN)
  
  Mt1=(c(round(rate,  digits=3)*100)) 
  Mt2=(c(round(rateU, digits=3)*100)) 
  Mt3=(c(round(1-rate-rateU, digits=3)*100)) 
  Mt4=(c(round(rateN, digits=2))) 
  LS=c(3,6,8,11,12,13) 
  print(KKK)
  print(t(cbind(Mt1[LS], Mt2[LS], Mt3[LS], Mt4[LS])))
  
  
}
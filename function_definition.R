
library(quantreg)
## compute Fisher's z-transform values and their covariance matrix

cov.fisher.z<-function(X){
  ## X: n *p matrix 
  dimx<-dim(X)
  n<-dimx[1]; p<-dimx[2]
  
  ## number of tests: N
  N<-p*(p-1)/2
  
  ## sample correlation matrix 
  rho<-cor(X)
  
  ## psi matrix & cov.z matrix : N*N
  psi<-matrix(0,N,N)
  cov.z<-diag(1,N)
  
  ## all combinations of elements 1:p taken two: p*(p-1)/2
  idx<-t(combn(1:p,2))
  
  ## Fisher's z transform
  rho.vec<-rho[lower.tri(rho)]
  Z<-0.5*log( (1+rho.vec)/(1-rho.vec))
  
  for (Ni in 1:(N-1)){
    ## j,k,h,m notations
    j<-idx[Ni,1]; k<-idx[Ni,2]     
    
    for (Nj in (Ni+1):N){
      
      ## j,k,h,m notations
      h<-idx[Nj,1]; m<-idx[Nj,2] 
      
      psi[Ni,Nj]<- 0.5* (
        (rho[j,h] -rho[j,k]*rho[k,h])*(rho[k,m] -rho[k,h]*rho[h,m])+
          (rho[j,m] -rho[j,h]*rho[h,m])*(rho[k,h] -rho[k,j]*rho[j,h])+ 
          (rho[j,h] -rho[j,m]*rho[m,h])*(rho[k,m] -rho[k,j]*rho[j,m])+ 
          (rho[j,m] -rho[j,k]*rho[k,m])*(rho[k,h] -rho[k,m]*rho[m,h])
      )
      cov.z[Ni,Nj]<-psi[Ni,Nj]/(1-rho[j,k]^2)/(1-rho[h,m]^2)
      
    }
  }
  
  cov.z<-cov.z/(n-3)
  cov.z<-cov.z+t(cov.z)
  diag(cov.z)<-1/(n-3)
  
  return(list(Z=Z, cov.z=cov.z, psi=psi))    
}



corr2.fisher.z<-function(X1, X2){
  n1<-dim(X1)[1]
  n2<-dim(X2)[1]
  
  
  temp1<-cov.fisher.z(X1)
  temp2<-cov.fisher.z(X2)
  
  
  a<-1/sqrt(1/(n1-3)+1/(n2-3))
  Z<-a*(temp1$Z-temp2$Z)
  
  Sigma<-a^2*(temp1$cov.z+temp2$cov.z)
  
  return(list(Z=Z,Sigma=Sigma))
  
}

comp_L1sum<-function(M){
  sum(abs(M))
}


FanFDP<-function(stat,corM,ecut=0.01, alpha = NULL, num.fix=NULL){
  
  len<-length(stat)
  eig_full<-eigen(corM)
  eigval<-eig_full$values
  
  # Fan's condition : this is bad for weak corr matrix
  #cutid<-len-max(which(sqrt(cumsum((rev(eigval))^2))/sum(eigval)<ecut)) 
  
  # more accurate condition
  
  cutid<-0
  L1sum<-comp_L1sum(corM)/len^2;
  if (L1sum<=ecut) {cutid<-0}
  
  if (is.null(num.fix)){
    if (L1sum>ecut) {M<-corM}
    while(L1sum>ecut){
      cutid<-cutid+1
      M<-M-eigval[cutid]*(eig_full$vectors[,cutid]%*%t(eig_full$vectors[,cutid]))
      L1sum<-comp_L1sum(M)/len^2
    }
  }
  if (is.null(num.fix)==FALSE){cutid<-num.fix}
  
  if (cutid==0){avec<-rep(1,len);eta<-rep(0,len)}
  
  if (cutid>0){
    
    L<- t(t(eig_full$vectors[,1:cutid])*c(sqrt(eigval[1:cutid])));
    
    A<-corM-L%*%t(L)
    
    ### computing W
    
    cut_zero<-quantile(abs(stat),0.90)
    zeroid<-which(abs(stat)<cut_zero)
    DM<-cbind(stat[zeroid],L[zeroid,])
    rq_res<-rq(DM[,1]~-1+DM[,-1])
    W<-coef(rq_res)
    avec<-apply(L,1,function(x) 1/sqrt(1-sum(x^2)))
    eta<-as.vector(L%*%W)
  }
  
  if (is.null(alpha)) {alpha<-10^(seq(-10, -1,,100))}
  
  
  nsim<-length(alpha); 
  mu<-stat
  FDP<-rep(0,nsim)
  RR<-rep(0,nsim)
  R.mat<-V.mat<-S.mat<-V0.mat<-S0.mat<-matrix(0,len, nsim) ## 
  
  for (j in 1:nsim){
    
    quant<-qnorm(alpha[j]/2)
    
    term1<-pnorm(avec*(quant+eta))+pnorm(avec*(quant-eta))
    term1.mu<-pnorm(avec*(quant+eta+mu))+pnorm(avec*(quant-eta-mu))
    
    
    term0<-2*pnorm(quant)
    term0.mu<-pnorm(quant-mu)+pnorm(quant+mu)
    
    
    V.mat[,j]<-term1; V0.mat[,j]<-term0
    S.mat[,j]<-term1.mu; S0.mat[,j]<-term0.mu
    
    R.mat[,j]<-abs(stat)>qnorm(1-alpha[j]/2)
    RR[j]<-sum(R.mat[,j])
    
    FDP[j]<-sum(V.mat[,j])/max(RR[j],1)
    
  }
  
  FDP<-rev(cummin(rev(FDP)))
  V<-FDP*RR
  #plot(alpha,FDP,ylim=c(0,1))
  
  return(list(alpha=alpha,V.mat=V.mat, S.mat=S.mat, 
              V0.mat=V0.mat, S0.mat=S0.mat, R.mat=R.mat, 
              V=V, R=RR, FDP=FDP))
  
}

FanFDP.control<-function(alpha, FDP, level=0.05){
  idx<-sum(FDP<level)  
  alpha.optimal<-alpha[idx]+(level- FDP[idx])/(FDP[idx+1]-FDP[idx])*(alpha[idx+1]-alpha[idx]) 
  return(alpha.optimal)
}

#pca버전 Fan 알고리즘
FanFDP_pca<-function(stat,corM,ecut=0.01, alpha = NULL, num.fix=NULL){
  
  len<-length(stat) #1225
  eig_full<-eigen(corM) 
  eigval<-eig_full$values  #eigenvalues 1225개. decresing order.
  
  
  
  # Fan's condition : this is bad for weak corr matrix
  #cutid<-len-max(which(sqrt(cumsum((rev(eigval))^2))/sum(eigval)<ecut)) 
  
  # more accurate condition
  
  #cutid<-0
  #L1sum<-comp_L1sum(corM)/len^2;
  #if (L1sum<=ecut) {cutid<-0}
  
  if (is.null(num.fix)){
    var_exp<-cumsum(eigval)/sum(eigval) #전체 분산에서 몇퍼센트 설명하는지
    cutid<-head(which(var_exp>0.8),1) #처음으로 설명력 80% 되는 지점
    
  }
  
  if (is.null(num.fix)==FALSE){cutid<-num.fix}
  
  #if (cutid==0){avec<-rep(1,len);eta<-rep(0,len)}
  
  if (cutid>0){
    
    L<- t(t(eig_full$vectors[,1:cutid])*c(sqrt(eigval[1:cutid]))); #차원축소된 L
    
    A<-corM-L%*%t(L)  
    
    ### computing W
    
    cut_zero<-quantile(abs(stat),0.90)
    zeroid<-which(abs(stat)<cut_zero)
    DM<-cbind(stat[zeroid],L[zeroid,])
    rq_res<-rq(DM[,1]~-1+DM[,-1])
    W<-coef(rq_res)
    avec<-apply(L,1,function(x) 1/sqrt(1-sum(x^2)))
    eta<-as.vector(L%*%W)
  }
  
  if (is.null(alpha)) {alpha<-10^(seq(-10, -1,length = 100))}
  
  
  nsim<-length(alpha); 
  mu<-stat
  FDP<-rep(0,nsim)
  RR<-rep(0,nsim)
  R.mat<-V.mat<-S.mat<-V0.mat<-S0.mat<-matrix(0,len, nsim) ## 
  
  for (j in 1:nsim){
    
    quant<-qnorm(alpha[j]/2)
    
    term1<-pnorm(avec*(quant+eta))+pnorm(avec*(quant-eta))
    term1.mu<-pnorm(avec*(quant+eta+mu))+pnorm(avec*(quant-eta-mu))
    
    
    term0<-2*pnorm(quant)
    term0.mu<-pnorm(quant-mu)+pnorm(quant+mu)
    
    
    V.mat[,j]<-term1; V0.mat[,j]<-term0
    S.mat[,j]<-term1.mu; S0.mat[,j]<-term0.mu
    
    R.mat[,j]<-abs(stat)>qnorm(1-alpha[j]/2)
    RR[j]<-sum(R.mat[,j])
    
    FDP[j]<-sum(V.mat[,j])/max(RR[j],1)
    
  }
  
  FDP<-rev(cummin(rev(FDP)))
  V<-FDP*RR
  #plot(alpha,FDP,ylim=c(0,1))
  
  return(list(alpha=alpha,V.mat=V.mat, S.mat=S.mat, 
              V0.mat=V0.mat, S0.mat=S0.mat, R.mat=R.mat, 
              V=V, R=RR, FDP=FDP))
  
}



#논문버전 Fan
FanFDP_t<-function(stat,corM,ecut=0.01, alpha = NULL, num.fix=NULL){
  
  len<-length(stat)
  eig_full<-eigen(corM)
  eigval<-eig_full$values
  
  # Fan's condition : this is bad for weak corr matrix
  cutid<-len-max(which(sqrt(cumsum((rev(eigval))^2))/sum(eigval)<ecut)) 
  
  # more accurate condition
  
  if (is.null(num.fix)==FALSE){cutid<-num.fix}
  
  if (cutid==0){avec<-rep(1,len);eta<-rep(0,len)}
  
  if (cutid>0){
    
    L<- t(t(eig_full$vectors[,1:cutid])*c(sqrt(eigval[1:cutid])));
    
    A<-corM-L%*%t(L)
    
    ### computing W
    
    cut_zero<-quantile(abs(stat),0.90)
    zeroid<-which(abs(stat)<cut_zero)
    DM<-cbind(stat[zeroid],L[zeroid,])
    rq_res<-rq(DM[,1]~-1+DM[,-1])
    W<-coef(rq_res)
    avec<-apply(L,1,function(x) 1/sqrt(1-sum(x^2)))
    eta<-as.vector(L%*%W)
  }
  
  if (is.null(alpha)) {alpha<-10^(seq(-10, -1,length = 100))}
  
  
  nsim<-length(alpha); 
  mu<-stat
  FDP<-rep(0,nsim)
  RR<-rep(0,nsim)
  R.mat<-V.mat<-S.mat<-V0.mat<-S0.mat<-matrix(0,len, nsim) ## 
  
  for (j in 1:nsim){
    
    quant<-qnorm(alpha[j]/2)
    
    term1<-pnorm(avec*(quant+eta))+pnorm(avec*(quant-eta))
    term1.mu<-pnorm(avec*(quant+eta+mu))+pnorm(avec*(quant-eta-mu))
    
    
    term0<-2*pnorm(quant)
    term0.mu<-pnorm(quant-mu)+pnorm(quant+mu)
    
    
    V.mat[,j]<-term1; V0.mat[,j]<-term0
    S.mat[,j]<-term1.mu; S0.mat[,j]<-term0.mu
    
    R.mat[,j]<-abs(stat)>qnorm(1-alpha[j]/2)
    RR[j]<-sum(R.mat[,j])
    
    FDP[j]<-sum(V.mat[,j])/max(RR[j],1)
    
  }
  
  FDP<-rev(cummin(rev(FDP)))
  V<-FDP*RR
  #plot(alpha,FDP,ylim=c(0,1))
  
  return(list(alpha=alpha,V.mat=V.mat, S.mat=S.mat, 
              V0.mat=V0.mat, S0.mat=S0.mat, R.mat=R.mat, 
              V=V, R=RR, FDP=FDP))
  
}


## Toni Cai's procedure
## input : Y is the data matrix (variable by sample)
##  g : group indicator



TCSstat<-function(Y,g,boot=type,nboots=nboots){ ## boot=0: only observed statistics, boot=1: perform bootstrap
  
  p<-dim(Y)[2]
  n<-dim(Y)[1]
  
  tabg<-table(g)
  
  gx1<-which(g==names(tabg)[1])
  gx2<-which(g==names(tabg)[2])
  
  nn1<-length(gx1)
  nn2<-length(gx2)
  
  #  COR1<-t(Y[gx1,])%*%Y[gx1,]/(nn1-1) ## some elements >1 or < -1
  #  COR2<-t(Y[gx2,])%*%Y[gx2,]/(nn2-1)
  
  COR1<-cor(Y[gx1,]) #####  변경 
  COR2<-cor(Y[gx2,]) #####  변경 
  
  vech1<-COR1[upper.tri(COR1)]  
  vech2<-COR2[upper.tri(COR2)]  
  
  vech_diff=vech1-vech2
  
  
  ave1<-apply(Y[gx1,],2,mean) 
  kappa1<-sum(apply((Y[gx1,]-ave1)^4,2,sum) * nn1 / ( apply((Y[gx1,]-ave1)^2,2,sum))^2)/(3*p) ###
  
  ave2<-apply(Y[gx2,],2,mean) 
  kappa2<-sum(apply((Y[gx2,]-ave2)^4,2,sum)*nn2/(apply((Y[gx2,]-ave2)^2,2,sum))^2)/(3*p)
  
  
  #tt<-sqrt((nn1-3)*(nn2-3)/(nn1+nn2-6))*(0.5*log((1+vech1)/(1-vech1))-0.5*log((1+vech2)/(1-vech2)))
  
  rhot1<-vech1*(abs(vech1)/sqrt(kappa1*(1-vech1^2)^2/nn1)>2*sqrt(log(p)/nn1))
  rhot2<-vech2*(abs(vech2)/sqrt(kappa2*(1-vech2^2)^2/nn2)>2*sqrt(log(p)/nn2))
  rhotilde2<-pmax(rhot1^2,rhot2^2)
  
  tt<-(vech1-vech2)/sqrt(kappa1*(1-rhotilde2)^2/nn1+kappa2*(1-rhotilde2)^2/nn2)
  
  vech_diff=vech1-vech2
  
  if (boot==1) {
    
    tt_boot<-matrix(0,nboots,length(vech_diff))
    
    
    for (k in 1:nboots){
      
      YY1<-Y[sample(gx1,replace=TRUE),]
      YY2<-Y[sample(gx2,replace=TRUE),]
      Ystar<-rbind(YY1,YY2)
      
      COR1<-cor(Ystar[gx1,])
      COR2<-cor(Ystar[gx2,])
      
      vech1<-COR1[upper.tri(COR1)]  
      vech2<-COR2[upper.tri(COR2)]  
      
      
      ave1<-apply(Ystar[gx1,],2,mean)
      kappa1<-sum(apply((Ystar[gx1,]-ave1)^4,2,sum)*nn1/(apply((Ystar[gx1,]-ave1)^2,2,sum))^2)/(3*p)
      
      ave2<-apply(Ystar[gx2,],2,mean)
      kappa2<-sum(apply((Ystar[gx2,]-ave2)^4,2,sum)*nn2/(apply((Ystar[gx2,]-ave2)^2,2,sum))^2)/(3*p)
      
      
      #tt<-sqrt((nn1-3)*(nn2-3)/(nn1+nn2-6))*(0.5*log((1+vech1)/(1-vech1))-0.5*log((1+vech2)/(1-vech2)))
      
      rhot1<-vech1*(abs(vech1)/sqrt(kappa1*(1-vech1^2)^2/nn1)>2*sqrt(log(p)/nn1))
      rhot2<-vech2*(abs(vech2)/sqrt(kappa2*(1-vech2^2)^2/nn2)>2*sqrt(log(p)/nn2))
      rhotilde2<-pmax(rhot1^2,rhot2^2)
      
      tt_boot[k,]<-(vech1-vech2-vech_diff)/sqrt(kappa1*(1-rhotilde2)^2/nn1+kappa2*(1-rhotilde2)^2/nn2)
      
    }
    
  }
  
  
  if (boot==0) {return(tt)}
  if (boot==1) {return(list(tt=tt,tt_boot=tt_boot))}
  
}



TCSstat_Fisher<-function(Y,g,boot=type,nboots=nboots){ ## boot=0: only observed statistics, boot=1: perform bootstrap
  
  
  p<-dim(Y)[2]
  n<-dim(Y)[1]
  
  tabg<-table(g)
  
  gx1<-which(g==names(tabg)[1]) 
  gx2<-which(g==names(tabg)[2])  
  
  nn1<-length(gx1)
  nn2<-length(gx2)
  
  #  COR1<-t(Y[gx1,])%*%Y[gx1,]/(nn1-1) ## some elements >1 or < -1
  #  COR2<-t(Y[gx2,])%*%Y[gx2,]/(nn2-1)
  
  COR1<-cor(Y[gx1,]) #####  변경 
  COR2<-cor(Y[gx2,]) #####  변경
  
  vech1<-COR1[upper.tri(COR1)]  
  vech2<-COR2[upper.tri(COR2)]  

  vech_diff=vech1-vech2
  
  tt<-sqrt((nn1-3)*(nn2-3)/(nn1+nn2-6))*(0.5*log((1+vech1)/(1-vech1))-0.5*log((1+vech2)/(1-vech2)))
  
  
  
  if (boot==1) {
    
    tt_boot<-matrix(0,nboots,length(vech_diff))
    
    
    for (k in 1:nboots){
      
      YY1<-Y[sample(gx1,replace=TRUE),]
      YY2<-Y[sample(gx2,replace=TRUE),]
      Ystar<-rbind(YY1,YY2)
      
      COR1<-cor(Ystar[gx1,])
      COR2<-cor(Ystar[gx2,])
      
      vech1<-COR1[upper.tri(COR1)]  
      vech2<-COR2[upper.tri(COR2)]  
      
      tt_boot[k,]<-sqrt((nn1-3)*(nn2-3)/(nn1+nn2-6))*(0.5*log((1+vech1)/(1-vech1))-0.5*log((1+vech2)/(1-vech2)))
      
    }
    
  }
  
  
  if (boot==0) {return(tt)}
  if (boot==1) {return(list(tt=tt,tt_boot=tt_boot))}
  
}


TCS<-function(Y,g,true=true,alpha=0.05,type=0,nboots=1000,emp=0,FisherZ=0){  ## alpha= FDR value to be controlled, # type=0: normal approximation, type 1: bootstrap
  ## emp = 0 : standard normal null, emp =1 : empirical normal null
  
  if (emp==1) {library(locfdr)}
  
  p<-dim(Y)[2] 
  tabg<-table(g)
  
  gx1<-which(g==names(tabg)[1])
  gx2<-which(g==names(tabg)[2])
  
  nn1<-length(gx1)
  nn2<-length(gx2)
  
  #  COR1<-t(Y[gx1,])%*%Y[gx1,]/(nn1-1) ## some elements >1 or < -1
  #  COR2<-t(Y[gx2,])%*%Y[gx2,]/(nn2-1)
  
  
  COR1<-cor(Y[gx1,])
  COR2<-cor(Y[gx2,]) 
  
  if (type==0 && FisherZ==0) {tt<-TCSstat(Y,g,type,nboots)}
  if (type==1 && FisherZ==0) {RR<-TCSstat(Y,g,type,nboots);tt<-RR$tt;tt_boot<-RR$tt_boot}
  
  if (type==0 && FisherZ==1) {tt<-TCSstat_Fisher(Y,g,type,nboots)}
  if (type==1 && FisherZ==1) {RR<-TCSstat_Fisher(Y,g,type,nboots);tt<-RR$tt;tt_boot<-RR$tt_boot}
  
  ## Now, FDR
  
  if (type==0){
    xx<-seq(0,sqrt(4*log(p)-2*log(log(p))),length=100)
    tpos<-xx
    
    for (i in 1:100){
      
      if (emp==0) {numer<-(2-2*pnorm(xx[i]))*p*(p-1)*0.5}
      if (emp==1) {empnull<-locfdr(tt,plot=0);  mle<-empnull$fp0[3,] ; mu<-mle[1]; sigma<-mle[2]; phi0<-mle[3]; numer<-(2-2*pnorm((xx[i]-mu)/sigma))*p*(p-1)*0.5 }
      denom<-  max(sum(abs(tt)>=xx[i]),1)
      ratio<-numer/denom
      
      tpos[i]<-(ratio<=alpha)
    }
    
    
    if (max(tpos)==1) {threshold<-xx[min(which(tpos==1))]}
    if (max(tpos)==0) {threshold<-sqrt(4*log(p))}
    
    numrej<-sum(abs(tt)>threshold)
    numrej_true<-sum(abs(tt)>threshold & true==0)  ######### V(t):차이가 없는데 기각된거. 즉 잘못 기각된거 
    notreject_False<-sum(abs(tt)<threshold & true==1) #### 차이가 있는데, 기각되지 않은거.
    
  }
  
  
  if (type==1){
    
    xx<-seq(0,sqrt(4*log(p)-2*log(log(p))),length=100)
    tpos<-xx
    
    for (i in 1:100){
      
      numer<-(2*sum(tt_boot>=xx[i])/(p*(p-1)*nboots))*p*(p-1)*0.5
      
      denom<-  max(sum(abs(tt)>=xx[i]),1)
      ratio<-numer/denom
      
      tpos[i]<-(ratio<=alpha)
    }
    
    
    if (max(tpos)==1) {threshold<-xx[min(which(tpos==1))]}
    if (max(tpos)==0) {threshold<-sqrt(4*log(p))}
    
    numrej<-sum(abs(tt)>threshold)
    numrej_true<-sum(abs(tt)>threshold & true==0)  ######### 수정 
    notreject_False<-sum(abs(tt)<threshold & true==1) #### 차이가 있는데, 기각되지 않은거.
    
  }
  
  delta<-as.numeric(abs(tt)>threshold)
  delta.mat<-diag(1,p)
  delta.mat[upper.tri(delta.mat)]<-delta
  for (i in 2:p){
    for (j in 1:i){
      delta.mat[i,j]<-delta.mat[j,i]
    }
  }
  
  return(list(tt=tt,threshold=threshold,numrej=numrej,numrej_true=numrej_true, notreject_False = notreject_False, delta=delta, delta.mat=delta.mat, COR1=COR1, COR2=COR2))
}

#true가 없는 경우 - 실제 데이터에 적용 

TCS<-function(Y,g,alpha=0.05,type=0,nboots=1000,emp=0,FisherZ=0){  ## alpha= FDR value to be controlled, # type=0: normal approximation, type 1: bootstrap
  ## emp = 0 : standard normal null, emp =1 : empirical normal null
  
  if (emp==1) {library(locfdr)}
  
  p<-dim(Y)[2] ####바꾼부분 
  tabg<-table(g)
  
  gx1<-which(g==names(tabg)[1])
  gx2<-which(g==names(tabg)[2])
  
  nn1<-length(gx1)
  nn2<-length(gx2)
  
  #  COR1<-t(Y[gx1,])%*%Y[gx1,]/(nn1-1) ## some elements >1 or < -1
  #  COR2<-t(Y[gx2,])%*%Y[gx2,]/(nn2-1)
  
  
  COR1<-cor(Y[gx1,]) ########### 바꾼부분  # 
  COR2<-cor(Y[gx2,]) ########### 바꾼부분 
  
  if (type==0 && FisherZ==0) {tt<-TCSstat(Y,g,type,nboots)}
  if (type==1 && FisherZ==0) {RR<-TCSstat(Y,g,type,nboots);tt<-RR$tt;tt_boot<-RR$tt_boot}
  
  if (type==0 && FisherZ==1) {tt<-TCSstat_Fisher(Y,g,type,nboots)}
  if (type==1 && FisherZ==1) {RR<-TCSstat_Fisher(Y,g,type,nboots);tt<-RR$tt;tt_boot<-RR$tt_boot}
  
  ## Now, FDR
  
  if (type==0){
    xx<-seq(0,sqrt(4*log(p)-2*log(log(p))),length=100)
    tpos<-xx
    
    for (i in 1:100){
      
      if (emp==0) {numer<-(2-2*pnorm(xx[i]))*p*(p-1)*0.5}
      if (emp==1) {empnull<-locfdr(tt,plot=0);  mle<-empnull$fp0[3,] ; mu<-mle[1]; sigma<-mle[2]; phi0<-mle[3]; numer<-(2-2*pnorm((xx[i]-mu)/sigma))*p*(p-1)*0.5 }
      denom<-  max(sum(abs(tt)>=xx[i]),1)
      ratio<-numer/denom
      
      tpos[i]<-(ratio<=alpha)
    }
    
    
    if (max(tpos)==1) {threshold<-xx[min(which(tpos==1))]}
    if (max(tpos)==0) {threshold<-sqrt(4*log(p))}
    
    numrej<-sum(abs(tt)>threshold)
#    numrej_true<-sum(abs(tt)>threshold & true==0)  ######### V(t):차이가 없는데 기각된거. 즉 잘못 기각된거 
#    notreject_False<-sum(abs(tt)<threshold & true==1) #### 차이가 있는데, 기각되지 않은거.
    
  }
  
  
  if (type==1){
    
    xx<-seq(0,sqrt(4*log(p)-2*log(log(p))),length=100)
    tpos<-xx
    
    for (i in 1:100){
      
      numer<-(2*sum(tt_boot>=xx[i])/(p*(p-1)*nboots))*p*(p-1)*0.5
      
      denom<-  max(sum(abs(tt)>=xx[i]),1)
      ratio<-numer/denom
      
      tpos[i]<-(ratio<=alpha)
    }
    
    
    if (max(tpos)==1) {threshold<-xx[min(which(tpos==1))]}
    if (max(tpos)==0) {threshold<-sqrt(4*log(p))}
    
    numrej<-sum(abs(tt)>threshold)
#    numrej_true<-sum(abs(tt)>threshold & true==0)  ######### 수정 
#    notreject_False<-sum(abs(tt)<threshold & true==1) #### 차이가 있는데, 기각되지 않은거.
    
  }
  
  delta<-as.numeric(abs(tt)>threshold)
  delta.mat<-diag(1,p)
  delta.mat[upper.tri(delta.mat)]<-delta
  for (i in 2:p){
    for (j in 1:i){
      delta.mat[i,j]<-delta.mat[j,i]
    }
  }
  
  return(list(tt=tt,threshold=threshold,numrej=numrej, delta=delta, delta.mat=delta.mat, COR1=COR1, COR2=COR2))
}

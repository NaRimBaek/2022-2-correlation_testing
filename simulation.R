

############################# case1: Equal Correlation ####################################
#서로 다른 ti들의 correlation 동일 
level<-0.05
rho<-0.5
n<-20
p<-50
sigma<- diag(1-rho, p, p)+rho 
## Alternatives 
eff<-2
p1<-1225
true<-c(rep(1,N-p1),rep(0,p1)) #두 그룹을 다르게 해서 차이가 나는 애를 앞으로 몰아넣기(200개) -> 30개 정도 
table(true)


########################### normal mixture distribution###########################
library(MASS)
nsim=150

R.cai = c()
R.cai_boot = c()
R.bh = c()
R.by = c()
R.proposed_procedure = c()

#cai
for (cnt in 1:nsim){
  U1<- runif(n, 0,1) ; U2<-runif(n,0,1)
  Z1<-mvrnorm(n, rep(0,p),sigma);Z2<-mvrnorm(n, rep(0,p),sigma)
  X<- U1*Z1; Y<- U2*Z2  # X: 그룹1, Y: 그룹2
  
  all=rbind(X,Y) #두 그룹 붙임 
  g=c(rep(0,20),rep(1,20)) #indicator 
  table(g)
  
  temp_nm = TCS(all,g,true,alpha=0.05,type=0,emp=0,FisherZ=0) #boostrap 없는 경우
  R.cai[cnt]<-temp_nm$numrej #기각된 개수
  
  temp_nm_boot = TCS(all,g,true,alpha=0.05,type=1,nboots=100,emp=0,FisherZ=0) #boostrap 있는 경우
  R.cai_boot[cnt]<-temp_nm_boot$numrej
  
  level=0.05
  Fz<-corr2.fisher.z(X, Y)
  stat=Fz$Z
  pval<-2*(1-pnorm(abs(stat)))
  
  p.BH<-p.adjust(pval,"BH")
  R.bh[cnt]<-sum(p.BH<level) #총 4개 기각 
  
  p.BY<-p.adjust(pval,"BY")
  R.by[cnt]<-sum(p.BY<level) #총 0 기각 
  

  #pca fan
  pval<-2*(1-pnorm(abs(stat)))
  res<-FanFDP_pca(stat=Fz$Z,corM=Fz$Sigma)
  alpha <-FanFDP.control(res$alpha, res$FDP, level = level)
  R.proposed_procedure[cnt]<-sum(pval<alpha) 
  
}

result=data.frame(R.cai, R.cai_boot,R.bh, R.by, R.proposed_procedure)
head(result)

#number of rejections Null hypothesis 
boxplot(result,ylim=c(0,N))
title(main="number of rejected Null hypothesis",sub="normal mixture distribution")

power = (N - result) / N
boxplot(power)
title(main="Power of Null hypothesis",sub="normal mixture distribution")


write.csv(result,"/Users/baeknarim/학회포스터발표_2022/correlation testing/case1/normal_mixture/case1_normalmixture_rejection", row.names = FALSE)

###################normal distribution#################

R.cai = c()
R.cai_boot = c()
R.bh = c()
R.by = c()
R.proposed_procedure = c()



for (cnt in 1:nsim){
  X<- mvrnorm(n,rep(0,p),sigma) ; Y<-mvrnorm(n,rep(0,p),sigma)
  
  all=rbind(X,Y) #두 그룹 붙임 
  g=c(rep(0,20),rep(1,20)) #indicator 
  table(g)

  temp_nm = TCS(all,g,true,alpha=0.05,type=0,emp=0,FisherZ=0) #boostrap 없는 경우
  R.cai[cnt]<-temp_nm$numrej #기각된 개수
  
  temp_nm_boot = TCS(all,g,true,alpha=0.05,type=1,nboots=100,emp=0,FisherZ=0) #boostrap 있는 경우
  R.cai_boot[cnt]<-temp_nm_boot$numrej
  
  level=0.05
  Fz<-corr2.fisher.z(X, Y)
  stat=Fz$Z
  pval<-2*(1-pnorm(abs(stat)))
  
  p.BH<-p.adjust(pval,"BH")
  R.bh[cnt]<-sum(p.BH<level) #총 4개 기각 
  
  p.BY<-p.adjust(pval,"BY")
  R.by[cnt]<-sum(p.BY<level) #총 0 기각 
  
  #pca fan
  pval<-2*(1-pnorm(abs(stat)))
  res<-FanFDP_pca(stat=Fz$Z,corM=Fz$Sigma)
  alpha <-FanFDP.control(res$alpha, res$FDP, level = level)
  
  R.proposed_procedure[cnt]<-sum(pval<alpha) 
}



result=data.frame(R.cai, R.cai_boot,R.bh, R.by,R.proposed_procedure )
head(result)


#number of rejections Null hypothesis 
boxplot(result,ylim=c(0,N))
title(main="number of rejected Null hypothesis",sub="normal distribution")

power = (N - result) / N
boxplot(power)
title(main="Power of Null hypothesis",sub="normal distribution")

write.csv(result,"/Users/baeknarim/학회포스터발표_2022/correlation testing/case1/normal/case1_normal_rejection", row.names = FALSE)


################### t-distribution#################
install.packages("LaplacesDemon")
library('LaplacesDemon')
R.cai = c()
R.cai_boot = c()
R.bh = c()
R.by = c()
R.proposed_procedure = c()


for (cnt in 1:nsim){
  X<-rmvt(n, rep(0,p), sigma, df=6) ;Y<-rmvt(n,rep(0,p),  sigma, df=6)

  all=rbind(X,Y) #두 그룹 붙임 
  g=c(rep(0,20),rep(1,20)) #indicator 
  table(g)
  
  temp_nm = TCS(all,g,true,alpha=0.05,type=0,emp=0,FisherZ=0) #boostrap 없는 경우
  R.cai[cnt]<-temp_nm$numrej #기각된 개수
  
  temp_nm_boot = TCS(all,g,true,alpha=0.05,type=1,nboots=100,emp=0,FisherZ=0) #boostrap 있는 경우
  R.cai_boot[cnt]<-temp_nm_boot$numrej
  
  level=0.05
  Fz<-corr2.fisher.z(X, Y)
  stat=Fz$Z
  pval<-2*(1-pnorm(abs(stat)))
  
  p.BH<-p.adjust(pval,"BH")
  R.bh[cnt]<-sum(p.BH<level) #총 4개 기각 
  
  p.BY<-p.adjust(pval,"BY")
  R.by[cnt]<-sum(p.BY<level) #총 0 기각 
  
  #pca fan
  pval<-2*(1-pnorm(abs(stat)))
  res<-FanFDP_pca(stat=Fz$Z,corM=Fz$Sigma)
  alpha <-FanFDP.control(res$alpha, res$FDP, level = level)
  R.proposed_procedure[cnt]<-sum(pval<alpha) 
  
}


result=data.frame(R.cai, R.cai_boot,R.bh, R.by,R.proposed_procedure)
head(result)

#number of rejections Null hypothesis 
boxplot(result,ylim=c(0,N))
title(main="number of rejected Null hypothesis",sub="t distribution")


power = (N - result) / N
boxplot(power)
title(main="Power of Null hypothesis",sub="t distribution")

write.csv(result,"/Users/baeknarim/학회포스터발표_2022/correlation testing/case1/t/case1_t_rejection", row.names = FALSE)



###################exponential distribution#################
library(LaplacesDemon)

R.cai = c()
R.cai_boot = c()
R.bh = c()
R.by = c()
R.proposed_procedure = c()


for (cnt in 1:nsim){
  X<-rmvpe(20, rep(0,p), sigma, kappa=6); Y<-rmvpe(20, rep(0,p), sigma, kappa=6)
  
  all=rbind(X,Y) #두 그룹 붙임 
  g=c(rep(0,20),rep(1,20)) #indicator 
  table(g)
  
  temp_nm = TCS(all,g,true,alpha=0.05,type=0,emp=0,FisherZ=0) #boostrap 없는 경우
  R.cai[cnt]<-temp_nm$numrej #기각된 개수
  
  temp_nm_boot = TCS(all,g,true,alpha=0.05,type=1,nboots=100,emp=0,FisherZ=0) #boostrap 있는 경우
  R.cai_boot[cnt]<-temp_nm_boot$numrej
  
  level=0.05
  Fz<-corr2.fisher.z(X, Y)
  stat=Fz$Z
  pval<-2*(1-pnorm(abs(stat)))
  
  p.BH<-p.adjust(pval,"BH")
  R.bh[cnt]<-sum(p.BH<level) #총 4개 기각 
  
  p.BY<-p.adjust(pval,"BY")
  R.by[cnt]<-sum(p.BY<level) #총 0 기각 
  
  #pca fan
  pval<-2*(1-pnorm(abs(stat)))
  res<-FanFDP_pca(stat=Fz$Z,corM=Fz$Sigma)
  alpha <-FanFDP.control(res$alpha, res$FDP, level = level)
  R.proposed_procedure[cnt]<-sum(pval<alpha) 
  
}


result=data.frame(R.cai, R.cai_boot,R.bh, R.by, R.proposed_procedure)
head(result)

#number of rejections Null hypothesis 
boxplot(result,ylim=c(0,N))
title(main="number of rejected Null hypothesis",sub="exponential distribution")


power = (N - result) / N
boxplot(power)
title(main="Power of Null hypothesis",sub="exponential distribution")

write.csv(result,"/Users/baeknarim/학회포스터발표_2022/correlation testing/case1/exp/case1_exponential_rejection", row.names = FALSE)



############################# case3: 두 집단 correlation 일부분 다름#############################
library(MASS)
library(lqmm)
library(matrixcalc)

level<-0.05
rho<-0.5
n<-20
p<-50

sigma1<- diag(1-rho, p, p)+rho #대각은 모두 1, 나머지는 0.5
sigma2<- diag(1-rho, p, p)+rho #대각은 모두 1, 나머지는 0.5

# 첫 4줄 corr 다르게. 49 + 48 + 47 + 46 = 190
difference1<-runif(49,0,0.1)  
difference2<-runif(49,0.8,0.9)
for (i in 1:49){
  sigma1[1,i+1] = difference1[i]
  sigma2[1,i+1] = difference2[i]
}

difference1<-runif(49,0,0.1)  
difference2<-runif(49,0.8,0.9)
for (i in 1:48){
  sigma1[2,i+2] = difference1[i]
  sigma2[2,i+2] = difference2[i]
}

difference1<-runif(49,0,0.1)  
difference2<-runif(49,0.8,0.9)
for (i in 1:47){
  sigma1[3,i+3] = difference1[i]
  sigma2[3,i+3] = difference2[i]
}

difference1<-runif(49,0,0.1)  
difference2<-runif(49,0.8,0.9)
for (i in 1:46){
  sigma1[4,i+4] = difference1[i]
  sigma2[4,i+4] = difference2[i]
}


sigma1[lower.tri(sigma1, diag=FALSE)]<-t(sigma1)[lower.tri(sigma1, diag=FALSE)]
sigma2[lower.tri(sigma2, diag=FALSE)]<-t(sigma2)[lower.tri(sigma2, diag=FALSE)]

#positive definite matrix 만들기
sigma1<-make.positive.definite(sigma1)
sigma2<-make.positive.definite(sigma2)
is.symmetric.matrix(sigma1)
is.symmetric.matrix(sigma2)

## Alternatives 
N<-1225
t=200
true<-c(rep(1,t),rep(0,N-t)) #두 그룹을 다르게 해서 차이가 나는 애를 앞으로 몰아넣기(190개) 
table(true)


######################### 1.normal mixture distribution  20*50#########################
nsim=150

R.cai = c()
V.cai =c() 
FDP.cai=c()
S.cai=c()
FNDP.cai=c()

R.cai_boot = c()
V.cai_boot =c() 
FDP.cai_boot=c()
S.cai_boot=c()
FNDP.cai_boot=c()

R.bh = c()
V.bh=c()
FDP.bh=c()
S.bh=c()
FNDP.bh=c()

R.by = c()
V.by=c()
FDP.by=c()
S.by=c()
FNDP.by=c()

R.proposed_procedure = c()
V.proposed_procedure=c()
FDP.proposed_procedure=c()
S.proposed_procedure=c()
FNDP.proposed_procedure=c()


#cai
for (cnt in 1:nsim){
  U1<- runif(n, 0,1) ; U2<-runif(n,0,1)
  Z1<-mvrnorm(n, rep(0,p),sigma1);Z2<-mvrnorm(n, rep(0,p),sigma2)
  X<- U1*Z1 ; Y<- U2*Z2  # X: 그룹1, Y: 그룹2
  
  all=rbind(X,Y) #두 그룹 붙임 
  g=c(rep(0,20),rep(1,20)) #indicator 
  table(g)
  
  #cai 
  temp_nm = TCS(all,g,true,alpha=0.05,type=0,emp=0,FisherZ=0) #boostrap 없는 경우
  R.cai[cnt]<-temp_nm$numrej #기각된 개수
  V.cai[cnt]<-temp_nm$numrej_true  #차이가 안나는데 잘못 기각되어버린 것들 
  FDP.cai[cnt]<-V.cai[cnt]/max(R.cai[cnt],1)
  S.cai[cnt] <-temp_nm$notreject_False
  FNDP.cai[cnt]<-S.cai[cnt]/ (N - max(R.cai[cnt],1) )
  
  #bootstrap 
  temp_nm_boot = TCS(all,g,true,alpha=0.05,type=1,nboots=100,emp=0,FisherZ=0) #boostrap 있는 경우
  R.cai_boot[cnt]<-temp_nm_boot$numrej #기각된 개수
  V.cai_boot[cnt]<-temp_nm_boot$numrej_true  #차이가 안나는데 잘못 기각되어버린 것들 
  FDP.cai_boot[cnt]<-V.cai_boot[cnt]/max(R.cai_boot[cnt],1)
  S.cai_boot[cnt] <-temp_nm_boot$notreject_False
  FNDP.cai_boot[cnt]<-S.cai_boot[cnt]/ (N - max(R.cai_boot[cnt],1) )
  
  level=0.05
  Fz<-corr2.fisher.z(X, Y)
  stat=Fz$Z
  pval<-2*(1-pnorm(abs(stat)))
  
  #bh
  p.BH<-p.adjust(pval,"BH")
  R.bh[cnt]<-sum(p.BH<level)   #기각된 개수
  V.bh[cnt]<-sum(p.BH<level & true==0) #차이가 안나는데 잘못 기각되어버린 것들 
  FDP.bh[cnt]<-V.bh[cnt]/max(R.bh[cnt],1)
  S.bh[cnt]<-sum(p.BH>=level & true==1) #차이가 나는데 기각 안된 것들 
  FNDP.bh[cnt]<- S.bh[cnt] / ( N- max(R.bh[cnt],1))
  
  #by
  p.BY<-p.adjust(pval,"BY")
  R.by[cnt]<-sum(p.BY<level) #총 0 기각 
  V.by[cnt]<-sum(p.BY<level & true==0)
  FDP.by[cnt]<-V.by[cnt]/max(R.by[cnt],1)
  S.by[cnt]<-sum(p.BY>=level & true==1) #차이가 나는데 기각 안된 것들 
  FNDP.by[cnt]<- S.by[cnt] / ( N- max(R.by[cnt],1))
  
  #pca fan
  res<-FanFDP_pca(stat=Fz$Z,corM=Fz$Sigma)
  alpha <-FanFDP.control(res$alpha, res$FDP, level = level)
  
  R.proposed_procedure[cnt]<-sum(pval<alpha) 
  V.proposed_procedure[cnt]<-sum(pval<alpha & true==0) # 기각이 되었는데 잘못 기각된거 (사실 차이가 없는애들)
  FDP.proposed_procedure[cnt]<-V.proposed_procedure[cnt]/max(R.proposed_procedure[cnt],1)
  S.proposed_procedure[cnt]<-sum(pval>=alpha & true==1) #차이가 나는데 기각 안된 것들 
  FNDP.proposed_procedure[cnt]<- S.proposed_procedure[cnt] / ( N- R.proposed_procedure[cnt])
  
}


result=data.frame(R.cai,V.cai,FDP.cai,FNDP.cai,
                  R.cai_boot,V.cai_boot,FDP.cai_boot,FNDP.cai_boot,
                  R.bh,V.bh,FDP.bh,FNDP.bh,
                  R.by,V.by,FDP.by,FNDP.by,
                  R.proposed_procedure, V.proposed_procedure,FDP.proposed_procedure,S.proposed_procedure,FNDP.proposed_procedure)

FDP = data.frame(FDP.cai,
                 FDP.cai_boot,
                 FDP.bh,
                 FDP.by,
                 FDP.proposed_procedure)

FNDP = data.frame(FNDP.cai,
                  FNDP.cai_boot,
                  FNDP.bh,
                  FNDP.by,
                  FNDP.proposed_procedure)
head(result)

#number of rejections Null hypothesis 

boxplot(FDP, ylim=c(0,1))
title(main="FDP",sub="normal mixture distribution")
abline(h=0.05,col="red",lty=3)

boxplot(FNDP, ylim=c(0,1))
title(main="FNDP",sub="normal mixture distribution")


write.csv(result,"/Users/baeknarim/학회포스터발표_2022/correlation testing/case3/normal_mixture/case3_normalmixture.csv", row.names = FALSE)



############################ 2.normal distribution###################### 
R.cai = c()
V.cai =c() 
FDP.cai=c()
S.cai=c()
FNDP.cai=c()

R.cai_boot = c()
V.cai_boot =c() 
FDP.cai_boot=c()
S.cai_boot=c()
FNDP.cai_boot=c()

R.bh = c()
V.bh=c()
FDP.bh=c()
S.bh=c()
FNDP.bh=c()

R.by = c()
V.by=c()
FDP.by=c()
S.by=c()
FNDP.by=c()

R.proposed_procedure = c()
V.proposed_procedure=c()
FDP.proposed_procedure=c()
S.proposed_procedure=c()
FNDP.proposed_procedure=c()



#cai
for (cnt in 1:nsim){
  X<- mvrnorm(n,rep(0,p),sigma1) ; Y<-mvrnorm(n,rep(0,p),sigma2)
  
  all=rbind(X,Y) #두 그룹 붙임 
  g=c(rep(0,20),rep(1,20)) #indicator 
  table(g)
  
  #cai
  temp_nm = TCS(all,g,true,alpha=0.05,type=0,emp=0,FisherZ=0) #boostrap 없는 경우
  R.cai[cnt]<-temp_nm$numrej #기각된 개수
  V.cai[cnt]<-temp_nm$numrej_true  #차이가 안나는데 잘못 기각되어버린 것들 
  FDP.cai[cnt]<-V.cai[cnt]/max(R.cai[cnt],1)
  S.cai[cnt] <-temp_nm$notreject_False
  FNDP.cai[cnt]<-S.cai[cnt]/ (N - max(R.cai[cnt],1) )
  
  #cai boostrap
  temp_nm_boot = TCS(all,g,true,alpha=0.05,type=1,nboots=100,emp=0,FisherZ=0) #boostrap 있는 경우 type1
  R.cai_boot[cnt]<-temp_nm_boot$numrej #기각된 개수
  V.cai_boot[cnt]<-temp_nm_boot$numrej_true  #차이가 안나는데 잘못 기각되어버린 것들 
  FDP.cai_boot[cnt]<-V.cai_boot[cnt]/max(R.cai_boot[cnt],1)
  S.cai_boot[cnt] <-temp_nm_boot$notreject_False
  FNDP.cai_boot[cnt]<-S.cai_boot[cnt]/ (N - max(R.cai_boot[cnt],1) )
  
  level=0.05
  Fz<-corr2.fisher.z(X, Y)
  stat=Fz$Z
  pval<-2*(1-pnorm(abs(stat)))
  
  #bh
  p.BH<-p.adjust(pval,"BH")
  R.bh[cnt]<-sum(p.BH<level)   #기각된 개수
  V.bh[cnt]<-sum(p.BH<level & true==0) #차이가 안나는데 잘못 기각되어버린 것들 
  FDP.bh[cnt]<-V.bh[cnt]/max(R.bh[cnt],1)
  S.bh[cnt]<-sum(p.BH>=level & true==1) #차이가 나는데 기각 안된 것들 
  FNDP.bh[cnt]<- S.bh[cnt] / ( N- max(R.bh[cnt],1))
  
  #by
  p.BY<-p.adjust(pval,"BY")
  R.by[cnt]<-sum(p.BY<level) #총 0 기각 
  V.by[cnt]<-sum(p.BY<level & true==0)
  FDP.by[cnt]<-V.by[cnt]/max(R.by[cnt],1)
  S.by[cnt]<-sum(p.BY>=level & true==1) #차이가 나는데 기각 안된 것들 
  FNDP.by[cnt]<- S.by[cnt] / ( N- max(R.by[cnt],1))
  
  #pca fan
  res<-FanFDP_pca(stat=Fz$Z,corM=Fz$Sigma)
  alpha <-FanFDP.control(res$alpha, res$FDP, level = level)
  
  R.proposed_procedure[cnt]<-sum(pval<alpha) 
  V.proposed_procedure[cnt]<-sum(pval<alpha & true==0) # 기각이 되었는데 잘못 기각된거 (사실 차이가 없는애들)
  FDP.proposed_procedure[cnt]<-V.proposed_procedure[cnt]/max(R.proposed_procedure[cnt],1)
  S.proposed_procedure[cnt]<-sum(pval>=alpha & true==1) #차이가 나는데 기각 안된 것들 
  FNDP.proposed_procedure[cnt]<- S.proposed_procedure[cnt] / ( N- R.proposed_procedure[cnt])

}


result=data.frame(R.cai,V.cai,FDP.cai,FNDP.cai,
                  R.cai_boot,V.cai_boot,FDP.cai_boot,FNDP.cai_boot,
                  R.bh,V.bh,FDP.bh,FNDP.bh,
                  R.by,V.by,FDP.by,FNDP.by,
                  R.proposed_procedure, V.proposed_procedure,FDP.proposed_procedure,S.proposed_procedure,FNDP.proposed_procedure)

FDP = data.frame(FDP.cai,
                 FDP.cai_boot,
                 FDP.bh,
                 FDP.by,
                 FDP.proposed_procedure)

FNDP = data.frame(FNDP.cai,
                  FNDP.cai_boot,
                  FNDP.bh,
                  FNDP.by,
                  FNDP.proposed_procedure)
head(result)

  

boxplot(FDP, ylim=c(0,1))
title(main="FDP",sub="normal distribution")
abline(h=0.05,col="red",lty=3)

boxplot(FNDP, ylim=c(0,1))
title(main="FNDP",sub="normal distribution")

write.csv(result,"/Users/baeknarim/학회포스터발표_2022/correlation testing/case3/normal/case3_normal.csv", row.names = FALSE)



############################3.t-distribution############################
library(LaplacesDemon)

R.cai = c()
V.cai =c() 
FDP.cai=c()
S.cai=c()
FNDP.cai=c()

R.cai_boot = c()
V.cai_boot =c() 
FDP.cai_boot=c()
S.cai_boot=c()
FNDP.cai_boot=c()

R.bh = c()
V.bh=c()
FDP.bh=c()
S.bh=c()
FNDP.bh=c()

R.by = c()
V.by=c()
FDP.by=c()
S.by=c()
FNDP.by=c()

R.proposed_procedure = c()
V.proposed_procedure=c()
FDP.proposed_procedure=c()
S.proposed_procedure=c()
FNDP.proposed_procedure=c()


#cai
for (cnt in 1:nsim){
  X<-rmvt(n, rep(0,p), sigma1, df=6); Y<-rmvt(n,rep(0,p),  sigma2, df=6)
  
  #cai
  all=rbind(X,Y) #두 그룹 붙임 
  g=c(rep(0,20),rep(1,20)) #indicator 
  table(g)
  
  temp_nm = TCS(all,g,true,alpha=0.05,type=0,emp=0,FisherZ=0) #boostrap 없는 경우
  R.cai[cnt]<-temp_nm$numrej #기각된 개수
  V.cai[cnt]<-temp_nm$numrej_true  #차이가 안나는데 잘못 기각되어버린 것들 
  FDP.cai[cnt]<-V.cai[cnt]/max(R.cai[cnt],1)
  S.cai[cnt] <-temp_nm$notreject_False
  FNDP.cai[cnt]<-S.cai[cnt]/ (N - max(R.cai[cnt],1) )
  
  #cai-boostrap
  temp_nm_boot = TCS(all,g,true,alpha=0.05,type=1,nboots=100,emp=0,FisherZ=0) #boostrap 있는 경우
  R.cai_boot[cnt]<-temp_nm_boot$numrej #기각된 개수
  V.cai_boot[cnt]<-temp_nm_boot$numrej_true  #차이가 안나는데 잘못 기각되어버린 것들 
  FDP.cai_boot[cnt]<-V.cai_boot[cnt]/max(R.cai_boot[cnt],1)
  S.cai_boot[cnt] <-temp_nm_boot$notreject_False
  FNDP.cai_boot[cnt]<-S.cai_boot[cnt]/ (N - max(R.cai_boot[cnt],1) )
  
  level=0.05
  Fz<-corr2.fisher.z(X, Y)
  stat=Fz$Z
  pval<-2*(1-pnorm(abs(stat)))
  
  #bh
  p.BH<-p.adjust(pval,"BH")
  R.bh[cnt]<-sum(p.BH<level)  
  V.bh[cnt]<-sum(p.BH<level & true==0)
  FDP.bh[cnt]<-V.bh[cnt]/max(R.bh[cnt],1)
  S.bh[cnt] <-sum(p.BH>=level & true==1)
  FNDP.bh[cnt]<-S.bh[cnt]/ (N - max(R.bh[cnt],1) )
  
  #by
  p.BY<-p.adjust(pval,"BY")
  R.by[cnt]<-sum(p.BY<level) #총 0 기각 
  V.by[cnt]<-sum(p.BY<level & true==0)
  FDP.by[cnt]<-V.by[cnt]/max(R.by[cnt],1)
  S.by[cnt]<-sum(p.BY>=level & true==1) #차이가 나는데 기각 안된 것들 
  FNDP.by[cnt]<- S.by[cnt] / ( N- max(R.by[cnt],1))
  
  #pca fan
  res<-FanFDP_pca(stat=Fz$Z,corM=Fz$Sigma)
  alpha <-FanFDP.control(res$alpha, res$FDP, level = level)
  
  R.proposed_procedure[cnt]<-sum(pval<alpha) 
  V.proposed_procedure[cnt]<-sum(pval<alpha & true==0) # 기각이 되었는데 잘못 기각된거 (사실 차이가 없는애들)
  FDP.proposed_procedure[cnt]<-V.proposed_procedure[cnt]/max(R.proposed_procedure[cnt],1)
  S.proposed_procedure[cnt]<-sum(pval>=alpha & true==1) #차이가 나는데 기각 안된 것들 
  FNDP.proposed_procedure[cnt]<- S.proposed_procedure[cnt] / ( N- R.proposed_procedure[cnt])
}


result=data.frame(R.cai,V.cai,FDP.cai,FNDP.cai,
                  R.cai_boot,V.cai_boot,FDP.cai_boot,FNDP.cai_boot,
                  R.bh,V.bh,FDP.bh,FNDP.bh,
                  R.by,V.by,FDP.by,FNDP.by,
                  R.proposed_procedure, V.proposed_procedure,FDP.proposed_procedure,S.proposed_procedure,FNDP.proposed_procedure)

FDP = data.frame(FDP.cai,
                 FDP.cai_boot,
                 FDP.bh,
                 FDP.by,
                 FDP.proposed_procedure)

FNDP = data.frame(FNDP.cai,
                  FNDP.cai_boot,
                  FNDP.bh,
                  FNDP.by,
                  FNDP.proposed_procedure)
head(result)


boxplot(FDP, ylim=c(0,1))
title(main="FDP",sub="t distribution")
abline(h=0.05,col="red",lty=3)

boxplot(FNDP, ylim=c(0,1))
title(main="FNDP",sub="t distribution")


write.csv(result,"/Users/baeknarim/학회포스터발표_2022/correlation testing/case3/t/case3_t.csv", row.names = FALSE)

##########################4.exponential distribution#############
R.cai = c()
V.cai =c() 
FDP.cai=c()
S.cai=c()
FNDP.cai=c()

R.cai_boot = c()
V.cai_boot =c() 
FDP.cai_boot=c()
S.cai_boot=c()
FNDP.cai_boot=c()

R.bh = c()
V.bh=c()
FDP.bh=c()
S.bh=c()
FNDP.bh=c()

R.by = c()
V.by=c()
FDP.by=c()
S.by=c()
FNDP.by=c()

R.proposed_procedure = c()
V.proposed_procedure=c()
FDP.proposed_procedure=c()
S.proposed_procedure=c()
FNDP.proposed_procedure=c()

for (cnt in 1:nsim){
  X<-rmvpe(20, rep(0,p), sigma1, kappa=6); Y<-rmvpe(20, rep(0,p), sigma2, kappa=6)
  
  #cai
  all=rbind(X,Y) #두 그룹 붙임 
  g=c(rep(0,20),rep(1,20)) #indicator 
  table(g)
  
  temp_nm = TCS(all,g,true,alpha=0.05,type=0,emp=0,FisherZ=0) #boostrap 없는 경우
  R.cai[cnt]<-temp_nm$numrej #기각된 개수
  V.cai[cnt]<-temp_nm$numrej_true  #차이가 안나는데 잘못 기각되어버린 것들 
  FDP.cai[cnt]<-V.cai[cnt]/max(R.cai[cnt],1)
  S.cai[cnt] <-temp_nm$notreject_False
  FNDP.cai[cnt]<-S.cai[cnt]/ (N - max(R.cai[cnt],1) )
  
  #cai-boostrap
  temp_nm_boot = TCS(all,g,true,alpha=0.05,type=1,nboots=100,emp=0,FisherZ=0) #boostrap 있는 경우
  R.cai_boot[cnt]<-temp_nm_boot$numrej #기각된 개수
  V.cai_boot[cnt]<-temp_nm_boot$numrej_true  #차이가 안나는데 잘못 기각되어버린 것들 
  FDP.cai_boot[cnt]<-V.cai_boot[cnt]/max(R.cai_boot[cnt],1)
  S.cai_boot[cnt] <-temp_nm_boot$notreject_False
  FNDP.cai_boot[cnt]<-S.cai_boot[cnt]/ (N - max(R.cai_boot[cnt],1) )
  
  level=0.05
  Fz<-corr2.fisher.z(X, Y)
  stat=Fz$Z
  pval<-2*(1-pnorm(abs(stat)))
  
  #bh
  p.BH<-p.adjust(pval,"BH")
  R.bh[cnt]<-sum(p.BH<level)  
  V.bh[cnt]<-sum(p.BH<level & true==0)
  FDP.bh[cnt]<-V.bh[cnt]/max(R.bh[cnt],1)
  S.bh[cnt] <-sum(p.BH>=level & true==1)
  FNDP.bh[cnt]<-S.bh[cnt]/ (N - max(R.bh[cnt],1) )
  
  #by
  p.BY<-p.adjust(pval,"BY")
  R.by[cnt]<-sum(p.BY<level) #총 0 기각 
  V.by[cnt]<-sum(p.BY<level & true==0)
  FDP.by[cnt]<-V.by[cnt]/max(R.by[cnt],1)
  S.by[cnt]<-sum(p.BY>=level & true==1) #차이가 나는데 기각 안된 것들 
  FNDP.by[cnt]<- S.by[cnt] / ( N- max(R.by[cnt],1))
  
  #pca fan
  res<-FanFDP_pca(stat=Fz$Z,corM=Fz$Sigma)
  alpha <-FanFDP.control(res$alpha, res$FDP, level = level)
  
  R.proposed_procedure[cnt]<-sum(pval<alpha) 
  V.proposed_procedure[cnt]<-sum(pval<alpha & true==0) # 기각이 되었는데 잘못 기각된거 (사실 차이가 없는애들)
  FDP.proposed_procedure[cnt]<-V.proposed_procedure[cnt]/max(R.proposed_procedure[cnt],1)
  S.proposed_procedure[cnt]<-sum(pval>=alpha & true==1) #차이가 나는데 기각 안된 것들 
  FNDP.proposed_procedure[cnt]<- S.proposed_procedure[cnt] / ( N- R.proposed_procedure[cnt])
  
}

result=data.frame(R.cai,V.cai,FDP.cai,FNDP.cai,
                  R.cai_boot,V.cai_boot,FDP.cai_boot,FNDP.cai_boot,
                  R.bh,V.bh,FDP.bh,FNDP.bh,
                  R.by,V.by,FDP.by,FNDP.by,
                  R.proposed_procedure, V.proposed_procedure,FDP.proposed_procedure,S.proposed_procedure,FNDP.proposed_procedure)

FDP = data.frame(FDP.cai,
                 FDP.cai_boot,
                 FDP.bh,
                 FDP.by,
                 FDP.proposed_procedure)

FNDP = data.frame(FNDP.cai,
                  FNDP.cai_boot,
                  FNDP.bh,
                  FNDP.by,
                  FNDP.proposed_procedure)
head(result)


boxplot(FDP, ylim=c(0,1))
title(main="FDP",sub="exponential distribution")
abline(h=0.05,col="red",lty=3)

boxplot(FNDP, ylim=c(0,1))
title(main="FNDP",sub="exponential distribution")


write.csv(result,"/Users/baeknarim/학회포스터발표_2022/correlation testing/case3/exp/case3_exp.csv", row.names = FALSE)



#plot
#case1
normal_mixture=read.csv("/Users/baeknarim/학회포스터발표_2022/correlation testing/case1/normal_mixture/case1_normalmixture_rejection")

result= data.frame(normal_mixture$R.cai,
                 normal_mixture$R.cai_boot,
                 normal_mixture$R.bh,
                 normal_mixture$R.by,
                 normal_mixture$R.proposed_procedure)
names(result) <- c("R.cai", "R.cai_boot", "R.bh", "R.by", "R.LEE")

boxplot(result,ylim=c(0,N))
title(main="number of rejected Null hypothesis",sub="normal mixture distribution")



normal=read.csv("/Users/baeknarim/학회포스터발표_2022/correlation testing/case1/normal/case1_normal_rejection")

result= data.frame(normal$R.cai,
                   normal$R.cai_boot,
                   normal$R.bh,
                   normal$R.by,
                   normal$R.proposed_procedure)
names(result) <- c("R.cai", "R.cai_boot", "R.bh", "R.by", "R.LEE")

boxplot(result,ylim=c(0,N))
title(main="number of rejected Null hypothesis",sub="normal distribution")



t=read.csv("/Users/baeknarim/학회포스터발표_2022/correlation testing/case1/t/case1_t_rejection")

result= data.frame(t$R.cai,
                   t$R.cai_boot,
                   t$R.bh,
                   t$R.by,
                   t$R.proposed_procedure)
names(result) <- c("R.cai", "R.cai_boot", "R.bh", "R.by", "R.LEE")

boxplot(result,ylim=c(0,N))
title(main="number of rejected Null hypothesis",sub="t distribution")

exp=read.csv("/Users/baeknarim/학회포스터발표_2022/correlation testing/case1/exp/case1_exp_rejection")

result= data.frame(exp$R.cai,
                   exp$R.cai_boot,
                   exp$R.bh,
                   exp$R.by,
                   exp$R.proposed_procedure)
names(result) <- c("R.cai", "R.cai_boot", "R.bh", "R.by", "R.LEE")

boxplot(result,ylim=c(0,N))
title(main="number of rejected Null hypothesis",sub="exp distribution")




#case3
normal_mixture=read.csv("/Users/baeknarim/학회포스터발표_2022/correlation testing/case3/normal_mixture/case3_normalmixture.csv")

FDP = data.frame(normal_mixture$FDP.cai,
                  normal_mixture$FDP.cai_boot,
                  normal_mixture$FDP.bh,
                  normal_mixture$FDP.by,
                  normal_mixture$FDP.proposed_procedure)

names(FDP) <- c("FDP.cai", "FDP.cai_boot", "FDP.bh", "FDP.by", "FDP.LEE")
boxplot(FDP, ylim=c(0,1))
abline(h=0.05,col="red",lty=3)
title(main="FDP",sub="normal_mixture distribution")


normal=read.csv("/Users/baeknarim/학회포스터발표_2022/correlation testing/case3/normal/case3_normal.csv")

FDP = data.frame(normal$FDP.cai,
                 normal$FDP.cai_boot,
                 normal$FDP.bh,
                 normal$FDP.by,
                 normal$FDP.proposed_procedure)

names(FDP) <- c("FDP.cai", "FDP.cai_boot", "FDP.bh", "FDP.by", "FDP.LEE")
boxplot(FDP, ylim=c(0,1))
abline(h=0.05,col="red",lty=3)
title(main="FDP",sub="normal distribution")


t=read.csv("/Users/baeknarim/학회포스터발표_2022/correlation testing/case3/t/case3_t.csv")
FDP = data.frame(t$FDP.cai,
                 t$FDP.cai_boot,
                 t$FDP.bh,
                 t$FDP.by,
                 t$FDP.proposed_procedure)

names(FDP) <- c("FDP.cai", "FDP.cai_boot", "FDP.bh", "FDP.by", "FDP.LEE")
boxplot(FDP, ylim=c(0,1))
abline(h=0.05,col="red",lty=3)
title(main="FDP",sub="t distribution")


exp=read.csv("/Users/baeknarim/학회포스터발표_2022/correlation testing/case3/exp/case3_exp.csv")
FDP = data.frame(exp$FDP.cai,
                 exp$FDP.cai_boot,
                 exp$FDP.bh,
                 exp$FDP.by,
                 exp$FDP.proposed_procedure)

names(FDP) <- c("FDP.cai", "FDP.cai_boot", "FDP.bh", "FDP.by", "FDP.LEE")
boxplot(FDP, ylim=c(0,1))
abline(h=0.05,col="red",lty=3)
title(main="FDP",sub="exp distribution")





normal_mixture=read.csv("/Users/baeknarim/학회포스터발표_2022/correlation testing/case3/normal_mixture/case3_normalmixture.csv")

FNDP = data.frame(normal_mixture$FNDP.cai,
                 normal_mixture$FNDP.cai_boot,
                 normal_mixture$FNDP.bh,
                 normal_mixture$FNDP.by,
                 normal_mixture$FNDP.proposed_procedure)
FNDP=1-FNDP
names(FNDP) <- c("FNDP.cai", "FNDP.cai_boot", "FNDP.bh", "FNDP.by", "FNDP.LEE")
boxplot(FNDP, ylim=c(0,1))

title(main="1-FNDP",sub="normal_mixture distribution")


normal=read.csv("/Users/baeknarim/학회포스터발표_2022/correlation testing/case3/normal/case3_normal.csv")

FNDP = data.frame(normal$FNDP.cai,
                  normal$FNDP.cai_boot,
                  normal$FNDP.bh,
                  normal$FNDP.by,
                  normal$FNDP.proposed_procedure)
FNDP=1-FNDP
names(FNDP) <- c("FNDP.cai", "FNDP.cai_boot", "FNDP.bh", "FNDP.by", "FNDP.LEE")
boxplot(FNDP, ylim=c(0,1))

title(main="1-FNDP",sub="normal distribution")


t=read.csv("/Users/baeknarim/학회포스터발표_2022/correlation testing/case3/t/case3_t.csv")
FNDP = data.frame(t$FNDP.cai,
                  t$FNDP.cai_boot,
                  t$FNDP.bh,
                  t$FNDP.by,
                  t$FNDP.proposed_procedure)
FNDP=1-FNDP
names(FNDP) <- c("FNDP.cai", "FNDP.cai_boot", "FNDP.bh", "FNDP.by", "FNDP.LEE")
boxplot(FNDP, ylim=c(0,1))

title(main="1-FNDP",sub="t distribution")


exp=read.csv("/Users/baeknarim/학회포스터발표_2022/correlation testing/case3/exp/case3_exp.csv")
FNDP = data.frame(exp$FNDP.cai,
                  exp$FNDP.cai_boot,
                  exp$FNDP.bh,
                  exp$FNDP.by,
                  exp$FNDP.proposed_procedure)
FNDP=1-FNDP
names(FNDP) <- c("FNDP.cai", "FNDP.cai_boot", "FNDP.bh", "FNDP.by", "FNDP.LEE")
boxplot(FNDP, ylim=c(0,1))

title(main="1-FNDP",sub="exp distribution")


#FDP : type 1 error(R/V) 
#FNDP : 기각 해야되는데, 기각 안된 애들의 비율  : type 2 error Power개념.. (reject안한것중에서 H1의 비율 /N-R) -얼마나 작은지! 

#새 방법론의 boxplot상자가 적던지, 0.05를 좀 더 잘 맞추던지. FNDP가 작던지... 그래야함 

#cai 일반적인 함수 찾아보기 

#1,2개 나올떄까지 level 높여서, 1개가 나오는 level이 어느정도고 그 1개는 뭐다~(아마 0.3정도)
#FAN같은거 돌렸을때 ranking을 보여줘서, 제일 유의한게 어떤애들이 있는지, ROI pair (TOP-5) 보여주기 
#알츠하이머면 sample size, p, 정도만 쓰고 data script는 금방
#brain 시각화, FNP,FDNP boxplot
#각 방법론들 가볍게 소개하기 
#PCA 넣을꺼면 알고리즘 테이블 넣어주기 
#simulation setting 어떻게 generate하는지 보여주기 




#### Libraries ####
library(ggplot2)
library(reshape2)

#### Data ####
flaws<-read.table("subwaydat.txt", header=TRUE)
flaws$rate<-flaws$flaws/flaws$length

############# Estimation #############

#### Homogeneous Model ####
## Likelihood ##
simple_lik<-function(lambda, y, n){
  sum(log((lambda*n)^y*exp(-lambda*n)))
}

## MLEs
lambda.est<-optim(.5, function(u) -simple_lik(u, y=flaws$flaws, n=flaws$length), hessian=TRUE)
ci<-c(lambda.est$par - qnorm(0.975)*1/sqrt(lambda.est$hessian), lambda.est$par + qnorm(0.975)*1/sqrt(lambda.est$hessian))

## Get Predictions 
flaws$exp.counts <- lambda.est$par*flaws$length
chi.test<-sum((flaws$flaws - flaws$exp.counts)^2/flaws$exp.counts)
chi.test

qchisq(0.975, 74)

ddply(data=flaws, .(length), summarize, flaws=sum(flaws), exp.counts=sum(exp.counts))


### Profile phi ###
## I think I need constraints on lambda1, lambda2, i.e. lambda1 > lambda2 ###
likelihood_pois<-function(par, phi, y, n){
  l1<-par[1]
  l2<-par[2]
  ifelse(l1 > l2 & l1 >0 & l2>0, sum(log(phi*(l1*n)^(y)*exp(-l1*n) + (1-phi)*(l2*n)^(y)*exp(-l2*n))), -1e50)
}

Lp_pois<-function(phi, y, n){
  mles <- optim(par0s, function(a){-likelihood_pois(a, phi=phi, y=y, n=n)})$par
  lp<-likelihood_pois(mles, y=y, phi=phi, n=n)
  par0s<<-mles
  return(c(l1=mles[1], l2=mles[2], phi=phi, lp=lp))
}

Lp_pois2<-function(phi, y, n){
  mles <- optim(par0s, function(a){-likelihood_pois(a, phi=phi, y=y, n=n)})$par
  lp<-likelihood_pois(mles, y=y, phi=phi, n=n)
  par0s<<-mles
  return(lp)
}

phi.est<-optimize(function(u) Lp_pois2(u, y=flaws$flaws, n=flaws$length), interval=c(0,1), maximum=TRUE)$max

optim(c(2,1), function(a) -likelihood_pois(a, phi=phi.est, y=flaws$flaws, n=flaws$length))
## phi = 0.4967

par0s<-c(5,1)
phi<-matrix(seq(0,1,.0001),ncol=1)
df<-as.data.frame(t(apply(phi, 1, function(u) Lp_pois(u, y=flaws$flaws, n=flaws$length))))  
df[which.max(df$lp),]
#       l1        l2    phi       lp
# 3.158086 0.5266318 0.4967 582.3936


## No profile ##
likelihood_pois2<-function(par, y, n){
  l1<-par[1]
  l2<-par[2]
  phi<-par[3]
  ifelse(l1 > l2 & l1 >0 & l2>0, sum(log(phi*(l1*n)^(y)*exp(-l1*n) + (1-phi)*(l2*n)^(y)*exp(-l2*n))), -1e50)
}

optim(c(2,1,.4), function(u) -likelihood_pois2(u, y=flaws$flaws, n=flaws$length))

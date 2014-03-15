###Zero-Inflated Poisson Model 
 source("code/library.R")
 source("code/data_format.R")
 source("code/helpers.R")
 source("code/newtraph_kaiser.R")

#### Function for newtraph ####
ders.zip <- function(ps, dat) {
  ##likelihood
  b0 <- ps[1]
  b1 <- ps[2]
  b2 <- ps[3]
  b3 <- ps[4]  
  
  x <- dat$price
  y <- dat$mvm
  
  T1 <- exp(b0 + b1*x)
  T2 <- exp(b2 + b3*x)
  
  p <- T1/(1+T1)
  lambda <- T2  
  
  delta <- as.numeric(y==0)
  
  loglik <- loglik.zip(ps, dat)
  
  ##gradient
  dldp <- delta*(exp(-lambda)-1)/(1-p + p*exp(-lambda)) + (1-delta)/p 
  dldlam <- delta*(-p*exp(-lambda))/(1-p + p*exp(-lambda)) + (1-delta)*(y/lambda-1)
  dpdb0 <- T1/(1+T1)^2
  dpdb1 <- x*T1/(1+T1)^2
  dpdb2 <- 0
  dpdb3 <- 0
  dlamdb0 <- 0
  dlamdb1 <- 0
  dlamdb2 <- T2
  dlamdb3 <- x*T2
   
  grad <- c(sum(dldp*dpdb0), sum(dldp*dpdb1), sum(dldlam*dlamdb2), sum(dldlam*dlamdb3))
   
  ##hessian
  d2ldp2 <- -delta*((exp(-lambda)-1)/(1-p + p*exp(-lambda)))^2 - (1-delta)/p^2 
  d2ldpdlam <- -delta*(exp(-lambda))/((1-p) + p*exp(-lambda))^2
  d2ldlam2 <- delta*p*(1-p)*exp(-lambda)/(1-p + p*exp(-lambda))^2 - (1-delta)*(y/lambda^2)
  d2pdb02 <- (1-T1)*T1/(1+T1)^3
  d2pdb12 <- x^2*(1-T1)*T1/(1+T1)^3
  d2pdb22 <- 0
  d2pdb32 <- 0  
  d2pdb0db1 <- x*(1-T1)*T1/(1+T1)^3
  d2pdb0db2 <- 0
  d2pdb0db3 <- 0
  d2pdb1db2 <- 0
  d2pdb1db3 <- 0
  d2pdb2db3 <- 0

  d2lamdb02 <- 0
  d2lamdb12 <- 0
  d2lamdb22 <- T2
  d2lamdb32 <- x^2*T2
  d2lamdb0db1 <- 0
  d2lamdb0db2 <- 0
  d2lamdb0db3 <- 0
  d2lamdb1db2 <- 0
  d2lamdb1db3 <- 0
  d2lamdb2db3 <- x*T2
  
  
  hess00 <- sum((d2ldp2*dpdb0 + d2ldpdlam*dlamdb0)*dpdb0 + d2pdb02*dldp + (d2ldpdlam*dpdb0 + d2ldlam2*dlamdb0)*dlamdb0 + d2lamdb02*dldlam) 
  hess01 <- sum((d2ldp2*dpdb0 + d2ldpdlam*dlamdb0)*dpdb1 + d2pdb0db1*dldp + (d2ldpdlam*dpdb0 + d2ldlam2*dlamdb0)*dlamdb1 + d2lamdb0db1*dldlam) 
  hess02 <- sum((d2ldp2*dpdb0 + d2ldpdlam*dlamdb0)*dpdb2 + d2pdb0db2*dldp + (d2ldpdlam*dpdb0 + d2ldlam2*dlamdb0)*dlamdb2 + d2lamdb0db2*dldlam) 
  hess03 <- sum((d2ldp2*dpdb0 + d2ldpdlam*dlamdb0)*dpdb3 + d2pdb0db3*dldp + (d2ldpdlam*dpdb0 + d2ldlam2*dlamdb0)*dlamdb3 + d2lamdb0db3*dldlam) 

  hess11 <- sum((d2ldp2*dpdb1 + d2ldpdlam*dlamdb1)*dpdb1 + d2pdb12*dldp + (d2ldpdlam*dpdb1 + d2ldlam2*dlamdb1)*dlamdb1 + d2lamdb12*dldlam) 
  hess12 <- sum((d2ldp2*dpdb1 + d2ldpdlam*dlamdb1)*dpdb2 + d2pdb1db2*dldp + (d2ldpdlam*dpdb1 + d2ldlam2*dlamdb1)*dlamdb2 + d2lamdb1db2*dldlam) 
  hess13 <- sum((d2ldp2*dpdb1 + d2ldpdlam*dlamdb1)*dpdb3 + d2pdb1db3*dldp + (d2ldpdlam*dpdb1 + d2ldlam2*dlamdb1)*dlamdb3 + d2lamdb1db3*dldlam) 

  hess22 <- sum((d2ldp2*dpdb2 + d2ldpdlam*dlamdb2)*dpdb2 + d2pdb22*dldp + (d2ldpdlam*dpdb2 + d2ldlam2*dlamdb2)*dlamdb2 + d2lamdb22*dldlam) 
  hess23 <- sum((d2ldp2*dpdb2 + d2ldpdlam*dlamdb2)*dpdb3 + d2pdb2db3*dldp + (d2ldpdlam*dpdb2 + d2ldlam2*dlamdb2)*dlamdb3 + d2lamdb2db3*dldlam) 

  hess33 <- sum((d2ldp2*dpdb3 + d2ldpdlam*dlamdb3)*dpdb3 + d2pdb32*dldp + (d2ldpdlam*dpdb3 + d2ldlam2*dlamdb3)*dlamdb3 + d2lamdb32*dldlam) 
  
  
  hess <- matrix(c(hess00, hess01, hess02, hess03, 
                   hess01, hess11, hess12, hess13,
                   hess02, hess12, hess22, hess23,
                   hess03, hess13, hess23, hess33), ncol = 4)
  
  return(list(loglik, grad, hess))  
}

## likelihood function
loglik.zip <- function(par, dat){
  b0<-par[1]
  b1<-par[2]
  b2<-par[3]
  b3<-par[4]
  
  x <- dat$price
  y <- dat$mvm
  
  T1 <- exp(b0 + b1*x)
  T2 <- exp(b2 + b3*x)  
  
  p <- T1/(1+T1)
  lambda <- T2
  
  delta<-rep(0,length(y))
  delta[y==0]<-1
  
  sum(delta*log(1-p + p*exp(-lambda)) + (1-delta)*(log(p) - lfactorial(y) + y*log(lambda) - lambda))
}

## get estimates
zpois.est <- dlply(dat, .(store), function(x) {
  x0 <- optim(c(0,0,0,0), function(u) -loglik.zip(u, dat=x))$par
  tryCatch(newtraph(ders=ders.zip, dat=x, x0=x0), 
           error=function(e) list(e))
})

zpois.mles <- ldply(zpois.est, function(x) if(length(x) == 3) t(x[[1]]))
zpois.mles.m <- melt(zpois.mles, id.vars="store")
dist.coef.zpois.plot <- qplot(data=zpois.mles.m, x=value, xlab="Coefficients", ylab="Count", fill=I("grey60"), colour=I("black")) + facet_wrap(~variable, scales="free_x")

### create plots for viewing fit
# plot_by_store <- function(store_id, data=dat) {
#   dat1 <- subset(data, store == store_id)
#   
#   est1 <- est[[as.character(dat1$store[1])]][[1]]
#   dat1$p <- exp(est1[1] + est1[2]*dat1$price)/(1 + exp(est1[1] + est1[2]*dat1$price))
#   dat1$lambda <- exp(est1[3] + est1[4]*dat1$price)
#   
#   distns <- ldply(mlply(unique(dat1$price), function(y) mdply(0:max(subset(dat1, price == y)$mvm), function(x) c(mvm=x, price=y, lik=exp(loglik.zip(est1, dat=list(price=y, mvm=x)))))), function(x) x)
#   
#   gg <- ggplot(dat1) + geom_histogram(aes(x=mvm, y=..density..)) + 
#     geom_line(aes(x=mvm, y=lik), data=distns, colour=I("red")) + facet_wrap(~price, scales="free")  +
#     ggtitle(sprintf("Store ID: %d", store_id))
#   
#   ##ggsave(sprintf("EA 1/figure/fit_store_%d.pdf", store_id), gg)
# }
# 
# mlply(1:length(est), function(x) {
#   if(length(est[[x]]) == 3) {
#     cat(sprintf("Counter: %d \r", x))
#     plot_by_store(as.numeric(names(est)[x]))     
#   }
# })

#### Generate Data ####
generate_data_zip<-function(store_dat, b, nsets=100){
 ys<-NULL
 x<-store_dat$price
 p <- exp(b[1] + b[2]*x)/(1+exp(b[1] + b[2]*x))
 lambda<-exp(b[3]+b[4]*x)
 for(i in 1:length(p)){
   z<-rbinom(nsets,1,p[i])
   y<-rpois(nsets, lambda[i])
   y[z==0]<-0
   ys<-rbind(ys, y)
   i<-i+1
 }
 return(ys)
} 
 
gens_zip<-list(NA) 
z<-unique(dat$store)
for(i in 1:length(z)){
  gens_zip[[i]]<-generate_data_zip(dat[dat$store==z[i],], as.numeric(zpois.mles[1,-1]), nsets=1000)  
cat(i, "\r")
}

save(gens_zip, "gens_zip.rda") 

## Proportion of 0's ##
prop0<-function(dat, gen_dat){
 orig<-sum(dat$mvm==0)/length(dat$mvm)
 gens<-apply(gen_dat, 2, function(u) sum(u==0)/length(u))
 pval=sum(gens > orig)/length(gens)
 return(pval)
} 

## Number of days with more than 40 sales ##
more40<-function(dat, gen_dat){
 orig<-sum(dat$mvm > 40)
 gens<-apply(gen_dat, 2, function(u) sum(u > 40))
 pval=sum(gens > orig)/length(gens)
 return(pval)
}

## Sample Variance ##
vars<-function(dat, gen_dat){
 orig<-var(dat$mvm)
 gens<-apply(gen_dat, 2, var)
 pval=sum(gens > orig)/length(gens)
 return(pval)
}

## Sample Variance ##
less4<-function(dat, gen_dat){
 orig<-sum(dat$mvm<=4)
 gens<-apply(gen_dat, 2, function(u) sum(u <= 4))
 pval=sum(gens > orig)/length(gens)
 return(pval)
}

 
### p-values for 4 criteria ###
pvals_prop0<-c()
 for(i in 1:length(z)){
   pvals_prop0<-c(pvals_prop0, prop0(dat[dat$store==z[i],], gens_zip[[i]]))
}

sum(pvals_prop0 > 0.05 & pvals_prop0 < 0.95)

pvals_more40<-c()
 for(i in 1:length(z)){
   pvals_more40<-c(pvals_more40, more40(dat[dat$store==z[i],], gens_zip[[i]]))
}

sum(pvals_more40 > 0.05 & pvals_more40 < 0.95)

pvals_var<-c()
 for(i in 1:length(z)){
   pvals_var<-c(pvals_var, vars(dat[dat$store==z[i],], gens_zip[[i]]))
}
 
sum(pvals_var > 0.05 & pvals_var < 0.95)

pvals_less4<-c()
 for(i in 1:length(z)){
   pvals_less4<-c(pvals_less4, less4(dat[dat$store==z[i],], gens_zip[[i]]))
} 
 
sum(pvals_less4 > 0.05 & pvals_less4 < 0.95)
 
crits<-cbind(pvals_prop0, pvals_more40, pvals_less4, pvals_var)

## 52 out of 161 stores can pass 1 criterion
sum(apply(crits, 1, function(u) any(u > 0.05 & u < 0.95)))

## 0 stores can pass all
sum(apply(crits, 1, function(u) all(u > 0.05 & u < 0.95)))
 

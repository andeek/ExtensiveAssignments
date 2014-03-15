#### Zero Inflated Poisson Gamma ####
# source("code/library.R")
# source("code/data_format.R")
# source("code/helpers.R")
# source("code/newtraph_kaiser.R")

## Function for newtraph
ders.zipg <- function(ps, dat, phi) {
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
  mu <- T2  
  
  T3 <- (1/(1+phi))^(mu/phi)
  T4 <- log(1/(1+phi))/phi
  
  delta <- as.numeric(y==0)
  
  loglik <- loglik(ps, dat, phi)
  
  ##gradient
  dldp <- delta*(T3 - 1)/(1-p + p*T3) + (1-delta)/p 
  dldmu <- delta*(p*T3*T4)/(1-p + p*T3) + (1-delta)*(-digamma(mu/phi) - log(phi) + digamma(y + mu/phi) + phi*T4)/phi
  dpdb0 <- T1/(1+T1)^2
  dpdb1 <- x*T1/(1+T1)^2
  dpdb2 <- 0
  dpdb3 <- 0
  dmudb0 <- 0
  dmudb1 <- 0
  dmudb2 <- T2
  dmudb3 <- x*T2
  
  grad <- c(sum(dldp*dpdb0), sum(dldp*dpdb1), sum(dldmu*dmudb2), sum(dldmu*dmudb3))
  
  ##hessian
  d2ldp2 <- -delta*((T3-1)/(1-p + p*T3))^2 - (1-delta)/p^2 
  d2ldpdmu <- -delta*(T3*T4)/((1-p) + p*T3)^2
  d2ldmu2 <- delta*p*(1-p)*T3*T4^2/(1-p + p*T3)^2 + (1-delta)/phi^2*(-trigamma(mu/phi) + trigamma(y + mu/phi))
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
  
  d2mudb02 <- 0
  d2mudb12 <- 0
  d2mudb22 <- T2
  d2mudb32 <- x^2*T2
  d2mudb0db1 <- 0
  d2mudb0db2 <- 0
  d2mudb0db3 <- 0
  d2mudb1db2 <- 0
  d2mudb1db3 <- 0
  d2mudb2db3 <- x*T2
  
  
  hess00 <- sum((d2ldp2*dpdb0 + d2ldpdmu*dmudb0)*dpdb0 + d2pdb02*dldp + (d2ldpdmu*dpdb0 + d2ldmu2*dmudb0)*dmudb0 + d2mudb02*dldmu) 
  hess01 <- sum((d2ldp2*dpdb0 + d2ldpdmu*dmudb0)*dpdb1 + d2pdb0db1*dldp + (d2ldpdmu*dpdb0 + d2ldmu2*dmudb0)*dmudb1 + d2mudb0db1*dldmu) 
  hess02 <- sum((d2ldp2*dpdb0 + d2ldpdmu*dmudb0)*dpdb2 + d2pdb0db2*dldp + (d2ldpdmu*dpdb0 + d2ldmu2*dmudb0)*dmudb2 + d2mudb0db2*dldmu) 
  hess03 <- sum((d2ldp2*dpdb0 + d2ldpdmu*dmudb0)*dpdb3 + d2pdb0db3*dldp + (d2ldpdmu*dpdb0 + d2ldmu2*dmudb0)*dmudb3 + d2mudb0db3*dldmu) 
  
  hess11 <- sum((d2ldp2*dpdb1 + d2ldpdmu*dmudb1)*dpdb1 + d2pdb12*dldp + (d2ldpdmu*dpdb1 + d2ldmu2*dmudb1)*dmudb1 + d2mudb12*dldmu) 
  hess12 <- sum((d2ldp2*dpdb1 + d2ldpdmu*dmudb1)*dpdb2 + d2pdb1db2*dldp + (d2ldpdmu*dpdb1 + d2ldmu2*dmudb1)*dmudb2 + d2mudb1db2*dldmu) 
  hess13 <- sum((d2ldp2*dpdb1 + d2ldpdmu*dmudb1)*dpdb3 + d2pdb1db3*dldp + (d2ldpdmu*dpdb1 + d2ldmu2*dmudb1)*dmudb3 + d2mudb1db3*dldmu) 
  
  hess22 <- sum((d2ldp2*dpdb2 + d2ldpdmu*dmudb2)*dpdb2 + d2pdb22*dldp + (d2ldpdmu*dpdb2 + d2ldmu2*dmudb2)*dmudb2 + d2mudb22*dldmu) 
  hess23 <- sum((d2ldp2*dpdb2 + d2ldpdmu*dmudb2)*dpdb3 + d2pdb2db3*dldp + (d2ldpdmu*dpdb2 + d2ldmu2*dmudb2)*dmudb3 + d2mudb2db3*dldmu) 
  
  hess33 <- sum((d2ldp2*dpdb3 + d2ldpdmu*dmudb3)*dpdb3 + d2pdb32*dldp + (d2ldpdmu*dpdb3 + d2ldmu2*dmudb3)*dmudb3 + d2mudb32*dldmu) 
  
  
  hess <- matrix(c(hess00, hess01, hess02, hess03, 
                   hess01, hess11, hess12, hess13,
                   hess02, hess12, hess22, hess23,
                   hess03, hess13, hess23, hess33), ncol = 4)
  
  return(list(loglik, grad, hess))  
}

## likelihood function
loglik <- function(par, dat, phi){
  b0<-par[1]
  b1<-par[2]
  b2<-par[3]
  b3<-par[4]
  
  x <- dat$price
  y <- dat$mvm
  
  T1 <- exp(b0 + b1*x)
  T2 <- exp(b2 + b3*x)  
  
  p <- T1/(1+T1)
  mu <- T2
  
  delta<-rep(0,length(y))
  delta[y==0]<-1
  
  sum(delta*log(1-p + p*(1/(1+phi))^(mu/phi)) + (1-delta)*(log(p) - lfactorial(y) - lgamma(mu/phi) - mu/phi*log(phi) + lgamma(y + mu/phi) + (y+mu/phi)*log(phi/(1+phi))))
}

## profile out phi
phis <- c(seq(0.001, 0.019, by=.001), seq(1.01, 1.05, by = .005))
profile.phi.results <- ldply(phis, function(phi) {
  cat(paste("phi:", phi))
  zipoisgamma_NR <- optim(c(0,0,0,0), function(u) -loglik(u, dat=dat, phi=phi)) 
  c(phi = phi, loglik = zipoisgamma_NR$value, b0 = zipoisgamma_NR$par[1], b1 = zipoisgamma_NR$par[2], b2 = zipoisgamma_NR$par[3], b3 = zipoisgamma_NR$par[4], converge=zipoisgamma_NR$convergence)
})


store1 <- samp_store(1, 1, dat)
store1_opt <- optim(c(0,0,0,0), function(u) -loglik(u, dat=store1, phi=1), hessian=TRUE)
ders.zipg(store1_opt$par, store1, 1)
newtraph(ders=ders.zipg, dat=store1, x0=, phi=1)

est <- dlply(dat, .(store), function(x) {
  x0 <- optim(c(0,0,0,0), function(u) -loglik(u, dat=x))$par
  tryCatch(newtraph(ders=ders.zipg, dat=x, x0=x0, phi=1), 
           error=function(e) list(e))
})



zigp.newtraph <- profile.alpha.results[which.max(profile.alpha.results$loglik), ]
par.zigp.newtraph <- as.numeric(zigp.newtraph[1,3:6])


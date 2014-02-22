
ders.newt <- function(ps, dat, theta) {
  ##likelihood
  b1 <- ps[1]
  b2 <- ps[2]
  b3 <- ps[3]
  sigma2 <- ps[4]  
  
  #   sigma2 <- gom.gls$sshat
  
  x <- dat$x
  y <- dat$y
  
  T1 <- exp(-exp(b2 - b3*x))
  T2 <- exp(b2 - b3*x)  
  u <- b1*T1
  v <- sigma2*u^(2*theta)
  
  loglik <- sum(-0.5*log(2*pi) - 0.5*log(v) - (y - u)^2/(2*v))
  
  
  ##gradient
  dldu <- (y - u)/v
  dldv <- -1/2*(1/v - ((y - u)/v)^2)
  dvdu <- sigma2*2*theta*u^(2*theta - 1)
  dudb1 <- T1
  dudb2 <- -b1*T1*T2
  dudb3 <- b1*x*T1*T2
  
  dvdsigma2 <- u^(2*theta)
  dldsigma2 <- dldv*dvdsigma2
  
  grad <- c(sum((dldu + dldv*dvdu)*dudb1), sum((dldu + dldv*dvdu)*dudb2), sum((dldu + dldv*dvdu)*dudb3), sum(dldsigma2))
  #   grad <- c(sum((dldu + dldv*dvdu)*dudb1), sum((dldu + dldv*dvdu)*dudb2), sum((dldu + dldv*dvdu)*dudb3))
  
  ##hessian
  d2ldu2 <- -1/v
  d2ldudv <- -(y - u)/(v^2)
  d2ldv2 <- 1/(2*v^2) - (y - u)^2/(v^3)
  d2vdu2 <- 2*theta*(2*theta - 1)*sigma2*u^(2*theta - 2)
  d2udb12 <- 0
  d2udb1db2 <- -T1*T2
  d2udb1db3 <- x*T1*T2
  d2udb22 <- b1*(T1*T2^2 - T1*T2)
  d2udb2db3 <- b1*x*(-T1*T2^2 + T1*T2)
  d2udb32 <- b1*x^2*(-T1*T2^2 + T1*T2)
  
  d2vdsigma2du <- 2*theta*u^(2*theta - 1)
  #d2ldsigma2dv <- u/(2*v^2) - u*(y-u)^2/(v^3)
  #d2ldsigma2du <- -.5*(1/v - (y-u)^2/(v^2)) - u*((y-u)/(v^2))
  
  hess11 <- sum((d2ldu2 + d2ldudv*dvdu)*dudb1*dudb1 + dldu*d2udb12 + (d2ldv2*dvdu + d2ldudv)*dudb1*dvdu*dudb1 + dldv*dvdu*d2udb12 + dldv*d2vdu2*dudb1*dudb1)
  hess12 <- sum((d2ldu2 + d2ldudv*dvdu)*dudb2*dudb1 + dldu*d2udb1db2 + (d2ldv2*dvdu + d2ldudv)*dudb2*dvdu*dudb1 + dldv*dvdu*d2udb1db2 + dldv*d2vdu2*dudb1*dudb2)
  hess13 <- sum((d2ldu2 + d2ldudv*dvdu)*dudb3*dudb1 + dldu*d2udb1db3 + (d2ldv2*dvdu + d2ldudv)*dudb3*dvdu*dudb1 + dldv*dvdu*d2udb1db3 + dldv*d2vdu2*dudb1*dudb3)
  
  hess22 <- sum((d2ldu2 + d2ldudv*dvdu)*dudb2*dudb2 + dldu*d2udb22 + (d2ldv2*dvdu + d2ldudv)*dudb2*dvdu*dudb2 + dldv*dvdu*d2udb22 + dldv*d2vdu2*dudb2*dudb2)
  hess23 <- sum((d2ldu2 + d2ldudv*dvdu)*dudb3*dudb2 + dldu*d2udb2db3 + (d2ldv2*dvdu + d2ldudv)*dudb3*dvdu*dudb2 + dldv*dvdu*d2udb2db3 + dldv*d2vdu2*dudb2*dudb3)
  
  hess33 <- sum((d2ldu2 + d2ldudv*dvdu)*dudb3*dudb3 + dldu*d2udb32 + (d2ldv2*dvdu + d2ldudv)*dudb3*dvdu*dudb3 + dldv*dvdu*d2udb32 + dldv*d2vdu2*dudb3*dudb3)
  
  hess14 <- sum((d2ldv2*dvdu + d2ldudv)*dudb1*dvdsigma2 + dldv*dudb1*d2vdsigma2du)
  hess24 <- sum((d2ldv2*dvdu + d2ldudv)*dudb2*dvdsigma2 + dldv*dudb2*d2vdsigma2du)
  hess34 <- sum((d2ldv2*dvdu + d2ldudv)*dudb3*dvdsigma2 + dldv*dudb3*d2vdsigma2du)
  hess44 <- sum(d2ldv2*(dvdsigma2)^2)
  
  #   hess<- matrix(c(hess11, hess12, hess13, 
  #                  hess12, hess22, hess23,
  #                  hess13, hess23, hess33), ncol = 3)
  
  hess <- matrix(c(hess11, hess12, hess13, hess14, 
                   hess12, hess22, hess23, hess24,
                   hess13, hess23, hess33, hess34,
                   hess14, hess24, hess34, hess44), ncol = 4)
  
  return(list(loglik, grad, hess))  
}
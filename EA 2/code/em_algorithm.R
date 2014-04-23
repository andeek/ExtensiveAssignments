phat_fun<-function(phi, l1, l2, y, n){
  (phi*(n*l1)^y*exp(-l1*n)/factorial(y))/(phi*(n*l1)^y*exp(-l1*n)/factorial(y) + (1-phi)*(n*l2)^y*exp(-l2*n)/factorial(y))
}

hessian<-function(phi, l1, l2, y, n){
  f<-function(l){
    exp(-l*n)*(n*l)^y/factorial(y)
  }
  h<-phi*f(l1) + (1-phi)*f(l2)

  dfdlambda<-function(l){
    f(l)*(y/l-n)
  }
  df2d2lambda<-function(l){
    f(l)*(y/l-n)^2 - f(l)*(y/l^2)
  }

  
  dhdl1<-phi*dfdlambda(l1)
  dhdl2<-(1-phi)*dfdlambda(l2)
  dhdphi<-f(l1) - f(l2)
  d2hd2l1<-phi*df2d2lambda(l1)
  d2hd2l2<-(1-phi)*df2d2lambda(l2)
  d2hd2phi<-0
  d2hdl1l2<-0
  d2dl1dphi<-dfdlambda(l1)
  d2dl2dphi<- -dfdlambda(l2)
  
  hess11<-sum(-dhdphi^2/h^2)
  hess22<-sum(d2hd2l1/h - dhdl1^2/h^2)
  hess33<-sum(d2hd2l2/h - dhdl2^2/h^2)
  hess12<-sum(d2dl1dphi/h - dhdl1*dhdphi/h^2)
  hess13<-sum(d2dl2dphi/h - dhdl2*dhdphi/h^2)
  hess23<-sum(d2hdl1l2/h - dhdl1*dhdl2/h^2)
  hessian<-matrix(c(hess11, hess12, hess13, hess12, hess22, hess23, hess13, hess23, hess33), ncol=3)
  return(-hessian)
}

em_algorithm<-function(phi, l1, l2, y, n, maxit=1000, conv.crit=10e-8){
  j<-1
  repeat{
    phat<-phat_fun(phi, l1, l2, y, n)
    phi2<-mean(phat)
    l1_2<-sum(y*phat)/sum(n*phat)
    l2_2<-sum(y*(1-phat))/sum(n*(1-phat))
    if(abs(phi2-phi) < conv.crit & abs(l1_2-l1) < conv.crit & abs(l2_2-l2) < conv.crit){
      cat("Convergence reached\n")
      break 
    } 
    if(j == maxit){
      cat("Maximum iterations reached\n")
      break
    }else{
      j<-j+1
      phi<-phi2
      l1<-l1_2
      l2<-l2_2
    }
  }
  return(c(phi=phi2, lambda1=l1, lambda2=l2, n.iter=j))
}

ests<-em_algorithm(.5, 5, 0.01, y=subway.dat$flaws, n=subway.dat$length)
## Match mles from optim

## So, how do we get probability each observation is in group A?
probs<-round(phat_fun(ests[1], ests[2], ests[3], y=subway.dat$flaws, n=subway.dat$length),3)
qplot(x=subway.dat$flaws, y=probs)
subway.dat$probs<-probs

hessian(ests[1], ests[2], ests[3], y=subway.dat$flaws, n=subway.dat$length)
# I think there must be a mistake somewhere, one positive two negative. Too tired, will check tomorrow.


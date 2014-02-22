newtraph <- function(ders, dat, x0) {
  ##function maximizes, for minization multiply likelihood by -1

  ##initialize critical values, initial values, etc.
  crit <- c(1e-10, 1e-06, 1e-06, 1000)
  stopped <- rep(0,3)
  curnt <- x0
  p<-length(x0)
  k <- 1
  
  while (sum(stopped) < 3 | k < crit[4]) {
    k <- k + 1
    int <- ders(curnt, dat)
    logL<-int[[1]]
    gi<-int[[2]]
    Gi <- int[[3]]
    
    step <- solve(Gi) %*% gi
    new <- curnt - step
    sc <- 1
    
    while(sc < 10) {
      sc <- sc + 1
      check <- ders(new, dat)
      if(check[[1]] < logL) {new <- curnt - (1/sc) * step}
      if(check[[1]] >= logL) {
        newL<-check[[1]]
        newg<-check[[2]]
        break
      }                        
    }
    if(sc == 10) stop("Step halving not effective, try new starting values.")
          
    dist1<-abs(newL-logL)
    dist2 <- sqrt(sum((new - curnt)^2))
    
    if(crit[1] > dist1){stopped[1] <- stopped[1] + 1} #convergence of loglik
    if(crit[2] > dist2){stopped[2] <- stopped[2] + 1} #convergence of estimates
    if(sum(crit[3] > abs(newg)) == p){ stopped[3] <- stopped[3] + 1} # convergence for sum of derivs
    curnt <- new
  }
  if(k == crit[4]) stop(sprintf("Convergence not reached in %x iterations.", crit[4]))
  
  final <- ders(new, dat)
  res <- list(pars = new, loglik = final[[1]], invNegHess = -1*solve(final[[3]]), grad = final[[2]])
  return(res)
}

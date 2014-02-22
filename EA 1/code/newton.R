newton <- function(ftn, x0, tol = 1e-9, max.iter = 100){
  x0 <- x0
  fx <- ftn(x0)
  iter <- 0
  while((abs(fx[1]) > tol) && (iter < max.iter)){
    x1 <- x0 - fx[1]/fx[2]
    fx <- ftn(x1)
    x0 <- x1
    iter <- iter+1
  }
  if(abs(fx[1]) > tol){
    cat("Algorithm failed to converge\n")
  }else {
    return(x0)
  }
}

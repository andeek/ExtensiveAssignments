cbindPad <- function(...){
  args <- list(...)
  n <- sapply(args,nrow)
  mx <- max(n)
  pad <- function(x, mx){
    if (nrow(x) < mx){
      nms <- colnames(x)
      padTemp <- matrix(NA,mx - nrow(x), ncol(x))
      colnames(padTemp) <- nms
      return(rbind(x,padTemp))
    }
    else{
      return(x)
    }
  }
  rs <- lapply(args,pad,mx)
  return(do.call(cbind,rs))
}
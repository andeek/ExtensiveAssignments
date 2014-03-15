### Sample Stores ###
samp_store<-function(n, seed, dat){
  set.seed(seed)
  s<-sample(unique(dat$store),n)
  dat[dat$store %in% s,]
}

### Root MSE ###
rmse<-function(x){
  sqrt(sum((x$y - fitted(x))^2/length(x$y)))
}








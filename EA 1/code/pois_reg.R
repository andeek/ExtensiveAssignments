#### Libraries ####
library(lubridate)
library(reshape2)
library(plyr)
library(ggplot2)

#### Data ####
dat<-read.table(file='./GitHub/ExtensiveAssignments/EA 1/data/greenbeandat.txt', header=TRUE)[,c(1,3,4,5)]
dat$date <-  ymd(as.character(dat$date))

#### Sample Stores ####
samp_store<-function(n, seed, dat){
  set.seed(seed)
  s<-sample(unique(dat$store),n)
  dat[dat$store %in% s,]
}

#### Poisson Regression ####
pois_fits <- dlply(dat, .(store), function(x) {
  glm(data=x, mvm ~ price, family=poisson, maxit=50)
})

rmse<-function(x){
  sqrt(sum((x$y - fitted(x))^2/length(x$y)))
}

rmse_store<-round(as.numeric(sapply(pois_fits, rmse)),3)













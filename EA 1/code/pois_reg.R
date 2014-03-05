#### Libraries ####
library(lubridate)
library(reshape)
library(plyr)
library(ggplot2)
theme_set(theme_bw(base_family="serif"))

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
## Stores 7022 and 7055 have constant price! Discard ##
pois_fits <- dlply(dat[!dat$store %in% c(7022, 7055),], .(store), function(x) {
  glm(data=x, mvm ~ price, family=poisson, maxit=50)
})

### Root MSE ###
rmse<-function(x){
  sqrt(sum((x$y - fitted(x))^2/length(x$y)))
}

rmse_store<-round(as.numeric(sapply(pois_fits, rmse)),3)

### Coefficients ###
ests<-t(sapply(pois_fits, function(u) u$coef))

df<-melt(ests)
qplot(data=df, x=value) + facet_wrap(~X2, scales="free_x")

### Which stores Price not significant ###
signif<-t(sapply(pois_fits, function(u) as.numeric(summary(u)$coef[,4])))

sum(apply(signif, 1, function(u) u[2] > 0.5))







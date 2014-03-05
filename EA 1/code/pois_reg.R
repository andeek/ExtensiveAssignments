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
## Stores 7030 and 7035 are the bad stores ##
which.min(ests[,1])
which.min(ests[,2])

## Why are 7030 and 7035 bad? very few sales ##
ndays_not0<-ddply(dat, .(store), function(x) sum(x$mvm!=0))
ndays_not0[order(ndays_not0[,2]),]
# apparently 7022 had no sales
# conclusion: remove 7022, 7030, 7035, 7055, at least

pois_fits <- dlply(dat[!dat$store %in% c(7022, 7030, 7035, 7055),], .(store), function(x) {
  glm(data=x, mvm ~ price, family=poisson, maxit=50)
})

rmse_store<-round(as.numeric(sapply(pois_fits, rmse)),3)

qplot(x=rmse_store, xlab="RMSE", ylab="Count", fill=I("grey60"), colour=I("black"), binwidth=2)

ests<-t(sapply(pois_fits, function(u) u$coef))
df<-melt(ests)
qplot(data=df, x=value, xlab="Coefficients", ylab="Count", fill=I("grey60"), colour=I("black")) + facet_wrap(~X2, scales="free_x")
## get rid of 7025, 7042, 1161? also the weird ones

ests<-t(sapply(pois_fits[-c(53, 186, 187)], function(u) u$coef))
df<-melt(ests)
qplot(data=df, x=value, xlab="Coefficients", ylab="Count", fill=I("grey60"), colour=I("black"), binwidth=.5) + facet_wrap(~X2, scales="free_x")

### Which stores Price not significant ###
signif<-t(sapply(pois_fits, function(u) as.numeric(summary(u)$coef[,4])))

sum(apply(signif, 1, function(u) u[2] > 0.5))
# All stores significant price coef - good.

### Final store list? ###
## Remove 1161, 7022, 7025, 7030, 7035, 7042, 7055
## Fails from newtraph? ##
# "1027" "1037" "1068" "1078" "1108" "1159" "1177" "1183" "1324"
# "1381" "1389" "1406" "1469" "1471" "1514" "1525" "1533" "1542"
# "1573" "1620" "1637" "1848" "1866" "7022" "7030" "7035" "7055"
rm_store<-c(1027, 1037, 1068, 1078, 1108, 1159, 1161, 1177, 1183, 1324, 1381, 1389, 1406, 1469, 1471, 1514,
            1525, 1533, 1542, 1573, 1620, 1637, 1848, 1866, 7022, 7025, 7030, 7035, 7042, 7055)

pois_fits <- dlply(dat[!dat$store %in% rm_store,], .(store), function(x) {
  glm(data=x, mvm ~ price, family=poisson, maxit=50)
})

rmse_store<-round(as.numeric(sapply(pois_fits, rmse)),3)

qplot(x=rmse_store, xlab="RMSE", ylab="Count", fill=I("grey60"), colour=I("black"), binwidth=2)

ests<-t(sapply(pois_fits, function(u) u$coef))
df<-melt(ests)
qplot(data=df, x=value, xlab="Coefficients", ylab="Count", fill=I("grey60"), colour=I("black"), binwidth=.5) + facet_wrap(~X2, scales="free_x")



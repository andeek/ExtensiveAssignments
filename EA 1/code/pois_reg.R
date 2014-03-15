#### Poisson Regression ####
# source("code/library.R")
# source("code/data_format.R")
# source("code/helpers.R")

fit.pois <- dlply(dat, .(store), function(x) {
  glm(data=x, mvm ~ price, family=poisson, maxit=50)
})


rmse.pois<-round(as.numeric(sapply(fit.pois, rmse)),3)

### Coefficients ###
est.pois<-t(sapply(fit.pois, function(u) u$coef))

est.pois.m<-melt(est.pois)
dist.coef.pois.plot <- qplot(data=est.pois.m, x=value, xlab="Coefficients", ylab="Count", fill=I("grey60"), colour=I("black")) + facet_wrap(~Var2, scales="free_x")
rmse.pois.plot <- qplot(x=rmse.pois, xlab="RMSE", ylab="Count", fill=I("grey60"), colour=I("black"), binwidth=2)

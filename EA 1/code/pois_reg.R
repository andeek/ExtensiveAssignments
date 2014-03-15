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

### Generating Data ###
generate_data_pois<-function(model.fit, nsets=100){
  x<-model.fit$data$price
  b<-model.fit$coef
  lambda<-exp(b[1]+b[2]*x)
  dats<-NULL
  for(i in 1:length(lambda)){
    dats<-rbind(dats, rpois(nsets, lambda[i]))
  }
  return(dats)
}

gens_pois<-lapply(fit.pois, generate_data_pois, nsets=1000)

save(gens_pois, "gens_pois.rda") 

## Proportion of 0's ##
prop0<-function(dat, gen_dat){
  orig<-sum(dat$mvm==0)/length(dat$mvm)
  gens<-apply(gen_dat, 2, function(u) sum(u==0)/length(u))
  pval=sum(gens > orig)/length(gens)
  return(pval)
}

more40<-function(dat, gen_dat){
  orig<-sum(dat$mvm > 40)
  gens<-apply(gen_dat, 2, function(u) sum(u > 40))
  pval=sum(gens > orig)/length(gens)
  return(pval)
}

## Sample Variance ##
vars<-function(dat, gen_dat){
  orig<-var(dat$mvm)
  gens<-apply(gen_dat, 2, var)
  pval=sum(gens > orig)/length(gens)
  return(pval)
}

## Sample Variance ##
less4<-function(dat, gen_dat){
  orig<-sum(dat$mvm<=4)
  gens<-apply(gen_dat, 2, function(u) sum(u <= 4))
  pval=sum(gens > orig)/length(gens)
  return(pval)
}

## p-values for 4 criteria ##
z<-unique(dat$store)
pvals_prop0<-c()
for(i in 1:length(z)){
  pvals_prop0<-c(pvals_prop0, prop0(fit.pois[[i]]$data, gens_pois[[i]]))
}

sum(pvals_prop0 > 0.05 & pvals_prop0 < 0.95)

pvals_more40<-c()
for(i in 1:length(z)){
  pvals_more40<-c(pvals_more40, more40(fit.pois[[i]]$data, gens_pois[[i]]))
}

sum(pvals_more40 > 0.05 & pvals_more40 < 0.95)

pvals_var<-c()
for(i in 1:length(z)){
  pvals_var<-c(pvals_var, vars(fit.pois[[i]]$data, gens_pois[[i]]))
}

sum(pvals_var > 0.05 & pvals_var < 0.95)

pvals_less4<-c()
for(i in 1:length(z)){
  pvals_less4<-c(pvals_less4, less4(fit.pois[[i]]$data, gens_pois[[i]]))
} 

sum(pvals_less4 > 0.05 & pvals_less4 < 0.95)

crits<-cbind(pvals_prop0, pvals_more40, pvals_less4, pvals_var)

## 52 out of 161 stores can pass 1 criterion
sum(apply(crits, 1, function(u) any(u > 0.05 & u < 0.95)))

## 0 stores can pass all
sum(apply(crits, 1, function(u) all(u > 0.05 & u < 0.95)))


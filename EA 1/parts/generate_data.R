#### Models as Data Generating Mechanisms ####
# source("code/library.R")
# source("code/data_format.R")
# source("code/helpers.R")
# load("data/zip_ests.rda")
# load("data/fit.pois.rda")


### Generating Data Functions ###
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

generate_data_zip<-function(store_dat, b, nsets=100){
  ys<-NULL
  x<-store_dat$price
  p <- exp(b[1] + b[2]*x)/(1+exp(b[1] + b[2]*x))
  lambda<-exp(b[3]+b[4]*x)
  for(i in 1:length(p)){
    z<-rbinom(nsets,1,p[i])
    y<-rpois(nsets, lambda[i])
    y[z==0]<-0
    ys<-rbind(ys, y)
    i<-i+1
  }
  return(ys)
} 

### Critera ###
## Proportion of 0's ##
prop0<-function(dat, gen_dat){
  orig<-sum(dat$mvm==0)/length(dat$mvm)
  gens<-apply(gen_dat, 2, function(u) sum(u==0)/length(u))
  pval=sum(gens > orig)/length(gens)
  return(pval)
}

## Number of days with more than 40 movement ##
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

## Number of days with less than 4 movement ##
less4<-function(dat, gen_dat){
  orig<-sum(dat$mvm<=4)
  gens<-apply(gen_dat, 2, function(u) sum(u <= 4))
  pval=sum(gens > orig)/length(gens)
  return(pval)
}

### Generate Data ###
# Poisson #
gens_pois<-lapply(fit.pois, generate_data_pois, nsets=500)

# ZIP #
gens_zip<-list(NA) 
z<-unique(dat$store)
for(i in 1:length(z)){
  gens_zip[[i]]<-generate_data_zip(dat[dat$store==z[i],], as.numeric(zpois.mles[1,-1]), nsets=500)  
cat(i, "\r")
}

# ## p-values for 4 criteria ##
# z<-unique(dat$store)
# pvals_prop0<-c()
# for(i in 1:length(z)){
#   pvals_prop0<-c(pvals_prop0, prop0(fit.pois[[i]]$data, gens_pois[[i]]))
# }
# 
# sum(pvals_prop0 > 0.05 & pvals_prop0 < 0.95)
# 
# pvals_more40<-c()
# for(i in 1:length(z)){
#   pvals_more40<-c(pvals_more40, more40(fit.pois[[i]]$data, gens_pois[[i]]))
# }
# 
# sum(pvals_more40 > 0.05 & pvals_more40 < 0.95)
# 
# pvals_var<-c()
# for(i in 1:length(z)){
#   pvals_var<-c(pvals_var, vars(fit.pois[[i]]$data, gens_pois[[i]]))
# }
# 
# sum(pvals_var > 0.05 & pvals_var < 0.95)
# 
# pvals_less4<-c()
# for(i in 1:length(z)){
#   pvals_less4<-c(pvals_less4, less4(fit.pois[[i]]$data, gens_pois[[i]]))
# } 
# 
# sum(pvals_less4 > 0.05 & pvals_less4 < 0.95)
# 
# crits<-cbind(pvals_prop0, pvals_more40, pvals_less4, pvals_var)
# 
# ## 52 out of 161 stores can pass 1 criterion
# sum(apply(crits, 1, function(u) any(u > 0.05 & u < 0.95)))
# 
# ## 0 stores can pass all
# sum(apply(crits, 1, function(u) all(u > 0.05 & u < 0.95)))
# 

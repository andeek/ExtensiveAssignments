#### Models as Data Generating Mechanisms ####
# source("code/library.R")
# source("code/data_format.R")
# source("code/helpers.R")
# load("data/zip_ests.rda")
# load("data/fit.pois.rda")
# load("data/gens_pois.rda")
# load("data/gens_zip.rda")
# load("data/gens_bayes.rda")
# load("data/pvals_gen.rda")
# load("data/res1.rda")
# load("data/res2.rda")
# load("data/res3.rda")

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
# ptm <- proc.time()
# gens_pois<-lapply(fit.pois, generate_data_pois, nsets=300)
# proc.time() - ptm

# ZIP #
# z<-unique(dat$store)
# gens_zip<-list(NA) 
# for(i in 1:length(z)){
#   gens_zip[[i]]<-generate_data_zip(dat[dat$store==z[i],], as.numeric(zpois.mles[1,-1]), nsets=300)  
# cat(i, "\r")
# }

## p-values for 4 criteria ##
# Poisson #
# pvals_pois<-c()
# for(i in 1:length(z)){
#   pvals1<-c(prop0(fit.pois[[i]]$data, gens_pois[[i]]), more40(fit.pois[[i]]$data, gens_pois[[i]]),
#             vars(fit.pois[[i]]$data, gens_pois[[i]]), less4(fit.pois[[i]]$data, gens_pois[[i]]))
#   pvals_pois<-as.data.frame(rbind(pvals_pois, pvals1))
#   cat(i, "\r")
# }
# row.names(pvals_pois)<-NULL
# names(pvals_pois)<-c("Prop 0", "More 40", "Var", "Less 4")

num_pass_pois<-apply(pvals_pois, 2, function(u) sum(u > 0.05 & u < 0.95))
prop_pass_pois<-num_pass_pois/161

# Poisson #
# pvals_zip<-c()
# for(i in 1:length(z)){
#   pvals1<-c(prop0(dat[dat$store==z[i],], gens_zip[[i]]), more40(dat[dat$store==z[i],], gens_zip[[i]]),
#             vars(dat[dat$store==z[i],], gens_zip[[i]]), less4(dat[dat$store==z[i],], gens_zip[[i]]))
#   pvals_zip<-as.data.frame(rbind(pvals_zip, pvals1))
#   cat(i, "\r")
# }
# row.names(pvals_zip)<-NULL
# names(pvals_zip)<-c("Prop 0", "More 40", "Var", "Less 4")

num_pass_zip<-apply(pvals_zip, 2, function(u) sum(u > 0.05 & u < 0.95))
prop_pass_zip<-num_pass_zip/161


#### Regional Model ####
samp<-sample(1:6000, 300)

beta0<-rbind(res1$beta0[1000:3000,], res2$beta0[1000:3000,], res3$beta0[1000:3000,])[samp,]
beta1<-rbind(res1$beta1[1000:3000,], res2$beta1[1000:3000,], res3$beta1[1000:3000,])[samp,]
beta2<-rbind(res1$beta2[1000:3000,], res2$beta2[1000:3000,], res3$beta2[1000:3000,])[samp,]
beta3<-rbind(res1$beta3[1000:3000,], res2$beta3[1000:3000,], res3$beta3[1000:3000,])[samp,]

generate_data_bayes<-function(store_dat, b0, b1, b2, b3){
  ys<-NULL
  x<-store_dat$price
  p <- apply(data.frame(b0, b1, b2, b3), 1, function(u) exp(u[1] + u[2]*x)/(1+exp(u[3] + u[4]*x)))
  lambda<-apply(data.frame(b2, b3), 1, function(u) exp(u[1]+u[2]*x))
  for(i in 1:nrow(p)){
    z<-apply(matrix(1:length(b0),ncol=1), 1, function(u) rbinom(1,1,min(p[i,u],1)))
    y<-apply(matrix(1:length(b0),ncol=1), 1, function(u) rpois(1,lambda[i,u]))
    y[z==0]<-0
    ys<-rbind(ys, y)
  }
  return(ys)
} 

gens_bayes<-list(NA) 
for(i in 1:length(z)){
  gens_bayes[[i]]<-generate_data_bayes(dat[dat$store==z[i],], beta0[,i], beta1[,i], beta2[,i], beta3[,i])  
  cat(i, "\r")
}

# Bayes #
pvals_bayes<-c()
for(i in 1:length(z)){
  pvals1<-c(prop0(dat[dat$store==z[i],], gens_bayes[[i]]), more40(dat[dat$store==z[i],], gens_bayes[[i]]),
            vars(dat[dat$store==z[i],], gens_bayes[[i]]), less4(dat[dat$store==z[i],], gens_bayes[[i]]))
  pvals_bayes<-as.data.frame(rbind(pvals_bayes, pvals1))
  cat(i, "\r")
}
row.names(pvals_bayes)<-NULL
names(pvals_bayes)<-c("Prop 0", "More 40", "Var", "Less 4")

num_pass_bayes<-apply(pvals_bayes, 2, function(u) sum(u > 0.05 & u < 0.95))
prop_pass_bayes<-num_pass_bayes/161

# source('../code/helpers.R')
# source('../code/library.R')
# source('../code/data_format.R')
load('../data/gens_pois.rda')
load('../data/gens_zip.rda')
load('../data/gens_bayes.rda')
load('../data/pvals_gen.rda')

num_pass_pois<-apply(pvals_pois, 2, function(u) sum(u > 0.05 & u < 0.95))
prop_pass_pois<-num_pass_pois/161

num_pass_zip<-apply(pvals_zip, 2, function(u) sum(u > 0.05 & u < 0.95))
prop_pass_zip<-num_pass_zip/161

num_pass_bayes<-apply(pvals_bayes, 2, function(u) sum(u > 0.05 & u < 0.95))
prop_pass_bayes<-num_pass_bayes/161

num_pass<-rbind(num_pass_pois, num_pass_zip, num_pass_bayes)
row.names(num_pass)<-c("Poisson", "ZIP", "Regional")

prop_pass<-rbind(prop_pass_pois, prop_pass_zip, prop_pass_bayes)
row.names(prop_pass)<-c("Poisson", "ZIP", "Regional")

#### Plot Stores ####
prop0_dat<-function(dat, gen_dat){
  orig<-sum(dat$mvm==0)/length(dat$mvm)
  gens<-apply(gen_dat, 2, function(u) sum(u==0)/length(u))
  return(list(orig, gens))
}

## Number of days with more than 40 movement ##
more40_dat<-function(dat, gen_dat){
  orig<-sum(dat$mvm > 40)
  gens<-apply(gen_dat, 2, function(u) sum(u > 40))
  return(list(orig, gens))
}

## Sample Variance ##
vars_dat<-function(dat, gen_dat){
  orig<-var(dat$mvm)
  gens<-apply(gen_dat, 2, var)
  return(list(orig, gens))
}

## Number of days with less than 4 movement ##
less4_dat<-function(dat, gen_dat){
  orig<-sum(dat$mvm<=4)
  gens<-apply(gen_dat, 2, function(u) sum(u <= 4))
  return(list(orig, gens))
}

store4<-samp_store(4,2,dat)

z<-unique(store4$store)

dat40_pois<-list(NA)
for(i in 1:length(z)){
  dat40_pois[[i]] <- more40_dat(store4[store4$store==z[i],], gens_pois[[i]])
}

dat40_zip<-list(NA)
for(i in 1:length(z)){
  dat40_zip[[i]] <- more40_dat(store4[store4$store==z[i],], gens_zip[[i]])
}

dat40_bayes<-list(NA)
for(i in 1:length(z)){
  dat40_bayes[[i]] <- more40_dat(store4[store4$store==z[i],], gens_bayes[[i]])
}

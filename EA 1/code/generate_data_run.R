#load('../data/gens_pois.rda')
#load('../data/gens_zip.rda')
#load('../data/gens_bayes.rda')
load('./data/pvals_gen.rda')
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




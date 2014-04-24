### Chi-square test ###
phi.est<-ests[1]
l1<-ests[2]
l2<-ests[3]
subway.dat$expected<-probs*l1*subway.dat$length + (1-probs)*l2*subway.dat$length

chi.test<-sum((subway.dat$flaws-subway.dat$expected)^2/subway.dat$expected)
chi.test

## Not valid ##
sums.obs.exp<-ddply(subway.dat, .(length), summarize, exp=sum(expected), obs=sum(flaws))
chi.test<-sum((sums.obs.exp$obs-sums.obs.exp$exp)^2/sums.obs.exp$exp)
chi.test

### Bootstrap CI ###

### Sample from model ###
n<-subway.dat$length
ests<-NULL
for(j in 1:1000){
  phis<-runif(75)
  y<-1:75
  for(i in 1:75){
    if(phis[i] < phi.est){
      y[i]<-rpois(1, l1*n[i])
    }else{
      y[i]<-rpois(1, l2*n[i])
    }
  }
  ests<-rbind(ests, c(em_algorithm(.5, 5, 0.01, y=y, n=subway.dat$length)[1:3]))
}
dim(ests)

boot_intervals<-t(apply(ests, 2, function(u) sort(u)[c(25, 975)]))
boot_intervals
### Plots MCMC for Regional Model ###
#source('/code/library.R')
#source('/code/data_format.R')
#source('/code/helpers.R')
#load('../data/res1.rda')
#load('../data/res2.rda')
#load('../data/res3.rda')
load('../data/bayes_ests.rda')

#### trace plots ####
stores4<-samp_store(4, 2, dat)
store1<-stores4[stores4$store==1418,]

#i=1
#beta0.list = mcmc.list(mcmc(res1$beta0[,i]),
#                       mcmc(res2$beta0[,i]),
#                       mcmc(res3$beta0[,i]))
#beta1.list = mcmc.list(mcmc(res1$beta1[,i]),
#                       mcmc(res2$beta1[,i]),
#                       mcmc(res3$beta1[,i]))
#beta2.list = mcmc.list(mcmc(res1$beta2[,i]),
#                       mcmc(res2$beta2[,i]),
#                       mcmc(res3$beta2[,i]))
#beta3.list = mcmc.list(mcmc(res1$beta3[,i]),
#                       mcmc(res2$beta3[,i]),
#                       mcmc(res3$beta3[,i]))

#par(mfrow=c(2,2))
#plot(beta0.list, smooth=F, density=F, ylim=c(3,12), auto.layout=F, main=expression(beta[paste("01")]), lwd=2, ylab="Estimate")
#plot(beta1.list, smooth=F, density=F, ylim=c(-16, -5), auto.layout=F, main=expression(beta[11]), lwd=2, ylab="Estimate")
#plot(beta2.list, smooth=F, density=F, auto.layout=F, main=expression(beta[21]), lwd=2, ylab="Estimate")
#plot(beta3.list, smooth=F, density=F, auto.layout=F, main=expression(beta[31]), lwd=2, ylab="Estimate")

#sigma0.list = mcmc.list(mcmc(res1$sigma[,1]),
#                       mcmc(res2$sigma[,1]),
#                       mcmc(res3$sigma[,1]))
#sigma1.list = mcmc.list(mcmc(res1$sigma[,2]),
#                        mcmc(res2$sigma[,2]),
#                        mcmc(res3$sigma[,2]))
#sigma2.list = mcmc.list(mcmc(res1$sigma[,3]),
#                        mcmc(res2$sigma[,3]),
#                        mcmc(res3$sigma[,3]))
#sigma3.list = mcmc.list(mcmc(res1$sigma[,4]),
#                        mcmc(res2$sigma[,4]),
#                        mcmc(res3$sigma[,4]))

#par(mfrow=c(2,2))
#plot(sigma0.list, smooth=F, density=F, auto.layout=F, main=expression(sigma[paste("0")]), lwd=2, ylab="Estimate")
#plot(sigma1.list, smooth=F, density=F, auto.layout=F, main=expression(sigma[1]), lwd=2, ylab="Estimate")
#plot(sigma2.list, smooth=F, density=F, auto.layout=F, main=expression(sigma[2]), lwd=2, ylab="Estimate")
#plot(sigma3.list, smooth=F, density=F, auto.layout=F, main=expression(sigma[3]), lwd=2, ylab="Estimate")

names(est.bayes)<-c("beta_0", "beta_1", "beta_2", "beta_3")
bayes.melt<-melt(est.bayes)
bayes.coef.dist<-qplot(data=bayes.melt, x=value, fill=I("grey40"), colour=I("black"), xlab="Coefficient", ylab="Count") + facet_wrap(~Var2, scales="free_x")
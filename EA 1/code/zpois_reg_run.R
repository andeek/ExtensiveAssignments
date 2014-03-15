load("data/fit.zpois.rda")

zpois.mles <- ldply(zpois.est, function(x) if(length(x) == 3) t(x[[1]]))
zpois.mles.m <- melt(zpois.mles, id.vars="store")
dist.coef.zpois.plot <- qplot(data=zpois.mles.m, x=value, xlab="Coefficients", ylab="Count", fill=I("grey60"), colour=I("black")) + facet_wrap(~variable, scales="free_x")


### create plots for viewing fit
# plot_by_store <- function(store_id, data=dat) {
#   dat1 <- subset(data, store == store_id)
#   
#   est1 <- est[[as.character(dat1$store[1])]][[1]]
#   dat1$p <- exp(est1[1] + est1[2]*dat1$price)/(1 + exp(est1[1] + est1[2]*dat1$price))
#   dat1$lambda <- exp(est1[3] + est1[4]*dat1$price)
#   
#   distns <- ldply(mlply(unique(dat1$price), function(y) mdply(0:max(subset(dat1, price == y)$mvm), function(x) c(mvm=x, price=y, lik=exp(loglik.zip(est1, dat=list(price=y, mvm=x)))))), function(x) x)
#   
#   gg <- ggplot(dat1) + geom_histogram(aes(x=mvm, y=..density..)) + 
#     geom_line(aes(x=mvm, y=lik), data=distns, colour=I("red")) + facet_wrap(~price, scales="free")  +
#     ggtitle(sprintf("Store ID: %d", store_id))
#   
#   ##ggsave(sprintf("EA 1/figure/fit_store_%d.pdf", store_id), gg)
# }
# 
# mlply(1:length(est), function(x) {
#   if(length(est[[x]]) == 3) {
#     cat(sprintf("Counter: %d \r", x))
#     plot_by_store(as.numeric(names(est)[x]))     
#   }
# })
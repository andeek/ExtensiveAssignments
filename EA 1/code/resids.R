source("code/pois_reg_run.R")
source("code/zipois_reg_run.R")

est.pois <- as.data.frame(est.pois)
est.pois$store <- as.character(rownames(est.pois))
names(est.pois) <- c("b_0", "b_1", "store")

### Poisson Deviance
input <- merge(dat, est.pois, all.x=TRUE, by="store")
input$lambda <-with(input, exp(b_0 + b_1*price ))

D_pois <- ddply(input, .(store), function(input) 2*sum(with(input, mvm*log(ifelse(mvm == 0, 1, mvm/lambda)) - (mvm - lambda))))

D_resids_pois <- dlply(input, .(store), function(input){
  with(input, sign(mvm - lambda)*sqrt(2*(mvm*log(ifelse(mvm == 0, 1, mvm/lambda)) - (mvm - lambda))))
})

### ZIP Deviance Resids
names(zpois.mles) <- c("store", "b_0", "b_1", "b_2", "b_3")
input.zip <- merge(dat, zpois.mles, all.x=TRUE, by="store")
input.zip$p <- with(input.zip, exp(b_0 + b_1*price)/(1+exp(b_0 + b_1*price)))
input.zip$lambda <-with(input.zip, exp(b_2 + b_3*price ))
input.zip$mu <- with(input.zip, lambda*p)
input.zip$delta <- as.numeric(input$mvm == 0)


D_resids_zip <- dlply(input.zip, .(store), function(input){
  with(input, sign(mvm - lambda)*sqrt(2*(delta*log(ifelse(lambda - mvm + mvm*exp(-lambda) <=0, 1, (lambda - mvm + mvm*exp(-lambda))/(lambda - mu + mu*exp(-lambda)))) + (1-delta)*log(ifelse(mvm == 0, 1, mvm/mu))*(log(ifelse(mvm == 0, 1, mvm/mu))>0))))
})

plot_D_resids <- function(store_id, deviance){
  store <- as.character(store_id)
  g <- qplot(input[input$store == store, "date"], deviance[[store]]) +
    xlab("Date") + ylab("Deviance Residuals") + ggtitle(sprintf("Store: %s", store))
  print(g)
}

plot_D_resids(1009, D_resids_pois)

qplot(price, (mvm-mu)/sqrt(1+lambda+mu), data=input.zip[input.zip$store==1481,]) +
geom_point(aes(price+.001, (mvm-lambda)/sqrt(lambda)), colour="red", data=input[input$store==1481,])


with(input.zip[input.zip$store==1009,], sum(((mvm-mu)/sqrt(1+lambda+mu))^2)/(length(mvm)-1))
with(input[input$store==1009,], sum(((mvm-lambda)/sqrt(lambda))^2)/(length(mvm)-1))


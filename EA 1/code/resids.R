### Residuals ###
# source("../code/pois_reg_run.R")
# source("../code/zipois_reg_run.R")

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

storeid_4 <- c(1065, 1075, 1416, 1481)
                    
stand_resid_plot_price <- ggplot() + 
  geom_jitter(aes(price, (mvm-lambda)/sqrt(lambda)), colour="grey", data=input[input$store%in% storeid_4,]) +
  geom_jitter(aes(price, (mvm-mu)/sqrt(1+lambda+mu)), colour="black", shape=17, data=input.zip[input.zip$store %in% storeid_4,]) +
  facet_wrap(~store, nrow=2) +
  xlab("Price") + ylab("Standardized Residuals")

stand_resid_plot_time <- ggplot() + 
  geom_jitter(aes(date, (mvm-lambda)/sqrt(lambda)), colour="grey", data=input[input$store%in% storeid_4,]) +
  geom_jitter(aes(date, (mvm-mu)/sqrt(1+lambda+mu)), colour="black", shape=17, data=input.zip[input.zip$store %in% storeid_4,]) +
  facet_wrap(~store, nrow=2) +
  xlab("Time") + ylab("Standardized Residuals")

MSE_standard <- cbind(
ddply(input.zip[input.zip$store%in% storeid_4,], .(store), summarise, MSE=sum(((mvm-mu)/sqrt(1+lambda+mu))^2)/(length(mvm)-1)),
ddply(input[input$store%in% storeid_4,], .(store), summarise, MSE=sum(((mvm-lambda)/sqrt(lambda))^2)/(length(mvm)-1))[,2])

names(MSE_standard) <- c("Store", "ZIP", "Pois")

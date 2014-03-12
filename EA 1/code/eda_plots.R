### Code for EDA ###
## Basic summaries
dat_store_s <- ddply(dat, .(store), summarise, 
                     count = length(price), 
                     num_sales = sum(mvm))

plot_nrow_hist <- ggplot(dat_store_s) + geom_histogram(aes(count)) + xlab("Number of observations") + 
  theme_bw()

plot_nsales_hist <- ggplot(dat_store_s) + geom_histogram(aes(num_sales)) + xlab("Number of units sold") + 
  theme_bw()

plot_density <- ggplot(aes(value, group=store), data=melt(dat[,-2], id.vars="store")) + 
  geom_line(alpha=.1, stat='density') + 
  facet_wrap(~variable, scales='free') + 
  theme_bw()


#### Sample Stores ####
samp_store<-function(n, seed, dat){
  set.seed(seed)
  s<-sample(unique(dat$store),n)
  dat[dat$store %in% s,]
}

#### Plot Price vs. Volume ####
## Sample 4 stores
stores4<-samp_store(4, 2, dat=dat)
stores4$store<-paste("Store ID:", stores4$store)

## Plot price vs. sales for 4 stores
stores4_m <- merge(x=stores4, y=ddply(stores4, .(store), summarise, max_mvm = max(mvm)), all.x = TRUE, by="store")
plot_four_stores <- ggplot(stores4_m, ) + 
  xlab("Date") + ylab("Scaled Movement vs. Price") + 
  geom_line(aes(date, mvm/max_mvm, colour="line")) + 
  geom_point(aes(date, price, colour="point")) +
  theme_bw() + facet_wrap(~store, nrow=1) +
  scale_colour_manual(values=c("blue", "red"),
                      labels=c("Scaled Movement", "Price"),
                      name="") +
  theme(legend.position="bottom")

plot_hist_4 <- ggplot(stores4) +
  geom_histogram(aes(mvm)) +
  facet_wrap(~store, nrow=2, scales="free") +
  theme_bw()

plot_price_mvm_4 <- ggplot(merge(x=stores4, y=ddply(stores4, .(store, price), summarise, med_mvm = median(mvm)), all.x = TRUE, by="store")) +
  geom_point(aes(x=price.x, y=mvm, colour=store)) + 
  facet_wrap(~store, nrow=1) +
  geom_line(aes(x=price.y, y=med_mvm, colour=store)) +
  theme_bw() +
  theme(legend.position="none")

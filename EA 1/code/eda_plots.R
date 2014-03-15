### Code for EDA ###
# source("code/library.R")
# source("code/data_format.R")

## Basic summaries
dat_store_s <- ddply(dat.orig, .(store), summarise, 
                     count = length(price), 
                     num_sales = sum(mvm))

plot_nrow_hist <- ggplot(dat_store_s) + geom_histogram(aes(count), fill=I("grey60"), colour=I("black")) + xlab("Number of observations")

plot_nsales_hist <- ggplot(dat_store_s) + geom_histogram(aes(num_sales), fill=I("grey60"), colour=I("black")) + xlab("Number of units sold")

plot_density <- ggplot(aes(value, group=store), data=melt(dat.orig[,-c(2, 3, 6, 7)], id.vars="store")) + 
  geom_line(alpha=.1, stat='density') + 
  facet_wrap(~variable, scales='free')

plot_density_zoom1 <- ggplot(aes(mvm, group=store), data=dat.orig[, c(1,4)]) + 
  geom_line(alpha=.1, stat='density') +
  xlim(c(0,1))
plot_density_zoom2 <- ggplot(aes(mvm, group=store), data=dat.orig[, c(1,4)]) + 
  geom_line(alpha=.1, stat='density') +
  xlim(c(1,25))

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
  facet_wrap(~store, nrow=1) +
  scale_colour_manual(values=c("blue", "red"),
                      labels=c("Scaled Movement", "Price"),
                      name="") +
  theme(legend.position="bottom")

plot_hist_4 <- ggplot(stores4) +
  geom_histogram(aes(mvm), fill=I("grey60"), colour=I("black")) +
  facet_wrap(~store, nrow=2, scales="free")

plot_price_mvm_4 <- ggplot(merge(x=stores4, y=ddply(stores4, .(store, price), summarise, med_mvm = median(mvm)), all.x = TRUE, by="store")) +
  geom_point(aes(x=price.x, y=mvm, colour=store)) + 
  facet_wrap(~store, nrow=1) +
  geom_line(aes(x=price.y, y=med_mvm, colour=store)) +
  theme(legend.position="none")

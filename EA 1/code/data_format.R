### Script for formatting data ###
## read file
dat.orig <- read.table("data/greenbeandat.txt", header=TRUE)

## remove unused columns
dat <- dat.orig[,-c(2, 6, 7)]

## force date to useful date format
dat$date <-  ymd(as.character(dat$date))

## rm unusable stores
rm_store<-c(1027, 1037, 1068, 1078, 1108, 1159, 1161, 1177, 1183, 1324, 1381, 1389, 1406, 1469, 1471, 1514,
            1525, 1533, 1542, 1573, 1620, 1637, 1848, 1866, 7022, 7025, 7030, 7035, 7042, 7055)

dat <- dat[!dat$store %in% rm_store,]
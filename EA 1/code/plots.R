#### Libraries ####
library(lubridate)
library(ggplot2)
library(reshape)
library(gtable)
library(grid)

setwd("./GitHub/ExtensiveAssignments")

#### Data ####
dat<-read.table(file='EA 1/data/greenbeandat.txt', header=TRUE)[,c(1,3,4,5)]
dat$date <-  ymd(as.character(dat$date))

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

# two plots
p1 <- ggplot(stores4, aes(date, mvm)) + xlab("Date") + ylab("Movement") + geom_line(colour = "blue") + theme_bw() + facet_wrap(~store, ncol=2)
p2 <- ggplot(stores4, aes(date, price)) + xlab("Date") + ylab("Movement") + geom_point(colour = "red", size=1) + facet_wrap(~store, ncol=2) + theme_bw()  %+replace% 
  theme(panel.background = element_rect(fill = NA))

# extract gtable
g1 <- ggplot_gtable(ggplot_build(p1))
g2 <- ggplot_gtable(ggplot_build(p2))

# overlap the panel of 2nd plot on that of 1st plot
pp <- c(subset(g1$layout, name == "panel-1", se = t:r))
g <- gtable_add_grob(g1, g2$grobs[which(g2$layout$name == "panel-1")], pp$t, 
                     pp$l, pp$b, pp$l)

# overlap the panel of 2nd plot on that of 1st plot
pp <- c(subset(g1$layout, name == "panel-2", se = t:r))
g <- gtable_add_grob(g, g2$grobs[which(g2$layout$name == "panel-2")], pp$t, 
                     pp$l, pp$b, pp$l)

# overlap the panel of 2nd plot on that of 1st plot
pp <- c(subset(g1$layout, name == "panel-3", se = t:r))
g <- gtable_add_grob(g, g2$grobs[which(g2$layout$name == "panel-3")], pp$t, 
                     pp$l, pp$b, pp$l)

# overlap the panel of 2nd plot on that of 1st plot
pp <- c(subset(g1$layout, name == "panel-4", se = t:r))
g <- gtable_add_grob(g, g2$grobs[which(g2$layout$name == "panel-4")], pp$t, 
                     pp$l, pp$b, pp$l)


ia <- which(g2$layout$name == "axis_l-1")
ga <- g2$grobs[[ia]]
ax <- ga$children[[2]]
ax$widths <- rev(ax$widths)
ax$grobs <- rev(ax$grobs)
ax$grobs[[1]]$x <- ax$grobs[[1]]$x - unit(1, "npc") + unit(0.15, "cm")
g <- gtable_add_cols(g, g2$widths[g2$layout[ia, ]$l], length(g$widths) - 1)
g <- gtable_add_grob(g, ax, pp$t, length(g$widths) - 1, pp$b)

# draw it
grid.draw(g)

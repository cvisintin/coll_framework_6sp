require(maptools)
require(raster)
require(xtable)
require(ggplot2)
require(ncf)
require(doMC)
require(data.table)

invcloglog <- function (x) {1-exp(-exp(x))}

plotPal <- c("#94d1c7", "#cccc2b", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#969696", "#bc80bd")

species.table <- read.delim("data/species_list.csv", header=T, sep=",")

load("data/coll_model_data")
load("output/coll_glm")

occ <- NULL
for (i in 1:nrow(species.table)) {
  data <- model.data[[i]]
  colnames(data) <- c("uid","occ","tvol","tspd","coll")
  model <- coll.glm[[i]]
  temp_df <- data.frame(x=data[,occ], y=invcloglog(cbind(1,log(data[,occ]),mean(log(data[,tvol])),mean((log(data[,tvol]))*(log(data[,tvol]))),mean(log(data[,tspd]))) %*% coef(model)), col=rep(toupper(species.table[i,2]), each=length(data[,occ])))
  occ <- rbind(occ,temp_df)
  rm(data)
  rm(model)
  rm(temp_df)
}  

tiff('figs/occ.tif', pointsize = 8, compression = "lzw", res=300, width = 1100, height = 900)
ggplot(occ,aes(x=x,y=y,group=col,colour=factor(col))) +
  geom_line(size=0.3) +
  ylab("Likelihood of Collision") +
  xlab("Likelihood of Species Occurrence") +
  labs(color = "Species") +
  theme_bw() +
  theme(legend.key = element_blank()) +
  theme(plot.margin=unit(c(.5,0,.1,.1),"cm")) +
  theme(axis.title.x = element_text(margin=unit(c(.3,0,0,0),"cm"))) +
  theme(axis.title.y = element_text(margin=unit(c(0,.3,0,0),"cm"))) +
  theme(panel.grid.major = element_line(size=0.1),panel.grid.minor = element_line(size=0.1)) +
  scale_colour_manual(values=plotPal) +
  theme(text = element_text(size = 8)) +
  scale_x_continuous(breaks=seq(0,1,by=.1), expand = c(0, 0), lim=c(0,1)) +
  scale_y_continuous(breaks=seq(0,1,by=.1), expand = c(0, 0), lim=c(0,1)) #+
  #guides(colour=FALSE)
dev.off()

tvol <- NULL
for (i in 1:nrow(species.table)) {
  data <- model.data[[i]]
  colnames(data) <- c("uid","occ","tvol","tspd","coll")
  model <- coll.glm[[i]]
  temp_df <- data.frame(x=data[,tvol], y=invcloglog(cbind(1,mean(log(data[,occ])),log(data[,tvol]),(log(data[,tvol]))*(log(data[,tvol])),mean(log(data[,tspd]))) %*% coef(model)), col=rep(toupper(species.table[i,2]), each=length(data[,tvol])))
  tvol <- rbind(tvol,temp_df)
  rm(data)
  rm(model)
  rm(temp_df)
}

tiff('figs/tvol.tif', pointsize = 8, compression = "lzw", res=300, width = 1100, height = 900)
ggplot(tvol,aes(x=x/1000,y=y,group=col,colour=factor(col))) +
  geom_line(size=0.3) +
  ylab("Likelihood of Collision") +
  xlab("Traffic Volume (1000 vehicles/day)") +
  labs(color = "Species") +
  theme_bw() +
  theme(legend.key = element_blank()) +
  theme(plot.margin=unit(c(.5,0,.1,.1),"cm")) +
  theme(axis.title.x = element_text(margin=unit(c(.3,0,0,0),"cm"))) +
  theme(axis.title.y = element_text(margin=unit(c(0,.3,0,0),"cm"))) +
  theme(panel.grid.major = element_line(size=0.1),panel.grid.minor = element_line(size=0.1)) +
  scale_colour_manual(values=plotPal) +
  theme(text = element_text(size = 8)) +
  scale_x_continuous(breaks=seq(0,25,by=5), expand = c(0, 0), lim=c(0,25)) +
  scale_y_continuous(breaks=seq(0,1,by=.1), expand = c(0, 0), lim=c(0,1)) #+
  #guides(colour=FALSE)
dev.off()

tspd <- NULL
for (i in 1:nrow(species.table)) {
  data <- model.data[[i]]
  colnames(data) <- c("uid","occ","tvol","tspd","coll")
  model <- coll.glm[[i]]
  temp_df <- data.frame(x=data[,tspd], y=invcloglog(cbind(1,mean(log(data[,occ])),mean(log(data[,tvol])),mean((log(data[,tvol]))*(log(data[,tvol]))),log(data[,tspd])) %*% coef(model)), col=rep(toupper(species.table[i,2]), each=length(data[,tspd])))
  tspd <- rbind(tspd,temp_df)
  rm(data)
  rm(model)
  rm(temp_df)
}  

tiff('figs/tspd.tif', pointsize = 8, compression = "lzw", res=300, width = 1100, height = 900)
ggplot(tspd,aes(x=x,y=y,group=col,colour=factor(col))) +
  geom_line(size=0.3) +
  ylab("Likelihood of Collision") +
  xlab("Traffic Speed (km/hour)") +
  labs(color = "Species") +
  theme_bw() +
  theme(legend.key = element_blank()) +
  theme(plot.margin=unit(c(.5,0,.1,.1),"cm")) +
  theme(axis.title.x = element_text(margin=unit(c(.3,0,0,0),"cm"))) +
  theme(axis.title.y = element_text(margin=unit(c(0,.3,0,0),"cm"))) +
  theme(panel.grid.major = element_line(size=0.1),panel.grid.minor = element_line(size=0.1)) +
  scale_colour_manual(values=plotPal) +
  theme(text = element_text(size = 8)) +
  scale_x_continuous(breaks=seq(40,110,by=10), expand = c(0, 0), lim=c(40,110)) +
  scale_y_continuous(breaks=seq(0,1,by=.1), expand = c(0, 0), lim=c(0,1)) #+
  #guides(colour=FALSE)
dev.off()


#calculate spatial autocorrelation across all species
registerDoMC(detectCores() - 1)

auto.resid <- foreach(i = 1:nrow(species.table), .packages = c("ncf"), .combine="rbind") %dopar% {
  data <- as.data.frame(model.data[[i]])
  model <- coll.glm[[i]]
  cor <- correlog(data[,2], data[,3], resid(model), increment=1000, resamp=0, latlon=FALSE)
  temp_df <- data.frame(x=as.numeric(names(cor$correlation[2:20])), y=cor$correlation[2:20], col=rep(toupper(species.table[i,2]), each=length(cor$correlation[2:20])))
  temp_df
}  

auto.coll <- foreach(i = 1:nrow(species.table), .packages = c("ncf"), .combine="rbind") %dopar% {
  data <- as.data.frame(model.data[[i]])
  cor <- correlog(data[,2], data[,3], data[,7], increment=1000, resamp=0, latlon=FALSE)
  temp_df <- data.frame(x=as.numeric(names(cor$correlation[2:20])), y=cor$correlation[2:20], col=rep(toupper(species.table[i,2]), each=length(cor$correlation[2:20])))
  temp_df
}  

shapes <- unlist(lapply(c("1", "2", "3", "4", "5", "6"), utf8ToInt))

ggplot(auto.resid,aes(x=x,y=y,group=col,shape=col)) + geom_line(colour=c("grey70"),size=.75) + geom_point(size=2.5) + ylab("Moran's I") + xlab("Distance (km)") + labs(shape = "Species") + theme_bw() + theme(legend.key = element_blank()) + scale_colour_manual(values=plotPal) + scale_shape_manual(values=shapes) + geom_hline(aes(yintercept=0), linetype=2) + scale_x_continuous(breaks=seq(1, 20, 1))

ggplot(auto.coll,aes(x=x,y=y,group=col,shape=col)) + geom_line(colour=c("grey70"),size=.75) + geom_point(size=2.5) + ylab("Moran's I") + xlab("Distance (km)") + labs(shape = "Species") + theme_bw() + theme(legend.key = element_blank()) + scale_colour_manual(values=plotPal) + scale_shape_manual(values=shapes) + geom_hline(aes(yintercept=0), linetype=2) + scale_x_continuous(breaks=seq(1, 20, 1))

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

species.names <- c("Eastern Grey Kangaroo","Common Brushtail Possum","Common Ringtail Possum","Black Swamp Wallaby","Common Wombat","Koala")

load("data/coll_model_data")
load("output/coll_glm")

occ <- NULL
for (i in 1:nrow(species.table)) {
  data <- model.data[[i]]
  colnames(data) <- c("uid","x","y","occ","tvol","tspd","coll")
  model <- coll.glm[[i]]
  temp_df <- data.frame(x=data[,occ], y=invcloglog(cbind(1,log(data[,occ]),mean(log(data[,tvol])),mean((log(data[,tvol]))*(log(data[,tvol]))),mean(log(data[,tspd]))) %*% coef(model)), name=rep(paste(species.names[i]), each=length(data[,occ])))
  occ <- rbind(occ,temp_df)
  rm(data)
  rm(model)
  rm(temp_df)
}  

tiff('figs/occ.tif', pointsize = 12, compression = "lzw", res=300, width = 1440, height = 900)
ggplot(occ,aes(x=x,y=y,group=name,colour=factor(name))) +
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
  theme(text = element_text(size = 10)) +
  scale_x_continuous(breaks=seq(0,1,by=.1), expand = c(0, 0), lim=c(0,1)) +
  scale_y_continuous(breaks=seq(0,1,by=.1), expand = c(0, 0), lim=c(0,1)) #+
  #guides(colour=FALSE)
dev.off()

tvol <- NULL
for (i in 1:nrow(species.table)) {
  data <- model.data[[i]]
  colnames(data) <- c("uid","x","y","occ","tvol","tspd","coll")
  model <- coll.glm[[i]]
  temp_df <- data.frame(x=data[,tvol], y=invcloglog(cbind(1,mean(log(data[,occ])),log(data[,tvol]),(log(data[,tvol]))*(log(data[,tvol])),mean(log(data[,tspd]))) %*% coef(model)), name=rep(paste(species.names[i]), each=length(data[,tvol])))
  tvol <- rbind(tvol,temp_df)
  rm(data)
  rm(model)
  rm(temp_df)
}

tiff('figs/tvol.tif', pointsize = 12, compression = "lzw", res=300, width = 900, height = 900)
ggplot(tvol,aes(x=x/1000,y=y,group=name,colour=factor(name))) +
  geom_line(size=0.3) +
  ylab("Likelihood of Collision") +
  xlab("Traffic Volume (1000 vehicles/day)") +
  labs(color = "Species") +
  theme_bw() +
  theme(legend.key = element_blank(), legend.position="none") +
  theme(plot.margin=unit(c(.5,.5,.1,.1),"cm")) +
  theme(axis.title.x = element_text(margin=unit(c(.3,0,0,0),"cm"))) +
  theme(axis.title.y = element_text(margin=unit(c(0,.3,0,0),"cm"))) +
  theme(panel.grid.major = element_line(size=0.1),panel.grid.minor = element_line(size=0.1)) +
  scale_colour_manual(values=plotPal) +
  theme(text = element_text(size = 10)) +
  scale_x_continuous(breaks=seq(0,25,by=5), expand = c(0, 0), lim=c(0,25)) +
  scale_y_continuous(breaks=seq(0,1,by=.1), expand = c(0, 0), lim=c(0,1)) #+
  #guides(colour=FALSE)
dev.off()

tspd <- NULL
for (i in 1:nrow(species.table)) {
  data <- model.data[[i]]
  colnames(data) <- c("uid","x","y","occ","tvol","tspd","coll")
  model <- coll.glm[[i]]
  temp_df <- data.frame(x=data[,tspd], y=invcloglog(cbind(1,mean(log(data[,occ])),mean(log(data[,tvol])),mean((log(data[,tvol]))*(log(data[,tvol]))),log(data[,tspd])) %*% coef(model)), name=rep(paste(species.names[i]), each=length(data[,tspd])))
  tspd <- rbind(tspd,temp_df)
  rm(data)
  rm(model)
  rm(temp_df)
}  

tiff('figs/tspd.tif', pointsize = 12, compression = "lzw", res=300, width = 900, height = 900)
ggplot(tspd,aes(x=x,y=y,group=name,colour=factor(name))) +
  geom_line(size=0.3) +
  ylab("Likelihood of Collision") +
  xlab("Traffic Speed (km/hour)") +
  labs(color = "Species") +
  theme_bw() +
  theme(legend.key = element_blank(), legend.position="none") +
  theme(plot.margin=unit(c(.5,.5,.1,.1),"cm")) +
  theme(axis.title.x = element_text(margin=unit(c(.3,0,0,0),"cm"))) +
  theme(axis.title.y = element_text(margin=unit(c(0,.3,0,0),"cm"))) +
  theme(panel.grid.major = element_line(size=0.1),panel.grid.minor = element_line(size=0.1)) +
  scale_colour_manual(values=plotPal) +
  theme(text = element_text(size = 10)) +
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
  temp_df <- data.frame(x=as.numeric(names(cor$correlation[2:20])), y=cor$correlation[2:20], name=rep(paste(species.names[i]), each=length(cor$correlation[2:20])))
  temp_df
}  

auto.coll <- foreach(i = 1:nrow(species.table), .packages = c("ncf"), .combine="rbind") %dopar% {
  data <- as.data.frame(model.data[[i]])
  cor <- correlog(data[,2], data[,3], data[,7], increment=1000, resamp=0, latlon=FALSE)
  temp_df <- data.frame(x=as.numeric(names(cor$correlation[2:20])), y=cor$correlation[2:20], name=rep(paste(species.names[i]), each=length(cor$correlation[2:20])))
  temp_df
}  

shapes <- unlist(lapply(c("1", "2", "3", "4", "5", "6"), utf8ToInt))

ggplot(auto.resid,aes(x=x,y=y,group=name,shape=name)) + geom_line(colour=c("grey70"),size=.75) + geom_point(size=2.5) + ylab("Moran's I") + xlab("Distance (km)") + labs(shape = "Species") + theme_bw() + theme(legend.key = element_blank()) + scale_colour_manual(values=plotPal) + scale_shape_manual(values=shapes) + geom_hline(aes(yintercept=0), linetype=2) + scale_x_continuous(breaks=seq(1, 20, 1))

ggplot(auto.coll,aes(x=x,y=y,group=name,shape=name)) + geom_line(colour=c("grey70"),size=.75) + geom_point(size=2.5) + ylab("Moran's I") + xlab("Distance (km)") + labs(shape = "Species") + theme_bw() + theme(legend.key = element_blank()) + scale_colour_manual(values=plotPal) + scale_shape_manual(values=shapes) + geom_hline(aes(yintercept=0), linetype=2) + scale_x_continuous(breaks=seq(1, 20, 1))


#plot marginal effects based on seasons
load("data/coll_model_data_sum")
load("data/coll_model_data_aut")
load("data/coll_model_data_win")
load("data/coll_model_data_spr")
load("output/coll_glm_summer")
load("output/coll_glm_autumn")
load("output/coll_glm_winter")
load("output/coll_glm_spring")

occ.summer <- NULL
for (i in 1:nrow(species.table)) {
  data <- model.data.summer[[i]]
  colnames(data) <- c("uid","x","y","occ","tvol","tspd","coll")
  model <- coll.glm.summer[[i]]
  temp_df <- data.frame(x=data[,occ], y=invcloglog(cbind(1,log(data[,occ]),mean(log(data[,tvol])),mean((log(data[,tvol]))*(log(data[,tvol]))),mean(log(data[,tspd]))) %*% coef(model)), name=rep(paste(species.names[i]), each=length(data[,occ])))
  occ.summer <- rbind(occ.summer,temp_df)
  rm(data)
  rm(model)
  rm(temp_df)
}
occ.summer <- cbind(occ.summer,season="summer")

occ.autumn <- NULL
for (i in 1:nrow(species.table)) {
  data <- model.data.autumn[[i]]
  colnames(data) <- c("uid","x","y","occ","tvol","tspd","coll")
  model <- coll.glm.autumn[[i]]
  temp_df <- data.frame(x=data[,occ], y=invcloglog(cbind(1,log(data[,occ]),mean(log(data[,tvol])),mean((log(data[,tvol]))*(log(data[,tvol]))),mean(log(data[,tspd]))) %*% coef(model)), name=rep(paste(species.names[i]), each=length(data[,occ])))
  occ.autumn <- rbind(occ.autumn,temp_df)
  rm(data)
  rm(model)
  rm(temp_df)
}
occ.autumn <- cbind(occ.autumn,season="autumn")


occ.winter <- NULL
for (i in 1:nrow(species.table)) {
  data <- model.data.winter[[i]]
  colnames(data) <- c("uid","x","y","occ","tvol","tspd","coll")
  model <- coll.glm.winter[[i]]
  temp_df <- data.frame(x=data[,occ], y=invcloglog(cbind(1,log(data[,occ]),mean(log(data[,tvol])),mean((log(data[,tvol]))*(log(data[,tvol]))),mean(log(data[,tspd]))) %*% coef(model)), name=rep(paste(species.names[i]), each=length(data[,occ])))
  occ.winter <- rbind(occ.winter,temp_df)
  rm(data)
  rm(model)
  rm(temp_df)
}
occ.winter <- cbind(occ.winter,season="winter")

occ.spring <- NULL
for (i in 1:nrow(species.table)) {
  data <- model.data.spring[[i]]
  colnames(data) <- c("uid","x","y","occ","tvol","tspd","coll")
  model <- coll.glm.spring[[i]]
  temp_df <- data.frame(x=data[,occ], y=invcloglog(cbind(1,log(data[,occ]),mean(log(data[,tvol])),mean((log(data[,tvol]))*(log(data[,tvol]))),mean(log(data[,tspd]))) %*% coef(model)), name=rep(paste(species.names[i]), each=length(data[,occ])))
  occ.spring <- rbind(occ.spring,temp_df)
  rm(data)
  rm(model)
  rm(temp_df)
}
occ.spring <- cbind(occ.spring,season="spring")

occ.all <- rbind(occ.summer,occ.autumn,occ.winter,occ.spring)

tiff('figs/occ_seasons.tif', pointsize = 12, compression = "lzw", res=300, width = 1415, height = 900)
ggplot(occ.all,aes(x=x,y=y,group=interaction(name, season),colour=factor(name),shape=factor(season))) +
  geom_line(size=0.3,aes(linetype=season)) +
  ylab("Likelihood of Collision") +
  xlab("Likelihood of Species Occurrence") +
  labs(colour = "Species", linetype = "Season") +
  theme_bw() +
  theme(legend.key = element_blank(),legend.key.size = unit(1.0, 'lines')) +
  theme(plot.margin=unit(c(.5,0,.1,.1),"cm")) +
  theme(axis.title.x = element_text(margin=unit(c(.3,0,0,0),"cm"))) +
  theme(axis.title.y = element_text(margin=unit(c(0,.3,0,0),"cm"))) +
  theme(panel.grid.major = element_line(size=0.1),panel.grid.minor = element_line(size=0.1)) +
  scale_colour_manual(values=plotPal) +
  theme(text = element_text(size = 10)) +
  scale_x_continuous(breaks=seq(0,1,by=.1), expand = c(0, 0), lim=c(0,1)) +
  scale_y_continuous(breaks=seq(0,1,by=.1), expand = c(0, 0), lim=c(0,1)) #+
#guides(colour=FALSE)
dev.off()


tvol.summer <- NULL
for (i in 1:nrow(species.table)) {
  data <- model.data.summer[[i]]
  colnames(data) <- c("uid","x","y","occ","tvol","tspd","coll")
  model <- coll.glm.summer[[i]]
  temp_df <- data.frame(x=data[,tvol], y=invcloglog(cbind(1,mean(log(data[,occ])),log(data[,tvol]),(log(data[,tvol]))*(log(data[,tvol])),mean(log(data[,tspd]))) %*% coef(model)), name=rep(paste(species.names[i]), each=length(data[,tvol])))
  tvol.summer <- rbind(tvol.summer,temp_df)
  rm(data)
  rm(model)
  rm(temp_df)
}
tvol.summer <- cbind(tvol.summer,season="summer")

tvol.autumn <- NULL
for (i in 1:nrow(species.table)) {
  data <- model.data.autumn[[i]]
  colnames(data) <- c("uid","x","y","occ","tvol","tspd","coll")
  model <- coll.glm.autumn[[i]]
  temp_df <- data.frame(x=data[,tvol], y=invcloglog(cbind(1,mean(log(data[,occ])),log(data[,tvol]),(log(data[,tvol]))*(log(data[,tvol])),mean(log(data[,tspd]))) %*% coef(model)), name=rep(paste(species.names[i]), each=length(data[,tvol])))
  tvol.autumn <- rbind(tvol.autumn,temp_df)
  rm(data)
  rm(model)
  rm(temp_df)
}
tvol.autumn <- cbind(tvol.autumn,season="autumn")

tvol.winter <- NULL
for (i in 1:nrow(species.table)) {
  data <- model.data.winter[[i]]
  colnames(data) <- c("uid","x","y","occ","tvol","tspd","coll")
  model <- coll.glm.winter[[i]]
  temp_df <- data.frame(x=data[,tvol], y=invcloglog(cbind(1,mean(log(data[,occ])),log(data[,tvol]),(log(data[,tvol]))*(log(data[,tvol])),mean(log(data[,tspd]))) %*% coef(model)), name=rep(paste(species.names[i]), each=length(data[,tvol])))
  tvol.winter <- rbind(tvol.winter,temp_df)
  rm(data)
  rm(model)
  rm(temp_df)
}
tvol.winter <- cbind(tvol.winter,season="winter")

tvol.spring <- NULL
for (i in 1:nrow(species.table)) {
  data <- model.data.spring[[i]]
  colnames(data) <- c("uid","x","y","occ","tvol","tspd","coll")
  model <- coll.glm.spring[[i]]
  temp_df <- data.frame(x=data[,tvol], y=invcloglog(cbind(1,mean(log(data[,occ])),log(data[,tvol]),(log(data[,tvol]))*(log(data[,tvol])),mean(log(data[,tspd]))) %*% coef(model)), name=rep(paste(species.names[i]), each=length(data[,tvol])))
  tvol.spring <- rbind(tvol.spring,temp_df)
  rm(data)
  rm(model)
  rm(temp_df)
}
tvol.spring <- cbind(tvol.spring,season="spring")

tvol.all <- rbind(tvol.summer,tvol.autumn,tvol.winter,tvol.spring)

tiff('figs/tvol_seasons.tif', pointsize = 12, compression = "lzw", res=300, width = 900, height = 900)
ggplot(tvol.all,aes(x=x/1000,y=y,group=interaction(name, season),colour=factor(name),shape=factor(season))) +
  geom_line(size=0.3,aes(linetype=season)) +
  ylab("Likelihood of Collision") +
  xlab("Traffic Volume (1000 vehicles/day)") +
  labs(color = "Species") +
  theme_bw() +
  theme(legend.key = element_blank(), legend.position="none") +
  theme(plot.margin=unit(c(.5,.5,.1,.1),"cm")) +
  theme(axis.title.x = element_text(margin=unit(c(.3,0,0,0),"cm"))) +
  theme(axis.title.y = element_text(margin=unit(c(0,.3,0,0),"cm"))) +
  theme(panel.grid.major = element_line(size=0.1),panel.grid.minor = element_line(size=0.1)) +
  scale_colour_manual(values=plotPal) +
  theme(text = element_text(size = 10)) +
  scale_x_continuous(breaks=seq(0,25,by=5), expand = c(0, 0), lim=c(0,25)) +
  scale_y_continuous(breaks=seq(0,1,by=.1), expand = c(0, 0), lim=c(0,1)) #+
#guides(colour=FALSE)
dev.off()


tspd.summer <- NULL
for (i in 1:nrow(species.table)) {
  data <- model.data.summer[[i]]
  colnames(data) <- c("uid","x","y","occ","tvol","tspd","coll")
  model <- coll.glm.summer[[i]]
  temp_df <- data.frame(x=data[,tspd], y=invcloglog(cbind(1,mean(log(data[,occ])),mean(log(data[,tvol])),mean((log(data[,tvol]))*(log(data[,tvol]))),log(data[,tspd])) %*% coef(model)), name=rep(paste(species.names[i]), each=length(data[,tspd])))
  tspd.summer <- rbind(tspd.summer,temp_df)
  rm(data)
  rm(model)
  rm(temp_df)
}
tspd.summer <- cbind(tspd.summer,season="summer")

tspd.autumn <- NULL
for (i in 1:nrow(species.table)) {
  data <- model.data.autumn[[i]]
  colnames(data) <- c("uid","x","y","occ","tvol","tspd","coll")
  model <- coll.glm.autumn[[i]]
  temp_df <- data.frame(x=data[,tspd], y=invcloglog(cbind(1,mean(log(data[,occ])),mean(log(data[,tvol])),mean((log(data[,tvol]))*(log(data[,tvol]))),log(data[,tspd])) %*% coef(model)), name=rep(paste(species.names[i]), each=length(data[,tspd])))
  tspd.autumn <- rbind(tspd.autumn,temp_df)
  rm(data)
  rm(model)
  rm(temp_df)
}
tspd.autumn <- cbind(tspd.autumn,season="autumn")

tspd.winter <- NULL
for (i in 1:nrow(species.table)) {
  data <- model.data.winter[[i]]
  colnames(data) <- c("uid","x","y","occ","tvol","tspd","coll")
  model <- coll.glm.winter[[i]]
  temp_df <- data.frame(x=data[,tspd], y=invcloglog(cbind(1,mean(log(data[,occ])),mean(log(data[,tvol])),mean((log(data[,tvol]))*(log(data[,tvol]))),log(data[,tspd])) %*% coef(model)), name=rep(paste(species.names[i]), each=length(data[,tspd])))
  tspd.winter <- rbind(tspd.winter,temp_df)
  rm(data)
  rm(model)
  rm(temp_df)
}
tspd.winter <- cbind(tspd.winter,season="winter")

tspd.spring <- NULL
for (i in 1:nrow(species.table)) {
  data <- model.data.spring[[i]]
  colnames(data) <- c("uid","x","y","occ","tvol","tspd","coll")
  model <- coll.glm.spring[[i]]
  temp_df <- data.frame(x=data[,tspd], y=invcloglog(cbind(1,mean(log(data[,occ])),mean(log(data[,tvol])),mean((log(data[,tvol]))*(log(data[,tvol]))),log(data[,tspd])) %*% coef(model)), name=rep(paste(species.names[i]), each=length(data[,tspd])))
  tspd.spring <- rbind(tspd.spring,temp_df)
  rm(data)
  rm(model)
  rm(temp_df)
}
tspd.spring <- cbind(tspd.spring,season="spring")

tspd.all <- rbind(tspd.summer,tspd.autumn,tspd.winter,tspd.spring)

tiff('figs/tspd_seasons.tif', pointsize = 12, compression = "lzw", res=300, width = 900, height = 900)
ggplot(tspd.all,aes(x=x,y=y,group=interaction(name, season),colour=factor(name),shape=factor(season))) +
  geom_line(size=0.3,aes(linetype=season)) +
  ylab("Likelihood of Collision") +
  xlab("Traffic Speed (km/hour)") +
  labs(color = "Species") +
  theme_bw() +
  theme(legend.key = element_blank(), legend.position="none") +
  theme(plot.margin=unit(c(.5,.5,.1,.1),"cm")) +
  theme(axis.title.x = element_text(margin=unit(c(.3,0,0,0),"cm"))) +
  theme(axis.title.y = element_text(margin=unit(c(0,.3,0,0),"cm"))) +
  theme(panel.grid.major = element_line(size=0.1),panel.grid.minor = element_line(size=0.1)) +
  scale_colour_manual(values=plotPal) +
  theme(text = element_text(size = 10)) +
  scale_x_continuous(breaks=seq(40,110,by=10), expand = c(0, 0), lim=c(40,110)) +
  scale_y_continuous(breaks=seq(0,1,by=.1), expand = c(0, 0), lim=c(0,1)) #+
#guides(colour=FALSE)
dev.off()

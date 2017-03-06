# require(maptools)
# require(raster)
# require(xtable)
require(ggplot2)
# require(ncf)
# require(doMC)
# require(data.table)
# require(plyr)

invcloglog <- function (x) {1-exp(-exp(x))}

plotPal <- c("#94d1c7", "#cccc2b", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#969696", "#bc80bd")

species.table <- read.delim("data/species_list.csv", header=T, sep=",")

species.names <- c("Eastern Grey Kangaroo","Common Brushtail Possum","Common Ringtail Possum","Black Swamp Wallaby","Common Wombat","Koala")

load("data/coll_model_data")
load("output/coll_glm")

occ <- NULL
for (i in 1:nrow(species.table)) {
  data <- model.data[[i]]
  colnames(data) <- c("uid","length","x","y","occ","tvol","tspd","coll")
  model <- coll.glm[[i]]
  temp_df <- data.frame(x=data[,occ], y=invcloglog(cbind(1,log(data[,occ]),mean(log(data[,tvol])),mean((log(data[,tvol]))*(log(data[,tvol]))),mean(log(data[,tspd]))) %*% coef(model)[1:5]), name=rep(paste(species.names[i]), each=length(data[,occ])))
  occ <- rbind(occ,temp_df)
  rm(data)
  rm(model)
  rm(temp_df)
}  

tiff('figs/occ.tif', pointsize = 12, compression = "lzw", res=300, width = 1440, height = 900)
ggplot(occ,aes(x=x,y=y,group=name,colour=factor(name))) +
  geom_line(size=0.3) +
  ylab("Relative Collision Risk") +
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
  scale_x_continuous(breaks=seq(0,1,by=.1), expand = c(0, 0), lim=c(0,1)) #+
  #scale_y_continuous(breaks=seq(0,1,by=.1), expand = c(0, 0), lim=c(0,1)) #+
  #guides(colour=FALSE)
dev.off()

# tiff('figs/occ_k.tif', pointsize = 12, compression = "lzw", res=300, width = 900, height = 900)
# ggplot(occ[occ$name==species.names[1],],aes(x=x,y=y)) +
#   geom_line(size=0.3, colour=plotPal[1]) +
#   ylab("Likelihood of Collision") +
#   xlab("Likelihood of Species Occurrence") +
#   theme_bw() +
#   theme(legend.key = element_blank()) +
#   theme(plot.margin=unit(c(.5,.5,.1,.1),"cm")) +
#   theme(axis.title.x = element_text(margin=unit(c(.3,0,0,0),"cm"))) +
#   theme(axis.title.y = element_text(margin=unit(c(0,.3,0,0),"cm"))) +
#   theme(panel.grid.major = element_line(size=0.1),panel.grid.minor = element_line(size=0.1)) +
#   theme(text = element_text(size = 10)) +
#   scale_x_continuous(breaks=seq(0,1,by=.1), expand = c(0, 0), lim=c(0,1)) +
#   scale_y_continuous(breaks=seq(0,1,by=.1), expand = c(0, 0), lim=c(0,1)) #+
# #guides(colour=FALSE)
# dev.off()

tvol <- NULL
for (i in 1:nrow(species.table)) {
  data <- model.data[[i]]
  colnames(data) <- c("uid","length","x","y","occ","tvol","tspd","coll")
  model <- coll.glm[[i]]
  temp_df <- data.frame(x=data[,tvol], y=invcloglog(cbind(1,mean(log(data[,occ])),log(data[,tvol]),(log(data[,tvol]))*(log(data[,tvol])),mean(log(data[,tspd]))) %*% coef(model)[1:5]), name=rep(paste(species.names[i]), each=length(data[,tvol])))
  tvol <- rbind(tvol,temp_df)
  rm(data)
  rm(model)
  rm(temp_df)
}

tiff('figs/tvol.tif', pointsize = 12, compression = "lzw", res=300, width = 900, height = 900)
ggplot(tvol,aes(x=x/1000,y=y,group=name,colour=factor(name))) +
  geom_line(size=0.3) +
  ylab("Relative Collision Risk") +
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
  scale_x_continuous(breaks=seq(0,25,by=5), expand = c(0, 0), lim=c(0,25)) #+
  #scale_y_continuous(breaks=seq(0,1,by=.1), expand = c(0, 0), lim=c(0,1)) #+
  #guides(colour=FALSE)
dev.off()

# tiff('figs/tvol_k.tif', pointsize = 12, compression = "lzw", res=300, width = 900, height = 900)
# ggplot(tvol[occ$name==species.names[1],],aes(x=x/1000,y=y)) +
#   geom_line(size=0.3, colour=plotPal[1]) +
#   ylab("Likelihood of Collision") +
#   xlab("Traffic Volume (1000 vehicles/day)") +
#   theme_bw() +
#   theme(legend.key = element_blank(), legend.position="none") +
#   theme(plot.margin=unit(c(.5,.5,.1,.1),"cm")) +
#   theme(axis.title.x = element_text(margin=unit(c(.3,0,0,0),"cm"))) +
#   theme(axis.title.y = element_text(margin=unit(c(0,.3,0,0),"cm"))) +
#   theme(panel.grid.major = element_line(size=0.1),panel.grid.minor = element_line(size=0.1)) +
#   theme(text = element_text(size = 10)) +
#   scale_x_continuous(breaks=seq(0,25,by=5), expand = c(0, 0), lim=c(0,25)) +
#   scale_y_continuous(breaks=seq(0,1,by=.1), expand = c(0, 0), lim=c(0,1)) #+
# #guides(colour=FALSE)
# dev.off()

tspd <- NULL
for (i in 1:nrow(species.table)) {
  data <- model.data[[i]]
  colnames(data) <- c("uid","length","x","y","occ","tvol","tspd","coll")
  model <- coll.glm[[i]]
  temp_df <- data.frame(x=data[,tspd], y=invcloglog(cbind(1,mean(log(data[,occ])),mean(log(data[,tvol])),mean((log(data[,tvol]))*(log(data[,tvol]))),log(data[,tspd])) %*% coef(model)[1:5]), name=rep(paste(species.names[i]), each=length(data[,tspd])))
  tspd <- rbind(tspd,temp_df)
  rm(data)
  rm(model)
  rm(temp_df)
}  

tiff('figs/tspd.tif', pointsize = 12, compression = "lzw", res=300, width = 900, height = 900)
ggplot(tspd,aes(x=x,y=y,group=name,colour=factor(name))) +
  geom_line(size=0.3) +
  ylab("Relative Collision Risk") +
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
  scale_x_continuous(breaks=seq(40,110,by=10), expand = c(0, 0), lim=c(40,110)) #+
  #scale_y_continuous(breaks=seq(0,1,by=.1), expand = c(0, 0), lim=c(0,1)) #+
  #guides(colour=FALSE)
dev.off()

# tiff('figs/tspd_k.tif', pointsize = 12, compression = "lzw", res=300, width = 900, height = 900)
# ggplot(tspd[occ$name==species.names[1],],aes(x=x,y=y)) +
#   geom_line(size=0.3, colour=plotPal[1]) +
#   ylab("Likelihood of Collision") +
#   xlab("Traffic Speed (km/hour)") +
#   theme_bw() +
#   theme(legend.key = element_blank(), legend.position="none") +
#   theme(plot.margin=unit(c(.5,.5,.1,.1),"cm")) +
#   theme(axis.title.x = element_text(margin=unit(c(.3,0,0,0),"cm"))) +
#   theme(axis.title.y = element_text(margin=unit(c(0,.3,0,0),"cm"))) +
#   theme(panel.grid.major = element_line(size=0.1),panel.grid.minor = element_line(size=0.1)) +
#   theme(text = element_text(size = 10)) +
#   scale_x_continuous(breaks=seq(40,110,by=10), expand = c(0, 0), lim=c(40,110)) +
#   scale_y_continuous(breaks=seq(0,1,by=.1), expand = c(0, 0), lim=c(0,1)) #+
# #guides(colour=FALSE)
# dev.off()

#calculate spatial autocorrelation across all species

#nb.list <- dnearneigh(as.matrix(data[,.(x,y)]), 0, 1000)
#nb.weights <- nb2listw(nb.list, zero.policy=TRUE)
registerDoMC(detectCores() - 1)

auto.resid <- foreach(i = 1:nrow(species.table), .packages = c("ncf"), .combine="rbind") %dopar% {
  data <- as.data.frame(model.data[[i]])
  model <- coll.glm[[i]]
  set.seed(123)
  sample.set <- sample(rownames(data[data$coll==0,]),2*nrow(data[data$coll==1,]))
  cor <- correlog(data[,2], data[,3], resid(model), increment=1000, resamp=0, latlon=FALSE)
 
  #cor <- correlog(data[,2], data[,3], resid(model), increment=1000, resamp=0, latlon=FALSE)
                                                                                                  
  temp_df <- data.frame(x=as.numeric(names(cor$correlation[1:20])), y=cor$correlation[1:20], name=rep(paste(species.names[i]), each=length(cor$correlation[1:20])))
  temp_df
}  

auto.coll <- foreach(i = 1:nrow(species.table), .packages = c("ncf"), .combine="rbind") %dopar% {
  data <- as.data.frame(model.data[[i]])
  cor <- correlog(data[,2], data[,3], data[,7], increment=1000, resamp=0, latlon=FALSE)
  temp_df <- data.frame(x=as.numeric(names(cor$correlation[1:20])), y=cor$correlation[1:20], name=rep(paste(species.names[i]), each=length(cor$correlation[1:20])))
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
  colnames(data) <- c("uid","length","x","y","occ","tvol","tspd","coll")
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
  colnames(data) <- c("uid","length","x","y","occ","tvol","tspd","coll")
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
  colnames(data) <- c("uid","length","x","y","occ","tvol","tspd","coll")
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
  colnames(data) <- c("uid","length","x","y","occ","tvol","tspd","coll")
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
  scale_x_continuous(breaks=seq(0,1,by=.1), expand = c(0, 0), lim=c(0,1)) #+
  #scale_y_continuous(breaks=seq(0,1,by=.1), expand = c(0, 0), lim=c(0,1)) #+
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

####################For Vic BioCon talk###############

tspd <- NULL
for (i in 1:nrow(species.table)) {
  data <- model.data[[i]]
  colnames(data) <- c("uid","x","y","occ","tvol","tspd","coll")
  model <- coll.glm[[i]]
  temp_df <- data.frame(x=seq(0.01,120,0.1), y=invcloglog(cbind(1,mean(log(data[,occ])),mean(log(data[,tvol])),mean((log(data[,tvol]))*(log(data[,tvol]))),log(seq(0.01,120,0.1))) %*% coef(model)[1:5]), name=rep(paste(species.names[i]), each=length(seq(0.01,120,0.1))))
  tspd <- rbind(tspd,temp_df)
  rm(data)
  rm(model)
  rm(temp_df)
} 

tspd <- NULL
for (i in 1:nrow(species.table)) {
  data <- model.data[[i]]
  colnames(data) <- c("uid","x","y","occ","tvol","tspd","coll")
  data2 <- data.frame(occ=mean(data$occ),tvol=mean(data$tvol),tspd=seq(0.01,120,0.1))
  colnames(data2) <- c(paste0(species.table[i,2]),"tvol","tspd")
  model <- coll.glm[[i]]
  tspd.fit <- predict.glm(model,data2,type="response",se.fit=TRUE)
  temp_df <- data.frame(x=seq(0.01,120,0.1),y=tspd.fit[["fit"]],ymin=tspd.fit[["fit"]]-1.96*tspd.fit[["se.fit"]],ymax=tspd.fit[["fit"]]+1.96*tspd.fit[["se.fit"]], name=rep(paste(species.names[i]), each=length(seq(0.01,120,0.1))))
  temp_df$y <- (temp_df$y/(data[coll==1,.N]/nrow(data)))/nrow(data)
  #temp_df$y <- (temp_df$y - min(temp_df$y)) / (max(temp_df$y) - min(temp_df$y))
  tspd <- rbind(tspd,temp_df)
  rm(data)
  rm(data2)
  rm(model)
  rm(temp_df)
} 

pdf('/home/casey/Research/Projects/VicBioConf/graphics/tspd.pdf', pointsize = 16)
ggplot(tspd,aes(x=x,y=y,group=name,colour=factor(name))) +
  geom_line(size=0.3) +
  ylab("Relative Collision Rate") +
  xlab("Traffic Speed (km/hour)") +
  labs(color = "Species") +
  theme_bw() +
  theme(legend.key = element_blank(), legend.position="none") +
  theme(plot.margin=unit(c(.5,.5,.1,.1),"cm")) +
  theme(axis.title.x = element_text(margin=unit(c(.3,0,0,0),"cm"))) +
  theme(axis.title.y = element_text(margin=unit(c(0,.3,0,0),"cm"))) +
  theme(panel.grid.major = element_line(size=0.1),panel.grid.minor = element_line(size=0.1)) +
  scale_colour_manual(values=plotPal) +
  theme(text = element_text(size = 16)) +
  scale_x_continuous(breaks=seq(0,120,by=20), expand = c(0, 0), lim=c(0,120)) #+
  #scale_y_continuous(breaks=seq(0,1,by=.1), expand = c(0, 0), lim=c(0,1)) #+
#guides(colour=FALSE)
dev.off()

###### alternative based on predictions for phd talk######

occ <- NULL
for (i in 1:nrow(species.table)) {
  data <- as.data.frame(model.data[[i]])
  occ.range <- seq(0,1,1/(nrow(data)+1))[-c(1,length(seq(0,1,1/(nrow(data)+1))))]
  model <- coll.glm[[i]]
  data2 <- data.frame(occ=occ.range,tvol=mean(data$tvol),tspd=mean(data$tspd),length=1)
  colnames(data2)[1] <- paste0(species.table[i,2])
  occ.fit <- predict.glm(model,data2,type="response",se.fit=TRUE)
  temp_df <- data.frame(x=occ.range,y=occ.fit[["fit"]],ymin=occ.fit[["fit"]]-1.96*occ.fit[["se.fit"]],ymax=occ.fit[["fit"]]+1.96*occ.fit[["se.fit"]],name=rep(paste(species.names[i]), each=nrow(data)))
  occ <- rbind(occ,temp_df)
  rm(data)
  rm(occ.range)
  rm(model)
  rm(data2)
  rm(occ.fit)
  rm(temp_df)
}

occ_k <- subset(occ,name=="Eastern Grey Kangaroo")

#pdf('/home/casey/Research/Projects/PhD_Thesis/completion_talk/graphics/6sp_occ_k.pdf', pointsize = 20)
png('/home/casey/Research/Projects/PhD_Thesis/completion_talk/graphics/6sp_occ_k.png', pointsize = 16, res=300, width = 1000, height = 900)
ggplot(occ_k,aes(x=x,y=y,ymin=ymin,ymax=ymax,group=name,colour=name)) +
  geom_line(size=0.8) +
  geom_ribbon(alpha=0.3, colour=NA) +
  ylab("RELATIVE COLLISION RATE") +
  xlab("LIKELIHOOD OF SPECIES OCCURRENCE") +
  theme_bw() +
  theme(legend.key = element_blank()) +
  theme(legend.position="none") +
  theme(plot.margin=unit(c(.5,.5,.1,.1),"cm")) +
  theme(axis.title.x = element_text(margin=unit(c(.3,0,0,0),"cm"))) +
  theme(axis.title.y = element_text(margin=unit(c(0,.3,0,0),"cm"))) +
  theme(panel.grid.major = element_line(size=0.1),panel.grid.minor = element_line(size=0.1)) +
  scale_colour_manual(values=plotPal) +
  theme(text = element_text(size = 8)) +
  scale_x_continuous(breaks=seq(0,1,by=.1), expand = c(0, 0), lim=c(0,1))
dev.off()

#pdf('/home/casey/Research/Projects/PhD_Thesis/completion_talk/graphics/6sp_occ.pdf', pointsize = 20)
png('/home/casey/Research/Projects/PhD_Thesis/completion_talk/graphics/6sp_occ.png', pointsize = 16, res=300, width = 1000, height = 900)
ggplot(occ,aes(x=x,y=y,group=name,colour=name)) +
  geom_line(size=0.8) +
  ylab("RELATIVE COLLISION RATE") +
  xlab("LIKELIHOOD OF SPECIES OCCURRENCE") +
  theme_bw() +
  theme(legend.key = element_blank()) +
  theme(legend.position="none") +
  theme(plot.margin=unit(c(.5,.5,.1,.1),"cm")) +
  theme(axis.title.x = element_text(margin=unit(c(.3,0,0,0),"cm"))) +
  theme(axis.title.y = element_text(margin=unit(c(0,.3,0,0),"cm"))) +
  theme(panel.grid.major = element_line(size=0.1),panel.grid.minor = element_line(size=0.1)) +
  scale_colour_manual(values=plotPal) +
  theme(text = element_text(size = 8)) +
  scale_x_continuous(breaks=seq(0,1,by=.1), expand = c(0, 0), lim=c(0,1))
dev.off()

lscale <- rep(2,nrow(species.table))

for (i in 1:nrow(species.table)) {
  plotPal.mod <- plotPal
  plotPal.mod[-i] <- "#c0c0c0"
  lscale.mod <- lscale
  lscale.mod[-i] <- 0.5
  occ.mod <- ddply(occ, "name",transform, y=(y-min(y))/(max(y)-min(y)))
  png(paste0('/home/casey/Research/Projects/PhD_Thesis/completion_talk/graphics/6sp_occ_',species.table[i,2],'.png'), pointsize = 16, res=300, width = 1000, height = 900)
  #pdf(paste0('/home/casey/Research/Projects/PhD_Thesis/completion_talk/graphics/6sp_occ_',species.table[i,2],'.pdf'), pointsize = 20)
  print(
  ggplot(occ.mod,aes(x=x,y=y,group=name,colour=name,size=name)) +
    geom_line() +
    ylab("") +
    xlab("") +
    theme_bw() +
    theme(legend.key = element_blank()) +
    theme(legend.position="none") +
    theme(plot.margin=unit(c(.1,.1,.1,.1),"cm")) +
    theme(axis.title.x = element_blank()) +
    theme(axis.title.y =element_blank()) +
    theme(panel.grid.major = element_line(size=0.1),panel.grid.minor = element_line(size=0.1)) +
    scale_colour_manual(values=plotPal.mod) +
    scale_size_manual(values=lscale.mod) + 
    scale_y_continuous(breaks=NULL, expand = c(0, 0), lim=c(0,1)) +
    scale_x_continuous(breaks=NULL, expand = c(0, 0), lim=c(0,1))
  )
  dev.off()
  rm(plotPal.mod)
  rm(lscale.mod)
  rm(occ.mod)
}


tvol <- NULL
for (i in 1:nrow(species.table)) {
  data <- as.data.frame(model.data[[i]])
  tvol.range <- seq(0,40000,40000/(nrow(data)+1))[-c(1,length(seq(0,40000,40000/(nrow(data)+1))))]
  model <- coll.glm[[i]]
  data2 <- data.frame(occ=mean(data[,5]),tvol=tvol.range,tspd=mean(data$tspd),length=1)
  colnames(data2)[1] <- paste0(species.table[i,2])
  tvol.fit <- predict.glm(model,data2,type="response",se.fit=TRUE)
  temp_df <- data.frame(x=tvol.range,y=tvol.fit[["fit"]],ymin=tvol.fit[["fit"]]-1.96*tvol.fit[["se.fit"]],ymax=tvol.fit[["fit"]]+1.96*tvol.fit[["se.fit"]],name=rep(paste(species.names[i]), each=nrow(data)))
  tvol <- rbind(tvol,temp_df)
  rm(data)
  rm(tvol.range)
  rm(model)
  rm(data2)
  rm(tvol.fit)
  rm(temp_df)
}

tvol_k <- subset(tvol,name=="Eastern Grey Kangaroo")

#pdf('/home/casey/Research/Projects/PhD_Thesis/completion_talk/graphics/6sp_tvol_k.pdf', pointsize = 20)
png('/home/casey/Research/Projects/PhD_Thesis/completion_talk/graphics/6sp_tvol_k.png', pointsize = 16, res=300, width = 1000, height = 900)
ggplot(tvol_k,aes(x=x/1000,y=y,ymin=ymin,ymax=ymax,group=name,colour=name)) +
  geom_line(size=0.8) +
  geom_ribbon(alpha=0.3, colour=NA) +
  ylab("RELATIVE COLLISION RATE") +
  xlab("TRAFFIC VOLUME (1000 VEHICLES/DAY)") +
  theme_bw() +
  theme(legend.key = element_blank()) +
  theme(legend.position="none") +
  theme(plot.margin=unit(c(.5,.5,.1,.1),"cm")) +
  theme(axis.title.x = element_text(margin=unit(c(.3,0,0,0),"cm"))) +
  theme(axis.title.y = element_text(margin=unit(c(0,.3,0,0),"cm"))) +
  theme(panel.grid.major = element_line(size=0.1),panel.grid.minor = element_line(size=0.1)) +
  scale_colour_manual(values=plotPal) +
  theme(text = element_text(size = 8)) +
  scale_x_continuous(breaks=seq(0,40,by=5), expand = c(0, 0), lim=c(0,40))
dev.off()

#pdf('/home/casey/Research/Projects/PhD_Thesis/completion_talk/graphics/6sp_tvol.pdf', pointsize = 20)
png('/home/casey/Research/Projects/PhD_Thesis/completion_talk/graphics/6sp_tvol.png', pointsize = 16, res=300, width = 1000, height = 900)
ggplot(tvol,aes(x=x/1000,y=y,group=name,colour=name)) +
  geom_line(size=0.8) +
  ylab("RELATIVE COLLISION RATE") +
  xlab("TRAFFIC VOLUME (1000 VEHICLES/DAY)") +
  theme_bw() +
  theme(legend.key = element_blank()) +
  theme(legend.position="none") +
  theme(plot.margin=unit(c(.5,.5,.1,.1),"cm")) +
  theme(axis.title.x = element_text(margin=unit(c(.3,0,0,0),"cm"))) +
  theme(axis.title.y = element_text(margin=unit(c(0,.3,0,0),"cm"))) +
  theme(panel.grid.major = element_line(size=0.1),panel.grid.minor = element_line(size=0.1)) +
  scale_colour_manual(values=plotPal) +
  theme(text = element_text(size = 8)) +
  scale_x_continuous(breaks=seq(0,40,by=5), expand = c(0, 0), lim=c(0,40))
dev.off()

lscale <- rep(2,nrow(species.table))

for (i in 1:nrow(species.table)) {
  plotPal.mod <- plotPal
  plotPal.mod[-i] <- "#c0c0c0"
  lscale.mod <- lscale
  lscale.mod[-i] <- 0.5
  tvol.mod <- ddply(tvol, "name",transform, y=(y-min(y))/(max(y)-min(y)))
  png(paste0('/home/casey/Research/Projects/PhD_Thesis/completion_talk/graphics/6sp_tvol_',species.table[i,2],'.png'), pointsize = 16, res=300, width = 1000, height = 900)
  #pdf(paste0('/home/casey/Research/Projects/PhD_Thesis/completion_talk/graphics/6sp_tvol_',species.table[i,2],'.pdf'), pointsize = 20)
  print(
  ggplot(tvol.mod,aes(x=x/1000,y=y,group=name,colour=name,size=name)) +
    geom_line() +
    ylab("") +
    xlab("") +
    theme_bw() +
    theme(legend.key = element_blank()) +
    theme(legend.position="none") +
    theme(plot.margin=unit(c(.1,.1,.1,.1),"cm")) +
    theme(axis.title.x = element_blank()) +
    theme(axis.title.y =element_blank()) +
    theme(panel.grid.major = element_line(size=0.1),panel.grid.minor = element_line(size=0.1)) +
    scale_colour_manual(values=plotPal.mod) +
    scale_size_manual(values=lscale.mod) + 
    scale_y_continuous(breaks=NULL, expand = c(0, 0), lim=c(0,1)) +
    scale_x_continuous(breaks=NULL, expand = c(0, 0), lim=c(0,40))
  )
  dev.off()
  rm(plotPal.mod)
  rm(lscale.mod)
  rm(tvol.mod)
}


tspd <- NULL
for (i in 1:nrow(species.table)) {
  data <- as.data.frame(model.data[[i]])
  tspd.range <- tspd.range <- seq(0,110,110/(nrow(data)+1))[-c(1,length(seq(0,110,110/(nrow(data)+1))))]
  model <- coll.glm[[i]]
  data2 <- data.frame(occ=mean(data[,5]),tvol=mean(data$tvol),tspd=tspd.range,length=1)
  colnames(data2)[1] <- paste0(species.table[i,2])
  tspd.fit <- predict.glm(model,data2,type="response",se.fit=TRUE)
  temp_df <- data.frame(x=tspd.range,y=tspd.fit[["fit"]],ymin=tspd.fit[["fit"]]-1.96*tspd.fit[["se.fit"]],ymax=tspd.fit[["fit"]]+1.96*tspd.fit[["se.fit"]],name=rep(paste(species.names[i]), each=nrow(data)))
  tspd <- rbind(tspd,temp_df)
  rm(data)
  rm(tspd.range)
  rm(model)
  rm(data2)
  rm(tspd.fit)
  rm(temp_df)
}

tspd_k <- subset(tspd,name=="Eastern Grey Kangaroo")

#pdf('/home/casey/Research/Projects/PhD_Thesis/completion_talk/graphics/6sp_tspd_k.pdf', pointsize = 20)
png('/home/casey/Research/Projects/PhD_Thesis/completion_talk/graphics/6sp_tspd_k.png', pointsize = 16, res=300, width = 1000, height = 900)
ggplot(tspd_k,aes(x=x,y=y,ymin=ymin,ymax=ymax,group=name,colour=name)) +
  geom_line(size=0.8) +
  geom_ribbon(alpha=0.3, colour=NA) +
  ylab("RELATIVE COLLISION RATE") +
  xlab("TRAFFIC SPEED (KM/HOUR)") +
  theme_bw() +
  theme(legend.key = element_blank()) +
  theme(legend.position="none") +
  theme(plot.margin=unit(c(.5,.5,.1,.1),"cm")) +
  theme(axis.title.x = element_text(margin=unit(c(.3,0,0,0),"cm"))) +
  theme(axis.title.y = element_text(margin=unit(c(0,.3,0,0),"cm"))) +
  theme(panel.grid.major = element_line(size=0.1),panel.grid.minor = element_line(size=0.1)) +
  scale_colour_manual(values=plotPal) +
  theme(text = element_text(size = 8)) +
  scale_x_continuous(breaks=seq(0,110,by=10), expand = c(0, 0), lim=c(0,110))
dev.off()

#pdf('/home/casey/Research/Projects/PhD_Thesis/completion_talk/graphics/6sp_tspd.pdf', pointsize = 20)
png('/home/casey/Research/Projects/PhD_Thesis/completion_talk/graphics/6sp_tspd.png', pointsize = 16, res=300, width = 1000, height = 900)
ggplot(tspd,aes(x=x,y=y,group=name,colour=name)) +
  geom_line(size=0.8) +
  ylab("RELATIVE COLLISION RATE") +
  xlab("TRAFFIC SPEED (KM/HOUR)") +
  theme_bw() +
  theme(legend.key = element_blank()) +
  theme(legend.position="none") +
  theme(plot.margin=unit(c(.5,.5,.1,.1),"cm")) +
  theme(axis.title.x = element_text(margin=unit(c(.3,0,0,0),"cm"))) +
  theme(axis.title.y = element_text(margin=unit(c(0,.3,0,0),"cm"))) +
  theme(panel.grid.major = element_line(size=0.1),panel.grid.minor = element_line(size=0.1)) +
  scale_colour_manual(values=plotPal) +
  theme(text = element_text(size = 8)) +
  scale_x_continuous(breaks=seq(0,110,by=10), expand = c(0, 0), lim=c(0,110))
dev.off()

lscale <- rep(2,nrow(species.table))

for (i in 1:nrow(species.table)) {
  plotPal.mod <- plotPal
  plotPal.mod[-i] <- "#c0c0c0"
  lscale.mod <- lscale
  lscale.mod[-i] <- 0.5
  tspd.mod <- ddply(tspd, "name",transform, y=(y-min(y))/(max(y)-min(y)))
  png(paste0('/home/casey/Research/Projects/PhD_Thesis/completion_talk/graphics/6sp_tspd_',species.table[i,2],'.png'), pointsize = 16, res=300, width = 1000, height = 900)
  #pdf(paste0('/home/casey/Research/Projects/PhD_Thesis/completion_talk/graphics/6sp_tspd_',species.table[i,2],'.pdf'), pointsize = 20)
  print(
  ggplot(tspd.mod,aes(x=x,y=y,group=name,colour=name,size=name)) +
    geom_line() +
    ylab("") +
    xlab("") +
    theme_bw() +
    theme(legend.key = element_blank()) +
    theme(legend.position="none") +
    theme(plot.margin=unit(c(.1,.1,.1,.1),"cm")) +
    theme(axis.title.x = element_blank()) +
    theme(axis.title.y =element_blank()) +
    theme(panel.grid.major = element_line(size=0.1),panel.grid.minor = element_line(size=0.1)) +
    scale_colour_manual(values=plotPal.mod) +
    scale_size_manual(values=lscale.mod) + 
    scale_y_continuous(breaks=NULL, expand = c(0, 0), lim=c(0,1)) +
    scale_x_continuous(breaks=NULL, expand = c(0, 0), lim=c(0,110))
  )
  dev.off()
  rm(plotPal.mod)
  rm(lscale.mod)
  rm(tspd.mod)
}
#############################################
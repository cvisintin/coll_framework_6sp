# require(maptools)
# require(raster)
# require(xtable)
require(ggplot2)
# require(ncf)
# require(doMC)
require(data.table)
require(plyr)

invcloglog <- function (x) {1-exp(-exp(x))}

plotPal <- c("#94d1c7", "#cccc2b", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#969696", "#bc80bd")

species.table <- read.delim("data/species_list.csv", header=T, sep=",", stringsAsFactors = F)

species.names <- c("Eastern Grey Kangaroo","Common Brushtail Possum","Common Ringtail Possum","Black Swamp Wallaby","Common Wombat","Koala")

#load("data/coll_model_data")
load("output/coll_glm")

occ <- NULL
for (i in 1:nrow(species.table)) {
  model <- coll.glm[[i]]
  occ.range <- seq(0,1,1/(nrow(model$data)+1))[-c(1,length(seq(0,1,1/(nrow(model$data)+1))))]
  data <- data.frame(occ=occ.range,tvol=mean(model$data$tvol),tspd=mean(model$data$tspd),length=1)
  colnames(data)[1] <- paste0(species.table[i,2])
  occ.fit <- predict.glm(model,data,type="response",se.fit=TRUE)
  temp_df <- data.frame(x=occ.range,y=occ.fit[["fit"]],ymin=occ.fit[["fit"]]-1.96*occ.fit[["se.fit"]],ymax=occ.fit[["fit"]]+1.96*occ.fit[["se.fit"]],name=rep(paste(species.names[i]), each=nrow(data)))
  occ <- rbind(occ,temp_df)
  rm(data)
  rm(occ.range)
  rm(model)
  rm(occ.fit)
  rm(temp_df)
}

tiff('figs/occ_a.tif', pointsize = 16, compression = "lzw", res=300, width = 900, height = 900)
ggplot(occ,aes(x=x,y=y,group=name,colour=factor(name))) +
  geom_line(size=0.5) +
  ylab("Relative Collision Rate") +
  xlab("Likelihood of Species Occurrence") +
  labs(color = "Species") +
  theme_bw() +
  theme(legend.key = element_blank(), legend.position="none") +
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
  tiff(paste0('figs/occ_a_',species.table[i,2],'.tif'), pointsize = 16, res=300, width = 900, height = 900)
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
  model <- coll.glm[[i]]
  tvol.range <- seq(0,40000,40000/(nrow(model$data)+1))[-c(1,length(seq(0,40000,40000/(nrow(model$data)+1))))]
  data <- data.frame(occ=mean(model$data[[5]]),tvol=tvol.range,tspd=mean(model$data$tspd),length=1)
  colnames(data)[1] <- paste0(species.table[i,2])
  tvol.fit <- predict.glm(model,data,type="response",se.fit=TRUE)
  temp_df <- data.frame(x=tvol.range,y=tvol.fit[["fit"]],ymin=tvol.fit[["fit"]]-1.96*tvol.fit[["se.fit"]],ymax=tvol.fit[["fit"]]+1.96*tvol.fit[["se.fit"]],name=rep(paste(species.names[i]), each=nrow(data)))
  tvol <- rbind(tvol,temp_df)
  rm(data)
  rm(tvol.range)
  rm(model)
  rm(tvol.fit)
  rm(temp_df)
}

tiff('figs/tvol_b.tif', pointsize = 16, compression = "lzw", res=300, width = 900, height = 900)
ggplot(tvol,aes(x=x/1000,y=y,group=name,colour=factor(name))) +
  geom_line(size=0.5) +
  ylab("Relative Collision Rate") +
  xlab("Traffic Volume (1000 vehicles/day)") +
  labs(color = "Species") +
  theme_bw() +
  theme(legend.key = element_blank(), legend.position="none") +
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
  tiff(paste0('figs/tvol_b_',species.table[i,2],'.tif'), pointsize = 16, res=300, width = 900, height = 900)
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
  model <- coll.glm[[i]]
  tspd.range <- tspd.range <- seq(0,110,110/(nrow(model$data)+1))[-c(1,length(seq(0,110,110/(nrow(model$data)+1))))]
  data <- data.frame(occ=mean(model$data[[5]]),tvol=mean(model$data$tvol),tspd=tspd.range,length=1)
  colnames(data)[1] <- paste0(species.table[i,2])
  tspd.fit <- predict.glm(model,data,type="response",se.fit=TRUE)
  temp_df <- data.frame(x=tspd.range,y=tspd.fit[["fit"]],ymin=tspd.fit[["fit"]]-1.96*tspd.fit[["se.fit"]],ymax=tspd.fit[["fit"]]+1.96*tspd.fit[["se.fit"]],name=rep(paste(species.names[i]), each=nrow(data)))
  tspd <- rbind(tspd,temp_df)
  rm(data)
  rm(tspd.range)
  rm(model)
  rm(tspd.fit)
  rm(temp_df)
} 

tiff('figs/tspd_c.tif', pointsize = 16, compression = "lzw", res=300, width = 900, height = 900)
ggplot(tspd,aes(x=x,y=y,group=name,colour=factor(name))) +
  geom_line(size=0.5) +
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
  tiff(paste0('figs/tspd_c_',species.table[i,2],'.tif'), pointsize = 16, res=300, width = 900, height = 900)
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

load("output/coll_glm_summer")
occ.summer <- NULL
for (i in 1:nrow(species.table)) {
  model <- coll.glm.summer[[i]]
  data <- data.frame("occ"=model$data[[5]],"tvol"=mean(model$data$tvol),"tspd"=mean(model$data$tspd),"length"=1)
  colnames(data)[1] <- species.table[i,2]
  temp_df <- data.frame(x=data[,1],
                        y=predict(model, data, type="response"),
                        name=rep(paste(species.names[i]), each=length(data[,1]))
                        )
  occ.summer <- rbind(occ.summer,temp_df)
  rm(data)
  rm(model)
  rm(temp_df)
}
rm(coll.glm.summer)
occ.summer <- cbind(occ.summer,season="summer")

load("output/coll_glm_autumn")
occ.autumn <- NULL
for (i in 1:nrow(species.table)) {
  model <- coll.glm.autumn[[i]]
  data <- data.frame("occ"=model$data[[5]],"tvol"=mean(model$data$tvol),"tspd"=mean(model$data$tspd),"length"=1)
  colnames(data)[1] <- species.table[i,2]
  temp_df <- data.frame(x=data[,1],
                        y=predict(model, data, type="response"),
                        name=rep(paste(species.names[i]), each=length(data[,1]))
  )
  occ.autumn <- rbind(occ.autumn,temp_df)
  rm(data)
  rm(model)
  rm(temp_df)
}
rm(coll.glm.autumn)
occ.autumn <- cbind(occ.autumn,season="autumn")

load("output/coll_glm_winter")
occ.winter <- NULL
for (i in 1:nrow(species.table)) {
  model <- coll.glm.winter[[i]]
  data <- data.frame("occ"=model$data[[5]],"tvol"=mean(model$data$tvol),"tspd"=mean(model$data$tspd),"length"=1)
  colnames(data)[1] <- species.table[i,2]
  temp_df <- data.frame(x=data[,1],
                        y=predict(model, data, type="response"),
                        name=rep(paste(species.names[i]), each=length(data[,1]))
  )
  occ.winter <- rbind(occ.winter,temp_df)
  rm(data)
  rm(model)
  rm(temp_df)
}
rm(coll.glm.winter)
occ.winter <- cbind(occ.winter,season="winter")

load("output/coll_glm_spring")
occ.spring <- NULL
for (i in 1:nrow(species.table)) {
  model <- coll.glm.spring[[i]]
  data <- data.frame("occ"=model$data[[5]],"tvol"=mean(model$data$tvol),"tspd"=mean(model$data$tspd),"length"=1)
  colnames(data)[1] <- species.table[i,2]
  temp_df <- data.frame(x=data[,1],
                        y=predict(model, data, type="response"),
                        name=rep(paste(species.names[i]), each=length(data[,1]))
  )
  occ.spring <- rbind(occ.spring,temp_df)
  rm(data)
  rm(model)
  rm(temp_df)
}
rm(coll.glm.spring)
occ.spring <- cbind(occ.spring,season="spring")

occ.all <- rbind(occ.summer,occ.autumn,occ.winter,occ.spring)

for(i in species.names) {
  tiff(paste0('figs/occ_seasons_',c('a','b','c','d','e','f')[which(species.names==i)],'.tif'), pointsize = 12, compression = "lzw", res=300, width = 1100, height = 900)
  print(
    ggplot(occ.all[occ.all$name==i,],aes(x=x,y=y,group=season,shape=factor(season))) +
      geom_line(size=0.3,aes(linetype=season)) +
      ylab("Relative Collision Rate") +
      xlab("Likelihood of Species Occurrence") +
      labs(linetype = "Season") +
      labs(title = i) +
      theme_bw() +
      theme(legend.key = element_blank(),legend.key.size = unit(1.0, 'lines')) +
      theme(plot.margin=unit(c(.5,0,.1,.1),"cm")) +
      theme(axis.title.x = element_text(margin=unit(c(.3,0,0,0),"cm"))) +
      theme(axis.title.y = element_text(margin=unit(c(0,.3,0,0),"cm"))) +
      theme(panel.grid.major = element_line(size=0.1),panel.grid.minor = element_line(size=0.1)) +
      theme(text = element_text(size = 10)) +
      scale_x_continuous(breaks=seq(0,1,by=.1), expand = c(0, 0), lim=c(0,1))
  )
  dev.off()
}


load("output/coll_glm_summer")
tvol.summer <- NULL
for (i in 1:nrow(species.table)) {
  model <- coll.glm.summer[[i]]
  data <- data.frame("occ"=mean(model$data[[5]]),"tvol"=model$data$tvol,"tspd"=mean(model$data$tspd),"length"=1)
  colnames(data)[1] <- species.table[i,2]
  temp_df <- data.frame(x=data$tvol,
                        y=predict(model, data, type="response"),
                        name=rep(paste(species.names[i]), each=length(data$tvol))
  )
  tvol.summer <- rbind(tvol.summer,temp_df)
  rm(data)
  rm(model)
  rm(temp_df)
}
rm(coll.glm.summer)
tvol.summer <- cbind(tvol.summer,season="summer")

load("output/coll_glm_autumn")
tvol.autumn <- NULL
for (i in 1:nrow(species.table)) {
  model <- coll.glm.autumn[[i]]
  data <- data.frame("occ"=mean(model$data[[5]]),"tvol"=model$data$tvol,"tspd"=mean(model$data$tspd),"length"=1)
  colnames(data)[1] <- species.table[i,2]
  temp_df <- data.frame(x=data$tvol,
                        y=predict(model, data, type="response"),
                        name=rep(paste(species.names[i]), each=length(data$tvol))
  )
  tvol.autumn <- rbind(tvol.autumn,temp_df)
  rm(data)
  rm(model)
  rm(temp_df)
}
rm(coll.glm.autumn)
tvol.autumn <- cbind(tvol.autumn,season="autumn")

load("output/coll_glm_winter")
tvol.winter <- NULL
for (i in 1:nrow(species.table)) {
  model <- coll.glm.winter[[i]]
  data <- data.frame("occ"=mean(model$data[[5]]),"tvol"=model$data$tvol,"tspd"=mean(model$data$tspd),"length"=1)
  colnames(data)[1] <- species.table[i,2]
  temp_df <- data.frame(x=data$tvol,
                        y=predict(model, data, type="response"),
                        name=rep(paste(species.names[i]), each=length(data$tvol))
  )
  tvol.winter <- rbind(tvol.winter,temp_df)
  rm(data)
  rm(model)
  rm(temp_df)
}
rm(coll.glm.winter)
tvol.winter <- cbind(tvol.winter,season="winter")

load("output/coll_glm_spring")
tvol.spring <- NULL
for (i in 1:nrow(species.table)) {
  model <- coll.glm.spring[[i]]
  data <- data.frame("occ"=mean(model$data[[5]]),"tvol"=model$data$tvol,"tspd"=mean(model$data$tspd),"length"=1)
  colnames(data)[1] <- species.table[i,2]
  temp_df <- data.frame(x=data$tvol,
                        y=predict(model, data, type="response"),
                        name=rep(paste(species.names[i]), each=length(data$tvol))
  )
  tvol.spring <- rbind(tvol.spring,temp_df)
  rm(data)
  rm(model)
  rm(temp_df)
}
rm(coll.glm.spring)
tvol.spring <- cbind(tvol.spring,season="spring")

tvol.all <- rbind(tvol.summer,tvol.autumn,tvol.winter,tvol.spring)

for(i in species.names) {
  tiff(paste0('figs/tvol_seasons_',c('a','b','c','d','e','f')[which(species.names==i)],'.tif'), pointsize = 12, compression = "lzw", res=300, width = 1100, height = 900)
  print(
    ggplot(tvol.all[tvol.all$name==i,],aes(x=x/1000,y=y,group=season,shape=factor(season))) +
      geom_line(size=0.3,aes(linetype=season)) +
      ylab("Relative Collision Rate") +
      xlab("Traffic Volume (1000 vehicles/day)") +
      labs(linetype = "Season") +
      labs(title = i) +
      theme_bw() +
      theme(legend.key = element_blank(),legend.key.size = unit(1.0, 'lines')) +
      theme(plot.margin=unit(c(.5,0,.1,.1),"cm")) +
      theme(axis.title.x = element_text(margin=unit(c(.3,0,0,0),"cm"))) +
      theme(axis.title.y = element_text(margin=unit(c(0,.3,0,0),"cm"))) +
      theme(panel.grid.major = element_line(size=0.1),panel.grid.minor = element_line(size=0.1)) +
      theme(text = element_text(size = 10)) +
      scale_x_continuous(breaks=seq(0,25,by=5), expand = c(0, 0), lim=c(0,25))
  )
  dev.off()
}


load("output/coll_glm_summer")
tspd.summer <- NULL
for (i in 1:nrow(species.table)) {
  model <- coll.glm.summer[[i]]
  data <- data.frame("occ"=mean(model$data[[5]]),"tvol"=mean(model$data$tvol),"tspd"=model$data$tspd,"length"=1)
  colnames(data)[1] <- species.table[i,2]
  temp_df <- data.frame(x=data$tspd,
                        y=predict(model, data, type="response"),
                        name=rep(paste(species.names[i]), each=length(data$tspd))
  )
  tspd.summer <- rbind(tspd.summer,temp_df)
  rm(data)
  rm(model)
  rm(temp_df)
}
rm(coll.glm.summer)
tspd.summer <- cbind(tspd.summer,season="summer")

load("output/coll_glm_autumn")
tspd.autumn <- NULL
for (i in 1:nrow(species.table)) {
  model <- coll.glm.autumn[[i]]
  data <- data.frame("occ"=mean(model$data[[5]]),"tvol"=mean(model$data$tvol),"tspd"=model$data$tspd,"length"=1)
  colnames(data)[1] <- species.table[i,2]
  temp_df <- data.frame(x=data$tspd,
                        y=predict(model, data, type="response"),
                        name=rep(paste(species.names[i]), each=length(data$tspd))
  )
  tspd.autumn <- rbind(tspd.autumn,temp_df)
  rm(data)
  rm(model)
  rm(temp_df)
}
rm(coll.glm.autumn)
tspd.autumn <- cbind(tspd.autumn,season="autumn")

load("output/coll_glm_winter")
tspd.winter <- NULL
for (i in 1:nrow(species.table)) {
  model <- coll.glm.winter[[i]]
  data <- data.frame("occ"=mean(model$data[[5]]),"tvol"=mean(model$data$tvol),"tspd"=model$data$tspd,"length"=1)
  colnames(data)[1] <- species.table[i,2]
  temp_df <- data.frame(x=data$tspd,
                        y=predict(model, data, type="response"),
                        name=rep(paste(species.names[i]), each=length(data$tspd))
  )
  tspd.winter <- rbind(tspd.winter,temp_df)
  rm(data)
  rm(model)
  rm(temp_df)
}
rm(coll.glm.winter)
tspd.winter <- cbind(tspd.winter,season="winter")

load("output/coll_glm_spring")
tspd.spring <- NULL
for (i in 1:nrow(species.table)) {
  model <- coll.glm.spring[[i]]
  data <- data.frame("occ"=mean(model$data[[5]]),"tvol"=mean(model$data$tvol),"tspd"=model$data$tspd,"length"=1)
  colnames(data)[1] <- species.table[i,2]
  temp_df <- data.frame(x=data$tspd,
                        y=predict(model, data, type="response"),
                        name=rep(paste(species.names[i]), each=length(data$tspd))
  )
  tspd.spring <- rbind(tspd.spring,temp_df)
  rm(data)
  rm(model)
  rm(temp_df)
}
rm(coll.glm.spring)
tspd.spring <- cbind(tspd.spring,season="spring")

tspd.all <- rbind(tspd.summer,tspd.autumn,tspd.winter,tspd.spring)

for(i in species.names) {
  tiff(paste0('figs/tspd_seasons_',c('a','b','c','d','e','f')[which(species.names==i)],'.tif'), pointsize = 12, compression = "lzw", res=300, width = 1100, height = 900)
  print(
    ggplot(tspd.all[tspd.all$name==i,],aes(x=x,y=y,group=season,shape=factor(season))) +
      geom_line(size=0.3,aes(linetype=season)) +
      ylab("Relative Collision Rate") +
      xlab("Traffic Speed (km/hour)") +
      labs(linetype = "Season") +
      labs(title = i) +
      theme_bw() +
      theme(legend.key = element_blank(),legend.key.size = unit(1.0, 'lines')) +
      theme(plot.margin=unit(c(.5,0,.1,.1),"cm")) +
      theme(axis.title.x = element_text(margin=unit(c(.3,0,0,0),"cm"))) +
      theme(axis.title.y = element_text(margin=unit(c(0,.3,0,0),"cm"))) +
      theme(panel.grid.major = element_line(size=0.1),panel.grid.minor = element_line(size=0.1)) +
      theme(text = element_text(size = 10)) +
      scale_x_continuous(breaks=seq(40,110,by=10), expand = c(0, 0), lim=c(40,110))
  )
  dev.off()
}
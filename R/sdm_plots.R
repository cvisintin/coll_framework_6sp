require(maptools)
require(ggplot2)
require(gbm)

species.table <- read.delim("data/species_list.csv", header=T, sep=",")

species.names <- c("Eastern Grey Kangaroo","Common Brushtail Possum","Common Ringtail Possum","Swamp Wallaby","Common Wombat","Koala")

sdm.colors = colorRampPalette(c("white","red"))
sdm.bw = colorRampPalette(c("white","black"))

plotPal <- c("#94d1c7", "#cccc2b", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#969696", "#bc80bd")

load(file = "data/brt_models_simp")
load(file = "output/sac_occ")

victoria <- readShapePoly("data/VIC_GDA9455_ADMIN_STATE.shp")

#create prediction plots for all species
for(i in 1:nrow(species.table)) {
  assign(paste(toupper(species.table[i,2]),sep=""),raster(paste("output/",toupper(species.table[i,2]),"_preds_brt.tif",sep="")))
  r <- get(paste(toupper(species.table[i,2])))
  png(paste("figs/",toupper(species.table[i,2]),"_preds_brt.png",sep=""), bg = "transparent", width = 1000, height = 700, pointsize = 24)
  par(mar=c(0,0,0,0)+0.0)
  plot(r, col=sdm.bw(100), axes=FALSE, box=FALSE, legend=FALSE, useRaster=FALSE, frame.plot=FALSE, legend.width=0, legend.mar=0)
  plot(victoria, add=TRUE, lwd=0.5)
  dev.off()
}

#plot spatial autocorrelation across all species
shapes <- unlist(lapply(c("1", "2", "3", "4", "5", "6"), utf8ToInt))

auto <- cbind(auto.occ,"name"=rep(species.names, each=20))

for(i in species.names) {
  tiff(paste0('figs/sac_',c('a','b','c','d','e','f')[which(species.names==i)],'.tif'), pointsize = 12, compression = "lzw", res=300, width = 900, height = 900)
  print(
    ggplot(auto[auto$name==i,],aes(x=x,y=y,group=name,shape=name)) + 
      geom_line(colour=c("grey70"),size=.5) + 
      geom_point(size=1) + 
      ylab("Moran's I") + 
      xlab("Distance (km)") +
      labs(title = i) +
      theme_bw() +
      theme(legend.key = element_blank(),legend.position = 'none') +
      theme(plot.margin=unit(c(.5,.1,.1,.1),"cm")) +
      theme(axis.title.x = element_text(margin=unit(c(.3,0,0,0),"cm"))) +
      theme(axis.title.y = element_text(margin=unit(c(0,.3,0,0),"cm"))) +
      theme(panel.grid.major = element_line(size=0.1),panel.grid.minor = element_line(size=0.1)) +
      theme(text = element_text(size = 8)) +
      geom_hline(aes(yintercept=0), linetype=2) + 
      scale_x_continuous(breaks=seq(1, 20, 1))
  )
  dev.off()
}


#plot effect of annual temperature range for all species
tempanrange <- NULL
for (i in 1:nrow(species.table)) {
  data <- read.delim(paste("data/",species.table[i,2],".data",sep=""), header=T, sep=",")
  model <- brt.models.simp[[i]]
  values <- plot.gbm(model, i.var="TEMPANRANGE",return.grid=TRUE, type="response")
  colnames(values) <- c("x","y")
  values$x <- values$x/10
  values$name <- as.factor(rep(paste(species.names[i]), each=length(values[,2])))
  tempanrange <- rbind(tempanrange,values)
  rm(data)
  rm(model)
  rm(values)
}  

tiff('figs/tempanrange.tif', pointsize = 11, compression = "lzw", res=300, width = 1500, height = 900)
ggplot(tempanrange,aes(x=x,y=y,group=name,colour=name)) + 
  geom_line(size=0.75) +  
  ylab("Occurence (Pr)") + 
  xlab(expression(paste("Annual Range of Temperature (",degree,"C)",sep=""))) + 
  labs(color = "Species") + 
  theme_bw() + 
  theme(legend.key = element_blank()) + 
  theme(plot.margin=unit(c(.5,.5,.1,.1),"cm")) +
  theme(axis.title.x = element_text(margin=unit(c(.3,0,0,0),"cm"))) +
  theme(axis.title.y = element_text(margin=unit(c(0,.3,0,0),"cm"))) +
  theme(panel.grid.major = element_line(size=0.1),panel.grid.minor = element_line(size=0.1)) +
  scale_colour_manual(values=plotPal) +
  theme(text = element_text(size = 11)) +
  scale_x_continuous(expand = c(0, 0), lim=c(15,30))
dev.off()


#plot effect of light for all species
light <- NULL
for (i in 1:nrow(species.table)) {
  data <- read.delim(paste("data/",species.table[i,2],".data",sep=""), header=T, sep=",")
  model <- brt.models.simp[[i]]
  values <- plot.gbm(model, i.var="LIGHT",return.grid=TRUE, type="response")
  colnames(values) <- c("x","y")
  values$name <- as.factor(rep(paste(species.names[i]), each=length(values[,2])))
  light <- rbind(light,values)
  rm(data)
  rm(model)
  rm(values)
}

tiff('figs/light.tif', pointsize = 11, compression = "lzw", res=300, width = 900, height = 900)
ggplot(light,aes(x=x,y=y,group=name,colour=name)) + 
  geom_line(size=0.75) + 
  ylab("Occurence (Pr)") + 
  xlab("Artificial Light (relative)") + 
  labs(color = "Species") + 
  theme_bw() + 
  theme(legend.key = element_blank(), legend.position="none") +
  theme(plot.margin=unit(c(.5,.5,.1,.1),"cm")) +
  theme(axis.title.x = element_text(margin=unit(c(.3,0,0,0),"cm"))) +
  theme(axis.title.y = element_text(margin=unit(c(0,.3,0,0),"cm"))) +
  theme(panel.grid.major = element_line(size=0.1),panel.grid.minor = element_line(size=0.1)) +
  scale_colour_manual(values=plotPal) +
  theme(text = element_text(size = 11)) +
  scale_x_continuous(expand = c(0, 0), lim=c(0,65))
dev.off()


#plot effect of precipitation of driest month for all species
precdm <- NULL
for (i in 1:nrow(species.table)) {
  data <- read.delim(paste("data/",species.table[i,2],".data",sep=""), header=T, sep=",")
  model <- brt.models.simp[[i]]
  values <- plot.gbm(model, i.var="PRECDM",return.grid=TRUE, type="response")
  colnames(values) <- c("x","y")
  values$name <- as.factor(rep(paste(species.names[i]), each=length(values[,2])))
  precdm <- rbind(precdm,values)
  rm(data)
  rm(model)
  rm(values)
}

tiff('figs/precdm.tif', pointsize = 11, compression = "lzw", res=300, width = 900, height = 900)
ggplot(precdm,aes(x=x,y=y,group=name,colour=name)) + 
  geom_line(size=0.75) + 
  ylab("Occurence (Pr)") + 
  xlab("Precipitation of Driest Month (mm)") + 
  labs(color = "Species") + 
  theme_bw() + 
  theme(legend.key = element_blank(), legend.position="none") + 
  theme(plot.margin=unit(c(.5,.5,.1,.1),"cm")) +
  theme(axis.title.x = element_text(margin=unit(c(.3,0,0,0),"cm"))) +
  theme(axis.title.y = element_text(margin=unit(c(0,.3,0,0),"cm"))) +
  theme(panel.grid.major = element_line(size=0.1),panel.grid.minor = element_line(size=0.1)) +
  scale_colour_manual(values=plotPal) +
  theme(text = element_text(size = 11)) +
  scale_x_continuous(expand = c(0, 0), lim=c(20,80))
dev.off()
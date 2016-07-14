require("ncf")
require("xtable")
require("ggplot2")

species.table <- read.delim("data/species_list.csv", header=T, sep=",")

sdm.colors = colorRampPalette(c("white","red"))
sdm.bw = colorRampPalette(c("white","black"))

load(file = "data/brt_models_simp")

#plot effect of greenness for all species
green <- NULL
for (i in 1:nrow(species.table)) {
  data <- read.delim(paste("data/",species.table[i,2],".data",sep=""), header=T, sep=",")
  model <- brt.models.simp[[i]]
  values <- plot.gbm(model, i.var="GREEN",return.grid=TRUE, type="response")
  colnames(values) <- c("x","y")
  values$col <- as.factor(rep(paste(toupper(species.table[i,2])), each=length(values[,2])))
  green <- rbind(green,values)
  rm(data)
  rm(model)
  rm(values)
}  
ggplot(green,aes(x=x,y=y,group=col,colour=col)) + 
  geom_line(size=1) + 
  ylab("Occurence (Pr)") + 
  xlab("Seasonal Change in Vegetation Greenness") + 
  labs(color = "Species") + 
  theme_bw() + 
  theme(legend.key = element_blank()) +
  theme(text = element_text(size = 20)) +
  scale_colour_manual(values=plotPal)


#plot effect of annual temperature range for all species
tempanrange <- NULL
for (i in 1:nrow(species.table)) {
  data <- read.delim(paste("data/",species.table[i,2],".data",sep=""), header=T, sep=",")
  model <- brt.models.simp[[i]]
  values <- plot.gbm(model, i.var="TEMPANRANGE",return.grid=TRUE, type="response")
  colnames(values) <- c("x","y")
  values$x <- values$x/10
  values$col <- as.factor(rep(paste(toupper(species.table[i,2])), each=length(values[,2])))
  tempanrange <- rbind(tempanrange,values)
  rm(data)
  rm(model)
  rm(values)
}  
ggplot(tempanrange,aes(x=x,y=y,group=col,colour=col)) + 
  geom_line(size=1) + 
  ylab("Occurence (Pr)") + 
  xlab(expression(paste("Annual Range of Temperature (",degree,"C)",sep=""))) + 
  labs(color = "Species") + 
  theme_bw() + 
  theme(legend.key = element_blank()) + 
  theme(text = element_text(size = 20)) +
  scale_colour_manual(values=plotPal)


#plot effect of elevation for all species
elev <- NULL
for (i in 1:nrow(species.table)) {
  data <- read.delim(paste("data/",species.table[i,2],".data",sep=""), header=T, sep=",")
  model <- brt.models.simp[[i]]
  values <- plot.gbm(model, i.var="ELEV",return.grid=TRUE, type="response")
  colnames(values) <- c("x","y")
  values$col <- as.factor(rep(paste(toupper(species.table[i,2])), each=length(values[,2])))
  elev <- rbind(elev,values)
  rm(data)
  rm(model)
  rm(values)
}  
ggplot(elev,aes(x=x,y=y,group=col,colour=col)) + 
  geom_line(size=1) + 
  ylab("Occurence (Pr)") + 
  xlab("Elevation (m above sea level)") + 
  labs(color = "Species") + 
  theme_bw() + 
  theme(legend.key = element_blank()) +
  theme(text = element_text(size = 20)) +
  scale_colour_manual(values=plotPal)


#plot effect of precipitation of driest month for all species
precdm <- NULL
for (i in 1:nrow(species.table)) {
  data <- read.delim(paste("data/",species.table[i,2],".data",sep=""), header=T, sep=",")
  model <- brt.models.simp[[i]]
  values <- plot.gbm(model, i.var="PRECDM",return.grid=TRUE, type="response")
  colnames(values) <- c("x","y")
  values$col <- as.factor(rep(paste(toupper(species.table[i,2])), each=length(values[,2])))
  precdm <- rbind(precdm,values)
  rm(data)
  rm(model)
  rm(values)
}  
ggplot(precdm,aes(x=x,y=y,group=col,colour=col)) + 
  geom_line(size=1) + 
  ylab("Occurence (Pr)") + 
  xlab("Precipitation of Driest Month (mm)") + 
  labs(color = "Species") + 
  theme_bw() + 
  theme(legend.key = element_blank()) + 
  theme(text = element_text(size = 20)) +
  scale_colour_manual(values=plotPal)


#plot effect of light for all species
light <- NULL
for (i in 1:nrow(species.table)) {
  data <- read.delim(paste("data/",species.table[i,2],".data",sep=""), header=T, sep=",")
  model <- brt.models.simp[[i]]
  values <- plot.gbm(model, i.var="LIGHT",return.grid=TRUE, type="response")
  colnames(values) <- c("x","y")
  values$col <- as.factor(rep(paste(toupper(species.table[i,2])), each=length(values[,2])))
  light <- rbind(light,values)
  rm(data)
  rm(model)
  rm(values)
}  
ggplot(light,aes(x=x,y=y,group=col,colour=col)) + 
  geom_line(size=1) + 
  ylab("Occurence (Pr)") + 
  xlab("Artificial Light (relative)") + 
  labs(color = "Species") + 
  theme_bw() + 
  theme(legend.key = element_blank()) +
  theme(text = element_text(size = 20)) +
  scale_colour_manual(values=plotPal)


#plot effect of tree density for all species
tree <- NULL
for (i in 1:nrow(species.table)) {
  data <- read.delim(paste("data/",species.table[i,2],".data",sep=""), header=T, sep=",")
  model <- brt.models.simp[[i]]
  values <- plot.gbm(model, i.var="TREEDENS",return.grid=TRUE, type="response")
  colnames(values) <- c("x","y")
  values$col <- as.factor(rep(paste(toupper(species.table[i,2])), each=length(values[,2])))
  tree <- rbind(tree,values)
  rm(data)
  rm(model)
  rm(values)
}  
ggplot(tree,aes(x=x,y=y,group=col,colour=col)) + 
  geom_line(size=1) + 
  ylab("Occurence (Pr)") + 
  xlab("Tree Density (percent coverage)") + 
  labs(color = "Species") + 
  theme_bw() + 
  theme(legend.key = element_blank()) +
  theme(text = element_text(size = 20)) +
  scale_colour_manual(values=plotPal)
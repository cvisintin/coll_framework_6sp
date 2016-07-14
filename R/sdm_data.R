require("maptools")
require("raster")
require("ncf")
require("xtable")
require("gbm")
require("ggplot2")
require("doMC")
require("dismo")
require("rgdal")
require("rgeos")
require("RPostgreSQL")

drv <- dbDriver("PostgreSQL")  #Specify a driver for postgreSQL type database
con <- dbConnect(drv, dbname="qaeco_spatial", user="qaeco", password="Qpostgres15", host="boab.qaeco.com", port="5432")  #Connection to database server on Boab

species.table <- read.delim("data/species_list.csv", header=T, sep=",")

plotPal <- c("#8dd3c7", "#e5e500", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#969696", "#bc80bd")

sdm.colors = colorRampPalette(c("white","red"))
sdm.bw = colorRampPalette(c("white","black"))

setwd('data/')

ascii.files <- list.files(pattern = "\\.asc$")

for(i in list(c('ncols',1),c('nrows',2),c('x.corner',3),c('y.corner',4))) {
  assign(i[1],as.numeric(scan(ascii.files[1],nlines=1,skip=as.numeric(i[2])-1,what="complex",quiet=T)[2]))
}

ascii.names <- unlist(strsplit(ascii.files,"\\."))[(1:(2*(length(ascii.files)))*2)-1][1:length(ascii.files)]

################

grid.files <- list.files(path='data/grids') #Create vector of filenames

grid.names <- unlist(strsplit(grid.files,"\\."))[(1:(2*(length(grid.files)))*2)-1][1:length(grid.files)] #Create vector of covariate names

vic.rst <- raster("data/VIC_GDA9455_GRID_STATE_1000.tif")

clip <- extent(-58000, 764000, 5661000, 6224000) #Define clipping extent of maps

#Read in grids, crop, and multiply with template to create consistent covariate maps
for (i in 1:length(grid.files)) {
  temp <- raster(paste0("data/grids/",grid.files[i]))
  temp <- crop(temp, clip)
  assign(grid.names[i],temp * vic.rst)
}
vars <- stack(mget(grid.names)) #Combine all maps to single stack

vars.cor <- layerStats(vars, 'pearson', na.rm=TRUE)

#write.csv(vars.cor[[1]], file = "/home/casey/Research/Projects/SDMs/Data/vars_cor.csv")

################

#extract environmental variables and write model data files
data0 <- read.csv("data/bg_data_pts.csv")
colnames(data0) <- c("x","y","occ")

for(i in 1:nrow(species.table)) {
  #data1 <- read.delim(paste("/home/casey/Research/Projects/SDMs/Data/",species.table[i,3],".csv",sep=""), header=T, sep=",")
  data1 <- dbGetQuery(con,paste0("
                                 SELECT
                                 ST_X(pts.geom) AS X, ST_Y(pts.geom) AS Y, CAST(1 AS INTEGER) AS OCC
                                 FROM
                                 gis_victoria.vic_gda9455_fauna_vba AS pts, gis_victoria.vic_gda9455_admin_state AS poly
                                 WHERE
                                 ST_contains(poly.geom, pts.geom)
                                 AND
                                 pts.start_year >= '2000'
                                 AND
                                 sci_name = '",species.table[i,3],"'
                                 GROUP BY
                                 pts.geom;
                                 "))
  data <- rbind(data1,data0)
  samples.df <- extract(vars,data[,1:2])
  data <- cbind(data,samples.df)
  data <- na.omit(data)
  write.csv(data, file = paste("data/",species.table[i,2],".data",sep=""), row.names=FALSE)
  rm(data1)
  rm(data)
}
rm(data0)


#construct metadata for all species
for(i in 1:nrow(species.table)) {
  assign(paste(species.table[i,3],".data",sep=""),read.delim(paste("data/",species.table[i,3],".data",sep=""), header=T, sep=","))
}

brt_data <- data.frame("SPP"=rep(NA,nrow(species.table)),"SO_N"=rep(NA,nrow(species.table)),"SO_P"=rep(NA,nrow(species.table)),"SO_A"=rep(NA,nrow(species.table)))
for(i in 1:nrow(species.table)) {
  brt_data[i,"SPP"] <- toupper(paste(species.table[i,3]))
  data <- get(paste(species.table[i,3],".data",sep=""))
  brt_data[i,"SO_N"] <- nrow(data)
  brt_data[i,"SO_P"] <- nrow(data[data$occ == 1, ])
  brt_data[i,"SO_A"] <- nrow(data[data$occ == 0, ])
  rm(data)
}
write.csv(brt_data, file = "data/brt_data_meta.csv", row.names=FALSE)


#fit BRT models for all species
registerDoMC(detectCores() - 1)

brt.models <- foreach(i = 1:nrow(species.table), .packages = c("gbm","dismo")) %dopar% {
  data <- read.delim(paste("data/",species.table[i,3],".data",sep=""), header=T, sep=",")
  set.seed(123)
  model <- gbm.step(data = data, gbm.x = c(4:21), gbm.y = 3, family = "bernoulli", tree.complexity = 5, learning.rate = 0.005, bag.fraction = 0.5)
  #model <- gbm.step(data = data, gbm.x = c(4:43), gbm.y = 3, family = "bernoulli", tree.complexity = 5, learning.rate = 0.005, bag.fraction = 0.5)
  model
}

#fit Random Forest models for all species
# registerDoMC(detectCores() - 1)
# 
# rf.models <- foreach(i = 1:nrow(species.table), .packages = c("gbm","dismo")) %dopar% {
#   data <- read.delim(paste("/home/casey/Research/Projects/SDMs/Data/",species.table[i,3],".data",sep=""), header=T, sep=",")
#   set.seed(123)
#   model2 <- randomForest(x=data[4:21], y=data[[3]], ntree=3850)
#   model2
# }


#fit Maxent models for all species
# registerDoMC(detectCores() - 1)
# 
#   max.models <- foreach(i = 1:nrow(species.table), .packages = c("dismo","rJava")) %dopar% {
#     data <- read.delim(paste("/home/casey/Research/Projects/SDMs/Data/",species.table[i,3],".data",sep=""), header=T, sep=",")
#     data.bg <- data.frame("X"=x0[,1], "Y"=x0[,2])
#     data.pres <- data.frame("X"=data$X[data[,3]==1], "Y"=data$Y[data[,3]==1])
#     model <- maxent(vars, data.pres, a=data.bg, args=c("-P", "noautofeature", "nothreshold", "noproduct", -t soil, "noaddsamplestobackground","noremoveduplicates"))
#     model
#   }

#egk.ll <- SpatialPoints(cbind(raw.data$decimallongitude, raw.data$decimallatitude), proj4string=CRS("+init=epsg:4283"))
#egk.UTM <- data.frame(spTransform(egk.ll, CRS("+init=epsg:28355")))
#names(egk.UTM) <- c('X','Y')
#egk.raster <- rasterize(egk.UTM, r)
#egk.raster[egk.raster>=1] <- 1
#egk1 <- rasterToPoints(egk.raster)
#egk1 <- data.frame("X"=egk1[,1], "Y"=egk1[,2])
#egk0 <- data.frame("X"=x0[,1], "Y"=x0[,2])

#make predictions for all species (Maxent)
# registerDoMC(detectCores() - 1)
# 
#   max.models.output <- foreach(i = 1:nrow(species.table), .packages = c("dismo","doParallel","raster","rJava")) %dopar% {
#     data <- read.delim(paste("/home/casey/Research/Projects/SDMs/Data/",species.table[i,3],".data",sep=""), header=T, sep=",")
#     data.bg <- data.frame("X"=x0[,1], "Y"=x0[,2])
#     data.pres <- data.frame("X"=data$X[data[,3]==1], "Y"=data$Y[data[,3]==1])
#     model <- max.models[[i]]
#     max.preds <- predict(model, vars)
#     writeRaster(max.preds, filename=paste("/home/casey/Research/Projects/SDMs/Data/Preds/",toupper(species.table[i,3]),".asc",sep=""), format="ascii", overwrite=TRUE)
#     #plot(max.preds, col=sdm.colors(100), axes=F, box=FALSE)
#     summary(model)
#   }

#simplify BRT models and build predictor lists for all species
registerDoMC(detectCores() - 1)

predlists <- foreach(i = 1:nrow(species.table), .packages = c("gbm","dismo"), .export=("gbm.step")) %dopar% {
  data <- read.delim(paste("/home/casey/Research/Projects/SDMs/Data/",species.table[i,3],".data",sep=""), header=T, sep=",")
  set.seed(123)
  model.simp <- gbm.simplify(brt.models[[i]])
  preds.no <- which.min(model.simp$deviance.summary$mean)
  model.simp[['pred.list']][[which.min(model.simp$deviance.summary$mean)]]
}
save(predlists, file = "/home/casey/Research/Projects/SDMs/Data/predlists")
#load(file = "/home/casey/Research/Projects/SDMs/Data/predlists")


# for(i in 1:nrow(species.table)) {
#   data <- read.delim(paste("/home/casey/Research/Projects/SDMs/Data/",species.table[i,3],".data",sep=""), header=T, sep=",")
#   cor <- cor(data[,predlists[[i]]])
#   write.csv(cor, file = paste("/home/casey/Research/Projects/SDMs/Data/",species.table[i,3],".cor",sep=""), row.names=FALSE)
# }


#fit final BRT models for all species
registerDoMC(detectCores() - 1)

brt.models.simp <- foreach(i = 1:nrow(species.table), .packages = c("gbm","dismo")) %dopar% {
  data <- read.delim(paste("/home/casey/Research/Projects/SDMs/Data/",species.table[i,3],".data",sep=""), header=T, sep=",")
  set.seed(123)
  model <- gbm.step(data = data, gbm.x = predlists[[i]], gbm.y = 3, family = "bernoulli", tree.complexity = 5, learning.rate = 0.005, bag.fraction = 0.5)
  model
}
save(brt.models.simp, file = "/home/casey/Research/Projects/SDMs/Data/brt_models")
#load(file = "/home/casey/Research/Projects/SDMs/Data/brt_models")

#make predictions for all species
registerDoMC(detectCores() - 1)

brt.models.output <- foreach(i = 1:nrow(species.table), .packages = c("gbm","dismo","raster")) %dopar% {
  data <- read.delim(paste("/home/casey/Research/Projects/SDMs/Data/",species.table[i,3],".data",sep=""), header=T, sep=",")
  model <- brt.models.simp[[i]]
  #model <- brt.models[[i]]
  model.preds <- predict(vars, model, n.trees=model[["gbm.call"]][["best.trees"]], type="response")
  writeRaster(model.preds, filename=paste("/home/casey/Research/Projects/SDMs/Data/Preds/",toupper(species.table[i,3]),".asc",sep=""), format="ascii", overwrite=TRUE)
  plot(model.preds, col=sdm.colors(100), axes=F, box=FALSE)
  #summary(model)
}

for(i in 1:nrow(species.table)) {
  assign(paste(toupper(species.table[i,3]),sep=""),raster(paste("/home/casey/Research/Projects/SDMs/Data/Preds/",toupper(species.table[i,3]),".asc",sep="")))
  png(paste("/home/casey/Research/Graphics_Repo/",toupper(species.table[i,3]),"_SDM2.png",sep=""), bg = "white", width = 1000, height = 700, pointsize = 24)
  par(mar=c(0,0,0,0)+0.0)
  plot(get(paste(toupper(species.table[i,3]))), col=sdm.colors(100), axes=F, box=FALSE)
  dev.off()
}

#summarize model output data
brt_deviance <- data.frame("SPP"=rep(NA,nrow(species.table)),"DEV"=rep(NA,nrow(species.table)),"ERR"=rep(NA,nrow(species.table)),"ROC"=rep(NA,nrow(species.table)),"ROCERR"=rep(NA,nrow(species.table)))
for(i in 1:nrow(species.table)) {
  brt_deviance[i,"SPP"] <- toupper(paste(species.table[i,3]))
  brt_deviance[i,"DEV"] <- ((brt.models.simp[i][[1]][["self.statistics"]][["mean.null"]] - brt.models.simp[i][[1]][["cv.statistics"]][["deviance.mean"]])/brt.models.simp[i][[1]][["self.statistics"]][["mean.null"]])*100
  #brt_deviance[i,"DEV"] <- ((brt.models[i][[1]][["self.statistics"]][["mean.null"]] - brt.models[i][[1]][["cv.statistics"]][["deviance.mean"]])/brt.models[i][[1]][["self.statistics"]][["mean.null"]])*100
  brt_deviance[i,"ERR"] <- (brt.models.simp[i][[1]][["cv.statistics"]][["deviance.se"]]/brt.models.simp[i][[1]][["self.statistics"]][["mean.null"]])*100
  #brt_deviance[i,"ERR"] <- (brt.models[i][[1]][["cv.statistics"]][["deviance.se"]]/brt.models[i][[1]][["self.statistics"]][["mean.null"]])*100
  brt_deviance[i,"ROC"] <- brt.models.simp[i][[1]][["cv.statistics"]][["discrimination.mean"]]
  #brt_deviance[i,"ROC"] <- brt.models[i][[1]][["cv.statistics"]][["discrimination.mean"]]
  brt_deviance[i,"ROCERR"] <- brt.models.simp[i][[1]][["cv.statistics"]][["discrimination.se"]]
  #brt_deviance[i,"ROCERR"] <- brt.models[i][[1]][["cv.statistics"]][["discrimination.se"]]
}
write.csv(brt_deviance, file = "/home/casey/Research/Projects/SDMs/Data/brt_devs.csv", row.names=FALSE)


brt_sums <- data.frame(rep(NA,length(ascii.names)))
colnames(brt_sums) <- c("Predictor")
brt_sums$Predictor <- ascii.names[1:(length(ascii.names))]
for(i in 1:nrow(species.table)) {
  x.var <- data.frame(brt.models.output[i],stringsAsFactors=FALSE)
  x.var[2] <- round(x.var[2],2)
  brt_sums <- merge(brt_sums,x.var,by.x="Predictor",by.y="var",all.x=TRUE)
  colnames(brt_sums)[i+1] <- toupper(paste(species.table[i,3]))
  brt_sums[is.na(brt_sums)] <- 0
  #brt_sums[brt_sums==max(na.omit(brt_sums[i+1]))] <- paste("\\B{",max(na.omit(brt_sums[i+1])),"}",sep="")
}

#check total influence of variables across all species
nums <- sapply(brt_sums, is.numeric)
brt_sums$Total <- as.integer(rep(0,length(ascii.names)))
brt_sums$wTotal <- as.integer(rep(0,length(ascii.names)))
count <- rowSums(brt_sums!=0)-1
for(i in 1:nrow(brt_sums)) {
  brt_sums[i,"Total"] <- sum(brt_sums[i, nums])
  brt_sums[i,"wTotal"] <- signif(brt_sums[i,"Total"]*count[i]/10,digits=4)
} 
brt_sums_wt <- brt_sums

brt_sums$Total <- NULL
brt_sums$wTotal <- NULL


for(i in 2:ncol(brt_sums)) {
  sums <- as.numeric(brt_sums[,i])
  sums[which(sums == max(sums), arr.ind=TRUE)] <- paste("\\B{",max(sums),"}",sep="")
  brt_sums[,i] <- sums
}

#construct LaTex table 
brt_sums[brt_sums==0] <- "---"
print(xtable(brt_sums), include.rownames=FALSE, sanitize.text.function=function(x){x}, floating=FALSE)  

write.csv(brt_sums, file = "/home/casey/Research/Projects/SDMs/Data/brt_sums.csv", row.names=FALSE)


#calculate spatial autocorrelation across all species

auto <- foreach(i = 1:nrow(species.table), .packages = c("ncf"), .combine="rbind") %do% {
  data <- read.delim(paste("/home/casey/Research/Projects/SDMs/Data/",species.table[i,3],".data",sep=""), header=T, sep=",")
  #model <- mget(paste(species.table[i,3],".brt",sep=""))[[paste(species.table[i,3],".brt",sep="")]]
  model <- brt.models.simp[[i]]
  #model <- brt.models[[i]]
  cor <- correlog(data[,1], data[,2], resid(model), increment=1000, resamp=0, latlon=FALSE)
  temp_df <- data.frame(x=as.numeric(names(cor$correlation[2:21])), y=cor$correlation[2:21], col=rep(toupper(species.table[i,3]), each=length(cor$correlation[2:21])))
  temp_df
}  

shapes <- unlist(lapply(c("1", "2", "3", "4", "5", "6", "7"), utf8ToInt))

ggplot(auto,aes(x=x,y=y,group=col,shape=col)) + 
  geom_line(colour=c("grey70"),size=.75) + 
  geom_point(size=3) + 
  ylab("Moran's I") + 
  xlab("Distance (km)") + 
  labs(shape = "Species") + 
  theme_bw() + 
  theme(legend.key = element_blank()) +
  theme(text = element_text(size = 16)) +
  scale_colour_manual(values=plotPal) + 
  scale_shape_manual(values=shapes) + 
  geom_hline(aes(yintercept=0), linetype=2) + 
  scale_x_continuous(breaks=seq(1, 20, 1))


#plot effect of greenness for all species
green <- NULL
for (i in 1:nrow(species.table)) {
  data <- read.delim(paste("/home/casey/Research/Projects/SDMs/Data/",species.table[i,3],".data",sep=""), header=T, sep=",")
  model <- brt.models.simp[[i]]
  values <- plot.gbm(model, i.var="GREEN",return.grid=TRUE, type="response")
  colnames(values) <- c("x","y")
  values$col <- as.factor(rep(paste(toupper(species.table[i,3])), each=length(values[,2])))
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
  data <- read.delim(paste("/home/casey/Research/Projects/SDMs/Data/",species.table[i,3],".data",sep=""), header=T, sep=",")
  model <- brt.models.simp[[i]]
  values <- plot.gbm(model, i.var="TEMPANRANGE",return.grid=TRUE, type="response")
  colnames(values) <- c("x","y")
  values$x <- values$x/10
  values$col <- as.factor(rep(paste(toupper(species.table[i,3])), each=length(values[,2])))
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
  data <- read.delim(paste("/home/casey/Research/Projects/SDMs/Data/",species.table[i,3],".data",sep=""), header=T, sep=",")
  model <- brt.models.simp[[i]]
  values <- plot.gbm(model, i.var="ELEV",return.grid=TRUE, type="response")
  colnames(values) <- c("x","y")
  values$col <- as.factor(rep(paste(toupper(species.table[i,3])), each=length(values[,2])))
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
  data <- read.delim(paste("/home/casey/Research/Projects/SDMs/Data/",species.table[i,3],".data",sep=""), header=T, sep=",")
  model <- brt.models.simp[[i]]
  values <- plot.gbm(model, i.var="PRECDM",return.grid=TRUE, type="response")
  colnames(values) <- c("x","y")
  values$col <- as.factor(rep(paste(toupper(species.table[i,3])), each=length(values[,2])))
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
  data <- read.delim(paste("/home/casey/Research/Projects/SDMs/Data/",species.table[i,3],".data",sep=""), header=T, sep=",")
  model <- brt.models.simp[[i]]
  values <- plot.gbm(model, i.var="LIGHT",return.grid=TRUE, type="response")
  colnames(values) <- c("x","y")
  values$col <- as.factor(rep(paste(toupper(species.table[i,3])), each=length(values[,2])))
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
  data <- read.delim(paste("/home/casey/Research/Projects/SDMs/Data/",species.table[i,3],".data",sep=""), header=T, sep=",")
  model <- brt.models.simp[[i]]
  values <- plot.gbm(model, i.var="TREEDENS",return.grid=TRUE, type="response")
  colnames(values) <- c("x","y")
  values$col <- as.factor(rep(paste(toupper(species.table[i,3])), each=length(values[,2])))
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


#model based on residual autocorrelation
rast_brt <- r

xy <- cbind(egk.data$x, egk.data$y)  

xy_res_brt <- cbind(xy, resid(brt.models[[1]]))

rast_brt[cellFromXY(rast_brt, xy_res_brt)] <- xy_res_brt[,3]

focal_rac_rast_brt <- focal(rast_brt, w=matrix(1,3,3), fun = mean, na.rm = TRUE)

RAC <- extract(focal_rac_rast_brt, xy)

egk.data.spatial <- cbind(egk.data,RAC)

set.seed(123)
egk.brt.spatial <- gbm.step(data = egk.data.spatial, gbm.x = c(predlists[[i]],22), gbm.y = 3, family = "bernoulli", tree.complexity = 5, learning.rate = 0.005, bag.fraction = 0.5)


require("spdep")

nb.list <- dnearneigh(as.matrix(egk.data[,c("x", "y")]), 0, 3000)

nb.weights <- nb2listw(nb.list, zero.policy=TRUE)

nbd <- nbdists(nb.list, xy)
gl <- lapply(nbd, function(x) 1/x)
lw <- nb2listw(nb.list, glist=gl, zero.policy=TRUE)
AC <- lag(lw, egk.data$occ)

AC[is.na(AC)] <- 0

egk.data.spatial <- cbind(egk.data.spatial,AC)


rast_brt <- r

xy <- cbind(egk.data$x, egk.data$y)  

xy_res_brt <- cbind(xy, egk.data$occ)

rast_brt[cellFromXY(rast_brt, xy_res_brt)] <- xy_res_brt[,3]

focal_rac_rast_brt <- focal(rast_brt, w=matrix(1,3,3), fun = sum, na.rm = TRUE)

AC <- extract(focal_rac_rast_brt, xy)

egk.data.spatial <- cbind(egk.data.spatial,AC)

set.seed(123)
egk.brt.spatial3 <- gbm.step(data = egk.data.spatial, gbm.x = c(predlists[[i]],22), gbm.y = 3, family = "bernoulli", tree.complexity = 5, learning.rate = 0.005, bag.fraction = 0.5)



egk.preds.spatial1 <- predict(stack(vars,"AC"=focal_rac_rast_brt), model, n.trees=model[["gbm.call"]][["best.trees"]], type="response")
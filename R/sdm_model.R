require("ncf")
require("xtable")
require("gbm")
require("doMC")
require("dismo")
require("ggplot2")

species.table <- read.delim("data/species_list.csv", header=T, sep=",")

sdm.colors = colorRampPalette(c("white","red"))
sdm.bw = colorRampPalette(c("white","black"))

#fit BRT models for all species
registerDoMC(detectCores() - 1)

  brt.models <- foreach(i = 1:nrow(species.table), .packages = c("gbm","dismo")) %dopar% {
    data <- read.delim(paste("data/",species.table[i,2],".data",sep=""), header=T, sep=",")
    set.seed(123)
    model <- gbm.step(data = data, gbm.x = c(4:21), gbm.y = 3, family = "bernoulli", tree.complexity = 5, learning.rate = 0.005, bag.fraction = 0.5)
    #model <- gbm.step(data = data, gbm.x = c(4:43), gbm.y = 3, family = "bernoulli", tree.complexity = 5, learning.rate = 0.005, bag.fraction = 0.5)
    model
  }
save(brt.models, file = "data/brt_models")
#load(file = "data/brt_models")
  
#simplify BRT models and build predictor lists for all species
registerDoMC(detectCores() - 1)

predlists <- foreach(i = 1:nrow(species.table), .packages = c("gbm","dismo"), .export=("gbm.step")) %dopar% {
  data <- read.delim(paste("data/",species.table[i,2],".data",sep=""), header=T, sep=",")
  set.seed(123)
  model.simp <- gbm.simplify(brt.models[[i]])
  preds.no <- which.min(model.simp$deviance.summary$mean)
  model.simp[['pred.list']][[which.min(model.simp$deviance.summary$mean)]]
}
save(predlists, file = "data/predlists")
#load(file = "data/predlists")


for(i in 1:nrow(species.table)) {
  data <- read.delim(paste("data/",species.table[i,2],".data",sep=""), header=T, sep=",")
  cor <- cor(data[,predlists[[i]]])
  write.csv(cor, file = paste("data/",species.table[i,2],".cor",sep=""), row.names=FALSE)
}


#fit final BRT models for all species
registerDoMC(detectCores() - 1)

brt.models.simp <- foreach(i = 1:nrow(species.table), .packages = c("gbm","dismo")) %dopar% {
  data <- read.delim(paste("data/",species.table[i,2],".data",sep=""), header=T, sep=",")
  set.seed(123)
  model <- gbm.step(data = data, gbm.x = predlists[[i]], gbm.y = 3, family = "bernoulli", tree.complexity = 5, learning.rate = 0.005, bag.fraction = 0.5)
  model
}
save(brt.models.simp, file = "data/brt_models_simp")
#load(file = "data/brt_models")


#make predictions for all species
registerDoMC(detectCores() - 1)

brt.models.output <- foreach(i = 1:nrow(species.table), .packages = c("gbm","dismo","raster")) %dopar% {
  data <- read.delim(paste("data/",species.table[i,2],".data",sep=""), header=T, sep=",")
  model <- brt.models.simp[[i]]
  #model <- brt.models[[i]]
  model.preds <- predict(vars, model, n.trees=model[["gbm.call"]][["best.trees"]], type="response")
  writeRaster(model.preds, filename=paste("output/",toupper(species.table[i,2]),"_preds.asc",sep=""), format="ascii", overwrite=TRUE)
  plot(model.preds, col=sdm.colors(100), axes=F, box=FALSE)
  #summary(model)
}

for(i in 1:nrow(species.table)) {
  assign(paste(toupper(species.table[i,2]),sep=""),raster(paste("output/",toupper(species.table[i,2]),".asc",sep="")))
  png(paste("figs/",toupper(species.table[i,2]),"_preds.png",sep=""), bg = "transparent", width = 1000, height = 700, pointsize = 24)
  par(mar=c(0,0,0,0)+0.0)
  plot(get(paste(toupper(species.table[i,2]))), col=sdm.colors(100), axes=F, box=FALSE)
  dev.off()
}

#summarize model output data
brt_deviance <- data.frame("SPP"=rep(NA,nrow(species.table)),"DEV"=rep(NA,nrow(species.table)),"ERR"=rep(NA,nrow(species.table)),"ROC"=rep(NA,nrow(species.table)),"ROCERR"=rep(NA,nrow(species.table)))
for(i in 1:nrow(species.table)) {
  brt_deviance[i,"SPP"] <- toupper(paste(species.table[i,2]))
  brt_deviance[i,"DEV"] <- ((brt.models.simp[i][[1]][["self.statistics"]][["mean.null"]] - brt.models.simp[i][[1]][["cv.statistics"]][["deviance.mean"]])/brt.models.simp[i][[1]][["self.statistics"]][["mean.null"]])*100
  #brt_deviance[i,"DEV"] <- ((brt.models[i][[1]][["self.statistics"]][["mean.null"]] - brt.models[i][[1]][["cv.statistics"]][["deviance.mean"]])/brt.models[i][[1]][["self.statistics"]][["mean.null"]])*100
  brt_deviance[i,"ERR"] <- (brt.models.simp[i][[1]][["cv.statistics"]][["deviance.se"]]/brt.models.simp[i][[1]][["self.statistics"]][["mean.null"]])*100
  #brt_deviance[i,"ERR"] <- (brt.models[i][[1]][["cv.statistics"]][["deviance.se"]]/brt.models[i][[1]][["self.statistics"]][["mean.null"]])*100
  brt_deviance[i,"ROC"] <- brt.models.simp[i][[1]][["cv.statistics"]][["discrimination.mean"]]
  #brt_deviance[i,"ROC"] <- brt.models[i][[1]][["cv.statistics"]][["discrimination.mean"]]
  brt_deviance[i,"ROCERR"] <- brt.models.simp[i][[1]][["cv.statistics"]][["discrimination.se"]]
  #brt_deviance[i,"ROCERR"] <- brt.models[i][[1]][["cv.statistics"]][["discrimination.se"]]
}
write.csv(brt_deviance, file = "data/brt_devs.csv", row.names=FALSE)


brt_sums <- data.frame(rep(NA,length(ascii.names)))
colnames(brt_sums) <- c("Predictor")
brt_sums$Predictor <- ascii.names[1:(length(ascii.names))]
for(i in 1:nrow(species.table)) {
  x.var <- data.frame(brt.models.output[i],stringsAsFactors=FALSE)
  x.var[2] <- round(x.var[2],2)
  brt_sums <- merge(brt_sums,x.var,by.x="Predictor",by.y="var",all.x=TRUE)
  colnames(brt_sums)[i+1] <- toupper(paste(species.table[i,2]))
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

#construct LaTex table 
for(i in 2:ncol(brt_sums)) {
  sums <- as.numeric(brt_sums[,i])
  sums[which(sums == max(sums), arr.ind=TRUE)] <- paste("\\B{",max(sums),"}",sep="")
  brt_sums[,i] <- sums
}

brt_sums[brt_sums==0] <- "---"
print(xtable(brt_sums), include.rownames=FALSE, sanitize.text.function=function(x){x}, floating=FALSE)  

write.csv(brt_sums, file = "data/brt_sums.csv", row.names=FALSE)


#calculate spatial autocorrelation across all species
auto <- foreach(i = 1:nrow(species.table), .packages = c("ncf"), .combine="rbind") %do% {
  data <- read.delim(paste("data/",species.table[i,2],".data",sep=""), header=T, sep=",")
  #model <- mget(paste(species.table[i,2],".brt",sep=""))[[paste(species.table[i,2],".brt",sep="")]]
  model <- brt.models.simp[[i]]
  #model <- brt.models[[i]]
  cor <- correlog(data[,1], data[,2], resid(model), increment=1000, resamp=0, latlon=FALSE)
  temp_df <- data.frame(x=as.numeric(names(cor$correlation[2:21])), y=cor$correlation[2:21], col=rep(toupper(species.table[i,2]), each=length(cor$correlation[2:21])))
  temp_df
} 

shapes <- unlist(lapply(c("1", "2", "3", "4", "5", "6"), utf8ToInt))

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
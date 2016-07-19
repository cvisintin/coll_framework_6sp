library("maptools")
library("raster")
library("ggplot2")
library("xtable")
library("R2jags")
library("ncf")
library("doMC")
library("data.table")

"roc" <- function (obsdat, preddat) {
    if (length(obsdat) != length(preddat)) 
      stop("obs and preds must be equal lengths")
    n.x <- length(obsdat[obsdat == 0])
    n.y <- length(obsdat[obsdat == 1])
    xy <- c(preddat[obsdat == 0], preddat[obsdat == 1])
    rnk <- rank(xy)
    wilc <- ((n.x * n.y) + ((n.x * (n.x + 1))/2) - sum(rnk[1:n.x]))/(n.x * n.y)
    return(round(wilc, 4))
}

plotPal <- c("#94d1c7", "#cccc2b", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#969696", "#bc80bd")

setwd('/home/casey/Research/Projects/SDMs/Data/Preds')

ascii.files <- list.files(pattern=".asc")

for(i in list(c('ncols',1),c('nrows',2),c('x.corner',3),c('y.corner',4))) {
  assign(i[1],as.numeric(scan(ascii.files[1],nlines=1,skip=as.numeric(i[2])-1,what="complex",quiet=T)[2]))
}

ascii.names <- unlist(strsplit(ascii.files,"\\."))[(1:(2*(length(ascii.files)))*2)-1][1:length(ascii.files)]

victoria <- readShapePoly("/home/casey/Research/GIS_Repo/VICTORIA/VIC_GDA9455_ADMIN_STATE.shp")

r <- raster(ncol=822, nrow=563, xmn=-58000, xmx=764000, ymn=5661000, ymx=6224000)

vic.rst <- rasterize(victoria, r, 'UFI')

clip <- extent(-58000, 764000, 5661000, 6224000)

for (i in 1:length(ascii.files)) {
  temp <- raster(ascii.files[i])
  temp <- crop(temp, clip)
  assign(ascii.names[i],temp * vic.rst)
}
vars <- stack(mget(ascii.names))

rm(list = ascii.names)
rm(vic.rst)

species.table <- read.delim("/home/casey/Research/Projects/Risk_Model/Data/species_list.csv", header=T, sep=",")

species.table <- species.table[1:6,]

species.list <- c("Eastern Grey Kangaroo","Common Brushtail Possum","Common Ringtail Possum","Black Swamp Wallaby","Common Wombat","Koala")
ascii.names <- ascii.names[match(toupper(species.table[,3]), ascii.names)]

# pb <- txtProgressBar(min = 0, max = nrow(species.table), style = 3)
# for (i in 1:length(ascii.files)){
#   sys.cmd.pres <- paste("spatialite -header -csv /home/casey/Research/Database_Repo/Spatialite/victoria.sqlite \"SELECT rds.UID AS UID, col.CASENUMBER AS CASENUMBER, col.YEAR AS YEAR, col.MONTH AS MONTH, col.DAY AS DAY, col.TIME AS TIME, X(Line_Interpolate_Point(rds.geometry, 0.5)) AS X, Y(Line_Interpolate_Point(rds.geometry, 0.5)) AS Y, '1' AS COLL, rds.TVOL AS TVOL, rds.TSPD AS TSPD FROM VIC_GDA9455_FAUNA_WV_20102014 AS col, VIC_GDA9455_ROADS_VICSTATE_GRID_CLIP_RDPREDS AS rds WHERE col.SPECIES LIKE '",species.table[i,1],"' AND col.CAUSE LIKE '%vehicle%' AND Intersects(col.geometry, Buffer(rds.geometry,10)) AND col.ROWID IN (SELECT ROWID FROM SpatialIndex WHERE f_table_name = 'VIC_GDA9455_FAUNA_WV_20102014' AND search_frame = Buffer(rds.geometry,10)) GROUP BY rds.UID; \" > /home/casey/Research/Projects/Risk_Model/Data/",species.table[i,3],"_1.csv",sep="")
#   system(sys.cmd.pres)
#   #sys.cmd.bg <- paste("spatialite -header -csv /home/casey/Research/Database_Repo/Spatialite/victoria.sqlite \"SELECT rds.UID AS UID, col.CASENUMBER AS CASENUMBER, col.YEAR AS YEAR, col.MONTH AS MONTH, col.DAY AS DAY, col.TIME AS TIME, X(Line_Interpolate_Point(rds.geometry, 0.5)) AS X, Y(Line_Interpolate_Point(rds.geometry, 0.5)) AS Y, '0' AS COLL, rds.TVOL AS TVOL, rds.TSPD AS TSPD FROM VIC_GDA9455_FAUNA_WV_20102014 AS col, VIC_GDA9455_ROADS_VICSTATE_GRID_CLIP_RDPREDS AS rds WHERE col.SPECIES NOT LIKE '",species.table[i,1],"' AND col.CAUSE LIKE '%vehicle%' AND Intersects(col.geometry, Buffer(rds.geometry,10)) AND col.ROWID IN (SELECT ROWID FROM SpatialIndex WHERE f_table_name = 'VIC_GDA9455_FAUNA_WV_20102014' AND search_frame = Buffer(rds.geometry,10)) GROUP BY rds.UID; \" > /home/casey/Research/Projects/Risk_Model/Data/",species.table[i,3],"_0.csv",sep="")
#   #system(sys.cmd.bg)
#   setTxtProgressBar(pb, i)
# }
# close(pb)

#sys.cmd.bg <- paste("spatialite -header -csv /home/casey/Research/Database_Repo/Spatialite/victoria.sqlite \"SELECT rds.UID AS UID, X(Line_Interpolate_Point(rds.geometry, 0.5)) AS X, Y(Line_Interpolate_Point(rds.geometry, 0.5)) AS Y, rds.TVOL AS TVOL, rds.TSPD AS TSPD FROM VIC_GDA9455_ROADS_VICSTATE_GRID_CLIP_RDPREDS AS rds; \" > /home/casey/Research/Projects/Risk_Model/Data/allroads_0.csv",sep="")
#system(sys.cmd.bg)

x0 <- read.delim(paste("/home/casey/Research/Projects/Risk_Model/Data/allroads_0.csv",sep=""), header=T, sep=",")
x0 <- cbind(x0[,1:3],"COLL"=rep(0,nrow(x0)),x0[,4:5])

#construct model datasets
for (i in 1:length(ascii.files)) {
  x1 <- read.delim(paste("/home/casey/Research/Projects/Risk_Model/Data/",tolower(ascii.names[i]),"_1.csv",sep=""), header=T, sep=",")
  x1 <- x1[,c(1,7:11)]
  #assign(paste(tolower(ascii.names[i]),".dist",sep=""), dist(as.matrix(cbind(x1$X,x1$Y)), method = "euclidean", diag = FALSE, upper = FALSE))
  set.seed(123)
  #x0 <- x0[sample(nrow(x0),3*nrow(x1)), ]
  #x0.dist <- mindistsamp(1000, x0[,2:3])
  #x0 <- merge(x0,x0.dist,by=c("X","Y"))
  x <- rbind(x1,x0)
  coord <- x[,2:3]
  samples.df <- extract(vars,coord)
  xx <- cbind(x,samples.df)
  xx <- xx[,c("X","Y","COLL","TVOL","TSPD",paste(ascii.names[i]))]
  assign(paste(tolower(ascii.names[i]),".data",sep=""), na.omit(xx))
  #write.csv(na.omit(xx), file = paste("/home/casey/Research/Projects/Risk_Model/Data/",tolower(ascii.names[i]),".data",sep=""), row.names=FALSE)
  rm(x)
  rm(x1)
  rm(xx)
  rm(samples.df)
  rm(coord)
}

#construct independent/prediction dataset
  coord <- x0[,2:3]
  samples.df <- extract(vars,coord)
  x <- cbind(x0,samples.df)
  xx <- na.omit(x[,c(-4)])
  xx[, "c.log.TVOL"] <- log(xx[,4])-mean(log(xx[,4]))
  xx[, "c.log.TSPD"] <- log(xx[,5])-mean(log(xx[,5]))
  for(i in 1:nrow(species.table)) {
    xx[, paste("c.log.",ascii.names[i],sep="")] <- log(xx[,i+5])-mean(log(xx[,i+5]))
  }
  assign(paste("all.ind.data",sep=""), xx)
  write.csv(na.omit(xx), file = paste("/home/casey/Research/Projects/Risk_Model/Data/all.ind.data",sep=""), row.names=FALSE)
  rm(x)  
  rm(xx)
  rm(samples.df)
  rm(coord)

#construct metadata for all species
coll_data <- data.frame("SPP"=rep(NA,nrow(species.table)),"CR_N"=rep(NA,nrow(species.table)),"CR_P"=rep(NA,nrow(species.table)),"CR_A"=rep(NA,nrow(species.table)))
for(i in 1:nrow(species.table)) {
  coll_data[i,"SPP"] <- toupper(paste(species.table[i,3]))
  colname <- "COLL"
  data <- get(paste(species.table[i,3],".data",sep=""))
  coll_data[i,"CR_N"] <- nrow(data)
  coll_data[i,"CR_P"] <- nrow(data[data[[colname]] == 1, ])
  coll_data[i,"CR_A"] <- nrow(data[data[[colname]] == 0, ])
  rm(colname)
  rm(data)
}
write.csv(coll_data, file = "/home/casey/Research/Projects/Risk_Model/Data/coll_data_meta.csv", row.names=FALSE)


#center and standardise all variables
for (i in 1:length(ascii.files)) {
  #x <- as.data.frame(mget(paste(tolower(ascii.names[i]),".data",sep=""))[[paste(tolower(ascii.names[i]),".data",sep="")]])
  x <- read.delim(paste("/home/casey/Research/Projects/Risk_Model/Data/",tolower(ascii.names[i]),".data",sep=""), header=T, sep=",")
  x[, "c.log.TVOL"] <- log(x[,4])-mean(log(x[,4]))
  x[, "c.log.TSPD"] <- log(x[,5])-mean(log(x[,5]))
  x[, paste("c.log.",ascii.names[i],sep="")] <- log(x[,6])-mean(log(x[,6]))
  assign(paste(tolower(ascii.names[i]),".data",sep=""), x)

  assign(paste("ml.TVOL.",i,sep=""), mean(log(x[,4])))
  assign(paste("ml.TSPD.",i,sep=""), mean(log(x[,5])))
  assign(paste("ml.",ascii.names[i],sep=""), mean(log(x[,6])))
  rm(x)
}

#model data using bayesian inference
# for (i in 1:nrow(species.table)) {
#   data <- as.data.frame(mget(paste(tolower(ascii.names[i]),".data",sep=""))[[paste(tolower(ascii.names[i]),".data",sep="")]]) 
#   coll <- data$COLL
#   occ <- data[,paste("c.log.",ascii.names[i],sep="")]
#   tvol <- data$c.log.TVOL
#   tspd <- data$c.log.TSPD
#   n <- nrow(data)
#   
#   model <- function() {
#     #priors
#     a ~ dnorm(0, 1.0E-6)
#     b[1] ~ dnorm(0, 1.0E-6)
#     b[2] ~ dnorm(0, 1.0E-6)
#     b[3] ~ dnorm(0, 1.0E-6)
#     
#     #likelihood
#     for (i in 1:n) {
#       coll[i] ~ dbern(p[i])
#       logit(p[i]) <- a + b[1] * occ[i] + b[2] * tvol[i] + b[3] * tspd[i] 
#     }
#   }
# 
# # bundle data for model
# jags.data <- list("coll", "occ", "tvol", "tspd", "n")
# 
# # initial values for model
# inits <- function() list(a=0, b=c(0,0,0))
# 
# # parameters to estimate
# parameters <- c("a", "b")
# 
# # start Gibbs sampling
# out <- jags.parallel(data = jags.data, inits = inits, parameters.to.save = parameters, model.file = model, n.chains=3, n.iter=10000)
# 
# assign(paste(tolower(ascii.names[i]),".bayes",sep=""), out)
# }


#model data with glm and summarize model output data
glm_sums <- data.frame(character(),character(),numeric(),numeric(),numeric(),numeric(),numeric(),stringsAsFactors=FALSE,row.names=NULL)
colnames(glm_sums) <- c("Species","Variable","Coefficient","Std. Error","$Z\\text{-value}$","$\\PRZ$","ANOVA")
for (i in 1:length(ascii.files)) {
  formula <- as.formula(paste("COLL ~ c.log.",ascii.names[i]," + c.log.TVOL + I(c.log.TVOL^2) + c.log.TSPD",sep=""))
  data <- as.data.frame(mget(paste(tolower(ascii.names[i]),".data",sep=""))[[paste(tolower(ascii.names[i]),".data",sep="")]])
  #model <- glmsimp(glm(formula = formula, family=binomial, data = data))
  #assign(paste(tolower(ascii.names[i]),"coll.glm",sep=""), model)
  #saveRDS(model,file=paste("/home/casey/Research/Projects/Risk_Model/Model/",tolower(ascii.names[i]),"coll.glm.rda",sep=""))
  #x <- sqrt(diag(vcov(glm(formula = formula, family=binomial, data = data))))
  #write.csv(x,paste("/home/casey/Research/Projects/Risk_Model/Model/",tolower(ascii.names[i]),"coll.glm.se",sep=""), row.names=FALSE)
  #model <- brglm(formula = formula, family = binomial(cloglog), data = data, method = "brglm.fit")
  model <- glm(formula = formula, family = binomial(link = "cloglog"), data = data)
  assign(paste(tolower(ascii.names[i]),"coll.glm",sep=""), model)
  assign(paste(tolower(ascii.names[i]),"coll.summary",sep=""),summary(model))
  assign(paste(tolower(ascii.names[i]),"coll.coef",sep=""),coef(model))
  x <- model
  assign(paste(tolower(ascii.names[i]),"coll.devexp",sep=""),round(((x$null.deviance - x$deviance)/x$null.deviance)*100,2))
  assign(paste(tolower(ascii.names[i]),"coll.roc",sep=""),roc(data$COLL,x$fitted.values))
  assign(paste(tolower(ascii.names[i]),"coll.resid",sep=""),resid(x))
  
  x.names <- c("Intercept", paste(ascii.names[i]), "TVOL", "TVOL$^2$", "TSPD")
  x.species <- c(paste(species.list[i],sep=""), NA, NA, NA, NA)
  x.coef <- signif(coef(summary(x))[,1],digits=4)
  x.se <- signif(coef(summary(x))[,2],digits=4)
  x.zvalue <- signif(coef(summary(x))[,3],digits=4)
  x.prz <- signif(coef(summary(x))[,4],digits=2)
  x.prz <- sapply(x.prz, function(x) ifelse(x < 2e-16, 2e-16, x))
  x.anova <- signif((anova(x)[,2]/sum(anova(x)[2:5,2]))*100,digits=4)
  x.anova[1] <- "---"
  x.all <- data.frame(cbind(x.species,x.names,x.coef,x.se,x.zvalue,x.prz,x.anova),stringsAsFactors=FALSE,row.names=NULL) 
  colnames(x.all) <- c("Species","Variable","Coefficient","Std. Error","$Z\\text{-value}$","$\\PRZ$","ANOVA")
  
  newrow = rep(NA,length(x.all))
  glm_sums <- rbind(glm_sums, x.all, newrow)
  
  rm(formula)
  rm(data)
  rm(model)
  rm(x)
  rm(x.names,x.coef,x.se,x.zvalue,x.prz,x.anova,x.species)
  rm(x.all)
  rm(newrow)
}

print(xtable(glm_sums), include.rownames=FALSE, sanitize.text.function=function(x){x}, floating=FALSE)

write.csv(glm_sums, file = "/home/casey/Research/Projects/Risk_Model/Data/glm_sums.csv", row.names=FALSE)


registerDoMC(detectCores() - 1)

glm.models <- foreach(i = 1:nrow(species.table)) %dopar% {
  formula <- as.formula(paste("COLL ~ c.log.",ascii.names[i]," + c.log.TVOL + I(c.log.TVOL^2) + c.log.TSPD",sep=""))
  data <- read.delim(paste("/home/casey/Research/Projects/Risk_Model/Data/",species.table[i,3],".data",sep=""), header=T, sep=",")
  data[, "c.log.TVOL"] <- log(data[,4])-mean(log(data[,4]))
  data[, "c.log.TSPD"] <- log(data[,5])-mean(log(data[,5]))
  data[, paste("c.log.",ascii.names[i],sep="")] <- log(data[,6])-mean(log(data[,6]))
  model <- glm(formula = formula, family=binomial(link = "cloglog"), data = data)
  model
}

#construct model performance table
coll_deviance <- data.frame("SPP"=rep(NA,nrow(species.table)),"DEV"=rep(NA,nrow(species.table)),"ROC"=rep(NA,nrow(species.table)))
for(i in 1:nrow(species.table)) {
  coll_deviance[i,"SPP"] <- toupper(paste(species.table[i,3]))
  coll_deviance[i,"DEV"] <- get(paste(species.table[i,3],"coll.devexp",sep=""))
  coll_deviance[i,"ROC"] <- get(paste(species.table[i,3],"coll.roc",sep=""))
}
write.csv(coll_deviance, file = "/home/casey/Research/Projects/Risk_Model/Data/coll_devs.csv", row.names=FALSE)

#make predictions based on models
#all.ind.data <- read.delim(paste("/home/casey/Research/Projects/Risk_Model/Data/all.ind.data",sep=""), header=T, sep=",")
glm.preds <- foreach(i = 1:nrow(species.table)) %dopar% {
  model <- glm.models[[i]]
  preds <- predict(model, all.ind.data, type="response")
  preds
}
coll_preds <- data.frame("UID"=all.ind.data$UID)
for(i in 1:nrow(species.table)) {
  x <- data.frame(glm.preds[[i]])
  rownames(x) <- NULL
  colnames(x) <- toupper(paste(species.table[i,3]))
  coll_preds <- cbind(coll_preds,x)
  rm(x)
}

#recode and classify road segment risk based on all species
coll_risk <- data.table(coll_preds)
setkey(coll_risk,UID,EGK,BTP,RTP,BSW,WOM,KOA)

coll_risk_tot <- coll_risk

coll_risk_tot[, TOTAL := rowSums(.SD>=0.90), .SDcols = c("EGK","BTP","RTP","BSW","WOM","KOA")]

coll_risk_tot[,.N,by=TOTAL]

coll_risk_95 <- coll_risk_tot[,.(UID,TOTAL)]

write.csv(coll_risk_95, file = "/home/casey/Research/Projects/Risk_Model/Data/coll_risk_95.csv", row.names=FALSE)

#plot effects
# occ <- NULL
# for (i in 1:length(ascii.files)) {
#   data <- as.data.frame(mget(paste(tolower(ascii.names[i]),".data",sep=""))[[paste(tolower(ascii.names[i]),".data",sep="")]])
#   model <- mget(paste(tolower(ascii.names[i]),"coll.glm",sep=""))[[paste(tolower(ascii.names[i]),"coll.glm",sep="")]]
#   coef <- mget(paste(tolower(ascii.names[i]),"coll.coef",sep=""))[[paste(tolower(ascii.names[i]),"coll.coef",sep="")]]
#   ml <- mget(paste("ml.",ascii.names[i],sep=""))[[paste("ml.",ascii.names[i],sep="")]]
#   temp_df <- data.frame(x=data[,6], y=invlogit(coef[1] + coef[2]*(log(data[,6])-ml)), col=rep(ascii.names[i], each=length(data[,6])))
#   occ <- rbind(occ,temp_df)
#   rm(data)
#   rm(model)
#   rm(coef)
#   rm(ml)
#   rm(temp_df)
# }
# tvol <- NULL
# for (i in 1:length(ascii.files)) {
#   data <- as.data.frame(mget(paste(tolower(ascii.names[i]),".data",sep=""))[[paste(tolower(ascii.names[i]),".data",sep="")]])
#   model <- mget(paste(tolower(ascii.names[i]),"coll.glm",sep=""))[[paste(tolower(ascii.names[i]),"coll.glm",sep="")]]
#   coef <- mget(paste(tolower(ascii.names[i]),"coll.coef",sep=""))[[paste(tolower(ascii.names[i]),"coll.coef",sep="")]]
#   ml <- mget(paste("ml.TVOL.",i,sep=""))[[paste("ml.TVOL.",i,sep="")]]
#   temp_df <- data.frame(x=data[,4], y=invlogit(coef[1] + coef[3]*(log(data[,4])-ml) + coef[4]*(log(data[,4])-ml)*(log(data[,4])-ml)), col=rep(ascii.names[i], each=length(data[,4])))
#   tvol <- rbind(tvol,temp_df)
#   rm(data)
#   rm(model)
#   rm(coef)
#   rm(ml)
#   rm(temp_df)
# }
# tspd <- NULL
# for (i in 1:length(ascii.files)) {
#   data <- as.data.frame(mget(paste(tolower(ascii.names[i]),".data",sep=""))[[paste(tolower(ascii.names[i]),".data",sep="")]])
#   model <- mget(paste(tolower(ascii.names[i]),"coll.glm",sep=""))[[paste(tolower(ascii.names[i]),"coll.glm",sep="")]]
#   coef <- mget(paste(tolower(ascii.names[i]),"coll.coef",sep=""))[[paste(tolower(ascii.names[i]),"coll.coef",sep="")]]
#   ml <- mget(paste("ml.TSPD.",i,sep=""))[[paste("ml.TSPD.",i,sep="")]]
#   temp_df <- data.frame(x=data[,5], y=invlogit(coef[1] + coef[5]*(log(data[,5])-ml)), col=rep(ascii.names[i], each=length(data[,5])))
#   tspd <- rbind(tspd,temp_df)
#   rm(data)
#   rm(model)
#   rm(coef)
#   rm(ml)
#   rm(temp_df)
# }

occ <- NULL
for (i in 1:length(ascii.files)) {
  data <- as.data.frame(mget(paste(tolower(ascii.names[i]),".data",sep=""))[[paste(tolower(ascii.names[i]),".data",sep="")]])
  model <- mget(paste(tolower(ascii.names[i]),"coll.glm",sep=""))[[paste(tolower(ascii.names[i]),"coll.glm",sep="")]]
  ml <- mget(paste("ml.",ascii.names[i],sep=""))[[paste("ml.",ascii.names[i],sep="")]]
  ml2 <- mget(paste("ml.TVOL.",i,sep=""))[[paste("ml.TVOL.",i,sep="")]]
  ml3 <- mget(paste("ml.TSPD.",i,sep=""))[[paste("ml.TSPD.",i,sep="")]]
  temp_df <- data.frame(x=data[,6], y=invlogit(cbind(1,log(data[,6])-ml,mean(log(data[,4])-ml2),mean((log(data[,4])-ml2)*(log(data[,4])-ml2)),mean(log(data[,5])-ml3)) %*% coef(model)), col=rep(ascii.names[i], each=length(data[,6])))
  occ <- rbind(occ,temp_df)
  rm(data)
  rm(model)
  rm(ml)
  rm(temp_df)
}  

ggplot(occ,aes(x=x,y=y,group=col,colour=factor(col))) +
  geom_line(size=1.5) +
  ylab("Likelihood of Collision\n") +
  xlab("\nLikelihood of Species Occurrence") +
  labs(color = "Species") +
  theme_bw() +
  theme(legend.key = element_blank()) +
  scale_colour_manual(values=plotPal) +
  theme(text = element_text(size = 20)) +
  scale_x_continuous(breaks=seq(0,1,by=.1), expand = c(0, 0), lim=c(0,1)) +
  scale_y_continuous(breaks=seq(0,1,by=.1), expand = c(0, 0), lim=c(0,1)) #+
  #guides(colour=FALSE)
 

tvol <- NULL
for (i in 1:length(ascii.files)) {
  data <- as.data.frame(mget(paste(tolower(ascii.names[i]),".data",sep=""))[[paste(tolower(ascii.names[i]),".data",sep="")]])
  model <- mget(paste(tolower(ascii.names[i]),"coll.glm",sep=""))[[paste(tolower(ascii.names[i]),"coll.glm",sep="")]]
  ml <- mget(paste("ml.",ascii.names[i],sep=""))[[paste("ml.",ascii.names[i],sep="")]]
  ml2 <- mget(paste("ml.TVOL.",i,sep=""))[[paste("ml.TVOL.",i,sep="")]]
  ml3 <- mget(paste("ml.TSPD.",i,sep=""))[[paste("ml.TSPD.",i,sep="")]]
  temp_df <- data.frame(x=data[,4], y=invlogit(cbind(1,mean(log(data[,6])-ml),log(data[,4])-ml2,(log(data[,4])-ml2)*(log(data[,4])-ml2),mean(log(data[,5])-ml3)) %*% coef(model)), col=rep(ascii.names[i], each=length(data[,6])))
  tvol <- rbind(tvol,temp_df)
  rm(data)
  rm(model)
  rm(ml)
  rm(temp_df)
}

ggplot(tvol,aes(x=x/1000,y=y,group=col,colour=factor(col))) +
  geom_line(size=1.5) +
  ylab("Likelihood of Collision\n") +
  xlab("\nTraffic Volume (1000 vehicles/day)") +
  labs(color = "Species") +
  theme_bw() +
  theme(legend.key = element_blank()) +
  scale_colour_manual(values=plotPal) +
  theme(text = element_text(size = 20)) +
  scale_x_continuous(breaks=seq(0,35,by=10), expand = c(0, 0), lim=c(0,35)) +
  scale_y_continuous(breaks=seq(0,1,by=.1), expand = c(0, 0), lim=c(0,1)) #+
  #guides(colour=FALSE)

tspd <- NULL
for (i in 1:length(ascii.files)) {
  data <- as.data.frame(mget(paste(tolower(ascii.names[i]),".data",sep=""))[[paste(tolower(ascii.names[i]),".data",sep="")]])
  model <- mget(paste(tolower(ascii.names[i]),"coll.glm",sep=""))[[paste(tolower(ascii.names[i]),"coll.glm",sep="")]]
  ml <- mget(paste("ml.",ascii.names[i],sep=""))[[paste("ml.",ascii.names[i],sep="")]]
  ml2 <- mget(paste("ml.TVOL.",i,sep=""))[[paste("ml.TVOL.",i,sep="")]]
  ml3 <- mget(paste("ml.TSPD.",i,sep=""))[[paste("ml.TSPD.",i,sep="")]]
  temp_df <- data.frame(x=data[,5], y=invlogit(cbind(1,mean(log(data[,6])-ml),mean(log(data[,4])-ml2),mean((log(data[,4])-ml2)*(log(data[,4])-ml2)),log(data[,5])-ml3) %*% coef(model)), col=rep(ascii.names[i], each=length(data[,6])))
  tspd <- rbind(tspd,temp_df)
  rm(data)
  rm(model)
  rm(ml)
  rm(temp_df)
}  

ggplot(tspd,aes(x=x,y=y,group=col,colour=factor(col))) +
  geom_line(size=1.5) +
  ylab("Likelihood of Collision\n") +
  xlab("\nTraffic Speed (km/hour)") +
  labs(color = "Species") +
  theme_bw() +
  theme(legend.key = element_blank()) +
  scale_colour_manual(values=plotPal) +
  theme(text = element_text(size = 20)) +
  scale_x_continuous(breaks=seq(30,110,by=10), expand = c(0, 0), lim=c(30,110)) +
  scale_y_continuous(breaks=seq(0,1,by=.1), expand = c(0, 0), lim=c(0,1)) #+
  #guides(colour=FALSE)


#calculate spatial autocorrelation across all species
registerDoMC(detectCores() - 1)

auto.resid <- foreach(i = 1:nrow(species.table), .packages = c("ncf"), .combine="rbind") %dopar% {
  data <- read.delim(paste("/home/casey/Research/Projects/Risk_Model/Data/",species.table[i,3],".data",sep=""), header=T, sep=",")
  data[, "c.log.TVOL"] <- log(data[,4])-mean(log(data[,4]))
  data[, "c.log.TSPD"] <- log(data[,5])-mean(log(data[,5]))
  data[, paste("c.log.",ascii.names[i],sep="")] <- log(data[,6])-mean(log(data[,6]))
  model <- glm.models[[i]]
  cor <- correlog(data[,1], data[,2], resid(model), increment=1000, resamp=0, latlon=FALSE)
  temp_df <- data.frame(x=as.numeric(names(cor$correlation[2:21])), y=cor$correlation[2:21], col=rep(toupper(species.table[i,3]), each=length(cor$correlation[2:21])))
  temp_df
}  

auto.coll <- foreach(i = 1:nrow(species.table), .packages = c("ncf"), .combine="rbind") %dopar% {
  data <- read.delim(paste("/home/casey/Research/Projects/Risk_Model/Data/",species.table[i,3],".data",sep=""), header=T, sep=",")
  cor <- correlog(data[,1], data[,2], data[,3], increment=1000, resamp=0, latlon=FALSE)
  temp_df <- data.frame(x=as.numeric(names(cor$correlation[2:21])), y=cor$correlation[2:21], col=rep(toupper(species.table[i,3]), each=length(cor$correlation[2:21])))
  temp_df
}  

shapes <- unlist(lapply(c("1", "2", "3", "4", "5", "6"), utf8ToInt))

ggplot(auto.resid,aes(x=x,y=y,group=col,shape=col)) + geom_line(colour=c("grey70"),size=.75) + geom_point(size=2.5) + ylab("Moran's I") + xlab("Distance (km)") + labs(shape = "Species") + theme_bw() + theme(legend.key = element_blank()) + scale_colour_manual(values=plotPal) + scale_shape_manual(values=shapes) + geom_hline(aes(yintercept=0), linetype=2) + scale_x_continuous(breaks=seq(1, 20, 1))

ggplot(auto.coll,aes(x=x,y=y,group=col,shape=col)) + geom_line(colour=c("grey70"),size=.75) + geom_point(size=2.5) + ylab("Moran's I") + xlab("Distance (km)") + labs(shape = "Species") + theme_bw() + theme(legend.key = element_blank()) + scale_colour_manual(values=plotPal) + scale_shape_manual(values=shapes) + geom_hline(aes(yintercept=0), linetype=2) + scale_x_continuous(breaks=seq(1, 20, 1))

rd.data <- read.delim("/home/casey/Research/Projects/Roads/Data/Pred/all_preds.csv", header=T, sep=",")


d.m <- dist(as.matrix(cbind(bsw.data$X,bsw.data$Y)), method = "euclidean", diag = FALSE, upper = FALSE)
require(maptools)
require(raster)
require(xtable)
require(R2jags)
require(ncf)
require(doMC)
require(data.table)

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

species.table <- read.delim("data/species_list.csv", header=T, sep=",")

load("data/coll_model_data")

registerDoMC(detectCores() - 1)
coll.glm <- foreach(i = 1:nrow(species.table)) %dopar% {
  formula <- as.formula(paste0("coll ~ log(",species.table[i,2],") + log(tvol) + I(log(tvol)^2) + log(tspd)"))
  model <- glm(formula = formula, family=binomial(link = "cloglog"), data = model.data[[i]])
  summary(model)
}

registerDoMC(detectCores() - 1)
coll.glm.deviance <- foreach(i = 1:nrow(species.table)) %dopar% {
  formula <- as.formula(paste0("coll ~ log(",species.table[i,2],") + log(tvol) + I(log(tvol)^2) + log(tspd)"))
  model <- glm(formula = formula, family=binomial(link = "cloglog"), data = model.data[[i]])
  paste("% Deviance Explained ",species.table[i,2],": ",round(((model$null.deviance - model$deviance)/model$null.deviance)*100,2),sep="")
}

##############
paste("% Deviance Explained: ",round(((coll.glm$null.deviance - coll.glm$deviance)/coll.glm$null.deviance)*100,2),sep="")  #Report reduction in deviance

write.csv(signif(summary(coll.glm)$coefficients, digits=4),"output/coll_coef.csv",row.names=FALSE)

write.csv(formatC(anova(coll.glm)[2:4,2]/sum(anova(coll.glm)[2:4,2]), format='f',digits=4),"output/coll_anova.csv",row.names=FALSE)

save(coll.glm,file="output/coll_glm")

save(model.data,file="output/coll_model_data")

coll.preds <- predict(coll.glm, cov.data, type="response")

coll.preds.df <- cbind("uid"=cov.data$uid,"collrisk"=coll.preds) #Combine predictions with unique IDs for all road segments
coll.preds.df <- na.omit(coll.preds.df)

write.csv(coll.preds.df, file = "output/coll_preds_glm.csv", row.names=FALSE)
###############

#model data with glm and summarize model output data (for word table)
glm_sums <- data.frame(character(),character(),numeric(),numeric(),numeric(),numeric(),numeric(),stringsAsFactors=FALSE,row.names=NULL)
colnames(glm_sums) <- c("Species","Variable","Coefficient","Std. Error","$Z\\text{-value}$","$\\PRZ$","ANOVA")
for (i in 1:nrow(species.table)) {
  formula <- as.formula(paste0("coll ~ log(",species.table[i,2],") + log(tvol) + I(log(tvol)^2) + log(tspd)"))
  data <- model.data[[i]]
  model <- glm(formula = formula, family = binomial(link = "cloglog"), data = data)
  assign(paste(tolower(species.table[i,2]),"coll.glm",sep=""), model)
  assign(paste(tolower(species.table[i,2]),"coll.summary",sep=""),summary(model))
  assign(paste(tolower(species.table[i,2]),"coll.coef",sep=""),coef(model))
  x <- model
  assign(paste(tolower(species.table[i,2]),"coll.devexp",sep=""),round(((x$null.deviance - x$deviance)/x$null.deviance)*100,2))
  assign(paste(tolower(species.table[i,2]),"coll.roc",sep=""),roc(data$coll,x$fitted.values))
  assign(paste(tolower(species.table[i,2]),"coll.resid",sep=""),resid(x))
  
  x.names <- c("Intercept", toupper(paste(species.table[i,2])), "TVOL", "TVOL2", "TSPD")
  x.species <- c(paste(species.list[i],sep=""), NA, NA, NA, NA)
  x.coef <- signif(coef(summary(x))[,1],digits=4)
  x.se <- signif(coef(summary(x))[,2],digits=4)
  x.zvalue <- signif(coef(summary(x))[,3],digits=4)
  x.prz <- signif(coef(summary(x))[,4],digits=2)
  x.prz <- sapply(x.prz, function(x) ifelse(x < 2e-16, 2e-16, x))
  x.anova <- signif((anova(x)[,2]/sum(anova(x)[2:5,2]))*100,digits=4)
  x.anova[1] <- "---"
  x.all <- data.frame(cbind(x.species,x.names,x.coef,x.se,x.zvalue,x.prz,x.anova),stringsAsFactors=FALSE,row.names=NULL) 
  colnames(x.all) <- c("Species","Variable","Coefficient","Std. Error","Z-value","Pr(Z)","ANOVA")
  
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

#model data with glm and summarize model output data (for latex table)
# glm_sums <- data.frame(character(),character(),numeric(),numeric(),numeric(),numeric(),numeric(),stringsAsFactors=FALSE,row.names=NULL)
# colnames(glm_sums) <- c("Species","Variable","Coefficient","Std. Error","$Z\\text{-value}$","$\\PRZ$","ANOVA")
# for (i in 1:nrow(species.table)) {
#   formula <- as.formula(paste0("coll ~ log(",species.table[i,2],") + log(tvol) + I(log(tvol)^2) + log(tspd)"))
#   data <- model.data[[i]]
#   model <- glm(formula = formula, family = binomial(link = "cloglog"), data = data)
#   assign(paste(tolower(species.table[i,2]),"coll.glm",sep=""), model)
#   assign(paste(tolower(species.table[i,2]),"coll.summary",sep=""),summary(model))
#   assign(paste(tolower(species.table[i,2]),"coll.coef",sep=""),coef(model))
#   x <- model
#   assign(paste(tolower(species.table[i,2]),"coll.devexp",sep=""),round(((x$null.deviance - x$deviance)/x$null.deviance)*100,2))
#   assign(paste(tolower(species.table[i,2]),"coll.roc",sep=""),roc(data$coll,x$fitted.values))
#   assign(paste(tolower(species.table[i,2]),"coll.resid",sep=""),resid(x))
#   
#   x.names <- c("Intercept", toupper(paste(species.table[i,2])), "TVOL", "TVOL$^2$", "TSPD")
#   x.species <- c(paste(species.list[i],sep=""), NA, NA, NA, NA)
#   x.coef <- signif(coef(summary(x))[,1],digits=4)
#   x.se <- signif(coef(summary(x))[,2],digits=4)
#   x.zvalue <- signif(coef(summary(x))[,3],digits=4)
#   x.prz <- signif(coef(summary(x))[,4],digits=2)
#   x.prz <- sapply(x.prz, function(x) ifelse(x < 2e-16, 2e-16, x))
#   x.anova <- signif((anova(x)[,2]/sum(anova(x)[2:5,2]))*100,digits=4)
#   x.anova[1] <- "---"
#   x.all <- data.frame(cbind(x.species,x.names,x.coef,x.se,x.zvalue,x.prz,x.anova),stringsAsFactors=FALSE,row.names=NULL) 
#   colnames(x.all) <- c("Species","Variable","Coefficient","Std. Error","$Z\\text{-value}$","$\\PRZ$","ANOVA")
#   
#   newrow = rep(NA,length(x.all))
#   glm_sums <- rbind(glm_sums, x.all, newrow)
#   
#   rm(formula)
#   rm(data)
#   rm(model)
#   rm(x)
#   rm(x.names,x.coef,x.se,x.zvalue,x.prz,x.anova,x.species)
#   rm(x.all)
#   rm(newrow)
# }
# print(xtable(glm_sums), include.rownames=FALSE, sanitize.text.function=function(x){x}, floating=FALSE)

write.csv(glm_sums, file = "output/glm_sums.csv", row.names=FALSE)


#construct model performance table
coll_deviance <- data.frame("SPP"=rep(NA,nrow(species.table)),"DEV"=rep(NA,nrow(species.table)),"ROC"=rep(NA,nrow(species.table)))
for(i in 1:nrow(species.table)) {
  coll_deviance[i,"SPP"] <- toupper(paste(species.table[i,2]))
  coll_deviance[i,"DEV"] <- get(paste(species.table[i,2],"coll.devexp",sep=""))
  coll_deviance[i,"ROC"] <- get(paste(species.table[i,2],"coll.roc",sep=""))
}
write.csv(coll_deviance, file = "output/coll_devs.csv", row.names=FALSE)

#make predictions based on models
glm.preds <- foreach(i = 1:nrow(species.table)) %dopar% {
  formula <- as.formula(paste0("coll ~ log(",species.table[i,2],") + log(tvol) + I(log(tvol)^2) + log(tspd)"))
  data <- model.data[[i]]
  model <- glm(formula = formula, family = binomial(link = "cloglog"), data = data)
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
# for (i in 1:nrow(species.table)) {
#   data <- as.data.frame(mget(paste(tolower(species.table[i,2]),".data",sep=""))[[paste(tolower(species.table[i,2]),".data",sep="")]])
#   model <- mget(paste(tolower(species.table[i,2]),"coll.glm",sep=""))[[paste(tolower(species.table[i,2]),"coll.glm",sep="")]]
#   coef <- mget(paste(tolower(species.table[i,2]),"coll.coef",sep=""))[[paste(tolower(species.table[i,2]),"coll.coef",sep="")]]
#   ml <- mget(paste("ml.",species.table[i,2],sep=""))[[paste("ml.",species.table[i,2],sep="")]]
#   temp_df <- data.frame(x=data[,6], y=invlogit(coef[1] + coef[2]*(log(data[,6])-ml)), col=rep(species.table[i,2], each=length(data[,6])))
#   occ <- rbind(occ,temp_df)
#   rm(data)
#   rm(model)
#   rm(coef)
#   rm(ml)
#   rm(temp_df)
# }
# tvol <- NULL
# for (i in 1:nrow(species.table)) {
#   data <- as.data.frame(mget(paste(tolower(species.table[i,2]),".data",sep=""))[[paste(tolower(species.table[i,2]),".data",sep="")]])
#   model <- mget(paste(tolower(species.table[i,2]),"coll.glm",sep=""))[[paste(tolower(species.table[i,2]),"coll.glm",sep="")]]
#   coef <- mget(paste(tolower(species.table[i,2]),"coll.coef",sep=""))[[paste(tolower(species.table[i,2]),"coll.coef",sep="")]]
#   ml <- mget(paste("ml.TVOL.",i,sep=""))[[paste("ml.TVOL.",i,sep="")]]
#   temp_df <- data.frame(x=data[,4], y=invlogit(coef[1] + coef[3]*(log(data[,4])-ml) + coef[4]*(log(data[,4])-ml)*(log(data[,4])-ml)), col=rep(species.table[i,2], each=length(data[,4])))
#   tvol <- rbind(tvol,temp_df)
#   rm(data)
#   rm(model)
#   rm(coef)
#   rm(ml)
#   rm(temp_df)
# }
# tspd <- NULL
# for (i in 1:nrow(species.table)) {
#   data <- as.data.frame(mget(paste(tolower(species.table[i,2]),".data",sep=""))[[paste(tolower(species.table[i,2]),".data",sep="")]])
#   model <- mget(paste(tolower(species.table[i,2]),"coll.glm",sep=""))[[paste(tolower(species.table[i,2]),"coll.glm",sep="")]]
#   coef <- mget(paste(tolower(species.table[i,2]),"coll.coef",sep=""))[[paste(tolower(species.table[i,2]),"coll.coef",sep="")]]
#   ml <- mget(paste("ml.TSPD.",i,sep=""))[[paste("ml.TSPD.",i,sep="")]]
#   temp_df <- data.frame(x=data[,5], y=invlogit(coef[1] + coef[5]*(log(data[,5])-ml)), col=rep(species.table[i,2], each=length(data[,5])))
#   tspd <- rbind(tspd,temp_df)
#   rm(data)
#   rm(model)
#   rm(coef)
#   rm(ml)
#   rm(temp_df)
# }

occ <- NULL
for (i in 1:nrow(species.table)) {
  data <- as.data.frame(mget(paste(tolower(species.table[i,2]),".data",sep=""))[[paste(tolower(species.table[i,2]),".data",sep="")]])
  model <- mget(paste(tolower(species.table[i,2]),"coll.glm",sep=""))[[paste(tolower(species.table[i,2]),"coll.glm",sep="")]]
  ml <- mget(paste("ml.",species.table[i,2],sep=""))[[paste("ml.",species.table[i,2],sep="")]]
  ml2 <- mget(paste("ml.TVOL.",i,sep=""))[[paste("ml.TVOL.",i,sep="")]]
  ml3 <- mget(paste("ml.TSPD.",i,sep=""))[[paste("ml.TSPD.",i,sep="")]]
  temp_df <- data.frame(x=data[,6], y=invlogit(cbind(1,log(data[,6])-ml,mean(log(data[,4])-ml2),mean((log(data[,4])-ml2)*(log(data[,4])-ml2)),mean(log(data[,5])-ml3)) %*% coef(model)), col=rep(species.table[i,2], each=length(data[,6])))
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
for (i in 1:nrow(species.table)) {
  data <- as.data.frame(mget(paste(tolower(species.table[i,2]),".data",sep=""))[[paste(tolower(species.table[i,2]),".data",sep="")]])
  model <- mget(paste(tolower(species.table[i,2]),"coll.glm",sep=""))[[paste(tolower(species.table[i,2]),"coll.glm",sep="")]]
  ml <- mget(paste("ml.",species.table[i,2],sep=""))[[paste("ml.",species.table[i,2],sep="")]]
  ml2 <- mget(paste("ml.TVOL.",i,sep=""))[[paste("ml.TVOL.",i,sep="")]]
  ml3 <- mget(paste("ml.TSPD.",i,sep=""))[[paste("ml.TSPD.",i,sep="")]]
  temp_df <- data.frame(x=data[,4], y=invlogit(cbind(1,mean(log(data[,6])-ml),log(data[,4])-ml2,(log(data[,4])-ml2)*(log(data[,4])-ml2),mean(log(data[,5])-ml3)) %*% coef(model)), col=rep(species.table[i,2], each=length(data[,6])))
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
for (i in 1:nrow(species.table)) {
  data <- as.data.frame(mget(paste(tolower(species.table[i,2]),".data",sep=""))[[paste(tolower(species.table[i,2]),".data",sep="")]])
  model <- mget(paste(tolower(species.table[i,2]),"coll.glm",sep=""))[[paste(tolower(species.table[i,2]),"coll.glm",sep="")]]
  ml <- mget(paste("ml.",species.table[i,2],sep=""))[[paste("ml.",species.table[i,2],sep="")]]
  ml2 <- mget(paste("ml.TVOL.",i,sep=""))[[paste("ml.TVOL.",i,sep="")]]
  ml3 <- mget(paste("ml.TSPD.",i,sep=""))[[paste("ml.TSPD.",i,sep="")]]
  temp_df <- data.frame(x=data[,5], y=invlogit(cbind(1,mean(log(data[,6])-ml),mean(log(data[,4])-ml2),mean((log(data[,4])-ml2)*(log(data[,4])-ml2)),log(data[,5])-ml3) %*% coef(model)), col=rep(species.table[i,2], each=length(data[,6])))
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
  data[, paste("c.log.",species.table[i,2],sep="")] <- log(data[,6])-mean(log(data[,6]))
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
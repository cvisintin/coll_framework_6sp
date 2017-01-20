require(maptools)
require(raster)
require(xtable)
require(R2jags)
require(ncf)
require(doMC)
require(data.table)
require(logistf)
require(RPostgreSQL)
require(rethinking)
require(rstan)

drv <- dbDriver("PostgreSQL")  #Specify a driver for postgreSQL type database
con <- dbConnect(drv, dbname="qaeco_spatial", user="qaeco", password="Qpostgres15", host="boab.qaeco.com", port="5432")  #Connection to database server on Boab


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

species.table <- read.delim("data/species_list.csv", header=T, sep=",")
species.list <- c("Eastern Grey Kangaroo","Common Brushtail Possum","Common Ringtail Possum","Black Swamp Wallaby","Common Wombat","Koala")

load("data/coll_model_data")

#model.data2 <- lapply(model.data, function(x) x[coll==0,AC:=0])

registerDoMC(detectCores() - 1)
coll.glm <- foreach(i = 1:nrow(species.table)) %dopar% {
  formula <- as.formula(paste0("coll ~ log(",species.table[i,2],") + log(tvol) + I(log(tvol)^2) + log(tspd) + AC"))
  model <- glm(formula = formula, family=binomial(link = "cloglog"), data = model.data[[i]])
}
save(coll.glm, file="output/coll_glm")

registerDoMC(detectCores() - 1)
coll.glm.deviance <- foreach(i = 1:nrow(species.table)) %dopar% {
  formula <- as.formula(paste0("coll ~ log(",species.table[i,2],") + log(tvol) + I(log(tvol)^2) + log(tspd) + AC"))
  model <- glm(formula = formula, family=binomial(link = "cloglog"), data = model.data[[i]])
  paste("% Deviance Explained ",species.table[i,2],": ",round(((model$null.deviance - model$deviance)/model$null.deviance)*100,2),sep="")
}

#combine all datasets for glm model
model.data2 <- model.data

for(i in 1:nrow(species.table)){
  colnames(model.data2[[i]]) <- c("uid","x","y","occ","tvol","tspd","coll")
}

model.data.all <- do.call(rbind,model.data2)

model <- glm(formula = coll ~ log(occ) + log(tvol) + I(log(tvol)^2) + log(tspd), family=binomial(link = "cloglog"), data = model.data.all[!duplicated(model.data.all$x),])
paste("% Deviance Explained: ",round(((model$null.deviance - model$deviance)/model$null.deviance)*100,2),sep="")

summary(model)

##############
# set.seed(123) 
# coll.model.map <- map(
#   alist(
#     y ~ dbinom(1,p), 
#     p <- 1 - exp(-exp(a + b1*x1 + b2*x2 + b3*x3 + b4*x4 + b5*x5)),
#     a ~ dnorm(0,1),
#     b1 ~ dnorm(0,1),
#     b2 ~ dnorm(0,1),
#     b3 ~ dnorm(0,1),
#     b4 ~ dnorm(0,1),
#     b5 ~ dnorm(0,100)
#   ),
#   start=list(a=0.5,b1=0.5,b2=0.5,b3=-0.5,b4=0.5,b5=0.5),
#   data=list(y=model.data[[i]][[7]],
#             x1=log(model.data[[i]][[4]]),
#             x2=log(model.data[[i]][[5]]/1000),
#             x3=log((model.data[[i]][[5]]/1000)*(model.data[[i]][[5]]/1000)),
#             x4=log(model.data[[i]][[6]]),
#             x5=model.data[[i]][[8]]
#             ),
#   #iter=1000,
#   #warmup=500,
#   #chains=3,
#   #cores=3
# )
# precis(coll.model.map)
###############

#model data with glm and summarize model output data (for word table)
glm_sums <- data.frame(character(),character(),numeric(),numeric(),numeric(),numeric(),numeric(),stringsAsFactors=FALSE,row.names=NULL)
colnames(glm_sums) <- c("Species","Variable","Coefficient","Std. Error","$Z\\text{-value}$","$\\PRZ$","ANOVA")
for (i in 1:nrow(species.table)) {
  formula <- as.formula(paste0("coll ~ log(",species.table[i,2],") + log(tvol) + I(log(tvol)^2) + log(tspd)"))
  data <- model.data[[i]]
  model <- glm(formula = formula, family = binomial(link = "cloglog"), data = data)
  #assign(paste(tolower(species.table[i,2]),"coll.glm",sep=""), model)
  #assign(paste(tolower(species.table[i,2]),"coll.summary",sep=""),summary(model))
  #assign(paste(tolower(species.table[i,2]),"coll.coef",sep=""),coef(model))
  x <- model
  #assign(paste(tolower(species.table[i,2]),"coll.devexp",sep=""),round(((x$null.deviance - x$deviance)/x$null.deviance)*100,2))
  #assign(paste(tolower(species.table[i,2]),"coll.roc",sep=""),roc(data$coll,x$fitted.values))
  #assign(paste(tolower(species.table[i,2]),"coll.resid",sep=""),resid(x))
  
  x.names <- c("Intercept", toupper(paste(species.table[i,2])), "TVOL", "TVOL2", "TSPD")
  x.species <- c(paste(species.list[i],sep=""), NA, NA, NA, NA)
  x.coef <- signif(coef(summary(x))[,1],digits=4)
  x.se <- signif(coef(summary(x))[,2],digits=4)
  x.zvalue <- signif(coef(summary(x))[,3],digits=4)
  x.prz <- signif(coef(summary(x))[,4],digits=2)
  x.prz <- sapply(x.prz, function(x) ifelse(x < 2e-16, "<2e-16", x))
  x.anova <- signif((anova(x)[,2]/sum(anova(x)[2:6,2]))*100,digits=4)
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


#Validate model fit with independent data
load("data/coll_val_data")

registerDoMC(detectCores() - 1)
coll.val.roc <- foreach(i = 1:nrow(species.table)) %dopar% {
  formula <- as.formula(paste0("coll ~ log(",species.table[i,2],") + log(tvol) + I(log(tvol)^2) + log(tspd)"))
  data <- as.data.frame(model.data[[i]])
  model <- glm(formula = formula, family=binomial(link = "cloglog"), data = data)
  val.pred.glm <- predict(model, val.data[[i]], type="response")  #Make predictions with regression model fit
  roc(val.data[[i]]$coll, val.pred.glm)  #Compare collision records to predictions using receiver operator characteristic (ROC) function and report value
}
write.csv(coll.val.roc, file = "output/coll_ind_roc.csv", row.names=FALSE)


#combine all validation data to assess all species model
val.data2 <- val.data

for(i in 1:nrow(species.table)){
  colnames(val.data2[[i]]) <- c("uid","x","y","occ","tvol","tspd","coll")
}

val.data.all <- do.call(rbind,val.data2)

val.pred.glm <- predict(model, val.data.all[!duplicated(val.data.all$x),], type="response")  #Make predictions with regression model fit
roc(val.data.all[!duplicated(val.data.all$x),coll], val.pred.glm)  #Compare collision records to predictions using receiver operator characteristic (ROC) function and report value


#make predictions based on models
load("data/cov_data")

registerDoMC(detectCores() - 1)
glm.preds <- foreach(i = 1:nrow(species.table)) %dopar% {
  formula <- as.formula(paste0("coll ~ log(",species.table[i,2],") + log(tvol) + I(log(tvol)^2) + log(tspd)"))
  data <- model.data[[i]]
  model <- glm(formula = formula, family = binomial(link = "cloglog"), data = data)
  preds <- predict(model, cov.data, type="response")
  preds
}

coll_preds <- data.frame("uid"=cov.data$uid)
for(i in 1:nrow(species.table)) {
  x <- data.frame(glm.preds[[i]])
  rownames(x) <- NULL
  colnames(x) <- paste(species.table[i,2])
  coll_preds <- cbind(coll_preds,x)
  rm(x)
}

dbWriteTable(con, c("gis_victoria", "vic_nogeom_roads_6spcollrisk"), value = coll_preds, row.names=FALSE)

#recode and classify road segment risk based on all species
coll_risk <- na.omit(data.table(coll_preds))
setkey(coll_risk,UID,EGK,BTP,RTP,BSW,WOM,KOA)

coll_risk_tot <- coll_risk

threshold <- 0.75

coll_risk_tot[, TOTAL := rowSums(.SD>=threshold), .SDcols = c("EGK","BTP","RTP","BSW","WOM","KOA")]

coll_risk_tot[,.N,by=TOTAL]

assign(paste0("coll_risk_",threshold*100), coll_risk_tot[,.(UID,TOTAL)])

write.csv(get(paste0("coll_risk_",threshold*100)), file = paste0("output/coll_risk_",threshold*100,".csv"), row.names=FALSE)


#model based on seasons
load("data/coll_model_data_sum")
load("data/coll_model_data_aut")
load("data/coll_model_data_win")
load("data/coll_model_data_spr")

# registerDoMC(detectCores() - 1)
# coll.glm.deviance.summer <- foreach(i = 1:nrow(species.table)) %dopar% {
#   formula <- as.formula(paste0("coll ~ log(",species.table[i,2],") + log(tvol) + I(log(tvol)^2) + log(tspd)"))
#   model <- logistf(formula = formula, family=binomial(link = "cloglog"), data = model.data.summer[[i]])
#   paste("% Deviance Explained ",species.table[i,2],": ",round(((model$null.deviance - model$deviance)/model$null.deviance)*100,2),sep="")
# }

registerDoMC(detectCores() - 1)
coll.glm.summer <- foreach(i = 1:nrow(species.table)) %dopar% {
  formula <- as.formula(paste0("coll ~ log(",species.table[i,2],") + log(tvol) + I(log(tvol)^2) + log(tspd)"))
  model <- logistf(formula = formula, family=binomial(link = "cloglog"), data = model.data.summer[[i]])
}
save(coll.glm.summer, file="output/coll_glm_summer")


# registerDoMC(detectCores() - 1)
# coll.glm.deviance.autumn <- foreach(i = 1:nrow(species.table)) %dopar% {
#   formula <- as.formula(paste0("coll ~ log(",species.table[i,2],") + log(tvol) + I(log(tvol)^2) + log(tspd)"))
#   model <- logistf(formula = formula, family=binomial(link = "cloglog"), data = model.data.autumn[[i]])
#   paste("% Deviance Explained ",species.table[i,2],": ",round(((model$null.deviance - model$deviance)/model$null.deviance)*100,2),sep="")
# }

registerDoMC(detectCores() - 1)
coll.glm.autumn <- foreach(i = 1:nrow(species.table)) %dopar% {
  formula <- as.formula(paste0("coll ~ log(",species.table[i,2],") + log(tvol) + I(log(tvol)^2) + log(tspd)"))
  model <- logistf(formula = formula, family=binomial(link = "cloglog"), data = model.data.summer[[i]])
}
save(coll.glm.autumn, file="output/coll_glm_autumn")


# registerDoMC(detectCores() - 1)
# coll.glm.deviance.winter <- foreach(i = 1:nrow(species.table)) %dopar% {
#   formula <- as.formula(paste0("coll ~ log(",species.table[i,2],") + log(tvol) + I(log(tvol)^2) + log(tspd)"))
#   model <- logistf(formula = formula, family=binomial(link = "cloglog"), data = model.data.winter[[i]])
#   paste("% Deviance Explained ",species.table[i,2],": ",round(((model$null.deviance - model$deviance)/model$null.deviance)*100,2),sep="")
# }

registerDoMC(detectCores() - 1)
coll.glm.winter <- foreach(i = 1:nrow(species.table)) %dopar% {
  formula <- as.formula(paste0("coll ~ log(",species.table[i,2],") + log(tvol) + I(log(tvol)^2) + log(tspd)"))
  model <- logistf(formula = formula, family=binomial(link = "cloglog"), data = model.data.summer[[i]])
}
save(coll.glm.winter, file="output/coll_glm_winter")


# registerDoMC(detectCores() - 1)
# coll.glm.deviance.spring <- foreach(i = 1:nrow(species.table)) %dopar% {
#   formula <- as.formula(paste0("coll ~ log(",species.table[i,2],") + log(tvol) + I(log(tvol)^2) + log(tspd)"))
#   model <- logistf(formula = formula, family=binomial(link = "cloglog"), data = model.data.spring[[i]])
#   paste("% Deviance Explained ",species.table[i,2],": ",round(((model$null.deviance - model$deviance)/model$null.deviance)*100,2),sep="")
# }

registerDoMC(detectCores() - 1)
coll.glm.spring <- foreach(i = 1:nrow(species.table)) %dopar% {
  formula <- as.formula(paste0("coll ~ log(",species.table[i,2],") + log(tvol) + I(log(tvol)^2) + log(tspd)"))
  model <- logistf(formula = formula, family=binomial(link = "cloglog"), data = model.data.summer[[i]])
}
save(coll.glm.spring, file="output/coll_glm_spring")

###########Spatial Generalized Linear Mixed Model###########

group <- factor(rep("a",nrow(model.data[[1]])))
model.data <- cbind(model.data[[1]], group)
attach(model.data)

data.coords <- data.frame("x"=model.data[,x], "y"=model.data[,y])
coordinates(data.coords) <- ~x+y
d <- gDistance(data.coords, byid=T)
min.d <- apply(d, 1, function(x) sort(x[x>0], decreasing=F)[1])


model.e <-glmmPQL(coll ~ log(egk) + log(tvol) + I(log(tvol)^2) + log(tspd), random=~1|group, data=model.data, correlation=corExp( form=~x+y), family=binomial(link = "cloglog"))

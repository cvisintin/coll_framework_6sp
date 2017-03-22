#require(maptools)
require(raster)
require(xtable)
#require(R2jags)
#require(ncf)
require(doMC)
require(data.table)
#require(logistf)
require(brglm)
require(RPostgreSQL)
require(PresenceAbsence)
#require(rethinking)
#require(rstan)

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
species.list <- c("Eastern Grey Kangaroo","Common Brushtail Possum","Common Ringtail Possum","Swamp Wallaby","Common Wombat","Koala")

load("data/coll_model_data")

#model.data2 <- lapply(model.data, function(x) x[coll==0,AC:=0])

registerDoMC(detectCores() - 1)
coll.glm <- foreach(i = 1:nrow(species.table)) %dopar% {
  formula <- as.formula(paste0("coll ~ log(",species.table[i,2],") + log(tvol) + I(log(tvol)^2) + log(tspd)"))# + AC"))
  model <- glm(formula = formula, offset=log(length*4), family=binomial(link = "cloglog"), data = model.data[[i]])
}
save(coll.glm, file="output/coll_glm")

registerDoMC(detectCores() - 1)
coll.glm.deviance <- foreach(i = 1:nrow(species.table)) %dopar% {
  formula <- as.formula(paste0("coll ~ log(",species.table[i,2],") + log(tvol) + I(log(tvol)^2) + log(tspd)"))# + AC"))
  model <- glm(formula = formula, offset=log(length*4), family=binomial(link = "cloglog"), data = model.data[[i]])
  paste("% Deviance Explained ",species.table[i,2],": ",round(((model$null.deviance - model$deviance)/model$null.deviance)*100,2),sep="")
}

#combine all datasets for glm model
# model.data2 <- model.data
# 
# for(i in 1:nrow(species.table)){
#   colnames(model.data2[[i]]) <- c("uid","length","x","y","occ","tvol","tspd","coll")
# }
# 
# model.data.all <- do.call(rbind,model.data2)
# 
# model <- glm(formula = coll ~ log(occ) + log(tvol) + I(log(tvol)^2) + log(tspd), offset=log(length*4), family=binomial(link = "cloglog"), data = model.data.all[!duplicated(model.data.all$x),])
# paste("% Deviance Explained: ",round(((model$null.deviance - model$deviance)/model$null.deviance)*100,2),sep="")
# 
# summary(model)

#model data with glm and summarize model output data
glm_sums <- data.frame(character(),character(),numeric(),numeric(),numeric(),numeric(),numeric(),stringsAsFactors=FALSE,row.names=NULL)
colnames(glm_sums) <- c("Species","Variable","Coefficient","Std. Error","$Z\\text{-value}$","$\\PRZ$","ANOVA")
for (i in 1:nrow(species.table)) {
  x <- coll.glm[[i]]
  x2 <- anova(x)
  x.names <- c("Intercept", toupper(paste(species.table[i,2])), "TVOL", "TVOL$^2$", "TSPD")
  #x.names <- c("Intercept", toupper(paste(species.table[i,2])), "TVOL", "TVOL2", "TSPD")
  x.species <- c(paste(species.list[i],sep=""), NA, NA, NA, NA)
  x.coef <- signif(coef(summary(x))[,1],digits=4)
  x.se <- signif(coef(summary(x))[,2],digits=4)
  x.zvalue <- signif(coef(summary(x))[,3],digits=4)
  x.prz <- signif(coef(summary(x))[,4],digits=2)
  x.prz <- sapply(x.prz, function(x) ifelse(x < .001, "$<$.001*", x))
  #x.prz <- sapply(x.prz, function(x) ifelse(x < .001, "<.001*", x))
  x.anova <- signif((x2[,2]/sum(x2[2:5,2]))*100,digits=4)
  x.anova[1] <- "---"
  x.all <- data.frame(cbind(x.species,x.names,x.coef,x.se,x.zvalue,x.prz,x.anova),stringsAsFactors=FALSE,row.names=NULL) 
  colnames(x.all) <- c("Species","Variable","Coefficient","Std. Error","Z-value","Pr(Z)","ANOVA")

  newrow = rep(NA,length(x.all))
  glm_sums <- rbind(glm_sums, x.all, newrow)

  rm(x)
  rm(x2)
  rm(x.names,x.species,x.coef,x.se,x.zvalue,x.prz,x.anova)
  rm(x.all)
  rm(newrow)
}

#print(xtable(glm_sums), include.rownames=FALSE, sanitize.text.function=function(x){x}, floating=FALSE)

write.csv(glm_sums, file = "output/glm_sums.csv", row.names=FALSE)


#construct model performance table
# coll_deviance <- data.frame("SPP"=rep(NA,nrow(species.table)),"DEV"=rep(NA,nrow(species.table)),"ROC"=rep(NA,nrow(species.table)))
# for(i in 1:nrow(species.table)) {
#   coll_deviance[i,"SPP"] <- toupper(paste(species.table[i,2]))
#   coll_deviance[i,"DEV"] <- get(paste(species.table[i,2],"coll.devexp",sep=""))
#   coll_deviance[i,"ROC"] <- get(paste(species.table[i,2],"coll.roc",sep=""))
# }
# write.csv(coll_deviance, file = "output/coll_devs.csv", row.names=FALSE)


#Validate model fit with independent data
load("data/coll_val_data")

registerDoMC(detectCores() - 1)
coll.val.roc <- foreach(i = 1:nrow(species.table)) %dopar% {
  #formula <- as.formula(paste0("coll ~ log(",species.table[i,2],") + log(tvol) + I(log(tvol)^2) + log(tspd)"))
  #data <- as.data.frame(model.data[[i]])
  #model <- glm(formula = formula, offset=log(length*4), family=binomial(link = "cloglog"), data = data)
  val.pred.glm <- predict(coll.glm[[i]], val.data[[i]], type="response")  #Make predictions with regression model fit
  roc(val.data[[i]]$coll, val.pred.glm)  #Compare collision records to predictions using receiver operator characteristic (ROC) function and report value
}
write.csv(coll.val.roc, file = "output/coll_ind_roc.csv", row.names=FALSE)

registerDoMC(detectCores() - 1)
coll.val.tss <- foreach(i = 1:nrow(species.table)) %dopar% {
  df <- cbind(val.data[[i]]$uid,val.data[[i]]$coll,predict(coll.glm[[i]], val.data[[i]], type="response"))
  opt <- optimal.thresholds(df,opt.methods="MaxSens+Spec")[,2]
  cfs <- cmx(df, threshold = opt)
  sens <- sensitivity(cfs, st.dev = TRUE)
  spec <- specificity(cfs, st.dev = TRUE)
  round(sens[,1] + spec[,1] - 1,2)
}
write.csv(coll.val.tss, file = "output/coll_ind_tss.csv", row.names=FALSE)

#combine all validation data to assess all species model
# val.data2 <- val.data
# 
# for(i in 1:nrow(species.table)){
#   colnames(val.data2[[i]]) <- c("uid","x","y","occ","tvol","tspd","coll")
# }
# 
# val.data.all <- do.call(rbind,val.data2)
# 
# val.pred.glm <- predict(model, val.data.all[!duplicated(val.data.all$x),], type="response")  #Make predictions with regression model fit
# roc(val.data.all[!duplicated(val.data.all$x),coll], val.pred.glm)  #Compare collision records to predictions using receiver operator characteristic (ROC) function and report value


#make predictions based on models
load("data/cov_data")

registerDoMC(detectCores() - 1)
glm.preds <- foreach(i = 1:nrow(species.table)) %dopar% {
  #formula <- as.formula(paste0("coll ~ log(",species.table[i,2],") + log(tvol) + I(log(tvol)^2) + log(tspd)"))
  #data <- model.data[[i]]
  #model <- coll.glm[[i]]
  preds <- predict(coll.glm[[i]], cov.data, type="response")
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

dbWriteTable(con, c("gis_victoria", "vic_nogeom_roads_6spcollrisk"), value = coll_preds, row.names=FALSE, overwrite=TRUE)

#calculate portions of roads that are within top 1% of risk for each species - same for all species
prop_risk <- coll_preds

n <- .1

for(i in 1:nrow(species.table)) {
  x <- coll_preds[,i+1]
  x[x >= quantile(x,prob=1-n/100, na.rm = TRUE)] <- 1
  x[!(x >= quantile(x,prob=1-n/100, na.rm = TRUE))] <- 0
  prop_risk[,i+1] <- x
}

prop_risk$sums <- rowSums(prop_risk[,-1])

sum(na.omit(prop_risk$sums[prop_risk$sums>0]))
sum(na.omit(prop_risk$sums[prop_risk$sums>1]))
sum(na.omit(prop_risk$sums[prop_risk$sums>2]))
sum(na.omit(prop_risk$sums[prop_risk$sums>3]))
sum(na.omit(prop_risk$sums[prop_risk$sums>4]))
sum(na.omit(prop_risk$sums[prop_risk$sums>5]))

######Models by season######

load("data/coll_model_data_sum")
registerDoMC(6)
coll.glm.summer <- foreach(i = 1:nrow(species.table)) %dopar% {
  formula <- as.formula(paste0("coll ~ log(",species.table[i,2],") + log(tvol) + I(log(tvol)^2) + log(tspd)"))
  model <- brglm(formula = formula, offset=log(length*4), family=binomial(link = "cloglog"), data = model.data.summer[[i]])
}
save(coll.glm.summer, file="output/coll_glm_summer")
rm(coll.glm.summer)
rm(model.data.summer)

load("data/coll_model_data_aut")
registerDoMC(6)
coll.glm.autumn <- foreach(i = 1:nrow(species.table)) %dopar% {
  formula <- as.formula(paste0("coll ~ log(",species.table[i,2],") + log(tvol) + I(log(tvol)^2) + log(tspd)"))
  model <- brglm(formula = formula, offset=log(length*4), family=binomial(link = "cloglog"), data = model.data.autumn[[i]])
}
save(coll.glm.autumn, file="output/coll_glm_autumn")
rm(coll.glm.autumn)
rm(model.data.autumn)

load("data/coll_model_data_win")
registerDoMC(6)
coll.glm.winter <- foreach(i = 1:nrow(species.table)) %dopar% {
  formula <- as.formula(paste0("coll ~ log(",species.table[i,2],") + log(tvol) + I(log(tvol)^2) + log(tspd)"))
  model <- brglm(formula = formula, offset=log(length*4), family=binomial(link = "cloglog"), data = model.data.winter[[i]])
}
save(coll.glm.winter, file="output/coll_glm_winter")
rm(coll.glm.winter)
rm(model.data.winter)

load("data/coll_model_data_spr")
registerDoMC(6)
coll.glm.spring <- foreach(i = 1:nrow(species.table)) %dopar% {
  formula <- as.formula(paste0("coll ~ log(",species.table[i,2],") + log(tvol) + I(log(tvol)^2) + log(tspd)"))
  model <- brglm(formula = formula, offset=log(length*4), family=binomial(link = "cloglog"), data = model.data.spring[[i]])
}
save(coll.glm.spring, file="output/coll_glm_spring")
rm(coll.glm.spring)
rm(model.data.spring)
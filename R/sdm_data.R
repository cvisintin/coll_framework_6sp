require(maptools)
require(raster)
require(ncf)
require(xtable)
require(rgdal)
require(rgeos)
require(RPostgreSQL)
require(doMC)
require(sp)
require(dismo)
require(fields)
require(spatstat)

drv <- dbDriver("PostgreSQL")  #Specify a driver for postgreSQL type database
con <- dbConnect(drv, dbname="qaeco_spatial", user="qaeco", password="Qpostgres15", host="boab.qaeco.com", port="5432")  #Connection to database server on Boab

species.table <- read.delim("data/species_list.csv", header=T, sep=",")


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
vars <- stack(c(mget(grid.names),"X"=X,"Y"=Y)) #Combine all maps to single stack
save(vars, file = "data/vars")

#vars.cor <- layerStats(vars, 'pearson', na.rm=TRUE)
#write.csv(vars.cor[[1]], file = "data/vars_cor.csv")

################Datasets using 10,000 random background points##################

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
  assign(paste(species.table[i,2],".data",sep=""),read.delim(paste("data/",species.table[i,2],".data",sep=""), header=T, sep=","))
}

brt_data <- data.frame("SPP"=rep(NA,nrow(species.table)),"SO_N"=rep(NA,nrow(species.table)),"SO_P"=rep(NA,nrow(species.table)),"SO_A"=rep(NA,nrow(species.table)))
for(i in 1:nrow(species.table)) {
  brt_data[i,"SPP"] <- toupper(paste(species.table[i,2]))
  data <- get(paste(species.table[i,2],".data",sep=""))
  brt_data[i,"SO_N"] <- nrow(data)
  brt_data[i,"SO_P"] <- nrow(data[data$occ == 1, ])
  brt_data[i,"SO_B"] <- nrow(data[data$occ == 0, ])
  brt_data[i,"ALL"] <- paste0(brt_data[i,"N"]," : ",brt_data[i,"SO_P"],"P / ",brt_data[i,"SO_B"],"B")
  rm(data)
}
write.csv(brt_data, file = "data/brt_data_meta.csv", row.names=FALSE)

###########Datasets using target group background##########
# registerDoMC(detectCores() - 1)
# 
# data <- foreach(i = 1:nrow(species.table)) %dopar% {
#   #data1 <- read.delim(paste("/home/casey/Research/Projects/SDMs/Data/",species.table[i,3],".csv",sep=""), header=T, sep=",")
#   drv <- dbDriver("PostgreSQL")  #Specify a driver for postgreSQL type database
#   con <- dbConnect(drv, dbname="qaeco_spatial", user="qaeco", password="Qpostgres15", host="boab.qaeco.com", port="5432")  #Connection to database server on Boab
#   
#   data1 <- dbGetQuery(con,paste0("
#                                  SELECT
#                                  ST_X(pts.geom) AS X, ST_Y(pts.geom) AS Y, CAST(1 AS INTEGER) AS OCC
#                                  FROM
#                                  gis_victoria.vic_gda9455_fauna_vba AS pts, gis_victoria.vic_gda9455_admin_state AS poly
#                                  WHERE
#                                  ST_contains(poly.geom, pts.geom)
#                                  AND
#                                  pts.start_year >= '2000'
#                                  AND
#                                  sci_name = '",species.table[i,3],"'
#                                  GROUP BY
#                                  pts.geom;
#                                  "))
#   
#   data0 <- dbGetQuery(con,paste0("
#                                  SELECT
#                                  ST_X(pts.geom) AS X, ST_Y(pts.geom) AS Y, CAST(0 AS INTEGER) AS OCC
#                                  FROM
#                                  gis_victoria.vic_gda9455_fauna_vba AS pts, gis_victoria.vic_gda9455_admin_state AS poly
#                                  WHERE
#                                  ST_contains(poly.geom, pts.geom)
#                                  AND
#                                  pts.start_year >= '2000'
#                                  GROUP BY
#                                  pts.geom;
#                                  "))
#   data <- rbind(data1,data0)
#   samples.df <- extract(vars,data[,1:2])
#   data <- cbind(data,samples.df)
#   data <- na.omit(data)
#   write.csv(data, file = paste("data/",species.table[i,2],"_tb.data",sep=""), row.names=FALSE)
#   data
# }
# 
# ###########Datasets with added spatial autocovariate##########
# 
# data0 <- read.csv("data/bg_data_pts.csv")
# colnames(data0) <- c("x","y","occ")
# 
# vic.rst0 <- vic.rst
# vic.rst0[vic.rst0 == 1] <- 0
# 
# for(i in 1:nrow(species.table)) {
#   #data1 <- read.delim(paste("/home/casey/Research/Projects/SDMs/Data/",species.table[i,3],".csv",sep=""), header=T, sep=",")
#   data1 <- dbGetQuery(con,paste0("
#                                  SELECT
#                                  ST_X(pts.geom) AS X, ST_Y(pts.geom) AS Y, CAST(1 AS INTEGER) AS OCC
#                                  FROM
#                                  gis_victoria.vic_gda9455_fauna_vba AS pts, gis_victoria.vic_gda9455_admin_state AS poly
#                                  WHERE
#                                  ST_contains(poly.geom, pts.geom)
#                                  AND
#                                  pts.start_year >= '2000'
#                                  AND
#                                  sci_name = '",species.table[i,3],"'
#                                  GROUP BY
#                                  pts.geom;
#                                  "))
#   data1.1 <- thin.algorithm(data1, 1000, 30)
#   data <- rbind(data1.1,data0)
#   #data <- rbind(data1,data0)
#   samples.df <- extract(vars,data[,1:2])
#   data <- cbind(data,samples.df)
#   data <- na.omit(data)
#   
#   #data.coords <- data.frame("x"=data[data$occ==1,1], "y"=data[data$occ==1,2])
#   #coordinates(data.coords) <- ~x+y
#   #d <- gDistance(data.coords, byid=T)
#   #min.d <- apply(d, 1, function(x) sort(x[x>0], decreasing=F)[1])
#   # radius <- ceiling(mean(min.d))
#   #weights.pr <- unlist(lapply(min.d, function(x) x/sum(min.d)))+1
#   #weights.bg <- rep(1,nrow(data[data$occ==0,]))
#   #weights.brt <- c(weights.pr, weights.bg)
#   
#   # rast_brt <- raster(ncol=822, nrow = 563, ymn = 5661000, ymx = 6224000, xmn = -58000, xmx = 764000)
#   # res(rast_brt) <- 1000
#   # xy_res_brt <- cbind(data[,1], data[,2], data[,3])
#   # rast_brt[cellFromXY(rast_brt, xy_res_brt)] <- xy_res_brt[,3]
#   # cf <- focalWeight(rast_brt, c(10,0.5), "Gaus")
#   # focal_ac_rast_brt <- focal(rast_brt, w=cf, fun = mean, na.rm = TRUE)
#   # s <- stack(focal_ac_rast_brt, vic.rst0)
#   # focal_ac_rast_brt <- sum(s, na.rm=TRUE)
#   # AUTOCOV <- extract(focal_ac_rast_brt, cbind(data[,1], data[,2]))
#   # data <- cbind(data,AUTOCOV)
#   
#   #data <- cbind(data, weights.brt)
#   
#   # writeRaster(focal_ac_rast_brt, filename=paste0("output/",toupper(species.table[i,2]),"_ac_rast.tif"), format="GTiff", overwrite=TRUE)
#   
#   write.csv(data, file = paste("data/",species.table[i,2],"_sp.data",sep=""), row.names=FALSE)
#   
#   # rm(data1,data,samples.df,data.coords,d,min.d,radius,rast_brt,xy_res_brt,cf,focal_ac_rast_brt,s,AUTOCOV)
#   #rm(data1,data,samples.df,data.coords,d,min.d,weights.pr,weights.bg,weights.brt)
#   
#   rm(data1,data,samples.df)
# }
# rm(data0)
# 
# # AC <- autocov_dist(rast_brt[rast_brt>=0], rasterToPoints(rast_brt, spatial=TRUE), nbs=75000, zero.policy=TRUE)
# #####################Address spatial autocorrelation###############################
# 
# #find maximum and mean distances between nearest neighbouring occurrence points
# data.coords <- data.frame("x"=data[data$occ==1,1], "y"=data[data$occ==1,2])
# coordinates(data.coords) <- ~x+y
# d <- gDistance(data.coords, byid=T)
# min.d <- apply(d, 1, function(x) sort(x[x>0], decreasing=F)[1])
# rm(d)
# max(min.d)  ###125227.6 metres for EGK
# mean(min.d)  ###3557.865 metres for EGK
# inv.dist2 <- 1/min.d^2
# sq.dist <- sqrt(min.d)
# 
# test <- glm(formula=formula, family="binomial", data=data)
# AC <- autocov_dist(data[,3], cbind(data[,1], data[,2]), nbs=75000, zero.policy=TRUE)
# formula2 <- as.formula(paste0(colnames(data)[3]," ~ ",paste(colnames(data)[4:21],collapse=" + ")," + AC"))
# test2 <- glm(formula=formula2, family="binomial", data=cbind(data,AC))
# 
# s <- stack(rast_brt, vic.rst0)
# rast_brt <- sum(s, na.rm=TRUE)
# nb <- cell2nb(nrow = nrow(rast_brt), ncol = ncol(rast_brt), type="queen")
# lwb <- nb2listw(nb, style = "B")
# 
# joincount.test(as.factor(rast_brt@data@values), lwb, alternative = "greater")
# joincount.mc(as.factor(rast_brt@data@values), lwb, nsim = 10, alternative = "greater")
# 
# for(i in 1:20){
#   nb <- knn2nb(knearneigh(cbind(data[,1], data[,2]), k=i))
#   test <- joincount.test(as.factor(data[,3]), listw=nb2listw(nb))
#   print(c(test[[2]]$estimate,test[[2]]$p.value))
# }
# 
# binom.krige(coords=cbind(data[,1], data[,2]), data=data[,3], krige = list(cov.pars = c(1,1), beta = 1), mcmc.input = mcmc.control(S.scale = 0.2, thin = 1), locations=)
# 
# grid.data <- as.data.frame(rasterToPoints(rast_brt))
# preds <- predict(test,grid.data,type="response")
# 
# rast_brt[cellFromXY(rast_brt, grid.data)] <- preds
# 
# formula <- as.formula(paste0(colnames(data)[3]," ~ ",paste(colnames(data)[4:21],collapse=" + ")))
# bw <- ggwr.sel(formula, data = data, coords=cbind(data[,1], data[,2]), family = binomial)
# temp.gwr <- ggwr(formula, data = data, coords=cbind(data[,1], data[,2]), bandwidth=100000, gweight = gwr.Gauss, fit.points=grid.data, family = binomial)
# 
# spdf_brt <- rasterToPoints(rast_brt, spatial=TRUE)
# proj4string(spdf_brt) <- CRS("+init=epsg:28355") 
# spobj <- spodt(occ ~ 1, data=spdf_brt)
# 
# rast_wvt <- raster(ncol=822, nrow = 822, ymn = 5542000, ymx = 6364000, xmn = -58000, xmx = 764000)
# res(rast_wvt) <- 1000
# rast_wvt[rast_wvt] <- 0
# xyz <- cbind(data[,1], data[,2], data[,3])
# rast_wvt[cellFromXY(rast_wvt, xyz[xyz[,3]==1,1:2])] <- 1
# mat_brt <- as.matrix(rast_wvt)
# test <- dwt.2d(mat_brt, "haar", 9)
# test2 <- idwt.2d(test)
# rast_brt <- raster(test2)
# 
# 
# rast_wvt <- raster(ncol=822, nrow = 563, ymn = 5661000, ymx = 6224000, xmn = -58000, xmx = 764000)
# res(rast_wvt) <- 1000
# rast_wvt[rast_wvt] <- 0
# xyz <- cbind(data[,1], data[,2], data[,3])
# rast_wvt[cellFromXY(rast_wvt, xyz[xyz[,3]==1,1:2])] <- 1
# mat_brt <- as.matrix(rast_wvt)
# test <- modwt.2d(mat_brt, "haar", 10, boundary = "periodic")
# test2 <- imodwt.2d(test)
# rast_brt <- raster(test2)
# 
# summary(gam.bino <- gam(as.formula(paste0(colnames(data)[3]," ~ ",paste(colnames(data)[4:21],collapse=" + "), " + s(x, y)")), data=data, family=binomial))
# 
# vic.mat <- as.matrix(vic.rst)
# vic.mat <- vic.mat[nrow(vic.mat):1, ]
# vic.mat[!is.finite(vic.mat)] <- 0
# vic.mat2 <- matrix(FALSE,nrow=563,ncol=822)
# vic.mat2[] <- as.logical(vic.mat)
# 
# vic.win <- owin(xrange=c(-58,764), yrange=c(5661,6224), units=c("metre","metres"), mask=vic.mat2)
# 
# vic.ppp <- ppp(data[,1]/1000,data[,2]/1000, window=vic.win, marks=data[,3])
# 
# vic.idw <- idw(vic.ppp)
# rast_idw <- raster(as.SpatialGridDataFrame.im(vic.idw))
# rast_idw <- setExtent(rast_idw, extent(-58000, 764000, 5661000, 6224000))
# proj4string(rast_idw) <- CRS("+init=epsg:28355")
# AUTOCOV <- extract(rast_idw, cbind(data[,1], data[,2]))

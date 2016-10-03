require(maptools)
require(raster)
require(ggplot2)
require(xtable)
require(R2jags)
require(ncf)
require(doMC)
require(data.table)
require(RPostgreSQL)

drv <- dbDriver("PostgreSQL")  #Specify a driver for postgreSQL type database
con <- dbConnect(drv, dbname="qaeco_spatial", user="qaeco", password="Qpostgres15", host="boab.qaeco.com", port="5432")  #Connection to database server on Boab

grid.files <- list.files(path="output/",pattern=".tif")
grid.names <- unlist(strsplit(grid.files,"\\_preds_brt."))[(1:(2*(length(grid.files)))*2)-1][1:length(grid.files)]

vic.rst <- raster("data/VIC_GDA9455_GRID_STATE_1000.tif")

clip <- extent(-58000, 764000, 5661000, 6224000)

for (i in 1:length(grid.files)) {
  temp <- raster(paste0("output/",grid.files[i]))
  temp <- crop(temp, clip)
  assign(grid.names[i],temp * vic.rst)
}
sdm.vars <- stack(mget(grid.names))

rm(list = grid.names)
rm(vic.rst)

species.table <- read.delim("data/species_list.csv", header=T, sep=",")
species.list <- c("Eastern Grey Kangaroo","Common Brushtail Possum","Common Ringtail Possum","Black Swamp Wallaby","Common Wombat","Koala")

roads <- as.data.table(dbGetQuery(con,"
  SELECT
    r.uid as uid, ST_X(r.geom) AS x, ST_Y(r.geom) AS y
  FROM
	  (SELECT
      uid, ST_ClosestPoint(geom, ST_Centroid(geom)) AS geom
		FROM
      gis_victoria.vic_gda9455_roads_state) AS r
  "))
setkey(roads,uid)

tvol.preds <- as.data.table(read.csv("../coll_framework_egk/output/tvol_preds_rf.csv"))  #Read in collision data training set (presences/absences of collisions and covariates)
tspd.preds <- as.data.table(read.csv("../coll_framework_egk/output/tspd_preds_rf.csv"))  #Read in collision data training set (presences/absences of collisions and covariates)

cov.data <- Reduce(function(x, y) merge(x, y, all=TRUE), list(roads,tvol.preds,tspd.preds))
cov.data$uid <- as.integer(cov.data$uid)
setkey(cov.data,uid)

registerDoMC(detectCores() - 1)
sdm.preds <- foreach(i = 1:nrow(species.table), .packages = c("raster")) %dopar% {
  raster(paste0("output/",toupper(species.table[i,2]),"_preds_brt.tif"))
}

for(i in 1:nrow(species.table)){
  #cov.data[,paste0(species.table[i,2]):=raster::extract(sdm.preds[[i]],cov.data[,.(x,y)])]
  set(cov.data, ,paste0(species.table[i,2]), raster::extract(sdm.preds[[i]],cov.data[,.(x,y)]))
}

cov.data$coll <- as.integer(0)

save(cov.data, file="data/cov_data")

registerDoMC(detectCores() - 1)
coll <- foreach(i = 1:nrow(species.table), .packages = c("RPostgreSQL")) %dopar% {
  drv <- dbDriver("PostgreSQL")  #Specify a driver for postgreSQL type database
  con <- dbConnect(drv, dbname="qaeco_spatial", user="qaeco", password="Qpostgres15", host="boab.qaeco.com", port="5432")  #Connection to database server on Boab
  data <- as.data.table(dbGetQuery(con,paste0("
    SELECT DISTINCT ON (p.id)
      r.uid AS uid, CAST(1 AS INTEGER) AS coll
	  FROM
      gis_victoria.vic_gda9455_roads_state as r,
        (SELECT
          id, geom
        FROM
          gis_victoria.vic_gda9455_fauna_wv
        WHERE
          species = '",species.table[i,1],"'
        AND
          cause = 'hit by vehicle'
        AND
          year < 2013) AS p
    WHERE ST_DWithin(p.geom,r.geom,100)
    ORDER BY p.id, ST_Distance(p.geom,r.geom)
    "))
  )
  setkey(data,uid)
  unique(data)
}


#construct modelling datasets using all road segments
# registerDoMC(detectCores() - 1)
# model.data <- foreach(i = 1:nrow(species.table), .packages = c("RPostgreSQL")) %dopar% {
#   data1 <- coll[[i]]
#   data <- cov.data
#   data[data1, coll := i.coll]
#   na.omit(data[,c(1,i+5,4,5,12),with=FALSE])
# }

#construct modelling datasets using 2x ncoll randomly sampled raods
registerDoMC(detectCores() - 1)
model.data <- foreach(i = 1:nrow(species.table), .packages = c("RPostgreSQL")) %dopar% {
  data.coll <- coll[[i]]
  set.seed(123)
  data0 <- cov.data[sample(seq(1:nrow(cov.data)),2*nrow(data.coll))]
  data.all <- cov.data
  data.all[data.coll, coll := i.coll]
  data1 <- data.all[coll==1]
  data <- rbind(data0,data1)
  data <- na.omit(data[,c(1:3,i+5,4,5,12),with=FALSE])
  #cols <- colnames(data)[4:6]
  #data[, (cols) := lapply(.SD, function(x){log(x)}), .SDcols=cols]
  #data[, (cols) := lapply(.SD, scale), .SDcols=cols]
  data
}
save(model.data, file="data/coll_model_data")


#construct metadata for all species
coll_data <- data.frame("SPP"=rep(NA,nrow(species.table)),"CR_N"=rep(NA,nrow(species.table)),"CR_P"=rep(NA,nrow(species.table)),"CR_A"=rep(NA,nrow(species.table)))
for(i in 1:nrow(species.table)) {
  coll_data[i,"SPP"] <- toupper(paste(species.table[i,2]))
  data <- model.data[[i]]
  coll_data[i,"CR_N"] <- data[,.N]
  coll_data[i,"CR_P"] <- data[coll==1,.N]
  coll_data[i,"CR_A"] <- data[coll==0,.N]
  rm(data)
}
write.csv(coll_data, file = "data/coll_data_meta.csv", row.names=FALSE)


registerDoMC(detectCores() - 1)
coll.ind <- foreach(i = 1:nrow(species.table), .packages = c("RPostgreSQL")) %dopar% {
  drv <- dbDriver("PostgreSQL")  #Specify a driver for postgreSQL type database
  con <- dbConnect(drv, dbname="qaeco_spatial", user="qaeco", password="Qpostgres15", host="boab.qaeco.com", port="5432")  #Connection to database server on Boab
  data <- as.data.table(dbGetQuery(con,paste0("
    SELECT DISTINCT ON (p.id)
      r.uid AS uid, CAST(1 AS INTEGER) AS coll
	  FROM
      gis_victoria.vic_gda9455_roads_state as r,
        (SELECT
          id, geom
        FROM
          gis_victoria.vic_gda9455_fauna_wv
        WHERE
          species = '",species.table[i,1],"'
        AND
          cause = 'hit by vehicle'
        AND
          year >= 2013
        AND
          year != 2014) AS p
    WHERE ST_DWithin(p.geom,r.geom,100)
    ORDER BY p.id, ST_Distance(p.geom,r.geom)
  ")) #~1 second query
  )
  setkey(data,uid)
  unique(data)
}

#construct validation datasets using 2x ncoll randomly sampled raods
registerDoMC(detectCores() - 1)
val.data <- foreach(i = 1:nrow(species.table), .packages = c("RPostgreSQL")) %dopar% {
  data.coll <- coll.ind[[i]]
  set.seed(123)
  data0 <- cov.data[sample(seq(1:nrow(cov.data)),2*nrow(data.coll))]
  data.all <- cov.data
  data.all[data.coll, coll := i.coll]
  data1 <- data.all[coll==1]
  
  data <- rbind(data0,data1)
  na.omit(data[,c(1:3,i+5,4,5,12),with=FALSE])
}

save(val.data, file="data/coll_val_data")

#construct seasonal model data sets
registerDoMC(detectCores() - 1)
coll.summer <- foreach(i = 1:nrow(species.table), .packages = c("RPostgreSQL")) %dopar% {
  drv <- dbDriver("PostgreSQL")  #Specify a driver for postgreSQL type database
  con <- dbConnect(drv, dbname="qaeco_spatial", user="qaeco", password="Qpostgres15", host="boab.qaeco.com", port="5432")  #Connection to database server on Boab
  data <- as.data.table(dbGetQuery(con,paste0("
    SELECT DISTINCT ON (p.id)
      r.uid AS uid, CAST(1 AS INTEGER) AS coll
	  FROM
      gis_victoria.vic_gda9455_roads_state as r,
        (SELECT
          id, geom
        FROM
          gis_victoria.vic_gda9455_fauna_wv
        WHERE
          species = '",species.table[i,1],"'
        AND
          cause = 'hit by vehicle'
        AND
          year < 2013
        AND
          (month = 12
          OR
          month = 01
          OR
          month = 02
          )) AS p
    WHERE ST_DWithin(p.geom,r.geom,100)
    ORDER BY p.id, ST_Distance(p.geom,r.geom)
    "))
  )
  setkey(data,uid)
  unique(data)
}

registerDoMC(detectCores() - 1)
model.data.summer <- foreach(i = 1:nrow(species.table), .packages = c("RPostgreSQL")) %dopar% {
  data.coll <- coll.summer[[i]]
  set.seed(123)
  data0 <- cov.data[sample(seq(1:nrow(cov.data)),2*nrow(data.coll))]
  data.all <- cov.data
  data.all[data.coll, coll := i.coll]
  data1 <- data.all[coll==1]
  
  data <- rbind(data0,data1)
  na.omit(data[,c(1:3,i+5,4,5,12),with=FALSE])
}
save(model.data.summer, file="data/coll_model_data_sum")


registerDoMC(detectCores() - 1)
coll.autumn <- foreach(i = 1:nrow(species.table), .packages = c("RPostgreSQL")) %dopar% {
  drv <- dbDriver("PostgreSQL")  #Specify a driver for postgreSQL type database
  con <- dbConnect(drv, dbname="qaeco_spatial", user="qaeco", password="Qpostgres15", host="boab.qaeco.com", port="5432")  #Connection to database server on Boab
  data <- as.data.table(dbGetQuery(con,paste0("
    SELECT DISTINCT ON (p.id)
      r.uid AS uid, CAST(1 AS INTEGER) AS coll
	  FROM
      gis_victoria.vic_gda9455_roads_state as r,
        (SELECT
          id, geom
        FROM
          gis_victoria.vic_gda9455_fauna_wv
        WHERE
          species = '",species.table[i,1],"'
        AND
          cause = 'hit by vehicle'
        AND
          year < 2013
        AND
          (month = 03
          OR
          month = 04
          OR
          month = 05
          )) AS p
    WHERE ST_DWithin(p.geom,r.geom,100)
    ORDER BY p.id, ST_Distance(p.geom,r.geom)
    "))
  )
  setkey(data,uid)
  unique(data)
}

registerDoMC(detectCores() - 1)
model.data.autumn <- foreach(i = 1:nrow(species.table), .packages = c("RPostgreSQL")) %dopar% {
  data.coll <- coll.autumn[[i]]
  set.seed(123)
  data0 <- cov.data[sample(seq(1:nrow(cov.data)),2*nrow(data.coll))]
  data.all <- cov.data
  data.all[data.coll, coll := i.coll]
  data1 <- data.all[coll==1]
  
  data <- rbind(data0,data1)
  na.omit(data[,c(1:3,i+5,4,5,12),with=FALSE])
}
save(model.data.autumn, file="data/coll_model_data_aut")


registerDoMC(detectCores() - 1)
coll.winter <- foreach(i = 1:nrow(species.table), .packages = c("RPostgreSQL")) %dopar% {
  drv <- dbDriver("PostgreSQL")  #Specify a driver for postgreSQL type database
  con <- dbConnect(drv, dbname="qaeco_spatial", user="qaeco", password="Qpostgres15", host="boab.qaeco.com", port="5432")  #Connection to database server on Boab
  data <- as.data.table(dbGetQuery(con,paste0("
    SELECT DISTINCT ON (p.id)
      r.uid AS uid, CAST(1 AS INTEGER) AS coll
	  FROM
      gis_victoria.vic_gda9455_roads_state as r,
        (SELECT
          id, geom
        FROM
          gis_victoria.vic_gda9455_fauna_wv
        WHERE
          species = '",species.table[i,1],"'
        AND
          cause = 'hit by vehicle'
        AND
          year < 2013
        AND
          (month = 06
          OR
          month = 07
          OR
          month = 08
          )) AS p
    WHERE ST_DWithin(p.geom,r.geom,100)
    ORDER BY p.id, ST_Distance(p.geom,r.geom)
    "))
  )
  setkey(data,uid)
  unique(data)
}

registerDoMC(detectCores() - 1)
model.data.winter <- foreach(i = 1:nrow(species.table), .packages = c("RPostgreSQL")) %dopar% {
  data.coll <- coll.winter[[i]]
  set.seed(123)
  data0 <- cov.data[sample(seq(1:nrow(cov.data)),2*nrow(data.coll))]
  data.all <- cov.data
  data.all[data.coll, coll := i.coll]
  data1 <- data.all[coll==1]
  
  data <- rbind(data0,data1)
  na.omit(data[,c(1:3,i+5,4,5,12),with=FALSE])
}
save(model.data.winter, file="data/coll_model_data_win")


registerDoMC(detectCores() - 1)
coll.spring <- foreach(i = 1:nrow(species.table), .packages = c("RPostgreSQL")) %dopar% {
  drv <- dbDriver("PostgreSQL")  #Specify a driver for postgreSQL type database
  con <- dbConnect(drv, dbname="qaeco_spatial", user="qaeco", password="Qpostgres15", host="boab.qaeco.com", port="5432")  #Connection to database server on Boab
  data <- as.data.table(dbGetQuery(con,paste0("
    SELECT DISTINCT ON (p.id)
      r.uid AS uid, CAST(1 AS INTEGER) AS coll
	  FROM
      gis_victoria.vic_gda9455_roads_state as r,
        (SELECT
          id, geom
        FROM
          gis_victoria.vic_gda9455_fauna_wv
        WHERE
          species = '",species.table[i,1],"'
        AND
          cause = 'hit by vehicle'
        AND
          year < 2013
        AND
          (month = 09
          OR
          month = 10
          OR
          month = 11
          )) AS p
    WHERE ST_DWithin(p.geom,r.geom,100)
    ORDER BY p.id, ST_Distance(p.geom,r.geom)
    "))
  )
  setkey(data,uid)
  unique(data)
}

registerDoMC(detectCores() - 1)
model.data.spring <- foreach(i = 1:nrow(species.table), .packages = c("RPostgreSQL")) %dopar% {
  data.coll <- coll.spring[[i]]
  set.seed(123)
  data0 <- cov.data[sample(seq(1:nrow(cov.data)),2*nrow(data.coll))]
  data.all <- cov.data
  data.all[data.coll, coll := i.coll]
  data1 <- data.all[coll==1]
  
  data <- rbind(data0,data1)
  na.omit(data[,c(1:3,i+5,4,5,12),with=FALSE])
}
save(model.data.spring, file="data/coll_model_data_spr")
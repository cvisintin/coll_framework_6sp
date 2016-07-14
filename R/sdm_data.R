require("maptools")
require("raster")
require("ncf")
require("xtable")
require("rgdal")
require("rgeos")
require("RPostgreSQL")

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
vars <- stack(mget(grid.names)) #Combine all maps to single stack

#vars.cor <- layerStats(vars, 'pearson', na.rm=TRUE)
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
  assign(paste(species.table[i,2],".data",sep=""),read.delim(paste("data/",species.table[i,2],".data",sep=""), header=T, sep=","))
}

brt_data <- data.frame("SPP"=rep(NA,nrow(species.table)),"SO_N"=rep(NA,nrow(species.table)),"SO_P"=rep(NA,nrow(species.table)),"SO_A"=rep(NA,nrow(species.table)))
for(i in 1:nrow(species.table)) {
  brt_data[i,"SPP"] <- toupper(paste(species.table[i,2]))
  data <- get(paste(species.table[i,2],".data",sep=""))
  brt_data[i,"SO_N"] <- nrow(data)
  brt_data[i,"SO_P"] <- nrow(data[data$occ == 1, ])
  brt_data[i,"SO_A"] <- nrow(data[data$occ == 0, ])
  rm(data)
}
write.csv(brt_data, file = "data/brt_data_meta.csv", row.names=FALSE)
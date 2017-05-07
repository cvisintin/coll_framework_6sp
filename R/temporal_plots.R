require(ggplot2)
require(doMC)
require(data.table)
require(RPostgreSQL)

plotPal <- c("#94d1c7", "#cccc2b", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#969696", "#bc80bd")

drv <- dbDriver("PostgreSQL")  #Specify a driver for postgreSQL type database
con <- dbConnect(drv, dbname="qaeco_spatial", user="qaeco", password="Qpostgres15", host="boab.qaeco.com", port="5432")  #Connection to database server on Boab

species.table <- read.delim("data/species_list.csv", header=T, sep=",")
species.list <- c("Eastern Grey Kangaroo","Common Brushtail Possum","Common Ringtail Possum","Swamp Wallaby","Common Wombat","Koala")

registerDoMC(detectCores() - 1)
coll.temp1 <- foreach(i = 1:nrow(species.table), .packages = c("RPostgreSQL"), .combine=rbind) %dopar% {
  drv <- dbDriver("PostgreSQL")  #Specify a driver for postgreSQL type database
  con <- dbConnect(drv, dbname="qaeco_spatial", user="qaeco", password="Qpostgres15", host="boab.qaeco.com", port="5432")  #Connection to database server on Boab
  data <- as.data.table(dbGetQuery(con,paste0("
                                              SELECT DISTINCT ON (p.id)
                                              r.uid AS uid, p.year AS year, p.month AS month, CAST(substring(p.time from 1 for 2) AS INTEGER) AS hour, CAST('",species.table[i,2],"' AS VARCHAR) AS name
                                              FROM
                                              gis_victoria.vic_gda9455_roads_state as r,
                                              (SELECT
                                              id, year, month, time, geom
                                              FROM
                                              gis_victoria.vic_gda9455_fauna_wv
                                              WHERE
                                              species = '",species.table[i,1],"'
                                              AND
                                              cause = 'hit by vehicle') AS p
                                              WHERE ST_DWithin(p.geom,r.geom,100)
                                              ORDER BY p.id, ST_Distance(p.geom,r.geom)
                                              "))
  )
  setkey(data,uid)
  unique(data)
}

coll.temp2 <- as.data.table(dbGetQuery(con,paste0("
                                              SELECT DISTINCT ON (p.id)
                                              r.uid AS uid, p.year AS year, p.month AS month, p.hour AS hour, CAST('",species.table[i,2],"' AS VARCHAR) AS name
                                              FROM
                                              gis_victoria.vic_gda9455_roads_state as r,
                                              (SELECT
                                              id, year, month, hour, geom
                                              FROM
                                              gis_victoria.vic_gda9455_fauna_wv_2015_egkcoll) AS p
                                              WHERE ST_DWithin(p.geom,r.geom,100)
                                              ORDER BY p.id, ST_Distance(p.geom,r.geom)
                                              "))
  )

coll.temporal <- rbind(coll.temp1,coll.temp2)
setkey(coll.temporal,name)

for (i in 1:nrow(species.table)) {
  png(paste0('/home/casey/Research/Projects/PhD_Thesis/completion_talk/graphics/',species.table[i,2],'_coll_hour.png'), pointsize = 16, res=300, width = 1000, height = 900)
  print(
    ggplot(coll.temporal[name==paste0(species.table[i,2]),.N,by='hour'],aes(x=hour,y=N)) +
      geom_bar(stat="identity", position="stack", colour=NA, fill=plotPal[i]) +
      ylab("TOTAL COLLISIONS") +
      xlab("HOUR") +
      theme_bw() +
      theme(legend.key = element_blank()) +
      theme(legend.position="none") +
      theme(plot.margin=unit(c(.5,.5,.1,.1),"cm")) +
      theme(axis.title.x = element_text(margin=unit(c(.3,0,0,0),"cm"))) +
      theme(axis.title.y = element_text(margin=unit(c(0,.3,0,0),"cm"))) +
      theme(panel.grid.major = element_line(size=0.1),panel.grid.minor = element_line(size=0.1)) +
      theme(text = element_text(size = 8)) +
      scale_x_continuous(breaks=seq(0,23,by=2), expand = c(0, 0), lim=c(-1,24))
  )
  dev.off()
}

for (i in 1:nrow(species.table)) {
  png(paste0('/home/casey/Research/Projects/PhD_Thesis/completion_talk/graphics/',species.table[i,2],'_coll_month.png'), pointsize = 16, res=300, width = 1000, height = 900)
  print(
  ggplot(coll.temporal[name==paste0(species.table[i,2]),.N,by='month'],aes(x=month,y=N)) +
    geom_bar(stat="identity", position="stack", colour=NA, fill=plotPal[i]) +
    ylab("TOTAL COLLISIONS") +
    xlab("MONTH") +
    theme_bw() +
    theme(legend.key = element_blank()) +
    theme(legend.position="none") +
    theme(plot.margin=unit(c(.5,.5,.1,.1),"cm")) +
    theme(axis.title.x = element_text(margin=unit(c(.3,0,0,0),"cm"))) +
    theme(axis.title.y = element_text(margin=unit(c(0,.3,0,0),"cm"))) +
    theme(panel.grid.major = element_line(size=0.1),panel.grid.minor = element_line(size=0.1)) +
    theme(text = element_text(size = 8)) +
    scale_x_continuous(breaks=seq(1,12,by=1), expand = c(0, 0), lim=c(0,13))
  )
  dev.off()
}

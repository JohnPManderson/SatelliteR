# code to get sat data for survey data
library(ncdf4)
library(raster)

# the sat data
# f <- nc_open("http://128.175.28.250:8080/thredds/dodsC/Aqua3DayAggregate.nc")
                                        # bins<-"3 day"
 # bins<-"1 day"
f <- nc_open("http://basin.ceoe.udel.edu/thredds/dodsC/Aqua1DayAggregate.nc")

lon = ncvar_get(f, "lon")
lat = ncvar_get(f, "lat")

#get and format the satellite times
SatTimes = ncvar_get(f, "time")
Sat.Times = as.POSIXct(SatTimes, origin = "1970-01-01", tz = "GMT")
plot(Sat.Times, rep(1,length(Sat.Times)))

# When is the earliest sattime
Start_here<-Sat.Times[1]-86400

# import the survey data
file<-"~/Dropbox/Mackerel_COOP research/Documents/data/NEFSC_bottom_trawl_data/MackerelNEFSCbt_01112017_w_maturity.csv"

# read it in and subset for season
data<-read.delim(file, header=T, sep=",")
data<-subset(data, SEASON=="SPRING")

# format the station times and convert from est to gmt
StationTimes.est<-as.POSIXct(strptime(data$EST_TOWDATE, "%Y-%m-%d %H:%M:%S",tz="EST") )
attributes(StationTimes.est)$tzone<-"GMT"
StationTimes.utc<-StationTimes.est
rm(StationTimes.est)

# bind gmt times so you can select statiosn that overlap with sat data
data$StationTimes.utc<-StationTimes.utc
data<-subset(data, StationTimes.utc>Start_here)
rm( StationTimes.utc)
dim(data)

StationTimes.est<-as.POSIXct(strptime(data$EST_TOWDATE, "%Y-%m-%d %H:%M:%S",tz="EST") )
attributes(StationTimes.est)$tzone<-"GMT"
StationTimes.utc<-StationTimes.est
rm(StationTimes.est)

#buffersize in meters around station location for the extraction
buffer.size=1000

results<-data.frame(matrix(nrow=dim(data)[1], ncol=8))

for(i in 1:length(StationTimes.utc)){

#results[i, 1]<-data$key[i]
#results[i, 2]<-data$LAT[i]
#results[i, 3]<-data$LON[i]
#results[i, 4]<-data$StationTimes.utc[i]

#calculate the absolute value of time time differences between the station time and all the sat times and pick the satellite scene with the least time difference
timediffs<-abs(difftime(StationTimes.utc[i], Sat.Times, units="hours"))
this.sat.scene<-which.min(timediffs)
#results[i, 5]<-Sat.Times[this.sat.scene]
results[i,1]<-i
results[i, 2]<-timediffs[this.sat.scene]

# pick the station location out and narrow the sattelite query to a
lon_lat.stat<-data[i,9:8]
projstatcoords<-SpatialPoints(lon_lat.stat,  proj4string = CRS("+proj=merc +lat_ts=37.6960626707359 +lon_0=-75.0 +datum=WGS84"))

lon_ind = which(lon >= lon_lat.stat[,1]-.02 & lon <= lon_lat.stat[,1]+.02)
lat_ind = which(lat >= lon_lat.stat[,2]-.02 & lat <= lon_lat.stat[,2]+.02)

var<-"sst"

var.a= ncvar_get(f, var, start = c(min(lon_ind),min(lat_ind),this.sat.scene), count = c(length(lon_ind), length(lat_ind),1))
var.a <- t(var.a[, ncol(var.a):1])

my.rast<-raster(var.a,  xmn=min(lon[lon_ind])-nudgelon, xmx=max(lon[lon_ind])+nudgelon, ymn=min(lat[lat_ind])-nudgelat, ymx=max(lat[lat_ind])+nudgelat, crs="+proj=merc +lat_ts=37.6960626707359 +lon_0=-75.0 +datum=WGS84")

results[i, 3]<-extract(my.rast, projstatcoords, buffer=buffer.size, fun=median)

var<-"chl_oc3" #Chlorophyll Concentration, OC3 Algorithm

var.a= ncvar_get(f, var, start = c(min(lon_ind),min(lat_ind),this.sat.scene), count = c(length(lon_ind), length(lat_ind),1))
var.a <- t(var.a[, ncol(var.a):1])

my.rast<-raster(var.a,  xmn=min(lon[lon_ind])-nudgelon, xmx=max(lon[lon_ind])+nudgelon, ymn=min(lat[lat_ind])-nudgelat, ymx=max(lat[lat_ind])+nudgelat, crs="+proj=merc +lat_ts=37.6960626707359 +lon_0=-75.0 +datum=WGS84")

results[i, 4]<-extract(my.rast, projstatcoords, buffer=buffer.size, fun=median)

var<-"poc"  #Particulate Organic Carbon, D. Stramski, 2007 (443/555 version)

var.a= ncvar_get(f, var, start = c(min(lon_ind),min(lat_ind),this.sat.scene), count = c(length(lon_ind), length(lat_ind),1))
var.a <- t(var.a[, ncol(var.a):1])

my.rast<-raster(var.a,  xmn=min(lon[lon_ind])-nudgelon, xmx=max(lon[lon_ind])+nudgelon, ymn=min(lat[lat_ind])-nudgelat, ymx=max(lat[lat_ind])+nudgelat, crs="+proj=merc +lat_ts=37.6960626707359 +lon_0=-75.0 +datum=WGS84")

results[i, 5]<-extract(my.rast, projstatcoords, buffer=buffer.size, fun=median)

var<-"a_443_qaa" # Total absorption at 443 nm, QAA algorithm  alias for CDOM

var.a= ncvar_get(f, var, start = c(min(lon_ind),min(lat_ind),this.sat.scene), count = c(length(lon_ind), length(lat_ind),1))
var.a <- t(var.a[, ncol(var.a):1])

my.rast<-raster(var.a,  xmn=min(lon[lon_ind])-nudgelon, xmx=max(lon[lon_ind])+nudgelon, ymn=min(lat[lat_ind])-nudgelat, ymx=max(lat[lat_ind])+nudgelat, crs="+proj=merc +lat_ts=37.6960626707359 +lon_0=-75.0 +datum=WGS84")

results[i, 6]<-extract(my.rast, projstatcoords, buffer=buffer.size, fun=median)

var<-"M_WK" # Water Mass Classifications using Wards and Kmeans clustering after Oliver and Irwin 2008

var.a= ncvar_get(f, var, start = c(min(lon_ind),min(lat_ind),this.sat.scene), count = c(length(lon_ind), length(lat_ind),1))
var.a <- t(var.a[, ncol(var.a):1])

my.rast<-raster(var.a,  xmn=min(lon[lon_ind])-nudgelon, xmx=max(lon[lon_ind])+nudgelon, ymn=min(lat[lat_ind])-nudgelat, ymx=max(lat[lat_ind])+nudgelat, crs="+proj=merc +lat_ts=37.6960626707359 +lon_0=-75.0 +datum=WGS84")

results[i, 7]<-extract(my.rast, projstatcoords)

var<-"M_WK_G" # Gradient Strengths Across Water Mass Classifications using Wards and Kmeans clustering after Oliver and Irwin 2008

var.a= ncvar_get(f, var, start = c(min(lon_ind),min(lat_ind),this.sat.scene), count = c(length(lon_ind), length(lat_ind),1))
var.a <- t(var.a[, ncol(var.a):1])

my.rast<-raster(var.a,  xmn=min(lon[lon_ind])-nudgelon, xmx=max(lon[lon_ind])+nudgelon, ymn=min(lat[lat_ind])-nudgelat, ymx=max(lat[lat_ind])+nudgelat, crs="+proj=merc +lat_ts=37.6960626707359 +lon_0=-75.0 +datum=WGS84")

results[i, 8]<-extract(my.rast, projstatcoords, buffer=buffer.size, fun=median)
}

names(results)<-c("row_index","timediff.hrs","sst","chl_oc3","poc","a_443_qaa","M_WK","M_WK_G")

data_w_sats<- cbind(data,results)

head(data_w_sats, 1)

file<-"~/Dropbox/Mackerel_COOP research/Documents/data/NEFSC_bottom_trawl_data/MackerelNEFSCbt_withSatObs.csv"

write.table(data_w_sats, file, row.names=F, col.names=T, sep=",")


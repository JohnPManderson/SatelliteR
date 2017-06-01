rm(list=ls(all=TRUE))
library(ncdf4)
library(raster)
library(maptools)
library(colorRamps)

# bins<-"1 day"
f <- nc_open("http://basin.ceoe.udel.edu/thredds/dodsC/Aqua3DayAggregate.nc")

lon = ncvar_get(f, "lon")
lat = ncvar_get(f, "lat")

ndays<-40

#get and format the satellite times
SatTimes = ncvar_get(f, "time")
Sat.Times = as.POSIXct(SatTimes, origin = "1970-01-01", tz = "GMT")

mytimes<-seq(length(Sat.Times)-(ndays), length(Sat.Times), 1)

lon = ncvar_get(f, "lon")
lat = ncvar_get(f, "lat")



#NEAUS
lat_ind = which(lat >= 36 & lat <= 44)
lon_ind = which(lon >= -76 & lon <= -68)

nudgelon<-0
nudgelat<-0

#SST
#var<-"chl_oc3"                                       #other Vars
#var<-"sst"
#var<-"TSS_gould" #:Total Suspended Particles (PIM+POM)
#var<-"PIM_gould" #:Particulate Inorganic Matter From Sediment and Detrital Absorption at 412nm, Gould (200711)
#var<-"POM_gould" #:Particulate Organic Matter From Phytoplankton Absorption at 443nm, Gould (200711)
#var<-"salinity" #:Salinity - based on CDOM, Ladner algorithm
#var<-"M_WK" #:Water Mass Classifications using Wards and Kmeans clustering after Oliver and Irwin 2008
var<-"M_WK_G" #:long_name: Gradient Strengths Across Water Mass Classifications using Wards and Kmeans clustering after Oliver and Irwin 2008


for(i in 1:length(mytimes))
{
var.a= ncvar_get(f, var, start = c(min(lon_ind),min(lat_ind),mytimes[i]), count = c(length(lon_ind), length(lat_ind),1))
var.a <- t(var.a[, ncol(var.a):1])

rast<-raster(var.a,  xmn=min(lon[lon_ind])-nudgelon, xmx=max(lon[lon_ind])+nudgelon, ymn=min(lat[lat_ind])-nudgelat, ymx=max(lat[lat_ind])+nudgelat, crs="+proj=merc +lat_ts=37.6960626707359 +lon_0=-75.0 +datum=WGS84")
#yrast<-log(rast+1)
myrast<-rast

png(file=paste("~/Dropbox/Mine/Graphics/GRAD/GRAD/",var,"_",i,".png", sep=""), width=650, height=650, pointsize=18, bg="white")

plot(myrast, las=1,col=matlab.like(50), zlim=c(0,15), xlab="Longitude",  ylab="Latitude")#
lines(read.table("~/Dropbox/Mine/data/spatial data/bigcoast.dat", header=F, sep=""),lwd=2)
contour(myrast, nlevels=20, add=T, col="dark grey", lwd=.5)#col=matlab.like(20)
plot(readShapeLines("~/Dropbox/Mine/data/spatial data/bottomdata/na_15_bathy_contour/na_15_contour_50m.shp"), add=T, col="black", lty=3, lwd=1)
plot(readShapeLines("~/Dropbox/Pessutti&Manderson/spatial data/bottomdata/na_15_bathy_contour/Meters/na15_100m.shp"), add=T, col="black", lty=4, lwd=1)
plot(readShapeLines("~/Dropbox/Mine/data/spatial data/bottomdata/na_15_bathy_contour/na_15_contour_200m.shp"), add=T, col="black", lty=2, lwd=1)
grid()
title(main=paste(var,": ",Sat.Times[mytimes[i]], sep=""))
legend("topleft", c("50M", "100M","200M"), lwd=2, lty=c(3,4, 2), bty="n")

dev.off()

}


file.names<-data.frame(matrix(NA,nrow=length(mytimes), ncol=1))

for (i in 1:length(mytimes))
{
    file.names[i,1]<-paste(var,"_",i,".png",sep="")
    
    }

write.table(file.names,"~/Dropbox/Mine/Graphics/GRAD/GRAD/filenames.txt", col.names=F, row.names=F, sep="")


#, zlim=c(cellStats(my.rast, "min"),cellStats(my.rast, "max")))





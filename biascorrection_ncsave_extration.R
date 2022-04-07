  library(ncdf4)
  library(date)
  library(raster)
  library(bigmemory)
  library(qmap)
  library(dplyr)
  library(arsenal)
  install.packages("ggplot")
  install.packages("rqdatatable")
  library(rqdatatable)
  library(easyNCDF)
rm(list=ls())

###################################################################################################################
###################################################################################################################

############Fit qmap for reference period############################################

pcp45n<-nc_open("E:/Proyecciones_mensuales ROM grilla regular/Escalados/PCP45_06-19.nc")
chrp<-nc_open("E:/Proyecciones_mensuales ROM grilla regular/Escalados/CHIRPS_06-19.nc")
pcp<-nc_open("E:/Proyecciones_mensuales ROM grilla regular/Escalados/PCP45.nc")

chb<-ncvar_get(chrp,varid= "prec" )
pcp45<-ncvar_get(pcp45,varid= "var143" )
pcpt<-ncvar_get(pcp,varid= "var143" )

mchrp<-array(chb,dim = c(68,59,168))
mpcp<-array(pcpt,dim = c(68,59,168))
mpcp45<-array(pcpt,dim = c(68,59,1128))

######Getting rid of NA´s#################################################



ind_non_na <- which(!is.na(mpcp45), TRUE)

arrayInd(ind_non_na)

pcp_non_na<-subset(mpcp,!is.na(pcp45[,,]))
  chb_non_na<-subset(mchrp,!is.na(chb[,,]))
    pcpt_non_na<-subset(mpcp45,!is.na(pcpt[,,]))
     
########### FIT map####################################################################

    mapQUANTparP <- fitQmapQUANT(chb_non_na, pcp_non_na, wet.day=TRUE,qstep=0.01) #include wet day correction for P
save(mapQUANTparP, file = "E:/Proyecciones_mensuales ROM grilla regular/Escalados/quantmap45.RData") #Save parameters


###########do qmap##################
QQ45 <- doQmapQUANT(pcpt_non_na, mapQUANTparP)
save(QQ45, file = "E:/Proyecciones_mensuales ROM grilla regular/Escalados/qq45.RData")
####Make a df ###########################################################

matrixQQ45<-array(QQ45,dim = c(68,59,1128))

####MAking matrixes for the data###########
mchrp<-array(chb_non_na,dim = c(68,59,168))
mpcp<-array(pcpt_non_na,dim = c(68,59,1128))
mpcp45<-array(pcp_non_na,dim = c(68,59,1128))
#####plot##################################################################
plot(sqrtquant(matrixQQ45[35,30,]),type="l",col='red',xlab="quantile", ylab="Precipitation(?)")
  lines(sqrtquant(mpcp[35,30,]),type = "l",col="blue")
  lines(sqrtquant(mchrp[35,30,]),type = "l",col="black")
  legend("top",
         legend=c("QQ45","PCP45","Hist"),
         lty=c(1,1,1),pch=c(NA,NA,NA),
         col=c("red","blue","black"))
  
  library(hydroGOF)
i<-1
y<-1
kgedist<-array(dim=c(68,59))
for(i in 1:68){
  for(y in 1:59){
kgedist<-KGE(sim=qqmatrix[20,20,1:168], obs = mchrp[20,20,1:168], method = "2012", out.type = "full")
 }  
    }

plot(kgedist) 
#}
##########!NA matrix################
  #Quantile chrips vrs Quantile remapping 
  
  plot(sqrtquant(matrixQQ45[40,40,1:168]),sqrtquant(mchrp[40,40,1:168]),type="p",xlab="Quantiles QQ45", ylab="Quantiles Chrips")
  abline(a=4,b=1,col=2,lwd=2)  
  plot(sqrtquant(mpcp[40,40,1:168]),sqrtquant(mchrp[40,40,1:168]),type="p",xlab="Quantiles PCP45 raw", ylab="Quantiles Chrips")
  abline(a=0,b=1,col=2,lwd=3)
#### NA grid   
  5.054474e-05
    
#######QQplots#############################################################
  sqrtquant <- function(x,qstep=0.01){
    qq <- quantile(x,prob=seq(0,1,by=qstep,na.rm=TRUE))
    sqrt(qq)
  }
#####Gridded data########
  mpcp<-array(pcpt,dim = c(4525536))
  qq_non_na<-array(dim=c(4525536))
  y<-1
  for(i in 1:4525536)
   {
          if( !is.na(mpcp[i])){
           qq_non_na[i]=QQ45[y]
            y=y+1
            
          }
   
  }
  correctedlist
#####Near 0 values to 0 #####          
qqmatrix<-array(qq_non_na,dim = c(68,59,1128)) 
for(i in 1:4525536)
{
  if( !is.na(qq_non_na[i]))
    {
      if(qq_non_na[i] < 0.001)
        {
         qq_non_na[i]=0
    
    }
  }
}
  
qqmatrix<-array(qq_non_na,dim = c(68,59,1128))

library(ncdf4)
#######extract points for comparison####
qqmatrix_sj<-qqmatrix[39:42,28:30,]  
qq_non_na_sj<-as.data.frame(qq_non_na)  
write.csv(qqmatrix[39:42,28:30,1:168],"E:/Proyecciones_mensuales ROM grilla regular/Escalados/cuenca481_2.csv", row.names = FALSE)

chrp_sj<-mchrp[37,34,] 
qq_non_na_sj<-as.data.frame(qq_non_na)  
write.csv(mchrp,"E:/Proyecciones_mensuales ROM grilla regular/Escalados/datoschrp.csv", row.names = FALSE)
plot(qqmatrix[37,34,1:168],mchrp[38,34,1:168])
plot(qqmatrix[,,])
###Save as a netcdf
##Define dims
dat1=list()
  dat1$x=c((pcp45n[["dim"]][["lon"]][["vals"]]))
  dat1$y=c(pcp45n[["dim"]][["lat"]][["vals"]])
dat1$t=c(seq(as.Date("2006-1-1"), as.Date("2099-12-1"), by = "months"))
###Add dims to list 

##units for time 
ntd <- 12
tunits3 <- "months since 1900-01-01"

londim <- ncdim_def("Lon","Longitude",as.double(dat1$x)) 
latdim <- ncdim_def("Lat","Latitude",as.double(dat1$y)) 
timedim <- ncdim_def("time",tunits3,as.double.POSIXlt(dat1$t))

##Define variables
fillvalue <- 1e31
dlname <- "PCP" #Precipitation
pcp_def <- ncvar_def("pcp","mm",list(londim,latdim,timedim),fillvalue,dlname,prec="single")
##NC file dir 
ncpath <- "E:/Proyecciones_mensuales ROM grilla regular/Escalados/"
ncname <- "PCP45_ds7"  
ncfname <- paste(ncpath, ncname, ".nc", sep="")
##Create the nc file ############
ncout <- nc_create(ncfname,pcp_def,force_v4=TRUE)
ncvar_put(ncout,pcp_def,qqmatrix)
ncatt_put(ncout,"Lon","axis","X") 
ncatt_put(ncout,"Lat","axis","Y")
ncatt_put(ncout,"time","axis","T")

nc_close(ncout)
nc45<-nc_open(ncfname)
nctest<-ncvar_get(nc45,varid= "pcp" )
###Load the nc file as a raster stack#########
nc_path<-"E:/Proyecciones_mensuales ROM grilla regular/Escalados/"
ncname<-"basin_3"
ncfname<-paste(nc_path,ncname,".nc",sep = "")
nctest_raster<- stack(ncfname)
plot(nctest_raster)
##asignar las fechas en el stack
##cut nc for test basins###########
library(sf)
shp_path<-"E:/Proyecciones_mensuales ROM grilla regular/Escalados/"
shp_name<-"basin_3"
shp_file<-paste(shp_path,shp_name,".shp",sep="")
basin_3_shp<-st_read(shp_file)
plot(basin_3_shp,axis=TRUE)
basin_3<-crop(nctest_raster,basin_3_shp)
plot(basin_3)
basin_3_brick<-brick(basin_3)
###average the basin precipitation index###########
average<-1:1128
for (i in 1:1128){
  average[i]<-cellStats(basin_3[[i]],"mean")
}  

  
dates<-c(seq(as.Date("2006-1-1"), as.Date("2099-12-1"), by = "months"))
##save as data frame (bc it can holds more than one type of data)##############
basin3_array<-data.frame(dates,average)
hist(basin3_array$average)



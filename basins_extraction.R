library(terra)
library(hydroTSM)


dir   <- "E:/Proyecciones_mensuales ROM grilla regular/Escalados/PCP45_ds7.nc"
shape <- "E:/Proyecciones_mensuales ROM grilla regular/Escalados/Catchments_CostaRica_1.geojson"
output <- "E:/Proyecciones_mensuales ROM grilla regular/Escalados/temporales"
start <- "2006-01-01"
end   <- "2099-12-01"

shape <- vect(shape)
prod  <- rast(dir)
nms   <- shape$subid
dates <- mip(start, end)


for(i in 1:length(nms)){
  
  
  s <- shape[i,]
  
  ts <- terra::extract(prod, s, fun = mean, na.rm = TRUE, weights = TRUE)
  ts <- zoo(as.numeric(ts), dates)
  
  n <- paste0(nms[i], "_Subcatchment.zoo")
  write.zoo(ts, file.path(output, n))
  
}

plot(ts[1:168])

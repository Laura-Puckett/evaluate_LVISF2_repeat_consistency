# purpose: Spatially join LVIS footprints from two different files where footprint centers
# are within a set distance of each other, and add store information for all pairs in a single spatial object. 


library(sp); library(sf); library(nngeo); library(dplyr); library(corrgram); library(raster)
# ABoVE_crs =  sp::CRS("+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")

# ## Interpret inputs from bash script
# args = commandArgs(trailingOnly = T)
# distance_threshold = as.numeric(args[1]) # in meters
# 
# input_dir1 = args[2]
# input_dir2 = args[3]
# filename1 = args[4]
# filename2 = args[5]
# arrayNum = as.numeric(args[6])
# outDir = "./output/"

input_dir1 = "C:/Users/lkp58/Downloads/"
input_dir2 = "C:/Users/lkp58/Downloads/"
filename1 = "LVISF2_ABoVE2019_0722_R2003_075347.shp"
filename2 = "LVISF2_ABoVE2019_0722_R2003_079096.shp"
distance_threshold = 1
outDir="//nau.froot.nau.edu/cirrus/scratch/lkp58/LVIS/accuracy_assessment/output/"

arrayLength = 20
out_filename = paste0(substr(filename1, 0, nchar(filename1)-4),
                        '_',
                        substr(filename2, 0, nchar(filename1)-4))

print(paste0("distance threshold: ", distance_threshold, "m"))
print(paste0("file1: ", input_dir1, filename1))
print(paste0("file2: ", input_dir2, filename2))

## read in data
pts1_orig = st_read(dsn=paste0(input_dir1, filename1))
pts2_orig = st_read(dsn=paste0(input_dir2, filename2))

## assign crs of the sf objects (for some reason this wasn't carrying over from the shapefiles)
st_crs(pts1_orig) <- sp::CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
st_crs(pts2_orig) <- sp::CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

# ## transform to projected coordinate system so distance calculations will be quicker and more transparent
pts1_orig = st_transform(pts1_orig, sp::CRS("+proj=utm +zone=6 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
pts2_orig = st_transform(pts2_orig, sp::CRS("+proj=utm +zone=6 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))


## figure out indices for subsetting for parallelization
# (splits into 500 even breaks)
# indices = floor(seq(from=0, to=nrow(pts1), length.out = 501))
# start_index=indices[arrayNum]+1
# end_index = indices[arrayNum+1]
# print(paste0("indices ", start_index, " to ", end_index))
## subset pts 1
#pts1 = pts1[start_index:end_index,]

## split up by even breaks in y dimension

overlap_extent = st_intersection(
  st_set_crs(st_as_sf(as(raster::extent(pts1_orig), "SpatialPolygons")), st_crs(pts1_orig)), # extent of pts1
  st_set_crs(st_as_sf(as(raster::extent(pts2_orig), "SpatialPolygons")), st_crs(pts1_orig))) # extent of pts2

old_extent = raster::extent(overlap_extent)
lat1 = old_extent[3]
lat2 = old_extent[4]
lat_seq = seq(from=lat1, to=lat2, length.out = arrayLength)

## to be used later for calculating differences for each metric within the loop
diff_fun <- function(.){  return(paste0(., "_2") - paste0(., "_1"))} 

metricNames = c( "ZG","RH10","RH15", "RH20", "RH25", "RH30", "RH35", "RH40", "RH45",
                 "RH50", "RH55", "RH60", "RH65", "RH70", "RH75", "RH80","RH85",
                 "RH90", "RH95", "RH96", "RH97", "RH98", "RH99","RH100")
## 
combined = list()

for(arrayNum in 1:(arrayLength-1)){
  print(arrayNum)
  new_lat_min = lat_seq[arrayNum]
  new_lat_max = lat_seq[arrayNum+1]
  
  new_extent = st_set_crs(st_as_sf(as(raster::extent(old_extent[1], old_extent[2], new_lat_min, new_lat_max), "SpatialPolygons")), st_crs(pts1_orig))
  
  pts1 = st_intersection(pts1_orig, new_extent)
  pts2 = st_intersection(pts2_orig, new_extent)
  
  print(pts1)
  print(pts2)
  
  
  ## get points from second dataset that are within 1m of points from first dataset
  if(nrow(pts1) >0 && nrow(pts2) > 0){
    join_result = nngeo::st_nn(pts1, pts2, returnDist = T, k = 1, progress = F, maxdist = distance_threshold, sparse = T)
    unlist(join_result$dist)
    
    ## only execute the rest of the code if matches are found
    if(is.null(unlist(join_result$dist))){
      print("No matches found")
    }else{
      
      ## make table of indices from join
      joined_indices = tibble(dist = join_result$dist,
                              ID1 = 1:nrow(pts1),
                              ID2 = join_result$nn) %>% filter(dist<distance_threshold)
      
      pts1 = st_drop_geometry(pts1)
      pts2 = st_drop_geometry(pts2)
      
      ## rename columns from each dataset
      colnames_1 = colnames(pts1)
      colnames_2 = colnames(pts2)
      colnames(pts1) = paste0(colnames_1, "_1")
      colnames(pts2) = paste0(colnames_2, "_2")
      
      # ## specify the new name of the geometry attribute of the objects
      # sf::st_geometry(pts1) <- "geometry_1"
      # sf::st_geometry(pts2) <- "geometry_2"
      
      ## subset original pts objects by indices from join
      pts1_from_join = pts1[joined_indices$ID1,]
      pts2_from_join = pts2[unlist(joined_indices$ID2),]
      
      ## combine the columns into one sf object
      # **note spatial attributes of the object will be for pts1 only.
      # **since joined points are close to each other, they're basically the same location anyway - no need to 
      # ** have separate spatial objects for each point in the pair
      pts_from_join = cbind(pts1_from_join, pts2_from_join)
      pts_from_join = cbind(pts_from_join, "distance" = unlist(join_result$dist))
      
      
      ## calculate differences for each metric
      
      for(metric in metricNames){
        pts_from_join[paste0(metric, "_diff")] = pts_from_join[paste0(metric, "_2")][[1]] - pts_from_join[paste0(metric, "_1")][[1]]
      }
      
      ## summarize metadata by limiting factor in each pair
      # (min, max, or whether or not values match)
      if("SENSITI_1" %in% colnames(pts_from_join) &&
         "SENSITI_2" %in% colnames(pts_from_join)){
        pts_from_join["SENS_min"] = min(pts_from_join$SENSITI_1, pts_from_join$SENSITI_2)
      }
      pts_from_join["AZIM_max"] = max(pts_from_join$AZIMUTH_1, pts_from_join$AZIMUTH_2)
      pts_from_join["INCI_max"] = max(pts_from_join$INCIDEN_1, pts_from_join$INCIDEN_2)
      pts_from_join["COMP_max"] = max(pts_from_join$COMPLEX_1, pts_from_join$COMPLEX_2)
      pts_from_join["CHAN_ZG"] = ifelse(pts_from_join$CHANNEL_ZG_1 == pts_from_join$CHANNEL_ZG_2, T, F)
      pts_from_join["CHAN_ZT"] = ifelse(pts_from_join$CHANNEL_ZT_1 == pts_from_join$CHANNEL_ZT_2, T, F)
      pts_from_join["CHAN_R"] = ifelse(pts_from_join$CHANNEL_R_1 == pts_from_join$CHANNEL_R_2, T, F)
      
      combined = rbind(combined,  pts_from_join)
      
      ## transform pts_from_join back to WGS84
      # pts_from_join = st_stransform(pts_from_join, sp::CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
      
      # just having it rewrite with every iteration because I'm tired of it failing before the end without producing any output
      write.csv(x = combined,
                file = paste0(outDir, out_filename, '.csv'))
      
      
      points = sp::SpatialPointsDataFrame(coords = cbind(combined$GLON_1, combined$GLAT_1),
                                          data = as.data.frame(combined),
                                          proj4string = sp::CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
      
      rgdal::writeOGR(obj = points,
                      dsn = paste0(outDir, out_filename, '.shp'),
                      layer = out_filename,
                      driver = 'ESRI Shapefile',
                      overwrite_layer = T,
                      verbose = T)
    }
  }
}




# 
# write.csv(x = combined,
#          file = paste0(outDir, out_filename, '.csv'))
# 
# 
# points = sp::SpatialPointsDataFrame(coords = cbind(combined$GLON_1, combined$GLAT_1),
#                                     data = as.data.frame(combined),
#                                     proj4string = sp::CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
# 
# rgdal::writeOGR(obj = points,
#          dsn = paste0(outDir, out_filename, '.shp'),
#          layer = out_filename,
#          driver = 'ESRI Shapefile',
#          overwrite_layer = T,
#          verbose = T)

pdf(paste0(outDir, 'plots.pdf'))

  hist(x = combined$ZG_diff,
       breaks = seq(from = min(combined$ZG_diff)-0.1,
                    to = max(combined$ZG_diff)+0.1,
                    by = 0.05),
       plot=T, axes=T,
       main = paste("Histogram of all ZG_diff"))
  
  hist(x = combined$ZG_diff,
       breaks = seq(from = min(combined$ZG_diff)-0.1,
                    to = max(combined$ZG_diff)+0.1,
                    by = 0.02),
       plot=T, 
       axes=T,
       xlim = c(-1,1),
       main = paste("Histogram of all ZG_diff, zoom to (-1,1)"))
  
  # plot(x = combined$distance,
  #      y = combined$ZG_diff,
  #      xlab = 'distance between footprint centers',
  #      ylab = 'difference in detected ground elevation',
  #      cex=0.75,
  #      pch=16)
  # 
  # plot(x = combined$SENS_min,
  #      y = combined$ZG_diff,
  #      xlab = 'minmum Sensitivity',
  #      ylab = 'difference in detected ground elevation',
  #      cex=0.75,
  #      pch=16)
  # 
  # plot(x = combined$CHAN_ZG,
  #      y = combined$ZG_diff,
  #      xlab = 'Same (1) or different (0) ZG channel',
  #      ylab = 'difference in detected ground elevation',
  #      cex=0.75,
  #      pch=16)
  # 
  # plot(x = combined$INCI_max,
  #      y = combined$ZG_diff,
  #      xlab = 'Max incidence angle',
  #      ylab = 'difference in detected ground elevation',
  #      cex=0.75,
  #      pch=16)
  # 
  # plot(x = combined$ZG_1,
  #      y = combined$ZG_diff,
  #      xlab = 'Elevation',
  #      ylab = 'difference in detected ground elevation',
  #      cex=0.75,
  #      pch=16)
  # 
  # plot(x = combined$RH98_1,
  #      y = combined$ZG_diff,
  #      xlab = 'RH 98',
  #      ylab = 'difference in detected ground elevation',
  #      cex=0.75,
  #      pch=16)
  
  corrgram(combined %>% select(ZG_diff, 
                               COMP_max,
                               AZIM_max,
                               INCI_max,
                               CHAN_ZG,
                               distance,
                               RH10_1,
                               RH50_1,
                               RH98_1,
                               RH98_diff),
           order=NULL,
           lower.panel=panel.shade,
           text.panel=panel.txt,
           upper.panel=panel.pts, 
           main="Correlogram")
dev.off()

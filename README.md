# evaluate_LVISF2_repeat_consistency

The goal of this project was to determine if repeat measurements of LVIS airborne waveform lidar are sufficiently consistent to be used for change detection of soil burning in boreal forest wildfires. The removal of thick soil organic layer in the boreal forest is ecologically relevant for forest regeneration trajectory and also useful for estinmating carbon release due to the combustion of organic soil. If differences in ground elevation estimates of the same location are generally greater than 5cm over a period of time where no changes would be expected (such as consecutive days) then it will not be possible to use these data for detecting burn depth. Ultimately this was the case, and the project shifted direction to assess the use of UAV structure-from-motion data, rather than airborne waveform lidar data, to study burn depth resulting from a wildfire near Fairbanks, Alaska. 

The main workflow was run on NAU's computing cluster using R Code. An alternate Google Earth Engine workflow implements some of the same steps.

Steps:  
R workflow:
  (1) Download LVIS Level 1B Geolocated waveform data from the [NSIDC webpage](https://lvis.gsfc.nasa.gov/Data/Maps/ABoVE2017Map.html). 
  (2) Create .csv with metadata for LVIS flight pairs for comparison
  (3) Run ./code/array_for_join_points.sh script which runs the R script ./code/join_points_from_txt_filelist.R over each member of a job array.
  (4) Write output plots using ./code/summary.R

GEE workflow:
  (1) Download LVIS Level 1B Geolocated waveform data from the [NSIDC webpage](https://lvis.gsfc.nasa.gov/Data/Maps/ABoVE2017Map.html).   
  (2) use R script to convert TXT data to shapefile  
  (3) Upload shapefiles as a Google Earth Engine asset  
  (4) In GEE: compare repeat measurements of elevation at the footprint level to evaluate consistency of LVIS data.  
  GEE script link: https://code.earthengine.google.com/8055502adbe4f95787cf035db4dd12c5  

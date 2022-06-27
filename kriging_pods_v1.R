#------------------------------SettingsSet--------------------------------------
#tell the code where to find your pod database text files
directory = "C:\\Users\\okorn\\Documents\\NASA Postdoc\\R Kriging\\Kriging Example Dataset"

#(Eventually) select which algorithm to use
#This code shows all 6 for comparative reasons
  #1 - Nearest Neighbor
  #2 - Inverse Distance Weighting
  #3 - Thin Plate Spline Regression
  #4 - Generalized Additive Model
  #5 - Triangular Irregular Surface
  #6 - Automatized Kriging
#----------------------------End settingsSet------------------------------------

#Based on tutorial from here: https://swilke-geoscience.net/post/spatial_interpolation/

# We will need some packages for (spatial) data processing
library(tidyverse) # wrangling tabular data and plotting
library(readr)
library(sf) # processing spatial vector data
library(sp) # another vector data package necessary for continuity
library(raster) # processing spatial raster data. !!!overwrites dplyr::select!!!

# And a lot of different packages to test their interpolation functions
library(gstat)  # inverse distance weighted, Kriging
library(fields) # Thin Plate Spline
library(interp) # Triangulation
library(mgcv)   # Spatial GAM
library(automap)# Automatic approach to Kriging

# Finally, some packages to make pretty plots
library(patchwork)
library(viridis)

#Other packages we need
library(magrittr) # needs to be run every time you start R and want to use %>%
library(dplyr)    # alternatively, this also loads %>%
require(readxl) #for reading data from file
library(data.table) #for manipulating data from text files


### Step 0: Have some data

#set the working directory from which the files will be read from
setwd(directory)
#create a list of the files from your target directory
fileList <- list.files(path=directory)
#get our pod list by removing "_Field.txt" from the list of names
podList = gsub('_Field.txt', '', fileList)

#loop through each of the files & extract only the columns we need:
  #podName, time, estimate, latitude, & longitude only

#get the correct number of entries to loop through
x <- c(1:length(fileList))


for(i in x) {
  #designate the currentPod for analysis
  currentPod = podList[i]
  #get the file directory for the current pod
  pod_dir = paste(directory,fileList[i],sep = "\\")
  #load the file in for the current pod
  temp_data <- fread(pod_dir, select = c("podName","time","estimate","latitude","longitude"),data.table = TRUE)
  
  #create a name using the currentPod name
  target <- paste(currentPod, "data", sep = "_")
  #use this new name to label the pod data
  assign(target,temp_data)
  
}

#get a list with the names of our new datasets
data_names = gsub('_Field.txt', '_data', fileList) 
  
#Create a timetable equivalent to put all our data together
  #We will need to look at 1 timestamp at a time

#see how many files we're working with
#(must be at least 2 in order to krig / interpolate between them)

if (length(x) == 2) {
  #only need to merge these 2 files - one step
  left_join(df1, df2, by = 'timestamp')
}

for (k in x-1) {
  #will need to merge 2 at a time, so repeat x-1 times


}

#Now plot the raw data to get a visual idea of it
point_plot <- ggplot(
  data = pts_NH4,
  mapping = aes(x = X, y = Y, color = NH4)) +
  geom_point(size = 3) +
  scale_color_gradientn(colors = c("blue", "yellow", "red"))

point_plot


### Step 1: Create a grid template

#Option a) The simple approach
# First let's define a bounding Box, a rectangle that contains all our data
# points. This is one of the simpler ways to do it

bbox <- c(
  "xmin" = min(pts_NH4$X),
  "ymin" = min(pts_NH4$Y),
  "xmax" = max(pts_NH4$X),
  "ymax" = max(pts_NH4$Y)
)

grd_template <- expand.grid(
  X = seq(from = bbox["xmin"], to = bbox["xmax"], by = 20),
  Y = seq(from = bbox["ymin"], to = bbox["ymax"], by = 20) # 20 m resolution
)

grid_plot <- ggplot() +
  geom_point(data = grd_template, aes(x = X, y = Y), size = 0.01) +
  geom_point(data = pts_NH4,
             mapping = aes(x = X, y = Y, color = NH4), size = 3) +
  scale_color_gradientn(colors = c("blue", "yellow", "red")) +
  coord_cartesian( #zooming in so we can actually see something
    xlim = c(408000, 409000), ylim = c(5815000, 5816000)) +
  theme_bw()

grid_plot

#Option b) classic approach
sf_NH4 <- st_as_sf(pts_NH4, coords = c("X", "Y"), crs = 25833)

alt_grd_template_sf <- sf_NH4 %>% 
  st_bbox() %>% 
  st_as_sfc() %>% 
  st_make_grid(
    cellsize = c(20, 20),
    what = "centers"
  ) %>%
  st_as_sf() %>%
  cbind(., st_coordinates(.)) %>% 
  st_drop_geometry() %>% 
  mutate(Z = 0)


### Step 1b: Rasterizing your grid template

#Will use multiple interpolation methods from different packages
#so create a copy of our grid template in a different data structure

# {raster} expects a PROJ.4 string, see https://epsg.io/25833
crs_raster_format <- "+proj=utm +zone=33 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

grd_template_raster <- grd_template %>% 
  dplyr::mutate(Z = 0) %>% 
  raster::rasterFromXYZ( 
    crs = crs_raster_format)

## Warning in showSRID(uprojargs, format = "PROJ", multiline = "NO"): Discarded datum Unknown based on GRS80 ellipsoid in CRS definition,
##  but +towgs84= values preserved

# Let's also carry with us the raster from the alternative approach
alt_grd_template_raster <- alt_grd_template_sf %>% 
  raster::rasterFromXYZ(
    crs = crs_raster_format
  )

## Warning in showSRID(uprojargs, format = "PROJ", multiline = "NO"): Discarded datum Unknown based on GRS80 ellipsoid in CRS definition,
##  but +towgs84= values preserved

##This version is to make only 1 - we need to loop through for each point in time to make a gif of all outcomes

#Step 2: Fit a model
# We start with functions that return model objects as this is the most 
# common case

# Nearest Neighbor
fit_NN <- gstat::gstat( # using package {gstat} 
  formula = NH4 ~ 1,    # The column `NH4` is what we are interested in
  data = as(sf_NH4, "Spatial"), # using {sf} and converting to {sp}, which is expected
  nmax = 10, nmin = 3 # Number of neighboring observations used for the fit
)

# Inverse Distance Weighting
fit_IDW <- gstat::gstat( # The setup here is quite similar to NN
  formula = NH4 ~ 1,
  data = as(sf_NH4, "Spatial"),
  nmax = 10, nmin = 3,
  set = list(idp = 0.5) # inverse distance power
)

# Thin Plate Spline Regression
fit_TPS <- fields::Tps( # using {fields}
  x = as.matrix(pts_NH4[, c("X", "Y")]), # accepts points but expects them as matrix
  Y = pts_NH4$NH4,  # the dependent variable
  miles = FALSE     # EPSG 25833 is based in meters
)

# Generalized Additive Model
fit_GAM <- mgcv::gam( # using {mgcv}
  NH4 ~ s(X, Y),      # here come our X/Y/Z data - straightforward enough
  data = pts_NH4      # specify in which object the data is stored
)

# Next we use a couple of functions that have a slightly different modus
# operandi as they in fact already return interpolated Z values.

# Triangular Irregular Surface
fit_TIN <- interp::interp( # using {interp}
  x = pts_NH4$X,           # the function actually accepts coordinate vectors
  y = pts_NH4$Y,
  z = pts_NH4$NH4,
  xo = grd_template$X,     # here we already define the target grid
  yo = grd_template$Y,
  output = "points"
) %>% bind_cols()

# Automatized Kriging  
fit_KRIG <- automap::autoKrige(      # using {automap}
  formula = NH4 ~ 1,                 # The interface is similar to {gstat} but
  input_data = as(sf_NH4, "Spatial") # {automap} makes a lot of assumptions for you
) %>% 
  .$krige_output %>%  # the function returns a complex object with lot's of metainfo
  as.data.frame() %>% # we keep only the data we are interested in
  dplyr::select(X = x1, Y = x2, Z = var1.pred) 

## [using ordinary kriging]

#Interpolate!

#Case 1: X, Y, and Z are provided by the algorithm --------------------

# Triangular Irregular Surface
interp_TIN <- raster::rasterFromXYZ(fit_TIN, crs = crs_raster_format)

## Warning in showSRID(uprojargs, format = "PROJ", multiline = "NO"): Discarded datum Unknown based on GRS80 ellipsoid in CRS definition,
##  but +towgs84= values preserved

# Automatized Kriging
interp_KRIG <- raster::rasterFromXYZ(fit_KRIG, crs = crs_raster_format)

## Warning in showSRID(uprojargs, format = "PROJ", multiline = "NO"): Discarded datum Unknown based on GRS80 ellipsoid in CRS definition,
##  but +towgs84= values preserved

# Note that you can inspect the result of your interpolation anytime with a 
# simple `plot()`, if it is a {raster*} object
plot(interp_KRIG)

#End of case 1 -------------------------------------------------------

#Case 2: algorithm only returns a model object, so we use raster::interpolate ---

print(fit_TPS)

# Nearest Neighbor
interp_NN <- interpolate(grd_template_raster, fit_NN)

## [inverse distance weighted interpolation]

# Inverse Distance Weighting
interp_IDW <- interpolate(grd_template_raster, fit_IDW)

## [inverse distance weighted interpolation]

# Thin Plate Spline Regression
interp_TPS <- interpolate(grd_template_raster, fit_TPS)

#End of case 2 --------------------------------------------------------------

#Case 3: using the official method stats::predict ----------------------

# Generalized Additive Model
interp_GAM <- grd_template %>% 
  mutate(Z = predict(fit_GAM, .)) %>% 
  rasterFromXYZ(crs = crs_raster_format)
## Warning in showSRID(uprojargs, format = "PROJ", multiline = "NO"): Discarded datum Unknown based on GRS80 ellipsoid in CRS definition,
##  but +towgs84= values preserved

#End of case 3 ---------------------------------------------------------


#Final step: Visualization
plot_my_rasters <- function(raster_object, raster_name){
  
  df <- rasterToPoints(raster_object) %>% as_tibble()
  colnames(df) <- c("X", "Y", "Z")
  
  ggplot(df, aes(x = X, y = Y, fill = Z)) +
    geom_raster() +
    ggtitle(label = raster_name) +
    scale_fill_viridis(option = "C") +
    theme_bw() +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank()
    )
}

rasterlist <- list(
  "Nearest Neighbor" = interp_NN, 
  "Inverse Distance Weighted" = interp_IDW, 
  "Kriging" = interp_KRIG, 
  "Thin Plate Spline Regression" = interp_TPS,
  "Triangular Irregular Surface" = interp_TIN, 
  "Generalized Additive Model" = interp_GAM
)

plotlist <- map2(
  rasterlist,
  names(rasterlist),
  plot_my_rasters
)

# Note that the next trick only works because of library(patchwork)
(plotlist[[1]] + plotlist[[2]]) /
  (plotlist[[3]] + plotlist[[4]]) /
  (plotlist[[5]] + plotlist[[6]])


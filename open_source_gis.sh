

### Examples of open source GIS tools using the command line - Sami Domisch, Feb17


# ###============================================================================
# ### Install tools
# sudo apt-get update && sudo apt-get upgrade

# # GDAL/OGR
# # www.gdal.org
# sudo apt-get install gdal-bin  

# # GRASS
# sudo apt-get install grass-core  grass-dev  grass-gui

# # pktools
# wget "http://download.savannah.gnu.org/releases/pktools/install_pktools.sh" $HOME
# sudo bash install_pktools.sh

# ### openforis
# ### http://www.openforis.org/tools/geospatial-toolkit.html
# ### Installation
# sudo apt-get install gcc g++ gdal-bin libgdal1-dev libgsl0-dev libgsl0ldbl libproj-dev python-gdal python-scipy python-tk python-qt4 perl
# wget wget http://foris.fao.org/static/geospatialtoolkit/releases/OpenForisToolkit.run
# sudo chmod u+x OpenForisToolkit.run
# sudo ./OpenForisToolkit.Run

# ### FW tools
# wget "http://fwtools.loskot.net/FWTools-linux-2.0.6.tar.gz"
# tar xzvf FWTools-linux-2.0.6.tar.gz
# cd FWTools-2.0.6
# sudo ./install.sh
# ### If you use Bash as your shell add the following to your startup script (ie. ~/.bash_profile):
# PATH=$PATH:$HOME/FWTools/bin_safe

# # OpenEV
# # http://openev.sourceforge.net/
# wget "https://sourceforge.net/projects/openev/files/OpenEV/1.8.0/openev-linux-180.tar.gz/download"  -O openev-linux-180.tar.gz
# tar xvf openev-linux-180.tar.gz
# cd openev
# sudo ./install linux   /usr/bin/openev
# # add lias openev to the ~/.bashrc
# echo "alias openev='/usr/bin/openev/bin/openev' "  >>  ~/.bashrc
# ###============================================================================



### Basic bash commands
pwd 					# print working directory
echo HOME   			# print text string in the console 
echo $HOME 				# print the variable that is assigned to the string
cd $HOME 				# enter your home directory
ll 						# or "ls -hl", displays items in folder
top						# check ongoing processes, quit with "q"
seq 1 10 				# print numbers 1 to 10
seq 1 10 > test.txt  	# print numbers 1 to 10 and write into a .txt file
echo 11 >> test.txt 	# append to the same .txt file
head test.txt   		# show first 10 lines
cat test.txt  			# show the entire file
rm test.txt   			# remove the file
# [ctrl + c]  			# stop running process

# Basic for loop
for i in $(seq 1 15); do
echo Printing number $i
done

### Use this interactive shell to learn more
http://www.learnshell.org/


# ###============================================================================



### Create working directory
DIR=$HOME/gis_intro # assign a variable
rm -rf $DIR 		# remove previous data
mkdir -p $DIR 		# create the folder
cd $DIR

### Download global world borders shapefile
wget -O world.zip   "http://thematicmapping.org/downloads/TM_WORLD_BORDERS-0.3.zip"

### Unzip the data
unzip -j   world.zip  -d  $DIR/world


###==================================================================================


### GDAL & openev 
# www.gdal.org
# http://openev.sourceforge.net/

### Plot data
openev $DIR/world/TM_WORLD_BORDERS-0.3.shp &
### Check metadata
ogrinfo $DIR/world/TM_WORLD_BORDERS-0.3.shp  -al -so

### Rasterize the shapefile- test different spatial grains and file sizes
# man gdal_rasterize
gdal_rasterize -h
# e.g to 0.1 degree
### Check EPSG code from www.spatialreference.org, e.g. 
gdal_rasterize  $DIR/world/TM_WORLD_BORDERS-0.3.shp   -l TM_WORLD_BORDERS-0.3  $DIR/world/world.tif   -a_srs EPSG:4326   -a_nodata  -9999  -tr 0.1  0.1  -a UN 
openev   $DIR/world/TM_WORLD_BORDERS-0.3.shp     $DIR/world/world.tif  &
# e.g to 0.0083 degree, aka 1km --> check file size
gdal_rasterize  $DIR/world/TM_WORLD_BORDERS-0.3.shp   -l TM_WORLD_BORDERS-0.3  $DIR/world/world.tif   -a_srs EPSG:4326   -a_nodata  -9999  -tr 0.0083  0.0083  -a UN 

### See the different data types for storing 
### https://grass.osgeo.org/grass72/manuals/r.out.gdal.html

# Ranges of GDAL data types
#   GDAL data type	      	   minimum  		maximum
# ========================================================
#   Byte  			   	    		 0  	 	    255
#   UInt16						     0  		 65,535
#   Int16, CInt16 	  	   	   -32,768  		 32,767
#   UInt32				    		 0    4,294,967,295
#   Int32, CInt32 		-2,147,483,648    2,147,483,647
#   Float32, CFloat32	       -3.4E38  		 3.4E38
#   Float64, CFloat64	     -1.79E308         1.79E308

# 1km with the correct output type --> check file size
gdal_rasterize $DIR/world/TM_WORLD_BORDERS-0.3.shp    -l TM_WORLD_BORDERS-0.3  $DIR/world/world.tif   -a_srs EPSG:4326 -ot Byte   -a_nodata  255  -tr 0.0083  0.0083 -a UN  -co COMPRESS=LZW
gdalinfo  $DIR/world/world.tif  -mm


### Create a global 1-km "distance to the sea" layer
### http://www.gdal.org/gdal_proximity.html
gdal_proximity.py  $DIR/world/world.tif   -values 255   $DIR/distance.tif  -co COMPRESS=LZW
openev   $DIR/distance.tif  &

###==================================================================================


### GRASS
### - download 1km elevation data for central Europe
### - run hydrological conditioning of the DEM
### - extract stream network and watersheds

### Create and enter the folder where the data will be stored
export DIR=/home/domisch/gis_intro
mkdir $DIR  # create new folder
cd $DIR


### Download and unzip a DEM from WorldClim (http://www.worldclim.org/tiles.php):
wget -O  $DIR/alt_16_tif.zip  "http://biogeo.ucdavis.edu/data/climate/worldclim/1_4/tiles/cur/alt_16_tif.zip"
unzip  -o $DIR/alt_16_tif.zip  -d  $DIR/dem
gdalinfo $DIR/dem/alt_16.tif        # check data

### Create the GRASS GIS data base and enter GRASS:
grass72  -text -c -e  $DIR/dem/alt_16.tif  $DIR/grass_location
grass72 -text $DIR/grass_location/PERMANENT  # enter GRASS

### Read data into GRASS
r.in.gdal input=$DIR/dem/alt_16.tif    output=elevation --overwrite
r.info elevation # check data

### Open GUI and visualize the layers
g.gui wxpython

### Run a hydrologic conditioning on the DEM, removing sinks and peaks
# g.extension  extension=r.hydrodem # download extension
r.hydrodem  input=elevation  output=elevation_cond  --overwrite

### Extract drainage direction and stream network
# r.watershed  --h  # see help regarding the options and flags
r.watershed  elevation=elevation_cond  drainage=drainage   stream=stream  accumulation=accumulation  threshold=100  --o

### Get drainage basins (last downstream segment: -l flag)
# g.extension  extension=r.stream.basins
r.stream.basins  direction=drainage  stream_rast=stream  basins=basins   --o  #-l

### Categorize the single basins:
r.clump -d input=basins  output=basins_cat  --o
r.info basins_cat

### Write files to disk:
r.out.gdal  input=stream   output=$DIR/stream.tif  type=Int32  nodata=-9999  --o  -c  -m    createopt="COMPRESS=LZW,ZLEVEL=9"
r.out.gdal  input=basins_cat   output=$DIR/basin.tif  type=Int32  nodata=-9999  --o  -c  -m  createopt="COMPRESS=LZW,ZLEVEL=9"

exit
###==================================================================================


### Fix NoData and alignment issue in Worldclim data (here: the basin.tif layer originating from the DEM)
gdalwarp  -ot Int32  -dstnodata -9999  -te 0 30 30 60  -tap   -tr  0.00833333333333  0.00833333333333  $DIR/basin.tif  $DIR/basin_al.tif 
gdalinfo $DIR/basin_al.tif | grep Size
gdalinfo $DIR/basin.tif | grep Size

### pktools
# http://pktools.nongnu.org/html/index.html
### Crop global distance layer to same extent as basins (central Europe)
# pkcrop  -i $DIR/distance.tif   -o $DIR/distance_crop1.tif  -align  -ulx 2  -uly 60  -lrx 20  -lry  35
pkinfo -i $DIR/basin_al.tif -bb -dx -dy
pkcrop -i $DIR/distance.tif  $(pkinfo -i $DIR/basin_al.tif  -bb -dx -dy)  -o $DIR/distance_mask.tif  -co COMPRESS=LZW -co ZLEVEL=9
openev $DIR/distance_mask.tif &
# gdalinfo $DIR/distance_mask.tif 


# Merge two raster layers
pkcrop  -i $DIR/distance.tif   -o $DIR/distance_crop1.tif  -align  -ulx 2  -uly 60  -lrx 20  -lry  35
pkcrop  -i $DIR/distance.tif  -o $DIR/distance_crop2.tif  -align  -ulx 15  -uly 40  -lrx 35  -lry  20
pkcomposite   -i $DIR/distance_crop1.tif   -i $DIR/distance_crop2.tif   -o $DIR/distance_merge.tif
openev $DIR/distance_crop1.tif    $DIR/distance_crop2.tif   $DIR/distance_merge.tif   &


### Create a polygon from the basin raster layer
# gdal_polygonize.py  -h
gdal_polygonize.py -8  -f "ESRI Shapefile"  $DIR/basin.tif  $DIR/basin.shp   basin  basin_id
ogrinfo   $DIR/basin.shp   -al -so
openev $DIR/basin.shp  &


###==================================================================================

### Open Foris Geospatial Toolkit
# http://www.openforis.org/OFwiki/index.php/Tools_%26_Exercises
### Calculate distance statistics for each basin
oft-stat  -i $DIR/distance_mask.tif  -o $DIR/stats.txt -um  $DIR/basin_al.tif   -mm 
head $DIR/stats.txt 
### Add the header
# echo basin_id pixel min max avg std > $DIR/stats_new.txt 
# cat $DIR/stats.txt  >> $DIR/stats_new.txt 
# cat $DIR/stats_new.txt | head -20  


### Attach the average "distance to the sea" to the polygon attribute table 
# First get the two columns as a lookup table:    ID | avg
awk  '{ print $1, $5 }' $DIR/stats.txt > $DIR/stats_tmp.txt
head $DIR/stats_tmp.txt
# Add data
oft-addattr-new.py  $DIR/basin.shp  basin_id  avg  Float  $DIR/stats_tmp.txt 
ogrinfo   $DIR/basin.shp   -al -so
openev $DIR/basin.shp  &


###==================================================================================

### R
### - load raster and polygon data
### - set projection
### - subset / crop layers to smaller extent
### - extract maximum elevation for each watershed

R --vanilla 

### Load packages and set working directory
# install.packages("raster")
library(raster)
library(sp)
library(rgdal)
library(maptools)
library(rgeos)


DIR="/home/domisch/gis_intro"
setwd(DIR)

### Read global shapefile and set projection
world <- readShapePoly("world/TM_WORLD_BORDERS-0.3.shp")
world
### See also www.spatialreference.org
proj4string(world) <- "+proj=longlat +ellps=WGS84"
x11(); plot(world)

### Load the elevation data and basin-polygons
dem <- raster("dem/alt_16.tif")
dem

basin <- readShapePoly("basin.shp")
basin

### Define projection for the basins
proj4string(basin) <- proj4string(dem)

### Dissolve global shapefile
world_dissolve <-  gUnaryUnion(world)
x11(20,10); plot(world_dissolve)

### Subset global shapefile by country
head(world)
world_sub <- subset(world, NAME== "Austria")
x11(); plot(world_sub)

### Which basins overlap with the selection?
basin_subset <- raster::intersect(world_sub, basin)
x11(); plot(basin_subset)
head(basin_subset) # both shapefiles are merged

### Crop elevation raster data
dem_crop <- crop(dem, extent(basin_subset))
x11(20,10); plot(dem_crop); plot(basin_subset, add=T)

### Alternatives to cropping the shapefiles to a smaller extent
# basin_subset <- subset(basin, avg > 700) 		# subset by distance to sea

# e <- extent(10, 13, 45, 48) 		# 3x3 degree window
# basin_subset <- crop(basin, e) 		# subset by extent
# basin_subset <- subset(basin, basin_id > 30000 & basin_id < 31000) 		# subset by the basin ID
# x11(); plot(dem); plot(basin_subset, add=T)


### Get the maximum elevation in each watershed

### Run extraction on (a subset of the) shapefile
beginCluster(8)
basin_subset_elevation <- extract(dem_crop, basin_subset, fun=max, sp=T) # only for a subset...
endCluster()

### Check data
head(basin_subset_elevation)

### Export the shapefile
writePolyShape(basin_subset_elevation, "basin_subset_elevation.shp")


### Run same analysis on raster layers
basin_r <- raster("basin.tif")
basin_r_crop <- crop(basin_r, dem_crop)

### Plot
x11(40,15); par(mfrow=c(1, 2))
plot(dem_crop); plot(basin_subset_elevation, add=T)
plot(basin_r_crop); plot(basin_subset_elevation, add=T)

### Get maximum elevation within each watershed
dem_max <- zonal(dem_crop, basin_r_crop, "max", na.rm=TRUE)
head(dem_max)

### Reclassify the basin -raster
dem_max <- as.data.frame(dem_max) 	 # convert matrix to dataframe 
names(dem_max) <- c("is", "becomes") # change headers, see ?reclassify
head(dem_max)
basin_maxElev <- reclassify(basin_r_crop, dem_max)

x11(); plot(basin_maxElev)
plot(basin_subset_elevation, add=T)

# ### Run using the cluster-function:
# beginCluster(4)
# myfun <- function(x) reclassify(x, dem_max)
# basin_maxElev <- clusterR(basin_r_crop, myfun, export="dem_max")
# endCluster()

### Export the raster
writeRaster(basin_maxElev, "elevation_max.tif", overwrite=T)
# graphics.off()

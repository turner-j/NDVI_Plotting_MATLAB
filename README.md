# Plotting NDVI from MODIS in MATLAB

These scripts will plot NDVI from a single MODIS file in MATLAB and create a time series of NDVI within a specific area.

## Description

These scripts were specifically made to process the Terra Moderate Resolution Imaging Spectroradiometer (MODIS) Vegetation Indices (MOD13Q1) Version 6 data. Composite images of NDVI are generated every 16 days at 250 meter (m) spatial resolution. However, the code can be edited to plot other variables or resolutions.

## Getting Started

### Dependencies

* Download ["kml2struct"](https://www.mathworks.com/matlabcentral/fileexchange/35642-kml2struct) and add its folder to the MATLAB path or copy the file into your working directory.

* Download the HDF4 files of interest from [USGS Earth Explorer](https://earthexplorer.usgs.gov/).

* Create a kml file of your area of interest. This can be done with Google Earth. Save it in your working directory as "MikeIslandPolygon" or alter the file name in the script.

### Installing

* Install the program into your working directory in MATLAB.
* Create a folder in this directory titled "MOD13Q1" containing the MODIS HDF4 files.

### Executing program

* The program can be ran from start to finish.

## Help

* You will need to change the sinusoidal projection constants for grid sizes other than 250 meters. These coordinates can be found in Appendix B of the [MODIS User Guide (Giglio et al, 2021)](https://modis-land.gsfc.nasa.gov/pdf/MODIS_C61_BA_User_Guide_1.0.pdf).

* Changes may need to be made to file paths, folder, and variable names.

* Issues plotting time series data will arise if more than one file has the same start or end date in its file name.


## Author

Jessica Richardson, Ph.D.


## Acknowledgments

This code was written with the help of ChatGPT.

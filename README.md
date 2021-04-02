# MODIS Greenland

This repository contains a script for georeferencing the Multilayer Greenland Ice Surface Temperature, Surface Albedo,
and Water Vapor from MODIS, Version 1 [dataset](https://nsidc.org/data/MODGRNLD/versions/1) for use
in QGIS or other GIS software.

## Steps:
### 1. Download/clone this repository

### 2. Create a conda environment

Once you have successfully cloned the repository, create a `conda` environment using the `environment.yml` file:

```
> conda env create -f environment.yml
```

### 3. Run the script

Navigate to the folder(s) where you have the NetCDF files stored. You can run the script like so:

```
python path_to_repository/copy_netcdf.py
```

The script will copy each of the `MODGRNLD*.nc` files to a new folder, `georeferenced`, while updating the georeferencing information for each file.

### 4. Next steps
When the script has finished, you can load the newly georeferenced files into QGIS, or your GIS software of choice.
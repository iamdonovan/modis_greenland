import os
from glob import glob
import errno
import time
from collections import OrderedDict
from osgeo import osr
import numpy as np
import xarray as xr


def parse_wkt(sref_wkt):
    split1 = sref_wkt.split(',PROJECTION')
    split2 = split1[1].split('],')
    split3 = split1[0].split(',GEOGCS[')[1]

    inverse_flattening = float(split3.split(',SPHEROID')[1].split(',')[2])
    return inverse_flattening


def mkdir_p(out_dir):
    """
    Add bash mkdir -p functionality to os.makedirs.

    :param out_dir: directory to create.
    """
    try:
        os.makedirs(out_dir)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(out_dir):
            pass
        else:
            raise


def create_crs_variable(ds):
    """
    Create a CRS variable for a netCDF file, given an input xarray Dataset with projection information

    :param ds: input xarray DataSet
    :type ds: xarray.core.dataset.DataSet

    :returns odict: An OrderedDictionary that can be used to create the crs variable
    """
    for varname in ds.keys():
        if 'PROJECTION_GRID_MAPPING_NAME' in ds[varname].attrs:
            break
        else:
            continue
        raise KeyError

    grid_mapping_name = 'polar_stereographic'
    straight_vertical_longitude_from_pole = ds[varname].attrs['PROJECTION_STRAIGHT_VERTICAL_LONGITUDE_FROM_POLE']
    false_easting = ds[varname].attrs['PROJECTION_FALSE_EASTING']
    false_northing = ds[varname].attrs['PROJECTION_FALSE_NORTHING']
    latitude_of_projection_origin = ds[varname].attrs['PROJECTION_LATITUDE_OF_PROJECTION_ORIGIN']
    standard_parallel = ds[varname].attrs['PROJECTION_STANDARD_PARALLEL']
    long_name = "CRS definition"
    longitude_of_prime_meridian = ds[varname].attrs['PROJECTION_LONGITUDE_OF_PROJECTION_ORIGIN']
    semi_major_axis = ds[varname].attrs['PROJECTION_SEMIMAJOR_RADIUS']
    # semi_minor_axis = ds[varname].attrs['PROJECTION_SEMIMINOR_RADIUS']
    # scale_factor = ds[varname].attrs['PROJECTION_SCALING_FACTOR']

    sref = osr.SpatialReference()
    sref.ImportFromEPSG(3413)

    sref_wkt = sref.ExportToWkt()

    inverse_flattening = parse_wkt(sref_wkt)

    odict = OrderedDict({'grid_mapping_name': grid_mapping_name,
                         'straight_vertical_longtiude_from_pole': straight_vertical_longitude_from_pole,
                         'false_easting': false_easting, 'false_northing': false_northing,
                         'latitude_of_projection_origin': latitude_of_projection_origin,
                         'standard_parallel': standard_parallel, 'long_name': long_name,
                         'longitude_of_prime_meridian': longitude_of_prime_meridian,
                         'semi_major_axis': semi_major_axis, 'inverse_flattening': inverse_flattening,
                         # 'semi_minor_axis': semi_minor_axis,
                         'spatial_ref': sref_wkt})

    return odict


def copy_nc(filename, outdir='georeferenced'):
    mkdir_p(outdir)
    outfilename = os.path.join(outdir, filename)

    # load the dataset to copy
    in_ds = xr.open_dataset(filename)
    out_ds = in_ds.copy()

    out_ds.attrs['Conventions'] = 'CF-1.6'
    out_ds.attrs['description'] = "Multilayer Greenland Ice Surface Temperature, Surface Albedo," + \
                                  " and Water Vapor from MODIS, Version 1"
    out_ds.attrs['history'] = "Created " + time.ctime(time.time())
    out_ds.attrs['source'] = "National Snow and Ice Data Center"
    out_ds.attrs['url'] = "https://nsidc.org/data/MODGRNLD/versions/1"
    out_ds.attrs['contributors'] = "Hall, D. K. and N. DiGirolamo."

    # make sure to tell netcdf how to map these variables
    for varname in in_ds.keys():
        out_ds[varname].attrs['grid_mapping'] = 'polar_stereographic'

    # create the new crs dataset
    out_ds['polar_stereographic'] = xr.Variable((), np.array(b'', dtype='|S1'))
    out_ds['polar_stereographic'].attrs = create_crs_variable(in_ds)

    pixel_size = out_ds['Albedo'].attrs['Pixel_Size']
    ulx, uly = in_ds['Albedo'].attrs['Proj_ul_xy']

    ncols, = out_ds['x'].shape
    nrows, = out_ds['y'].shape

    # set x,y to the correct mapped values
    out_ds['x'] = ulx + pixel_size * (np.arange(ncols) + 0.5)
    out_ds['y'] = uly + -pixel_size * (np.arange(nrows) + 0.5)

    out_ds.assign_coords({'x': out_ds['x'],
                          'y': out_ds['y']})

    out_ds['x'].attrs = {'standard_name': 'projection_x_coordinate',
                         'long_name': 'x coordinate of projection',
                         'units': 'm'}

    out_ds['y'].attrs = {'standard_name': 'projection_y_coordinate',
                         'long_name': 'y coordinate of projection',
                         'units': 'm'}

    # set the geotransform for completeness
    geotrans = [ulx, pixel_size, 0, uly, 0, -pixel_size]
    out_ds['polar_stereographic'].attrs['GeoTransform'] = ' '.join([str(i) for i in geotrans])

    # remove the extra variables that we no longer need
    for_removal = ['Pixel_Size', 'Add_Offset', 'Scale_Factor', 'Proj_Type', 'Proj_Params',
                   'Proj_ul_xy', 'Proj_lr_xy']
    for varname in in_ds.keys():
        for attr in list(out_ds[varname].attrs.keys()):
            if attr in for_removal or 'PROJECTION' in attr:
                del out_ds[varname].attrs[attr]

    # save the file
    out_ds.to_netcdf(outfilename)


# -----------------------------
# get a list of MODGRNLD netCDF files
filelist = glob('MODGRNLD*.nc')
filelist.sort()

# default copies files to a new folder, georeferenced
print('copying {} files to georeferenced/'.format(len(filelist)))

for fn in filelist:
    print(fn)
    copy_nc(fn)

print('finished.')

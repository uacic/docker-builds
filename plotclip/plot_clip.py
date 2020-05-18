#!/usr/bin/env python3
"""
Author : Emmanuel Gonzalez
References: TERRA-REF & AgPipeline
Date   : 2020-05-17
Purpose: Clip plots using a shapefile/GeoJSON
"""
import argparse
from terrautils.spatial import clip_raster, \
geometry_to_geojson, convert_json_geometry, convert_geometry
from terrautils.imagefile import image_get_geobounds, get_epsg
from terrautils.betydb import get_site_boundaries
import geopandas as gpd
from osgeo import ogr, gdal
import numpy as np
import liblas
import subprocess
import re
import osr
import logging
import yaml
import json
import glob
import datetime
import os
from typing import Optional
import copy
#mport configuration

# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description='Rock the Casbah',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('dir',
                        metavar='dir',
                        type=str,
                        help='Input directory')

    parser.add_argument('-sen',
                        '--sensor',
                        metavar='sensor',
                        type=str,
                        required=True,
                        help='Sensor used to collect images')

    parser.add_argument('-e',
                        '--epsg',
                        metavar='epsg',
                        type=str,
                        default='4326',
                        help='Coordinate system of images')

    parser.add_argument('-shp',
                        '--shape',
                        metavar='GeoJSON',
                        type=str,
                        required=True,
                        help='GeoJSON with plot polygons')

    parser.add_argument('-o',
                        '--outdir',
                        help='Output directory',
                        metavar='outdir',
                        type=str,
                        default='plotclip_out')

    return parser.parse_args()


# --------------------------------------------------
def get_files_to_process(file_list: list, sensor: str, default_epsg: int = None) -> dict:
    """Returns a dictionary of georeferenced files to process
    Arguments:
        file_list: the list of file paths to process
        sensor: the name of the sensor associated with the files
        default_epsg: the default EPSG value to use if a file is missing one
    Return:
        Returns a dictionary with the file names as keys. Each key's value is another dictionary containing
        the file path, file bounds (as GeoJSON), and the sensor name
    """
    files_to_process = {}
    for one_file in file_list:
        filename = os.path.basename(one_file)
        if filename in files_to_process:
            continue
        if not os.path.exists(one_file):
            logging.warning("Skipping file that does not exist: '%s'", one_file)
            continue

        if one_file.endswith('.tif'):
            files_to_process[filename] = {
                'path': one_file,
                'bounds': get_image_bounds_json(one_file, default_epsg),
                'sensor_name': sensor
            }
        elif one_file.endswith(".las"):
            files_to_process[filename] = {
                'path': one_file,
                'bounds': get_las_extents(one_file, default_epsg),
                'sensor_name': sensor
            }
    return files_to_process


# --------------------------------------------------
def get_las_epsg_from_header(header: liblas.header.Header) -> str:
    """Returns the found EPSG code from the LAS header
    Arguments:
        header: the loaded LAS header to find the SRID in
    Return:
        Returns the SRID as a string if found, None is returned otherwise
    """
    epsg = None
    search_terms_ordered = ['DATUM', 'AUTHORITY', '"EPSG"', ',']
    try:
        # Get the WKT from the header, find the DATUM, then finally the EPSG code
        srs = header.get_srs()
        wkt = srs.get_wkt().decode('UTF-8')
        idx = -1
        for term in search_terms_ordered:
            idx = wkt.find(term)
            if idx < 0:
                break
        if idx >= 0:
            epsg = re.search(r'\d+', wkt[idx:])[0]
    except Exception as ex:
        logging.debug("Unable to find EPSG in LAS file header")
        logging.debug("    exception caught: %s", str(ex))

    return epsg


# --------------------------------------------------
def get_las_extents(file_path: str, default_epsg: int = None) -> Optional[str]:
    """Calculate the extent of the given las file and return as GeoJSON.
    Arguments:
        file_path: path to the file from which to load the bounds
        default_epsg: the default EPSG to assume if a file has a boundary but not a coordinate system
    Return:
        Returns the JSON representing the image boundary, or None if the
        bounds could not be loaded
    Notes:
        If a file doesn't have a coordinate system and a default epsg is specified, the
        return JSON will use the default_epsg.
        If a file doesn't have a coordinate system and there isn't a default epsg specified, the boundary
        of the image is not returned (None) and a warning is logged.
    """
    # Get the bounds and the EPSG code
    las_info = liblas.file.File(file_path, mode='r')
    min_bound = las_info.header.min
    max_bound = las_info.header.max
    epsg = get_las_epsg_from_header(las_info.header)
    if epsg is None:
        if default_epsg is not None:
            epsg = default_epsg
        else:
            logging.warning("Unable to find EPSG and not default is specified for file '%s'", file_path)
            return None

    ring = ogr.Geometry(ogr.wkbLinearRing)
    ring.AddPoint(min_bound[1], min_bound[0])  # Upper left
    ring.AddPoint(min_bound[1], max_bound[0])  # Upper right
    ring.AddPoint(max_bound[1], max_bound[0])  # lower right
    ring.AddPoint(max_bound[1], min_bound[0])  # lower left
    ring.AddPoint(min_bound[1], min_bound[0])  # Closing the polygon

    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)

    ref_sys = osr.SpatialReference()
    if ref_sys.ImportFromEPSG(int(epsg)) == ogr.OGRERR_NONE:
        poly.AssignSpatialReference(ref_sys)
        return geometry_to_geojson(poly)

    logging.error("Failed to import EPSG %s for las file %s", str(epsg), file_path)
    return None


# --------------------------------------------------
def clip_las(las_path: str, clip_tuple: tuple, out_path: str) -> None:
    """Clip LAS file to polygon.
    Arguments:
        las_path: path to point cloud file
        clip_tuple: tuple containing (minX, maxX, minY, maxY) of clip bounds
        out_path: output file to write
    Notes:
        The clip_tuple is assumed to be in the correct coordinate system for the point cloud file
    """
    bounds_str = "([%s, %s], [%s, %s])" % (clip_tuple[0], clip_tuple[1], clip_tuple[2], clip_tuple[3])

    pdal_dtm = out_path.replace(".las", "_dtm.json")
    with open(pdal_dtm, 'w') as dtm:
        dtm_data = """{
            "pipeline": [
                "%s",
                {
                    "type": "filters.crop",
                    "bounds": "%s"
                },
                {
                    "type": "writers.las",
                    "filename": "%s"
                }
            ]
        }""" % (las_path, bounds_str, out_path)
        logging.debug("Writing dtm file contents: %s", str(dtm_data))
        dtm.write(dtm_data)

    cmd = 'pdal pipeline "%s"' % pdal_dtm
    logging.debug("Running pipeline command: %s", cmd)
    subprocess.call([cmd], shell=True)
    os.remove(pdal_dtm)


# --------------------------------------------------
def get_image_bounds_json(file_path: str, default_epsg: int = None) -> Optional[str]:
    """Loads the boundaries of the image file and returns the GeoJSON
       representing the bounds (including EPSG code)
    Arguments:
        file_path: path to the file from which to load the bounds
        default_epsg: the default EPSG to assume if a file has a boundary but not a coordinate system
    Return:
        Returns the JSON representing the image boundary, or None if the
        bounds could not be loaded
    Notes:
        If a file doesn't have a coordinate system and a default epsg is specified, the
        return JSON will use the default_epsg.
        If a file doesn't have a coordinate system and there isn't a default epsg specified, the boundary
        of the image is not returned (None) and a warning is logged.
    """
    # Get the bounds (if they exist)
    bounds = image_get_geobounds(file_path)
    if bounds[0] == np.nan:
        return None

    epsg = get_epsg(file_path)
    if epsg is None:
        if default_epsg:
            epsg = default_epsg
        else:
            logging.warning("Files does not have a coordinate system defined and no default was specified: '%s'",
                            file_path)
            return None

    ring = ogr.Geometry(ogr.wkbLinearRing)
    ring.AddPoint(bounds[2], bounds[1])  # Upper left
    ring.AddPoint(bounds[3], bounds[1])  # Upper right
    ring.AddPoint(bounds[3], bounds[0])  # lower right
    ring.AddPoint(bounds[2], bounds[0])  # lower left
    ring.AddPoint(bounds[2], bounds[1])  # Closing the polygon

    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)

    ref_sys = osr.SpatialReference()
    if ref_sys.ImportFromEPSG(int(epsg)) == ogr.OGRERR_NONE:
        poly.AssignSpatialReference(ref_sys)
        return geometry_to_geojson(poly)

    logging.error("Failed to import EPSG %s for image file %s", str(epsg), file_path)
    return None


# --------------------------------------------------
def find_plots_intersect_boundingbox(bounding_box, all_plots, fullmac=True):
    """Take a list of plots from BETY and return only those overlapping bounding box.

    fullmac -- only include full plots (omit KSU, omit E W partial plots)

    """
    bbox_poly = ogr.CreateGeometryFromJson(str(bounding_box))
    bb_sr = bbox_poly.GetSpatialReference()
    intersecting_plots = dict()

#     for plotname in all_plots:
# #         if fullmac and (plotname.find("KSU") > -1 or plotname.endswith(" E") or plotname.endswith(" W")):
# #             continue
    for index, row in all_plots.iterrows():
        bounds = str(row.geometry)
        #bounds = all_plots[plotname]

        yaml_bounds = yaml.safe_load(bounds)
        #current_poly = ogr.CreateGeometryFromJson(json.dumps(yaml_bounds))
        current_poly = ogr.CreateGeometryFromWkt(yaml_bounds)
        # Check for a need to convert coordinate systems
        check_poly = current_poly
        if bb_sr:
            poly_sr = current_poly.GetSpatialReference()
            if poly_sr and not bb_sr.IsSame(poly_sr):
                # We need to convert to the same coordinate system before an intersection
                check_poly = convert_geometry(current_poly, bb_sr)
                transform = osr.CreateCoordinateTransformation(poly_sr, bb_sr)
                new_poly = current_poly.Clone()
                if new_poly:
                    new_poly.Transform(transform)
                    check_poly = new_poly

        intersection_with_bounding_box = bbox_poly.Intersection(check_poly)

        if intersection_with_bounding_box is not None:
            intersection = json.loads(intersection_with_bounding_box.ExportToJson())
            if 'coordinates' in intersection and len(intersection['coordinates']) > 0:
                intersecting_plots[row.ID] = bounds

    return intersecting_plots



# --------------------------------------------------
def get_spatial_reference_from_json(geojson: str) -> Optional[osr.SpatialReference]:
    """Returns the spatial reference embedded in the geojson.
    Args:
        geojson(str): the geojson to get the spatial reference from
    Return:
        The osr.SpatialReference that represents the geographic coordinate system
        in the geojson. None is returned if a spatial reference isn't found
    """
    yaml_geom = yaml.safe_load(geojson)
    current_geom = ogr.CreateGeometryFromJson(json.dumps(yaml_geom))

    if current_geom:
        return current_geom.GetSpatialReference()
    return None


# --------------------------------------------------
def calculate_overlap_percent(check_bounds: str, bounding_box: str) -> float:
    """Calculates and returns the percentage overlap between the two boundaries.
       The calculation determines the overlap shape between the two parameters and
       then calculates the percentage by dividing the overlap area by the bounding
       box area, and returns that value.
    Args:
        check_bounds: GeoJSON of boundary to check
        bounding_box: GeoJSON of boundary to check against
    Return:
        The calculated overlap percent (0.0 - 1.0) or 0.0 if there is no overlap.
        If an exception is detected, a warning message is logged and 0.0 is returned.
    """
    try:
        check_poly = ogr.CreateGeometryFromJson(str(check_bounds))
        bbox_poly = ogr.CreateGeometryFromJson(str(bounding_box))

        if check_poly and bbox_poly:
            intersection = bbox_poly.Intersection(check_poly)
            if intersection:
                return intersection.Area() / check_poly.Area()
    except Exception as ex:
        logging.warning("Exception caught while calculating shape overlap: %s", str(ex))

    return 0.0


# --------------------------------------------------
def geojson_to_tuples(bounding_box: str) -> tuple:
    """Returns the bounds of the shape
    Arguments:
        bounding_box: the JSON of the geometry
    Return:
        A tuple containing the bounds in (min Y, max Y, min X, max X) order
    """

    #yaml_geom = yaml.safe_load(bounding_box)
    current_geom = ogr.CreateGeometryFromJson(bounding_box)
    current_env = current_geom.GetEnvelope()

    return current_env[2], current_env[3], current_env[0], current_env[1]


# --------------------------------------------------
def to_tuples(bounding_box: str) -> tuple:
    """Returns the bounds of the shape
    Arguments:
        bounding_box: the JSON of the geometry
    Return:
        A tuple containing the bounds in (min Y, max Y, min X, max X) order
    """

    #yaml_geom = yaml.safe_load(bounding_box)
    current_geom = ogr.CreateGeometryFromWkt(bounding_box)
    current_env = current_geom.GetEnvelope()

    return current_env[2], current_env[3], current_env[0], current_env[1]


# --------------------------------------------------
def cleanup_request_md(source_md: dict) -> dict:
    """Makes a copy of the source metadata and cleans it up for use as plot-level information
    Arguments:
        source_md: the source metadata to clone and clean up
    Returns:
        returns the cleaned up metadata
    """
    if not source_md:
        return {}

    new_md = copy.deepcopy(source_md)
    new_md.pop('list_files', None)
    new_md.pop('context_md', None)
    new_md.pop('working_folder', None)


# --------------------------------------------------
def clip_raster_intersection(file_path: str, file_bounds: str, plot_bounds: str, out_file: str) -> Optional[int]:
    """Clips the raster to the intersection of the file bounds and plot bounds
    Arguments:
        file_path: the path to the source file
        file_bounds: the geometric boundary of the source file as JSON
        plot_bounds: the geometric boundary of the plot to clip to as JSON
        out_file: the path to store the clipped image
    Return:
        The number of pixels in the new image, or None if no pixels were saved
    Notes:
        Assumes the boundaries are in the same coordinate system
    Exceptions:
        Raises RuntimeError if the polygons are invalid
    """
    logging.debug("Clip to intersect of plot boundary: File: '%s' '%s' Plot: '%s'", file_path, str(file_bounds), str(plot_bounds))
    try:
        file_poly = ogr.CreateGeometryFromJson(str(file_bounds))
        plot_poly = ogr.CreateGeometryFromWkt(str(plot_bounds))

        if not file_poly or not plot_poly:
            logging.error("Invalid polygon specified for clip_raster_intersection: File: '%s' plot: '%s'",
                          str(file_bounds), str(plot_bounds))
            raise RuntimeError("One or more invalid polygons specified when clipping raster")

        intersection = file_poly.Intersection(plot_poly)
        if not intersection or not intersection.Area():
            logging.info("File does not intersect plot boundary: %s", file_path)
            return None

        # Make sure we pass a multipolygon down to the tuple converter
        if intersection.GetGeometryName().startswith('MULTI'):
            multi_polygon = intersection
        else:
            multi_polygon = ogr.Geometry(ogr.wkbMultiPolygon)
            multi_polygon.AddGeometry(intersection)

        # Proceed to clip to the intersection
        tuples = geojson_to_tuples(geometry_to_geojson(multi_polygon))
        return clip_raster(file_path, tuples, out_path=out_file, compress=True)

    except Exception as ex:
        logging.exception("Exception caught while clipping image to plot intersection")
        raise ex


# --------------------------------------------------
def main():
    args = get_args()
    #file_list = glob.glob('/mnt/c/Users/emman/Documents/psII_dimen_corr_test/*.tif', recursive=True)
    file_list = glob.glob(f'{args.dir}/*.tif', recursive=True) + glob.glob(f'{args.dir}/*.las', recursive=True)
    print(file_list)
    #sensor = 'ps2Top'
    sensor = args.sensor
    #default_epsg = '4326'
    default_epsg = args.epsg

    processed_files = 0
    processed_plots = 0
    start_timestamp = datetime.datetime.now()
    files_to_process = get_files_to_process(file_list, sensor, default_epsg)
    logging.info("Found %s files to process", str(len(files_to_process)))

    if files_to_process:

        # Get all the possible plots
        all_plots = gpd.read_file(args.shape)

        print(f'Have {len(all_plots):,} plots for site.')

        for filename in files_to_process:
                processed_files += 1
                file_path = files_to_process[filename]['path']
                file_bounds = files_to_process[filename]['bounds']
                sensor = files_to_process[filename]['sensor_name']
                #print(file_bounds)

                overlap_plots = find_plots_intersect_boundingbox(file_bounds, all_plots)
                file_spatial_ref = get_spatial_reference_from_json(file_bounds)

                for plot_name in overlap_plots:
                    print(plot_name)
                    processed_plots += 1
                    plot_bounds = convert_json_geometry(overlap_plots[plot_name], file_spatial_ref)
                    tuples = to_tuples(plot_bounds)

                    if filename.endswith('.tif'):
                        out_dir = args.outdir
                        out_path = os.path.join(out_dir, plot_name)

                        if not os.path.isdir(out_path):
                            os.makedirs(out_path)

                        out_file = os.path.join(out_path, filename)
                        clip_raster_intersection(file_path, file_bounds, plot_bounds, out_file)

                    elif filename.endswith('.las'):
                        out_dir = args.outdir
                        out_path = os.path.join(out_dir, plot_name)

                        if not os.path.isdir(out_path):
                            os.makedirs(out_path)
                        out_file = os.path.join(out_path, filename)

                        clip_las(file_path, tuples, out_path=out_file)


# --------------------------------------------------
if __name__ == '__main__':
    main()
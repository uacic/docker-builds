#!/usr/bin/env python3
"""
Author : emmanuelgonzalez
Date   : 2020-04-03
Purpose: Add latitude and longitude (from GDAL info) to EXIF metadata
"""

import argparse
import os
import sys
import gdal
import pprint
from GPSPhoto import gpsphoto
import glob


# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description='Rock the Casbah',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('dir',
                        metavar='str',
                        type=str,
                        help='Directory where images are located')

    parser.add_argument('-o',
                        '--outdir',
                        help='Output filename',
                        metavar='str',
                        type=str,
                        default='addexif_out/')

    return parser.parse_args()


# --------------------------------------------------
def main():
    """Get center coordinate using GDAL, then add to EXIF metadata"""

    args = get_args()

    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)

    # Scan for tif images in input directory
    images = glob.glob(args.dir + "*.tif", recursive=True)
    
    num = 0
    for i in images:
        num += 1
        ds = gdal.Open(i)
        meta = gdal.Info(ds)
        coord_list = []
        lines = meta.splitlines()

        for line in lines:
            if 'Center' in line:
                location = ' '.join(line.split()[:1]).strip('()')
                lat_dec = ' '.join(line.split()[2:3]).strip('()')
                long_dec =  ' '.join(line.split()[1:2]).strip('(),')
                print(f'{num}: ' + i + "\n" + 'Lat, Long: ' + f'({lat_dec}, {long_dec})' + "\n")

        filename = os.path.splitext(os.path.basename(i))[0]
        photo = gpsphoto.GPSPhoto(i)
        info = gpsphoto.GPSInfo((float(lat_dec), float(long_dec)))
        photo.modGPSData(info, args.outdir + filename + '_exif.tif')

    print(f'Done, images saved in {args.outdir}')


# --------------------------------------------------
if __name__ == '__main__':
    main()

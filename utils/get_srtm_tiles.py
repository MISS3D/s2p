#!/usr/bin/env python
# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>
# Copyright (C) 2018, David Youssefi <david.youssefi@cnes.fr>
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import argparse
import json
import os
import zipfile
import copy
import numpy as np

from osgeo import gdal

import s2p
from s2plib import common
from s2plib import rpc_utils
from s2plib import rpc_model

# This function is here as a workaround to python bug #24313 When
# using python3, json does not know how to serialize numpy.int64 on
# some platform numpy also decides to go for int64 when numpy.arange
# is called. This results in our json not being serializable anymore
# Calling json.dump(..,default=workaround_json_int64) fixes this
# https://bugs.python.org/issue24313
def workaround_json_int64(o):
    if isinstance(o,np.integer) : return int(o)
    raise TypeError

def download_srtm_tile(srtm_tile, out_dir,
                       srtm_url='http://data_public:GDdci@data.cgiar-csi.org/srtm/tiles/GeoTIFF'):
    """
    Downloads and extract an srtm tile from the internet.

    Args:
    srtm_tile: string following the pattern 'srtm_%02d_%02d', identifying
    the desired strm tile
    out_dir: directory where to store and extract the srtm tiles
    """
    # check if the tile is already there
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    if os.path.exists(os.path.join(out_dir, '%s.tif' % srtm_tile)):
        return

    # download the zip file
    srtm_tile_url = '%s/%s.zip' % (srtm_url, srtm_tile)
    zip_path = os.path.join(out_dir, '%s.zip' % srtm_tile)
    common.download(zip_path, srtm_tile_url)

    # extract the tif file
    if zipfile.is_zipfile(zip_path):
        z = zipfile.ZipFile(zip_path, 'r')
        z.extract('%s.tif' % srtm_tile, out_dir)
    else:
        print("%s not available" % srtm_tile)

    # remove the zip file
    os.remove(zip_path)
    return os.path.join(out_dir, '%s.tif' % srtm_tile)


def list_srtm_tiles(rpcfile, x, y, w, h):
    """
    Tells which srtm tiles are needed to cover a given region of interest.

    Args:
    rpcfile: path to the xml file describing the rpc model.
    x, y, w, h: four integers defining a rectangular region of interest
    (ROI) in the image. (x, y) is the top-left corner, and (w, h) are
    the dimensions of the rectangle.

    Returns:
    the set of srtm tiles corresponding to the input ROI.
    """
    rpc = rpc_model.RPCModel(rpcfile)
    lon_min, lon_max, lat_min, lat_max = rpc_utils.geodesic_bounding_box(rpc,
                                                                         x, y,
                                                                         w, h)
    srtm_list = []

    for lon_inc in np.arange(lon_min, lon_max, 5):
        for lat_inc in np.arange(lat_min, lat_max, 5):
            lat = min(lat_inc, 60)
            lat = max(lat_inc, -60)

            # tiles longitude indexes go from 1 to 72, covering the range from -180 to
            # +180
            tlon = (1+np.floor((lon_inc+180)/5)) % 72
            tlon = 72 if tlon == 0 else tlon

            # tiles latitude indexes go from 1 to 24, covering the range from 60 to
            # -60
            tlat = 1+np.floor((60-lat)/5)
            tlat = 24 if tlat == 25 else tlat
            srtm_list.append("srtm_%02d_%02d" % (tlon, tlat))

    srtm_list = set(srtm_list)
    print("Needed srtm tiles: ", srtm_list)
    return srtm_list


def main(cfg_in, srtm_dir):
    """
    Get SRTM tiles from the reference image given in a json file.

    Args:
        cfg_in: user config dictionary
        srtm_dir: path to the srtm directory
    """
    img = cfg_in['images'][0]['img']
    rpc = cfg_in['images'][0]['rpc']

    dataset = gdal.Open(img)
    x, y, w, h = 0, 0, dataset.RasterXSize, dataset.RasterYSize

    srtm_list = list_srtm_tiles(rpc, x, y, w, h)
    files_list = []

    for s in srtm_list:
        files_list.append(download_srtm_tile(s, srtm_dir))

    cfg_out = copy.deepcopy(cfg_in)
    vrt = os.path.join(srtm_dir, 'srtm.vrt')
    gdal.BuildVRT(vrt, files_list)
    cfg_out['exogenous_dem'] = vrt
    return cfg_out


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=('Get SRTM tiles tool'))
    
    parser.add_argument('config_in', metavar='config_in.json',
                        help=('path to a json file containing the paths to '
                              'input and output files and the algorithm '
                              'parameters'))

    parser.add_argument('srtm_dir', metavar='srtm_dir',
                        help=('path to the srtm directory'))

    parser.add_argument('config_out', metavar='config_out.json',
                        help=('path to a json output file (adding exogenous_dem location)'))

    args = parser.parse_args()

    cfg_in = s2p.read_config_file(args.config_in)
    cfg_out = main(cfg_in, args.srtm_dir)

    # store a json dump of the config.cfg dictionary
    with open(args.config_out, 'w') as f:
        json.dump(cfg_out, f, indent=2, default=workaround_json_int64)

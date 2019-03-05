#!/usr/bin/env python

# Copyright (C) 2017, Julien Michel <julien.michel@cnes.fr>
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published                                                                                                                                        # by the Free Software Foundation, either version 3 of the License, or
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
import sys
import os
import tempfile

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from s2plib import rpc_model
from s2plib import common


class ImageRPCModifier():
    def __init__(self):
        self.available_filters = ['near', 'bilinear',
                     'cubic', 'cubicspline', 'lanczos', 'average']


    def scale_and_crop(self, image_in, rpc_in, write_folder, scaling_filter, scale_x, crop, scale_first=False, scale_y=None):
        # Either scales, crops or does both an image and a RPC file

        if not crop and not scaling_filter:
            print("Neither crop nor scaling filter was provided, nothing was done")
            return(0)

        if not os.path.exists(write_folder):
            os.makedirs(write_folder)

        if crop and not scaling_filter:
            print("Only croping.")
            self.crop_image_and_rpc(image_in, rpc_in, write_folder, crop)
            return(0)

        if scaling_filter and not crop:
            print("Only scaling.")
            self.scale_image_and_rpc(image_in, rpc_in, write_folder, scaling_filter, scale_x, scale_y)
            return(0)

        else:
            if scale_first:
                print("Scaling then croping")
                self.scale_image_and_rpc(image_in, rpc_in, write_folder, scaling_filter, scale_x, scale_y)
                new_image_in = self._get_output_name(write_folder, image_in, False)
                new_rpc_in = self._get_output_name(write_folder, rpc_in, False)
                self.crop_image_and_rpc(new_image_in, new_rpc_in, write_folder, crop)
                os.remove(new_image_in)
                os.remove(new_rpc_in)
                return(0)
            else:
                print("Croping then scaling")
                self.crop_image_and_rpc(image_in, rpc_in, image_out, rpc_out, crop)
                new_image_in = self._get_output_name(write_folder, image_in)
                new_rpc_in = self._get_output_name(write_folder, rpc_in)
                self.scale_image_and_rpc(new_image_in, new_rpc_in, write_folder, scaling_filter, scale_x, scale_y)
                os.remove(new_image_in)
                os.remove(new_rpc_in)
                return(0)

    def scale_image_and_rpc(self, image_in, rpc_in, write_folder, scaling_filter, scale_x, scale_y=None):
        scale_x = float(scale_x)
        if scale_y is None:
            scale_y = scale_x
        else:
            scale_y = float(scale_y)
        assert scaling_filter in self.available_filters

        self.scale_image(image_in, write_folder, scaling_filter, scale_x, scale_y)
        rpc_model = self.scale_rpc(rpc_in, scale_x, scale_y)
        rpc_out = self._get_output_name(write_folder, rpc_in, False)
        self.write_rpc(rpc_model, rpc_out)
        return(0)


    def scale_image(self, image_in, write_folder, scaling_filter, scale_x, scale_y):

        image_out = self._get_output_name(write_folder, image_in, False)

        # generate image
        print("Generating {} ...".format(image_out))

        # First get input image size
        w, h, _ = common.image_size_gdal(image_in)

        # Generate a temporary vrt file to have the proper geotransform
        fd, tmp_vrt = tempfile.mkstemp(suffix='.vrt',
                                       dir=write_folder)

        os.close(fd)


        common.run('gdal_translate -of VRT -a_ullr 0 0 {} {} {} {}'.format(w, h, image_in, tmp_vrt))

        common.run('gdalwarp -co RPB=NO -co PROFILE=GeoTIFF -r {} -co "BIGTIFF=IF_NEEDED" -co "TILED=YES" -ovr NONE \
                     -to SRC_METHOD=NO_GEOTRANSFORM -to DST_METHOD=NO_GEOTRANSFORM -tr \
                    {} {} {} {}'.format(scaling_filter, scale_x, scale_y, tmp_vrt, image_out))

        try:
            # Remove aux files if any
            os.remove(image_out + ".aux.xml")
        except OSError:
            pass

        # Clean tmp vrt file
        os.remove(tmp_vrt)

        print("Done")

        return(0)


    def scale_rpc(self, rpc_in, scale_x, scale_y):


        r = rpc_model.RPCModel(rpc_in)
        r.linScale /= scale_y
        r.linOff /= scale_y
        r.colScale /= scale_x
        r.colOff /= scale_x

        print("Done")
        return(r)

    def write_rpc(self, rpc_model, rpc_out):
        rpc_model.write(rpc_out)


    def crop_image_and_rpc(self, image_in, rpc_in, write_folder, crop):
        # crop should be x_min, y_min, width, height

        assert len(crop) == 4

        self.crop_image(image_in, write_folder, crop)
        rpc_model = self.crop_rpc(rpc_in, crop)
        rpc_out = self._get_output_name(write_folder, rpc_in)
        self.write_rpc(rpc_model, rpc_out)
        return(0)

    def crop_image(self, image_in, write_folder, crop):
        # crop should be x_min, y_min, width, height

        x_min, y_min, width, height = crop

        image_out = self._get_output_name(write_folder, image_in)
        # generate image
        print("Croping to image {} ...".format(image_out))

        common.run('gdal_translate -srcwin {} {} {} {} {} {}'.format(x_min, y_min, width, height, image_in, image_out))

        print("Done")

        return(0)


    def crop_rpc(self, rpc_in, crop):
        # crop should be x_min, y_min, width, height

        x_min, y_min, width, height = crop

        # generate rpc file
        print("croping rpc {} ...".format(rpc_in))

        r = rpc_model.RPCModel(rpc_in)

        r.linOff += y_min
        r.colOff += x_min

        print("Done")
        return(r)


    def _get_output_name(self, write_folder, file_path, crop=True):
        _, file_name = os.path.split(file_path) 
        if crop:
            new_file_name = os.path.join(write_folder, file_name.replace(".", "_cropped."))
        else:
            new_file_name = os.path.join(write_folder, file_name.replace(".", "_scaled."))
        return new_file_name


    def _verify_rpc_model_scaled(self, rpc_1, rpc_2, scale_x, scale_y):
        #sanity check of model
        lon_gt, lat_gt, alt_gt = rpc_1.direct_estimate(scale_x*5000, scale_y*10000, 0)
        lon, lat, alt = rpc_2.direct_estimate(5000, 10000, 0)

        print(abs(lon_gt - lon), abs(lat_gt - lat), abs(alt_gt - alt))
        assert abs(lon_gt - lon) < 10**(-8)
        assert abs(lat_gt - lat) < 10**(-8)
        assert abs(alt_gt - alt) < 10**(-8)

        print("models match")

        return(0)

    def _verify_rpc_model_crop(self, rpc_1, rpc_2):
        #sanity check of model

        rpc_model_before_croping = rpc_model.RPCModel(rpc_1)
        rpc_model_after_after = rpc_model.RPCModel(rpc_2)

        offset_w, offset_h = [100, 100]

        lon_gt, lat_gt, alt_gt = rpc_model_before_scaling.direct_estimate(200, 300, 0)
        lon, lat, alt = rpc_model_after_croping.direct_estimate(200 + offset_w, 300 + offset_h, 0)

        assert abs(lon_gt - lon) < 10**(-8)
        assert abs(lat_gt - lat) < 10**(-8)
        assert abs(alt_gt - alt) < 10**(-8)

        return(0)

if __name__ == "__main__":




    parser = argparse.ArgumentParser(description='Rescale and/or crop an image and the RPC model.')
    parser.add_argument('image_in', type=str,
                        help='path to input image')
    parser.add_argument('rpc_in', type=str,
                        help='path to RPC model')
    parser.add_argument('write_folder', type=str,
                        help='path to write output')
    parser.add_argument('--scaling_filter', type=str, nargs='?',
                        help='filter for rescaling, can be one of [near, bilinear, cubic, cubicspline, lanczos, average]')
    parser.add_argument('--scale_x', type=float, nargs='?',
                        help='scaling factor along x axis')
    parser.add_argument('--scale_y', type=float, nargs='?',
                        help='scaling factor along y axis, if not provided scale_x will be used.')
    parser.add_argument('--crop', type=int, nargs='*',
                        help='crop information: x_min, y_min, width, height.')
    parser.add_argument('--scale_first', type=bool, default=True,
                        help='Defautl is True: in case of croping and scaling, wheter to scale first or not.')


    args = parser.parse_args()
    args = vars(args)

    modifier = ImageRPCModifier()
    modifier.scale_and_crop(**args)
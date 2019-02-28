#!/usr/bin/env python

import importlib
import os
import numpy as np
import s2plib.common as common
from abc import ABCMeta, abstractmethod


# link to ease circular import implied by class factory
this_file_realpath = os.path.dirname(os.path.realpath(__file__))
this_package_path = this_file_realpath[this_file_realpath.find('s2plib'):].replace(os.path.sep, '.')


# list of supported stereo matching algotirhms
stereo_methods_avail = {'mgm': {'module': 'mgm', 'class': 'mgmMatching'},
                        'msm3': {'module': 'msmw', 'class': 'msmwMatching'}}


# Metaclass Abstract Factory
class StereoMatching(object):
    __metaclass__ = ABCMeta

    def __new__(cls, stereo_method='Unknown'):
        if cls is StereoMatching:
            try:
                stereo_method_desc = stereo_methods_avail[stereo_method]
                stereo_module = importlib.import_module('.'.join([this_package_path, stereo_method_desc['module']]))
                stereo_class = getattr(stereo_module, stereo_method_desc['class'])
                return super(StereoMatching, cls).__new__(stereo_class)
            except KeyError:
                raise KeyError('No stereo matching algorithm named {} supported'.format(stereo_method))
        else:
            return super(StereoMatching, cls).__new__(cls, stereo_method)

    @abstractmethod
    def desc(self):
        print('StereoMatching description')

    @abstractmethod
    def compute_disparity_map(self, im_ref, im_sec, out_disp_path, out_mask_path,
                              disp_min=None, disp_max=None, **kwargs):
        """
        Stereo matching

        :param im_ref: str, path to input ref image
        :param im_sec: str, path to input sec image
        :param out_disp_path: str, path to output disparity map
        :param out_mask_path: str, path to output mask map
        :param disp_min: float, min disp range
        :param disp_max: float, max disp range
        :param kwargs: optional arguments
        :return:
        """

    @abstractmethod
    def rectify_secondary_tile_only(self):
        pass

    def _check_disp_range(self, im_ref, disp_min, disp_max):
        """
        Constrain disp range edges according to input ref image size

        :param im_ref
        :param disp_min:
        :param disp_max:
        :return: disp_min, disp_max
        """

        if not self.rectify_secondary_tile_only():
            disp_min = [disp_min]
            disp_max = [disp_max]

        # limit disparity bounds
        np.alltrue(len(disp_min) == len(disp_max))
        for dim in range(len(disp_min)):
            if disp_min[dim] is not None and disp_max[dim] is not None:
                image_size = common.image_size_gdal(im_ref)
                if disp_max[dim] - disp_min[dim] > image_size[dim]:
                    center = 0.5 * (disp_min[dim] + disp_max[dim])
                    disp_min[dim] = int(center - 0.5 * image_size[dim])
                    disp_max[dim] = int(center + 0.5 * image_size[dim])

            # round disparity bounds
            if disp_min[dim] is not None:
                disp_min[dim] = int(np.floor(disp_min[dim]))
            if disp_max is not None:
                disp_max[dim] = int(np.ceil(disp_max[dim]))

        if not self.rectify_secondary_tile_only():
            disp_min = disp_min[0]
            disp_max = disp_max[0]

        return disp_min, disp_max


if __name__ == '__main__':
    #nothing = StereoMatching()
    mgm = StereoMatching('mgm')
    msmw = StereoMatching('msmw')

    #nothing.name()
    mgm.desc()
    msmw.desc()

    #msmw.a_method()

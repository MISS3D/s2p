#!/usr/bin/env python

import numpy as np
import s2plib.common as common
from abc import ABCMeta, abstractmethod


# StereoMatching is an abstract metaclass
#   - abstract because it was created from ABCMeta metaclass (and not type itself), and ABCMeta creates abstract classes
#   - metaclass because we overwrite __new__, which is a constructor, and a constructor creates a class,
#     and what creates a class is a metaclass
class StereoMatching(object):
    __metaclass__ = ABCMeta

    stereo_methods_avail = {}

    def __new__(cls, stereo_method='Unknown'):
        if cls is StereoMatching:
            if type(stereo_method) is str:
                try:
                    return super(StereoMatching, cls).__new__(cls.stereo_methods_avail[stereo_method])
                except KeyError:
                    raise KeyError('No stereo matching algorithm named {} supported'.format(stereo_method))
            else:
                if type(stereo_method) is unicode:
                    try:
                        return super(StereoMatching, cls).__new__(cls.stereo_methods_avail[stereo_method.encode('utf-8')])
                    except KeyError:
                        raise KeyError('No stereo matching algorithm named {} supported'.format(stereo_method))
                else:
                    try:
                        return super(StereoMatching, cls).__new__(stereo_method)
                    except:
                        raise
        else:
            return super(StereoMatching, cls).__new__(cls)

    @classmethod
    def register_subclass(cls, short_name):
        def decorator(subclass):
            cls.stereo_methods_avail[short_name] = subclass
            return subclass

        return decorator

    @abstractmethod
    def desc(self):
        print('StereoMatching description')

    @abstractmethod
    def compute_disparity_map(self, im_ref, im_sec, out_disp_path, out_mask_path,
                              min_disp_range=None, max_disp_range=None, **kwargs):
        """
        Stereo matching

        :param im_ref: str, path to input ref image
        :param im_sec: str, path to input sec image
        :param out_disp_path: str, path to output disparity map
        :param out_mask_path: str, path to output mask map
        :param min_disp_range: float, min disp range
        :param max_disp_range: float, max disp range
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
        image_size = common.image_size_gdal(im_ref)
        for dim in range(len(disp_min)):
            if disp_min[dim] is not None and disp_max[dim] is not None:
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
    """
    >>> import stereo_matching
    >>> stereo_matching.StereoMatching.__subclasses__()
    []
    >>> import mgm
    >>> stereo_matching.StereoMatching.__subclasses__()
    [<class 'mgm.mgmMatching'>, <class 'mgm.mgm_multiMatching'>]
    >>> stereo_matching.StereoMatching(mgm.mgmMatching)
    <mgm.mgmMatching object at 0x2adca3c8feb8>
    >>> stereo_matching.StereoMatching('mgm')
    <mgm.mgmMatching object at 0x2adc99c78ef0>
    >>> stereo_matching.StereoMatching.stereo_methods_avail
    {'mgm': <class 'mgm.mgmMatching'>}
    """


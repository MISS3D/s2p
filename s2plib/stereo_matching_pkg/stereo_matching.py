#!/usr/bin/env python

import importlib
import os
from abc import ABCMeta, abstractmethod


# link to ease circular import implied by class factory
this_file_realpath = os.path.dirname(os.path.realpath(__file__))
this_package_path = this_file_realpath[this_file_realpath.find('s2plib'):].replace(os.path.sep, '.')


# list of supported stereo matching algotirhms
stereo_methods_avail = {'mgm': {'module': 'mgm', 'class': 'mgmMatching'}}


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

    def a_method(self):
        print('a_method')


if __name__ == '__main__':
    #nothing = StereoMatching()
    mgm = StereoMatching('mgm')
    msmw = StereoMatching('msmw')

    #nothing.name()
    mgm.desc()
    msmw.desc()

    msmw.a_method()

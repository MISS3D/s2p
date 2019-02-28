#!/usr/bin/env python


from s2plib.stereo_matching_pkg.mgm import mgmMatching
from s2plib.stereo_matching_pkg.msmw import msmwMatching

from abc import ABCMeta, abstractmethod


stereo_methods_avail = {'mgm': mgmMatching,
                        'msmw': msmwMatching}


class StereoMatching(object):
    __metaclass__ = ABCMeta

    def __new__(cls, stereo_method='Unknown'):
        if cls is StereoMatching:
            try:
                print('StereoMatching __new__ if {}'.format(stereo_method))
                return super(StereoMatching, cls).__new__(stereo_methods_avail[stereo_method])
            except KeyError:
                raise KeyError('No stereo matching algorithm named {} supported'.format(stereo_method))
        else:
            print('StereoMatching __new__ else')
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

import os
import s2plib.common as common
try:
    import stereo_matching
except ImportError:
    from . import stereo_matching


@stereo_matching.StereoMatching.register_subclass('mgm')
class mgmMatching(stereo_matching.StereoMatching):

    _omp_num_threads = 1
    _median = 1
    _census_ncc_win = 5
    _tsgm = 3

    def desc(self):
        print('mgmMatching algorithm')

    def rectify_secondary_tile_only(self):
        return False

    def compute_disparity_map(self, im_ref, im_sec, out_disp_path, out_mask_path,
                              min_disp_range=None, max_disp_range=None, **kwargs):
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

        self._check_disp_range(im_ref, min_disp_range, max_disp_range)

        # define environment variables
        if 'omp_num_threads' in kwargs:
            self._omp_num_threads = str(kwargs['omp_num_threads'])
        if 'census_ncc_win' in kwargs:
            self._census_ncc_win = str(kwargs['census_ncc_win'])
        varenv = dict()
        varenv['OMP_NUM_THREADS'] = self._omp_num_threads
        varenv['MEDIAN'] = self._median
        varenv['CENSUS_NCC_WIN'] = self._census_ncc_win
        varenv['TSGM'] = self._tsgm
        env = os.environ.copy()
        self._export_environment(env, varenv)

        conf = '{}_confidence.tif'.format(os.path.splitext(out_disp_path)[0])
        common.run('{0} -r {1} -R {2} -s vfit -t census -O 8 {3} {4} {5} -confidence_consensusL {6}'.format('mgm',
                                                                                                            min_disp_range,
                                                                                                            max_disp_range,
                                                                                                            im_ref,
                                                                                                            im_sec,
                                                                                                            out_disp_path,
                                                                                                            conf),
                   env)

        # produce the mask: rejected pixels are marked with nan of inf in disp
        # map
        common.run('plambda {0} "isfinite" -o {1}'.format(out_disp_path, out_mask_path))


@stereo_matching.StereoMatching.register_subclass('mgm_multi')
class mgm_multiMatching(stereo_matching.StereoMatching):

    _omp_num_threads = 1
    _stereo_speckle_filter = 25
    _mindiff = 1
    _census_ncc_win = 5
    _subpix = 2

    def desc(self):
        print('mgm_multiMatching algorithm')

    def rectify_secondary_tile_only(self):
        return False

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

        self._check_disp_range(im_ref, min_disp_range, max_disp_range)

        # define environment variables
        if 'omp_num_threads' in kwargs:
            self._omp_num_threads = str(kwargs['omp_num_threads'])
        if 'stereo_speckle_filter' in kwargs:
            self._stereo_speckle_filter = str(kwargs['stereo_speckle_filter'])
        if 'census_ncc_win' in kwargs:
            self._census_ncc_win = str(kwargs['census_ncc_win'])

        varenv = dict()
        varenv['OMP_NUM_THREADS'] = self._omp_num_threads
        varenv['REMOVESMALLCC'] = self._stereo_speckle_filter
        varenv['MINDIFF']= self._mindiff
        varenv['CENSUS_NCC_WIN'] = self._census_ncc_win
        varenv['SUBPIX'] = self._subpix
        env = os.environ.copy()
        self._export_environment(env, varenv)

        # it is required that p2 > p1. The larger p1, p2, the smoother the disparity
        # TODO
        #   next three lines are commented since P1 & P2 do not seem to be used
        #regularity_multiplier = kwargs['stereo_regularity_multiplier']
        #P1 = 8*regularity_multiplier   # penalizes disparity changes of 1 between neighbor pixels
        #P2 = 32*regularity_multiplier  # penalizes disparity changes of more than 1
        conf = '{}_confidence.tif'.format(os.path.splitext(out_disp_path)[0])
        common.run('{0} -r {1} -R {2} -S 6 -s vfit '
                   '-t census {3} {4} {5} -confidence_consensusL {6}'.format('mgm_multi',
                                                                             min_disp_range,
                                                                             max_disp_range,
                                                                             im_ref, im_sec,
                                                                             out_disp_path, conf),
                   env)

        # produce the mask: rejected pixels are marked with nan of inf in disp
        # map
        common.run('plambda {0} "isfinite" -o {1}'.format(out_disp_path, out_mask_path))

import os
import s2plib.common as common
from . import stereo_matching


class mgmMatching(stereo_matching.StereoMatching):

    def desc(self):
        print('mgmMatching class')

    def rectify_secondary_tile_only(self):
        return False

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

        # define environment variables
        env = os.environ.copy()
        env['OMP_NUM_THREADS'] = str(kwargs['omp_num_threads'])
        env['MEDIAN'] = '1'
        env['CENSUS_NCC_WIN'] = str(kwargs['census_ncc_win'])
        env['TSGM'] = '3'

        conf = '{}_confidence.tif'.format(os.path.splitext(out_disp_path)[0])
        common.run('{0} -r {1} -R {2} -s vfit -t census -O 8 {3} {4} {5} -confidence_consensusL {6}'.format('mgm',
                                                                                                            disp_min,
                                                                                                            disp_max,
                                                                                                            im_ref,
                                                                                                            im_sec,
                                                                                                            out_disp_path,
                                                                                                            conf),
                   env)

        # produce the mask: rejected pixels are marked with nan of inf in disp
        # map
        common.run('plambda {0} "isfinite" -o {1}'.format(out_disp_path, out_mask_path))

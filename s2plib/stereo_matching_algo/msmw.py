import os
import s2plib.common as common
try:
    # case when the package stereo matching lives on his own
    import stereo_matching
except ImportError:
    # case when the package stereo matching is used from another package (say s2plib)
    from . import stereo_matching


@stereo_matching.StereoMatching.register_subclass('msmw3')
class msmwMatching(stereo_matching.StereoMatching):

    def desc(self):
        print('msmwMatching algorithm')

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
        env = os.environ.copy()
        env['OMP_NUM_THREADS'] = str(kwargs['omp_num_threads'])

        bm_binary = 'msmw'
        common.run('{0} -m {1} -M {2} -il {3} -ir {4} -dl {5} -kl {6}'.format(
                bm_binary, min_disp_range, max_disp_range, im_ref, im_sec, out_disp_path, out_mask_path))


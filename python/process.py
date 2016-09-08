# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>
# Copyright (C) 2015, Julien Michel <julien.michel@cnes.fr>


import numpy as np
from config import cfg
import os
import sys
import multiprocessing

from python import common
from python import geographiclib
from python import triangulation
from python import fusion
from python import block_matching
from python import rectification
from python import masking


def color_crop_ref(tile_info, clr=None):
    """
    Colorizations of a crop_ref (for a given tile)

    Args:
        tile_info: a dictionary that provides all you need to process a tile
        clr (optional): if crop_ref is a pan image, will perform the pansharpening with the color image clr

        If clr is None then:
            case 1: tile is an RGBI image, so removes I channel, and perform rescaling of the remaining channels
            case 2: tile is already an RGB image, so just perform rescaling

        Note that if rescaling is already performed, then the file applied_minmax.txt exists:
            if applied_minmax.txt exists and cfg['skip_existing'] is True, rescaling won't be performed again
            if applied_minmax.txt exists and is different from global_minmax, rescaling will be compulsorily performed (can occur if a new tile is added)
    """
    # get info
    x, y, w, h = tile_info['coordinates']
    tile_dir = tile_info['directory']
    z = cfg['subsampling_factor']

    # paths
    crop_ref = os.path.join(tile_dir , 'roi_ref_crop.tif')
    global_minmax = os.path.join(cfg['out_dir'] , 'global_minmax.txt')
    applied_minmax = os.path.join(tile_dir , 'applied_minmax.txt')
    
    global_minmax_arr = np.loadtxt(global_minmax)

    if cfg['color_ply']:

        doProcess = False
        if not os.path.exists(applied_minmax):
            doProcess = True
            applied_minmax_arr = global_minmax_arr
        else:
            applied_minmax_arr = np.loadtxt(applied_minmax)

            if (applied_minmax_arr[0] != global_minmax_arr[0]) or (applied_minmax_arr[1] != global_minmax_arr[1]):
                doProcess = True
                applied_minmax_arr = global_minmax_arr

        if not doProcess and cfg['skip_existing']:
            print 'Rescaling of tile %s already done, skip' % tile_dir
        else:

            crop_color = os.path.join(tile_dir , 'roi_color_ref.tif')
            if clr is not None:
                triangulation.colorize(crop_ref, clr, x, y, z, crop_color,
                                       applied_minmax_arr[0],
                                       applied_minmax_arr[1])
            else:  # use of image_rescaleintensities
                np.savetxt(applied_minmax, applied_minmax_arr)

                if common.image_pix_dim_tiffinfo(crop_ref) == 4:
                    print 'the image is pansharpened fusioned'
                    tmp = common.rgbi_to_rgb(crop_ref, out=None, tilewise=True)
                    #common.image_qauto(tmp, crop_color, tilewise=False)
                    common.image_rescaleintensities(tmp, crop_color,
                                                    applied_minmax_arr[0],
                                                    applied_minmax_arr[1])
                else:
                    print 'no color data'
                    #common.image_qauto(crop_ref, crop_color, tilewise=False)
                    common.image_rescaleintensities(crop_ref, crop_color,
                                                    applied_minmax_arr[0],
                                                    applied_minmax_arr[1])


def generate_cloud(tile_info, do_offset=False, utm_zone=None):
    """
    Args:
        tile_info: a dictionary that provides all you need to process a tile
        do_offset (optional, default: False): boolean flag to decide wether the
            x, y coordinates of points in the ply file will be translated or not
            (translated to be close to 0, to avoid precision loss due to huge
            numbers)
    """
    print "\nComputing point cloud..."

    # get info
    tile_dir = tile_info['directory']
    x, y, w, h = tile_info['coordinates']
    img1, rpc1 = cfg['images'][0]['img'], cfg['images'][0]['rpc']

    height_map = os.path.join(tile_dir , 'height_map_crop.tif')
    crop_color = os.path.join(tile_dir , 'roi_color_ref.tif')
    if not os.path.exists(crop_color):
        crop_color = ''

    z = cfg['subsampling_factor']

    # Compute the homography transforming the coordinates system of the
    # original full size image into the coordinates system
    # of the crop we are dealing with
    # col and row have been divided by z inside 'finalize_tile' for
    # convinience; col*z and row*z allow to get the initial values back.
    A = common.matrix_translation(-x * z, -y * z)
    f = 1.0 / z
    Z = np.diag([f, f, 1])
    A = np.dot(Z, A)
    trans = os.path.join(tile_dir , 'trans.txt')
    np.savetxt(trans, A, fmt='%9.3f')

    # compute coordinates (offsets) of the point we want to use as origin in
    # the local coordinate system of the computed cloud
    if do_offset:
        r = rpc_model.RPCModel(rpc1)
        lat = r.latOff
        lon = r.lonOff
        off_x, off_y = geographiclib.geodetic_to_utm(lat, lon)[0:2]
    else:
        off_x, off_y = 0, 0

    # output
    cloud = os.path.join(tile_dir , 'cloud.ply')

    triangulation.compute_point_cloud(cloud, height_map, rpc1, trans, crop_color,
                                      off_x, off_y, utm_zone=utm_zone)

    common.garbage_cleanup()


def finalize_tile(tile_info, utm_zone=None):
    """
    Finalize the processing of a tile.

    remove overlapping areas, get the
    colors from a XS image and use it to color and generate a ply file
    (colorization is not mandatory)

    Args:
        tile_info: a dictionary that provides all you need to process a tile
    """
    ## get info
    tile_dir = tile_info['directory']
    x, y, w, h = tile_info['coordinates']
    ov = tile_info['overlap']
    pos = tile_info['position_type']
    img1, rpc1 = cfg['images'][0]['img'], cfg['images'][0]['rpc']
    nb_pairs = tile_info['number_of_pairs']

    ## remove overlapping areas
    height_map = os.path.join(tile_dir , 'height_map.tif')
    height_map_crop = os.path.join(tile_dir , 'height_map_crop.tif')
    rpc_err_all = os.path.join(tile_dir , 'rpc_err_all.tif')
    rpc_err_all_crop = os.path.join(tile_dir , 'rpc_err_all_crop.tif')
    crop_ref = os.path.join(tile_dir , 'roi_ref.tif')
    crop_ref_crop = os.path.join(tile_dir , 'roi_ref_crop.tif')
    
    if cfg['full_vrt']:
        nb_views = os.path.join(tile_dir , 'nb_views.tif')
        nb_views_crop = os.path.join(tile_dir , 'nb_views_crop.tif')

    dicoPos = {}
    dicoPos['M'] = [ov / 2, ov / 2, -ov, -ov]
    dicoPos['L'] = [0, ov / 2, -ov / 2, -ov]
    dicoPos['R'] = [ov / 2, ov / 2, -ov / 2, -ov]
    dicoPos['U'] = [ov / 2, 0, -ov, -ov / 2]
    dicoPos['B'] = [ov / 2, ov / 2, -ov, -ov / 2]
    dicoPos['UL'] = [0, 0, -ov / 2, -ov / 2]
    dicoPos['UR'] = [ov / 2, 0, -ov / 2, -ov / 2]
    dicoPos['BR'] = [ov / 2, ov / 2, -ov / 2, -ov / 2]
    dicoPos['BL'] = [0, ov / 2, -ov / 2, -ov / 2]
    dicoPos['Single'] = [0, 0, 0, 0]

    z = cfg['subsampling_factor']
    newcol, newrow, difftw, diffth = np.array(dicoPos[pos]) / z
    x = x / z + newcol
    y = y / z + newrow
    w = w / z + difftw
    h = h / z + diffth
    tile_info['coordinates'] = (x, y, w, h)
    
    # z=1 because height_map, crop_ref (and so forth) have
    # already been zoomed. So don't zoom again to crop these images.
    if not (os.path.isfile(height_map_crop) and cfg['skip_existing']):
        common.cropImage(height_map, height_map_crop,
                         newcol, newrow, w, h)
    if not (os.path.isfile(rpc_err_all_crop) and cfg['skip_existing']):
        common.cropImage(rpc_err_all, rpc_err_all_crop,
                         newcol, newrow, w, h)
    
    if cfg['full_vrt']:
        
        if not (os.path.isfile(nb_views_crop) and cfg['skip_existing']):
            common.cropImage(nb_views, nb_views_crop, newcol, newrow, w, h)
            
        for img_id in xrange(1,len(cfg['images'])+1): 
            #selected sights
            selected_sighti = os.path.join(tile_dir ,
                            'selected_sight%d.tif' % img_id)
            selected_sighti_crop = os.path.join(tile_dir ,
                                'selected_sight%d_crop.tif' % img_id) 
            if not (os.path.isfile(selected_sighti_crop) and cfg['skip_existing']):
                common.cropImage(selected_sighti, selected_sighti_crop,
                             newcol, newrow, w, h)
            
            # err by sight                   
            rpc_err_sighti = os.path.join(tile_dir , 
                                'rpc_err_sight%d.tif' % img_id)
            rpc_err_sighti_crop = os.path.join(tile_dir , 
                                'rpc_err_sight%d_crop.tif' % img_id)
            if not (os.path.isfile(rpc_err_sighti_crop) and cfg['skip_existing']):
                common.cropImage(rpc_err_sighti, rpc_err_sighti_crop,
                             newcol, newrow, w, h) 
                             
            # err vector by sight (from opt point to a given sight)                                   
            rpc_err_veci = os.path.join(tile_dir , 
                                'rpc_err_vec%d.tif' % img_id)
            rpc_err_veci_crop = os.path.join(tile_dir , 
                                'rpc_err_vec%d_crop.tif' % img_id)
            if not (os.path.isfile(rpc_err_veci_crop) and cfg['skip_existing']):
                common.cropImage(rpc_err_veci, rpc_err_veci_crop,
                             newcol, newrow, w, h) 
                                       
            # reprojected err vector by sight (from opt point to a given sight)                                   
            rpc_err_vec_rpji = os.path.join(tile_dir , 
                                'rpc_err_vec_rpj%d.tif' % img_id)
            rpc_err_vec_rpji_crop = os.path.join(tile_dir , 
                                'rpc_err_vec_rpj%d_crop.tif' % img_id)
            if not (os.path.isfile(rpc_err_vec_rpji_crop) and cfg['skip_existing']):
                common.cropImage(rpc_err_vec_rpji, rpc_err_vec_rpji_crop,
                             newcol, newrow, w, h)           
                             
        for pair_id in xrange(1,nb_pairs+1):
            # 2D disparities (if originaly computed in epipolar geometry)
            disp2Di = os.path.join(tile_dir ,
                            'pair_%d/disp2D.tif' % pair_id)
            disp2Di_crop = os.path.join(tile_dir ,
                            'pair_%d/disp2D_crop.tif' % pair_id)
            if not (os.path.isfile(disp2Di_crop) and cfg['skip_existing']):
                common.cropImage(disp2Di, disp2Di_crop,
                             newcol, newrow, w, h)
        
    
    # ref                     
    if not (os.path.isfile(crop_ref_crop) and cfg['skip_existing']):
        common.cropImage(crop_ref, crop_ref_crop, newcol, newrow, w, h)

    # colors
    color_crop_ref(tile_info, cfg['images'][0]['clr'])
    
    ## Generate cloud
    generate_cloud(tile_info, cfg['offset_ply'], utm_zone)


def rectify(tile_dir, A_global, img1, rpc1, img2, rpc2, x=None, y=None,
            w=None, h=None, prv1=None):
    """
    Computes rectifications, without tiling

    Args:
        tile_dir: path to the output directory
        A_global: global pointing corrections
        img1: path to the reference image.
        rpc1: paths to the xml file containing the rpc coefficients of the
            reference image
        img2: path to the secondary image.
        rpc2: paths to the xml file containing the rpc coefficients of the
            secondary image
        x, y, w, h: four integers defining the rectangular ROI in the reference
            image. (x, y) is the top-left corner, and (w, h) are the dimensions
            of the rectangle.
        prv1 (optional): path to a preview of the reference image.
    """
    # output files
    rect1 = os.path.join(tile_dir,'rectified_ref.tif')
    rect2 = os.path.join(tile_dir,'rectified_sec.tif')
    disp = os.path.join(tile_dir,'rectified_disp.tif')
    mask = os.path.join(tile_dir,'rectified_mask.png')
    subsampling = os.path.join(tile_dir,'subsampling.txt')
    pointing = os.path.join(tile_dir,'pointing.txt')
    center = os.path.join(tile_dir,'center_keypts_sec.txt')
    sift_matches = os.path.join(tile_dir,'sift_matches.txt')
    sift_matches_plot = os.path.join(tile_dir,'sift_matches_plot.png')
    H_ref = os.path.join(tile_dir,'H_ref.txt')
    H_sec = os.path.join(tile_dir,'H_sec.txt')
    disp_min_max = os.path.join(tile_dir,'disp_min_max.txt')
    config = os.path.join(tile_dir,'config.json')

    A, m = None, None

    if os.path.isfile(os.path.join(tile_dir,'pointing.txt')):
        A = np.loadtxt(os.path.join(tile_dir,'pointing.txt'))
    else:
        A = A_global
    if os.path.isfile(os.path.join(tile_dir,'sift_matches.txt')):
        m = np.loadtxt(os.path.join(tile_dir,'sift_matches.txt'))

    # rectification
    H1, H2, disp_min, disp_max = rectification.rectify_pair(img1, img2, rpc1,
                                                            rpc2, x, y, w, h,
                                                            rect1, rect2, A, m)

    z = cfg['subsampling_factor']
    np.savetxt(subsampling, np.array([z]), fmt='%.1f')
    np.savetxt(H_ref, H1, fmt='%12.6f')
    np.savetxt(H_sec, H2, fmt='%12.6f')
    np.savetxt(disp_min_max, np.array([disp_min, disp_max]), fmt='%3.1f')


def disparity(tile_dir, img1, rpc1, img2, rpc2, x=None, y=None,
              w=None, h=None, prv1=None):
    """
    Computes a disparity map from a Pair of Pleiades images, without tiling

    Args:
        tile_dir: path to the output directory
        img1: path to the reference image.
        rpc1: paths to the xml file containing the rpc coefficients of the
            reference image
        img2: path to the secondary image.
        rpc2: paths to the xml file containing the rpc coefficients of the
            secondary image
        x, y, w, h: four integers defining the rectangular ROI in the reference
            image. (x, y) is the top-left corner, and (w, h) are the dimensions
            of the rectangle.
        prv1 (optional): path to a preview of the reference image
        cld_msk (optional): path to a gml file containing a cloud mask
        roi_msk (optional): path to a gml file containing a mask defining the
            area contained in the full image
        wat_msk (optional): path to a tiff file containing a water mask.
    """
    # output files
    rect1 = os.path.join(tile_dir,'rectified_ref.tif')
    rect2 = os.path.join(tile_dir,'rectified_sec.tif')
    disp = os.path.join(tile_dir,'rectified_disp.tif')
    mask = os.path.join(tile_dir,'rectified_mask.png')
    cwid_msk = os.path.join(tile_dir,'cloud_water_image_domain_mask.png')
    cwid_msk_rect = os.path.join(tile_dir,'rectified_cloud_water_image_domain_mask.png')
    subsampling = os.path.join(tile_dir,'subsampling.txt')
    pointing = os.path.join(tile_dir,'pointing.txt')
    center = os.path.join(tile_dir,'center_keypts_sec.txt')
    sift_matches = os.path.join(tile_dir,'sift_matches.txt')
    sift_matches_plot = os.path.join(tile_dir,'sift_matches_plot.png')
    H_ref = os.path.join(tile_dir,'H_ref.txt')
    H_sec = os.path.join(tile_dir,'H_sec.txt')
    disp_min_max = os.path.join(tile_dir,'disp_min_max.txt')
    config = os.path.join(tile_dir,'config.json')

    # disparity (block-matching)
    disp_min, disp_max = np.loadtxt(disp_min_max)

    if cfg['disp_min'] is not None:
        disp_min = cfg['disp_min']
    if cfg['disp_max'] is not None:
        disp_max = cfg['disp_max']
    block_matching.compute_disparity_map(rect1, rect2, disp, mask,
                                         cfg['matching_algorithm'], disp_min,
                                         disp_max)

    # intersect mask with the cloud_water_image_domain mask (recomputed here to
    # get to be sampled on the epipolar grid)
    ww, hh = common.image_size(rect1)
    H1 = np.loadtxt(H_ref)
    H_inv = np.array([[1, 0, x], [0, 1, y], [0, 0, 1]])
    common.image_apply_homography(cwid_msk_rect, cwid_msk, np.dot(H1,H_inv), ww, hh)

    try:
        masking.intersection(mask, mask, cwid_msk_rect)
        masking.erosion(mask, mask, cfg['msk_erosion'])
    except OSError:
        print "file %s not produced" % mask


def triangulate(tile_info, prv1=None, A=None):
    """
    Computes triangulations, without tiling

    Args:
        tile_info : a dictionary that provides all you need to process a tile
        prv1 (optional): path to a preview of the reference image
        #A (optional, default None): pointing correction matrix.

    """
    # get info
    col, row, tw, th = tile_info['coordinates']
    z = cfg['subsampling_factor']
    global_out_dir = cfg['out_dir']
       
    rpc_list=[]
    images = cfg['images']
    for i in xrange(len(images)):
        rpc_list.append(images[i]['rpc'])

    # triangulation
    triangulation.compute_height_map(global_out_dir, 
                              col, row, tw, th, z,
                              rpc_list)

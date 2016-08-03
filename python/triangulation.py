# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>
# Copyright (C) 2015, Julien Michel <julien.michel@cnes.fr>

import os
import sys
import numpy as np
import multiprocessing

import common
from config import cfg


def compute_height_map(global_out_dir, 
                       col, row, tw, th, z, rpc_list):
    """
    Computes a height map from a disparity map, using a list of rpc.

    Args:
        global_out_dir : global s2p output directory
        col, row, tw, th: four integers defining the rectangular ROI in the original
            image. (col, row) is the top-left corner, and (tw, th) are the dimensions
            of the tile.
        z: zoom factor (usually 1, 2 or 4) used to produce the input disparity
            map
        rpc_list : list of rpc for each involved image
    """
    #TODO
    #if A is not None:
        #HH2 = common.tmpfile('.txt')
        #np.savetxt(HH2, np.dot(np.loadtxt(H2), np.linalg.inv(A)))
    #else:
        #HH2 = H2
        
    #bool_A = (A is None)
    
    rpc_list_str=''
    for rpc in rpc_list:
        rpc_list_str += rpc + ' '
        
    cmd = "disp_to_heights %s %s %s %s %s %s %s " % (global_out_dir, 
                       col, row, tw, th, z, 
                       rpc_list_str)

    #print cmd
    common.run(cmd)

    return


def transfer_map(in_map, H, x, y, w, h, zoom, out_map):
    """
    Transfer the heights computed on the rectified grid to the original
    Pleiades image grid.

    Args:
        in_map: path to the input map, usually a height map or a mask, sampled
            on the rectified grid
        H: path to txt file containing a numpy 3x3 array representing the
            rectifying homography
        x, y, w, h: four integers defining the rectangular ROI in the original
            image. (x, y) is the top-left corner, and (w, h) are the dimensions
            of the rectangle.
        zoom: zoom factor (usually 1, 2 or 4) used to produce the input height
            map
        out_map: path to the output map
    """
    # write the inverse of the resampling transform matrix. In brief it is:
    # homography * translation * zoom
    # This matrix transports the coordinates of the original cropped and
    # zoomed grid (the one desired for out_height) to the rectified cropped and
    # zoomed grid (the one we have for height)
    Z = np.diag([zoom, zoom, 1])
    A = common.matrix_translation(x, y)
    HH = np.dot(np.loadtxt(H), np.dot(A, Z))

    # apply the homography
    # write the 9 coefficients of the homography to a string, then call synflow
    # to produce the flow, then backflow to apply it
    # zero:256x256 is the iio way to create a 256x256 image filled with zeros
    hij = ' '.join(['%r' % num for num in HH.flatten()])
    common.run('synflow hom "%s" zero:%dx%d /dev/null - | BILINEAR=1 backflow - %s %s' % (
        hij, w/zoom, h/zoom, in_map, out_map))

    # w and h, provided by the s2p.process_pair_single_tile function, are
    # always multiples of zoom.

    # replace the -inf with nan in the heights map, because colormesh filter
    # out nans but not infs
    # implements: if isinf(x) then nan, else x
    # common.run('plambda %s "x isinf nan x if" > %s' % (tmp_h, out_height))


def compute_dem(global_out_dir, 
                col, row, tw, th, z, rpc_list):
    """
    Computes an altitude map, on the grid of the original reference image, from
    a disparity map given on the grid of the rectified reference image.

    Args:
        global_out_dir : global s2p output directory
        col, row, tw, th: four integers defining the rectangular ROI in the original
            image. (col, row) is the top-left corner, and (tw, th) are the dimensions
            of the tile.
        z: zoom factor (usually 1, 2 or 4) used to produce the input disparity
            map
        rpc_list : list of rpc for each involved image
    """

    compute_height_map(global_out_dir, 
                       col, row, tw, th, z, rpc_list)


def colorize(crop_panchro, im_color, x, y, zoom, out_colorized, rmin,rmax):
    """
    Colorizes a Pleiades gray crop using low-resolution color information.

    Args:
        crop_panchro: path to the panchro (ie gray) crop
        im_color: path to the full color image (tiff or jp2)
        x, y: coordinates of the top-left corner of crop_panchro, in the full
            Pleiade image frame.
        zoom: subsampling zoom-factor that was used to generate crop_panchro
        out_colorized: path to the output file
    """
    # get a translated and zoomed crop from the color image. It has to be
    # sampled on exactly the same grid as the panchro crop.
    # To do that we compose the translation + zoom transformation with a 4x
    # zoom (because color pleiades images have 4x lower resolution).  There is
    # also a small horizontal translation (4 pixels at the panchro resolution)
    # The resulting transformation is the composition of:
    #   translation (-1 - x/4, -y/4)
    #   zoom 4/z
    w, h = common.image_size_tiffinfo(crop_panchro)
    xx = np.floor(x / 4.0) 
    yy = np.floor(y / 4.0)
    ww = np.ceil((x + w * zoom) / 4.0) - xx 
    hh = np.ceil((y + h * zoom) / 4.0) - yy
    crop_ms = common.image_crop_tif(im_color, xx, yy, ww, hh)
    crop_ms = common.image_zoom_gdal(crop_ms, zoom/4.0)
    # crop_ms = common.image_safe_zoom_fft(crop_ms, zoom/4.0)

    # crop the crop_ms image to remove the extra-pixels due to the integer crop
    # followed by zoom
    x0 = max(0,x - 4*xx)
    y0 = max(0,y - 4*yy)
    crop_ms = common.image_crop_tif(crop_ms, x0, y0, w, h)
    assert(common.image_size_tiffinfo(crop_panchro) ==
           common.image_size_tiffinfo(crop_ms))

    # convert rgbi to rgb
    rgb = common.rgbi_to_rgb(crop_ms, out=None, tilewise=True)

    # blend intensity and color to obtain the result
    # each channel value r, g or b is multiplied by 3*y / (r+g+b), where y
    # denotes the panchro intensity
    tmp = common.tmpfile('.tif')
    pcmd = "dup split + + / * 3 *"
    os.environ['TMPDIR'] = os.path.join(cfg['temporary_dir'], 'meta/')
    cmd = 'tiffu meta \"plambda ^ ^1 \\\"%s\\\" -o @\" %s %s -- %s' % (pcmd,
                                                                      crop_panchro,
                                                                      rgb, tmp)
    common.run(cmd)
    if w * h > 25e6:  # image larger than 5000 x 5000 pixels
        common.image_qauto_otb(out_colorized, tmp)
    else:
        #common.image_qauto(tmp, out_colorized)
        common.image_rescaleintensities(tmp, out_colorized, rmin, rmax)
    return


def compute_point_cloud(cloud, heights, rpc, H=None, crop_colorized='',
                        off_x=None, off_y=None, ascii_ply=False,
                        with_normals=False, utm_zone=None):
    """
    Computes a color point cloud from a height map.

    Args:
        cloud: path to the output points cloud (ply format)
        heights: height map, sampled on the same grid as the crop_colorized
            image. In particular, its size is the same as crop_colorized.
        rpc: path to xml file containing RPC data for the current Pleiade image
        H (optional, default None): path to the file containing the coefficients
            of the homography transforming the coordinates system of the
            original full size image into the coordinates system of the crop we
            are dealing with.
        crop_colorized (optional, default ''): path to a colorized crop of a
            Pleiades image
        off_{x,y} (optional, default None): coordinates of the point we want to
            use as origin in the local coordinate system of the computed cloud
        ascii_ply (optional, default false): boolean flag to tell if the output
            ply file should be encoded in plain text (ascii).
        utm_zone (optional, default None):
    """
    hij = " ".join([str(x) for x in np.loadtxt(H).flatten()]) if H else ""
    asc = "--ascii" if ascii_ply else ""
    nrm = "--with-normals" if with_normals else ""
    utm = "--utm-zone %s" % utm_zone if utm_zone else ""
    if not (os.path.isfile(cloud) and cfg['skip_existing']):
        command = "colormesh %s %s %s %s -h \"%s\" %s %s %s" % (cloud, heights, rpc,
                                                                crop_colorized, hij,
                                                                asc, nrm, utm)
        if off_x:
            command += " --offset_x %d" % off_x
        if off_y:
            command += " --offset_y %d" % off_y
        common.run(command)

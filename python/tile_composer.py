# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>
# Copyright (C) 2015, Julien Michel <julien.michel@cnes.fr>

import gc
import sys
import os.path
import numpy as np

import piio
import common
from config import cfg


def mosaic_gdal(fout, w, h, list_tiles, tw, th, ov):
    """
    Compose several tiles of the same size into a bigger image (using gdal vrt)

    Args:
        fout: path to the output image
        w, h: output image dimensions
        list_tiles: list containing paths to the input tiles
        tw, th: dimensions of a tile (they must all have the same dimensions)
        ov: overlap between tiles (in pixels)

    Returns:
        nothing
    """
    N = len(list_tiles)
    ntx = np.ceil(float(w - ov) / (tw - ov)).astype(int)
    nty = np.ceil(float(h - ov) / (th - ov)).astype(int)
    assert(ntx * nty == N)

    vrtfilename = fout+'.vrt'

    vrtfile = open(vrtfilename, 'w')

    vrtfile.write("<VRTDataset rasterXSize=\"%i\" rasterYSize=\"%i\">\n" % (w,
                                                                            h))
    vrtfile.write("\t<VRTRasterBand dataType=\"Float32\" band=\"1\">\n")
    vrtfile.write("\t\t<ColorInterp>Gray</ColorInterp>\n")

    # loop over all the tiles
    for j in range(nty):
        for i in range(ntx):
            x0 = i * (tw - ov)
            y0 = j * (th - ov)
            x1 = min(x0 + tw, w)
            y1 = min(y0 + th, h)
            f = list_tiles[j * ntx + i]
            if os.path.isfile(f):
                # remove first dir name from path
                tile_fname = os.path.join(os.path.split(os.path.dirname(f))[1],
                                                        os.path.basename(f))
                vrtfile.write("\t\t<SimpleSource>\n")
                vrtfile.write("\t\t\t<SourceFilename relativeToVRT=\"1\">%s</SourceFilename>\n" % tile_fname)
                vrtfile.write("\t\t\t<SourceBand>1</SourceBand>\n")
                vrtfile.write("\t\t\t<SrcRect xOff=\"%i\" yOff=\"%i\" xSize=\"%i\" ySize=\"%i\"/>\n" % (0, 0, x1-x0, y1-y0))
                vrtfile.write("\t\t\t<DstRect xOff=\"%i\" yOff=\"%i\" xSize=\"%i\" ySize=\"%i\"/>\n" % (x0, y0, x1-x0, y1-y0))
                vrtfile.write("\t\t</SimpleSource>\n")

    vrtfile.write("\t</VRTRasterBand>\n")
    vrtfile.write("</VRTDataset>\n")
    vrtfile.close()

    common.run('gdal_translate %s %s' % (vrtfilename, fout))

    return


def mosaic_gdal2(fout, tiles_full_info, filename, w,h,z=1):
    """
    Compose several tiles of differents sizes into a bigger image (using gdal vrt)

    Args:
        fout: path to the output image
        fullInfo : all that you need to process a tile:
            col,row,tw,th,ov,i,j,pos,images=tiles_full_info[tile_dir]

    Returns:
        nothing
    """
    
    vrt_row = {}

    tw = 0
    th = 0

    files_to_remove = []
    
    for tile_dir in tiles_full_info:
        
        col,row,tw,th=tiles_full_info[tile_dir]

        vrt_row.setdefault(row,{'vrt_body' : "",'vrt_dir' : tile_dir.split('/')[0], 'th' : th})
        height_map = os.path.join(tile_dir,filename)
        s = height_map.split("/")
        height_map = os.path.join(*s[1:])

        if os.path.isfile(os.path.join(cfg['out_dir'],vrt_row[row]['vrt_dir'],height_map)):
            files_to_remove.append(os.path.join(cfg['out_dir'],vrt_row[row]['vrt_dir'],height_map))
            vrt_row[row]['vrt_body']+="\t\t<SimpleSource>\n"
            vrt_row[row]['vrt_body']+="\t\t\t<SourceFilename relativeToVRT=\"1\">%s</SourceFilename>\n" % height_map
            vrt_row[row]['vrt_body']+="\t\t\t<SourceBand>1</SourceBand>\n"
            vrt_row[row]['vrt_body']+="\t\t\t<SrcRect xOff=\"%i\" yOff=\"%i\" xSize=\"%i\" ySize=\"%i\"/>\n" % (0, 0, tw/z, th/z)
            vrt_row[row]['vrt_body']+="\t\t\t<DstRect xOff=\"%i\" yOff=\"%i\" xSize=\"%i\" ySize=\"%i\"/>\n" % (col/z, 0, tw/z, th/z)
            vrt_row[row]['vrt_body']+="\t\t</SimpleSource>\n"

    vrtfilename = fout
    vrtfile = open(vrtfilename, 'w')
    vrtfile.write("<VRTDataset rasterXSize=\"%i\" rasterYSize=\"%i\">\n" % (w/z,h/z))
    vrtfile.write("\t<VRTRasterBand dataType=\"Float32\" band=\"1\">\n")
    vrtfile.write("\t\t<ColorInterp>Gray</ColorInterp>\n")
    
    
    for row,vrt_data in vrt_row.iteritems():
        th = vrt_data['th']
        # First, write row vrt file
        col_vrt_filename = os.path.join(cfg['out_dir'],vrt_data['vrt_dir'],os.path.basename(vrtfilename))
        files_to_remove.append(col_vrt_filename)
        tmp_vrt_file = open(col_vrt_filename,'w')
        tmp_vrt_file.write("<VRTDataset rasterXSize=\"%i\" rasterYSize=\"%i\">\n" % (w/z,
                                                                            th/z))
        tmp_vrt_file.write("\t<VRTRasterBand dataType=\"Float32\" band=\"1\">\n")
        tmp_vrt_file.write("\t\t<ColorInterp>Gray</ColorInterp>\n")
        tmp_vrt_file.write(vrt_data['vrt_body'])
        tmp_vrt_file.write("\t</VRTRasterBand>\n")
        tmp_vrt_file.write("</VRTDataset>\n")
        tmp_vrt_file.close()
                
        # Next, write entry in final vrt file
        vrtfile.write("\t\t<SimpleSource>\n")
        vrtfile.write("\t\t\t<SourceFilename relativeToVRT=\"1\">%s</SourceFilename>\n" % os.path.join(vrt_data['vrt_dir'],os.path.basename(vrtfilename)))
        vrtfile.write("\t\t\t<SourceBand>1</SourceBand>\n")
        vrtfile.write("\t\t\t<SrcRect xOff=\"%i\" yOff=\"%i\" xSize=\"%i\" ySize=\"%i\"/>\n" % (0, 0, w/z, th/z))
        vrtfile.write("\t\t\t<DstRect xOff=\"%i\" yOff=\"%i\" xSize=\"%i\" ySize=\"%i\"/>\n" % (0, row/z, w/z, th/z))
        vrtfile.write("\t\t</SimpleSource>\n")

    vrtfile.write("\t</VRTRasterBand>\n")
    vrtfile.write("</VRTDataset>\n")
    vrtfile.close()

    if cfg['vrt_to_tiff']:
        common.run('gdal_translate %s %s' % (vrtfilename, os.path.splitext(vrtfilename)[0]+".tif"))

        if cfg['clean_intermediate']:
            for f in files_to_remove:
                common.remove_if_exists(f)
    return


def mosaic(fout, w, h, list_tiles, tw, th, ov):
    """
    Compose several tiles of the same size into a bigger image.

    Args:
        fout: path to the output image
        w, h: output image dimensions
        list_tiles: list containing paths to the input tiles
        tw, th: dimensions of a tile (they must all have the same dimensions)
        ov: overlap between tiles (in pixels)

    Returns:
        nothing
    """
    N = len(list_tiles)
    ntx = np.ceil(float(w - ov) / (tw - ov)).astype(int)
    nty = np.ceil(float(h - ov) / (th - ov)).astype(int)
    assert(ntx * nty == N)

    # default numpy datatype is float64, useless as the ouput file will be
    # stored with float32
    out = np.zeros([h, w], dtype=np.float32)
    count = np.zeros([h, w], dtype=np.uint8)

    # loop over all the tiles
    for j in range(nty):
        for i in range(ntx):
            sys.stdout.write("\tPasting tile %02d %02d\r" % (j, i))
            sys.stdout.flush()
            # top-left and bottom-right corners of the tile in the output full
            # image
            x0 = i * (tw - ov)
            y0 = j * (th - ov)
            x1 = min(x0 + tw, w)
            y1 = min(y0 + th, h)

            # read the tile with piio. If the tile has not been produced,
            # nothing needs to be done. The corresponding pixels will get the
            # value 'nan' in the output full image.
            tile_fname = list_tiles[j * ntx + i]
            if os.path.isfile(tile_fname):
                tile = piio.read(tile_fname).astype(np.float32)[:, :, 0]
                assert(np.shape(tile) == (th, tw))

                # count the pixels different from nan and inf
                ind = np.isfinite(tile)
                count[y0:y1, x0:x1] += ind[:y1 - y0, :x1 - x0]

                # replace nan and inf with zeros, then add the tile to the
                # output. ~ind is the negation of ind
                tile[~ind] = 0
                out[y0:y1, x0:x1] += tile[:y1 - y0, :x1 - x0]

    # free mem
    if 'tile' in locals():
        del tile
    if 'ind' in locals():
        del ind
    gc.collect()

    sys.stdout.write('\n')
    # put nan where count is zero, and take the average where count is nonzero.
    sys.stdout.write('\tCounting...\n')
    sys.stdout.flush()
    ind = (count > 0)

    sys.stdout.write('\tAveraging...\n')
    sys.stdout.flush()
    out[ind] /= count[ind]

    sys.stdout.write('\tPutting nans on empty pixels...\n')
    sys.stdout.flush()
    out[~ind] = np.nan

    del count

    # saving the 'out' numpy array in TIFF with piio requires to much memory
    # (something like twice de file size) because tiled tiff image writing is
    # not implemented yet in iio.
    # As an alternative, the numpy array is stored in raw and the libtiff util
    # 'raw2tiff' is used to produce a tiff file from it.
    sys.stdout.write('\twriting raw data to disk...\n')
    raw_file = common.tmpfile('')
    out.tofile(raw_file)
    common.run('raw2tiff -w %d -l %d -d float -c zip %s %s' % (w, h, raw_file,
                                                               fout))


    # sys.stdout.write('Writing output file...')
    # sys.stdout.flush()
    # piio.write(fout, out)



# How prepare the list of tiles to launch this function from ipython:
# import glob
# tiles = glob.glob('tile_*')
# tiles.sort()
# ind = []
# line = range(0, (nx-1)*ny+1, ny)
# for y in range(ny): ind = ind + [y+x for x in line]
# tiles_sorted = [os.path.join(tiles[i], 'height_map.tif') for i in ind]

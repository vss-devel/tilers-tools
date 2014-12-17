#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
# Copyright (c) 2011-2013 Vadim Shlyakhov
#
#  Permission is hereby granted, free of charge, to any person obtaining a
#  copy of this software and associated documentation files (the "Software"),
#  to deal in the Software without restriction, including without limitation
#  the rights to use, copy, modify, merge, publish, distribute, sublicense,
#  and/or sell copies of the Software, and to permit persons to whom the
#  Software is furnished to do so, subject to the following conditions:
#
#  The above copyright notice and this permission notice shall be included
#  in all copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
#  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
#  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#  DEALINGS IN THE SOFTWARE.
###############################################################################

from __future__ import print_function
import os
import os.path
import shutil
import math
import cgi
from PIL import Image

try:
    from osgeo import gdal
    from osgeo import osr
    from osgeo import ogr
    from osgeo.gdalconst import *
#    gdal.TermProgress = gdal.TermProgress_nocb
except ImportError:
    import gdal
    import osr
    import ogr
    from gdalconst import *

from tiler_functions import *
import map2gdal

profile_map = []

resampling_map = {
    'near':     Image.NEAREST,
    'nearest':  Image.NEAREST,
    'bilinear': Image.BILINEAR,
    'bicubic':  Image.BICUBIC,
    'antialias':Image.ANTIALIAS,
    }
def resampling_lst():
    return resampling_map.keys()

base_resampling_map = {
    'near':         'NearestNeighbour',
    'nearest':      'NearestNeighbour',
    'bilinear':     'Bilinear',
    'cubic':        'Cubic',
    'cubicspline':  'CubicSpline',
    'lanczos':      'Lanczos',
    }
def base_resampling_lst():
    return base_resampling_map.keys()

#############################

class TilingScheme(object):

#############################
    tile_origin_corner = 'tl'
    tile_size = (256, 256)
    axis_inv = (1, -1) # y goes downwards

    #----------------------------

    def normalize_tile(self, tile):
        'normalize according to the tile grid'
    #----------------------------
        z, x, y = tile
        ntiles_x, ntiles_y = self.n_tiles_xy(z)
        return (z, x % ntiles_x, y)

    #----------------------------

    def tile_path(self, tile):
        'relative path to a tile'
    #----------------------------
        z, x, y = self.normalize_tile(tile)
        return '%i/%i/%i%s' % (z, x, y, self.tile_ext)

class XYZtiling(TilingScheme):
    pass

class ZYXtiling(XYZtiling):

    #----------------------------

    def tile_path(self, tile):
        'relative path to a tile'
    #----------------------------
        z, x, y = self.normalize_tile(tile)
        return 'z%i/%i/%i%s' % (z, y, x, self.tile_ext)

class TMStiling(TilingScheme):
    axis_inv = (1, 1)
    tile_origin_corner = 'bl'

    #----------------------------

    def normalize_tile(self, tile):
        'normalize according to the tile grid'
    #----------------------------
        z, x, y = super(TMStiling, self).normalize_tile(tile)
        ntiles_x, ntiles_y = self.n_tiles_xy(zoom)
        return (z, x, ntiles_y - 1 - y)

#############################

class BaseImg(object):
    '''Tile feeder for a base zoom level'''
#############################

    def __init__(self, dataset, tl_offsets, transparency=None):
        self.ds = dataset
        self.tl_offsets = tl_offsets
        self.transparency = transparency

        self.size = self.ds.RasterXSize, self.ds.RasterYSize
        self.bands = [self.ds.GetRasterBand(i + 1) for i in range(self.ds.RasterCount)]

    def __del__(self):
        del self.bands
        del self.ds

    def get_tile(self, corners):
        '''crop raster as per pair of world pixel coordinates'''

        tl = [corners[0][c] - self.tl_offsets[c] for c in (0, 1)]
        sz = [corners[1][c] - corners[0][c] for c in (0, 1)]

        tile_bands = [bnd.ReadRaster(tl[0], tl[1], sz[0], sz[1], sz[0], sz[1], GDT_Byte)
                    for bnd in self.bands]
        n_bands = len(self.bands)
        if n_bands == 1:
            opacity = 1
            mode = 'L'
            if self.transparency is not None:
                if chr(self.transparency) in tile_bands[0]:
                    colorset = set(tile_bands[0])
                    if len(colorset) == 1:  # fully transparent
                        return None, 0
                    else:                   # semi-transparent
                        opacity = -1
            img = Image.frombuffer('L', sz, tile_bands[0], 'raw', 'L', 0, 1)
        else:
            aplpha = tile_bands[-1]
            if min(aplpha) == '\xFF':       # fully opaque
                opacity = 1
                tile_bands = tile_bands[:-1]
                mode = 'RGB' if n_bands > 2 else 'L'
            elif max(aplpha) == '\x00':     # fully transparent
                return None, 0
            else:                           # semi-transparent
                opacity = -1
                mode = 'RGBA' if n_bands > 2 else 'LA'
            img = Image.merge(mode, [Image.frombuffer('L', sz, bnd, 'raw', 'L', 0, 1) for bnd in tile_bands])
        return img, opacity
# BaseImg


#############################

class Pyramid(object):
    '''Tile pyramid generator and utilities'''
#############################

    zoom0_tiles = [1, 1] # tiles at zoom 0, default value

    palette = None
    transparency = None
    zoom_range = None
    min_res = None
    max_extent = None
    max_raster_origin = None

    default_zoom_range = (0, 22)

    #----------------------------

    def __init__(self, src=None, dest=None, options=None):

    #----------------------------
        gdal.UseExceptions()

        self.temp_files = []
        self.src = src
        self.dest = dest
        ld('src dest',src, dest)
        self.options = LooseDict(options)
        self.name = self.options.name
        self.tile_ext = self.options.tile_ext
        self.description = ''

        # self.proj_srs may be changed later, for example, to avoid crossing longitude 180
        self.proj_srs = txt2proj4(self.srs)
        self.geog_srs = proj_cs2geog_cs(self.proj_srs)
        ld('proj, longlat', self.proj_srs, self.geog_srs)

        self.proj2geog = GdalTransformer(SRC_SRS=self.proj_srs, DST_SRS=self.geog_srs)

        self.init_parameters()

    #----------------------------

    def __del__(self):

    #----------------------------
        try:
            if self.options.verbose < 2:
                for f in self.temp_files:
                    os.remove(f)
        except: pass

    #----------------------------

    def init_parameters(self):
        # init tile grid parameters
    #----------------------------

        self.tile_origin = self.get_corner_coords(self.tile_origin_corner)
        self.max_raster_origin = self.get_corner_coords('tl')

        # default raster corners to max_extent
        self.raster_corners = (
            self.get_corner_coords('tl'),
            self.get_corner_coords('br')
            )

        ld('zoom0_tiles', self.zoom0_tiles, 'tile_size', self.tile_size, 'max_raster_origin', self.max_raster_origin, 'tile_origin', self.tile_origin,)

    #----------------------------

    def open_source_dataset(self):
        ''
    #----------------------------
        self.src_dir, src_f = os.path.split(self.src)
        self.base = os.path.splitext(src_f)[0]

        if self.options.delete_src:
            self.temp_files.append(self.src)

        if os.path.isdir(self.dest):
            if self.options.noclobber and os.path.exists(self.dest):
                raise RuntimeError('Target already exists: skipping')
            else:
                shutil.rmtree(self.dest, ignore_errors=True)

        self.base_resampling = base_resampling_map[self.options.base_resampling]
        self.resampling = resampling_map[self.options.overview_resampling]

        self.src_path = self.src
        if os.path.exists(self.src):
            self.src_path = os.path.abspath(self.src)
            #~ pf('')
            ld('self.src_path',self.src_path, self.src)

        # check for source raster type
        self.src_ds = gdal.Open(self.src_path, GA_ReadOnly)
        self.description = self.src_ds.GetMetadataItem('DESCRIPTION')

        # source is successfully opened, then create destination dir
        os.makedirs(self.dest)

        self.modify_src_raster()

    #----------------------------

    def modify_src_raster(self):
        'convert to RGB(A) if required'
    #----------------------------

        override_srs = self.options.srs

        src_bands = self.src_ds.RasterCount
        band1 = self.src_ds.GetRasterBand(1)
        if src_bands == 1 and band1.GetColorInterpretation() == GCI_PaletteIndex : # source is a paletted raster
            transparency = None
            if self.base_resampling == 'NearestNeighbour' and self.resampling == Image.NEAREST :
                # check if src can be rendered in paletted mode
                color_table = band1.GetColorTable()
                ncolors = color_table.GetCount()
                palette = [color_table.GetColorEntry(i) for i in range(ncolors)]
                r, g, b, a = zip(*palette)
                pil_palette = flatten(zip(r, g, b)) # PIL doesn't support RGBA palettes
                if self.options.dst_nodata is not None:
                    transparency = int(self.options.dst_nodata.split(',')[0])
                elif min(a) == 0:
                    transparency = a.index(0)
                elif ncolors < 256:
                    pil_palette += [0, 0, 0] # the last color index is a transparency
                    transparency = len(pil_palette)/3-1

            ld('transparency', transparency)
            if transparency is not None: # render in paletted mode
                self.transparency = transparency
                self.palette = pil_palette
                ld('self.palette', self.palette)
            else: # convert src to rgb VRT

                src_geotr = self.src_ds.GetGeoTransform()
                src_proj = txt2proj4(self.src_ds.GetProjection())
                gcps = self.src_ds.GetGCPs()
                if gcps:
                    ld('src GCPsToGeoTransform', gdal.GCPsToGeoTransform(gcps))

                if not src_proj and gcps :
                    src_proj = txt2proj4(self.src_ds.GetGCPProjection())

                if override_srs is not None:
                    src_proj = txt2proj4(override_srs)
                    override_srs = None

                ld('src_proj', src_proj, 'src geotr', src_geotr)
                assert src_proj, 'The source does not have a spatial reference system assigned'

                if not src_geotr or src_geotr == (0.0, 1.0, 0.0, 0.0, 0.0, 1.0):
                    geotr_txt = ''
                else:
                    geotr_txt = geotr_templ % src_geotr

                gcplst_txt = ''
                if gcps:
                    gcp_lst = '\n'.join((gcp_templ % (g.Id, g.GCPPixel, g.GCPLine, g.GCPX, g.GCPY, g.GCPZ)
                                        for g in gcps))
                    gcp_proj = txt2proj4(self.src_ds.GetGCPProjection()) if override_srs is None else src_proj
                    gcplst_txt = gcplst_templ % (gcp_proj, gcp_lst)

                metadata = self.src_ds.GetMetadata()
                ld('metadata', metadata)
                if metadata:
                    mtd_lst = [xml_txt('MDI', metadata[mdkey], 4, key=mdkey) for mdkey in metadata]
                    meta_txt = meta_templ % '\n'.join(mtd_lst)
                else:
                    meta_txt = ''

                xsize, ysize = (self.src_ds.RasterXSize, self.src_ds.RasterYSize)
                blxsize, blysize = band1.GetBlockSize()

                band_lst = ''.join((band_templ % {
                    'band':     band,
                    'color':    color,
                    'src':      cgi.escape(self.src_path, quote=True),
                    'srcband':  1,
                    'xsize':    xsize,
                    'ysize':    ysize,
                    'blxsize':  blxsize,
                    'blysize':  blysize,
                    } for band, color in ((1, 'Red'), (2, 'Green'), (3, 'Blue'))))

                vrt_txt = vrt_templ % {
                    'xsize':    xsize,
                    'ysize':    ysize,
                    'metadata': meta_txt,
                    'srs':      (srs_templ % src_proj) if src_proj else '',
                    'geotr':    geotr_txt,
                    'gcp_list': gcplst_txt,
                    'band_list':band_lst,
                    }

                src_vrt = os.path.abspath(os.path.join(self.dest, self.base+'.src.vrt')) # auxilary VRT file

                self.temp_files.append(src_vrt)
                self.src_path = src_vrt
                with open(src_vrt, 'w') as f:
                    f.write(vrt_txt.encode('utf-8'))

                self.src_ds = gdal.Open(src_vrt, GA_ReadOnly)
                # finished with a paletted raster

        if override_srs is not None: # src SRS needs to be relpaced
            src_vrt = os.path.join(self.dest, self.base+'.src.vrt') # auxilary VRT file
            self.temp_files.append(src_vrt)
            self.src_path = src_vrt

            vrt_drv = gdal.GetDriverByName('VRT')
            self.src_ds = vrt_drv.CreateCopy(src_vrt, src_ds) # replace src dataset

            ld('override_srs', override_srs, 'txt2wkt(override_srs)', txt2wkt(override_srs))
            self.src_ds.SetProjection(txt2wkt(override_srs)) # replace source SRS
            gcps = self.src_ds.GetGCPs()
            if gcps :
                self.src_ds.SetGCPs(gcps, txt2wkt(override_srs))

    #----------------------------

    def init_output(self):
        'initialize geo-parameters and generate base zoom level'
    #----------------------------
        self.tiles_prefix = self.options.tiles_prefix

        #~ if self.options.verbose > 0:
            #~ print('\n%s -> %s '%(self.src, self.dest), end='')
        logging.info(' %s -> %s '%(self.src, self.dest))

        # get corners at the target SRS
        target_corners = self.auto_warp_corners()

        ld('max raster', self.raster_corners, self.proj2geog.transform(self.raster_corners))
        # self.raster_corners were set to the max raster, now clip them to the max tileset area
        self.raster_corners = (
            (max(self.raster_corners[0][0], target_corners[0][0]),
                min(self.raster_corners[0][1], target_corners[0][1])),
            (min(self.raster_corners[1][0], target_corners[1][0]),
                max(self.raster_corners[1][1], target_corners[1][1]))
            )

        ld('target_corners', target_corners, self.proj2geog.transform(target_corners))
        ld('target raster', self.raster_corners, self.proj2geog.transform(self.raster_corners))

        self.calc_zoom(target_corners)

    #----------------------------

    def calc_zoom(self, corners):
        'determine and set a list of zoom levels to generate'
    #----------------------------

        # check raster parameters to find default zoom range
        ld('automatic zoom levels')

        auto_warp_res = self.auto_warp_res(corners)
        ld('auto_warp_corners', corners, 'auto_warp_res', auto_warp_res)

        max_zoom = max(self.res2zoom_xy(auto_warp_res))

        tl = corners[0]
        br = corners[1]
        wh = (br[0] - tl[0], tl[1] - br[1])
        ld('calc_zoom tl, br, wh', tl, br, wh)
        min_zoom = min(self.res2zoom_xy([wh[i] / self.tile_size[i] for i in (0, 1)]))

        self.set_zoom_range(self.options.zoom, (min_zoom, max_zoom))

    #----------------------------

    def auto_warp_corners(self):

    #----------------------------

        width = self.src_ds.RasterXSize
        height = self.src_ds.RasterYSize

        n_samples = 100

        def chunks(length):
            return (int(round(i * length / n_samples)) for i in range(n_samples + 1))

        top_line = ((j, 0) for j in chunks(width))
        bottom_line = ((j, height) for j in chunks(width))

        left_line = ((0, j) for j in chunks(height))
        right_line = ((width, j) for j in chunks(height))

        def iter_transformer(*args):
            transformer = GdalTransformer(
                self.src_ds,
                DST_SRS=self.proj_srs,
                )
            for p in itertools.chain(*args):
                try:
                    yield transformer.transform_point(p)
                except RuntimeError:
                    yield None

        out_pts = filter(None, iter_transformer(top_line, bottom_line, left_line, right_line))
        #~ ld('out_pts', out_pts)

        xx, yy = zip(*out_pts)
        left = min(xx)
        right = max(xx)
        bottom = min(yy)
        top = max(yy)

        return ((left, top), (right, bottom))

    #----------------------------

    def auto_warp_res(self, corners):

    #----------------------------
        left, top = corners[0]
        right, bottom = corners[1]
        diag = math.sqrt((right - left) ** 2 + (top - bottom) ** 2)

        pix_diag = math.sqrt(self.src_ds.RasterXSize ** 2 + self.src_ds.RasterYSize ** 2)

        r = diag / pix_diag
        return (r, r)

    #----------------------------

    def create_target_dataset(self):

    #----------------------------

        # adjust raster extents to tile boundaries
        tile_tl, tile_br = self.corner_tiles(self.max_zoom)
        ld('base_raster')
        ld('tile_tl', tile_tl, 'tile_br', tile_br)
        tl_c = self.tile_corners(tile_tl)[0]
        br_c = self.tile_corners(tile_br)[1]

        # warp base raster
        base_ds = self.create_warped_vrt((tl_c, br_c), self.zoom2res(self.max_zoom))

        # close source dataset
        del self.src_ds

        tl_pix = self.tile_pixcorners(tile_tl)[0]
        # create base_image raster
        self.base_img = BaseImg(base_ds, tl_pix, self.transparency)

    #----------------------------

    def create_warped_vrt(self, corners=None, res=None):

    #----------------------------

        # generate warp transform
        src_geotr = self.src_ds.GetGeoTransform()
        src_proj = txt2proj4(self.src_ds.GetProjection())
        gcp_proj = None

        if not self.options.tps and src_geotr and src_geotr != (0.0, 1.0, 0.0, 0.0, 0.0, 1.0):
            ok, src_igeotr = gdal.InvGeoTransform(src_geotr)
            assert ok
            src_transform = '%s\n%s' % (warp_src_geotr % src_geotr, warp_src_igeotr % src_igeotr)
        else:
            gcps = self.src_ds.GetGCPs()
            assert gcps, 'Neither geotransform, nor gpcs are in the source file %s' % self.src

            gcp_lst = [(g.Id, g.GCPPixel, g.GCPLine, g.GCPX, g.GCPY, g.GCPZ) for g in gcps]
            ld('src_proj', self.src_ds.GetProjection(), 'gcp_proj', self.src_ds.GetGCPProjection())
            gcp_proj = txt2proj4(self.src_ds.GetGCPProjection())
            if src_proj and gcp_proj != src_proj:
                coords = GdalTransformer(
                    SRC_SRS=gcp_proj,
                    DST_SRS=src_proj
                    ).transform([g[3:6] for g in gcp_lst])

                gcp_lst = [tuple(p[:3] + c) for p, c in zip(gcp_lst, coords)]

            gcp_txt = '\n'.join((gcp_templ % g for g in gcp_lst))
            #src_transform = warp_src_gcp_transformer % (0, gcp_txt)
            src_transform = warp_src_tps_transformer % gcp_txt

        # base zoom level raster size
        tl_c, br_c = corners
        left, top  = tl_c
        right, bottom = br_c
        ld('right, left, top, bottom', right, left, top, bottom)
        size = (abs((right - left) / res[0]), abs((top - bottom) / res[0]))

        #tl_ll, br_ll = self.coords2longlat([tl_c, br_c])
        ld('create_target_dataset', 'max_zoom', self.max_zoom, 'size', size[0], size[1], '-tr', res[0], res[1], '-te', tl_c[0], br_c[1], br_c[0], tl_c[1], '-t_srs', self.proj_srs)

        dst_geotr = ( left, res[0], 0.0,
                      top, 0.0, -res[1] )
        ok, dst_igeotr = gdal.InvGeoTransform(dst_geotr)
        assert ok
        dst_transform = '%s\n%s' % (warp_dst_geotr % dst_geotr, warp_dst_igeotr % dst_igeotr)

        # generate warp options
        warp_options = []
        def w_option(name, value): # warp options template
            return '    <Option name="%s">%s</Option>' % (name, value)

        warp_options.append(w_option('INIT_DEST', 'NO_DATA'))

        # generate cut line
        if self.options.cut or self.options.cutline:
            cut_wkt = self.get_cutline()
        else:
            cut_wkt = None
        if cut_wkt:
            warp_options.append(w_option('CUTLINE', cut_wkt))
            if self.options.blend_dist:
                warp_options.append(w_option('CUTLINE_BLEND_DIST', self.options.blend_dist))

        src_bands = self.src_ds.RasterCount
        ld('src_bands', src_bands)

        # process nodata info
        src_nodata = None
        if self.options.src_nodata:
            src_nodata = map(int, self.options.src_nodata.split(','))
            assert len(src_nodata) == src_bands, 'Nodata must match the number of bands'
            if src_bands > 1:
                warp_options.append(w_option('UNIFIED_SRC_NODATA', 'YES'))
        dst_nodata = None
        if self.palette is not None:
            dst_nodata = [self.transparency]
        ld('nodata', src_nodata, dst_nodata)

        # src raster bands mapping
        vrt_bands = []
        wo_BandList = []
        for i in range(src_bands):
            vrt_bands.append(warp_band % (i + 1, '/'))
            if src_nodata or dst_nodata:
                band_mapping_info = warp_band_mapping_nodata % (
                        warp_band_src_nodata % (src_nodata[i], 0) if src_nodata else '',
                        warp_band_dst_nodata % (dst_nodata[i], 0) if dst_nodata else '')
            else:
                band_mapping_info = '/'
            wo_BandList.append(warp_band_mapping % (i + 1, i + 1, band_mapping_info))

        if src_bands < 4 and self.palette is None:
            vrt_bands.append(warp_band % (src_bands + 1, warp_band_color % 'Alpha'))

        vrt_text = warp_vrt % {
            'xsize':            size[0],
            'ysize':            size[1],
            'srs':              self.proj_srs,
            'geotr':            geotr_templ % dst_geotr,
            'band_list':        '\n'.join(vrt_bands),
            'blxsize':          self.tile_size[0],
            'blysize':          self.tile_size[1],
            'wo_ResampleAlg':   self.base_resampling,
            'wo_src_path':      cgi.escape(self.src_path, quote=True),
            'warp_options':     '\n'.join(warp_options),
            'wo_src_srs':       gcp_proj if gcp_proj else src_proj,
            'wo_dst_srs':       self.proj_srs,
            'wo_src_transform': src_transform,
            'wo_dst_transform': dst_transform,
            'wo_BandList':      '\n'.join(wo_BandList),
            'wo_DstAlphaBand':  warp_dst_alpha_band % (src_bands + 1) if src_bands < 4  and self.palette is None else '',
            'wo_Cutline':       (warp_cutline % cut_wkt) if cut_wkt else '',
            }

        temp_vrt = os.path.join(self.dest, self.base + '.tmp.vrt') # auxilary VRT file
        self.temp_files.append(temp_vrt)
        with open(temp_vrt, 'w') as f:
            f.write(vrt_text.encode('utf-8'))

        return gdal.Open(vrt_text, GA_ReadOnly)

    #----------------------------

    def get_cutline(self):

    #----------------------------
        cutline = self.src_ds.GetMetadataItem('CUTLINE')
        ld('cutline', cutline)
        if cutline and not self.options.cutline:
            return cutline

        # try to find an external cut line
        if self.options.cutline:
            cut_file = self.options.cutline
        else: # try to find a file with a cut shape
            for ext in ('.gmt', '.shp', '.kml'):
                cut_file = os.path.join(self.src_dir, self.base+ext)
                if os.path.exists(cut_file):
                    break
            else:
                return None

        feature_name = self.base if self.options.cutline_match_name else None
        return shape2cutline(cut_file, self.src_ds, feature_name)

    #----------------------------

    def generate_tiles(self):
        'generate tiles'
    #----------------------------

        # connect to src dataset
        try:
            self.open_source_dataset()
        except RuntimeError as exc:
            if self.options.skip_invalid:
                logging.error(exc.message)
                return
            else:
                raise

        self.init_output()

        # create a raster source for a base zoom
        self.create_target_dataset()

        if not self.name:
            self.name = os.path.basename(self.dest)

        ld('generate tiles')

        self.progress()

        top_results = filter(None, itertools.imap(self.make_tile_raster, self.get_top_tiles()))

        self.progress(finished=True)

        children, images, opacities = zip(*top_results)
        # write top-level metadata (html/kml)
        self.write_metadata(None, children)

        # cache back tiles transparency
        transparency = dict((
            (self.tile_path(tile), opc)
                for tile, opc in itertools.chain(*opacities)
            ))
        write_transparency(self.dest, transparency)

    #----------------------------

    def get_top_tiles(self):

    #----------------------------
        min_zoom = self.zoom_range[-1]
        tile_tl, tile_br = self.corner_tiles(min_zoom)
        xx = (tile_tl[1], tile_br[1])
        yy = (tile_tl[2], tile_br[2])
        return ((min_zoom, x, y) for y in range(min(yy), max(yy)+1) for x in range(min(xx), max(xx)+1))

    #----------------------------

    def make_tile_raster(self, tile):

    #----------------------------

        if not self.in_range(tile, check_zoom=False):
            return

        zoom, x, y = tile
        if zoom == self.max_zoom: # get from the base image
            tile_img, opacity = self.base_img.get_tile(self.tile_pixcorners(tile))
            opacity_lst = [(tile, opacity)]
        else: # merge children
            tile_img, opacity_lst = self.assemble_tile(tile)

        #~ ld('make_tile_raster', tile, tile_img, opacity)
        if tile_img is not None and self.zoom_in_range(zoom):
            if self.palette:
                tile_img.putpalette(self.palette)

            self.write_tile(tile, tile_img)

            # write tile-level metadata (html/kml)
            children = [ch for ch, opacity in opacity_lst[1:]]
            self.write_metadata(tile, children)

            return tile, tile_img, opacity_lst

    #----------------------------

    def assemble_tile(self, tile):

    #----------------------------

        zoom, x, y = tile
        ch_zoom = self.zoom_range[self.zoom_range.index(zoom) - 1] # children's zoom

        # map children locations inside the parent raster
        len_xy = int(2 ** (ch_zoom - zoom))
        children_map = dict(
            (((ch_zoom, x * len_xy + i, y * len_xy + j), # child tile
                (self.tile_size[0] * i, self.tile_size[1] * j)) # offset inside the parent tile
            for i in range(len_xy) for j in range(len_xy))
            )
        #ld(tile, ch_mozaic)

        ch_results = filter(None, itertools.imap(self.make_tile_raster, children_map.keys()))
        #~ ld('tile', tile, 'children', children, 'ch_results', ch_results)

        opacity = 0
        # combine into the parent tile
        if (len(ch_results) == len_xy * len_xy and
                all([opacity_data[0][1] == 1 for ch, img, opacity_data in ch_results])
            ):
            opacity = 1
            mode_opacity = ''
        else:
            opacity = -1
            mode_opacity = 'A'

        tile_img = None
        opacity_lst = []
        for ch, ch_img, opacity_data in ch_results:
            opacity_lst.extend(opacity_data)

            ch_mask = ch_img.split()[-1] if 'A' in ch_img.mode else None

            if tile_img is None:
                if 'P' in ch_img.mode:
                    tile_mode = 'P'
                elif 'L' in ch_img.mode:
                    tile_mode = 'L' + mode_opacity
                else:
                    tile_mode = 'RGB' + mode_opacity

                img_size = [i * len_xy for i in self.tile_size]
                if self.transparency is not None:
                    tile_img = Image.new(tile_mode, img_size, self.transparency)
                else:
                    tile_img = Image.new(tile_mode, img_size)

            tile_img.paste(ch_img, children_map[ch], ch_mask)

        opacity_lst.insert(0, (tile, opacity))

        return None if not tile_img else tile_img.resize(self.tile_size, self.resampling), opacity_lst

    #----------------------------

    def write_tile(self, tile, tile_img):

    #----------------------------
        rel_path = self.tile_path(tile)
        full_path = os.path.join(self.dest, rel_path)
        try:
            os.makedirs(os.path.dirname(full_path))
        except:
            pass

        tile_format = self.options.tile_format
        if self.options.paletted and tile_format == 'png':
            try:
                tile_img = tile_img.convert('P', palette=Image.ADAPTIVE, colors=255)
            except ValueError:
                #ld('tile_img.mode', tile_img.mode)
                pass
        elif tile_img.mode == 'P' and tile_format in ('jpeg', 'webp'):
            mode = 'RGB' # + 'A' if self.transparency else ''
            try:
                tile_img = tile_img.convert(mode)
            except ValueError:
                #ld('tile_img.mode', tile_img.mode)
                pass

        if self.transparency is not None:
            tile_img.save(full_path, transparency=self.transparency)
        else:
            tile_img.save(full_path)

        self.progress()

    #----------------------------

    def write_metadata(self, tile=None, children=[]):

    #----------------------------
        if tile is None:
            self.write_tilemap()

    #----------------------------

    def write_tilemap(self):
        '''Generate JSON for a tileset description'''
    #----------------------------

        # reproject extents back to the unshifted SRS
        bbox = GdalTransformer(SRC_SRS=self.proj_srs, DST_SRS=self.srs).transform(self.raster_corners)

        tile_mime = mime_from_ext(self.tile_ext)
        tilemap = {
            'type': 'TileMap',
            'properties': {
                'title':        self.name,
                'description':  self.description,
                },
            'tiles': {
                'size':         self.tile_size,
                'inversion':    self.axis_inv,
                'ext':          self.tile_ext[1:],
                'mime':         tile_mime,
                'origin':       self.tile_origin,
                'max_extent':   self.max_extent
                },
            'bbox': (
                bbox[0][0],
                bbox[1][1],
                bbox[1][0],
                bbox[0][1]),
            'crs': {
                "type": "name",
                "properties": {
                    "name": self.tilemap_crs,
                    }
                },
            'tilesets': dict([
                (zoom,
                    {"href": 'z%d' % zoom,
                    "units_per_pixel": self.zoom2res(zoom)[0]})
                for zoom in reversed(self.zoom_range)]),
            }


        write_tilemap(self.dest, tilemap)
        #~ ld('tilemap', tilemap)

    #----------------------------
    #
    # utility functions
    #
    #----------------------------

    @staticmethod
    def profile_class(profile_name):
        for cls in profile_map:
            if cls.profile == profile_name:
                return cls
        else:
            raise Exception("Invalid profile: %s" % profile_name)

    @staticmethod
    def profile_lst(tty=False):
        if not tty:
            return [c.profile for c in profile_map]
        print('\nOutput profiles and compatibility:\n')
        [print('%10s - %s' % (c.profile, c.__doc__)) for c in profile_map]
        print()

    def get_corner_coords(self, corner):
        corner_map = {
            'l': 0,
            'b': 1,
            'r': 2,
            't': 3,
            }
        return [self.max_extent[corner_map[corner[a]]] for a in (1, 0)]

    def zoom2res(self, zoom):
        return [self.min_res[c]/2**zoom for c in (0, 1)]

    def res2zoom_xy(self, res):
        '''resolution to zoom levels (separate for x and y)'''
        z = [int(math.floor(math.log(self.min_res[c]/res[c], 2))) for c in (0, 1)]
        return [v if v>0 else 0 for v in z]

    def coord2tile(self, zoom, coord):
        '''cartesian coordinates to tile numbers'''
        return self.pix2tile(zoom, self.coord2pix(zoom, coord))

    def tile_pixcorners(self, tile):
        '''pixel coordinates of a tile'''
        tl_br = [
            [(tile[1+c] + t) * self.tile_size[c]
                for c in (0, 1)]
                    for t in (0, 1)]
        return tl_br

    def tile_corners(self, tile):
        '''cartesian coordinates of a tile's corners'''
        zoom = tile[0]
        coord_corners = [self.pix2coord(zoom, pc) for pc in self.tile_pixcorners(tile)]
        return coord_corners

    def coord2pix(self, zoom, xy):
        '''cartesian coordinates to pixel coordinates'''
        res = self.zoom2res(zoom)
        return [int(round(
            (xy[c] - self.max_raster_origin[c]) / res[c] * (-1 if c else 1)
            )) for c in (0, 1)]

    def pix2coord(self, zoom, pix_coord):
        res = self.zoom2res(zoom)
        return [
            self.max_raster_origin[c] + pix_coord[c] * res[c] * (-1 if c else 1)
                for c in (0, 1)
            ]

    def pix2tile(self, zoom, pix):
        '''pixel coordinates to tile (z, x, y)'''
        tile = [int(pix[c] // self.tile_size[c]) for c in (0, 1)]
        tile.insert(0, zoom)
        return tile

    def n_tiles_xy(self, zoom):
        '''number of tiles along X and Y axes'''
        return map(lambda v: v * 2**zoom, self.zoom0_tiles)

    def coords2longlat(self, coords):
        longlat = [i[:2] for i in self.proj2geog.transform(coords)]
        return longlat

    def corners_lst2longlat(self, box_lst):
        deg_lst = self.coords2longlat(flatten(box_lst))
        tl_lst = deg_lst[0::2]
        br_lst = deg_lst[1::2]
        return [[
            (tl[0] if tl[0] <  180 else tl[0] - 360, tl[1]),
            (br[0] if br[0] > -180 else br[0] + 360, br[1]),
            ] for tl, br in zip(tl_lst, br_lst)]

    def corner_tiles(self, zoom):
        p_tl = self.coord2pix(zoom, self.raster_corners[0])
        t_tl = self.pix2tile(zoom, (p_tl[0], p_tl[1]))

        p_br = self.coord2pix(zoom, self.raster_corners[1])
        t_br = self.pix2tile(zoom, (p_br[0], p_br[1]))

        box_tl, box_br = [self.tile_corners(t) for t in (t_tl, t_br)]
        return t_tl, t_br

    def set_zoom_range(self, zoom_parm, default_range=None):
        'set a list of zoom levels from a parameter list'

        if not default_range:
            default_range = self.default_zoom_range
        if not zoom_parm:
            zoom_parm = '%d:%d' % default_range

        zchunk_lst = [z.split(':') for z in zoom_parm.split(',')]
        zlist = []
        for zchunk in zchunk_lst:
            if len(zchunk) == 1:
                zlist.append(int(zchunk[0]))
            else:
                # calculate zoom range
                zrange = []
                for n, d in zip(zchunk, default_range):
                    if n == '':              # set to default
                        z = d
                    elif n.startswith('-'): # set to default - n
                        z = d-int(n[1:])
                    elif n.startswith('+'): # set to default + n
                        z = d+int(n[1:])
                    else:                   # set to n
                        z = int(n)
                    zrange.append(z)

                # update range list
                zlist += range(min(zrange), max(zrange)+1)

        self.zoom_range = list(reversed(sorted(set(zlist))))
        self.max_zoom = self.zoom_range[0]
        ld('zoom_range', self.zoom_range, default_range)

    def zoom_in_range(self, zoom):
        return not self.zoom_range or zoom in self.zoom_range

    def in_range(self, tl_tile, br_tile=None, check_zoom=True):
        if not tl_tile:
            return False
        if not br_tile:
             br_tile = tl_tile

        # Y axis goes downwards
        zoom, tile_xmin, tile_ymin = tl_tile
        zoom, tile_xmax, tile_ymax = br_tile

        if check_zoom and not self.zoom_in_range(zoom):
            return False

        tl_zoom, br_zoom = self.corner_tiles(zoom)

        # Y axis goes downwards
        z, zoom_xmin, zoom_ymin = tl_zoom
        z, zoom_xmax, zoom_ymax = br_zoom

        out = not (
            tile_xmin > zoom_xmax or tile_xmax < zoom_xmin or
            tile_ymin > zoom_ymax or tile_ymax < zoom_ymin
            )

        #~ ld('in_range zoom', tl_zoom, br_zoom)
        #~ ld('in_range tile', tl_tile, br_tile, res)

        return out

    def set_region(self, point_lst, source_srs=None):
        if source_srs and source_srs != self.proj_srs:
            point_lst = GdalTransformer(SRC_SRS=source_srs, DST_SRS=self.proj_srs).transform(point_lst)

        x_coords, y_coords = zip(*point_lst)[0:2]
        top_left = min(x_coords), max(y_coords)
        bottom_right = max(x_coords), min(y_coords)
        self.raster_corners = [top_left, bottom_right]

    def load_region(self, datasource):
        if not datasource:
            return
        point_lst = flatten(shape2mpointlst(datasource, self.proj_srs))
        #~ ld(datasource, point_lst)
        self.set_region(point_lst)

    # progress display
    tick_rate = 50
    count = 0
    def progress(self, finished=False):
        #~ pf('+', end='')
        #~ return
        if self.options.verbose == 0:
            pass
        elif finished:
            pf('')
        elif self.count % self.tick_rate == 0:
            pf('.', end='')
        self.count += 1

# Pyramid

#############################

class MercatorPyramid(Pyramid):
    '''Mercator specifics for the Pyramid class'''
#############################

    #----------------------------

    def init_parameters(self):

    #----------------------------

        max_x = self.proj2geog.transform_point((180, 0), inv=True)[0] # Equator's half length
        ld('max_x', max_x)

        # pixel resolution at the top of the pyramid
        min_res = max_x * 2 / (self.zoom0_tiles[0] * self.tile_size[0])
        self.min_res = [min_res, min_res] # pixel 'y' goes downwards

        # top left corner of a world raster
        tl = (-max_x, self.min_res[1] * self.tile_size[1] * self.zoom0_tiles[1] / 2)
        br = [-tl[0], -tl[1]]

        self.max_extent = (tl[0], br[1], br[0], tl[1])
        ld('max extent', self.max_extent)

        super(MercatorPyramid, self).init_parameters()

    #----------------------------

    def init_output(self):
        'initialize geo-parameters and generate base zoom level'
    #----------------------------
        self.modify_srs()
        super(MercatorPyramid, self).init_output()

    #----------------------------

    def modify_srs(self):
        'change prime meridian to allow charts crossing 180 meridian'
    #----------------------------

        def srs_replace(proj4_txt, parm_lst):
            parm_dict = dict(((p.split('=')[0], p) for p in parm_lst))
            proj_lst = proj4_txt.split()
            proj_new = filter(None,
                (i if i.split('=')[0] not in parm_dict else None for i in proj_lst)
                )
            proj_new.extend(parm_lst)
            return ' '.join(proj_new)
        # srs_replace

        tl, br = GdalTransformer(
            self.src_ds,
            DST_SRS=self.geog_srs
            ).transform([
            (0, 0),
            (self.src_ds.RasterXSize, self.src_ds.RasterYSize)])
        ld('shift_srs tl', tl, 'br', br)

        l_lon = tl[0]
        r_lon = br[0]
        if r_lon <= 180 and l_lon >= -180 and l_lon < r_lon:
            return

        while l_lon > r_lon:
            r_lon += 360
        lon_0 = l_lon + (r_lon - l_lon) / 2
        while lon_0 > 180:
            lon_0 -= 360
        new_parms = ['+lon_0=%r' % lon_0]

        new_srs = srs_replace(self.proj_srs, new_parms)
        ld( 'lon_0', lon_0, 'l_lon', l_lon, 'r_lon', r_lon)

        old2new = GdalTransformer(
            SRC_SRS=self.proj_srs,
            DST_SRS=new_srs
            )
        shift_x = old2new.transform_point((0, 0))[0]

        if shift_x != 0:
            self.proj_srs = new_srs
            self.proj2geog = GdalTransformer(
                SRC_SRS=self.proj_srs,
                DST_SRS=self.geog_srs
                )
            self.max_raster_origin = (
                self.max_raster_origin[0] + shift_x,
                self.max_raster_origin[1]
                )
            ld('new_srs', new_srs, 'shift_x', shift_x, 'max_raster_origin', self.max_raster_origin)

#----------------------------
#
# templates for VRT XML
#
#----------------------------

def xml_txt(name, value=None, indent=0, **attr_dict):
    attr_txt = ''.join((' %s="%s"' % (key, attr_dict[key]) for key in attr_dict))
    val_txt = ('>%s</%s' % (cgi.escape(value, quote=True), name)) if value else '/'
    return '%s<%s%s%s>' % (' '*indent, name, attr_txt, val_txt)

warp_vrt = '''<VRTDataset rasterXSize="%(xsize)d" rasterYSize="%(ysize)d" subClass="VRTWarpedDataset">
  <SRS>%(srs)s</SRS>
%(geotr)s%(band_list)s
  <BlockXSize>%(blxsize)d</BlockXSize>
  <BlockYSize>%(blysize)d</BlockYSize>
  <GDALWarpOptions>
    <!-- <WarpMemoryLimit>6.71089e+07</WarpMemoryLimit> -->
    <ResampleAlg>%(wo_ResampleAlg)s</ResampleAlg>
    <WorkingDataType>Byte</WorkingDataType>
    <SourceDataset relativeToVRT="0">%(wo_src_path)s</SourceDataset>
%(warp_options)s
    <Transformer>
      <ApproxTransformer>
        <MaxError>0.125</MaxError>
        <BaseTransformer>
          <GenImgProjTransformer>
%(wo_src_transform)s
%(wo_dst_transform)s
            <ReprojectTransformer>
              <ReprojectionTransformer>
                <SourceSRS>%(wo_src_srs)s</SourceSRS>
                <TargetSRS>%(wo_dst_srs)s</TargetSRS>
              </ReprojectionTransformer>
            </ReprojectTransformer>
          </GenImgProjTransformer>
        </BaseTransformer>
      </ApproxTransformer>
    </Transformer>
    <BandList>
%(wo_BandList)s
    </BandList>
%(wo_DstAlphaBand)s%(wo_Cutline)s  </GDALWarpOptions>
</VRTDataset>
'''
warp_band = '  <VRTRasterBand dataType="Byte" band="%d" subClass="VRTWarpedRasterBand"%s>'
warp_band_color = '>\n    <ColorInterp>%s</ColorInterp>\n  </VRTRasterBand'
warp_dst_alpha_band = '    <DstAlphaBand>%d</DstAlphaBand>\n'
warp_cutline = '    <Cutline>%s</Cutline>\n'
warp_dst_geotr = '            <DstGeoTransform> %r, %r, %r, %r, %r, %r</DstGeoTransform>'
warp_dst_igeotr = '            <DstInvGeoTransform> %r, %r, %r, %r, %r, %r</DstInvGeoTransform>'
warp_src_geotr = '            <SrcGeoTransform> %r, %r, %r, %r, %r, %r</SrcGeoTransform>'
warp_src_igeotr = '            <SrcInvGeoTransform> %r, %r, %r, %r, %r, %r</SrcInvGeoTransform>'
warp_band_mapping = '      <BandMapping src="%d" dst="%d"%s>'
warp_band_src_nodata = '''
        <SrcNoDataReal>%d</SrcNoDataReal>
        <SrcNoDataImag>%d</SrcNoDataImag>'''
warp_band_dst_nodata = '''
        <DstNoDataReal>%d</DstNoDataReal>
        <DstNoDataImag>%d</DstNoDataImag>'''
warp_band_mapping_nodata = '''>%s%s
      </BandMapping'''
warp_src_gcp_transformer = '''            <SrcGCPTransformer>
              <GCPTransformer>
                <Order>%d</Order>
                <Reversed>0</Reversed>
                <GCPList>
%s
                </GCPList>
              </GCPTransformer>
            </SrcGCPTransformer>'''
warp_src_tps_transformer = '''            <SrcTPSTransformer>
              <TPSTransformer>
                <Reversed>0</Reversed>
                <GCPList>
%s
                </GCPList>
              </TPSTransformer>
            </SrcTPSTransformer>'''

gcp_templ = '    <GCP Id="%s" Pixel="%r" Line="%r" X="%r" Y="%r" Z="%r"/>'
gcplst_templ = '  <GCPList Projection="%s">\n%s\n  </GCPList>\n'
geotr_templ = '  <GeoTransform> %r, %r, %r, %r, %r, %r</GeoTransform>\n'
meta_templ = '  <Metadata>\n%s\n  </Metadata>\n'
band_templ = '''  <VRTRasterBand dataType="Byte" band="%(band)d">
    <ColorInterp>%(color)s</ColorInterp>
    <ComplexSource>
      <SourceFilename relativeToVRT="0">%(src)s</SourceFilename>
      <SourceBand>%(srcband)d</SourceBand>
      <SourceProperties RasterXSize="%(xsize)d" RasterYSize="%(ysize)d" DataType="Byte" BlockXSize="%(blxsize)d" BlockYSize="%(blysize)d"/>
      <SrcRect xOff="0" yOff="0" xSize="%(xsize)d" ySize="%(ysize)d"/>
      <DstRect xOff="0" yOff="0" xSize="%(xsize)d" ySize="%(ysize)d"/>
      <ColorTableComponent>%(band)d</ColorTableComponent>
    </ComplexSource>
  </VRTRasterBand>
'''
srs_templ = '  <SRS>%s</SRS>\n'
vrt_templ = '''<VRTDataset rasterXSize="%(xsize)d" rasterYSize="%(ysize)d">
%(metadata)s%(srs)s%(geotr)s%(gcp_list)s%(band_list)s</VRTDataset>
'''

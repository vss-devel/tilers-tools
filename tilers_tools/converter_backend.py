#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
# Copyright (c) 2010, 2013 Vadim Shlyakhov
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

import os
import os.path
import glob
import shutil
import json
from struct import Struct

from tiler_functions import *
from gdal_tiler import Pyramid

#############################

class Tile(object):

#############################
    def __init__(self, coord):
        self._coord = coord

    def coord(self):
        return self._coord

    def get_mime(self):
        return mime_from_ext(self.get_ext())

#############################

class FileTile(Tile):

#############################
    def __init__(self, coord, path):
        super(FileTile, self).__init__(coord)
        self.path = path

    def data(self):
        return open(self.path, 'rb').read()

    def get_ext(self):
        return os.path.splitext(self.path)[1]

    def copy2file(self, dst, link=False):
        if link and os.name == 'posix':
            dst_dir = os.path.split(dst)[0]
            src = os.path.relpath(self.path, dst_dir)
            os.symlink(src, dst)
        else:
            shutil.copy(self.path, dst)

#############################

class FileTileNoExt(FileTile):

#############################
    def get_ext(self):
        return ext_from_file(self.path)

#############################

class PixBufTile(Tile):

#############################
    def __init__(self, coord, pixbuf, key):
        super(PixBufTile, self).__init__(coord)
        self.pixbuf = pixbuf
        self.path = repr(key) # only for debugging

    def data(self):
        return self.pixbuf

    def get_ext(self):
        try:
            ext = ext_from_buffer(self.pixbuf)
        except KeyError:
            error('PixBufTile: wrong data', self.coord())
            raise
        return ext

    def copy2file(self, dest_path, link=False):
        open(dest_path, 'wb').write(self.pixbuf)

tile_profiles = []

#############################

class TileSet(object):

#############################
    def __init__(self, root=None, options=None, src=None):
        options = LooseDict(options)
        options.write = src

        self.root = root
        self.options = options
        self.src = src

        self.srs = self.options.proj4def or self.options.tiles_srs
        self.tilemap_crs = self.options.tiles_srs or self.tilemap_crs
        self.options.tiles_srs = self.srs

        self.zoom_levels = {}
        self.pyramid = Pyramid.profile_class('generic')(options=options)

        if not self.options.write:
            assert os.path.exists(root), 'No file or directory found: %s' % root
            if self.options.region:
                self.pyramid.set_zoom_range(self.options.zoom)
                self.pyramid.load_region(self.options.region)

        else:
            basename = os.path.splitext(os.path.basename(self.root or src.root))[0]
            df_name = os.path.splitext(basename)[0]
            if self.options.region:
                df_name += '-' + os.path.splitext(self.options.region)[0]
            self.name = self.options.name or df_name

            if not self.root:
                self.root = os.path.join(options.dst_dir, self.name + self.ext)

            if os.path.exists(self.root):
                if self.options.remove_dest:
                    if os.path.isdir(self.root):
                        shutil.rmtree(self.root, ignore_errors=True)
                    else:
                        os.remove(self.root)
                else:
                    assert self.options.append, 'Destination already exists: %s' % root

    @staticmethod
    def get_class(profile, write=False):
        for cls in tile_profiles:
            if profile == cls.format and ((not write and cls.input) or (write and cls.output)):
                return cls
        else:
            raise Exception('Invalid format: %s' % profile)

    @staticmethod
    def list_profiles():
        for cl in tile_profiles:
            print '%10s\t%s%s\t%s' % (
                cl.format,
                'r' if cl.input else ' ',
                'w' if cl.output else ' ',
                cl.__doc__
                )

    def in_range(self, ul_coords, lr_coords=None):
        if not ul_coords:
            return False
        if not self.pyramid:
            return True
        return self.pyramid.in_range(ul_coords, lr_coords)

    def __del__(self):
        log('self.count', self.count)

    def __iter__(self): # to be defined by a child
        raise Exception('Unimplemented!')

    def convert(self):
        pf('%s -> %s ' % (self.src.root, self.root), end='')
        map(self.process_tile, self.src)
        self.finalize_pyramid()
        self.finalize_tileset()
        pf('')

    def process_tile(self, tile):
        #log('process_tile', tile)
        self.store_tile(tile)

        # collect min max values for tiles processed
        zxy = list(tile.coord())
        z = zxy[0]

        min_max = self.zoom_levels.get(z, []) # min, max
        zzz, xxx, yyy = zip(*(min_max+[zxy]))
        self.zoom_levels[z] = [[z, min(xxx), min(yyy)], [z, max(xxx), max(yyy)]]

    def finalize_pyramid(self):
        log('self.zoom_levels', self.zoom_levels)

        # compute "effective" covered area
        prev_sq = 0
        for z in reversed(sorted(self.zoom_levels)):
            ul_zxy, lr_zxy = self.zoom_levels[z]
            ul_c = self.pyramid.tile_bounds(ul_zxy)[0]
            lr_c = self.pyramid.tile_bounds(lr_zxy)[1]
            sq = (lr_c[0]-ul_c[0])*(ul_c[1]-lr_c[1])
            area_diff = round(prev_sq/sq, 5)
            log('ul_c, lr_c', z, ul_c, lr_c, sq, area_diff)
            if area_diff == 0.25:
                break # this must be an exact zoom of a previous level
            area_coords = [ul_c, lr_c]
            prev_sq = sq

        self.pyramid.set_region(area_coords)
        self.pyramid.set_zoom_range(','.join(map(str, self.zoom_levels.keys())))

        self.pyramid.name = self.name

    def finalize_tileset(self):
        pass

    count = 0
    tick_rate = 100
    def counter(self):
        self.count+=1
        if self.count % self.tick_rate == 0:
            pf('.', end='')
            return True
        else:
            return False

# TileSet

#############################

class TileDir(TileSet):

#############################
    tile_class = FileTile

    def __init__(self, *args, **kw_args):
        super(TileDir, self).__init__(*args, **kw_args)

        if self.options.write:
            try:
                os.makedirs(self.root)
            except os.error: pass

    def __iter__(self):
        for f in glob.iglob(os.path.join(self.root, self.dir_pattern)):
            coord = self.path2coord(f)
            if not self.in_range(coord):
                continue
            self.counter()
            yield self.tile_class(coord, f)

    def path2coord(self, tile_path):
        raise Exception('Unimplemented!')

    def coord2path(self, z, x, y):
        raise Exception('Unimplemented!')

    def dest_ext(self, tile):
        return tile.get_ext()

    def store_tile(self, tile):
        try:
            tile_ext = self.dest_ext(tile)
            self.tile_ext = tile_ext
        except KeyError:
            tile_ext = '.xxx' # invalid file type
        dest_path = os.path.join(self.root, self.coord2path(*tile.coord())) + tile_ext
        log('%s -> %s' % (tile.path, dest_path))
        try:
            os.makedirs(os.path.split(dest_path)[0])
        except os.error: pass
        tile.copy2file(dest_path, self.options.link)
        self.counter()
# TileDir

#############################

class TileMapDir(TileDir):

#############################

    def finalize_tileset(self):
        self.pyramid.tile_ext = self.tile_ext
        self.pyramid.dest = self.root
        self.pyramid.write_metadata()

# TileMapDir

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
        ext = ext_from_buffer(self.pixbuf)
        return ext

    def copy2file(self, dest_path, link=False):
        open(dest_path, 'wb').write(self.pixbuf)

tile_formats = []

#############################

class TileSet(object):

#############################
    def __init__(self, root, options=None):
        options = options or LooseDict()
        log('root', root, options)
        self.root = root
        self.options = options
        self.zoom_levels = {}
        self.pyramid = None

        if not self.options.write:
            assert os.path.exists(root), 'No file or directory found: %s' % root
            if self.options.region:
                self.pyramid = Pyramid.profile_class('zyx')()
                self.pyramid.set_zoom_range(self.options.zoom)
                self.pyramid.load_region(self.options.region)
        else:
            if os.path.exists(self.root):
                if self.options.remove_dest:
                    if os.path.isdir(self.root):
                        shutil.rmtree(self.root, ignore_errors=True)
                    else:
                        os.remove(self.root)
                else:
                    assert self.options.append, 'Destination already exists: %s' % root

    def in_range(self, coords, range_lr_coords=None):
        if not coords:
            return False
        if not self.pyramid:
            return True
        return self.pyramid.in_range(coords, range_lr_coords)

    def __del__(self):
        log('self.count', self.count)

    def __iter__(self): # to be defined by a child
        raise Exception('Unimplemented!')

    def load_from(self, src_tiles):
        log('load_from', (src_tiles.root, self.root))
        map(self.process_tile, src_tiles)

    def process_tile(self, tile):
        #log('process_tile', tile)
        self.store_tile(tile)

        # collect min max values for tiles processed
        zxy = list(tile.coord())
        z = zxy[0]
        min_max = self.zoom_levels.get(z)
        if min_max is None:
            self.zoom_levels[z] = [zxy, zxy] # min, max
        else:
            zz, xx, yy = zip(*(min_max+[zxy]))
            self.zoom_levels[z] = [[z, min(xx), min(yy)], [z, max(xx), max(yy)]]

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


    def __init__(self, root, options=None):
        super(TileDir, self).__init__(root, options)

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
        self.tile_ext = self.dest_ext(tile)
        dest_path = os.path.join(self.root, self.coord2path(*tile.coord())) + self.tile_ext
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
    def __del__(self):
        if self.options.write:
            self.store_metadata()

    def store_metadata(self):
        prm = self.init_pyramid()
        prm.write_tilemap()
        prm.write_metadata()

    def init_pyramid(self):
        log('self.zoom_levels', self.zoom_levels)
        prm = Pyramid.profile_class(self.format)(
            dest = self.root,
            options = dict(
                name=os.path.split(self.root)[1],
                tile_format=self.tile_ext[1:]
                )
            )
        # compute "effective" covered area
        prev_sq = 0
        for z in reversed(sorted(self.zoom_levels)):
            ul_zxy, lr_zxy = self.zoom_levels[z]
            ul_c = prm.tile_bounds(ul_zxy)[0]
            lr_c = prm.tile_bounds(lr_zxy)[1]
            sq = (lr_c[0]-ul_c[0])*(ul_c[1]-lr_c[1])
            area_diff = round(prev_sq/sq, 5)
            log('ul_c, lr_c', z, ul_c, lr_c, sq, area_diff)
            if area_diff == 0.25:
                break # this must be an exact zoom of a previous level
            area_coords = [ul_c, lr_c]
            prev_sq = sq

        prm.set_region(area_coords)
        prm.set_zoom_range(','.join(map(str, self.zoom_levels.keys())))
        return prm

# TileMapDir

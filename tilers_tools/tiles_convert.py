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

import sys
import os
import os.path
import glob
import shutil
import logging
import optparse
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
                self.pyramid = Pyramid.profile_class('zxy')()
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

#############################

class TMStiles(TileMapDir): # see TileMap Diagram at http://wiki.osgeo.org/wiki/Tile_Map_Service_Specification
    'TMS tiles'
#############################
    format, ext, input, output = 'tms', '.tms', True, True
    dir_pattern = '[0-9]*/*/*.*'

    def path2coord(self, tile_path):
        z, x, y = map(int, path2list(tile_path)[-4:-1])
        return (z, x, 2**z-y-1)

    def coord2path(self, z, x, y):
        return '%d/%d/%d' % (z, x, 2**z-y-1)

tile_formats.append(TMStiles)

#############################

class ZXYtiles(TileMapDir): # http://code.google.com/apis/maps/documentation/javascript/v2/overlays.html#Google_Maps_Coordinates
    'Popular ZXY aka XYZ format (Google Maps, OSM, mappero-compatible)'
#############################
    format, ext, input, output = 'zxy', '.zxy', True, True
    dir_pattern = '[0-9]*/*/*.*'

    def path2coord(self, tile_path):
        return map(int, path2list(tile_path)[-4:-1])

    def coord2path(self, z, x, y):
        return '%d/%d/%d' % (z, x, y)

tile_formats.append(ZXYtiles)

#############################

class MapNav(TileDir): # http://mapnav.spb.ru/site/e107_plugins/forum/forum_viewtopic.php?29047.post
    'MapNav (Global Mapper - compatible)'
#############################
    format, ext, input, output = 'mapnav', '.mapnav', True, True
    dir_pattern = 'Z[0-9]*/*/*.pic'
    tile_class = FileTileNoExt

    def dest_ext(self, tile):
        return '.pic'

    def path2coord(self, tile_path):
        z, y, x = path2list(tile_path)[-4:-1]
        return map(int, (z[1:], x, y))

    def coord2path(self, z, x, y):
        return 'Z%d/%d/%d' % (z, y, x)

tile_formats.append(MapNav)

#############################

class SASPlanet(TileDir): # http://sasgis.ru/forum/viewtopic.php?f=2&t=24
    'SASPlanet cache'
#############################
    format, ext, input, output = 'sasplanet', '.sasplanet', True, True
    dir_pattern = 'z[0-9]*/*/x[0-9]*/*/y[0-9]*.*'

    def path2coord(self, tile_path):
        z, dx, x, dy, y = path2list(tile_path)[-6:-1]
        z, x, y = map(int, (z[1:], x[1:], y[1:]))
        return (z-1, x, y)

    def coord2path(self, z, x, y):
        return 'z%d/%d/x%d/%d/y%d' % (z+1, x//1024, x, y//1024, y)

tile_formats.append(SASPlanet)

#############################

class SASGoogle(TileDir):
    'SASPlanet google maps cache'
#############################
    format, ext, input, output = 'sasgoogle', '.sasgoogle', True, True
    dir_pattern = 'z[0-9]*/*/*.*'

    def path2coord(self, tile_path):
        z, y, x = path2list(tile_path)[-4:-1]
        return map(int, (z[1:], x, y))

    def coord2path(self, z, x, y):
        return 'z%d/%d/%d' % (z, x, y)

tile_formats.append(SASGoogle)

#############################

class MapperSQLite(TileSet):
    'maemo-mapper SQLite cache'
#############################
    format, ext, input, output = 'mapper', '.db', True, True
    max_zoom = 20

    def __init__(self, root, options=None):
        super(MapperSQLite, self).__init__(root, options)

        import sqlite3

        self.db = sqlite3.connect(self.root)
        self.dbc = self.db.cursor()
        if self.options.write:
            try:
                self.dbc.execute (
                    'CREATE TABLE maps ('
                        'zoom INTEGER, '
                        'tilex INTEGER, '
                        'tiley INTEGER, '
                        'pixbuf BLOB, '
                        'PRIMARY KEY (zoom, tilex, tiley));'
                    )
            except:
                pass

    def __del__(self):
        self.db.commit()
        self.db.close()
        TileSet.__del__(self)

    def __iter__(self):
        self.dbc.execute('SELECT * FROM maps')
        for z, x, y, pixbuf in self.dbc:
            coord = self.max_zoom+1-z, x, y
            if not self.in_range(coord):
                continue
            self.counter()
            yield PixBufTile(coord, str(pixbuf), (z, x, y))

    def store_tile(self, tile):
        z, x, y = tile.coord()
        # convert to maemo-mapper coords
        z = self.max_zoom+1-z
        log('%s -> SQLite %d, %d, %d' % (tile.path, z, x, y))
        self.dbc.execute('INSERT OR REPLACE INTO maps (zoom, tilex, tiley, pixbuf) VALUES (?, ?, ?, ?);',
            (z, x, y, buffer(tile.data())))
        self.counter()

tile_formats.append(MapperSQLite)

# MapperSQLite

#############################

class MapperGDBM(TileSet): # due to GDBM weirdness on ARM this only works if run on the tablet itself
    'maemo-mapper GDBM cache (works only on Nokia tablet)'
#############################
    format, ext, input, output = 'gdbm', '.gdbm', True, True
    max_zoom = 20

    def __init__(self, root, options=None):

        super(MapperGDBM, self).__init__(root, options)
        #print self.root

        import platform
        assert platform.machine().startswith('arm'), 'This convertion works only on a Nokia tablet'

        import gdbm
        self.db = gdbm.open(self.root, 'cf' if write else 'r')

        self.key = Struct('>III')

    def __del__(self):
        self.db.sync()
        self.db.close()
        TileSet.__del__(self)

    def __iter__(self):
        key = self.db.firstkey()
        while key:
            z, x, y = self.key.unpack(key)
            coord = self.max_zoom+1-z, x, y
            if not self.in_range(coord):
                continue
            self.counter()
            yield PixBufTile(coord, self.db[key], (z, x, y))
            key = self.db.nextkey(key)

    def store_tile(self, tile):
        z, x, y = tile.coord()
        # convert to maemo-mapper coords
        z = self.max_zoom+1-z
        log('%s -> GDBM %d, %d, %d' % (tile.path, z, x, y))
        key = self.key.pack(z, x, y)
        self.db[key] = tile.data()
        self.counter()

tile_formats.append(MapperGDBM)

# MapperGDBM

#############################

class SASBerkeley(TileDir):
    'SASPlanet Berkeley DB'
#############################
    format, ext, input, output = 'sdb', '.sdb', True, False
    dir_pattern = 'z[0-9]*/[0-9]*/[0-9]*/*.sdb'

    def __init__(self, root, options=None):
        super(SASBerkeley, self).__init__(root, options)

        from bsddb3 import db
        self.db = db

        #  0  4b 4s // magic
        #  4  4b I, // crc32
        #  8  4b I, // tile size
        #  12 8b d, // tile date
        #  20 2c string // tile version
        #     2c string // tile content-type
        #     BLOB // tile data
        self.header = Struct('<3sBIId')
        self.key = Struct('>Q') # 64 bit, swap bytes

    def __iter__(self):
        log('__iter__', os.path.join(self.root, self.dir_pattern), glob.iglob(os.path.join(self.root, self.dir_pattern)))
        for db_file in glob.iglob(os.path.join(self.root, self.dir_pattern)):
            log('db_file', db_file)
            for coord, tile, path in self.iter_tiles(db_file):
                self.counter()
                yield PixBufTile(coord, tile, path)

    def iter_tiles(self, db_path):
        zoom = self.get_zoom(db_path) # also checks data in range
        if not zoom:
            return
        d = self.db.DB()
        d.open(db_path, '', self.db.DB_BTREE, self.db.DB_RDONLY)
        c = d.cursor()
        item = c.first(dlen=0, doff=0)
        while item:
            key = item[0]
            coord = self.get_coord(zoom, key)
            if self.in_range(coord):
                data = c.current()[1]
                tile = self.get_tile(data)
                if tile:
                    yield coord, tile, [db_path, key]
            item = c.next(dlen=0, doff=0)
        d.close()

    def get_zoom(self, db_path): # u_TileFileNameBerkeleyDB
        z, x10, y10, xy8 = path2list(db_path)[-5:-1]
        zoom = int(z[1:]) - 1
        x_min, y_min = [int(d) << 8 for d in xy8.split('.')]
        x_max, y_max = [d | 0xFF for d in x_min, y_min]

        if self.in_range((zoom, x_min, y_min), (zoom, x_max, y_max)):
            return zoom
        else:
            return None

    def get_coord(self, zoom, key): # u_BerkeleyDBKey.pas TBerkeleyDBKey.PointToKey
        if key == '\xff\xff\xff\xff\xff\xff\xff\xff':
            return None
        kxy = self.key.unpack(key)[0] # swaps bytes
        xy = [0, 0]
        for bit_n in range(64): # bits for x and y are interleaved in the key
            x0y1 = bit_n % 2 # x, y
            xy[x0y1] += (kxy >> bit_n & 1) << (bit_n - x0y1) / 2

        coord = [zoom] + xy
        log('get_coord', coord, zoom, key, hex(kxy), hex(xy[0]), hex(xy[1]))
        return coord

    def get_tile(self, data): # u_BerkeleyDBValue

        magic, magic_v, crc32, tile_size, tile_date = self.header.unpack_from(data)
        if magic != 'TLD' or magic_v != 3:
            log('get_tile', 'wrong magic', magic, magic_v)
            return None

        strings = []
        i = start = self.header.size
        while True:
            if data[i] == '\x00':
                strings.append(data[start: i: 2])
                start = i + 2
                if len(strings) == 2:
                    break
            i += 2

        tile_version, content_type = strings
        tile_data = data[start: start + tile_size]

        log('get_tile', self.header.size, magic, magic_v, tile_version, content_type, tile_size, tile_data[:20])#, data[4:60:2])
        return tile_data

tile_formats.append(SASBerkeley)

# SASBerkeley

#############################

class Mmap(TileSet):
    'mmap'
#############################
    format, ext, input, output = 'mmap', '.sqlite', True, True
    max_zoom = 20

    def __init__(self, root, options=None):

        path, name = os.path.split(root)
        #~ root = (self.format if path in ['', '.'] else path) + self.ext

        super(Mmap, self).__init__(root, options)
        self.name = self.options.name or os.path.splitext(name)[0]

        import sqlite3
        import base64
        self.b64encode = base64.b64encode
        self.b64decode = base64.b64decode

        self.db = sqlite3.connect(self.root)
        self.dbc = self.db.cursor()
        if self.options.write:
            try:
                self.dbc.execute ('PRAGMA auto_vacuum = INCREMENTAL;')
            finally:
                pass
            self.dbc.execute (
                'CREATE TABLE IF NOT EXISTS __WebKitDatabaseInfoTable__ ('
                    'key TEXT NOT NULL ON CONFLICT FAIL UNIQUE ON CONFLICT REPLACE,'
                    'value TEXT NOT NULL ON CONFLICT FAIL'
                    ');'
                )
            self.dbc.execute (
                "INSERT OR REPLACE INTO __WebKitDatabaseInfoTable__ VALUES('WebKitDatabaseVersionKey','');"
                )
            statements = [
            ]
            for table in ['layers', 'tiles']:
                self.dbc.execute (
                    'CREATE TABLE IF NOT EXISTS "%(table)s" ('
                        'id INTEGER PRIMARY KEY,'
                        'rank INTEGER,'
                        'xmin INTEGER,'
                        'xmax INTEGER,'
                        'ymin INTEGER,'
                        'ymax INTEGER,'
                        '"group" INTEGER,'
                        'geometry TEXT,'
                        'properties TEXT,'
                        'data TEXT'
                    ');'
                    % {'table': table}
                    )
            self.db.commit()

            properties = json.dumps(dict(
                name = self.name,
                url = self.options.url,
                isBaseLayer = not self.options.overlay,
                description = self.options.description
            ))

            self.dbc.execute(
                'INSERT OR REPLACE INTO layers (properties) VALUES (?);',
                (properties, )
                )
            self.layer_id = self.dbc.lastrowid
            log('self.layer', self.layer_id)

    def __del__(self):
        if self.options.write:
            for table in ['layers', 'tiles']:
                self.dbc.execute (
                    'CREATE INDEX IF NOT EXISTS "%(table)s_rank_bbox" ON "%(table)s"'
                        '(rank, xmin, xmax, ymin, ymax);'
                    % {'table': table}
                    )
        self.db.commit()
        self.db.close()
        super(Mmap, self).__del__()

    def __iter__(self):
        yield None
        #~ self.dbc.execute('SELECT * FROM tiles')
        #~ for z, x, y, data in self.dbc:
            #~ coord = z, x, y
            #~ if not self.in_range(coord):
                #~ continue
            #~ self.counter()
            #~ yield PixBufTile(coord, self.b64decode(data), (self.root, coord))

    def store_tile(self, tile):
        z, x, y = tile.coord()
        log('%s -> %s:%d, %d, %d' % (tile.path, self.name, z, x, y))

        self.dbc.execute(
            'INSERT INTO tiles'
                '(rank, xmin, xmax, ymin, ymax, "group", geometry, properties, data)'
                'VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?);',
            (z, x, x, y, y, self.layer_id, None, None,
                'data:' + tile.get_mime() + ';base64,' + self.b64encode(tile.data())))

        self.counter()

tile_formats.append(Mmap)

# Mmap

#----------------------------

def list_formats():

#----------------------------
    for cl in tile_formats:
        print '%10s\t%s%s\t%s' % (
            cl.format,
            'r' if cl.input else ' ',
            'w' if cl.output else ' ',
            cl.__doc__
            )

#----------------------------

def tiles_convert(src_lst, options):

#----------------------------
    for in_class in tile_formats:
        if in_class.input and options.in_fmt == in_class.format:
            break
    else:
        raise Exception('Invalid input format: %s' % options.in_fmt)
    for out_class in tile_formats:
        if out_class.output and options.out_fmt == out_class.format:
            break
    else:
        raise Exception('Invalid output format: %s' % options.out_fmt)

    for src in src_lst:
        dest = os.path.join(options.dst_dir, os.path.splitext(os.path.split(src)[1])[0]+out_class.ext)
        pf('%s -> %s ' % (src, dest), end='')
        out_class(dest, LooseDict(options, write=True)).load_from(in_class(src, options))
        pf('')

#----------------------------

def main(argv):

#----------------------------
    parser = optparse.OptionParser(
        usage='usage: %prog [<options>...] <source>...',
        version=version,
        description='copies map tiles from one structure to another')
    parser.add_option('--from', dest='in_fmt', default='zxy',
        help='input tiles format (default: zxy)')
    parser.add_option('--to', dest='out_fmt', default='mmap',
        help='output tiles format (default: mmap)')
    parser.add_option('--formats', action='store_true', dest='list_formats',
        help='list available formats')
    parser.add_option('-a', '--append', action='store_true', dest='append',
        help='append tiles to an existing destination')
    parser.add_option('-r', '--remove-dest', action='store_true',dest='remove_dest',
        help='delete destination directory before merging')
    parser.add_option('-t', '--dest-dir', default='.', dest='dst_dir',
        help='destination directory (default: current)')
    parser.add_option('--name', default=None,
        help='map name (default: derived from the source)')
    parser.add_option('--description', default='',
        help='map decription (default: None)')
    parser.add_option('--overlay', action='store_true',
        help='non-base layer (default: False)')
    parser.add_option('--url', default=None,
        help='URL template (default: None)')
    parser.add_option('-l', '--link', action='store_true', dest='link',
        help='make links to source tiles instead of copying if possible')
    parser.add_option('--region', default=None, metavar='DATASOURCE',
        help='region to process (OGR shape)')
    parser.add_option('-z', '--zoom', default=None,metavar='ZOOM_LIST',
        help='list of zoom ranges to process')
    parser.add_option('-d', '--debug', action='store_true', dest='debug')
    parser.add_option('-q', '--quiet', action='store_true', dest='quiet')

    global options
    (options, args) = parser.parse_args(argv[1:])

    logging.basicConfig(level=logging.DEBUG if options.debug else
        (logging.ERROR if options.quiet else logging.INFO))
    log(options.__dict__)

    if options.list_formats:
        list_formats()
        sys.exit(0)

    src_lst=args

    tiles_convert(src_lst, LooseDict(options))

# main()

if __name__ == '__main__':

    main(sys.argv)

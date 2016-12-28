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

import binascii
import time
import os.path

from converter_backend import *

#############################

class SASPlanet(TileDir): # http://www.sasgis.org/forum/viewtopic.php?f=2&t=24
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

tileset_profiles.append(SASPlanet)

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
        self.header = struct.Struct('<3sBIId')
        self.key = struct.Struct('>Q') # 64 bit, swap bytes

    def __iter__(self):
        log('__iter__', os.path.join(self.root, self.dir_pattern), glob.iglob(os.path.join(self.root, self.dir_pattern)))
        for db_file in glob.iglob(os.path.join(self.root, self.dir_pattern)):
            log('db_file', db_file)
            for coord, tile, path in self.iter_tiles(db_file):
                #~ log('db tile', coord, tile[:20], path)
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
                tile = self.get_image(data)
                if tile:
                    log('tile', coord)
                    yield coord, tile, [db_path, key]
            item = c.next(dlen=0, doff=0)
        d.close()

    def get_zoom(self, db_path): # u_TileFileNameBerkeleyDB
        z, x10, y10, xy8 = path2list(db_path)[-5:-1]
        zoom = int(z[1:]) - 1
        x_min, y_min = [int(d) << 8 for d in xy8.split('.')]
        x_max, y_max = [d | 0xFF for d in x_min, y_min]

        if not self.in_range((zoom, x_min, y_min), (zoom, x_max, y_max)):
            return None
            pass
        log('get_zoom', zoom, x_min, x_max, y_min, y_max,db_path)
        return zoom

    def get_coord(self, zoom, key): # u_BerkeleyDBKey.pas TBerkeleyDBKey.PointToKey
        if key == '\xff\xff\xff\xff\xff\xff\xff\xff':
            return None
        kxy = self.key.unpack(key)[0] # swaps bytes
        xy = [0, 0]
        for bit_n in range(64): # bits for x and y are interleaved in the key
            x0y1 = bit_n % 2 # x, y
            xy[x0y1] += (kxy >> bit_n & 1) << (bit_n - x0y1) / 2

        coord = [zoom] + xy
        #~ log('get_coord', coord, zoom, key, hex(kxy), hex(xy[0]), hex(xy[1]))
        return coord

    def get_image(self, data): # u_BerkeleyDBValue

        magic, magic_v, crc32, tile_size, tile_date = self.header.unpack_from(data)
        if magic != 'TLD' or magic_v != 3:
            log('get_image', 'wrong magic', magic, magic_v)
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

        #~ log('get_image', self.header.size, magic, magic_v, tile_version, content_type, tile_size, tile_data[:20])#, data[4:60:2])
        return tile_data

tileset_profiles.append(SASBerkeley)

# SASBerkeley

#############################

class SASSQLite(TileDir):
    'SASPlanet SQLite DB'
#############################
    format, ext, input, output = 'slite', '.sqlitedb', True, True
    dir_pattern = 'z[0-9]*/[0-9]*/[0-9]*/*.sqlitedb'

    # from u_TileStorageSQLiteHolder.pas
    create_sql = [
        'PRAGMA synchronous = OFF;',
        'PRAGMA page_size = 16384;',
        'CREATE TABLE IF NOT EXISTS t ('+
            'x INTEGER NOT NULL,'+
            'y INTEGER NOT NULL,'+
            'v INTEGER DEFAULT 0 NOT NULL,'+ # version
            'c TEXT,'+                       # content_type
            's INTEGER DEFAULT 0 NOT NULL,'+ # size
            'h INTEGER DEFAULT 0 NOT NULL,'+ # crc32
            'd INTEGER NOT NULL,'+           # date as unix seconds DEFAULT (strftime(''%s'',''now'')))
            'b BLOB,'+                       # body
            'constraint PK_TB primary key (x,y,v));',
        'CREATE INDEX IF NOT EXISTS t_v_idx on t (v);'
    ]

    db_path = None
    db = None

    def __init__(self, *args, **kw_args):
        super(SASSQLite, self).__init__(*args, **kw_args)

        import sqlite3
        self.sqlite3 = sqlite3

    def __iter__(self):
        log('__iter__', os.path.join(self.root, self.dir_pattern), glob.iglob(os.path.join(self.root, self.dir_pattern)))
        for db_file in glob.iglob(os.path.join(self.root, self.dir_pattern)):
            log('db_file', db_file)
            for tile in self.iter_db_file(db_file):
                yield tile

    def iter_db_file(self, db_path):
        zoom = self.get_zoom(db_path) # also checks data in range
        if not zoom:
            return

        db = self.sqlite3.connect(db_path)
        dbc = db.cursor()

        dbc.execute('SELECT x, y, MAX(v), b FROM t GROUP BY x,y;')
        for x, y, version, data in dbc:
            if data:
                coord = [zoom, x, y]
                #~ log('db tile', coord, tile[:20], path)
                #~ log('tile', coord, data)
                yield PixBufTile(coord, data, key=(db_path, coord))

        db.close()

    def get_zoom(self, db_path): # u_TileFileNameSQLite.pas
        z, x10, y10, xy8 = path2list(db_path)[-5:-1]
        zoom = int(z[1:]) - 1
        x_min, y_min = [int(d) << 8 for d in xy8.split('.')]
        x_max, y_max = [d | 0xFF for d in x_min, y_min]

        if not self.in_range((zoom, x_min, y_min), (zoom, x_max, y_max)):
            return None
            pass
        log('get_zoom', zoom, x_min, x_max, y_min, y_max,db_path)
        return zoom

    def store_tile(self, tile):
        z, x, y = tile.coord()
        data = buffer(tile.data())
        log('%s -> %s:%d, %d, %d' % (tile.path, self.name, z, x, y))

        try:
            data_type = tile.get_mime()
        except KeyError:
            return

        self.set_db(tile)

        timestamp = int(time.time())
        self.dbc.execute(
            'INSERT OR REPLACE INTO t '
                '(x, y, s, h, d, b) '
                'VALUES (?, ?, ?, ?, ?, ?);',
            (x, y, len(data), binascii.crc32(data) % (1 << 32), timestamp, data)
        )

    def set_db(self, tile):
        z, x, y = tile.coord()
        db_dir = os.path.join(self.root, 'z' + str(z + 1), str(x >> 10), str(y >> 10))
        db_name = '%d.%d%s' % (x >> 8, y >> 8, self.ext)

        db_path = os.path.join(db_dir, db_name)
        log(db_path)
        if self.db_path == db_path:
            return

        if self.db:
            self.db.commit()
            self.db.close()

        if not os.path.exists(db_dir):
            os.makedirs(db_dir)

        self.db_path = db_path
        self.db = self.sqlite3.connect(db_path)
        self.dbc = self.db.cursor()

        for stmt in self.create_sql:
            self.dbc.execute(stmt)
        self.db.commit()

    def finalize_tileset(self):
        if self.db:
            self.db.commit()
            self.db.close()
            self.db = None

tileset_profiles.append(SASSQLite)

# SASSQLite

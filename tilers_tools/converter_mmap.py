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

from converter_backend import *

#############################

class Mmap(TileSet):
    'mmap'
#############################
    format, ext, input, output = 'mmap', '.sqlite', True, True
    max_zoom = 20

    def __init__(self, *args, **kw_args):

        super(Mmap, self).__init__(*args, **kw_args)

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
                projection = self.srs,
                description = self.options.description
            ))

            self.dbc.execute(
                'INSERT OR REPLACE INTO layers (properties) VALUES (?);',
                (properties, )
                )
            self.layer_id = self.dbc.lastrowid
            log('self.layer', self.layer_id)

    def __del__(self):
        self.db.commit()
        if self.options.write:
            for table in ['layers', 'tiles']:
                stmt = (
                    'CREATE INDEX IF NOT EXISTS "%(table)s_rank_bbox" ON "%(table)s" '
                        '(rank, xmin, xmax, ymin, ymax);'
                    % {'table': table}
                    )
                self.dbc.execute (stmt)

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

        try:
            mime_type = tile.get_mime()
        except KeyError:
            return
        b64_data = self.b64encode(tile.data())
        data_url = 'data:' + mime_type + ';base64,' + b64_data
        self.dbc.execute(
            'INSERT INTO tiles'
                '(rank, xmin, xmax, ymin, ymax, "group", geometry, properties, data)'
                'VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?);',
            (z, x, x, y, y, self.layer_id, None, None, data_url)
        )

        self.counter()

tile_formats.append(Mmap)

# Mmap

#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
# Copyright (c) 2010,2011 Vadim Shlyakhov
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

from tiler_functions import *
from gdal_tiler import Pyramid

#############################

class Tile(object):

#############################
    def __init__(self,coord):
        self._coord=coord

    def coord(self):
        return self._coord

    def get_mime(self):
        return mime_from_ext(self.get_ext())

#############################

class FileTile(Tile):

#############################
    def __init__(self,coord,path):
        super(FileTile, self).__init__(coord)
        self.path=path

    def data(self):
        return open(self.path,'rb').read()

    def get_ext(self):
        return os.path.splitext(self.path)[1]

    def copy2file(self,dst,link=False):
        if link and os.name == 'posix':
            dst_dir=os.path.split(dst)[0]
            src=os.path.relpath(self.path,dst_dir)
            os.symlink(src,dst)
        else:
            shutil.copy(self.path,dst)

#############################

class FileTileNoExt(FileTile):

#############################
    def get_ext(self):
        return ext_from_file(self.path)

#############################

class PixBufTile(Tile):

#############################
    def __init__(self,coord,pixbuf,key):
        super(PixBufTile, self).__init__(coord)
        self.pixbuf=pixbuf
        self.path=repr(key) # only for debugging

    def data(self):
        return self.pixbuf

    def get_ext(self):
        ext=ext_from_buffer(self.pixbuf)
        return ext

    def copy2file(self,dest_path,link=False):
        open(dest_path,'wb').write(self.pixbuf)

#############################

class TileSet(object):

#############################
    def __init__(self,root,options=None):
        options=options or LooseOptions()
        log('root',root,options)
        self.root=root
        self.options=options
        self.write=False
        self.write=options.write
        self.zoom_levels={}

        if not self.write:
            assert os.path.exists(root), "No file or directory found: %s" % root
        else:
            if os.path.exists(self.root):
                if self.options.remove_dest:
                    if os.path.isdir(self.root):
                        shutil.rmtree(self.root,ignore_errors=True)
                    else:
                        os.remove(self.root)
                else:
                    assert self.options.append, "Destination already exists: %s" % root
            if self.options.region:
                prm=Pyramid.profile_class('zxy')()
                prm.set_zoom_range(self.options.zoom)
                prm.load_region(self.options.region)
                self.my_tile=lambda tile: prm.belongs_to(tile.coord())

    def my_tile(self, tile):
        return True

    def __del__(self):
        log('self.count',self.count)

    def __iter__(self): # to be defined by a child
        raise Exception("Unimplemented!")

    def load_from(self,src_tiles):
        log('load_from',(src_tiles.root, self.root))
        map(self.process_tile,src_tiles)

    def process_tile(self, tile):
        #log("process_tile",tile)
        if not self.my_tile(tile):
            return
        self.store_tile(tile)

        # collect min max values for tiles processed
        zxy=list(tile.coord())
        z=zxy[0]
        min_max=self.zoom_levels.get(z)
        if min_max is None:
            self.zoom_levels[z]=[zxy,zxy] # min,max
        else:
            zz,xx,yy=zip(*(min_max+[zxy]))
            self.zoom_levels[z]=[[z,min(xx),min(yy)],[z,max(xx),max(yy)]]

    count=0
    tick_rate=100
    def counter(self):
        self.count+=1
        if self.count % self.tick_rate == 0:
            pf('.',end='')
            return True
        else:
            return False

# TileSet

#############################

class TileDir(TileSet):

#############################
    tile_class = FileTile


    def __init__(self,root,options=None,write=False):
        super(TileDir, self).__init__(root,options)

        if self.write:
            try:
                os.makedirs(self.root)
            except os.error: pass

    def __iter__(self):
        for f in glob.iglob(os.path.join(self.root,self.dir_pattern)):
            self.counter()
            yield self.tile_class(self.path2coord(f),f)

    def path2coord(self,tile_path):
        raise Exception("Unimplemented!")

    def coord2path(self,z,x,y):
        raise Exception("Unimplemented!")

    def dest_ext(self, tile):
        return tile.get_ext()

    def store_tile(self, tile):
        self.tile_ext=self.dest_ext(tile)
        dest_path=os.path.join(self.root,self.coord2path(*tile.coord())) + self.tile_ext
        log('%s -> %s' % (tile.path,dest_path))
        try:
            os.makedirs(os.path.split(dest_path)[0])
        except os.error: pass
        tile.copy2file(dest_path,self.options.link)
        self.counter()
# TileDir

#############################

class TileMapDir(TileDir):

#############################
    def __del__(self):
        if self.write:
            self.store_metadata()

    def store_metadata(self):
        prm=self.init_pyramid()
        prm.write_tilemap()
        prm.write_html()

    def init_pyramid(self):
        log('self.zoom_levels',self.zoom_levels)
        prm=Pyramid.profile_class(self.format)(
            dest=self.root,
            options=dict(
                name=os.path.split(self.root)[1],
                tile_format=self.tile_ext[1:]
                )
            )
        # compute "effective" covered area
        prev_sq=0
        for z in reversed(sorted(self.zoom_levels)):
            ul_zxy,lr_zxy=self.zoom_levels[z]
            ul_c=prm.tile_bounds(ul_zxy)[0]
            lr_c=prm.tile_bounds(lr_zxy)[1]
            sq=(lr_c[0]-ul_c[0])*(ul_c[1]-lr_c[1])
            area_diff=round(prev_sq/sq,5)
            log('ul_c,lr_c',z,ul_c,lr_c,sq,area_diff)
            if area_diff == 0.25:
                break # this must be an exact zoom of a previous level
            area_coords=[ul_c,lr_c]
            prev_sq=sq

        prm.set_region(area_coords)
        prm.set_zoom_range(','.join(map(str,self.zoom_levels.keys())))
        return prm

#############################

class TMStiles(TileMapDir): # see TileMap Diagram at http://wiki.osgeo.org/wiki/Tile_Map_Service_Specification
    'TMS tiles'
#############################
    format,ext,input,output='tms','.tms',True,True
    dir_pattern='[0-9]*/*/*.*'

    def path2coord(self,tile_path):
        z,x,y=map(int,path2list(tile_path)[-4:-1])
        return (z,x,2**z-y-1)

    def coord2path(self,z,x,y):
        return '%d/%d/%d' % (z,x,2**z-y-1)

#############################

class ZXYtiles(TileMapDir): # http://code.google.com/apis/maps/documentation/javascript/v2/overlays.html#Google_Maps_Coordinates
    'Popular ZXY aka XYZ format (Google Maps, OSM, mappero-compatible)'
#############################
    format,ext,input,output='zxy','.zxy',True,True
    dir_pattern='[0-9]*/*/*.*'

    def path2coord(self,tile_path):
        return map(int,path2list(tile_path)[-4:-1])

    def coord2path(self,z,x,y):
        return '%d/%d/%d' % (z,x,y)

#############################

class MapNav(TileDir): # http://mapnav.spb.ru/site/e107_plugins/forum/forum_viewtopic.php?29047.post
    'MapNav (Global Mapper - compatible)'
#############################
    format,ext,input,output='mapnav','.mapnav',True,True
    dir_pattern='Z[0-9]*/*/*.pic'
    tile_class = FileTileNoExt

    def dest_ext(self, tile):
        return '.pic'

    def path2coord(self,tile_path):
        z,y,x=path2list(tile_path)[-4:-1]
        return map(int,(z[1:],x,y))

    def coord2path(self,z,x,y):
        return 'Z%d/%d/%d' % (z,y,x)

#############################

class SASPlanet(TileDir): # http://sasgis.ru/forum/viewtopic.php?f=2&t=24
    'SASPlanet cache'
#############################
    format,ext,input,output='sasplanet','.sasplanet',True,True
    dir_pattern='z[0-9]*/*/x[0-9]*/*/y[0-9]*.*'

    def path2coord(self,tile_path):
        z,dx,x,dy,y=path2list(tile_path)[-6:-1]
        z,x,y=map(int,(z[1:],x[1:],y[1:]))
        return (z-1,x,y)

    def coord2path(self,z,x,y):
        return 'z%d/%d/x%d/%d/y%d' % (z+1, x//1024, x, y//1024, y)

#############################

class SASGoogle(TileDir):
    'SASPlanet google maps cache'
#############################
    format,ext,input,output='sasgoogle','.sasgoogle',True,True
    dir_pattern='z[0-9]*/*/*.*'

    def path2coord(self,tile_path):
        z,y,x=path2list(tile_path)[-4:-1]
        return map(int,(z[1:],x,y))

    def coord2path(self,z,x,y):
        return 'z%d/%d/%d' % (z,x,y)

#############################

class WebSQL(TileSet):
    'WebSQL'
#############################
    format,ext,input,output='websql','.sqlite',True,True
    max_zoom=20

    def __init__(self,root,options=None,write=False):

        path,name=os.path.split(root)
        root= (self.format if path in ['','.'] else path)+self.ext

        super(WebSQL, self).__init__(root,options)
        self.name= self.options.name or os.path.splitext(name)[0]

        import sqlite3
        import base64
        self.b64encode=base64.b64encode
        self.b64decode=base64.b64decode

        self.db=sqlite3.connect(self.root)
        self.dbc = self.db.cursor()
        if self.write:
            self.dbc.execute (
                "CREATE TABLE IF NOT EXISTS __WebKitDatabaseInfoTable__ ("
                    "key TEXT NOT NULL ON CONFLICT FAIL UNIQUE ON CONFLICT REPLACE,"
                    "value TEXT NOT NULL ON CONFLICT FAIL"
                    ");"
                )
            self.dbc.execute (
                "INSERT OR REPLACE INTO __WebKitDatabaseInfoTable__ VALUES('WebKitDatabaseVersionKey','');"
                )
            self.dbc.execute (
                "CREATE TABLE IF NOT EXISTS layers ("
                    "layer INTEGER PRIMARY KEY NOT NULL,"
                    "name TEXT UNIQUE NOT NULL,"
                    "url TEXT,"
                    "overlay BOOLEAN,"
                    "description TEXT"
                    ")"
                )
            self.dbc.execute (
                "CREATE TABLE IF NOT EXISTS tiles ("
                    "tile INTEGER PRIMARY KEY,"
                    "layer INTEGER NOT NULL REFERENCES layers,"
                    "z INTEGER,"
                    "x INTEGER,"
                    "y INTEGER,"
                    "data TEXT,"
                    "UNIQUE (layer, z, x, y)"
                    ")"
                )
#                    "raster INTEGER REFERENCES rasters,"
            #~ self.dbc.execute (
                #~ "CREATE TABLE IF NOT EXISTS rasters ("
                    #~ "raster INTEGER PRIMARY KEY,"
                    #~ "tile INTEGER NOT NULL REFERENCES tiles,"
                    #~ "layer INTEGER NOT NULL REFERENCES layers,"
                    #~ "data TEXT"
                    #~ ")"
                #~ )
#                    "mime_type TEXT,"
            self.db.commit()

            self.dbc.execute("INSERT OR IGNORE INTO layers (name) VALUES (?);",
                (self.name,))
            self.dbc.execute("SELECT layer FROM layers WHERE name=?",(self.name,))
            row=self.dbc.fetchone()
            self.layer_id=row[0]
            log("self.layer",self.layer_id)

            self.dbc.execute("UPDATE layers SET url=?,overlay=?,description=? WHERE name=?;",
                (self.options.url,self.options.overlay,self.options.description,self.name))

    def __del__(self):
        self.db.commit()
        self.db.close()
        super(WebSQL, self).__del__()

    def __iter__(self):
        self.dbc.execute("SELECT * FROM tiles")
        for z,x,y,data in self.dbc:
            self.counter()
            yield PixBufTile((z,x,y),self.b64decode(data),(z,x,y))

    def store_tile(self, tile):
        z,x,y=tile.coord()
        log("%s -> WebSQL %s:%d,%d,%d" % (tile.path,self.name, z, x, y))
        #~ self.dbc.execute(
            #~ "SELECT tile,raster "
            #~ "FROM tiles "
            #~ "WHERE layer=? AND z=? AND x=? AND y=?",
            #~ (self.layer_id,z,x,y))
        #~ row=self.dbc.fetchone()
        #~ if row:
            #~ tile_id,raster_id=row
            #~ self.dbc.execute("DELETE FROM rasters WHERE raster=?",(raster_id,))
            #~ self.dbc.execute("DELETE FROM tiles WHERE tile=?",(tile_id,))

        self.dbc.execute("INSERT INTO tiles (layer,z,x,y,data) VALUES (?,?,?,?,?);",
            (self.layer_id,z,x,y,'data:' + tile.get_mime() + ';base64,' + self.b64encode(tile.data())))

#~ #        self.dbc.execute("INSERT INTO tiles (layer,z,x,y) VALUES (?,?,?,?);",
#~ #            (self.layer_id,z,x,y))
        #~ tile_id=self.dbc.lastrowid
        #~ self.dbc.execute(
#~ #            "INSERT INTO rasters (tile,layer,mime_type,data) VALUES (?,?,?,?);",
#~ #            (tile_id,self.layer_id,tile.get_mime(),self.b64encode(tile.data()))
            #~ "INSERT INTO rasters (tile,layer,data) VALUES (?,?,?);",
            #~ (tile_id,self.layer_id,'data:' + tile.get_mime() + ';base64,' + self.b64encode(tile.data()))
            #~ )
        #~ raster_id=self.dbc.lastrowid
        #~ self.dbc.execute("UPDATE tiles SET raster=? WHERE tile=?;",(raster_id,tile_id))

        self.counter()
# tst

#############################

class MapperSQLite(TileSet):
    'maemo-mapper SQLite cache'
#############################
    format,ext,input,output='sqlite','.db',True,True
    max_zoom=20

    def __init__(self,root,options=None,write=False):
        super(MapperSQLite, self).__init__(root,options)

        import sqlite3

        self.db=sqlite3.connect(self.root)
        self.dbc = self.db.cursor()
        if self.write:
            try:
                self.dbc.execute (
                    "CREATE TABLE maps ("
                        "zoom INTEGER,"
                        "tilex INTEGER,"
                        "tiley INTEGER,"
                        "pixbuf BLOB,"
                        "PRIMARY KEY (zoom, tilex, tiley));"
                    )
            except:
                pass

    def __del__(self):
        self.db.commit()
        self.db.close()
        TileSet.__del__(self)

    def __iter__(self):
        self.dbc.execute("SELECT * FROM maps")
        for z,x,y,pixbuf in self.dbc:
            self.counter()
            yield PixBufTile((self.max_zoom+1-z,x,y),str(pixbuf),(z,x,y))

    def store_tile(self, tile):
        z,x,y=tile.coord()
        # convert to maemo-mapper coords
        z=self.max_zoom+1-z
        log('%s -> SQLite %d,%d,%d' % (tile.path, z, x, y))
        self.dbc.execute("INSERT OR REPLACE INTO maps (zoom,tilex,tiley,pixbuf) VALUES (?,?,?,?);",
            (z,x,y,buffer(tile.data())))
        self.counter()
# MapperSQLite

#############################

class MapperGDBM(TileSet): # due to GDBM weirdness on ARM this only works if run on the tablet itself
    'maemo-mapper GDBM cache (works only on Nokia tablet)'
#############################
    format,ext,input,output='gdbm','.gdbm',True,True
    max_zoom=20

    def __init__(self,root,options=None,write=False):

        super(MapperGDBM, self).__init__(root,options)

        import platform
        assert platform.machine().startswith('arm'), 'This convertion works only on a Nokia tablet'
        import gdbm
        import struct

        self.pack=struct.pack
        self.unpack=struct.unpack
        print self.root
        self.db=gdbm.open(self.root, 'cf' if write else 'r')

    def __del__(self):
        self.db.sync()
        self.db.close()
        TileSet.__del__(self)

    def __iter__(self):
        key=self.db.firstkey()
        while key:
            z,x,y=self.unpack('>III',key)
            self.counter()
            yield PixBufTile((self.max_zoom+1-z,x,y),self.db[key],(z,x,y))
            key=self.db.nextkey(key)

    def store_tile(self, tile):
        z,x,y=tile.coord()
        # convert to maemo-mapper coords
        z=self.max_zoom+1-z
        log('%s -> GDBM %d,%d,%d' % (tile.path, z, x, y))
        key=self.pack('>III',z,x,y)
        self.db[key]=tile.data()
        self.counter()
# MapperGDBM

tile_formats=(
    TMStiles,
    MapperSQLite,
    MapperGDBM,
    ZXYtiles,
    MapNav,
    SASPlanet,
    SASGoogle,
    WebSQL
    )

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

def tiles_convert(src_lst,options):

#----------------------------
    for in_class in tile_formats:
        if in_class.input and options.in_fmt == in_class.format:
            break
    else:
        raise Exception("Invalid input format: %s" % options.in_fmt)
    for out_class in tile_formats:
        if out_class.output and options.out_fmt == out_class.format:
            break
    else:
        raise Exception("Invalid output format: %s" % options.out_fmt)

    for src in src_lst:
        dest=os.path.join(options.dst_dir,os.path.splitext(os.path.split(src)[1])[0]+out_class.ext)
        pf('%s -> %s ' % (src,dest),end='')
        out_options=LooseOptions(options)
        out_options.write=True
        out_class(dest,out_options).load_from(in_class(src))
        pf('')

#----------------------------

def main(argv):

#----------------------------
    parser = optparse.OptionParser(
        usage="usage: %prog [<options>...] <source>...",
        version=version,
        description="copies map tiles from one structure to another")
    parser.add_option("--from", dest="in_fmt", default='zxy',
        help='input tiles format (default: zxy)')
    parser.add_option("--to", dest="out_fmt", default='websql',
        help='output tiles format (default: websql)')
    parser.add_option("--formats", action="store_true", dest="list_formats",
        help='list available formats')
    parser.add_option("-a", "--append", action="store_true", dest="append",
        help="append tiles to an existing destination")
    parser.add_option("-r", "--remove-dest", action="store_true",dest="remove_dest",
        help='delete destination directory before merging')
    parser.add_option("-t", "--dest-dir", default='.', dest="dst_dir",
        help='destination directory (default: current)')
    parser.add_option("--name", default=None,
        help='map name (default: derived from the source)')
    parser.add_option("--description", default='',
        help='map decription (default: None)')
    parser.add_option("--overlay", action="store_true",
        help='non-base layer (default: False)')
    parser.add_option("--url", default=None,
        help='URL template (default: None)')
    parser.add_option("-l", "--link", action="store_true", dest="link",
        help='make links to source tiles instead of copying if possible')
    parser.add_option("--region", default=None, metavar="DATASOURCE",
        help='region to process (OGR shape)')
    parser.add_option("-z", "--zoom", default=None,metavar="ZOOM_LIST",
        help='list of zoom ranges to process')
    parser.add_option("-d", "--debug", action="store_true", dest="debug")
    parser.add_option("-q", "--quiet", action="store_true", dest="quiet")

    global options
    (options, args) = parser.parse_args(argv[1:])

    logging.basicConfig(level=logging.DEBUG if options.debug else
        (logging.ERROR if options.quiet else logging.INFO))
    log(options.__dict__)

    if options.list_formats:
        list_formats()
        sys.exit(0)

    src_lst=args

    tiles_convert(src_lst,LooseOptions(options))

# main()

if __name__=='__main__':

    main(sys.argv)


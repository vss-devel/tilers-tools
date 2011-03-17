#!/usr/bin/env python
# -*- coding: utf-8 -*-

# 2011-01-27 11:44:37 

###############################################################################
# Copyright (c) 2010, Vadim Shlyakhov
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

def path2list(path):
    head,ext=os.path.splitext(path)
    split=[ext]
    while head:
        head,p=os.path.split(head)
        split.append(p)
    split.reverse()
    return split

class Tile:
    def __init__(self,coord):
        self._coord=coord
    def coord(self):
        return self._coord

class FileTile(Tile):
    def __init__(self,coord,path):
        Tile.__init__(self,coord)
        self.path=path
    def data(self):
        return open(self.path,'rb').read()
    def copy2file(self,dest_path):
        # dest_path is w/o extension!
        tile_ext=os.path.splitext(self.path)[1]
        shutil.copy(self.path,dest_path+tile_ext)

class PixBufTile(Tile):
    def __init__(self,coord,pixbuf,key):
        Tile.__init__(self,coord)
        self.pixbuf=pixbuf
        self.path=repr(key) # only for debugging
    def data(self):
        return self.pixbuf
    def copy2file(self,dest_path):
        # dest_path is w/o extension!
        tile_ext='.png' if self.pixbuf.startswith('\x89PNG\x0D\x0A\x1A\x0A') else '.jpg'
        open(dest_path+tile_ext,'wb').write(self.pixbuf)

class TileSet:
    def __init__(self,root,write=False,append=False):
        ld(root)
        self.root=root
        self.count=0
        self.tick_rate=100
        if write:
            if not append and os.path.exists(self.root):
                if os.path.isdir(self.root):
                    shutil.rmtree(self.root,ignore_errors=True)
                else:
                    os.remove(self.root)
        else:
            assert os.path.exists(root), "No file or directory found: %s" % root

    def __del__(self):
        ld(self.count)

#    def store_tile(self, tile):
#        raise Exception("Unimplemented!")

    def counter(self):
        self.count+=1
        if self.count % self.tick_rate == 0:
            pf('.',end='')
            return True
        else:
            return False

    def load_from(self,src_tiles):
        ld((src_tiles.root, self.root))
        map(self.store_tile,src_tiles)
# TileSet

class TileDir(TileSet):
    def __init__(self,root,write=False,append=False):
        TileSet.__init__(self,root,write,append)
        if write:
            try:
                os.makedirs(root)
            except os.error: pass
        
    def __iter__(self):
        try:
            wrk_dir=os.getcwd()
            os.chdir(self.root)
            file_lst=glob.glob(self.dir_pattern)
        finally:
            os.chdir(wrk_dir)
        for f in file_lst:
            self.counter()
            yield FileTile(self.path2coord(f),os.path.join(self.root,f))

#    def path2coord(self,tile_path):
#        raise Exception("Unimplemented!")

#    def coord2path(self,z,x,y):
#        raise Exception("Unimplemented!")

    def store_tile(self, tile):
        # dest_path is w/o extension!
        dest_path=os.path.join(self.root,self.coord2path(*tile.coord()))
        ld('%s -> %s' % (tile.path,dest_path))
        try:
            os.makedirs(os.path.split(dest_path)[0])
        except os.error: pass
        tile.copy2file(dest_path)
        self.counter()
# TileDir
                            
class TMStiles(TileDir): # see TileMap Diagram at http://wiki.osgeo.org/wiki/Tile_Map_Service_Specification
    'TMS tiles'
    format,ext,input,output='tms','.tms',True,True
    dir_pattern='[0-9]*/*/*.*'

    def path2coord(self,tile_path):
        z,x,y=map(int,path2list(tile_path)[0:3])
        return (z,x,2**z-1-y)

    def coord2path(self,z,x,y):
        return '%d/%d/%d' % (z,x,2**z-1-y)

class Gmaps(TileDir): # http://code.google.com/apis/maps/documentation/javascript/v2/overlays.html#Google_Maps_Coordinates
    'Google map tiles (mappero-compatible)'
    format,ext,input,output='gmaps','.gmaps',True,True
    dir_pattern='[0-9]*/*/*.*'

    def path2coord(self,tile_path):
        return map(int,path2list(tile_path)[0:3])

    def coord2path(self,z,x,y):
        return '%d/%d/%d' % (z,x,y)

class SASPlanet(TileDir): # http://sasgis.ru/forum/viewtopic.php?f=2&t=24
    'SASPlanet cache'
    format,ext,input,output='sasplanet','.sasplanet',True,True
    dir_pattern='z[0-9]*/*/x[0-9]*/*/y[0-9]*.*'

    def path2coord(self,tile_path):
        z,dx,x,dy,y,ext=path2list(tile_path)
        z,x,y=map(int,(z[1:],x[1:],y[1:]))
        return (z-1,x,y)

    def coord2path(self,z,x,y):
        return 'z%d/%d/x%d/%d/y%d%s' % (z+1,x % 1024,x,y % 1024,y)

class SASGoogle(TileDir):
    'SASPlanet google maps cache'
    format,ext,input,output='sasgoogle','.sasgoogle',True,True
    dir_pattern='z[0-9]*/*/*.*'

    def path2coord(self,tile_path):
        return map(int,path2list(tile_path[1:])[0:3])

    def coord2path(self,z,x,y):
        return 'z%d/%d/%d' % (z,x,y)

class MapperSQLite(TileSet):
    'maemo-mapper SQLite cache'
    format,ext,input,output='sqlite','.db',True,True
    max_zoom=20
    
    def __init__(self,root,write=False,append=False):
        import sqlite3
        TileSet.__init__(self,root,write,append)
        self.db=sqlite3.connect(self.root)
        self.dbc = self.db.cursor()
        if write:
            try:
                self.dbc.execute ('''
                    create table maps (
                        zoom integer,
                        tilex integer,
                        tiley integer,
                        pixbuf blob,
                        primary key (zoom, tilex, tiley));
                    ''')
            except:
                pass

    def __del__(self):
        self.db.commit()
        self.db.close()
        TileSet.__del__(self)

    def __iter__(self):
        self.dbc.execute('select * from maps')
        for z,x,y,pixbuf in self.dbc:
            self.counter()
            yield PixBufTile((self.max_zoom+1-z,x,y),str(pixbuf),(z,x,y))
    
    def store_tile(self, tile):
        z,x,y=tile.coord()
        # convert to maemo-mapper coords
        z=self.max_zoom+1-z
        ld('%s -> SQLite %d,%d,%d' % (tile.path, z, x, y))
        self.dbc.execute('insert or replace into maps (zoom,tilex,tiley,pixbuf) values (?,?,?,?);',
            (z,x,y,buffer(tile.data())))
        self.counter()
# MapperSQLite

class MapperGDBM(TileSet): # due to GDBM weirdness on ARM this only works if run on the tablet itself
    'maemo-mapper GDBM cache (works only on Nokia tablet)'
    format,ext,input,output='gdbm','.gdbm',True,True
    max_zoom=20
    
    def __init__(self,root,write=False,append=False):
        import platform
        assert platform.machine().startswith('arm'), 'This convertion works only on a Nokia tablet'
        import gdbm
        import struct
        TileSet.__init__(self,root,write,append)
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
        ld('%s -> GDBM %d,%d,%d' % (tile.path, z, x, y))
        key=self.pack('>III',z,x,y)
        self.db[key]=tile.data()
        self.counter()
# MapperGDBM

tile_formats=(
    TMStiles,
    MapperSQLite,
    MapperGDBM,
    Gmaps,
    SASPlanet,
    SASGoogle,
    )

def list_formats():
    for cl in tile_formats:
        print '%10s\t%s%s\t%s' % (
            cl.format,
            'r' if cl.input else ' ',
            'w' if cl.output else ' ',
            cl.__doc__
            )

def tiles_convert(in_fmt,src,out_fmt,dest,append):
    for in_class in tile_formats:
        if in_class.input and in_fmt == in_class.format:
            break
    else:
        raise Exception("Invalid input format: %s" % in_fmt)
    for out_class in tile_formats:
        if out_class.output and out_fmt == out_class.format:
            break
    else:
        raise Exception("Invalid output format: %s" % out_fmt)
    dest=dest if dest else src+out_class.ext
    out_class(dest,write=True,append=append).load_from(in_class(src))
      
if __name__=='__main__':
    parser = optparse.OptionParser(
        usage="usage: %prog  <source> [<target>]",
        description="copies map tiles from one structure to another")
    parser.add_option("--from", dest="in_fmt", default='gmaps',
        help='input tiles format (default: gmaps)')
    parser.add_option("--to", dest="out_fmt", default='sqlite',
        help='output tiles format (default: sqlite)')
    parser.add_option("-l", "--list-formats", action="store_true", dest="list_formats",
        help='list available formats')
    parser.add_option("-a", "--append", action="store_true", dest="append",
        help="don't delete destination, append to it")
    parser.add_option("-d", "--debug", action="store_true", dest="debug")
    parser.add_option("-q", "--quiet", action="store_true", dest="quiet")

    (options, args) = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG if options.debug else 
        (logging.ERROR if options.quiet else logging.INFO))
    ld(options)

    if options.list_formats:
        list_formats()
        sys.exit(0)
        
    try:
        src=args[0]
    except:
        raise Exception("No source specified")
    try:
        dest=args[1]
    except:
        dest=None

    tiles_convert(options.in_fmt,src,options.out_fmt,dest,options.append)
    pf('')


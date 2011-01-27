#!/usr/bin/env python

# 2010-10-26 17:29:06 

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
#******************************************************************************

import sys
import os
import shutil
import glob
import logging
import re
import optparse
from PIL import Image

try:
    import multiprocessing # available in python 2.6 and above
except:
    multiprocessing=None

def parallel_map(func,iterable):
    if not multiprocessing:
        res=map(func,iterable)
    else:
        # process files in parallel
        mp_pool = multiprocessing.Pool() # multiprocessing pool
        res=mp_pool.map(func,iterable)
        # wait for threads to finish
        mp_pool.close()
        mp_pool.join()
    return res

def l_d(smth):
    logging.debug(str(smth))
    
def p_f(smth, nl=True):
    #logging.debug(str(smth))
    s=str(smth)
    if nl: s+='\n'
    sys.stdout.write(s)
    sys.stdout.flush()
    
def re_sub_file(fname, substitutions):
    'stream edit file using reg exp substitution list supplied'
    f=open(fname+'.new', 'w')
    for l in open(fname):
        for (pattern,repl) in substitutions:
            l=re.sub(pattern,repl,string=l)
        f.write(l)
    f.close()
    shutil.move(fname+'.new',fname) # mind Windows

class ZoomSet:
    def __init__(self,tiles_root,curr_zoom):
        self.curr_zoom=curr_zoom
        self.zoom=curr_zoom-1
        p_f('%i'%self.zoom,False)
        os.chdir(os.path.join(tiles_root,'%i' % curr_zoom))
        self.src_lst=set([tuple(map(eval,os.path.split(os.path.splitext(i)[0])))
                          for i in glob.glob('*/*.png')])
        l_d(self.src_lst)
        if len(self.src_lst) == 0:
            raise Exception("No tiles in %s" % os.getcwd())
        os.chdir(tiles_root)
        self.dest_lst=set([(src_x/2,src_y/2) for (src_x,src_y) in self.src_lst])
        l_d(self.dest_lst)
        shutil.rmtree('%i' % self.zoom,ignore_errors=True)
        for i in set([x for (x,y) in self.dest_lst]):
            os.makedirs('%i/%i' % (self.zoom,i))
        self.zoom_out()
            
    def __call__(self,dest_xy):
        p_f('.',False)
        (x,y)=dest_xy
        im = Image.new("RGBA",(256,256),(0,0,0,0))
        tiles_map=[(0,128), (128,128),
                    (0,0), (128,0)]
        tiles_in=[(x*2,y*2),(x*2+1,y*2),
                    (x*2,y*2+1),(x*2+1,y*2+1)]
        for (src_xy,out_loc) in zip(tiles_in,tiles_map):
            if src_xy in self.src_lst:
                src_path='%i/%i/%i.png' % (self.curr_zoom,src_xy[0],src_xy[1])
                im.paste(Image.open(src_path).resize((128,128),Image.ANTIALIAS),out_loc)
        dst_path='%i/%i/%i.png' % (self.zoom,x,y)
        im.save(dst_path)

    def zoom_out(self):
        parallel_map(self,self.dest_lst)

# ZoomSet end

if __name__=='__main__':
    usage = "usage: %prog tiles_dir ..."
    parser = optparse.OptionParser(usage=usage, description="")
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose")
    parser.add_option("-z", "--zoom", dest="zoom", type='int', 
        help='target zoom level)')
        
    (options, args) = parser.parse_args()
    if options.zoom == None:
        parser.error('No target zoom specified')

    logging.basicConfig(level=logging.DEBUG if options.verbose else logging.INFO)

    start_dir=os.getcwd()
    for tiles_dir in args if len(args)>0 else ['.']:
        p_f('%s '%tiles_dir,False)    
        os.chdir(tiles_dir)
        tiles_root=os.getcwd()
        min_zoom=options.zoom
        max_zoom=min([int(i) for i in glob.glob('[0-9]*')])

        for zoom in range(max_zoom,min_zoom,-1):
            ZoomSet(tiles_root,zoom)

        os.chdir(tiles_root)
        p_f('')
        # modify googlemaps.html and openlayers.html
        re_sub_file(os.path.join(tiles_root,'googlemaps.html'),
            [('(var mapMinZoom =).*;','\\1 %i;' % min_zoom)])
        re_sub_file(os.path.join(tiles_root,'openlayers.html'),
            [('(var mapMinZoom =).*;','\\1 %i;' % min_zoom)])
        os.chdir(start_dir)

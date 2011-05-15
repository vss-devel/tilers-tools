#!/usr/bin/env python

# 2011-05-11 11:09:08 

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
import stat
import shutil
import logging
import optparse
from PIL import Image

from tiler_functions import *

def find_by_ext(path, ext):
    return flatten([os.path.join(path, name) for name in files if name.endswith(ext)] 
                    for path, dirs, files in os.walk(src_dir))

tick_rate=50
tick_count=0

def counter():
    global tick_count
    tick_count+=1
    if tick_count % tick_rate == 0:
        pf('.',end='')
        return True
    else:
        return False

class KeyboardInterruptError(Exception): pass

def optimize_png(fl):
    'optimize png using pngnq utility'
    counter()
    try:
        command(['pngnq','-fn', options.colors, fl])
    except KeyboardInterrupt: # http://jessenoller.com/2009/01/08/multiprocessingpool-and-keyboardinterrupt/
        pf('got KeyboardInterrupt')
        raise KeyboardInterruptError()

def to_jpeg(fl):
    'convert to jpeg'
    counter()
    try:
        dst_fl=os.path.splitext(fl)[0]+'.jpg'
        img=Image.open(fl)
        img.save(dst_fl,optimize=True,quality=options.quality)
    except KeyboardInterrupt: # http://jessenoller.com/2009/01/08/multiprocessingpool-and-keyboardinterrupt/
        pf('got KeyboardInterrupt')
        raise KeyboardInterruptError()
   
if __name__=='__main__':
    logging.basicConfig(level=logging.INFO)
    #logging.basicConfig(level=logging.DEBUG)

    parser = optparse.OptionParser()
    usage = "usage: %prog [options] arg"
    parser.add_option("-q", "--quiet", action="store_true")
    parser.add_option("-n", "--colors", dest="colors", default='256',
        help='Specifies  the  number  of colors to quantize to (default: 256)')
    parser.add_option("--jpeg", action="store_true",
        help='convert tiles to JPEG')
    parser.add_option("--quality", dest="quality", type="int", default=75,
        help='JPEG quality (default: 75)')
        
    (options, args) = parser.parse_args()

    if not args:
        parser.error('No input directory(s) specified')

    for src_dir in args:
        pf(src_dir,end='')
        # delete rouge *-nq8.png if any
        map(os.remove, find_by_ext(src_dir, '-nq8.png'))
        # find *.png
        src_lst=filter(lambda fl: not stat.S_ISLNK(os.lstat(fl)[stat.ST_MODE]), # skip symlinks
                    find_by_ext(src_dir, '.png'))

        if options.jpeg:                
            parallel_map(to_jpeg,src_lst)            
            map(os.remove,src_lst)
            
            gmaps=os.path.join(src_dir,'gmaps.html')
            if os.path.exists(gmaps):
                re_sub_file(gmaps,[('tile_ext = .*;','tile_ext = ".jpg";')])
        else:        
            parallel_map(optimize_png,src_lst)

            # 'pngnq' creates *-nq8.png files, so rename files back to original names
            map(lambda f: os.rename(f,f[:-len('-nq8.png')]+'.png'), find_by_ext(src_dir, '-nq8.png'))

        pf('')


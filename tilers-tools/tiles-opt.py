#!/usr/bin/env python

# 2010-10-28 17:18:32 

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
from subprocess import *
import re

try:
    import multiprocessing # available in python 2.6 and above
except:
    multiprocessing=None

def p_f(smth, nl=True):
    if options.quiet:
        return
    s=str(smth)
    if nl: s+='\n'
    sys.stdout.write(s)
    sys.stdout.flush()

def flatten(listOfLists):
    out=[]
    for l in listOfLists:
        out += list(l)
    return out

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

def find_byext(path, ext):
    return flatten([os.path.join(path, name) for name in files if name.endswith(ext)] 
                    for path, dirs, files in os.walk(src_dir))

def command(params,stdin=None):
    process=Popen(params,stdin=PIPE, stdout=PIPE, stderr=PIPE, universal_newlines=True)
    out=process.communicate(stdin)
    if out[1]: p_f(out[1])
    if process.returncode != 0: raise Exception("*** External program failed: %s" % params[0])
    return out[0]

class KeyboardInterruptError(Exception): pass

def proc_file(fl):
    try:
        p_f( '.', False)
        command(['nice','pngnq','-fn', options.colors, fl])
    except KeyboardInterrupt: # http://jessenoller.com/2009/01/08/multiprocessingpool-and-keyboardinterrupt/
        p_f('got KeyboardInterrupt')
        raise KeyboardInterruptError()

if __name__=='__main__':
    #logging.basicConfig(level=logging.INFO)
    logging.basicConfig(level=logging.DEBUG)

    parser = optparse.OptionParser()
    usage = "usage: %prog [options] arg"
    parser.add_option("-n", "--colors", dest="colors", default='256',
        help='Specifies  the  number  of colors to quantize to. Defaults to 256 which is the maximum.  The minimum here is 2')
    parser.add_option("-q", "--quiet", action="store_true")
        
    (options, args) = parser.parse_args()

    if not args: args=['./']
    for src_dir in args:
        print src_dir
        # delete rouge *-nq8.png if any
        map(os.remove, find_byext(src_dir, '-nq8.png'))
        # find *.png
        src_lst=filter(lambda fl: not stat.S_ISLNK(os.lstat(fl)[stat.ST_MODE]), # skip symlinks
                    find_byext(src_dir, '.png'))
        parallel_map(proc_file,src_lst)
        p_f('')

        # 'pngnq' creates *-nq8.png files, so rename files back to original names
        map(lambda f: os.rename(f,f[:-len('-nq8.png')]+'.png'), find_byext(src_dir, '-nq8.png'))


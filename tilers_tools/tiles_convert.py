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
#~ import os
#~ import os.path
#~ import glob
#~ import shutil
import logging
import optparse

from tiler_functions import *

import converter_backend
import converter_xyz
import converter_maemomapper
import converter_sasplanet
import converter_mmap

#----------------------------

def list_formats():

#----------------------------
    for cl in converter_backend.tile_formats:
        print '%10s\t%s%s\t%s' % (
            cl.format,
            'r' if cl.input else ' ',
            'w' if cl.output else ' ',
            cl.__doc__
            )

#----------------------------

def tiles_convert(src_lst, options):

#----------------------------
    for in_class in converter_backend.tile_formats:
        if in_class.input and options.in_fmt == in_class.format:
            break
    else:
        raise Exception('Invalid input format: %s' % options.in_fmt)
    for out_class in converter_backend.tile_formats:
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
    parser.add_option('--from', dest='in_fmt', default='zyx',
        help='input tiles format (default: zyx)')
    parser.add_option('--to', dest='out_fmt', default='mmap',
        help='output tiles format (default: mmap)')
    parser.add_option('-f', '--formats', action='store_true', dest='list_formats',
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

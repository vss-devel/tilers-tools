#!/usr/bin/env python
# -*- coding: utf-8 -*-

# 2011-01-27 11:39:40 

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

import os
import sys
import logging
import re
from optparse import OptionParser

from tiler_functions import *

def kml_parm(hdr,name,lst=False):
    l=re.split('</?%s>' % name,hdr)
    # return only even elements as they are inside <name> </name> 
    return [i.strip() for i in l[1::2]] if lst else l[1].strip()

def find_image(img_path, map_dir):
    imp_path_slashed=img_path.replace('\\','/') # get rid of windows separators
    imp_path_lst=imp_path_slashed.split('/')
    img_patt=imp_path_lst[-1].lower()
    match=[i for i in os.listdir(map_dir if map_dir else '.') if i.lower() == img_patt]
    try:
        return os.path.join(map_dir, match[0])
    except IndexError: raise Exception("*** Image file not found: %s" % img_path)

def overlay2vrt(ol,map_dir):
    ld(ol)
    img_file=kml_parm(ol,'href')
    ld(img_file)
    img_path=find_image(img_file,map_dir)
    base=os.path.splitext(img_path)[0]
    out_vrt= base + '.vrt'        # output VRT file
#    pf(out_vrt)
    if os.path.exists(out_vrt): os.remove(out_vrt)

    # http://trac.osgeo.org/proj/wiki/FAQ#ChangingEllipsoidWhycantIconvertfromWGS84toGoogleEarthVirtualGlobeMercator
    out_srs="+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext +no_defs"
    #"+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +no_defs"

    info_out=command(['gdalinfo',img_path]).splitlines()
    width, height=eval([i for i in info_out if 'Size is' in i][0].replace('Size is',''))
    ld((width, height))
    points=[(0,height), (width,height), (width,0),(0,0)]
    if 'gx:LatLonQuad' in ol:
        refs=[i.split(',') for i in kml_parm(ol,'coordinates').split()]
    else: # assume LatLonBox
        if '<rotation>' in ol: 
            raise Exception("*** rotation in <LatLonBox> is not supported. Suggestion: convert to LatLonQuad")
        n=kml_parm(ol,'north')
        s=kml_parm(ol,'south')
        e=kml_parm(ol,'east')
        w=kml_parm(ol,'west')
        refs=[(w,s),(e,s),(e,n),(w,n)]
    latlong=''.join(['%s %s\n' % (ref[0],ref[1]) for ref in refs])
    ld(latlong)
    refs_proj=[ i.split() for i in 
        command(['proj'] + out_srs.split(), latlong).splitlines()]
    ld(refs_proj)
    gcps=flatten([('-gcp', str(i[0][0]),str(i[0][1]),i[1][0],i[1][1]) for i in zip(points, refs_proj)])
    transl_cmd=['gdal_translate','-of','VRT',img_path, out_vrt,'-a_srs', out_srs]
    if 'Color Table' in info_out:
        transl_cmd=transl_cmd+['-expand','rgb']
    transl_out=command(transl_cmd + gcps)
    logging.info( '"%s" %s' % (img_path,transl_out.strip()))

def kml2vrt(map_path):
    map_dir, map_fname=os.path.split(map_path)
    f=open(map_path, 'r').read()
    if '<GroundOverlay>' not in f: 
        raise Exception("*** Incorrect file: <GroundOverlay> required")
    overlay_lst=kml_parm(f,'GroundOverlay', lst=True) # get list of <GroundOverlay> content
    for ol in overlay_lst:
        overlay2vrt(ol,map_dir)

if __name__=='__main__':
    usage = "usage: %prog [--cut] [--dest-dir=DEST_DIR] MAP_file..."
    parser = OptionParser(usage=usage,
        description="simple KML converter into GDAL .VRT format")
    parser.add_option("-d", "--debug", action="store_true", dest="debug")
    parser.add_option("-t", "--dest-dir", dest="dest_dir", default='',
        help='destination directory (default: current)')

    options, args = parser.parse_args()
    if not args:
        parser.error('No input file(s) specified')
    logging.basicConfig(level=logging.DEBUG if options.debug else logging.INFO)

    for f in args:
        kml2vrt(f)

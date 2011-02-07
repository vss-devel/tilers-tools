#!/usr/bin/env python
# -*- coding: utf-8 -*-

# 2011-02-07 13:33:49 

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

import os
import logging
import locale

from optparse import OptionParser

from tiler_functions import *
from translate2gdal import *

datum_map={    
    'WGS84':                '+datum=WGS84',
    'NAD83':                '+datum=NAD83',
    'ED50':                 '+towgs84=-84.0000,-97.0000,-117.0000 +ellps=intl',
    'POTSDAM':              '+datum=potsdam',
    'EUROPEAN 1950':        '+towgs84=-84.0000,-97.0000,-117.0000 +ellps=intl',
    'EUROPEAN 1950 (NORWAY FINLAND)':
        '+towgs84=-85,-95,-120 +ellps=intl', #http://earth-info.nga.mil/GandG/coordsys/onlinedatum/CountryEuropeTable.html
    'ROMA DATUM 1940':      '+towgs84=-104.1,-49.1,-9.9,0.971,-2.917,0.714,-11.68 +ellps=intl',
    'ROMA 1940':            '+towgs84=-104.1,-49.1,-9.9,0.971,-2.917,0.714,-11.68 +ellps=intl',
    'HERMANSKOGEL DATUM':   '+datum=hermannskogel',
    'OSGB36':               '+towgs84=446.448,-125.157,542.060,0.1502,0.2470,0.8421,-20.4894 +ellps=airy',
    'RT90 (SWEDEN)': 
        '+towgs84=414.0978567149,41.3381489658,603.0627177516,-0.8550434314,2.1413465185,-7.0227209516,0 +ellps=bessel', # http://sv.wikipedia.org/wiki/RT_90
    
    # 'LOCAL DATUM'
    # 'LOCAL DATUM UNKNOWN'
    }

proj_map={
    'MERCATOR':                                 '+proj=merc',
    }

class GeoNosMap(MapTranslator):

    def get_header(self): 
        'read map header'
        with open(self.map_file, 'rU') as f:
            hdr=[[i.strip() for i in l.decode('cp1252','ignore').split('=')] for l in f]
        if not (hdr and hdr[0][0] == '[MainChart]'): 
            raise Exception(" Invalid file: %s" % self.map_file)
        ld(hdr)
        return hdr

    def hdr_parms(self, patt): 
        'filter header for params starting with "patt"'
        return [i[1] for i in self.header if i[0].startswith(patt)]

    def hdr_parms2list(self, patt):
        return [s.split() for s in self.hdr_parms(patt)]
        
    def get_dtm(self):
        'get DTM northing, easting'
        dtm_parm=options.dtm_shift
        if dtm_parm is None:
            try:
                dtm_parm=[self.hdr_parms(i)[0] for i in ('Longitude Offset','Latitude Offset')]
                ld('DTM',dtm_parm)
            except IndexError: # DTM not found
                ld('DTM not found')
                dtm_parm=[0,0]
        dtm=[float(s) for s in dtm_parm]
        return dtm if dtm != [0,0] else None

    def get_refs(self):
        'get a list of geo refs in tuples'
        refs=[(
            (int(i[3]),int(i[2])),                  # pixel
            (float(i[0]),float(i[1]))               # lat/long
            ) for i in self.hdr_parms2list('Point')]
        ld('refs',refs)
        return refs

    def get_plys(self):
        'boundary polygon'
        plys_ll=[(float(i[1]),float(i[0])) for i in self.hdr_parms2list('Vertex')]
        return [((),i) for i in plys_ll]
        
    def get_srs(self):
        options=self.options
        refs=self.refs
        dtm=None
        if options.srs:
            return(self.options.srs,refs)
        if options.proj:
            proj=options.proj
        else:
            proj_id=self.hdr_parms('Projection')[0]
            #parm_lst=self.hdr_parms('Projection Setup')[0]
            try:
                proj=[proj_map[proj_id]]
            except KeyError: 
                raise Exception("*** Unsupported projection (%s)" % proj_id)
            if '+proj=' in proj[0]: # overwise assume it already has a full data defined
                parms=[]
                # setup a central meridian artificialy to allow charts crossing meridian 180
                leftmost=min(refs,key=lambda r: r[0][0])
                rightmost=max(refs,key=lambda r: r[0][0])
                ld('leftmost',leftmost,'rightmost',rightmost)
                if leftmost[1][0] > rightmost[1][0] and '+lon_0=' not in proj[0]:
                    parms.append('+lon_0=%i' % int(leftmost[1][0]))
                if parms:
                    proj.extend(parms)
        datum_id=self.hdr_parms('Datum')[0]
        logging.info(' %s, %s' % (datum_id,proj_id))
        if options.datum: 
            datum=options.datum
        elif options.force_dtm or options.dtm_shift:
            datum='+datum=WGS84'
            dtm=self.get_dtm() # get northing, easting to WGS84 if any
        else:
            try:
                datum=datum_map[datum_id] 
            except KeyError: 
                dtm=self.get_dtm() # get northing, easting to WGS84 if any
                if dtm: 
                    logging.warning(' Unknown datum %s, trying WGS 84 with DTM shifts' % datum_id)
                    datum='+datum=WGS84'
                else: # assume DTM is 0,0
                    logging.warning(' Unknown datum %s, trying WGS 84' % datum_id)
                    datum='+datum=WGS84'
        srs=' '.join(proj)+' '+datum+' +nodefs'
        ld(srs)
        return srs,dtm

    def get_raster(self):
        name_patt=self.hdr_parms('Bitmap')[0].lower()
        map_dir,map_fname=os.path.split(self.map_file)
        dir_lst=[i.decode(locale.getpreferredencoding(),'ignore') 
                    for i in os.listdir(map_dir if map_dir else '.')]
        match=[i for i in dir_lst if i.lower() == name_patt]
        try:
            fn=match[0]
            ld(map_dir, fn)
            img_file=os.path.join(map_dir, fn)
        except:
            raise Exception("*** Image file not found: %s" % img_path)
        return img_file

    def get_size(self):
        return None

    def get_name(self):
        return self.hdr_parms('Name')[0]
# GeoNosMap

class Opt(object):
    def __init__(self,**dictionary):
        self.dict=dictionary
    def __getattr__(self, name):
        return self.dict.setdefault(name,None)

def proc_src(src):
    GeoNosMap(src,options=options).convert()

if __name__=='__main__':
    usage = "usage: %prog <options>... KAP_file..."
    parser = OptionParser(usage=usage,
        description="Extends builtin GDAL's support for BSB charts: WGS84 northing/easting, more projections, border line. "
        "The script converts .kap file with into GDAL .vrt, optionally clipping it out accroding to BSB border line.")
    parser.add_option("-d", "--debug", action="store_true", dest="debug")
    parser.add_option("-q", "--quiet", action="store_true", dest="quiet")
    parser.add_option("-t", "--dest-dir", default=None,
        help='destination directory (default: current)')
    parser.add_option("-l", "--long-names", action="store_true", 
        help='give an output file a long name')
    parser.add_option("--get-cutline", action="store_true", 
        help='print cutline polygon from KAP file then exit')
    parser.add_option("--expand", choices=('gray','rgb','rgba'),
        help='expose a dataset with 1 band with a color table as a dataset with 3 (RGB) or 4 (RGBA) bands')
    parser.add_option("--no-cut-file", action="store_true", 
        help='do not create a file with a cutline polygon from KAP file')
    parser.add_option("--force-dtm", action="store_true", 
        help='force using BSB datum shift to WGS84 instead of native BSB datum')
    parser.add_option("--dtm-shift",dest="dtm_shift",default=None,metavar="SHIFT_LAT,SHIFT_LON",
        help='override DTM: BSB northing, easting (in seconds!)')
    parser.add_option("--srs", default=None,
        help='override full chart with PROJ.4 definition of the spatial reference system')
    parser.add_option("--datum", default=None,
        help="override chart's datum (PROJ.4 definition)")
    parser.add_option("--proj", default=None,
        help="override chart's projection (BSB definition)")
    parser.add_option("--bsb-datum", default=None,dest="datum_id",
        help="override chart's datum (BSB definition)")
    parser.add_option("--bsb-proj", default=None,dest="proj_id",
        help="override chart's projection (BSB definition)")
    parser.add_option("--last-column-bug", action="store_true", 
        help='some BSB files are missing value for last column, here is a workaround')
    parser.add_option("--broken-raster", action="store_true", 
        help='try to workaround some BSB broken rasters (requires "convert" from ImageMagick)')

    (options, args) = parser.parse_args()
    
    if not args:
        parser.error('No input file(s) specified')

    logging.basicConfig(level=logging.DEBUG if options.debug else 
        (logging.ERROR if options.quiet else logging.INFO))

    ld(os.name)
    ld(options)

    map(proc_src,args)


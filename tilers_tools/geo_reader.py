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

from __future__ import with_statement

import os
import logging
import locale

from optparse import OptionParser

from tiler_functions import *
from base_reader import *

datum_map={
    'WGS84':                '+datum=WGS84',
    'NAD83':                '+datum=NAD83',
    'ED50':                 '+towgs84=-84.0000,-97.0000,-117.0000 +ellps=intl',
    'EUROPEAN':             '+towgs84=-84.0000,-97.0000,-117.0000 +ellps=intl',
    'WGS72':                '+ellps=WGS72 +towgs84=0,0,4.5,0,0,0.554,0.219',
    'WGS1972':              '+ellps=WGS72 +towgs84=0,0,4.5,0,0,0.554,0.219',
    'MERCHICH':             '+ellps=clrk80 +towgs84=31,146,47,0,0,0,0',
    # 'LOCAL DATUM'
    # 'LOCAL DATUM UNKNOWN'
    }

proj_map={
    'MERCATOR':                                 '+proj=merc',
    }

class GeoNosMap(MapTranslator):
    magic='[MainChart]'

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
        denominator=3600 # seconds if options.dtm_shift
        if dtm_parm is None:
            denominator=1 # degrees otherwise
            try:
                dtm_parm=[self.hdr_parms(i)[0] for i in ('Longitude Offset','Latitude Offset')]
                ld('DTM',dtm_parm)
            except IndexError: # DTM not found
                ld('DTM not found')
                dtm_parm=[0,0]
        dtm=[float(s)/denominator for s in dtm_parm]
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
        with open(self.img_file) as img:
            hdr=img.readline()
        assert hdr.startswith('NOS/')
        patt='RA='
        sz=hdr[hdr.index(patt)+len(patt):].split(',')[2:4]
        return map(int,sz)

    def get_name(self):
        return self.hdr_parms('Name')[0]
# GeoNosMap

if __name__=='__main__':

    print('\nPlease use translate2gdal.py\n')
    sys.exit(1)


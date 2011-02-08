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
from translate2gdal import *

datum_map={    
    'WGS84':                '+datum=WGS84',
    'NAD83':                '+datum=NAD83',
    'POTSDAM':              '+datum=potsdam',
    'EUROPEAN DATUM 1950' : '+towgs84=-84.0000,-97.0000,-117.0000 +ellps=intl',
    'EUROPEAN 1950':        '+towgs84=-84.0000,-97.0000,-117.0000 +ellps=intl',
    'EUROPEAN 1950 (NORWAY FINLAND)':
        '+towgs84=-85,-95,-120 +ellps=intl', #http://earth-info.nga.mil/GandG/coordsys/onlinedatum/CountryEuropeTable.html
    'ROMA DATUM 1940':      '+towgs84=-104.1,-49.1,-9.9,0.971,-2.917,0.714,-11.68 +ellps=intl',
    'ROMA 1940':            '+towgs84=-104.1,-49.1,-9.9,0.971,-2.917,0.714,-11.68 +ellps=intl',
    'HERMANSKOGEL DATUM':   '+datum=hermannskogel',
    'OSGB36':               '+towgs84=446.448,-125.157,542.060,0.1502,0.2470,0.8421,-20.4894 +ellps=airy',
    'RT90 (SWEDEN)': 
        '+towgs84=414.0978567149,41.3381489658,603.0627177516,-0.8550434314,2.1413465185,-7.0227209516,0 +ellps=bessel', # http://sv.wikipedia.org/wiki/RT_90
    u'SYSTEM KÜSTE':        '+datum=potsdam', # ???
    
    # 'LOCAL DATUM'
    # 'LOCAL DATUM UNKNOWN'
    }

datum_guess={ # guess the datum by a comment/copyright string pattern
    'Croatia':
        # http://spatial-analyst.net/wiki/index.php?title=MGI_/_Balkans_coordinate_systems
        '+ellps=bessel +towgs84=550.499,164.116,475.142,5.80967,2.07902,-11.62386,0.99999445824',
        # http://earth-info.nga.mil/GandG/coordsys/onlinedatum/DatumTable.html
        #'+datum=hermannskogel',
        #'+ellps=bessel +towgs84=682,-203,480',
    }

proj_knp={ # projection parameters
    'MERCATOR':
        {'PROJ': '+proj=merc',  'PP': '+lat_ts='},
    'TRANSVERSE MERCATOR':
        {'PROJ': '+proj=tmerc', 'PP': '+lon_0=', 'P1': '+lat_0=', 'P2': '+k=', 'P3': '+y_0=', 'P4': '+x_0='},
    'UNIVERSAL TRANSVERSE MERCATOR':
        {'PROJ': '+proj=utm', 'PP': '+lon_0='},
    'GNOMONIC':
        {'PROJ': '+proj=gnom',  'PP': '+lon_0=', 'P1': '+lat_0='},
    'LAMBERT CONFORMAL CONIC':
        {'PROJ': '+proj=lcc', 'PP': '+lon_0='},
    'POLYCONIC':
        {'PROJ': '+proj=poly',  'PP': '+lon_0='},
    'SWEDISH GRID':
        {'PROJ': '+proj=tmerc +lon_0=15.808277777778 +x_0=1500000 +y_0=0'},
    }

proj_knq={ # extra projection parameters for BSB v. 3.xx
    'TRANSVERSE MERCATOR':
        {'P1': '+lon_0=', 'P2': '+k=', 'P3': '+lat_0='}, # P3 - guess
    'LAMBERT CONFORMAL CONIC':
        {'P1': '+lon_0=', 'P2': '+lat_1=', 'P3': '+lat_2='},
    'POLYCONIC':
        {'P1': '+lon_0=', 'P2': '+lat_0='}, # P2 - guess
    }

class BsbKapMap(MapTranslator):

    def get_header(self): 
        'read map header'
        header=[]
        with open(self.map_file,'rU') as f:
            for l in f:
                if '\x1A' in l:
                    break
                l=l.decode('cp1252','ignore')
                if l.startswith((' ','\t')):
                    header[-1] += ','+l.strip()
                else:
                    header.append(l.strip())
        ld(header)
        if not (header and any(((s.startswith('BSB/') or s.startswith('KNP/')) for s in header))): 
            raise Exception(" Invalid file: %s" % self.map_file)
        return header

    def hdr_parms(self, patt): 
        'filter header for params starting with "patt/", if knd is empty then return comment lines'
        if patt != '!': 
            patt += '/'
        return [i[len(patt):] for i in self.header if i.startswith(patt)]

    def hdr_parms2list(self, knd):
        return [i.split(',') for i in self.hdr_parms(knd)]

    def hdr_parm2dict(self, knd):
        out={}
        for i in self.hdr_parms2list(knd)[0]:
            if '=' in i:
                (key,val)=i.split('=')
                out[key]=val
            else:
                out[key] += ','+i
        return out

    def assemble_parms(self,parm_map,parm_info):    
        check_parm=lambda s: (s not in ['NOT_APPLICABLE','UNKNOWN']) and s.replace('0','').replace('.','')
        res=' '.join([parm_map[i]+parm_info[i] for i in parm_map
                        if  i in parm_info and check_parm(parm_info[i])])
        return ' '+res if res else ''
        
    def get_dtm(self):
        'get DTM northing, easting'
        dtm_parm=options.dtm_shift
        if dtm_parm is None:
            try:
                dtm_parm=self.hdr_parms2list('DTM')[0]
                ld('DTM',dtm_parm)
            except IndexError: # DTM not found
                ld('DTM not found')
                dtm_parm=[0,0]
        dtm=[float(s)/3600 for s in reversed(dtm_parm)]
        return dtm if dtm != [0,0] else None

    def get_refs(self):
        'get a list of geo refs in tuples'
        refs=[(
            (int(i[1]),int(i[2])),                  # pixel
            (float(i[4]),float(i[3]))               # lat/long
            ) for i in self.hdr_parms2list('REF')]
        ld('refs',refs)
        return refs

    def get_plys(self):
        'boundary polygon'
        plys_ll=[(float(i[2]),float(i[1])) for i in self.hdr_parms2list('PLY')]
        return [((),i) for i in plys_ll]
        
    def get_srs(self):
        'returns srs for the BSB chart projection and srs for the REF points'
        options=self.options
        refs=self.refs
        dtm=None
        # Get a list of geo refs in tuples
        if options.srs:
            return options.srs
        # evaluate chart's projection
        proj=options.proj
        if not proj:
            knp_info=self.hdr_parm2dict('KNP')
            ld(knp_info)
            proj_id=if_set(options.proj_id,knp_info['PR'])
            try:            
                knp_parm=proj_knp[proj_id.upper()]
            except KeyError: 
                raise Exception(' Unsupported projection %s' % proj_id)
            # get projection and parameters
            proj=knp_parm['PROJ']
            if '+proj=utm' in proj[0]: # UTM
                # GDAL 1.7.2 doesn't seem to make use of lon_0 with UTM, but BSB doesn't use zones
                northing='10000000' if refs[0][1][1] < 0 else '0' # Southern hemisphere?
                proj='+proj=tmerc +k=0.9996 +x_0=500000 +y_0=%s +lon_0=%s' % (northing,knp_info['PP'])
            else:
                try: # extra projection parameters for BSB 3.xx, put them before KNP parms
                    knq_info=self.hdr_parm2dict('KNQ')
                    ld(knq_info)
                    knq_parm=proj_knq[proj_id.upper()]
                    proj+=assemble_parms(knq_parm,knq_info)
                except IndexError:  # No KNQ
                    pass
                except KeyError:    # No such proj in KNQ map
                    pass
                proj+=self.assemble_parms(knp_parm,knp_info)
        # setup a central meridian artificialy to allow charts crossing meridian 180
        leftmost=min(refs,key=lambda r: r[0][0])
        rightmost=max(refs,key=lambda r: r[0][0])
        ld('leftmost',leftmost,'rightmost',rightmost)
        if leftmost[1][0] > rightmost[1][0] and '+lon_0=' not in proj:
            proj+=' +lon_0=%i' % int(leftmost[1][0])
        
        # evaluate chart's datum
        datum_id=if_set(options.datum_id,knp_info['GD'])
        logging.info(' %s, %s' % (datum_id,proj_id))
        datum=options.datum
        if datum:
            pass
        elif options.force_dtm or options.dtm_shift:
            datum='+datum=WGS84'
            dtm=self.get_dtm() # get northing, easting to WGS84 if any
        elif not '+proj=' in proj: 
            datum='' # assume datum is defined already
        else:
            try:
                datum=datum_map[datum_id.upper()]
            except KeyError: 
                # try to guess the datum by comment and copyright string(s)
                crr=' '.join(self.hdr_parms('!')+self.hdr_parms('CRR'))
                try:
                    datum=[datum_guess[crr_patt] 
                        for crr_patt in datum_guess if crr_patt in crr][0]
                    logging.warning(' Unknown datum "%s", guessed as "%s"' % (datum_id,datum))
                except IndexError:
                    # datum still not found
                    dtm=self.get_dtm() # get northing, easting to WGS84 if any
                    if dtm: 
                        logging.warning(' Unknown datum %s, trying WGS 84 with DTM shifts' % datum_id)
                        datum='+datum=WGS84'
                    else: # assume DTM is 0,0
                        logging.warning(' Unknown datum %s, trying WGS 84' % datum_id)
                        datum='+datum=WGS84'
        srs=proj+' '+datum+' +nodefs'
        ld(srs)
        return srs,dtm

    def get_raster(self):
        return self.map_file

    def get_size(self):
        bsb_info=self.hdr_parm2dict('BSB') # general BSB parameters
        return map(int,bsb_info['RA'].split(','))

    def get_name(self):
        bsb_info=self.hdr_parm2dict('BSB') # general BSB parameters
        bsb_name=bsb_info['NA']
        return bsb_name
# BsbKapMap

class Opt(object):
    def __init__(self,**dictionary):
        self.dict=dictionary
    def __getattr__(self, name):
        return self.dict.setdefault(name,None)

def proc_src(src):
    BsbKapMap(src,options=options).convert()

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


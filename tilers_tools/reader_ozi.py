#!/usr/bin/env python
# -*- coding: utf-8 -*-

# 2011-02-08 17:38:41 

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

from __future__ import with_statement

import os
import logging
import locale

from optparse import OptionParser

from tiler_functions import *
from reader_backend import *

class OziMap(MapTranslator):
    magic='OziExplorer Map Data File'
    data_file='reader_ozi_data.csv'

    proj_parms=(
        '+lat_0=', # 1. Latitude Origin
        '+lon_0=', # 2. Longitude Origin
        '+k=',     # 3. K Factor
        '+x_0=',   # 4. False Easting
        '+y_0=',   # 5. False Northing
        '+lat_1=', # 6. Latitude 1
        '+lat_2=', # 7. Latitude 2
        '+h=',     # 8. Height - used in the Vertical Near-Sided Perspective Projection
                   # 9. Sat - not used
                   #10. Path - not used
        )
        
    def load_data(self):
        'load datum definitions, ellipses, projections from a file'
        self.datum_map={}
        self.ellps_map={}
        self.proj_map={}
        csv_map={
            'datum': (self.datum_map,self.ini_lst),
            'ellps': (self.ellps_map,self.ini_lst),
            'proj': (self.proj_map,self.ini_lst),
            }
        self.load_csv(self.data_file,csv_map)

    def get_header(self): 
        'read map header'
        with open(self.map_file, 'rU') as f:
            hdr=[[i.strip() for i in l.split(',')] for l in f]
        if not (hdr and hdr[0][0].startswith('OziExplorer Map Data File')): 
            raise Exception(" Invalid file: %s" % self.map_file)
        ld(hdr)
        return hdr

    def hdr_parms(self, patt): 
        'filter header for params starting with "patt"'
        return [i for i in self.header if i[0].startswith(patt)]
        
    def get_refs(self):
        'get a list of geo refs in tuples'
        points=[i for i in self.hdr_parms('Point') if i[2] != ''] # Get a list of geo refs
        if points[0][14] != '': # refs are cartesian
            refs=RefPoints(self,[(
                    i[0],                                   # id
                    (int(i[2]),int(i[3])),                  # pixel
                    (float(i[14]),float(i[15])),            # cartesian coords
                    ) for i in points],
                cartesian=True,
                extra=[(i[13],i[16]) for i in points]   # zone, hemisphere
                )
        else:
            refs=RefPoints(self,[(
                i[0],                                   # id
                (int(i[2]),int(i[3])),                  # pixel
                (dms2dec(*i[9:12]), dms2dec(*i[6:9])),  # lat/long
                ) for i in points])
        return refs

    def get_plys(self):
        'boundary polygon'
        ply_pix=[(int(i[2]),int(i[3])) for i in self.hdr_parms('MMPXY')]    # Moving Map border pixels
        ply_ll=[(float(i[2]),float(i[3])) for i in self.hdr_parms('MMPLL')] # Moving Map border lat,lon
        ids=[i[0] for i in self.hdr_parms('MMPXY')]    # Moving Map border pixels
        plys=RefPoints(self,ids=ids,pixels=ply_pix,coords=ply_ll)
        return plys

    def get_dtm(self):
        'get DTM northing, easting'
        dtm=[float(s)/3600 for s in self.header[4][2:4]]
        return dtm if dtm != [0,0] else None

    def get_proj_id(self):
        return self.hdr_parms('Map Projection')[0][1]
    
    def get_proj(self):
        proj_id=self.get_proj_id()
        parm_lst=self.hdr_parms('Projection Setup')[0]
        try:
            proj=self.proj_map[proj_id][0:1]
        except KeyError: 
            raise Exception("*** Unsupported projection (%s)" % proj_id)
        if '+proj=' in proj[0]: # overwise assume it already has a full data defined
            # get projection parameters
            proj.extend([ i[0]+i[1] for i in zip(self.proj_parms,parm_lst[1:]) 
                            if i[1].translate(None,'0.')])
            if '+proj=utm' in proj[0]:
                if not refs[0][1]: # refs are cartesian with a zone defined
                    hemisphere=refs[0][3]
                    utm_zone=int(refs[0][4])
                    proj.append('+zone=%i' % utm_zone)
                    if hemisphere != 'N': 
                        proj.append('+south')
                else: # refs are lat/long, then find zone, hemisphere
                    # doesn't seem to have central meridian for UTM
                    lon,lat=refs[0][1]
                    zone=(lon + 3 % 360) // 6 + 30
                    proj.append('+zone=%d' % zone)
                    if lat < 0: 
                        proj.append('+south')
        return proj

    def get_datum_id(self):
        return self.header[4][0]

    def get_datum(self):
        datum_id=self.get_datum_id()
        try:
            datum_def=self.datum_map[datum_id]
            if datum_def[5]: # PROJ4 datum defined ?
                datum=[datum_def[5]]
            else:
                datum=['+towgs84=%s,%s,%s' % tuple(datum_def[2:5])]               
                ellps_id=datum_def[1]
                ellps_def=self.ellps_map[ellps_id]
                ellps=if_set(ellps_def[2])
                if ellps:
                    datum.append(ellps)
                else:
                    datum.append('+a=%s',ellps_def[0])
                    datum.append('+rf=%s',ellps_def[1])                        
        except KeyError: 
            raise Exception("*** Unsupported datum (%s)" % datum_id)
        return datum

    try_encodings=(locale.getpreferredencoding(),'utf_8','cp1251','cp1252')

    def get_raster(self):
        img_path=self.header[2][0]
        img_path_slashed=img_path.replace('\\','/') # get rid of windows separators
        img_path_lst=os.path.split(img_path_slashed)
        img_fname=img_path_lst[-1]
        map_dir,map_fname=os.path.split(self.map_file)
        dir_lst=[i.decode(locale.getpreferredencoding(),'ignore') 
                    for i in os.listdir(map_dir if map_dir else '.')]
        # try a few encodings
        for enc in self.try_encodings:
            name_patt=img_fname.decode(enc,'ignore').lower()
            match=[i for i in dir_lst if i.lower() == name_patt]
            if match:
                fn=match[0]
                ld(map_dir, fn)
                img_file=os.path.join(map_dir, fn)
                break
        else:
            raise Exception("*** Image file not found: %s" % img_path)
        return img_file

    def get_size(self):
        return map(int,self.hdr_parms( 'IWH')[0][2:])

    def get_name(self):
        ozi_name=self.header[1][0]
        # try a few encodings
        for enc in self.try_encodings:
            try:
                if enc == 'cp1251' and any([ # ascii chars ?
                        ((c >= '\x41') and (c <= '\x5A')) or 
                        ((c >= '\x61') and (c <= '\x7A')) 
                            for c in ozi_name]):
                    continue # cp1251 name shouldn't have any ascii
                ozi_name=ozi_name.decode(enc)
                break
            except:
                pass
        ld('ozi_name',ozi_name)
        return ozi_name
                
# OziMap

if __name__=='__main__':

    print('\nPlease use translate2gdal.py\n')
    sys.exit(1)


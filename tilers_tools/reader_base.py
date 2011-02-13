#!/usr/bin/env python
# -*- coding: utf-8 -*-

# 2011-02-07 11:41:56 

###############################################################################
# Copyright (c) 2011, Vadim Shlyakhov
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
import csv

from tiler_functions import *

try:
    from osgeo import gdal
    from osgeo import osr
    from osgeo.gdalconst import *
    gdal.TermProgress = gdal.TermProgress_nocb
except ImportError:
    import osr
    import gdal
    from gdalconst import *

def dms2dec(degs='0',mins='0',ne='E',sec='0'):
    return (float(degs)+float(mins)/60+float(sec)/3600)*(-1 if ne in ('W','S') else 1 )

def dest_path(src,dest_dir,ext='',template='%s'):
    src_dir,src_file=os.path.split(src)
    base,sext=os.path.splitext(src_file)
    dest=(template % base)+ext
    if not dest_dir:
        dest_dir=src_dir
    if dest_dir:
        dest='%s/%s' % (dest_dir,dest)
    ld(base,dest)
    return dest

class Opt(object):
    def __init__(self,**dictionary):
        self.dict=dictionary
    def __getattr__(self, name):
        return self.dict.setdefault(name,None)

class RefPoints(object):
    'source geo-reference points and polygons'
    def __init__(self,owner,ref_lst=None,ids=None,pixels=None,lonlat=None,cartesian=None,extra=None):
        self.owner=owner
        self._ids=ids
        self._pixels=pixels
        self._lonlat=lonlat
        self._cartesian=cartesian
        self._extra=extra
        if ref_lst:
            transposed=[i if any(i) else None for i in zip(*ref_lst)]
            pad=[None]*(4-len(transposed))
            transposed.extend(pad)
            self._ids,self._pixels,self._lonlat,self._cartesian=transposed

        items=len(filter(None,(self._pixels,self._lonlat,self._cartesian))[0])
        if not self._ids:
            self._ids=map(str,range(1,items+1))

        if items == 2:
            logging.warning(' Only 2 reference points: assuming the chart is north alligned')
            for i in (self._pixels,self._lonlat,self._cartesian):
                if i is not None:
                    i.append((i[0][0],i[1][1]))
            self._ids.append('3')

        self._ids=[s.encode('utf-8') for s in self._ids]

    def ids(self):
        return self._ids

    def pixels(self,geo_transform=None):
        if self._pixels is None: # in a case of a mere polygon
            assert geo_transform is not None # must have geo_transform
            inv_gt=gdal.InvGeoTransform(geo_transform)
            assert inv_gt[0]
            self._pixels=[gdal.ApplyGeoTransform(inv_gt[1],x,y) for x,y,z in self.projected()]
        return self._pixels
        
    def projected(self):
        srs=self.owner.srs
        if srs.startswith('+proj=latlong'):
            return self.longlat()
        if self._cartesian is None:
            srs_proj = osr.SpatialReference()
            srs_proj.ImportFromProj4(srs)
            srs_geo = osr.SpatialReference()
            srs_geo.ImportFromProj4('+proj=latlong')
            srs_geo.CopyGeogCSFrom(srs_proj)

            tr=osr.CoordinateTransformation(srs_geo,srs_proj)
            
            geo_coords=self.longlat()
            self._cartesian=tr.TransformPoints(geo_coords)
            ld('geo_coords',geo_coords,'cartesian',self._cartesian)
        return self._cartesian

    def longlat(self):
        dtm=self.owner.dtm
        lonlat=self._lonlat
        if not dtm:
            return lonlat
        else: # alter refs as per DTM values
            return [(lon+dtm[0],lat+dtm[1]) for lon,lat in self._lonlat]

    def over_180(self):
        lonlat=self._lonlat
        if lonlat: # refs are lat/long
            leftmost=min(lonlat,key=lambda r: r[0])
            rightmost=max(lonlat,key=lambda r: r[0])
            ld('leftmost',leftmost,'rightmost',rightmost)
            if leftmost[0] > rightmost[0]:
                return leftmost[0]
        return None
        
class MapTranslator(object):
        
    def __init__(self,src_file,options=None):
        self.options=options

        self.load_data() # load datum definitions, ellipses, projections
        self.map_file=src_file.decode(locale.getpreferredencoding(),'ignore')
        self.header=self.get_header()       # Read map header
        self.img_file=self.get_raster()
        self.name=self.get_name()
        logging.info(' %s : %s (%s)' % (self.map_file,self.name,self.img_file))

        self.refs=self.get_refs()           # fetch reference points
        self.srs,self.dtm=self.get_srs()    # estimate SRS

    def load_csv(self,csv_file,csv_map):
        'load datum definitions, ellipses, projections from a file'
        data_dir=sys.path[0]
        ld('data_dir',data_dir)
        csv.register_dialect('strip', skipinitialspace=True)
        with open(os.path.join(data_dir,csv_file),'rb') as data_f:
            data_csv=csv.reader(data_f,'strip')
            for row in data_csv:
                row=[s.decode('utf-8') for s in row]
                ld(row)
                try:
                    dct,unpack=csv_map[row[0]]
                    unpack(dct,row)
                except IndexError:
                    pass
                except KeyError:
                    pass
        for dct,func in csv_map.values():
            ld(dct)
            
    def ini_lst(self,dct,row):
        dct[row[1]]=row[2:]

    def ini_map(self,dct,row):
        vlst=row[2:]
        keys=vlst[0::2]
        vals=vlst[1::2]
        dct[row[1]]=dict(zip(keys,vals))

    def get_srs(self):
        'returns srs for the map, and DTM shifts if any'
        options=self.options
        dtm=None
        srs=[]
        datum_id,proj_id='',''
        if options.srs:
            return(self.options.srs,None)
        # evaluate chart's projection
        if options.proj:
            srs.append(options.proj)
        else:
            srs,proj_id=self.get_proj()

        # setup a central meridian artificialy to allow charts crossing meridian 180
        leftmost=self.refs.over_180()
        if leftmost and '+lon_0=' not in proj:
            srs.append(' +lon_0=%i' % int(leftmost))
        
        # evaluate chart's datum
        if options.datum: 
            srs.append(options.datum)
        elif options.force_dtm or options.dtm_shift:
            dtm=self.get_dtm() # get northing, easting to WGS84 if any
            srs.append('+datum=WGS84')
        elif not '+proj=' in srs[0]: 
            pass # assume datum is defined already
        else:
            datum,datum_id=self.get_datum()
            srs.extend(datum)
        srs.extend(['+wktext','+nodefs'])
        ld(srs)
        logging.info(' %s, %s' % (datum_id,proj_id))
        return str(' '.join(srs)),dtm

    def convert(self,dest=None):
        options=self.options
        
        gdal.UseExceptions()

        if dest:
            base=os.path.split(dest)[0]
        else:
            name_patt=self.map_file
            if options.as_image:
                name_patt=self.img_file
            base=dest_path(name_patt,self.options.dest_dir)
            if options.long_name:
                base+=' - '+self.name
        dest_dir=os.path.split(base)[0]
        out_format='VRT'
        ext='.'+out_format.lower()
        dest_file= os.path.basename(base+ext).encode('utf-8') # output file
        img_path=os.path.relpath(self.img_file,dest_dir).encode('utf-8')

        try:
            cdir=os.getcwd()
            if dest_dir:
                os.chdir(dest_dir)

            src_ds = gdal.Open(img_path,GA_ReadOnly)
            dest_drv = gdal.GetDriverByName(out_format)
            dest_ds = dest_drv.CreateCopy(dest_file,src_ds,0)
            dest_ds.SetProjection(self.srs)

            refs=self.refs
            #double x = 0.0, double y = 0.0, double z = 0.0, double pixel = 0.0, 
            #double line = 0.0, char info = "", char id = ""
            gcps=[gdal.GCP(c[0],c[1],0,p[0],p[1], '',i)
                    for i,p,c in zip(refs.ids(),refs.pixels(),refs.projected())]
            dest_ds.SetGCPs(gcps,self.srs)
            dest_geotr=gdal.GCPsToGeoTransform(gcps)
            dest_ds.SetGeoTransform(dest_geotr)

            poly,gmt_data=self.cut_poly(dest_ds)
            del dest_ds
        finally:
            os.chdir(cdir)

        if self.options.get_cutline: # print cutline then return
            print poly
            return
        if gmt_data and not self.options.no_cut_file: # create shapefile with a cut polygon
            with open(base+'.gmt','w+') as f:
                f.write(gmt_data)

    gmt_templ='''# @VGMT1.0 @GPOLYGON
# @Jp"%s"
# FEATURE_DATA
>
# @P
%s
'''

    def cut_poly(self,dest_ds):
        plys=self.get_plys()
        if not plys:
            return '',''

        pix_lst=plys.pixels(dest_ds.GetGeoTransform())

        # check if the raster really needs cutting
        width=dest_ds.RasterXSize
        height=dest_ds.RasterYSize
        inside=[i for i in pix_lst # check if the polygon is inside the image border
            if (i[0] > 0 or i[0] < width) or (i[1] > 0 or i[1] < height)]
        if not inside:
            return '',''

        # Create cutline
        poly_proj=plys.projected()
        poly_shape=self.gmt_templ % (self.srs,'\n'.join(['%r %r %r' % tuple(i) for i in poly_proj]))
        poly_wkt='POLYGON((%s))' % ','.join(['%r %r' % tuple(i) for i in pix_lst]) # Create cutline
        return poly_wkt,poly_shape

# MapTranslator


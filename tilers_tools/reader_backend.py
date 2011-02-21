#!/usr/bin/env python
# -*- coding: utf-8 -*-

# 2011-02-16 18:11:18 

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

def dst_path(src,dst_dir,ext='',template='%s'):
    src_dir,src_file=os.path.split(src)
    base,sext=os.path.splitext(src_file)
    dest=(template % base)+ext
    if not dst_dir:
        dst_dir=src_dir
    if dst_dir:
        dest='%s/%s' % (dst_dir,dest)
    ld(base,dest)
    return dest

class Opt(object):
    def __init__(self,**dictionary):
        self.dict=dictionary
    def __getattr__(self, name):
        return self.dict.setdefault(name,None)

def srs2srs(proj4_src,proj4_dst):
    srs_src = osr.SpatialReference()
    srs_src.ImportFromProj4(proj4_src)
    srs_dst = osr.SpatialReference()
    srs_dst.ImportFromProj4(proj4_dst)
    return osr.CoordinateTransformation(srs_src,srs_dst)

class RefPoints(object):
    'source geo-reference points and polygons'
    def __init__(self,owner,ref_lst=None,ids=None,pixels=None,coords=None,cartesian=False,extra=None):
        self.owner=owner
        self._ids=ids
        self._pixels=pixels
        self._coords=coords
        self._cartesian=cartesian
        self._extra=extra
        if ref_lst:
            transposed=[i if any(i) else None for i in zip(*ref_lst)]
            self._ids,self._pixels,self._coords=transposed

        items=len(filter(None,(self._pixels,self._coords))[0])
        if not self._ids:
            self._ids=map(str,range(1,items+1))

        if items == 2:
            logging.warning(' Only 2 reference points: assuming the chart is north alligned')
            for i in (self._pixels,self._coords):
                if i is not None:
                    i.append((i[0][0],i[1][1]))
            self._ids.append('3')
        self._ids=[s.encode('utf-8') for s in self._ids]

    def srs(self):
        return self.owner.srs #if self._cartesian else self.geo_srs()

    def geo_srs(self):
        srs_proj = osr.SpatialReference()
        srs_proj.ImportFromProj4(self.owner.srs)
        srs_geo = osr.SpatialReference()
        srs_geo.CopyGeogCSFrom(srs_proj)
        return srs_geo.ExportToProj4()

    def __iter__(self):
        for i in zip(self._ids,self.pixels(),self.coords()):
            yield i
    
    def pixels(self,dataset=None):
        if self._pixels:
            return self._pixels
#        srs_tr=srs2srs(self.srs(),self.owner.srs)
#        p_dst=srs_tr.TransformPoints(self.coords())
        p_dst=self.coords()
        ld(p_dst)
        pix_tr=gdal.Transformer(dataset,None,['METHOD=GCP_TPS'])
        p_pix,ok=pix_tr.TransformPoints(True,p_dst)
        ld(p_pix)
        assert all(ok)
        return [(p[0],p[1]) for p in p_pix]

    def coords(self):
        if self._cartesian:
            return self._coords
        dtm=self.owner.dtm
        if not dtm:
            dtm=[0,0]
        latlon=[(lon+dtm[0],lat+dtm[1]) for lon,lat in self._coords]
#        return latlon        
        srs_tr=srs2srs(self.geo_srs(),self.owner.srs)
        coords=srs_tr.TransformPoints(latlon)
        return coords      

    def over_180(self):
        if not self._cartesian: # refs are lat/long
            leftmost=min(zip(self._pixels,self._coords),key=lambda r: r[0][0])
            rightmost=max(zip(self._pixels,self._coords),key=lambda r: r[0][0])
            ld('leftmost',leftmost,'rightmost',rightmost)
            if leftmost[1][0] > rightmost[1][0]:
                return leftmost[1][0]
        return None
        
class MapTranslator(object):
        
    def __init__(self,src_file,options=None):
        self.options=options
        gdal.UseExceptions()

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
                #ld(row)
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
        proj4=[]
        logging.info(' %s, %s' % (self.get_datum_id(),self.get_proj_id()))
        if options.srs:
            return(self.options.srs,None)
        
        # evaluate chart's projection
        if options.proj:
            proj4.append(options.proj)
        else:
            proj4=self.get_proj()

        # setup a central meridian artificialy to allow charts crossing meridian 180
        leftmost=self.refs.over_180()
        if leftmost and '+lon_0=' not in proj4[0]:
            proj4.append(' +lon_0=%i' % int(leftmost))
        
        # evaluate chart's datum
        if options.datum: 
            proj4.append(options.datum)
        elif options.force_dtm or options.dtm_shift:
            dtm=self.get_dtm() # get northing, easting to WGS84 if any
            proj4.append('+datum=WGS84')
        elif not '+proj=' in proj4[0]: 
            pass # assume datum is defined already
        else:
            datum=self.get_datum()
            proj4.extend(datum)
        proj4.extend(['+nodefs']) # '+wktext',
        ld(proj4)
        return ' '.join(proj4).encode('utf-8'),dtm

    def convert(self,dest=None):
        options=self.options
        
        if dest:
            base=os.path.split(dest)[0]
        else:
            name_patt=self.map_file
            if options.as_image:
                name_patt=self.img_file
            base=dst_path(name_patt,self.options.dst_dir)
            if options.long_name:
                base+=' - '+self.name
        dst_dir=os.path.split(base)[0]
        out_format='VRT'
        ext='.'+out_format.lower()
        dst_file= os.path.basename(base+ext).encode('utf-8') # output file
        img_path=os.path.relpath(self.img_file,dst_dir).encode('utf-8')

        try:
            cdir=os.getcwd()
            if dst_dir:
                os.chdir(dst_dir)

            src_ds = gdal.Open(img_path,GA_ReadOnly)
            dst_drv = gdal.GetDriverByName(out_format)
            dst_ds = dst_drv.CreateCopy(dst_file,src_ds,0)
            dst_ds.SetProjection(self.srs)

            #double x = 0.0, double y = 0.0, double z = 0.0, double pixel = 0.0, 
            #double line = 0.0, char info = "", char id = ""
            gcps=[gdal.GCP(c[0],c[1],0,p[0],p[1],'',i) for i,p,c in self.refs]
            dst_ds.SetGCPs(gcps,self.refs.srs())
            dst_geotr=gdal.GCPsToGeoTransform(gcps)
            dst_ds.SetGeoTransform(dst_geotr)
            poly,gmt_data=self.cut_poly(dst_ds)
            if poly:
                dst_ds.SetMetadataItem('CUTLINE',poly)
            if self.name:
                dst_ds.SetMetadataItem('DESCRIPTION',self.name.encode('utf-8'))

            del dst_ds # close dataset
            re_sub_file(dst_file, [
                    ('^.*<GeoTransform>.*\n',''),
                    ('^.*<SRS>.*\n','')
                    ])
        finally:
            os.chdir(cdir)

        if self.options.get_cutline: # print cutline then return
            print poly
            return
        if gmt_data and self.options.cut_file: # create shapefile with a cut polygon
            with open(base+'.gmt','w+') as f:
                f.write(gmt_data)

    gmt_templ='''# @VGMT1.0 @GPOLYGON
# @Jp"%s"
# FEATURE_DATA
>
# @P
%s
'''

    def cut_poly(self,dst_ds):
        plys=self.get_plys()
        if not plys:
            return '',''

        pix_lst=plys.pixels(dst_ds)

        # check if the raster really needs cutting
        width=dst_ds.RasterXSize
        height=dst_ds.RasterYSize
        inside=[i for i in pix_lst # check if the polygon is inside the image border
            if (i[0] > 0 or i[0] < width) or (i[1] > 0 or i[1] < height)]
        if not inside:
            return '',''

        # Create cutline
        poly_shape=self.gmt_templ % (self.refs.srs(),'\n'.join(['%r %r' % (i[0],i[1]) for i in plys.coords()]))
        poly_wkt='MULTIPOLYGON(((%s)))' % ','.join(['%r %r' % tuple(i) for i in pix_lst]) # Create cutline
        return poly_wkt,poly_shape

# MapTranslator


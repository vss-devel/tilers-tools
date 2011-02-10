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

    gmt_templ='''# @VGMT1.0 @GPOLYGON
# @Jp"%s"
# FEATURE_DATA
>
# @P
%s'''

    def cut_poly(self,dest_ds,lonlat2proj):
        plys=self.shift_lonlat(self.get_plys(),self.dtm)   # as per dtm value
        if not plys:
            return '',''

        # Create cutline
        dest_geotr=dest_ds.GetGeoTransform()

#        lonlat=''.join(['%r %r\n' % i[1] for i in plys])
        lonlat=[i[1] for i in plys]
        if not plys[0][0]: # convert cutline coordinates to pixel xy using GDAL's navive srs for this raster
#            pix_lines=command(['gdaltransform','-tps','-i','-t_srs','+proj=longlat',out_dataset],
#                                lonlat).splitlines()
#            pix_lst=[(int(i[0]),int(i[1])) for i in pix_lines]
            inv_gt=gdal.InvGeoTransform(dest_geotr)
            pix_lst=[gdal.ApplyGeoTransform(inv_gt,lon,lat) for lon,lat in lonlat]
        else:
            pix_lst=[i[0] for i in plys]
            pix_lines=['%d %d' % i for i in pix_lst]
        poly='POLYGON((%s))' % ','.join(pix_lines) # Create cutline

        width=dest_ds.GetRasterXSize
        height=dest_ds.GetRasterYSize
        inside=[i for i in pix_lst # check if the polygon is inside the image border
            if (i[0] > 0 or i[0] < width) or (i[1] > 0 or i[1] < height)]
        if not inside:
            return '',''

        # convert cutline geo coordinates to the chart's srs
        poly_proj=lonlat2proj.TransformPoints(lonlat)
#        poly_xy=command(['gdaltransform','-tps','-s_srs','+proj=longlat','-t_srs',self.srs],lonlat)
        return poly,self.gmt_templ % (self.srs,poly_xy)

    def shift_lonlat(self,refs,dtm):
        if not dtm:
            return refs
        else:
            # alter refs as per DTM values
            split=zip(*refs) # split refs
            split[1]=[(lonlat[0]+dtm[0],lonlat[1]+dtm[1]) for lonlat in split[1]]
            return zip(*split) # repack refs

    def get_srs(self):
        'returns srs for the map, and DTM shifts if any'
        options=self.options
        refs=self.refs
        dtm=None
        srs=[]
        datum_id,proj_id='',''
        # Get a list of geo refs in tuples
        if options.srs:
            return(self.options.srs,None)
        # evaluate chart's projection
        if options.proj:
            srs.append(options.proj)
        else:
            srs,proj_id=self.get_proj()
        # setup a central meridian artificialy to allow charts crossing meridian 180
        if refs[0][1]: # refs are lat/long
            leftmost=min(refs,key=lambda r: r[0][0])
            rightmost=max(refs,key=lambda r: r[0][0])
            ld('leftmost',leftmost,'rightmost',rightmost)
            if leftmost[1][0] > rightmost[1][0] and '+lon_0=' not in proj:
                srs.append(' +lon_0=%i' % int(leftmost[1][0]))
        
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

        refs=self.shift_lonlat(self.refs,self.dtm)   # as per dtm value
        if not refs[0][1]: # refs are cartesian with a zone defined
            refs_proj=[(i[0],i[2]) for i in refs]
        else: # refs are lat/long
            if self.srs.startswith('+proj=latlong'):
                refs_proj=refs
            else: # reproject coordinates
#                ll = '\n'.join(['%r %r' % i[1] for i in refs])
#                refs_out=command(['gdaltransform','-s_srs','+proj=longlat','-t_srs',self.srs], ll)
#                coord_proj=[map(float,i.split()[:2]) for i in refs_out.splitlines()]

                srs_proj = osr.SpatialReference()
                srs_proj.ImportFromProj4(self.srs)
                srs_geo = osr.SpatialReference()
                srs_geo.ImportFromProj4('+proj=latlong')
                srs_geo.CopyGeogCSFrom(srs_proj)

                lonlat2proj=osr.CoordinateTransformation(srs_geo,srs_proj)
                lonlat=[i[1] for i in refs]
                coord_proj=lonlat2proj.TransformPoints(lonlat)
                ld('lonlat',lonlat,'\ncoord_proj',coord_proj)
                refs_proj=[(ref[0],coord) for ref,coord in zip(refs,coord_proj)]
        if len(refs) == 2:
            logging.warning(' Only 2 reference points: assuming the chart is north alligned')
            refs_proj.append(((refs_proj[0][0][0],refs_proj[1][0][1]),
                                (refs_proj[0][1][0],refs_proj[1][1][1])))
        ld('refs_proj',refs_proj)
        
        #double x = 0.0, double y = 0.0, double z = 0.0, double pixel = 0.0, 
        #double line = 0.0, char info = "", char id = ""
        gcps=[gdal.GCP(r[1][0],r[1][1],0,r[0][0],r[0][1], '','%d'%i)
                for r,i in zip(refs_proj,range(1,len(refs_proj)+1))]
                
#        transl_cmd=['gdal_translate','-of',format,img_path,out_dataset,'-a_srs', self.srs]
#        if self.options.expand:
#            transl_cmd=transl_cmd+['-expand',self.options.expand]

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
            dest_ds.SetGCPs(gcps,self.srs)
            dest_geotr=gdal.GCPsToGeoTransform(gcps)
            dest_ds.SetGeoTransform(dest_geotr)
#            poly,gmt_data=self.cut_poly(dest_ds,lonlat2proj)
            del dest_ds
        finally:
            os.chdir(cdir)

#        if self.options.get_cutline: # print cutline then return
#            print poly
#            return
#        if gmt_data and not self.options.no_cut_file: # create shapefile with a cut polygon
#            with open(base+'.gmt','w+') as f:
#                f.write(gmt_data)
# MapTranslator


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

from tiler_functions import *

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
        self.map_file=src_file.decode(locale.getpreferredencoding(),'ignore')

        self.header=self.get_header()       # Read map header
        self.img_file=self.get_raster()

        self.name=self.get_name()
        logging.info(' %s : %s (%s)' % (self.map_file,self.name,self.img_file))

        self.refs=self.get_refs()           # fetch reference points
        self.srs,self.dtm=self.get_srs()    # estimate SRS

    gmt_templ='''# @VGMT1.0 @GPOLYGON
# @Jp"%s"
# FEATURE_DATA
>
# @P
%s'''

    def cut_poly(self,out_dataset):
        plys=self.shift_lonlat(self.get_plys(),self.dtm)   # as per dtm value
        if not plys:
            return '',''

        # Create cutline
        lonlat=''.join(['%r %r\n' % i[1] for i in plys])
        if not plys[0][0]: # convert cutline coordinates to pixel xy using GDAL's navive srs for this raster
            pix_lines=command(['gdaltransform','-tps','-i','-t_srs','+proj=longlat',out_dataset],
                                lonlat).splitlines()
            pix_lst=[(int(i[0]),int(i[1])) for i in pix_lines]
        else:
            pix_lst=[i[0] for i in plys]
            pix_lines=['%d %d' % i for i in pix_lst]
        poly='POLYGON((%s))' % ','.join(pix_lines) # Create cutline

        size=self.get_size()
        if size:
            width,height=size
            inside=[i for i in pix_lst # check if the polygon is inside the image border
                if (i[0] > 0 or i[0] < width) or (i[1] > 0 or i[1] < height)]
            if not inside:
                return '',''

        # convert cutline geo coordinates to the chart's srs
        poly_xy=command(['gdaltransform','-tps','-s_srs','+proj=longlat','-t_srs',self.srs],lonlat)
        return poly,self.gmt_templ % (self.srs,poly_xy)

    def shift_lonlat(self,refs,dtm):
        if not dtm:
            return refs
        else:
            # alter refs as per DTM values
            split=zip(*refs) # split refs
            split[1]=[(lonlat[0]+dtm[0],lonlat[1]+dtm[1]) for lonlat in split[1]]
            return zip(*split) # repack refs

    def convert(self,dest=None):
        if dest:
            base=os.path.split(dest)[0]
        else:
            base=dest_path(self.map_file,self.options.dest_dir)

        dest_dir=os.path.split(base)[0]
        img_path=os.path.relpath(self.img_file,dest_dir)
        out_dataset= os.path.basename(base+'.vrt') # output VRT file    

        refs=self.shift_lonlat(self.refs,self.dtm)   # as per dtm value
        if not refs[0][1]: # refs are cartesian with a zone defined
            refs_proj=[(i[0],i[2]) for i in refs]
        else: # refs are lat/long
            if self.srs.startswith('+proj=latlong'):
                refs_proj=refs
            else: # reproject coordinates
                ll = '\n'.join(['%r %r' % i[1] for i in refs])
                refs_out=command(['gdaltransform','-s_srs','+proj=longlat','-t_srs',self.srs], ll)
                coord_proj=[map(float,i.split()[:2]) for i in refs_out.splitlines()]
                refs_proj=[(ref[0],coord) for ref,coord in zip(refs,coord_proj)]
        if len(refs) == 2:
            logging.warning(' Only 2 reference points: assuming the chart is north alligned')
            refs_proj.append(((refs_proj[0][0][0],refs_proj[1][0][1]),
                                (refs_proj[0][1][0],refs_proj[1][1][1])))
        ld('refs_proj',refs_proj)
        gcps=flatten([['-gcp']+map(repr, pix)+map(repr, coord) for pix,coord in refs_proj])
        transl_cmd=['gdal_translate','-of','VRT',img_path,out_dataset,'-a_srs', self.srs]
        if self.options.expand:
            transl_cmd=transl_cmd+['-expand',self.options.expand]
        try:
            cdir=os.getcwd()
            if dest_dir:
                os.chdir(dest_dir)
            command(transl_cmd + gcps)
            poly,gmt_data=self.cut_poly(out_dataset)
        finally:
            os.chdir(cdir)

        if self.options.get_cutline: # print cutline then return
            print poly
            return
        if gmt_data and not self.options.no_cut_file: # create shapefile with a cut polygon
            with open(base+'.gmt','w+') as f:
                f.write(gmt_data)
# MapTranslator


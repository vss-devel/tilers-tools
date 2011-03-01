#!/usr/bin/env python
# -*- coding: utf-8 -*-

# 2011-03-01 14:51:28 

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
###############################################################################

from __future__ import print_function
import sys
import os
import os.path
import logging
import shutil
from optparse import OptionParser
import math
from PIL import Image
import pickle
import mmap
import operator
import struct

try:
    from osgeo import gdal
    from osgeo import osr
    from osgeo import ogr
    from osgeo.gdalconst import *
#    gdal.TermProgress = gdal.TermProgress_nocb
except ImportError:
    import gdal
    import osr
    import ogr
    from gdalconst import *

from tiler_functions import *

#############################
#
# class TiledTiff
#
#############################

class TiledTiff(object):
    def __init__(self,img_fname,tile_first,tile_last):
        self.fname=img_fname
        self.tile_first=tile_first
        self.tile_last=tile_last
        self.file=open(img_fname,'r+b')
        self.mmap=mmap.mmap(self.file.fileno(),0,access=mmap.ACCESS_READ)

        self.read_hdr() # in child class
        self.IFD0,next_IFD=self.read_IFD(self.IFD0_ofs)

        assert self.tag_data('PlanarConfiguration')[0] == 2
        self.size=self.tag_data('ImageWidth')[0],self.tag_data('ImageLength')[0]
        self.tile_sz=self.tag_data('TileWidth')[0],self.tag_data('TileLength')[0]
        self.tile_range=map(lambda sz,tsz: (sz-1)//tsz+1,self.size,self.tile_sz)
        self.ntiles=self.tile_range[0]*self.tile_range[1]

        self.tile_ofs=self.tag_data('TileOffsets')
        self.tile_lengths=self.tag_data('TileByteCounts')

    def __del__(self):
        self.mmap.close()
        self.file.close()
        if options.verbose < 2:
            os.remove(self.fname)
                
    def tile(self,tile_x,tile_y):
        tile_sz=self.tile_sz
        ntiles=self.ntiles
        src=self.mmap
        ofs=self.tile_ofs
        lns=self.tile_lengths
        d_x,d_y=self.tile_first
        idx=tile_x-d_x+(tile_y-d_y)*self.size[0]//tile_sz[0]
        assert idx >= 0 and idx < ntiles, 'tile: %s range: %s' % ((tile_x,tile_y),self.tile_range)
        
        rgba=[ src[ofs[idx+i*ntiles] : ofs[idx+i*ntiles] + lns[idx+i*ntiles]] 
                for i in range(4)]
        if min(rgba[3]) == '\xFF': # fully opaque
            opacity=1
            mode='RGB'
            bands=rgba[:3]
        elif max(rgba[3]) == '\x00': # fully transparent
            return None,0
        else: # semi-transparent
            opacity=-1
            mode='RGBA'
            bands=rgba
        img=Image.merge(mode,[Image.frombuffer('L',tile_sz,b,'raw','L',0,1) for b in bands])
        return img,opacity

    def unpack(self,fmt,start,len,src=None):
        if not src: 
            src=self.mmap
        #ld(self.order+fmt,src[start:start+len])
        r=struct.unpack_from(self.order+fmt,src,start)
        #ld(r)
        return r

    def tag_in(self,name,tag_dict):
        return self.tag_map[name] in tag_dict

    tag_map={
        'ImageWidth':                   (256, 'LONG'),      # SHORT or LONG
        'ImageLength':                  (257, 'LONG'),      # SHORT or LONG
        'BitsPerSample':                (258, 'SHORT'),     # 4 or 8
        'Compression':                  (259, 'SHORT'),     # 1 or 32773
        'PhotometricInterpretation':    (262, 'SHORT'),     # 3
        'StripOffsets':                 (273, 'LONG'),
        'SamplesPerPixel':              (277, 'SHORT'),
        'RowsPerStrip':                 (278, 'LONG'),      # SHORT or LONG
        'StripByteCounts':              (279, 'LONG'),      # SHORT or LONG
        'XResolution':                  (282, 'RATIONAL'),
        'YResolution':                  (283, 'RATIONAL'),  
        'PlanarConfiguration':          (284, 'SHORT'),
        'ResolutionUnit':               (296, 'SHORT'),     # 1 or 2 or 3
        'ColorMap':                     (320, 'SHORT'),
        'TileWidth':                    (322, 'LONG'),      # SHORT or LONG
        'TileLength':                   (323, 'LONG'),      # SHORT or LONG
        'TileOffsets':                  (324, 'LONG'),
        'TileByteCounts':               (325, 'LONG'),      # SHORT or LONG
        'SampleFormat':                 (339, 'SHORT'),
        }
    tag_types={
        1: (1,'B'), # BYTE
        2: (1,'c'), # ASCII
        3: (2,'H'), # SHORT
        4: (4,'I'), # LONG
        5: (8,'II'), # RATIONAL
        # ...
        16:(8,"Q"), # TIFF_LONG8
        17:(8,"q"), # TIFF_SLONG8
        18:(8,"Q"), # TIFF_IFD8        
        }

    def tag_data(self,name):
        tag_id=self.tag_map[name][0]
        try:
            tag,datatype,data_cnt,data_ofs=self.IFD0[tag_id]
        except IndexError:
            raise Exception('Tiff tag not found: %s' % name)        
        type_len,code=self.tag_types[datatype]
        nbytes=type_len*data_cnt       
        if nbytes <= self.ptr_len:
            data_src=struct.pack(self.order+self.ptr_code,data_ofs)[0:nbytes]
            data_ofs=0
        else:
            data_src=self.mmap
        data=self.unpack(str(data_cnt)+code,data_ofs,nbytes,data_src)
        return data
                
    def read_IFD(self,IFD_ofs):
        tag_count=self.unpack(self.IFD_cnt_code,IFD_ofs,self.IFD_cnt_len)[0]
        IFD_ofs+=self.IFD_cnt_len
        tags_end=IFD_ofs+self.tag_len*tag_count
        next_IFD=self.unpack(self.ptr_code,tags_end,self.ptr_len)[0]
        tag_lst=[self.unpack('HH'+2*self.ptr_code,i,self.tag_len) 
            for i in range(IFD_ofs,tags_end,self.tag_len)] # tag,datatype,data_cnt,data_ofs
        tag_ids=zip(*tag_lst)[0] 
        tag_dict=dict(zip(tag_ids,tag_lst)) # tag: tag, datatype,data_cnt,data_ofs
        ld(tag_dict)
        return tag_dict, next_IFD

# TiledTiff

class TiffLE(TiledTiff):
    signature='II*\x00'
    order  = '<' # little endian
    tag_len= 12
    ptr_len = 4
    ptr_code='L'
    IFD_cnt_len=2
    IFD_cnt_code='H'
    
    def read_hdr(self):
        self.IFD0_ofs=self.unpack('2sH'+self.ptr_code,0,4+self.ptr_len)[-1]

class TiffBE(TiffLE):
    signature='MM\x00*'
    order  = '>' # big endian

class BigTiffLE(TiledTiff):
    signature='II+\x00'
    order  = '<' # little endian
    tag_len= 20
    ptr_len = 8
    ptr_code='Q'
    IFD_cnt_len=8
    IFD_cnt_code='Q'
    
    def read_hdr(self):
        endian,version,plen,resr,self.IFD0_ofs=self.unpack('2sHHH'+self.ptr_code,0,8+self.ptr_len)
        assert plen==self.ptr_len and resr==0

class BigTiffBE(BigTiffLE):
    signature='MM\x00+'
    order  = '>' # big endian

img_type_map=(
    TiffLE,
    TiffBE,
    BigTiffLE,
    BigTiffBE,
    )

def BaseImg(img_fname,tile_ul,tile_lr):
    f=open(img_fname,'r+b')
    magic=f.read(4)
    f.close()
    ld(magic)
    for img_cls in img_type_map:
        if magic == img_cls.signature:
            break
    else:
        raise Exception('Invalid base image')

    return img_cls(img_fname,tile_ul,tile_lr)

resampling_map={
    'near':     Image.NEAREST,
    'nearest':  Image.NEAREST,
    'bilinear': Image.BILINEAR,
    'bicubic':  Image.BICUBIC,
    'antialias':Image.ANTIALIAS,
    }
def resampling_lst(): return resampling_map.keys()
    
base_resampling_map={
    'near':         'NearestNeighbour', 
    'nearest':      'NearestNeighbour', 
    'bilinear':     'Bilinear',
    'cubic':        'Cubic',
    'cubicspline':  'CubicSpline',
    'lanczos':      'Lanczos',
    }
def base_resampling_lst(): return base_resampling_map.keys()
    
#############################
#
# templates for VRT XML
#
#############################

def xml_txt(elm_name,elm_value=None,elm_indent=0,**attr_dict):
    attr_txt=''.join((' %s="%s"' % (key,attr_dict[key]) for key in attr_dict))
    val_txt=('>%s</%s' % (elm_value,elm_name)) if elm_value else '/'
    return '%s<%s%s%s>' % (' '*elm_indent,elm_name,attr_txt,val_txt)

#    <WarpMemoryLimit>2.68435e+08</WarpMemoryLimit>
warp_vrt='''<VRTDataset rasterXSize="%(xsize)d" rasterYSize="%(ysize)d" subClass="VRTWarpedDataset">
  <SRS>%(srs)s</SRS>
%(geotr)s%(band_list)s
  <BlockXSize>%(blxsize)d</BlockXSize>
  <BlockYSize>%(blysize)d</BlockYSize>
  <GDALWarpOptions>
    <ResampleAlg>%(wo_ResampleAlg)s</ResampleAlg>
    <WorkingDataType>Byte</WorkingDataType>
    <SourceDataset relativeToVRT="0">%(wo_src_path)s</SourceDataset>
%(warp_options)s
    <Transformer>
      <ApproxTransformer>
        <MaxError>0.125</MaxError>
        <BaseTransformer>
          <GenImgProjTransformer>
%(wo_src_transform)s
%(wo_dst_transform)s
            <ReprojectTransformer>
              <ReprojectionTransformer>
                <SourceSRS>%(wo_src_srs)s</SourceSRS>
                <TargetSRS>%(wo_dst_srs)s</TargetSRS>
              </ReprojectionTransformer>
            </ReprojectTransformer>
          </GenImgProjTransformer>
        </BaseTransformer>
      </ApproxTransformer>
    </Transformer>
    <BandList>
%(wo_BandList)s
    </BandList>
%(wo_DstAlphaBand)s%(wo_Cutline)s  </GDALWarpOptions>
</VRTDataset>
'''
warp_band='  <VRTRasterBand dataType="Byte" band="%d" subClass="VRTWarpedRasterBand"%s>'
warp_band_color='>\n    <ColorInterp>%s</ColorInterp>\n  </VRTRasterBand'
warp_dst_alpha_band='    <DstAlphaBand>4</DstAlphaBand>\n'
warp_cutline='    <Cutline>%s</Cutline>\n'
warp_dst_geotr= '            <DstGeoTransform> %r, %r, %r, %r, %r, %r</DstGeoTransform>'
warp_dst_igeotr='            <DstInvGeoTransform> %r, %r, %r, %r, %r, %r</DstInvGeoTransform>'
warp_src_geotr= '            <SrcGeoTransform> %r, %r, %r, %r, %r, %r</SrcGeoTransform>'
warp_src_igeotr='            <SrcInvGeoTransform> %r, %r, %r, %r, %r, %r</SrcInvGeoTransform>'
warp_band_mapping='      <BandMapping src="%d" dst="%d"%s>'
warp_band_mapping_nodata='''>
        <SrcNoDataReal>%d</SrcNoDataReal>
        <SrcNoDataImag>%d</SrcNoDataImag>
      </BandMapping'''
warp_src_gcp_transformer='''            <SrcGCPTransformer>
              <GCPTransformer>
                <Order>%d</Order>
                <Reversed>0</Reversed>
                <GCPList>
%s
                </GCPList>
              </GCPTransformer>
            </SrcGCPTransformer>'''
warp_src_tps_transformer='''            <SrcTPSTransformer>
              <TPSTransformer>
                <Reversed>0</Reversed>
                <GCPList>
%s
                </GCPList>
              </TPSTransformer>
            </SrcTPSTransformer>'''

gcp_templ='    <GCP Id="%s" Pixel="%r" Line="%r" X="%r" Y="%r" Z="%r"/>'
gcplst_templ='  <GCPList Projection="%s">\n%s\n  </GCPList>\n'
geotr_templ='  <GeoTransform> %r, %r, %r, %r, %r, %r</GeoTransform>\n'
meta_templ='  <Metadata>\n%s\n  </Metadata>\n'
band_templ='''  <VRTRasterBand dataType="Byte" band="%(band)d">
    <ColorInterp>%(color)s</ColorInterp>
    <ComplexSource>
      <SourceFilename relativeToVRT="0">%(src)s</SourceFilename>
      <SourceBand>1</SourceBand>
      <SourceProperties RasterXSize="%(xsize)d" RasterYSize="%(ysize)d" DataType="Byte" BlockXSize="%(blxsize)d" BlockYSize="%(blysize)d"/>
      <SrcRect xOff="0" yOff="0" xSize="%(xsize)d" ySize="%(ysize)d"/>
      <DstRect xOff="0" yOff="0" xSize="%(xsize)d" ySize="%(ysize)d"/>
      <ColorTableComponent>%(band)d</ColorTableComponent>
    </ComplexSource>
  </VRTRasterBand>
'''
srs_templ='  <SRS>%s</SRS>\n'
vrt_templ='''<VRTDataset rasterXSize="%(xsize)d" rasterYSize="%(ysize)d">
%(metadata)s%(srs)s%(geotr)s%(gcp_list)s%(band_list)s</VRTDataset>
'''
    
#############################
#
# class Pyramid
#
#############################

class Pyramid(object):
    def __init__(self,src,dest,options):
        gdal.UseExceptions()
        self.src=src
        self.dest=dest
        self.options=options
        self.src_path=self.src
        self.tiles_prefix=options.tiles_prefix
        self.init_proj()
        self.tile_ext='.'+options.tile_format.lower()
        self.tile_sz=tuple(map(int,options.tile_size.split(',')))
        self.src_dir,src_f=os.path.split(src)
        self.base=os.path.splitext(src_f)[0]
        if os.path.isdir(self.dest):
            shutil.rmtree(self.dest,ignore_errors=True)
        os.makedirs(self.dest)
        self.base_resampling=base_resampling_map[options.base_resampling]
        self.resampling=resampling_map[options.overview_resampling]
        self.init_map(options.zoom)

    #############################
    #
    # initialize geo-parameters and generate base zoom level
    #
    #############################

    def init_map(self,zoom_parm):
        options=self.options
        src_vrt=os.path.join(self.dest,self.base+'.src.vrt') # auxilary VRT file
        temp_vrt=os.path.join(self.dest,self.base+'.tmp.vrt') # auxilary VRT file
        self.src_ds=self.get_src_ds(src_vrt)

        # calculate zoom range
        pf('\n%s -> %s '%(self.src,self.dest),end='')
        self.zoom_range=self.calc_zoom(zoom_parm,temp_vrt)
             
        # reproject to base zoom
        zoom=self.zoom_range[0]
        #
        # extend the raster to cover full tiles
        ld("extend the raster")
        shifted_srs=self.shift_srs(zoom)

        # get origin and corners at the target SRS
        t_ds=gdal.AutoCreateWarpedVRT(self.src_ds,None,proj4wkt(shifted_srs))
        tr=MyTransformer(t_ds)
        self.origin,self.extent=tr.transform([(0,0),(t_ds.RasterXSize,t_ds.RasterYSize)])

        tr_srs=MyTransformer(SRC_SRS=shifted_srs,DST_SRS=self.proj)
        shift_x,y,z=tr_srs.transform_pt((0,0))
        ld('new_srs',shifted_srs,'shift_x',shift_x,'coord_offset',self.coord_offset)
        if shift_x < 0:
            shift_x+=self.zoom0_tile_dim[0]*2
        self.coord_offset[0]+=shift_x
        self.shift_x=shift_x
        self.proj=shifted_srs

        # adjust raster extents to tile boundaries
        tile_ul,tile_lr=self.corner_tiles(zoom)
        ld('tile_ul',tile_ul,'tile_lr',tile_lr)
        ul_c=self.tile2coord_box(tile_ul)[0]
        lr_c=self.tile2coord_box(tile_lr)[1]
        ul_pix=self.tile_corners(tile_ul)[0]
        lr_pix=self.tile_corners(tile_lr)[1]
        xsize=lr_pix[0]-ul_pix[0]
        ysize=lr_pix[1]-ul_pix[1]
        res=self.zoom2res(zoom)
        geotr=( ul_c[0], res[0],     0.0,
                ul_c[1],    0.0, -res[1] )
        ok,igeotr=gdal.InvGeoTransform(geotr)
        assert ok
        
        ld('Upper Left ',self.origin,ul_c)
        ld('Lower Right',self.extent,lr_c)
        ld('coord_offset',self.coord_offset,'ul_c+c_off',map(operator.add,ul_c,self.coord_offset))

        src_proj=wkt2proj4(self.src_ds.GetProjection())
        src_geotr=self.src_ds.GetGeoTransform()
        gcp_proj=None
        if src_geotr and src_geotr != (0.0, 1.0, 0.0, 0.0, 0.0, 1.0):
            ok,src_igeotr=gdal.InvGeoTransform(src_geotr)
            assert ok
            src_transform='%s\n%s' % (warp_src_geotr % src_geotr,warp_src_igeotr % src_igeotr)
        else:
            gcps=self.src_ds.GetGCPs()
            assert gcps, ' Neither geotransform, nor gpcs are in the source file %s' % self.src

            gcp_lst=[(g.Id,g.GCPPixel,g.GCPLine,g.GCPX,g.GCPY,g.GCPZ) for g in gcps]
            ld('src_proj',self.src_ds.GetProjection())
            ld('gcp_proj',self.src_ds.GetGCPProjection())
            gcp_proj=wkt2proj4(self.src_ds.GetGCPProjection())
            if src_proj and gcp_proj != src_proj:
                gcp2src=MyTransformer(SRC_SRS=gcp_proj,DST_SRS=srs_proj)
                coords=gcp2src.transform([g[3:6] for g in gcp_lst])
                gcp_lst=[tuple(p[:3]+c) for p,c in zip(gcp_lst,coords)]

            gcp_txt='\n'.join((gcp_templ % g for g in gcp_lst))
            #src_transform=warp_src_gcp_transformer % (0,gcp_txt)
            src_transform=warp_src_tps_transformer % gcp_txt

        dst_transform='%s\n%s' % (warp_dst_geotr % geotr,warp_dst_igeotr % igeotr)

        def w_option(name,value):
            return '    <Option name="%s">%s</Option>' % (name,value)
            
        warp_options=[]
        warp_options.append(w_option('INIT_DEST','0'))

        if options.cut:
            cut_wkt=self.cut_line()
        else:
            cut_wkt=None
        if cut_wkt:
            warp_options.append(w_option('CUTLINE',cut_wkt))
            if options.blend_dist:
                warp_options.append(w_option('CUTLINE_BLEND_DIST',options.blend_dist))

        nodata=None
        if options.no_data:
            nodata=options.no_data.split(',')
            warp_options.append(w_option('UNIFIED_SRC_NODATA','YES'))

        src_bands=self.src_ds.RasterCount
        assert src_bands in (3,4)
        vrt_bands=[]
        wo_BandList=[]
        for i in range(1,src_bands+1):
            vrt_bands.append(warp_band % (i,'/'))
            wo_BandList.append(warp_band_mapping % (i,i,
                (warp_band_mapping_nodata % (nodata[i],0)) if nodata else '/'))

        if src_bands < 4:
            vrt_bands.append(warp_band % (i+1,warp_band_color % 'Alpha'))

        block_sz=self.src_ds.GetRasterBand(1).GetBlockSize()

        vrt_text=warp_vrt % {
            'xsize':            xsize,
            'ysize':            ysize,
            'srs':              self.proj,
            'geotr':            geotr_templ % geotr,
            'band_list':        '\n'.join(vrt_bands),
            'blxsize':          self.tile_sz[0],
            'blysize':          self.tile_sz[1],
            'wo_ResampleAlg':   self.base_resampling,
            'wo_src_path':      self.src_path,
            'warp_options':     '\n'.join(warp_options),
            'wo_src_srs':       gcp_proj if gcp_proj else src_proj,
            'wo_dst_srs':       self.proj,
            'wo_src_transform': src_transform,
            'wo_dst_transform': dst_transform,
            'wo_BandList':      '\n'.join(wo_BandList),
            'wo_DstAlphaBand':  warp_dst_alpha_band if src_bands < 4 else '',
            'wo_Cutline':       (warp_cutline % cut_wkt) if cut_wkt else '',
            }

        with open(temp_vrt,'w') as f:
            f.write(vrt_text)

        t_ds = gdal.Open(vrt_text,GA_ReadOnly)
        # create Gtiff
        pf('...',end='')
        dst_drv = gdal.GetDriverByName('Gtiff')
        temp_tif=os.path.join(self.dest,self.base+'.tmp_%i.tiff' % zoom) # img for the base zoom
            
        dst_ds = dst_drv.CreateCopy(temp_tif,t_ds,0,[
			'TILED=YES',
            'INTERLEAVE=BAND',
            'BLOCKXSIZE=%i' % self.tile_sz[0],
            'BLOCKYSIZE=%i' % self.tile_sz[1],
            ])#, gdal.TermProgress)
        del dst_ds
        del t_ds
        del self.src_ds
        pf('.',end='')

        # create base_image raster
        self.base_img=BaseImg(temp_tif,tile_ul[1:],tile_lr[1:])

        if options.verbose < 2:
            try:
                os.remove(temp_vrt)
                os.remove(src_vrt)
            except: pass

    def get_src_ds(self,src_vrt):
        'get src dataset, convert to RGB(A) if required'
        if os.path.exists(self.src):
            self.src_path=os.path.abspath(self.src)

        # check for source raster type
        src_ds=gdal.Open(self.src_path,GA_ReadOnly)
        src_bands=src_ds.RasterCount
        if src_bands >= 3:
            return src_ds

        # convert to rgb VRT
        assert src_bands == 1
        xsize,ysize=(src_ds.RasterXSize,src_ds.RasterYSize)

        srs_proj=wkt2proj4(src_ds.GetProjection())
        src_geotr=src_ds.GetGeoTransform()
        ld('src geotr',src_geotr)
        if not src_geotr or src_geotr == (0.0, 1.0, 0.0, 0.0, 0.0, 1.0):
            geotr_txt=''
        else:
            geotr_txt=geotr_templ % src_geotr

        gcplst_txt=''
        gcps=src_ds.GetGCPs()
        if gcps:
            gcp_lst='\n'.join((gcp_templ % (g.Id,g.GCPPixel,g.GCPLine,g.GCPX,g.GCPY,g.GCPZ) 
                                for g in gcps))
            gcp_proj=wkt2proj4(src_ds.GetGCPProjection())
            gcplst_txt=gcplst_templ % (gcp_proj,gcp_lst)

        metadt=src_ds.GetMetadata()
        if metadt:
            mtd_lst=[xml_txt('MDI',metadt[mdkey],4,key=mdkey) for mdkey in metadt]
            meta_txt=meta_templ % '\n'.join(mtd_lst)
        else:
            meta_txt=''

        block_sz=src_ds.GetRasterBand(1).GetBlockSize()
        band_lst=''.join((band_templ % {
            'band':     band,
            'color':    color,
            'src':      self.src_path,
            'xsize':    xsize,
            'ysize':    ysize,
            'blxsize':    block_sz[0],
            'blysize':    block_sz[1],
            } for band,color in ((1,'Red'),(2,'Green'),(3,'Blue'))))
        vrt_txt=vrt_templ % {
            'xsize':    xsize,
            'ysize':    ysize,
            'metadata': meta_txt,
            'srs':      (srs_templ % srs_proj) if srs_proj else '',
            'geotr':    geotr_txt,
            'gcp_list': gcplst_txt,
            'band_list':band_lst,
            }
        with open(src_vrt,'w') as f:
            f.write(vrt_txt)
        self.src_path=src_vrt

        src_ds=gdal.Open(vrt_txt,GA_ReadOnly)
        return src_ds

    def shift_srs(self,zoom=None):
        'change prime meridian to allow charts crossing 180 meridian'
        t_ds=self.src_ds
#        tr_pix_ll=gdal.Transformer(t_ds,None,['METHOD=GCP_TPS','DST_SRS='+proj4wkt(self.latlong)])
        tr_pix_ll=MyTransformer(t_ds,DST_SRS=self.latlong)
        ll=tr_pix_ll.transform([(0,0,0),(t_ds.RasterXSize,t_ds.RasterYSize,0)])
        ul,lr=ll
        ld('ul',ul,'lr',lr)
        if ul[0] < lr[0]:
            return self.proj

        left_lon=int(math.floor(ul[0]))
        tr_srs=MyTransformer(SRC_SRS=self.latlong,DST_SRS=self.proj)
        left_xy=tr_srs.transform_pt((left_lon,0))
        if zoom is not None: # adjust to a tile boundary
            left_xy=self.tile2coord_box(self.coord2tile(zoom,left_xy))[0]

#        tr_srs=get_transformer(SRC_SRS=self.proj,DST_SRS=self.latlong)
        new_pm=int(math.floor(
                tr_srs.transform_pt(left_xy,inv=True)[0]
                ))#-180
        ld('left_xy',left_xy,'new_pm',new_pm)
        return '%s +lon_0=%d' % (self.proj,new_pm)

    def calc_zoom(self,zoom_parm,temp_vrt):
        'determite zoom levels to generate'
        res=None
        if not zoom_parm: # calculate "automatic" zoom levels
            # check raster parameters to find default zoom range
            # modify target srs to allow charts crossing meridian 180
            ld('"automatic" zoom levels')
            shifted_srs=self.shift_srs()

            t_ds=gdal.AutoCreateWarpedVRT(self.src_ds,None,proj4wkt(shifted_srs))
            geotr=t_ds.GetGeoTransform()
            res=(geotr[1], -geotr[5])
            max_zoom=max(self.res2zoom_xy(res))

            # calculate min_zoom
            ul_c=(geotr[0], geotr[3])
            lr_c=gdal.ApplyGeoTransform(geotr,t_ds.RasterXSize,t_ds.RasterYSize)
            wh=(lr_c[0]-ul_c[0],ul_c[1]-lr_c[1])
            ld(ul_c,lr_c,wh)
            min_zoom=min(self.res2zoom_xy([wh[i]/self.tile_sz[i]for i in (0,1)]))
            zoom_parm='%d-%d'%(min_zoom,max_zoom)
        zchunks=[map(int,z.split('-')) for z in zoom_parm.split(',')]
        zrange=[]
        for z in zchunks:
            if len(z) == 1:
                zrange+=z
            else:
                zrange+=range(min(z),max(z)+1)
        zoom_range=list(reversed(sorted(set(zrange))))
        ld(('res',res,'zoom_range',zoom_range,'z0 (0,0)',self.coord2pix(0,(0,0))))
        return zoom_range

    def cut_line(self):
        src_ds=self.src_ds
        src_proj=wkt2proj4(src_ds.GetProjection())
        cutline=src_ds.GetMetadataItem('CUTLINE')
        ld('cutline',cutline)
        if options.cutline:
            cut_file=self.options.cutline
        else: # try to find a file with a cut shape
            if cutline:
                return cutline
            for ext in ('.gmt','.shp'):
                cut_file=os.path.join(self.src_dir,self.base+ext)
                if os.path.exists(cut_file):
                    break
            else:
                cut_file=None
        if cut_file:
            cut_ds=ogr.Open(cut_file)
            assert cut_ds, ' Invalid cutline file %s' % cut_file
            layer=cut_ds.GetLayer()
            poly_srs=layer.GetSpatialRef()
            if poly_srs:
                poly_proj=poly_srs.ExportToProj4()
            else:
                poly_proj=src_proj
            feature=layer.GetFeature(0)
            geom=feature.GetGeometryRef()
            gname=geom.GetGeometryName()
            poly_iter={
                'MULTIPOLYGON':(geom.GetGeometryRef(i) for i in range(geom.GetGeometryCount())),
                'POLYGON': (geom,),
                }[gname]
            mpoly=[]
            srs_tr=MyTransformer(SRC_SRS=poly_proj,DST_SRS=src_proj)
            if poly_proj == src_proj:
                srs_tr.transform=lambda x:x
            pix_tr=MyTransformer(src_ds)
            for pl in poly_iter:
                assert pl.GetGeometryName() == 'POLYGON'
                for ln in (pl.GetGeometryRef(j) for j in range(pl.GetGeometryCount())):
                    assert ln.GetGeometryName() == 'LINEARRING'
                    points=[ln.GetPoint(n) for n in range(ln.GetPointCount())]
                    p_src=srs_tr.transform(points)
                    p_pix=pix_tr.transform_pt(p_src,inv=True)
                    mpoly.append(','.join(['%r %r' % (p[0],p[1]) for p in p_pix]))
            cutline='MULTIPOLYGON(%s)' % ','.join(['((%s))' % poly for poly in mpoly])
        return cutline

    #############################
    #
    # utility finctions
    #
    #############################

    def zoom2res(self,zoom):
        return map(lambda z0,ts:
            z0/(ts*2**zoom), self.zoom0_tile_dim,self.tile_sz)

    def res2zoom_xy(self,res):
        'resolution to zoom levels (separate for x and y)'
        z=map(lambda z0,ts,r:
            int(math.floor(math.log(z0/(ts*r),2))), self.zoom0_tile_dim,self.tile_sz,res)
        return [v if v>0 else 0 for v in z]

    def pix2tile(self,zoom,pix):
        'pixel coordinates to tile (z,x,y)'
        zoom_tiles=map(lambda v: v*2**zoom,self.zoom0_tiles)
        return [zoom]+map(lambda p,ts:p/ts,pix,self.tile_sz)

    def tile2pix(self,tile):
        'tile to pixel coordinates'
        return map(operator.mul,self.tile_sz,tile[1:])

    def coord2tile(self,zoom,coord):
        'cartesian coordinates to tile numbers'
        return self.pix2tile(zoom,self.coord2pix(zoom,coord))

    def tile_corners(self,tile):
        'pixel coordinates of a tile'
        z,x,y=tile
        return [self.tile2pix((z,x,y)),self.tile2pix((z,x+1,y+1))]

    def tile2coord_box(self,tile):
        'cartesian coordinates of tile corners'
        z=tile[0]
        return map(self.pix2coord,(z,z),self.tile_corners(tile))

    def coord2pix(self,zoom,coord):
        'cartesian coordinates to pixel coordinates'
        ul_c=(coord[0]+self.coord_offset[0],self.coord_offset[1]-coord[1])
        out=tuple([int(c/r) for c,r in zip(ul_c,self.zoom2res(zoom))])
        return out

    def pix2coord(self,zoom,pix_coord):
        pix00_ofs=map(operator.mul,pix_coord,self.zoom2res(zoom))
        return (pix00_ofs[0]-self.coord_offset[0],-(pix00_ofs[1]-self.coord_offset[1]))

    def tile_path(self,tile):
        'relative path to a tile'
        z,x,y=self.tile_norm(tile)
        return '%i/%i/%i%s' % (z,x,y,self.tile_ext)

    def tile_norm(self,tile):
        'x,y of a tile do not exceed max values'
        z,x,y=tile
        zoom_tiles=map(lambda v: v*2**z,self.zoom0_tiles)
        return [z,x%zoom_tiles[0],y%zoom_tiles[1]]

    def coords2latlong(self, coords): # redefined in PlateCarree
        try:
            self.proj2latlong
        except: # first time call
            self.proj2latlong=MyTransformer(SRC_SRS=self.proj,DST_SRS=self.latlong)
        return [i[:2] for i in self.proj2latlong.transform(coords)]

    def boxes2latlong(self,box_lst):
        out=self.coords2latlong(flatten(box_lst))
        deg_lst=[[(j+180)%360-180 for j in i] for i in out]
        ul_lst=deg_lst[0::2]
        lr_lst=[((x if x != -180 else 180),y) for x,y in deg_lst[1::2]]
        return zip(ul_lst,lr_lst)

    def corner_tiles(self,zoom):
        p_ul=self.coord2pix(zoom,self.origin)
        t_ul=[zoom]+map(lambda p,ts:(p+1)//ts,p_ul,self.tile_sz)

        p_lr=self.coord2pix(zoom,self.extent)
        t_lr=[zoom]+map(lambda p,ts:(p-1)//ts,p_lr,self.tile_sz)

        ld('zoom',zoom,'p_ul',p_ul,'p_lr',p_lr,'t_ul',t_ul,'t_lr',t_lr,
            'zoom tiles',map(lambda zt:zt*2**zoom,self.zoom0_tiles),
            'zoom size',map(lambda zt,ts:ts*zt*2**zoom,self.zoom0_tiles,self.tile_sz))
        return t_ul,t_lr

    def info(self,info_lst,pattern):
        info_line=[i for i in info_lst if pattern in i][0]
        if '(' in info_line:
            out=info_line[info_line.find('(')+1:info_line.find(')')]
        else:
            out=info_line[info_line.find(pattern)+len(pattern):]
        return out.strip()
        
    tick_rate=50
    count=0

    def counter(self):
        self.count+=1
        if self.count % self.tick_rate == 0:
            pf('.',end='')
            return True
        else:
            return False

    #############################
    #
    # generate pyramid
    #
    #############################

    def walk_pyramid(self):
        self.make_googlemaps()
        tiles=[]
        for zoom in self.zoom_range:
            tile_ul,tile_lr=self.corner_tiles(zoom)
            zoom_tiles=flatten([[(zoom,x,y) for x in range(tile_ul[1],tile_lr[1]+1)] 
                                           for y in range(tile_ul[2],tile_lr[2]+1)])
            tiles.extend(zoom_tiles)
        ld('min_zoom',zoom,'tile_ul',tile_ul,'tile_lr',tile_lr,'zoom_tiles',zoom_tiles)
        self.all_tiles=set(tiles)
        top_tiles=filter(None,map(self.proc_tile,zoom_tiles))
        # write top kml
        self.make_kml(None,[ch for img,ch,opacities in top_tiles])
        pf('')
        
        # cache back tiles opacity
        file_opacities=[(self.tile_path(tile),opc)
            for tile,opc in flatten([opacities for img,ch,opacities in top_tiles])]
        try:
            pickle.dump(dict(file_opacities),open(os.path.join(self.dest, 'merge-cache'),'w'))
        except:
            logging.warning("opacity cache save failed")

    def proc_tile(self,tile):
        #ld(tile)
        ch_opacities=[]
        ch_tiles=[]
        zoom,x,y=tile
        if zoom==self.zoom_range[0]: # crop from the base image
            tile_img,opacity=self.base_img.tile(*tile[1:])
        else: # mosaic from children
            opacity=0
            cz=self.zoom_range[self.zoom_range.index(zoom)-1] # child's zoom
            dz=int(2**(cz-zoom))
            children=dict(flatten([[((cz,x*dz+dx,y*dz+dy),(dx*self.tile_sz[0]//dz,dy*self.tile_sz[1]//dz))
                               for dx in range(dz)]
                                   for dy in range(dz)]))
            ch_tiles=filter(None,map(self.proc_tile,self.all_tiles & set(children)))
            if ch_tiles:
                if len(ch_tiles) == 4 and all([opacities[0][1]==1 for img,ch,opacities in ch_tiles]):
                    mode="RGB"
                    opacity=1
                else:
                    mode="RGBA"
                    opacity=-1
                tile_img=Image.new(mode,self.tile_sz,(0,0,0,0))
                for img,ch,opacities in ch_tiles:
                    #ld((tile,dz,ch,children[ch]))
                    img=img.resize([i//dz for i in img.size],self.resampling)
                    mask=img if img.mode =='RGBA' else None
                    tile_img.paste(img,children[ch],mask)
                    ch_opacities+=opacities
        if opacity:
            self.write_tile(tile,tile_img)
            self.make_kml(tile,[ch for img,ch,opacities in ch_tiles])
            return tile_img,tile,[(tile,opacity)]+ch_opacities

    def write_tile(self,tile,img):
        rel_path=self.tile_path(tile)
        full_path=os.path.join(self.dest,rel_path)
        try:
            os.makedirs(os.path.dirname(full_path))
        except: pass
        if options.to_palette and self.tile_ext == '.png':
            try:
                img=img.convert('P', palette=Image.ADAPTIVE, colors=255)
            except ValueError:
                #ld('img.mode',img.mode)
                pass
        img.save(full_path)
        self.counter()

    def make_kml(self,tile,children=[]): # 'virtual'
        pass

    def make_googlemaps(self): # 'virtual'
        pass        

# Pyramid        

class PlateCarree(Pyramid):
    'Google Earth (plate carrée), Google tile numbering'
    format='earth'
    defaul_ext='.earth'

    # http://earth.google.com/support/bin/static.py?page=guide.cs&guide=22373&topic=23750
    # "Google Earth uses Simple Cylindrical projection for its imagery base. This is a simple map 
    # projection where the meridians and parallels are equidistant, straight lines, with the two sets 
    # crossing at right angles. This projection is also known as Lat/Lon WGS84"    

    proj='+proj=eqc +datum=WGS84 +ellps=WGS84' # Equirectangular (aka plate carrée, aka Simple Cylindrical)
    latlong='+proj=latlong +datum=WGS84'

    def init_proj(self):
        tr=MyTransformer(SRC_SRS=self.latlong,DST_SRS=self.proj)
        semi_circ,semi_meridian,z=tr.transform_pt((180,90))
        self.zoom0_tiles=[2,1] # tiles at zoom 0
        self.zoom0_tile_dim=[semi_circ,semi_meridian*2] # dimentions of a tile at zoom 0
        self.coord_offset=[semi_circ,semi_meridian]

    def coords2latlong(self, coords):
        out=[map(lambda v,dim,c_off: (v+c_off)/dim*180,
                coord,self.zoom0_tile_dim,(self.shift_x,0)) for coord in coords]
        return [(lon,lat) for lon,lat in out]

    def kml_child_links(self,children,parent=None,path_prefix=''):
        kml_links=[]
        # convert tiles to degree boxes
        latlong_lst=self.boxes2latlong([self.tile2coord_box(t) for t in children])
        
        for tile,latlong in zip(children,latlong_lst):
            #ld(tile,latlong)
            w,n,e,s=['%.11f'%v for v in flatten(latlong)]
            name=os.path.splitext(self.tile_path(tile))[0]
            # fill in kml link template
            kml_links.append( kml_link_templ % { 
                'name':    name,
                'href':    path_prefix+'%s.kml' % name,
                'west':    w, 'north':    n,
                'east':    e, 'south':    s,
                'min_lod': 128,
                'max_lod': 2048 if parent else -1,
                })
        return ''.join(kml_links)

    def write_kml(self,rel_path,name,links='',overlay=''):
        kml= kml_templ % {
            'name':      name,
            'links':     links,
            'overlay':   overlay,
            'dbg_start': '' if options.verbose < 2 else '    <!--\n',
            'dbg_end':   '' if options.verbose < 2 else '      -->\n',
            }
        open(os.path.join(self.dest,rel_path+'.kml'),'w+').write(kml)

    def make_kml(self,tile,children=[]): #
        if not tile: # create top level kml
            self.write_kml(os.path.basename(self.base),os.path.basename(self.base),self.kml_child_links(children))
            return
        # fill in kml templates
        rel_path=self.tile_path(tile)
        name=os.path.splitext(rel_path)[0]
        kml_links=self.kml_child_links(children,tile,'../../')
        tile_box=self.boxes2latlong([self.tile2coord_box(tile)])[0]
        w,n,e,s=['%.11f'%v for v in flatten(tile_box)]
        kml_overlay = kml_overlay_templ % {
            'name':    name,
            'href':    os.path.basename(rel_path),
            'min_lod': 128,
            'max_lod': 2048 if kml_links else -1,
            'order':   tile[0],
            'west':    w, 'north':    n,
            'east':    e, 'south':    s,
            }
        self.write_kml(name,name,kml_links,kml_overlay)
            
# PlateCarree

class PlateCarreeTMS(PlateCarree):
    'Google Earth (plate carrée), TMS tile numbering'
    format='earth-tms'
    defaul_ext='.earth-tms'

    def tile_path(self,tile):
        z,x,y=self.tile_norm(tile)
        y_tiles=self.zoom0_tiles[1]*2**z
        return '%i/%i/%i%s' % (z,x,y_tiles-y-1,self.tile_ext)

kml_templ='''<?xml version="1.0" encoding="utf-8"?>
<kml xmlns="http://earth.google.com/kml/2.1">
    <Document>
    <!-- Generated by gdal_tiler.py (http://code.google.com/p/tilers-tools/) -->
%(dbg_start)s        <Style> 
            <ListStyle id="hideChildren"> <listItemType>checkHideChildren</listItemType> </ListStyle>
        </Style>
%(dbg_end)s        <name>%(name)s</name>%(overlay)s%(links)s
    </Document>
</kml>
'''

kml_overlay_templ='''
        <Region> 
            <Lod> 
                <minLodPixels>%(min_lod)s</minLodPixels> 
                <maxLodPixels>%(max_lod)s</maxLodPixels>
            </Lod>
            <LatLonAltBox>
            	<west>%(west)s</west> <north>%(north)s</north>
            	<east>%(east)s</east> <south>%(south)s</south>
            </LatLonAltBox>
        </Region>
        <GroundOverlay>
            <name>%(name)s</name>
            <drawOrder>%(order)s</drawOrder>
            <Icon> <href>%(href)s</href> </Icon>
            <LatLonBox>
                <west>%(west)s</west> <north>%(north)s</north>
                <east>%(east)s</east> <south>%(south)s</south>
            </LatLonBox>
        </GroundOverlay>'''

kml_link_templ='''
        <NetworkLink>
            <name>%(name)s</name>
            <Region> 
                <Lod> 
                    <minLodPixels>%(min_lod)s</minLodPixels> 
                    <maxLodPixels>%(max_lod)s</maxLodPixels>
                </Lod>
                <LatLonAltBox>
                    <west>%(west)s</west> <north>%(north)s</north>
                    <east>%(east)s</east> <south>%(south)s</south>
                </LatLonAltBox>
            </Region>
            <Link> <viewRefreshMode>onRegion</viewRefreshMode>
                <href>%(href)s</href>
            </Link>
        </NetworkLink>'''

class Gmaps(Pyramid):
    'Google Maps (Global Mercator), native tile numbering'
    format='gmaps'
    defaul_ext='.gmaps'

    #proj='epsg:900913' # Google Maps Global Mercator
    #proj='epsg:3857' # Google Maps Global Mercator
    proj='+proj=merc +a=6378137 +b=6378137 +nadgrids=@null +wktext' # Google Maps Global Mercator
    latlong='+proj=latlong +a=6378137 +b=6378137  +nadgrids=@null +wktext'
    
    def init_proj(self):
        # Half Equator length in meters
        tr=MyTransformer(SRC_SRS=self.latlong,DST_SRS=self.proj)
        semi_circ,y,z=tr.transform_pt((180,0))
        ld(semi_circ)
        self.zoom0_tiles=[1,1] # tiles at zoom 0
        self.zoom0_tile_dim=[semi_circ*2,semi_circ*2]  # dimentions of a tile at zoom 0
        self.coord_offset=[semi_circ,semi_circ]

    def make_googlemaps(self):
        ul,lr=self.boxes2latlong([(self.origin,self.extent)])[0]
        googlemaps = google_templ % dict(
            title=      os.path.basename(self.dest),
            longlat_ll= '%s, %s' % (lr[1],ul[0]),
            longlat_ur= '%s, %s' % (ul[1],lr[0]),
            minzoom=    self.zoom_range[-1],
            maxzoom=    self.zoom_range[0],
            header=     os.path.basename(self.dest), 
            tms_tiles=  'true' if self.defaul_ext == '.tms' else 'false',
            map_type=   'ROADMAP',
            tile_ext=   self.tile_ext,
            tile_size=  '%s, %s' % self.tile_sz,
            tiles_root= self.tiles_prefix,
            )
        open(os.path.join(self.dest,'gmaps.html'),'w').write(googlemaps)

# GMaps

class GMercatorTMS(Gmaps):
    'Google Maps (Global Mercator), TMS tile numbering'
    defaul_ext='.tms'
    format='tms'

    def tile_path(self,tile):
        z,x,y=self.tile_norm(tile)
        y_tiles=self.zoom0_tiles[1]*2**z
        return '%i/%i/%i%s' % (z,x,y_tiles-y-1,self.tile_ext)
    
google_templ='''<!DOCTYPE html>

<!-- Generated by gdal_tiler.py (http://code.google.com/p/tilers-tools/) -->

<html>
<head>
<title>%(title)s</title>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8">
<!-- <meta name="viewport" content="initial-scale=1.0, user-scalable=no" /> -->
<style type="text/css">
  html { height: 100%% }
  body { height: 100%%; margin: 0px; padding: 0px }
  #map_canvas { height: 100%% }
</style>
<script type="text/javascript"
    src="http://maps.google.com/maps/api/js?sensor=false">
</script>
<script type="text/javascript">
    var G=google.maps; // use G instead of google.maps

    var mapBounds = new G.LatLngBounds(new G.LatLng(%(longlat_ll)s), new G.LatLng(%(longlat_ur)s));
    var mapMinZoom = %(minzoom)d;
    var mapMaxZoom = %(maxzoom)d;
    var tile_ext = "%(tile_ext)s";
    var tile_size = new G.Size(%(tile_size)s);
    var tiles_root = "%(tiles_root)s";
    var tms_tiles = %(tms_tiles)s;
    var map_type = G.MapTypeId.%(map_type)s;

    var opacity = 0.5;

    function log(msg) {
        setTimeout(function() {
            throw new Error(msg);
        }, 0);
    }

    function map_overlay(){
        return new G.ImageMapType({
            getTileUrl: function(coord, zoom) {
                max_x=1<<zoom
                max_y=1<<zoom
                x=coord.x %% max_x;
                if (x < 0)
                    x=max_x+x
                y=coord.y;
                if (tms_tiles) y=(1<<zoom)-coord.y-1;
                //log(tiles_root+zoom+"/"+x+"/"+y+tile_ext)
                return tiles_root+zoom+"/"+x+"/"+y+tile_ext;
                },
            tileSize: tile_size,
            opacity: opacity,
            isPng: (tile_ext == ".png")
            })
        }

    function opacity_str(opacity){
        s=String(Math.round(opacity*100));
        while (s.length < 3) s='+'+s;
        return '<+++'+s+'%%+++>';
        }

    function opacity_control(map,overlay_index) {
        controlDiv=document.createElement('DIV');
        // Set CSS styles for the DIV containing the control
        // Setting padding to 5 px will offset the control
        // from the edge of the map
        controlDiv.style.padding = '5px';
        controlDiv.id = 'op-control-div';

        // Set CSS for the control border
        var controlUI = document.createElement('DIV');
        controlUI.style.backgroundColor = 'white';
        controlUI.style.borderStyle = 'solid';
        controlUI.style.borderWidth = '2px';
        controlUI.style.cursor = 'pointer';
        controlUI.style.textAlign = 'center';
        controlUI.title = 'Click to set opacity of the overlay';
        controlUI.id = 'op-control';
        controlDiv.appendChild(controlUI);

        // Set CSS for the control interior
        var controlText = document.createElement('DIV');
        controlText.style.fontFamily = 'Arial,sans-serif';
        controlText.style.fontile_sz = '10px';
        controlText.style.paddingLeft = '4px';
        controlText.style.paddingRight = '4px';
        controlText.innerHTML = opacity_str(opacity);
        controlText.id = 'op-control-txt';
        controlUI.appendChild(controlText);

        G.event.addDomListener(controlText, 'click', function(event) {
            rect=controlText.getBoundingClientRect();
            margin=7
            w=rect.right-rect.left+1-margin*2;
            offx=Math.round(event.clientX-rect.left-margin);
            opacity=offx/w; // global
            if (opacity < 0) opacity=0;
            if (opacity > 1) opacity=1;
            controlText.innerHTML = opacity_str(opacity);
            map.overlayMapTypes.removeAt(overlay_index);
            map.overlayMapTypes.insertAt(overlay_index,map_overlay());
            });
        return controlDiv;
        }

    function initialize() {
        var map = new G.Map(document.getElementById("map_canvas"));
        map.fitBounds(mapBounds);
        map.setMapTypeId(map_type);

        var overlay_index=map.overlayMapTypes.push(map_overlay())-1;

        map.controls[G.ControlPosition.TOP_RIGHT].push(opacity_control(map,overlay_index));
        }
</script>
</head>
<body onload="initialize()">
<!--    <div id="header"><h1>%(header)s</h1></div> -->
    <div id="map_canvas" style="width:100%%; height:100%%"></div>
</body>
</html>
'''

format_map=(
    Gmaps,
    GMercatorTMS,
    PlateCarree,
    PlateCarreeTMS,
    )
def formats_lst(tty=False):
    if not tty:
        return [c.format for c in format_map]    
    print('\nOutput formats and compatibility:\n')
    [print('%10s - %s' % (c.format,c.__doc__)) for c in format_map]
    print()

def proc_src(src):
    for cls in format_map:
        if cls.format == options.format:
            break
    else:
        raise Exception("Invalid format: %s" % format)

    dest=dest_path(src,options.dest_dir,cls.defaul_ext)
    cls(src,dest,options).walk_pyramid()
    
profile_map={
    'none':     None,
    'draft':    ('nearest','nearest'),
    'release':  ('antialias','bilinear'),
    }
def profile_lst(): return profile_map.keys()

def main(argv):
    
    usage = "usage: %prog <options>... input_file..."
    parser = OptionParser(usage=usage,
        description='Tile cutter for GDAL-compatible raster maps')
    parser.add_option('-f',"--to",dest="format",metavar='FORMAT',
        default='gmaps',choices=formats_lst(),
        help='output tiles format (default: gmaps)')
    parser.add_option("-l", "--list-formats", action="store_true",
        help='list output formats')
    parser.add_option("-z", "--zoom", default=None,metavar="ZOOM_LIST",
        help='list of zoom ranges to generate')
    parser.add_option('-r','--overview-resampling', default='nearest',metavar="METHOD1",
        choices=resampling_lst(),
        help='overview tiles resampling method (default: nearest)')
    parser.add_option('--base-resampling', default='nearest',metavar="METHOD2",
        choices=base_resampling_lst(),
        help='base image resampling method (default: nearest)')
    parser.add_option('-p','--profile', default='none',
        choices=profile_lst(),
        help='resampling profile (default: draft)')
    parser.add_option("-c", "--cut", action="store_true", 
        help='cut the raster as per cutline provided')
    parser.add_option("--cutline", default=None, metavar="DATASOURCE",
        help='cutline data: OGR datasource')
    parser.add_option("--cutline-blend", dest="blend_dist",default=None,metavar="N",
        help='CUTLINE_BLEND_DIST in pixels')
    parser.add_option("--no-data", dest="no_data", metavar='N[,N]...',
        help='Nodata masking values for input bands')
    parser.add_option("--tiles-prefix", default='',metavar="URL",
        help='prefix for tile URLs at googlemaps.hml')
    parser.add_option("--tile-format", default='png',metavar="FMT",
        help='tile image format (default: PNG)')
    parser.add_option("--to-palette", action="store_true", 
        help='convert tiles to paletted format (8 bit/pixel)')
    parser.add_option("--tile-size", default='256,256',metavar="SIZE_X,SIZE_Y",
        help='tile size (default: 256,256)')
    parser.add_option("-t", "--dest-dir", dest="dest_dir", default=None,
        help='destination directory (default: source)')
    parser.add_option("-q", "--quiet", action="store_const", 
        const=0, default=1, dest="verbose")
    parser.add_option("-d", "--debug", action="store_const", 
        const=2, dest="verbose")

    global options
    (options, args) = parser.parse_args(argv[1:])
    
    logging.basicConfig(level=logging.DEBUG if options.verbose==2 else 
        (logging.ERROR if options.verbose==0 else logging.INFO))

    ld(os.name)
    ld(options)
    
    if options.list_formats:
        formats_lst(tty=True)
        sys.exit(0)

    if options.profile != 'none':
        options.overview_resampling,options.base_resampling=profile_map[options.profile]

    if not args:
        parser.error('No input file(s) specified')
    try:
        sources=args
    except:
        raise Exception("No sources specified")

    parallel_map(proc_src,sources)

# main()

if __name__=='__main__':

    main(sys.argv)


#!/usr/bin/env python
# -*- coding: utf-8 -*-

# 2011-06-30 13:42:13 

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

def sasplanet_hlg2ogr(fname):
    with open(fname) as f:
        lines=f.readlines(4096)
        if not lines[0].startswith('[HIGHLIGHTING]'):
            return None
        coords=[[],[]]
        for l in lines[2:]:
            val=float(l.split('=')[1].replace(',','.'))
            coords[1 if 'Lat' in l else 0].append(val)
        points=zip(*coords)
        ld('points',points)

    ring = ogr.Geometry(ogr.wkbLinearRing)
    for p in points:
        ring.AddPoint(*p)
    polygon = ogr.Geometry(ogr.wkbPolygon)
    polygon.AddGeometry(ring)

    ds = ogr.GetDriverByName('Memory').CreateDataSource( 'wrk' )
    assert ds is not None, 'Unable to create datasource'

    src_srs = osr.SpatialReference()
    src_srs.ImportFromProj4('+proj=latlong +a=6378137 +b=6378137 +nadgrids=@null +wktext')

    layer = ds.CreateLayer('sasplanet_hlg',srs=src_srs)

    feature = ogr.Feature(layer.GetLayerDefn())
    feature.SetGeometry(polygon)
    layer.CreateFeature(feature)
    
    del feature
    del polygon
    del ring
    
    return ds

#############################

class TiledTiff(object):
    '''Tile feeder for a base zoom level'''
#############################

    def __init__(self,img_fname,tile_first,tile_last,transparency=None):
        self.fname=img_fname
        self.tile_first=tile_first
        self.tile_last=tile_last
        self.transparency=transparency
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
        self.samples_pp=self.tag_data('SamplesPerPixel')[0]

    def __del__(self):
        self.mmap.close()
        self.file.close()
                
    def tile(self,tile_x,tile_y):
        tile_sz=self.tile_sz
        ntiles=self.ntiles
        src=self.mmap
        ofs=self.tile_ofs
        lns=self.tile_lengths
        d_x,d_y=self.tile_first
        idx=tile_x-d_x+(tile_y-d_y)*self.size[0]//tile_sz[0]
        assert idx >= 0 and idx < ntiles, 'tile: %s range: %s' % ((tile_x,tile_y),self.tile_range)
        
        bands=[ src[ofs[idx+i*ntiles] : ofs[idx+i*ntiles] + lns[idx+i*ntiles]] 
                for i in range(self.samples_pp)]
        if self.samples_pp == 1:
            opacity=1
            mode='L'
            if self.transparency is not None:
                if chr(self.transparency) in bands[0]:
                    colorset=set(bands[0])
                    if len(colorset) == 1:  # fully transparent
                        return None,0
                    else:                   # semi-transparent
                        opacity=-1
            img=Image.frombuffer('L',tile_sz,bands[0],'raw','L',0,1)
        else:
            aplpha=bands[-1]
            if min(aplpha) == '\xFF':       # fully opaque
                opacity=1
                bands=bands[:-1]
                mode='RGB' if self.samples_pp > 2 else 'L'
            elif max(aplpha) == '\x00':     # fully transparent
                return None,0
            else:                           # semi-transparent
                opacity=-1
                mode='RGBA' if self.samples_pp > 2 else 'LA'
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

#############################

def BaseImg(img_fname,tile_ul,tile_lr,transparency_color=None):

#############################

    with open(img_fname,'r+b') as f:
        magic=f.read(4)
    ld('tiff magic',magic)
    for img_cls in img_type_map:
        if magic == img_cls.signature:
            break
    else:
        raise Exception('Invalid base image')

    return img_cls(img_fname,tile_ul,tile_lr,transparency_color)
    
#############################

# templates for VRT XML

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
    <!-- <WarpMemoryLimit>6.71089e+07</WarpMemoryLimit> -->
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
warp_dst_alpha_band='    <DstAlphaBand>%d</DstAlphaBand>\n'
warp_cutline='    <Cutline>%s</Cutline>\n'
warp_dst_geotr= '            <DstGeoTransform> %r, %r, %r, %r, %r, %r</DstGeoTransform>'
warp_dst_igeotr='            <DstInvGeoTransform> %r, %r, %r, %r, %r, %r</DstInvGeoTransform>'
warp_src_geotr= '            <SrcGeoTransform> %r, %r, %r, %r, %r, %r</SrcGeoTransform>'
warp_src_igeotr='            <SrcInvGeoTransform> %r, %r, %r, %r, %r, %r</SrcInvGeoTransform>'
warp_band_mapping='      <BandMapping src="%d" dst="%d"%s>'
warp_band_src_nodata='''
        <SrcNoDataReal>%d</SrcNoDataReal>
        <SrcNoDataImag>%d</SrcNoDataImag>'''
warp_band_dst_nodata='''
        <DstNoDataReal>%d</DstNoDataReal>
        <DstNoDataImag>%d</DstNoDataImag>'''
warp_band_mapping_nodata='''>%s%s
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
      <SourceBand>%(srcband)d</SourceBand>
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

class Pyramid(object):
    '''Tile pyramid generator and utilities'''
#############################

    tile_sz=(256,256)

    #############################

    def __init__(self,src=None,dest=None,options=None):

    #############################
        gdal.UseExceptions()

        self.temp_files=[]
        self.palette=None
        self.transparency=None
        self.zoom_range=None
        self.shift_x=0

        self.longlat=proj_cs2geog_cs(self.proj)
        ld('proj,longlat',self.proj,self.longlat)
        self.proj2geog=MyTransformer(SRC_SRS=self.proj,DST_SRS=self.longlat)

        self.init_grid()
        
        # init to maximum extent
        self.origin=self.pix2coord(0,(0,0))
        self.extent=self.pix2coord(0,(self.zoom0_tiles[0]*self.tile_sz[0],
                                      self.zoom0_tiles[1]*self.tile_sz[1]))
        self.src=src
        self.dest=dest
        self.options=options

    #############################

    def __del__(self):

    #############################
        if self.options and self.options.verbose < 2:
            try:
                for f in self.temp_files:
                    os.remove(f)
            except: pass

    #############################

    def init_grid(self):
        # init tile grid parameters
    #############################
        semi_circ,foo=self.proj2geog.transform_point((180,0),inv=True) # Half Equator length
        ld('semi_circ',semi_circ)
        # dimentions of a tile at zoom 0
        res0=semi_circ*2/(self.zoom0_tiles[0]*self.tile_sz[0])
        self.zoom0_res=[res0,res0]
        # offset from the upper left corner to the origin of the coord system
        self.coord_offset=[semi_circ,
                           self.zoom0_res[1]*self.tile_sz[1]*self.zoom0_tiles[1]/2]
        zoom0_tile_dim=[self.zoom0_res[0]*self.tile_sz[0],
                        self.zoom0_res[1]*self.tile_sz[1]]
        ld('zoom0_tiles',self.zoom0_tiles,'zoom0_tile_dim',zoom0_tile_dim,'coord_offset',self.coord_offset)
    
    #############################

    def walk_pyramid(self):
        'generate pyramid'
    #############################

        if not self.init_map(options.zoom):
            return

        # reproject to base zoom
        self.make_base_raster()

        tiles=[]
        for zoom in self.zoom_range:
            tile_ul,tile_lr=self.corner_tiles(zoom)
            zoom_tiles=flatten([[(zoom,x,y) for x in range(tile_ul[1],tile_lr[1]+1)] 
                                           for y in range(tile_ul[2],tile_lr[2]+1)])
            tiles.extend(zoom_tiles)
        ld('min_zoom',zoom,'tile_ul',tile_ul,'tile_lr',tile_lr,'zoom_tiles',zoom_tiles)
        self.all_tiles=frozenset(tiles)
        top_tiles=filter(None,map(self.proc_tile,zoom_tiles))

        # write top-level metadata (html/kml)
        self.write_metadata(None,[ch for img,ch,opacities in top_tiles])
        
        # cache back tiles opacity
        file_opacities=[(self.tile_path(tile),opc)
            for tile,opc in flatten([opacities for img,ch,opacities in top_tiles])]
        try:
            pickle.dump(dict(file_opacities),open(os.path.join(self.dest, 'merge-cache'),'w'))
        except:
            logging.warning("opacity cache save failed")

    #############################

    def init_map(self,zoom_parm):
        'initialize geo-parameters and generate base zoom level'
    #############################

        self.src_path=self.src
        self.tiles_prefix=options.tiles_prefix
        self.tile_ext='.'+options.tile_format.lower()
        self.src_dir,src_f=os.path.split(self.src)
        self.base=os.path.splitext(src_f)[0]
        self.base_resampling=base_resampling_map[options.base_resampling]
        self.resampling=resampling_map[options.overview_resampling]

        pf('\n%s -> %s '%(self.src,self.dest),end='')

        if os.path.isdir(self.dest):
            if options.noclobber and os.path.exists(os.path.join(self.dest,'merge-cache')):
                pf('*** Pyramid already exists: skipping',end='')
                return False
            else:
                shutil.rmtree(self.dest,ignore_errors=True)

        self.get_src_ds()
        # calculate zoom range
        self.calc_zoom(zoom_parm)
        self.base_zoom=self.zoom_range[0]        
             
        # shift target SRS to avoid crossing 180 meridian
        shifted_srs=self.shift_srs(self.base_zoom)
        shift_x,y=MyTransformer(SRC_SRS=shifted_srs,DST_SRS=self.proj).transform_point((0,0))
        if shift_x != 0:
            ld('new_srs',shifted_srs,'shift_x',shift_x,'coord_offset',self.coord_offset)
            self.coord_offset[0]+=shift_x
            self.shift_x=shift_x
            self.proj=shifted_srs
            self.proj2geog=MyTransformer(SRC_SRS=self.proj,DST_SRS=self.longlat)

        # get corners at the target SRS
        target_ds=gdal.AutoCreateWarpedVRT(self.src_ds,None,proj4wkt(shifted_srs))
        target_origin,target_extent=MyTransformer(target_ds).transform([(0,0),(target_ds.RasterXSize,target_ds.RasterYSize)])

        # clip to the max tileset area (set at the __init__)
        ld('target raster')
        ld('Upper Left ',self.origin,target_origin,self.proj2geog.transform([self.origin,target_origin]))
        ld('Lower Right',self.extent,target_extent,self.proj2geog.transform([self.extent,target_extent]))

        self.origin[0]=max(self.origin[0],target_origin[0])
        self.origin[1]=min(self.origin[1],target_origin[1])
        self.extent[0]=min(self.extent[0],target_extent[0])
        self.extent[1]=max(self.extent[1],target_extent[1])
        
        return True
        
    #############################

    def get_src_ds(self):
        'get src dataset, convert to RGB(A) if required'
    #############################
        override_srs=self.options.srs
        
        if os.path.exists(self.src):
            self.src_path=os.path.abspath(self.src)

        # check for source raster type
        src_ds=gdal.Open(self.src_path,GA_ReadOnly)
        self.src_ds=src_ds

        # source is successfully opened, then create destination dir
        os.makedirs(self.dest)

        src_geotr=src_ds.GetGeoTransform()
        src_proj=wkt2proj4(src_ds.GetProjection())
        gcps=src_ds.GetGCPs()
        if gcps:
            ld('src GCPsToGeoTransform',gdal.GCPsToGeoTransform(gcps))

        if not src_proj and gcps :
            src_proj=wkt2proj4(src_ds.GetGCPProjection())

        if self.options.srs is not None:
            src_proj=self.options.srs

        ld('src_proj',src_proj,'src geotr',src_geotr)
        assert src_proj, 'The source does not have a spatial reference system assigned'

        src_bands=src_ds.RasterCount
        band1=src_ds.GetRasterBand(1)
        if src_bands == 1 and band1.GetColorInterpretation() == GCI_PaletteIndex : # source is a paletted raster
            transparency=None
            if self.base_resampling == 'NearestNeighbour' and self.resampling == Image.NEAREST :
                # check if src can be rendered in paletted mode
                color_table=band1.GetColorTable()
                ncolors=color_table.GetCount()
                palette=[color_table.GetColorEntry(i) for i in range(ncolors)]
                r,g,b,a=zip(*palette)
                pil_palette=flatten(zip(r,g,b))             # PIL doesn't support RGBA palettes
                if self.options.dst_nodata is not None:
                    transparency=int(self.options.dst_nodata.split(',')[0])
                elif min(a) == 0:
                    transparency=a.index(0)
                elif ncolors < 256:
                    pil_palette+=[0,0,0]                   # the last color added is for transparency
                    transparency=len(pil_palette)/3-1

            ld('transparency',transparency)
            if transparency is not None: # render in paletted mode
                self.transparency=transparency
                self.palette=pil_palette
                ld('self.palette',self.palette)

            else: # convert src to rgb VRT
                if not src_geotr or src_geotr == (0.0, 1.0, 0.0, 0.0, 0.0, 1.0):
                    geotr_txt=''
                else:
                    geotr_txt=geotr_templ % src_geotr

                gcplst_txt=''
                if gcps:
                    gcp_lst='\n'.join((gcp_templ % (g.Id,g.GCPPixel,g.GCPLine,g.GCPX,g.GCPY,g.GCPZ) 
                                        for g in gcps))
                    if self.options.srs is None:
                        gcp_proj=wkt2proj4(src_ds.GetGCPProjection())
                    else:
                        gcp_proj=src_proj
                    gcplst_txt=gcplst_templ % (gcp_proj,gcp_lst)

                metadata=src_ds.GetMetadata()
                if metadata:
                    mtd_lst=[xml_txt('MDI',metadata[mdkey],4,key=mdkey) for mdkey in metadata]
                    meta_txt=meta_templ % '\n'.join(mtd_lst)
                else:
                    meta_txt=''

                xsize,ysize=(src_ds.RasterXSize,src_ds.RasterYSize)
                blxsize,blysize=band1.GetBlockSize()

                band_lst=''.join((band_templ % {
                    'band':     band,
                    'color':    color,
                    'src':      self.src_path,
                    'srcband':  1,
                    'xsize':    xsize,
                    'ysize':    ysize,
                    'blxsize':    blxsize,
                    'blysize':    blysize,
                    } for band,color in ((1,'Red'),(2,'Green'),(3,'Blue'))))
                vrt_txt=vrt_templ % {
                    'xsize':    xsize,
                    'ysize':    ysize,
                    'metadata': meta_txt,
                    'srs':      (srs_templ % src_proj) if src_proj else '',
                    'geotr':    geotr_txt,
                    'gcp_list': gcplst_txt,
                    'band_list':band_lst,
                    }

                src_vrt=os.path.join(self.dest,self.base+'.src.vrt') # auxilary VRT file
                self.temp_files.append(src_vrt)
                self.src_path=src_vrt
                with open(src_vrt,'w') as f:
                    f.write(vrt_txt)

                self.src_ds=gdal.Open(src_vrt,GA_ReadOnly)
                return # rgb VRT created
            # finished with a paletted raster

        if override_srs is not None: # src SRS needs to be relpaced
            src_vrt=os.path.join(self.dest,self.base+'.src.vrt') # auxilary VRT file
            self.temp_files.append(src_vrt)
            self.src_path=src_vrt

            vrt_drv = gdal.GetDriverByName('VRT')
            self.src_ds = vrt_drv.CreateCopy(src_vrt,src_ds) # replace src dataset

            self.src_ds.SetProjection(proj4wkt(override_srs)) # replace source SRS
            gcps=self.src_ds.GetGCPs()
            if gcps :
                self.src_ds.SetGCPs(gcps,proj4wkt(override_srs))

        # debug print
        ld('source_raster')
        src_origin,src_extent=MyTransformer(src_ds).transform([(0,0),(src_ds.RasterXSize,src_ds.RasterYSize)])
        src_proj=wkt2proj4(src_ds.GetProjection())
        src_proj2geog=MyTransformer(SRC_SRS=src_proj,DST_SRS=proj_cs2geog_cs(src_proj))
        ld('Upper Left ',src_origin,src_proj2geog.transform([src_origin]))
        ld('Lower Right',src_extent,src_proj2geog.transform([src_extent]))        

    #############################

    def shift_srs(self,zoom=None):
        'change prime meridian to allow charts crossing 180 meridian'
    #############################
        ul,lr=MyTransformer(self.src_ds,DST_SRS=self.longlat).transform([(0,0),(self.src_ds.RasterXSize,self.src_ds.RasterYSize)])
        ld('shift_srs ul',ul,'lr',lr)
        if lr[0] <= 180 and ul[0] >=-180 and ul[0] < lr[0]:
            return self.proj

        left_lon=int(math.floor(ul[0]))
        left_xy=self.proj2geog.transform_point((left_lon,0),inv=True)
        if zoom is not None: # adjust to a tile boundary
            left_xy=self.tile2coord_box(self.coord2tile(zoom,left_xy))[0]
            left_lon=int(math.floor(self.proj2geog.transform_point(left_xy)[0]))
        lon_0=left_lon+180
        ld('left_lon',left_lon,'left_xy',left_xy,'lon_0',lon_0)
        return '%s +lon_0=%d' % (self.proj,lon_0)

    #############################

    def calc_zoom(self,zoom_parm):
        'determine and set a list of zoom levels to generate'
    #############################
        res=None
        if not zoom_parm: # calculate "automatic" zoom levels
            # check raster parameters to find default zoom range
            # modify target srs to allow charts crossing meridian 180
            ld('automatic zoom levels')
            shifted_srs=self.shift_srs()

            t_ds=gdal.AutoCreateWarpedVRT(self.src_ds,None,proj4wkt(shifted_srs))
            geotr=t_ds.GetGeoTransform()
            res=(geotr[1], -geotr[5])
            max_zoom=max(self.res2zoom_xy(res))

            # calculate min_zoom
            ul_c=(geotr[0], geotr[3])
            lr_c=gdal.ApplyGeoTransform(geotr,t_ds.RasterXSize,t_ds.RasterYSize)
            wh=(lr_c[0]-ul_c[0],ul_c[1]-lr_c[1])
            ld('ul_c,lr_c,wh',ul_c,lr_c,wh)
            min_zoom=min(self.res2zoom_xy([wh[i]/self.tile_sz[i]for i in (0,1)]))
            zoom_parm='%d-%d'%(min_zoom,max_zoom)

        self.set_zoom_range(zoom_parm)
        ld(('res',res,'zoom_range',self.zoom_range,'z0 (0,0)',self.coord2pix(0,(0,0))))

    #############################

    def make_base_raster(self):

    #############################

        # adjust raster extents to tile boundaries
        tile_ul,tile_lr=self.corner_tiles(self.base_zoom)
        ld('base_raster')
        ld('tile_ul',tile_ul,'tile_lr',tile_lr)
        ul_c=self.tile2coord_box(tile_ul)[0]
        lr_c=self.tile2coord_box(tile_lr)[1]
        ul_pix=self.tile_corners(tile_ul)[0]
        lr_pix=self.tile_corners(tile_lr)[1]

        # base zoom level raster size
        dst_xsize=lr_pix[0]-ul_pix[0] 
        dst_ysize=lr_pix[1]-ul_pix[1]        

        ld('Upper Left ',self.origin,ul_c,self.proj2geog.transform([self.origin,ul_c]))
        ld('Lower Right',self.extent,lr_c,self.proj2geog.transform([self.extent,lr_c]))
        ld('coord_offset',self.coord_offset,'ul_c+c_off',map(operator.add,ul_c,self.coord_offset))

        # create VRT for base image warp

        # generate warp transform
        src_geotr=self.src_ds.GetGeoTransform()
        src_proj=wkt2proj4(self.src_ds.GetProjection())
        gcp_proj=None

        if src_geotr and src_geotr != (0.0, 1.0, 0.0, 0.0, 0.0, 1.0):
            ok,src_igeotr=gdal.InvGeoTransform(src_geotr)
            assert ok
            src_transform='%s\n%s' % (warp_src_geotr % src_geotr,warp_src_igeotr % src_igeotr)
        else:
            gcps=self.src_ds.GetGCPs()
            assert gcps, 'Neither geotransform, nor gpcs are in the source file %s' % self.src

            gcp_lst=[(g.Id,g.GCPPixel,g.GCPLine,g.GCPX,g.GCPY,g.GCPZ) for g in gcps]
            ld('src_proj',self.src_ds.GetProjection())
            ld('gcp_proj',self.src_ds.GetGCPProjection())
            gcp_proj=wkt2proj4(self.src_ds.GetGCPProjection())
            if src_proj and gcp_proj != src_proj:
                coords=MyTransformer(SRC_SRS=gcp_proj,DST_SRS=src_proj).transform([g[3:6] for g in gcp_lst])
                gcp_lst=[tuple(p[:3]+c) for p,c in zip(gcp_lst,coords)]

            gcp_txt='\n'.join((gcp_templ % g for g in gcp_lst))
            #src_transform=warp_src_gcp_transformer % (0,gcp_txt)
            src_transform=warp_src_tps_transformer % gcp_txt

        res=self.zoom2res(self.base_zoom)
        ul_ll,lr_ll=self.coords2longlat([ul_c,lr_c])
        ld('base_zoom',self.base_zoom,'size',dst_xsize,dst_ysize,'-tr',res[0],res[1],'-te',ul_c[0],lr_c[1],lr_c[0],ul_c[1])
        dst_geotr=( ul_c[0], res[0],     0.0,
                    ul_c[1],    0.0, -res[1] )
        ok,dst_igeotr=gdal.InvGeoTransform(dst_geotr)
        assert ok
        dst_transform='%s\n%s' % (warp_dst_geotr % dst_geotr,warp_dst_igeotr % dst_igeotr)

        # generate warp options
        warp_options=[]
        def w_option(name,value): # warp options template
            return '    <Option name="%s">%s</Option>' % (name,value)

        warp_options.append(w_option('INIT_DEST','NO_DATA'))

        # generate cut line
        if self.options.cut:
            cut_wkt=self.get_cutline()
        else:
            cut_wkt=None
        if cut_wkt:
            warp_options.append(w_option('CUTLINE',cut_wkt))
            if self.options.blend_dist:
                warp_options.append(w_option('CUTLINE_BLEND_DIST',self.options.blend_dist))

        src_bands=self.src_ds.RasterCount
        ld('src_bands',src_bands)

        # process nodata info
        src_nodata=None
        if self.options.src_nodata:
            src_nodata=map(int,options.src_nodata.split(','))
            assert len(src_nodata) == src_bands, 'Nodata must match the number of bands'
            if src_bands > 1:
                warp_options.append(w_option('UNIFIED_SRC_NODATA','YES'))
        dst_nodata=None
        if self.palette is not None:
            dst_nodata=[self.transparency]
        ld('nodata',src_nodata,dst_nodata)
        
        # src raster bands mapping
        vrt_bands=[]
        wo_BandList=[]
        for i in range(src_bands):
            vrt_bands.append(warp_band % (i+1,'/'))
            if src_nodata or dst_nodata:
                band_mapping_info=warp_band_mapping_nodata % (
                        warp_band_src_nodata % (src_nodata[i],0) if src_nodata else '',
                        warp_band_dst_nodata % (dst_nodata[i],0) if dst_nodata else '')
            else:
                band_mapping_info='/'
            wo_BandList.append(warp_band_mapping % (i+1,i+1,band_mapping_info))

        if src_bands < 4 and self.palette is None:
            vrt_bands.append(warp_band % (src_bands+1,warp_band_color % 'Alpha'))

        block_sz=self.tile_sz

        vrt_text=warp_vrt % {
            'xsize':            dst_xsize,
            'ysize':            dst_ysize,
            'srs':              self.proj,
            'geotr':            geotr_templ % dst_geotr,
            'band_list':        '\n'.join(vrt_bands),
            'blxsize':          block_sz[0],
            'blysize':          block_sz[1],
            'wo_ResampleAlg':   self.base_resampling,
            'wo_src_path':      self.src_path,
            'warp_options':     '\n'.join(warp_options),
            'wo_src_srs':       gcp_proj if gcp_proj else src_proj,
            'wo_dst_srs':       self.proj,
            'wo_src_transform': src_transform,
            'wo_dst_transform': dst_transform,
            'wo_BandList':      '\n'.join(wo_BandList),
            'wo_DstAlphaBand':  warp_dst_alpha_band % (src_bands+1) if src_bands < 4  and self.palette is None else '',
            'wo_Cutline':       (warp_cutline % cut_wkt) if cut_wkt else '',
            }

        temp_vrt=os.path.join(self.dest,self.base+'.tmp.vrt') # auxilary VRT file
        self.temp_files.append(temp_vrt)
        with open(temp_vrt,'w') as f:
            f.write(vrt_text)

        # warp base raster
        tmp_ds = gdal.Open(vrt_text,GA_ReadOnly)
        dst_drv = gdal.GetDriverByName('Gtiff')
        base_tiff=os.path.join(self.dest,self.base+'.tmp_%i.tiff' % self.base_zoom) # img for the base zoom
        self.temp_files.append(base_tiff)
            
        pf('...',end='')
        dst_ds = dst_drv.CreateCopy(base_tiff,tmp_ds,0,[
			'TILED=YES',
            'INTERLEAVE=BAND',
            'BLOCKXSIZE=%i' % self.tile_sz[0],
            'BLOCKYSIZE=%i' % self.tile_sz[1],
            ])#, gdal.TermProgress)
        pf('.',end='')
        
        # close datasets in a proper order
        del dst_ds
        del tmp_ds
        del self.src_ds

        # create base_image raster
        self.base_img=BaseImg(base_tiff,tile_ul[1:],tile_lr[1:],self.transparency)

    #############################

    def get_cutline(self):

    #############################
        src_ds=self.src_ds
        src_proj=wkt2proj4(src_ds.GetProjection())
        cutline=src_ds.GetMetadataItem('CUTLINE')
        ld('cutline',cutline)
        if self.options.cutline:
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
            multipoint_lst=self.shape2mpointlst(cut_file,src_proj)

            p_pix=MyTransformer(src_ds).transform_point(point_lst,inv=True)

            mpoly=[]
            for points in multipoint_lst:
                p_pix=pix_tr.transform_point(points,inv=True)
                mpoly.append(','.join(['%r %r' % (p[0],p[1]) for p in p_pix]))
            cutline='MULTIPOLYGON(%s)' % ','.join(['((%s))' % poly for poly in mpoly])
        return cutline

    #############################

    def proc_tile(self,tile):

    #############################

        ch_opacities=[]
        ch_tiles=[]
        zoom,x,y=tile
        if zoom==self.zoom_range[0]: # get from the base image
            tile_img,opacity=self.base_img.tile(*tile[1:])
            if tile_img and self.palette:
                tile_img.putpalette(self.palette)
        else: # merge children
            opacity=0
            cz=self.zoom_range[self.zoom_range.index(zoom)-1] # child's zoom
            dz=int(2**(cz-zoom))

            children_ofs=dict(flatten(
                [[((cz,x*dz+dx,y*dz+dy),(dx*self.tile_sz[0]//dz,dy*self.tile_sz[1]//dz))
                               for dx in range(dz)]
                                   for dy in range(dz)]))

            ch_tiles=filter(None,map(self.proc_tile,self.all_tiles & set(children_ofs)))
            if len(ch_tiles) == 4 and all([opacities[0][1]==1 for img,ch,opacities in ch_tiles]):
                opacity=1
                mode_opacity=''
            else:
                opacity=-1
                mode_opacity='A'

            tile_img=None
            for img,ch,opacity_lst in ch_tiles:
                ch_img=img.resize([i//dz for i in img.size],self.resampling)
                ch_mask=ch_img.split()[-1] if 'A' in ch_img.mode else None
                
                if tile_img is None:
                    if 'P' in ch_img.mode:
                        tile_mode='P'
                    else:
                        tile_mode=('L' if 'L' in ch_img.mode else 'RGB')+mode_opacity
                    
                    if self.transparency is not None:
                        tile_img=Image.new(tile_mode,self.tile_sz,self.transparency)
                    else:
                        tile_img=Image.new(tile_mode,self.tile_sz)
                    if self.palette is not None:
                        tile_img.putpalette(self.palette)

                tile_img.paste(ch_img,children_ofs[ch],ch_mask)
                ch_opacities.extend(opacity_lst)

        if tile_img is not None and opacity != 0:
            self.write_tile(tile,tile_img)
            
            # write tile-level metadata (html/kml)            
            self.write_metadata(tile,[ch for img,ch,opacities in ch_tiles])
            return tile_img,tile,[(tile,opacity)]+ch_opacities

    #############################

    def write_tile(self,tile,tile_img):

    #############################
        rel_path=self.tile_path(tile)
        full_path=os.path.join(self.dest,rel_path)
        try:
            os.makedirs(os.path.dirname(full_path))
        except: pass

        if self.options.paletted and self.tile_ext == '.png':
            try:
                tile_img=tile_img.convert('P', palette=Image.ADAPTIVE, colors=255)
            except ValueError:
                #ld('tile_img.mode',tile_img.mode)
                pass

        if self.transparency is not None:
            tile_img.save(full_path,transparency=self.transparency)
        else:
            tile_img.save(full_path)
        
        self.counter()

    #############################

    def tile_path(self,tile):
        'relative path to a tile'
    #############################
        z,x,y=self.tile_norm(tile)
        return '%i/%i/%i%s' % (z,x,y,self.tile_ext)

    #############################

    def tms_tile_path(self,tile):
        'relative path to a tile, TMS style'
    #############################
        z,x,y=self.tile_norm(tile)
        y_tiles=self.zoom0_tiles[1]*2**z
        return '%i/%i/%i%s' % (z,x,y_tiles-y-1,self.tile_ext)

    #############################

    def write_metadata(self,tile,children=[]): 

    #############################
        pass # 'virtual'


    #############################
    #
    # utility functions
    #    
    #############################

    @staticmethod
    def profile_class(profile_name):
        for cls in profile_map:
            if cls.profile == profile_name:
                return cls
        else:
            raise Exception("Invalid profile: %s" % profile_name)

    @staticmethod
    def profile_lst(tty=False):
        if not tty:
            return [c.profile for c in profile_map]    
        print('\nOutput profiles and compatibility:\n')
        [print('%10s - %s' % (c.profile,c.__doc__)) for c in profile_map]
        print()

    def zoom2res(self,zoom):
        return map(lambda res: res/2**zoom, self.zoom0_res)

    def res2zoom_xy(self,res):
        'resolution to zoom levels (separate for x and y)'
        z=map(lambda z0r,r: int(math.floor(math.log(z0r/r,2))), self.zoom0_res,res)
        return [v if v>0 else 0 for v in z]

    def pix2tile(self,zoom,pix_coord):
        'pixel coordinates to tile (z,x,y)'
        return [zoom]+map(lambda pc,ts: pc//ts, pix_coord,self.tile_sz)

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
        ul_distance=(coord[0]+self.coord_offset[0],self.coord_offset[1]-coord[1])
        out=tuple([int(dist/res) for dist,res in zip(ul_distance,self.zoom2res(zoom))])
        return out

    def pix2coord(self,zoom,pix_coord):
        pix00_ofs=map(operator.mul,pix_coord,self.zoom2res(zoom))
        return [pix00_ofs[0]-self.coord_offset[0],-(pix00_ofs[1]-self.coord_offset[1])]

    def zoom_tiles(self,zoom):
        return map(lambda v: v*2**zoom,self.zoom0_tiles)

    def tile_norm(self,tile):
        'x,y of a tile do not exceed max values'
        z,x,y=tile
        nztiles=self.zoom_tiles(z)
        return [z,x%nztiles[0],y%nztiles[1]]

    def coords2longlat(self, coords): # redefined in PlateCarree
        longlat=[i[:2] for i in self.proj2geog.transform(coords)]
        #ld('coords',coords)
        #ld('longlat',longlat)
        return longlat

    def boxes2longlat(self,box_lst):
        deg_lst=self.coords2longlat(flatten(box_lst))
        ul_lst=deg_lst[0::2]
        lr_lst=deg_lst[1::2]
        res=[[
            (ul[0] if ul[0] <  180 else ul[0]-360,ul[1]),
            (lr[0] if lr[0] > -180 else lr[0]+360,lr[1]),
            ] for ul,lr in zip(ul_lst,lr_lst)]
        return res

    def corner_tiles(self,zoom):
        p_ul=self.coord2pix(zoom,self.origin)
        t_ul=self.pix2tile(zoom,(p_ul[0]+1,p_ul[1]+1))

        p_lr=self.coord2pix(zoom,self.extent)
        t_lr=self.pix2tile(zoom,(p_lr[0]-1,p_lr[1]-1))

        nztiles=self.zoom_tiles(zoom)
        box_ul,box_lr=[self.tile2coord_box(t) for t in (t_ul,t_lr)]
        ld('corner_tiles zoom',zoom,
            'zoom tiles',nztiles,
            'zoom pixels',map(lambda zt,ts:zt*ts,nztiles,self.tile_sz),
            'p_ul',p_ul,'p_lr',p_lr,'t_ul',t_ul,'t_lr',t_lr,
            'longlat', self.coords2longlat([box_ul[0],box_lr[1]])
            )
        return t_ul,t_lr

    def belongs_to(self,tile):
        zoom,x,y=tile
        if self.zoom_range and zoom not in self.zoom_range:
            return False
        t_ul,t_lr=self.corner_tiles(zoom)
        return x>=t_ul[1] and y>=t_ul[2] and x<=t_lr[1] and y<=t_lr[2]

    def shape2mpointlst(self,datasource,target_srs):
        ds=ogr.Open(datasource)
        if not ds:
            ds=sasplanet_hlg2ogr(datasource)
        assert ds, 'Invalid datasource %s' % datasource

        layer=ds.GetLayer()
        feature=layer.GetFeature(0)
        geom=feature.GetGeometryRef()
        geom_name=geom.GetGeometryName()
        geom_lst={
            'MULTIPOLYGON':(geom.GetGeometryRef(i) for i in range(geom.GetGeometryCount())),
            'POLYGON': (geom,),
            }[geom_name]

        l_srs=layer.GetSpatialRef()
        if l_srs:
            layer_proj=l_srs.ExportToProj4()
        else:
            layer_proj=target_srs
        srs_tr=MyTransformer(SRC_SRS=layer_proj,DST_SRS=target_srs)
        if layer_proj == target_srs:
            srs_tr.transform=lambda x:x

        mpointlst=[]
        for pl in geom_lst:
            assert pl.GetGeometryName() == 'POLYGON'
            for ln in (pl.GetGeometryRef(j) for j in range(pl.GetGeometryCount())):
                assert ln.GetGeometryName() == 'LINEARRING'
                points=[ln.GetPoint(n) for n in range(ln.GetPointCount())]
                transformed_points=srs_tr.transform(points)
                mpointlst.append(transformed_points)
                
        ld('mpointlst',mpointlst)
        return mpointlst

    def set_region(self,point_lst,source_srs=None):
        if source_srs and source_srs != self.proj:
            point_lst=MyTransformer(SRC_SRS=source_srs,DST_SRS=self.proj).transform(point_lst)

        x_coords,y_coords=zip(*point_lst)[0:2]
        upper_left=min(x_coords),max(y_coords)
        lower_right=max(x_coords),min(y_coords)
        self.origin,self.extent=upper_left,lower_right

    def load_region(self,datasource):
        if not datasource:
            return
        point_lst=flatten(self.shape2mpointlst(datasource,self.proj))
        self.set_region(point_lst)

    def set_zoom_range(self,zoom_parm):
        'set a list of zoom levels from a parameter list'

        if not zoom_parm:
            return

        zchunks=[map(int,z.split('-')) for z in zoom_parm.split(',')]
        zrange=[]
        for z in zchunks:
            if len(z) == 1:
                zrange+=z
            else:
                zrange+=range(min(z),max(z)+1)
        self.zoom_range=list(reversed(sorted(set(zrange))))
    
    tick_rate=50
    count=0

    def counter(self):
        self.count+=1
        if self.count % self.tick_rate == 0:
            pf('.',end='')
            return True
        else:
            return False

# Pyramid        

#############################

class GenericMap(Pyramid):
    'full profile options are to be specified'
#############################
    profile='generic'
    defaul_ext='.generic'
    
    def __init__(self,src=None,dest=None,options=None):
        self.proj=options.t_srs
        assert self.proj, 'Target SRS is not specified'
        self.tile_sz=tuple(map(int,options.tile_size.split(',')))
        self.zoom0_tiles=map(int,options.zoom0_tiles.split(','))
        if options.tms:
            tile_path=Pyramid.tms_tile_path # method shortcut

        super(GenericMap, self).__init__(src,dest,options)

#############################

class PlateCarree(Pyramid):
    'Google Earth (plate carrée), Google tile numbering'
#############################
    profile='earth'
    defaul_ext='.earth'
    zoom0_tiles=[2,1] # tiles at zoom 0

    # http://earth.google.com/support/bin/static.py?page=guide.cs&guide=22373&topic=23750
    # "Google Earth uses Simple Cylindrical projection for its imagery base. This is a simple map 
    # projection where the meridians and parallels are equidistant, straight lines, with the two sets 
    # crossing at right angles. This projection is also known as Lat/Lon WGS84"    

    # Equirectangular (epsg:32662 aka plate carrée, aka Simple Cylindrical)
    proj='+proj=eqc +datum=WGS84 +ellps=WGS84'

    def coords2longlat(self, coords):
        out=[map(lambda c,res0,tsz,c_off: (((c+c_off)/(res0*tsz)*180)+180)%360-180,
                coord,self.zoom0_res,self.tile_sz,(self.shift_x,0)) for coord in coords]
        return [(lon,lat) for lon,lat in out]

    def kml_child_links(self,children,parent=None,path_prefix=''):
        kml_links=[]
        # convert tiles to degree boxes
        longlat_lst=self.boxes2longlat([self.tile2coord_box(t) for t in children])
        
        for tile,longlat in zip(children,longlat_lst):
            #ld(tile,longlat)
            w,n,e,s=['%.11f'%v for v in flatten(longlat)]
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

    def write_metadata(self,tile,children=[]): #
        if not tile: # create top level kml
            self.write_kml(os.path.basename(self.base),os.path.basename(self.base),self.kml_child_links(children))
            return
        # fill in kml templates
        rel_path=self.tile_path(tile)
        name=os.path.splitext(rel_path)[0]
        kml_links=self.kml_child_links(children,tile,'../../')
        tile_box=self.boxes2longlat([self.tile2coord_box(tile)])[0]
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

#############################

class PlateCarreeTMS(PlateCarree):
    'Google Earth (plate carrée), TMS tile numbering'
#############################
    profile='earth-tms'
    defaul_ext='.earth-tms'
    tile_path=Pyramid.tms_tile_path # method shortcut

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

#############################

class Yandex(Pyramid):
    'Yandex Maps (WGS 84 / World Mercator, epsg:3395)'
#############################
    profile='yandex'
    defaul_ext='.yandex'
    proj='+proj=merc +datum=WGS84 +ellps=WGS84'
# Yandex

#############################

class Gmaps(Pyramid):
    'Google Maps (Global Mercator), native tile numbering'
#############################
    profile='gmaps'
    defaul_ext='.gmaps'
    zoom0_tiles=[1,1] # tiles at zoom 0

    # Google Maps Global Mercator (epsg:3857)
    proj='+proj=merc +a=6378137 +b=6378137 +nadgrids=@null +wktext'

    def write_metadata(self,tile,children=[]): 
        if not tile: # create top level html
            self.write_html_maps()

    def write_html_maps(self):
        ul,lr=self.boxes2longlat([(self.origin,self.extent)])[0]
        googlemaps = google_templ % dict(
            title=      os.path.basename(self.dest),
            longlat_ll= '%s, %s' % (lr[1],ul[0]),
            longlat_ur= '%s, %s' % (ul[1],lr[0]),
            minzoom=    self.zoom_range[-1],
            maxzoom=    self.zoom_range[0],
            header=     os.path.basename(self.dest), 
            tms_tiles=  'true' if self.defaul_ext == '.tms' else 'false',
            map_type=   'SATELLITE', #'ROADMAP',
            tile_ext=   self.tile_ext,
            tile_size=  '%s, %s' % self.tile_sz,
            tiles_root= self.tiles_prefix,
            )
        open(os.path.join(self.dest,'gmaps.html'),'w').write(googlemaps)
# GMaps

#############################

class GMercatorTMS(Gmaps):
    'Google Maps (Global Mercator), TMS tile numbering'
#############################
    defaul_ext='.tms'
    profile='tms'
    tile_path=Pyramid.tms_tile_path # method shortcut

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
    var G=google.maps; // use G. instead of google.maps.

    var mapBounds = new G.LatLngBounds(new G.LatLng(%(longlat_ll)s), new G.LatLng(%(longlat_ur)s));
    var mapMinZoom = %(minzoom)d;
    var mapMaxZoom = %(maxzoom)d;
    var tile_ext = "%(tile_ext)s";
    var tile_size = new G.Size(%(tile_size)s);
    var tiles_root = "%(tiles_root)s";
    var tms_tiles = %(tms_tiles)s;
    var map_type = G.MapTypeId.%(map_type)s;
    var opacity = 0.5;
    var transparent_url='http://maps.gstatic.com/mapfiles/transparent.png';

    function log(msg) {
        setTimeout(function() {
            throw new Error(msg);
        }, 0);
    }

    function map_overlay(){
        return new G.ImageMapType({
            getTileUrl: function(coord, zoom) {
                max_x=1<<zoom;
                max_y=1<<zoom;
                y=coord.y;
                if(y >= max_y || y < 0)
                    return transparent_url
                x=coord.x %% max_x;
                if (x < 0)
                    x=max_x+x;
                if (tms_tiles) y=(1<<zoom)-coord.y-1;
                url=tiles_root+zoom+"/"+x+"/"+y+tile_ext;
                //log(url);
                return url;
                },
            tileSize: tile_size,
            opacity: opacity,
            isPng: (tile_ext == ".png")
            })
        }

    function opacity_str(opacity){
        var s = String(Math.round(opacity*100));
        while (s.length < 3) s = '+' + s;
        return '<++' + s + '%%+++>';
        }

    function opacity_control(map,overlay_index) {
        var controlDiv=document.createElement('DIV');
        // Set CSS styles for the DIV containing the control
        // Setting padding will offset the control from the edge of the map
        controlDiv.style.padding = '7px';
        controlDiv.id = 'op-control-div';

        // Set CSS for the control border
        var controlUI = document.createElement('DIV');
        controlUI.style.backgroundColor = 'white';
        controlUI.style.borderStyle = 'solid';
        controlUI.style.borderWidth = '1px';
        controlUI.style.cursor = 'pointer';
        controlUI.style.textAlign = 'center';
        controlUI.title = 'Click to set opacity of the overlay';
        controlUI.id = 'op-control';
        controlDiv.appendChild(controlUI);

        // Set CSS for the control interior
        var controlText = document.createElement('DIV');
        controlText.style.fontFamily = 'monospace';
        controlText.style.fontSize = '10px';
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
# google_templ

profile_map=(
    Gmaps,
    GMercatorTMS,
    PlateCarree,
    PlateCarreeTMS,
    Yandex,
    GenericMap,
    )

def proc_src(src):
    cls=Pyramid.profile_class(options.profile)
    ext= cls.defaul_ext if options.strip_dest_ext is None else ''
    dest=dest_path(src,options.dest_dir,ext)
    #
    cls(src,dest,options).walk_pyramid()

#############################

def main(argv):

#############################
    
    parser = OptionParser(
        usage = "usage: %prog <options>... input_file...",
        version=version,
        description='Tile cutter for GDAL-compatible raster maps')
    parser.add_option('-p','--profile','--to',dest="profile",metavar='PROFILE',
        default='gmaps',choices=Pyramid.profile_lst(),
        help='output tiles profile (default: gmaps)')
    parser.add_option("-l", "--list-profiles", action="store_true",
        help='list tile profiles')
    parser.add_option("--t-srs", default=None,metavar="TARGET_SRS",
        help='generic profile: PROJ.4 definition for target srs (default: None)')
    parser.add_option("--tile-size", default='256,256',metavar="SIZE_X,SIZE_Y",
        help='generic profile: tile size (default: 256,256)')
    parser.add_option("--zoom0-tiles", default='1,1',metavar="NTILES_X,NTILES_Y",
        help='generic profile: number of tiles along the axis at the zoom 0 (default: 1,1)')
    parser.add_option("--tms", action="store_true", 
        help='generic profile: generate TMS tiles (default: google)')
    parser.add_option("--srs", default=None,metavar="PROJ4_SRS",
        help="override source's spatial reference system with PROJ.4 definition")
    parser.add_option("-z", "--zoom", default=None,metavar="ZOOM_LIST",
        help='list of zoom ranges to generate')
    parser.add_option('--overview-resampling', default='nearest',metavar="METHOD1",
        choices=resampling_lst(),
        help='overview tiles resampling method (default: nearest)')
    parser.add_option('--base-resampling', default='nearest',metavar="METHOD2",
        choices=base_resampling_lst(),
        help='base image resampling method (default: nearest)')
    parser.add_option('-r','--release', action="store_true",
        help='set resampling options to (antialias,bilinear)')
    parser.add_option("-c", "--cut", action="store_true", 
        help='cut the raster as per cutline provided')
    parser.add_option("--cutline", default=None, metavar="DATASOURCE",
        help='cutline data: OGR datasource')
    parser.add_option("--cutline-blend", dest="blend_dist",default=None,metavar="N",
        help='CUTLINE_BLEND_DIST in pixels')
    parser.add_option("--src-nodata", dest="src_nodata", metavar='N[,N]...',
        help='Nodata values for input bands')
    parser.add_option("--dst-nodata", dest="dst_nodata", metavar='N',
        help='Assign nodata value for output paletted band')
    parser.add_option("--tiles-prefix", default='',metavar="URL",
        help='prefix for tile URLs at googlemaps.hml')
    parser.add_option("--tile-format", default='png',metavar="FMT",
        help='tile image format (default: PNG)')
    parser.add_option("--paletted", action="store_true", 
        help='convert tiles to paletted format (8 bit/pixel)')
    parser.add_option("-t", "--dest-dir", dest="dest_dir", default=None,
        help='destination directory (default: source)')
    parser.add_option("--noclobber", action="store_true", 
        help='skip processing if the target pyramyd already exists')
    parser.add_option("-s", "--strip-dest-ext", action="store_true",
        help='do not add a default extension suffix from a destination directory')
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
    
    if options.list_profiles:
        Pyramid.profile_lst(tty=True)
        sys.exit(0)

    if options.release:
        options.overview_resampling,options.base_resampling=('antialias','bilinear')

    if not args:
        parser.error('No input file(s) specified')
    try:
        sources=args
    except:
        raise Exception("No sources specified")

    parallel_map(proc_src,sources)
    pf('')

# main()

if __name__=='__main__':

    main(sys.argv)


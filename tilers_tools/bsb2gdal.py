#!/usr/bin/env python
# -*- coding: utf-8 -*-

# 2011-01-27 11:26:03 

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

def hdr_parms(header, patt): 
    'filter header for params starting with "patt/", if knd is empty then return comment lines'
    if patt != '!': 
        patt += '/'
    return [i[len(patt):] for i in header if i.startswith(patt)]

def hdr_parms2list(header, knd):
    return [i.split(',') for i in hdr_parms(header, knd)]

def hdr_parm2dict(header, knd):
    out={}
    for i in hdr_parms2list(header, knd)[0]:
        if '=' in i:
            (key,val)=i.split('=')
            out[key]=val
        else:
            out[key] += ','+i
    return out
    
def hdr_read(kap_name): 
    'read KAP header'
    header=[]
    f=open(kap_name, 'rU')
    for l in f:
        if '\x1A' in l:
            break
        l=l.decode('iso-8859-1','ignore')
        if l.startswith((' ','\t')):
            header[-1] += ','+l.strip()
        else:
            header.append(l.strip())
    if not (header and header[0].startswith('!') and hdr_parms(header,'BSB') and hdr_parms(header,'KNP')): 
        raise Exception("*** invalid file: %s" % kap_name)
    ld(header)
    f.close()
    return header

def assemble_parms(parm_map,parm_info):    
    check_parm=lambda s: (s not in ['NOT_APPLICABLE','UNKNOWN']) and s.replace('0','').replace('.','')
    res=' '.join([parm_map[i]+parm_info[i] for i in parm_map
                    if  i in parm_info and check_parm(parm_info[i])])
    return ' '+res if res else ''
    
def get_dtm(header):
    'get DTM northing, easting'
    dtm_parm=options.dtm_shift
    if dtm_parm is None:
        try:
            dtm_parm=hdr_parms2list(header,'DTM')[0]
            ld('DTM',dtm_parm)
            dtm=[float(s)/3600 for s in reversed(dtm_parm)]
        except IndexError: # DTM not found
            ld('DTM not found')
            dtm=[0,0]
    return dtm

def get_refs(header,dtm):
    'get a list of geo refs in tuples'
    refs=[((int(i[1]),int(i[2])),(float(i[4]),float(i[3]))) for i in hdr_parms2list(header,'REF')]
    ld('refs',refs)
    if options.dtm or options.dtm_shift:
        # alter refs as per DTM values
        refs=[(pix,(ll[0]+dtm[0],ll[1]+dtm[1])) for pix,ll in refs]
    return refs

def get_plys(header,dtm):
    'boundary polygon'
    plys=[(float(i[2]),float(i[1])) for i in hdr_parms2list(header, 'PLY')]
    if options.dtm or options.dtm_shift:
        # alter points as per DTM values
        plys=[(i[0]+dtm[0],i[1]+dtm[1]) for i in plys]
    return plys,None

def srs_refs(header, refs, options):
    'returns srs for the BSB chart projection and srs for the REF points'
    # Get a list of geo refs in tuples
    if options.srs:
        return options.srs
    # evaluate chart's projection
    knp_info=hdr_parm2dict(header, 'KNP')
    ld(knp_info)
    proj=options.proj
    if not proj:
        ld(dir(options))
        proj_id=if_set(options.proj_id,knp_info['PR'])
        try:            
            knp_parm=proj_knp[proj_id.upper()]
        except KeyError: raise Exception('*** Unsupported projection %s' % proj_id)
        # get projection and parameters
        proj=knp_parm['PROJ']
        if '+proj=utm' in proj[0]: # UTM
            # GDAL 1.7.2 doesn't seem to make use of lon_0 with UTM, but BSB doesn't use zones
            northing='10000000' if refs[0][1][1] < 0 else '0' # Southern hemisphere?
            proj='+proj=tmerc +k=0.9996 +x_0=500000 +y_0=%s +lon_0=%s' % (northing,knp_info['PP'])
        else:
            try: # extra projection parameters for BSB 3.xx, put them before KNP parms
                knq_info=hdr_parm2dict(header, 'KNQ')
                ld(knq_info)
                knq_parm=proj_knq[proj_id.upper()]
                proj+=assemble_parms(knq_parm,knq_info)
            except IndexError:  # No KNQ
                pass
            except KeyError:    # No such proj in KNQ map
                pass
            proj+=assemble_parms(knp_parm,knp_info)
    # setup a central meridian artificialy to allow charts crossing meridian 180
    leftmost=min(refs,key=lambda r: r[0][0])
    rightmost=max(refs,key=lambda r: r[0][0])
    ld('leftmost',leftmost,'rightmost',rightmost)
    if leftmost[1][0] > rightmost[1][0] and '+lon_0=' not in proj:
        proj+=' +lon_0=%i' % int(leftmost[1][0])
    
    # evaluate chart's datum
    datum=options.datum
    if datum:
        pass
    elif options.dtm or options.dtm_shift:
        datum='+datum=WGS84'
    elif not '+proj=' in proj: 
        datum='' # assume it already has a full data defined
    else:
        datum_id=if_set(options.datum_id,knp_info['GD'])
        logging.info('\t%s, %s' % (datum_id,proj_id))
        try:
            datum=datum_map[datum_id.upper()]
        except KeyError: 
            # try to guess the datum by comment and copyright string(s)
            crr=' '.join(hdr_parms(header, '!')+hdr_parms(header, 'CRR'))
            try:
                datum=[datum_guess[crr_patt] 
                    for crr_patt in datum_guess if crr_patt in crr][0]
            except IndexError:
                if dtm == [0.0,0.0]: 
                    datum='+datum=WGS84'
                else: # datum still not found
                    raise Exception('*** Unsupported or unknown datum %s. Try with "--use-dtm"' % datum_id)
            logging.warning('*** Unknown datum "%s", guessed as "%s"' % (datum_id,datum))
    srs=proj+' '+datum+' +nodefs'
    ld(srs)
    return srs

gmt_templ='''# @VGMT1.0 @GPOLYGON
# @Jp"%s"
# FEATURE_DATA
>
# @P
%s'''

def cut_poly(hdr,kap,out_srs,dtm,raster_size):
    ply_ll,ply_pix=get_plys(hdr,dtm)
    if not ply_ll:
        return '',''
    # Create cutline
    ll=''.join(['%r %r\n' % i for i in ply_ll])
    if not ply_pix: # convert cutline geo coordinates to pixel xy using GDAL's navive srs for this KAP
        pix=command(['gdaltransform','-tps','-i','-t_srs','+proj=longlat', kap],ll).splitlines()
        ply_pix=[(int(i[0]),int(i[1])) for i in pix]
    else:
        pix=['%d %d' % i for i in ply_pix]

    width,height=raster_size
    inside=[i for i in ply_pix # check if the polygon is inside the image border
        if (i[0] > 0 or i[0] < width) or (i[1] > 0 or i[1] < height)]

    if not inside:
        return '',''
    poly='POLYGON((%s))' % ','.join(pix) # Create cutline

    # convert cutline geo coordinates to the chart's srs
    poly_xy=command(['gdaltransform','-tps','-s_srs','+proj=longlat','-t_srs',out_srs],ll)
    return poly,gmt_templ % (out_srs,poly_xy)

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
        
def kap2vrt(kap,dest=None,options=None):
    if not options:
        options=Opt()

    kap=kap.decode(locale.getpreferredencoding(),'ignore')
    hdr=hdr_read(kap)                 # Read BSB header
    
    bsb_info=hdr_parm2dict(hdr,'BSB') # general BSB parameters
    bsb_name=bsb_info['NA']
    bsb_num=bsb_info['NU']
    raster_size=eval(bsb_info['RA'])
    ld(kap,bsb_name)

    if dest:
        base=os.path.split(dest)[0]
    else:
        base=dest_path(kap,options.dest_dir,
            template=('%s - '+bsb_name if options.long_names else '%s'))
    dest_dir=os.path.split(base)[0]
    kap_path=os.path.relpath(kap,dest_dir)
    out_dataset= os.path.basename(base+'.vrt') # output VRT file
    logging.info(' %s : %s -> %s' % (kap,bsb_name,out_dataset))

    dtm=get_dtm(hdr)                    # DTM shifts
    refs=get_refs(hdr,dtm)              # reference points
    out_srs=srs_refs(hdr,refs,options)  # estimate SRS

    # get cut polygon
    poly,gmt_data=cut_poly(hdr,kap,out_srs,dtm,raster_size)
    if options.get_cutline: # print cutline and finish
        print poly
        return
    if gmt_data and not options.no_cut_file: # create shapefile with a cut polygon
        f=open(base+'.gmt','w+')
        f.write(gmt_data)
        f.close()
        
    # convert latlong gcps to projected coordinates
    latlong='\n'.join(['%f %f' % ll for pix,ll in refs])
    refs_out=command(['gdaltransform','-tps','-s_srs','+proj=longlat','-t_srs',out_srs], latlong)
    refs_proj=[ i.split() for i in refs_out.splitlines()]
    # create vrt
    gcps=flatten([['-gcp']+map(repr, rp[0])+rc for rp,rc in zip(refs, refs_proj)])
    transl_cmd=['gdal_translate','-of','VRT','-a_srs',out_srs]
    if options.last_column_bug: # see http://trac.osgeo.org/gdal/ticket/3777
        transl_cmd+=['-srcwin', '0', '0', str(raster_size[0]-1), str(raster_size[1]) ]
    if options.broken_raster: # see http://trac.osgeo.org/gdal/ticket/3790
        idx_png= base+'.idx.png'
        rgb_png= base+'.rgb.png'
        command(['gdal_translate','-of','png',kap,idx_png])
        try:
            command(['convert','-type','TrueColor',idx_png,rgb_png])
        except: pass
        os.remove(idx_png)
        transl_cmd+=[rgb_png,out_dataset]
    else: 
        transl_cmd+=[kap_path,out_dataset]
        if options.expand:
            transl_cmd=transl_cmd+['-expand',options.expand]
    
    try:
        cdir=os.getcwd()
        if dest_dir:
            os.chdir(dest_dir)
        command(transl_cmd + gcps) # gdal_translate
    finally:
        os.chdir(cdir)

def proc_src(src):
    kap2vrt(src,options=options)

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
    parser.add_option("--dtm", action="store_true", 
        help='use BSB datum shifts record to convert to WGS84')
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


#!/usr/bin/env python

# 2010-12-28 15:22:23 

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
from optparse import OptionParser
from subprocess import Popen, PIPE
import itertools
imap=itertools.imap

knp_map={ # projection parameters
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
    }
knq_map={ # extra projection parameters for BSB v. 3.xx
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
    'ROMA DATUM 1940':      '+towgs84=-104.1,-49.1,-9.9,0.971,-2.917,0.714,-11.68 +ellps=intl',
    'ROMA 1940':            '+towgs84=-104.1,-49.1,-9.9,0.971,-2.917,0.714,-11.68 +ellps=intl',
    'HERMANSKOGEL DATUM':   '+datum=hermannskogel',
    'OSGB36':               '+towgs84=446.448,-125.157,542.060,0.1502,0.2470,0.8421,-20.4894 +ellps=airy',
    # 'LOCAL DATUM'
    # 'LOCAL DATUM UNKNOWN'
    }
guessed_datum_map={ # guess the datum by a comment/copyright string pattern
    'Croatia':
        # http://spatial-analyst.net/wiki/index.php?title=MGI_/_Balkans_coordinate_systems
        # '+ellps=bessel +towgs84=550.499,164.116,475.142,5.80967,2.07902,-11.62386,0.99999445824',
        # http://earth-info.nga.mil/GandG/coordsys/onlinedatum/DatumTable.html
        '+datum=hermannskogel',
    }
    
def ld(*parms):
    logging.debug(' '.join(imap(repr,parms)))

def pf(*parms,**kparms):
    try:
        if not options.verbose:
            return
    except:
        return
    end=kparms['end'] if 'end' in kparms else '\n'
    sys.stdout.write(' '.join(imap(str,parms))+end)
    sys.stdout.flush()
    
try:
    import win32pipe
except:
    win32pipe=None
    
def command(params,child_in=None):
    cmd_str=' '.join(('"%s"' % i if ' ' in i else i for i in params))
    ld((cmd_str,child_in))
    if win32pipe:
        (stdin,stdout,stderr)=win32pipe.popen3(cmd_str,'t')
        if child_in:
            stdin.write(child_in)
        stdin.close()
        child_out=stdout.read()
        child_err=stderr.read()
        if child_err:
            logging.warning(child_err)
    else:
        process=Popen(params,stdin=PIPE, stdout=PIPE, stderr=PIPE, universal_newlines=True)
        (child_out,child_err)=process.communicate(child_in)
        if process.returncode != 0: 
            raise Exception("*** External program failed: %s\n%s" % (cmd_str,child_err))
    ld((child_out,child_err))
    return child_out

def flatten(two_level_list): 
    return list(itertools.chain(*two_level_list))

def hdr_parms(header, patt): 
    'filter header for params starting with "patt/", if knd is empty then return comment lines'
    if patt != '!': patt += '/'
    return [i[len(patt):] for i in header if i.startswith(patt)]

def hdr_parm2dict(header, knd):
    out={}
    for i in hdr_parms(header, knd)[0].split(','):
        if '=' in i:
            (key,val)=i.split('=')
            out[key]=val
        else:
            out[key] += ','+i
    return out

def hdr_read(kap_name): 
    'read KAP header'
    header=[]
    for l in open(kap_name, 'rU'):
        if '\x1A' in l:
            break
        if l.startswith((' ','\t')):
            header[-1] += ','+l.strip()
        else:
            header.append(l.strip())
    if not (header and header[0].startswith('!') and hdr_parms(header,'BSB') and hdr_parms(header,'KNP')): 
        raise Exception("*** invalid file: %s" % kap_name)
    ld(header)
    return header

def assemble_parms(parm_map,parm_info):    
    check_parm=lambda s: (s not in ['NOT_APPLICABLE','UNKNOWN']) and s.translate(None,'0.')
    res=' '.join([parm_map[i]+parm_info[i] for i in parm_map
                    if  i in parm_info and check_parm(parm_info[i])])
    return ' '+res if res else ''
    
def srs_refs(bsb_header, options):
    'returns srs for the BSB chart projection and srs for the REF points'
    # Get a list of geo refs in tuples
    refs=map(eval,hdr_parms(bsb_header, 'REF'))
    ld(refs)
    if options.srs:
        return(options.srs,refs)
    # evaluate chart's projection
    knp_info=hdr_parm2dict(bsb_header, 'KNP')
    ld(knp_info)
    bsb_proj=knp_info['PR']
    bsb_datum=knp_info['GD']
    logging.info('\t%s, %s' % (bsb_datum,bsb_proj))
    if options.bsb_proj:
        bsb_proj=options.bsb_proj
    proj=options.proj
    if not proj:
        try:            
            knp_parm=knp_map[bsb_proj.upper()]
        except KeyError: raise Exception('*** Unsupported projection %s' % bsb_proj)
        # get projection and parameters
        proj=knp_parm['PROJ']
        if '+proj=utm' in proj[0]: # UTM
            # GDAL 1.7.2 doesn't seem to make use of lon_0 with UTM, but BSB doesn't use zones
            #zone=(eval(knp_info['PP']) + 3 % 360) // 6 + 30
            #proj.append('+zone='+str(zone))
            #if refs[0][3] < 0:
            #    proj.append('+south')
            northing='10000000' if refs[0][3] < 0 else '0' # Southern hemisphere?
            proj='+proj=tmerc +k=0.9996 +x_0=500000 +y_0=%s +lon_0=%s' % (northing,knp_info['PP'])
        else:
            try: # extra projection parameters for BSB 3.xx, put them before KNP parms
                knq_info=hdr_parm2dict(bsb_header, 'KNQ')
                ld(knq_info)
                knq_parm=knq_map[bsb_proj.upper()]
                proj+=assemble_parms(knq_parm,knq_info)
            except: pass
            proj+=assemble_parms(knp_parm,knp_info)
    # setup a central meridian artificialy to allow charts crossing meridian 180
    if '+lon_0=' not in proj:
        leftmost_ref=min(refs,key=lambda r: r[1])
        ld(leftmost_ref)
        proj+=' +lon_0=%i' % int(leftmost_ref[4])
    
    # evaluate chart's datum
    if options.bsb_datum:
        bsb_datum=options.bsb_datum
    datum=options.datum
    try: # get DTM northing, easting
        dtm_parm=hdr_parms(bsb_header, 'DTM')[0] if not options.dtm else options.dtm
        dtm=eval(dtm_parm)
        ld(dtm)
    except:
        dtm=None
    if options.use_dtm or options.dtm:
        # Use DTM northing/easting correctiongs towards WGS84
        datum='+datum=WGS84'
        # alter refs as per dtm values
        refs=[ref[0:3]+(ref[3]+dtm[0]/3600,ref[4]+dtm[1]/3600) for ref in refs]
    if not datum:
        try:
            datum=datum_map[bsb_datum.upper()]
        except KeyError: 
            # try to guess the datum by comment and copyright string(s)
            crr=' '.join(hdr_parms(bsb_header, '!')+hdr_parms(bsb_header, 'CRR'))
            try:
                datum=[guessed_datum_map[crr_patt] 
                    for crr_patt in guessed_datum_map if crr_patt in crr][0]
            except:
                if dtm == (0,0): 
                    datum='+datum=WGS84'
                else: # datum still not found
                    raise Exception('*** Unsupported or unknown datum %s. Try with "--use-dtm"' % bsb_datum)
            logging.warning('*** Unknown datum "%s", guessed as "%s"' % (bsb_datum,datum))
    srs=proj+' '+datum+' +nodefs'
    ld(srs)
    return (srs,refs)

gmt_templ='''# @VGMT1.0 @GPOLYGON
# @Jp"%s"
# FEATURE_DATA
>
# @P
%s'''

def cut_poly(kap,hdr,out_srs):
    # Gather border polygon coordinates early
    plys=map(eval,hdr_parms(hdr, 'PLY'))
    if not plys:
        return '',''
    # Create cutline
    poly_latlong=''.join(['%f %f\n' % (i[2],i[1]) for i in plys])
    # convert cutline geo coordinates to pixel xy using GDAL's navive srs for this KAP
    poly_pixels=command(['gdaltransform','-tps','-i','-t_srs','+proj=longlat', kap], poly_latlong).splitlines()
    poly='POLYGON((%s))' % ','.join(poly_pixels )

    # convert cutline geo coordinates to the chart's srs
    out=command(['gdaltransform','-tps','-s_srs','+proj=longlat','-t_srs',out_srs], poly_latlong)
    return poly,gmt_templ % (out_srs,out)

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

    hdr=hdr_read(kap)           # Read BSB header
    
    bsb_info=hdr_parm2dict(hdr, 'BSB') # general BSB parameters
    bsb_name=bsb_info['NA']
    bsb_num=bsb_info['NU']
    kap_size=eval(bsb_info['RA'])
    logging.info(' %s - %s' % (kap,bsb_name))

    if dest:
        base=os.path.split(dest)[0]
    else:
        base=dest_path(kap,options.dest_dir,
            template=('%s - '+bsb_name if options.long_names else '%s'))
    dest_dir=os.path.split(base)[0]
    kap_path=os.path.relpath(kap,dest_dir)
    out_dataset= os.path.basename(base+'.vrt') # output VRT file

    # estimate SRS
    (out_srs,refs)=srs_refs(hdr, options)

    # get cut polygon
    poly,gmt_data=cut_poly(kap,hdr,out_srs)
    if options.get_cutline: # print cutline then return
        print poly
        return
    if gmt_data and not options.no_cut_file: # create shapefile with a cut polygon
        f=open(base+'.gmt','w+')
        f.write(gmt)
        f.close()
        
    # convert latlong gcps to projected coordinates
    latlong=''.join(['%f %f\n' % (ref[4],ref[3]) for ref in refs])
    refs_out=command(['gdaltransform','-tps','-s_srs','+proj=longlat','-t_srs',out_srs], latlong)
    refs_proj=[ i.split() for i in refs_out.splitlines()]
    # create vrt
    gcps=flatten([('-gcp', `i[0][1]`,`i[0][2]`,i[1][0],i[1][1]) for i in zip(refs, refs_proj)])
    transl_cmd=['gdal_translate','-of','VRT','-a_srs',out_srs]
    if options.last_column_bug: # see http://trac.osgeo.org/gdal/ticket/3777
        transl_cmd+=['-srcwin', '0', '0', str(kap_size[0]-1), str(kap_size[1]) ]
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
    
    try:
        cdir=os.getcwd()
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
    parser.add_option("--get-cutline", action="store_true", 
        help='print cutline polygon from KAP file then exit')
    parser.add_option("--no-cut-file", action="store_true", 
        help='do not create a file with a cutline polygon from KAP file')
    parser.add_option("-u", "--use-dtm", action="store_true", 
        help='use BSB datum shifts record to convert to WGS84')
    parser.add_option("--dtm", dest="dtm",default=None,
        help='specify BSB datum shifts to convert SRS to WGS84')
    parser.add_option("--srs", default=None,
        help='override full chart with PROJ.4 definition of the spatial reference system')
    parser.add_option("--datum", default=None,
        help="override chart's datum (PROJ.4 definition)")
    parser.add_option("--bsb-datum", default=None,
        help="override chart's datum (BSB definition)")
    parser.add_option("--proj", default=None,
        help="override chart's projection (BSB definition)")
    parser.add_option("--bsb-proj", default=None,
        help="override chart's projection (BSB definition)")
    parser.add_option("-l", "--long-names", action="store_true", 
        help='give an output file a long name')
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

    if options.use_dtm:
        options.datum='+datum=WGS84'

    map(proc_src,args)


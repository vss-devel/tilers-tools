#!/usr/bin/env python
# -*- coding: utf-8 -*-

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

import os
import sys
import re
import logging
import optparse
import xml.dom.minidom
import sqlite3
import base64
import htmlentitydefs
import csv
import json
from PIL import Image

from tiler_functions import *

def re_subs(sub_list,l):
    for (pattern,repl) in sub_list:
        l=re.sub(pattern,repl,l)
    return l

htmlentitydefs.name2codepoint['apos']=27

def strip_html(text):
    'Removes HTML markup from a text string. http://effbot.org/zone/re-sub.htm#strip-html'

    def replace(match): # pattern replacement function
        text = match.group(0)
        if text == '<br>':
            return '\n'
        if text[0] == '<':
            return '' # ignore tags
        if text[0] == '&':
            if text[1] == '#':
                try:
                    if text[2] == 'x':
                        return unichr(int(text[3:-1], 16))
                    else:
                        return unichr(int(text[2:-1]))
                except ValueError:
                    pass
            else:
                return unichr(htmlentitydefs.name2codepoint[text[1:-1]])
        return text # leave as is
        # fixup end

    return re.sub('(?s)<[^>]*>|&#?\w+;', replace, text)

def attr_update(self,**updates):
        self.__dict__.update(updates)

class Category (object):
    def __init__(self,label):
        attr_update(self,label=label,enabled=1,desc='',cat_id=None,icons={})

    def update(self,enabled=None,desc=None,cat_id=None,icon=None, url=None):
        if icon:
            self.icons[icon]=url
        if desc:
            self.desc=desc
        if cat_id:
            self.cat_id=cat_id
        if enabled is not None:
            self.enabled= 1 if enabled else 0;

class Poi (object):
    def __init__(self,label=None,lat=None,lon=None,desc='',categ=None):
        if desc is None:
            desc=''
        attr_update(self,label=label,desc=desc,lat=lat,lon=lon,categ=categ.lower())

class Poi2Db (object):
    ''' Generic class '''
    def __init__ (self,src,dest_db):
        attr_update(self,src=src,categories={},styles={},icons={},pois=[])

        self.toGMercator=MyTransformer(SRC_SRS='+init=epsg:4326',DST_SRS='+init=epsg:3857')

        if dest_db:
            self.base=os.path.splitext(dest_db)[0]
            self.dest_db=dest_db
        else:
            self.base=os.path.splitext(src[0])[0]
            self.dest_db=self.base+'.db'
        self.base_dir=os.path.split(self.base)[0]

        if os.path.exists(self.dest_db):
            if options.remove_dest:
                os.remove(self.dest_db)
            else:
                raise Exception('Destination already exists: %s' % self.dest_db)

    def categ_add_update(self,label=None,enabled=1,desc=None,cat_id=None,icon=None, url=None):
        if not icon:
            icon=label+'.jpg'
        if not label:
            label=re.sub(r'\.[^.]*$','',icon)
        categ=label.lower()
        ic_id=icon.lower()
        if ic_id not in self.icons:
            self.icons[ic_id]=categ
        else:
            categ=self.icons[ic_id]
        if categ not in self.categories:
            self.categories[categ]=Category(label)
        self.categories[categ].update(enabled=enabled,desc=desc,cat_id=cat_id,icon=icon, url=url)
        return categ

    def load_categ(self,src):
        path=os.path.splitext(src)[0] + '.categories'
        if os.path.exists(path):
            cats_lst=[[unicode(i.strip(),'utf-8') for i in l.split(',',4)]
                        for l in open(path)]
            for d in cats_lst:
                try:
                    (enabled,icon,categ,desc) = d + [None for i in range(len(d),4)]
                    log(enabled,icon,categ,desc)
                    if enabled not in ('0','1') or not categ:
                        continue
                    self.categ_add_update(categ,int(enabled),icon=icon,desc=desc)
                except: pass

    def read_csv(self,path):
        col_id={
            'name': None,
            'desc': None,
            'lat':  None,
            'lon':  None,
            'categ':None,
            }
        csv.register_dialect('strip', skipinitialspace=True)
        with open(path,'rb') as data_f:
            data_csv=csv.reader(data_f,'strip')
            header=[s.decode('utf-8').lower() for s in data_csv.next()]

            for col in range(len(header)): # find relevant colunm numbers
                for id in col_id:
                    if header[col].startswith(id):
                        col_id[id]=col

            cat_ids={}
            for row in data_csv:
                row=[s.decode('utf-8') for s in row]
                poi_parms={}
                for col in col_id:
                    try:
                        poi_parms[col]=row[col_id[col]]
                    except:
                        poi_parms[col]=''
                if poi_parms['categ']:
                    icon=poi_parms['categ'].lower()+'.jpg'
                else:
                    icon='__undefined__.jpg'

                categ=self.categ_add_update(icon=icon)
                self.pois.append(Poi(
                    poi_parms['name'],
                    categ=categ,
                    lat=poi_parms['lat'],
                    lon=poi_parms['lon'],
                    desc=poi_parms['desc']
                    ))

    def handleStyle(self,elm):
        url=None
        style_id=elm.getAttribute('id')
        log(style_id)
        icon=u'__%s__.jpg' % style_id
        if elm.getElementsByTagName('IconStyle') != []:
            try:
                url=elm.getElementsByTagName('href')[0].firstChild.data
                icon=re.sub('^.*/','',url)
            except: pass
        elif elm.getElementsByTagName('PolyStyle') != []:
            icon=u'__polygon__.jpg'
        elif elm.getElementsByTagName('LineStyle') != []:
            icon=u'__line__.jpg'
        return (style_id, self.categ_add_update(None,icon=icon,url=url))

    def get_coords(self,elm):
        coords=elm.getElementsByTagName('coordinates')[0].firstChild.data.split()
        return [map(float,c.split(',')) for c in coords]

    def handlePlacemark(self,pm):
        point=pm.getElementsByTagName('Point')
        if point == []:
            return None
        coords=self.get_coords(point[0])
        (lon,lat)=coords[0][0:2]

        label=pm.getElementsByTagName('name')[0].firstChild.data
        style_id=pm.getElementsByTagName('styleUrl')[0].firstChild.data[1:]
        style=self.styles[style_id]
        log((label,style_id,style))
        if style.startswith('__') and style.endswith('__'):
            logging.warning(' No icon for "%s"' % label)
        desc=None
        try:
            desc_elm=pm.getElementsByTagName('description')[0]
            if desc_elm.firstChild:
                cdata=(desc_elm.firstChild.nodeType == self.doc.CDATA_SECTION_NODE)
                desc=desc_elm.firstChild.data
                if cdata:
                    desc=strip_html(desc)
        except IndexError:
            pass
        return Poi(label,lat=lat,lon=lon,desc=desc,categ=self.styles[style_id])

    def write_aux(self):
        cat_list=['# enabled, icon, category, desc']
        icon_urls=[]
        icon_aliases=[] #"ln -s '%s.db' 'poi.db'"  % self.base]
        icon_aliase_templ="ln -s '%s' '%s'"
        for (c_key,c) in self.categories.iteritems():
            for icon_name in c.icons:
                cat_list.append('%i, %s, %s%s' % (c.enabled,icon_name,c.label,((', '+c.desc) if c.desc else '')))
                url=c.icons[icon_name]
                if url:
                    icon_urls.append("wget -nc '%s'" % url)
                if c_key+'.jpg' != icon_name:
                    icon_aliases.append(icon_aliase_templ % (icon_name, c_key+'.jpg'))
        with open(self.base+'.categories.gen','w') as f:
            for s in cat_list:
                print >>f, s
        with open(self.base+'.sh','w') as f:
            for ls in [icon_urls,icon_aliases]:
                for s in ls:
                    print >>f, s

    def proc_src(self,src):
        log(src)
        self.load_categ(src)
        try: # to open as kml file
            self.doc = xml.dom.minidom.parse(src)
            self.name=[n for n in self.doc.getElementsByTagName('Document')[0].childNodes
                        if n.nodeType == n.ELEMENT_NODE and n.tagName == 'name'][0].firstChild.data

            self.styles=dict(map(self.handleStyle,self.doc.getElementsByTagName('Style')))
            self.pois+=filter(None,map(self.handlePlacemark,self.doc.getElementsByTagName('Placemark')))
            self.doc.unlink()
        except IOError:
            logging.warning(' No input file: %s' % src)
        except xml.parsers.expat.ExpatError:
            try: # to open as db
                self.read_db(src)
            except sqlite3.DatabaseError:
                try: # to open as csv
                    self.read_csv(src)
                except csv.Error:
                    raise Exception('Invalid input file: %s' % src)

    def proc_icon(self, c):
        pass

    def proc_all(self):
        map(self.proc_src, self.src)

        self.create_db()

        map(self.proc_category,self.categories.itervalues())
        map(self.proc_poi,self.pois)
        #map(self.proc_icon,self.pois)

        self.write_aux()
        self.close_db()
#
# Poi2Db

class Poi2Mapper (Poi2Db):
    def read_db(self,path):
        cat_ids={}
        if os.path.exists(path):
            db=sqlite3.connect(path)
            dbc=db.cursor()
            dbc.execute('SELECT * FROM category')
            for (cat_id,label,desc,enabled) in dbc:
                self.categ_add_update(label,enabled,desc=desc)
                cat_ids[cat_id]=label

            dbc.execute('SELECT * FROM poi')
            for row in dbc:
                (poi_id,lat,lon,name,desc,cat_id)=row[:6]
                self.pois.append(Poi(name,lat=lat,lon=lon,desc=desc,categ=cat_ids[cat_id]))
            db.close()

    def proc_category_icon(self, c):
        #log('c.icons',c.icons)
        icon_file=os.path.join(self.base_dir,c.label+'.jpg')
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!! Make case insensitive !!!!!!!!!!!!!!!!!!!!!!!!!!
        if not os.path.exists(icon_file):
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!! Make case insensitive !!!!!!!!!!!!!!!!!!!!!!!!!!
            logging.warning('No icon image for %s' % c.label)
            return
        with open(icon_file) as f:
            icon_data = f.read()
        mime_type = mime_from_ext(ext_from_buffer(icon_data))
        img = Image.open(icon_file)
        width, height = img.size
        x_offset, y_offset = (width/2, 0)
        #log('icon',c.label,img.format,width,height)
        self.dbc.execute(
            'INSERT INTO icon '
                '(cat_id,mime_type,data,width,height,x_offset,y_offset) '
                'VALUES (?,?,?,?,?,?,?);',
            (c.cat_id,mime_type,base64.b64encode(icon_data),width,height,x_offset,y_offset))

    def proc_category(self,c):
        self.dbc.execute('INSERT INTO category (label, desc, enabled) VALUES (?,?,?);',
            (c.label,c.desc,c.enabled))
        c.update(cat_id=self.dbc.lastrowid)
        self.proc_category_icon(c)

    def proc_poi(self,p):
        log('poi',p.label,p.lon,p.lat)
        x,y=self.toGMercator.transform_point([p.lon,p.lat])
        self.dbc.execute('INSERT INTO poi (x,y,lon,lat,label,desc,cat_id) VALUES (?,?,?,?,?,?,?);',
            (x,y,p.lon,p.lat,p.label,p.desc,self.categories[p.categ].cat_id))

    def create_db(self):

        self.db=sqlite3.connect(self.dest_db)
        self.dbc = self.db.cursor()
        try:
            self.db.execute(
                'CREATE TABLE IF NOT EXISTS category ('
                    'cat_id INTEGER PRIMARY KEY,'
                    'label TEXT,'
                    'desc TEXT,'
                    'enabled INTEGER'
                    ')'
                )
            self.db.execute(
                'CREATE TABLE IF NOT EXISTS poi ('
                    'poi_id INTEGER PRIMARY KEY,'
                    'lat REAL,'
                    'lon REAL,'
                    'label TEXT,'
                    'desc TEXT,'
                    'cat_id INTEGER,'
                    'x REAL,'
                    'y REAL'
                    ')'
                )
            self.db.execute(
                'CREATE TABLE IF NOT EXISTS icon ('
                    'icon_id INTEGER PRIMARY KEY,'
                    'cat_id INTEGER,'
                    'width INTEGER,'
                    'height INTEGER,'
                    'x_offset INTEGER,'
                    'y_offset INTEGER,'
                    'mime_type TEXT,'
                    'data TEXT'
                    ')'
                )
        except:
            pass

    def close_db(self):
        self.db.commit()
        self.db.close()
#
# class Poi2Mapper

class Poi2Mmap (Poi2Mapper):

    def __init__ (self, src, dest_db):
        super(Poi2Mmap, self).__init__(src,dest_db)
        self.inserted_icons = set()

    def read_db(self,path):
        #cat_ids={}
        if os.path.exists(path):
            db=sqlite3.connect(path)
            db.row_factory = sqlite3.Row
            dbc=db.cursor()
            dbc.execute('SELECT * FROM categories')
            for row in dbc:
                log('row',row)
                props = json.loads(row['properties'])
                self.categ_add_update(props['label'], props['enabled'], desc=props['description'])
                #cat_ids[cat_id] = props['label']

            dbc.execute('SELECT * FROM pois')
            for row in dbc:
                props = json.loads(row['properties'])
                self.pois.append(Poi(props['label'], lat=props['lat'], lon=props['lon'], desc=props['description'], categ=props['category']))
            db.close()

    def proc_icon(self, i):
        if i in self.inserted_icons:
            return
        self.inserted_icons.add(i)

        icon_file = os.path.join(self.base_dir, i)
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!! Make case insensitive !!!!!!!!!!!!!!!!!!!!!!!!!!
        if not os.path.exists(icon_file):
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!! Make case insensitive !!!!!!!!!!!!!!!!!!!!!!!!!!
            logging.warning('No icon image %s' % i)
            return

        with open(icon_file) as f:
            icon_data = f.read()

        mime_type = mime_from_ext(ext_from_buffer(icon_data))
        img = Image.open(icon_file)
        width, height = img.size
        x_offset, y_offset = (width/2, 0)
        #log('icon',c.label,img.format,width,height)

        props = json.dumps({
            'width': width,
            'height': height,
            'x_offset': x_offset,
            'y_offset': y_offset,
            'label': i,
            })
        data = 'data:' + mime_type + ';base64,' + base64.b64encode(icon_data)

        self.dbc.execute('INSERT INTO icons (properties, data) VALUES (?,?);',(props,data))

    def proc_category(self, c):
        icon = c.label+'.jpg'
        props = json.dumps({
            'description': c.desc,
            'icon': icon,
            'enabled': c.enabled,
            'label': c.label,
            })

        self.dbc.execute('INSERT INTO categories (properties) VALUES (?);',(props,))
        c.update(cat_id=self.dbc.lastrowid)

        log('c.icons',c.icons)
        map(self.proc_icon, [icon] + c.icons.keys())

    def proc_poi(self, p):
        log('poi',p.label,p.lon,p.lat)
        x, y = self.toGMercator.transform_point([p.lon, p.lat])
        props = json.dumps({
            'label': p.label,
            'description': p.desc,
            'lon': p.lon,
            'lat': p.lat,
            'category': self.categories[p.categ].label
            })
        geom = '{"type":"Point","coordinates":[%f,%f]}' % (x,y)
        self.dbc.execute(
            'INSERT INTO pois (xmin,xmax,ymin,ymax,geometry,properties) VALUES (?,?,?,?,?,?);',
            (x,y,x,y,geom,props))

    def create_db(self):

        self.db=sqlite3.connect(self.dest_db)
        self.dbc = self.db.cursor()
        try:
            self.dbc.execute ('PRAGMA auto_vacuum = INCREMENTAL;')
        finally:
            pass
        try:
            self.dbc.execute (
                'CREATE TABLE IF NOT EXISTS __WebKitDatabaseInfoTable__ ('
                    'key TEXT NOT NULL ON CONFLICT FAIL UNIQUE ON CONFLICT REPLACE,'
                    'value TEXT NOT NULL ON CONFLICT FAIL'
                    ');'
                )
            self.dbc.execute (
                "INSERT OR REPLACE INTO __WebKitDatabaseInfoTable__ VALUES('WebKitDatabaseVersionKey','');"
                )
            for table in ('categories', 'icons', 'pois'):
                self.dbc.execute((
                    'CREATE TABLE IF NOT EXISTS "%s" ('+
                        'id INTEGER PRIMARY KEY,'+
                        '"group" INTEGER,'+
                        'rank INTEGER,'+
                        'xmin FLOAT,'+
                        'xmax FLOAT,'+
                        'ymin FLOAT,'+
                        'ymax FLOAT,'+
                        'geometry TEXT,'+
                        'properties TEXT,'+
                        'data TEXT'+
                    ')') % table )
        #~ except:
        finally:
            pass
#
# class Poi2Mmap

if __name__=='__main__':
    parser = optparse.OptionParser(
        usage='usage: %prog [-o <output_db>] [<kml_file>]... [<input_db>]...',
        version=version,
        description='makes maemo-mapper POI db from a kml file(s)')
    parser.add_option('-d', '--debug', action='store_true', dest='debug')
    parser.add_option('-q', '--quiet', action='store_true', dest='quiet')
    parser.add_option('-r', '--remove-dest', action='store_true',
        help='delete destination before processing')
    parser.add_option('-o', '--output',dest='dest_db',
                      type='string',help='output POIs db file')

    (options, args) = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG if options.debug else
        (logging.ERROR if options.quiet else logging.INFO))

    if args == []:
        raise Exception('No source specified')

    #Poi2Mapper(args,options.dest_db).proc_all()
    Poi2Mmap(args, options.dest_db).proc_all()

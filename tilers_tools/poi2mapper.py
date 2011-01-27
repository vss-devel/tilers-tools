#!/usr/bin/env python

# 2010-10-28 17:49:52 

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

def l_d(smth, nl=True):
    logging.debug(repr(smth))

def p_f(smth, nl=True):
#    return
    #logging.debug(str(smth))
    s=str(smth)
    if nl: s+='\n'
    sys.stdout.write(s)
    sys.stdout.flush()

def re_subs(sub_list,l):
    for (pattern,repl) in sub_list:
        l=re.sub(pattern,repl,l)
    return l

cdata_subs=[
    ('</?span[^>]*>',''),
    ('(<br>)+',      '\n'),
    ('<[^>]*>',      ' '),
    ('&amp;',        '&'),
    ('&quot;',       '"'),
    ('&apos;',       "'"),
    ('&lt;',         '<'),
    ('&gt;',         '>'),
    ('&#39;',        "'"),
    ('^ *',          ''),
    ('\n *$',        ''),
]
def cdata_cleanup(l):
    return re_subs(cdata_subs, l)

def attr_update(self,**updates):
        self.__dict__.update(updates)

class Category:
    def __init__(self,label,enabled=1,desc=None,cat_id=None):
        attr_update(self,label=label,desc=desc,enabled=enabled,cat_id=cat_id,icons={})

    def update(self,enabled=None,desc=None,cat_id=None,icon=None, url=None):
        if icon:
            self.icons[icon]=url
        if desc and desc != '':
            self.desc=desc
        if cat_id:
            self.cat_id=cat_id
        if enabled:
            self.enabled=enabled

class Poi:
    def __init__(self,label,lat=None,lon=None,desc=None,categ=None):
        attr_update(self,label=label,desc=desc,lat=lat,lon=lon,categ=categ.lower())

class Poi2Mapper:
    def __init__ (self,src,dest_db,input_db):
        attr_update(self,src=src,input_db=input_db,categories={},styles={},icons={},pois=[])
        if dest_db:
            self.base=os.path.splitext(dest_db)[0]
            self.dest_db=dest_db
        else:
            self.base=os.path.splitext(src[0])[0]
            self.dest_db=self.base+'.db'
        # do the whole thing
        self.process_all()

    def categ_add_update(self,label=None,enabled=1,desc=None,cat_id=None,icon=None, url=None):
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

    def load_categ_from(self,path):
        if os.path.exists(path):
            cats_def=[[unicode(i.strip(),'utf-8') for i in l.split(',',4)] 
                        for l in open(path)]
            for d in cats_def:
                try:
                    (type,icon,categ,desc) = d + [None for i in range(len(d),4)]
                    l_d((type[-1] ,type,icon,categ,desc))
                    if type[-1] not in ('0','1') or not categ:
                        continue
                    self.categ_add_update(categ,int(type[-1]),icon=icon,desc=desc)
                except: pass

    def input_db_from(self,db_path):
        cat_ids={}
        if os.path.exists(db_path):
            db=sqlite3.connect(db_path)
            dbc=db.cursor()
            dbc.execute('select * from category')
            for (cat_id,label,desc,enabled) in dbc:
                icon=label+'.jpg'
                self.categ_add_update(label,enabled,icon=icon,desc=desc)
                cat_ids[cat_id]=label

            dbc.execute('select * from poi')
            for (poi_id,lat,lon,poi,desc, cat_id) in dbc:
                self.pois.append(Poi(poi,lat=lat,lon=lon,desc=desc,categ=cat_ids[cat_id]))
            db.close()

    def handleStyle(self,elm):
        url=None
        style_id=elm.getAttribute('id')
        l_d(style_id)
        icon=u'__%s__.jpg' % style_id
        if elm.getElementsByTagName("IconStyle") != []:
            try:
                url=elm.getElementsByTagName("href")[0].firstChild.data
                icon=re.sub('^.*/','',url)
            except: pass
        elif elm.getElementsByTagName("PolyStyle") != []:
            icon=u'__polygon__.jpg'
        elif elm.getElementsByTagName("LineStyle") != []:
            icon=u'__line__.jpg'
        return (style_id, self.categ_add_update(None,icon=icon,url=url))
                    
    def get_coords(self,elm):
        coords=elm.getElementsByTagName("coordinates")[0].firstChild.data.split()
        return [map(float,c.split(',')) for c in coords]
        
    def handlePlacemark(self,pm):
        point=pm.getElementsByTagName("Point")
        if point == []:
            return None
        coords=self.get_coords(point[0])
        (lon,lat)=coords[0][0:2]

        label=pm.getElementsByTagName("name")[0].firstChild.data
        style_id=pm.getElementsByTagName("styleUrl")[0].firstChild.data[1:]
        style=self.styles[style_id]
        l_d((label,style_id,style))
        if style.startswith('__') and style.endswith('__'):
            logging.warning(' No icon for "%s"' % label)
        desc=None
        try:
            desc_elm=pm.getElementsByTagName("description")[0]
            if desc_elm.firstChild:
                cdata=(desc_elm.firstChild.nodeType == self.doc.CDATA_SECTION_NODE)
                desc=desc_elm.firstChild.data
                if cdata:
                    desc=cdata_cleanup(desc)
        except IndexError:
            pass
        return Poi(label,lat=lat,lon=lon,desc=desc,categ=self.styles[style_id])

    def write_aux(self):
        cat_list=['# enabled, icon, category, desc']
        icon_urls=[]
        icon_aliases=["ln -s '%s.db' 'poi.db'"  % self.base]
        for (c_key,c) in self.categories.iteritems():
            for i_key in c.icons:
                cat_list.append('# %i, %s, %s%s' % (c.enabled,i_key,c.label,((', '+c.desc) if c.desc else '')))
                url=c.icons[i_key]
                if url:
                    icon_urls.append("wget -nc '%s'" % url)
                if c_key+'.jpg' != i_key:
                    icon_aliases.append("ln -s '%s' '%s'" % (i_key, c_key+'.jpg'))
        f=open(self.base+'.sh.gen','w')
        for ls in [cat_list,icon_urls,icon_aliases]:
            for s in ls:
                print >>f, s
            print >>f

    def proc_category(self,c):
        self.dbc.execute('insert into category (label, desc, enabled) values (?,?,?);',
            (c.label,c.desc,c.enabled))
        c.update(cat_id=self.dbc.lastrowid)
            
    def proc_poi(self,p):
        self.dbc.execute('insert into poi (lat, lon, label, desc, cat_id) values (?,?,?,?,?);',
            (p.lat,p.lon,p.label,p.desc,self.categories[p.categ].cat_id))

    def process_all(self):
        self.load_categ_from(self.base+'.sh')
        if self.input_db:
            map(self.input_db_from, self.input_db)

        for kml in self.src:
            self.doc = xml.dom.minidom.parse(kml)
            self.name=[n for n in self.doc.getElementsByTagName("Document")[0].childNodes 
                        if n.nodeType == n.ELEMENT_NODE and n.tagName == 'name'][0].firstChild.data

            self.styles=dict(map(self.handleStyle,self.doc.getElementsByTagName("Style")))
            self.pois+=filter(None,map(self.handlePlacemark,self.doc.getElementsByTagName("Placemark")))

        if os.path.exists(self.dest_db):
            os.remove(self.dest_db)
        self.db=sqlite3.connect(self.dest_db)
        self.dbc = self.db.cursor()
        try:
            self.db.execute ('''
                create table category (cat_id integer PRIMARY KEY, label text, desc text, enabled integer);
                ''')
            self.db.execute ('''
                create table poi (poi_id integer PRIMARY KEY, lat real, lon real, label text, desc text, cat_id integer);
                ''')
        except:
            pass

        map(self.proc_category,self.categories.itervalues())
        map(self.proc_poi,self.pois)
        self.db.commit()
        self.db.close()
        self.write_aux()

if __name__=='__main__':
    parser = optparse.OptionParser(
        usage="usage: %prog [-i <input_db>]... [-o <output_db>] [<kml_file>]...",
        description="makes maemo-mapper POI db from a kml file(s)")
    parser.add_option("-d", "--debug", action="store_true", dest="debug")
    parser.add_option("-q", "--quiet", action="store_true", dest="quiet")
    parser.add_option("-o", "--output",dest="dest_db", 
                      type="string",help="output POIs db file")
    parser.add_option("-i", "--input-db", action="append",dest="input_db", 
                      type="string",help="read POIs from another db")

    (options, args) = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG if options.debug else 
        (logging.ERROR if options.quiet else logging.INFO))

    if args == [] and not options.input_db:
        raise Exception("No source specified")

    Poi2Mapper(args,options.dest_db,options.input_db)

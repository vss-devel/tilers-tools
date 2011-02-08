#!/usr/bin/env python
# -*- coding: utf-8 -*-

# 2011-01-27 10:28:47 

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

from __future__ import with_statement
from __future__ import print_function

import sys
import os
import os.path
import logging
from subprocess import *
import itertools
import re
#from optparse import OptionParser

try:
    import multiprocessing # available in python 2.6 and above

    class KeyboardInterruptError(Exception): 
        pass
except:
    multiprocessing=None

def parallel_map(func,iterable):
    if not multiprocessing or len(iterable) < 2:
        return map(func,iterable)
    else:
        # map in parallel
        mp_pool = multiprocessing.Pool() # multiprocessing pool
        res=mp_pool.map(func,iterable)
        # wait for threads to finish
        mp_pool.close()
        mp_pool.join()
    return res

def ld(*parms):
    logging.debug(' '.join(itertools.imap(repr,parms)))

def ld_nothing(*parms):
    return

def pf(*parms,**kparms):
    end=kparms['end'] if 'end' in kparms else '\n'
    sys.stdout.write(' '.join(itertools.imap(str,parms))+end)
    sys.stdout.flush()

def pf_nothing(*parms,**kparms):
    return

def flatten(two_level_list): 
    return list(itertools.chain(*two_level_list))

try:
    import win32pipe 
except:
    win32pipe=None

def if_set(x,default=None):
    return x if x is not None else default

def command(params,child_in=None):
    cmd_str=' '.join(('"%s"' % i if ' ' in i else i for i in params))
    ld('>',cmd_str,child_in)
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
    ld('<',child_out,child_err)
    return child_out

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
    
def re_sub_file(fname, substitutions):
    'stream edit file using reg exp substitution list supplied'
    with open(fname+'.new', 'w') as out:
        for l in open(fname, 'r'):
            for (pattern,repl) in substitutions:
                l=re.sub(pattern,repl,string=l)
            out.write(l)
    shutil.move(fname+'.new',fname)


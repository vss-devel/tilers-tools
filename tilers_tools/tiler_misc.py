#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
# Copyright (c) 2011-2013 Vadim Shlyakhov
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


from tiler_functions import *
from tiler_backend import *

#############################

class GenericMap(Pyramid):
    'full profile options are to be specified'
#############################
    profile = 'generic'
    defaul_ext = '.generic'

    def __init__(self, src=None, dest=None, options=None):
        super(GenericMap, self).__init__(src, dest, options)

        self.srs = self.options.t_srs
        assert self.proj_srs, 'Target SRS is not specified'
        self.zoom0_tiles = map(int, self.options.zoom0_tiles.split(','))
        self.tile_dim = tuple(map(int, self.options.tile_dim.split(',')))
#
profile_map.append(GenericMap)
#

#~ class Yandex(Pyramid):
    #~ 'Yandex Maps (WGS 84 / World Mercator, epsg:3395)'
#~ ##############################
    #~ profile = 'yandex'
    #~ defaul_ext = '.yandex'
    #~ srs = '+proj=merc +datum=WGS84 +ellps=WGS84'
#~ #
#~ profile_map.append(Yandex)
#~ #


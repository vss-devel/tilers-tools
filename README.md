A few scripts for creating and handling a tile sets from digital raster maps. The scripts are based on GDAL tools.

Download from https://github.com/vadp/tilers-tools/
----
`tiler.py` \-- converts a [GDAL](http://www.gdal.org/)-compatible map file (dataset) into a set of zoom-levelled tile directories (a pyramid). A few output pyramid structure/projections (profiles) are supported: compatible with Google Maps (native and TMS-comatible), Google Earth and generic. The script is relatively fast, especially when processing a paletted source in the "draft mode" or rendering a few datasets simultaneously. It is also less picky with dataset formats and projections. In particular, it can cope with maps crossing the 180° meridian.

`tiler.py` should read from standard GDAL datasets, but for some of them it uses `map2gdal.py` internally to implements more accurate metadata translation: BSB (`.kap`), Geo/Nos (`.geo`), Ozi (`.map`) and KML with ground overlays (raster images). GDAL itself has some support for these formats, but the script extracts somewhat more geographical data from these source formats. For example, for BSB charts it supports more data and projections than GDAL does natively (currently the only datum for BSB charts is WGS84), it also makes use of BSB's DTM northing/easting data.

`tiles_merge.py` \-- merges a few separate tile sets created by `tiler.py` into a single one, so to cover a larger area and/or more zoom levels.

`tiles_convert.py` \-- converts tile sets between a few tile set structure and tile image formats.

Some less related scripts:

`ozf_decoder.py` \-- converts `ozf2` and `ozfx3` files into `tiff` format. If locally installed GDAL tools do not have an ozf reader built in, this one would help. As an advantage this script is also able to cope with some broken files.

`hdr_pcx_merge.py` \-- assembles a set of HDR-PCX chart tiles into a single PNG file;

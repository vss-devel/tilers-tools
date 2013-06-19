A few scripts for creating and handling a tile sets from digital raster maps. The scripts are based on GDAL tools.

To get the scripts go to https://code.google.com/p/gdal-tiler/source/browse/#hg%2Fgdal-tiler
----
 * `gdal-tiler.py` -- creates a tile set tree directory from a GDAL dataset;

 * `tiles-merge.py` -- sequentially merges a few tile sets in a single one to cover the area required;
 * `tiles2mapper.py` -- converts tile sets between a different tile structures: TMS, Google map-compatible (maemo mappero), SASPlanet cache, maemo-mapper sqlite3 and gmdb databases;

 * `tiles-opt.py` -- optimizes png tiles into a palleted form using pngnq tool;
 * `tiles-scale.py`

 * `bsb2gdal.py` -- creates geo-referenced GDAL .vrt file from BSB chart;
 * `ozi2gdal.py` -- creates geo-referenced GDAL .vrt file from Ozi map;

 * `ozf-decoder.py` -- converts .ozf2 or .ozfx3 file into .tiff
 * `hdr-pcx-merge.py` -- converts hdr-pcx chart image into .png

 * `poi2mapper.py`-- converts POI downloaded from the Google My Maps into maemo-mapper format

<wiki:comment>
 * gdal2kmz.py
 * kml2gdal.py
 * kml2gpx.sh
 * mk-merge-order.sh
 * vrt-bsb-cut.sh
</wiki:comment>


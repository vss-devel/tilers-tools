# Prerequisites for Linux

The prerequisites: Python (2.6 - 2.7), GDAL's Python bindings (python-gdal), Python imaging library (PIL), GDAL binaries (`gdal-bin` 1.7+) and, optionally, `pngnq` utility.

You may need to upgrade GDAL binaries to version 1.7+. For some Ubuntu distributions you may get it from the Ubuntu GIS Repository (see <http://trac.osgeo.org/ubuntugis/wiki/UbuntuGISRepository> ):

    sudo apt-get install python-software-properties
    sudo add-apt-repository ppa:ubuntugis/ubuntugis-unstable

or for stable PPA:

    sudo add-apt-repository ppa:ubuntugis/ppa

For other distributions, perhaps it's worth checking <http://trac.osgeo.org/gdal/wiki/DownloadingGdalBinaries>.

Then take `tilers-tools` from the Downloads page and put them into some directory. Make sure `tilers-tools` directory is included into the `PATH` environment variable.

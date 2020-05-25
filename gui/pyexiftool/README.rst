PyExifTool
==========

PyExifTool is a Python library to communicate with an instance of Phil
Harvey's excellent ExifTool_ command-line application.  The library
provides the class ``exiftool.ExifTool`` that runs the command-line
tool in batch mode and features methods to send commands to that
program, including methods to extract meta-information from one or
more image files.  Since ``exiftool`` is run in batch mode, only a
single instance needs to be launched and can be reused for many
queries.  This is much more efficient than launching a separate
process for every single query.

.. _ExifTool: http://www.sno.phy.queensu.ca/~phil/exiftool/

Getting PyExifTool
------------------

The source code can be checked out from the github repository with

::

    git clone git://github.com/smarnach/pyexiftool.git

Alternatively, you can download a tarball_.  There haven't been any
releases yet.

.. _tarball: https://github.com/smarnach/pyexiftool/tarball/master

Installation
------------

PyExifTool runs on Python 2.6 and above, including 3.x.  It has been
tested on Windows and Linux, and probably also runs on other Unix-like
platforms.

You need an installation of the ``exiftool`` command-line tool.  The
code has been tested with version 8.60, but should work with version
8.40 or above (which was the first production version of exiftool
featuring the ``-stay_open`` option for batch mode).

PyExifTool currently only consists of a single module, so you can
simply copy or link this module to a place where Python finds it, or
you can call

::

    python setup.py install [--user|--prefix=<installation-prefix]

to automatically install that module.

Documentation
-------------

The documentation is available at
http://smarnach.github.com/pyexiftool/.

Licence
-------

PyExifTool is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the licence, or
(at your option) any later version, or the BSD licence.

PyExifTool is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See COPYING.GPL or COPYING.BSD for more details.

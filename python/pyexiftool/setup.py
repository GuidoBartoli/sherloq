# PyExifTool <http://github.com/smarnach/pyexiftool>
# Copyright 2012 Sven Marnach

# This file is part of PyExifTool.
#
# PyExifTool is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the licence, or
# (at your option) any later version, or the BSD licence.
#
# PyExifTool is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# See COPYING.GPL or COPYING.BSD for more details.

from distutils.core import setup

setup(name="PyExifTool",
      version="0.1",
      description="Python wrapper for exiftool",
      license="GPLv3+",
      author="Sven Marnach",
      author_email="sven@marnach.net",
      url="http://github.com/smarnach/pyexiftool",
      classifiers=[
          "Development Status :: 3 - Alpha",
          "Intended Audience :: Developers",
          "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
          "Programming Language :: Python :: 2.6",
          "Programming Language :: Python :: 2.7",
          "Programming Language :: Python :: 3",
          "Topic :: Multimedia"],
      py_modules=["exiftool"])

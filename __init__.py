# -*- coding: utf-8; -*-
# stacker, Python module for stacking of interferometric data.
# Copyright (C) 2014  Lukas Lindroos
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor,
# Boston, MA 02110-1301, USA.

"""
    Library to stack interferometric images.
"""

import math
import os
from ctypes import cdll

__author__ = 'Lukas Lindroos'
__copyright__ = 'Copyright 2014'
__license__ = 'GPL'
__version__ = '1.0.3'
__maintainer__ = 'Lukas Lindroos'
__email__ = 'lindroos@chalmers.se'

PB_CONST = 0
PB_MS = 1
PB_FITS = 2

FILETYPENAME = {}
FILE_TYPE_NONE = 0
FILETYPENAME[FILE_TYPE_NONE] = 'none'
FILE_TYPE_MS = 1
FILETYPENAME[FILE_TYPE_MS] = 'ms'
FILE_TYPE_FITS = 2
FILETYPENAME[FILE_TYPE_FITS] = 'fits'

MS_DATACOLUMN_DATA = 1
MS_MODELCOLUMN_DATA = 2


clib_path = os.path.join(os.path.abspath(__path__[0]),
                         'stacker_clib')


def _casa_svnversion(casapath):
    import re
    import os.path
    try:
        casapyinfofile = open(os.path.join(casapath, 'casapyinfo'))
        svnversion = None
        for line in casapyinfofile:
            match = re.match('SVNVERSION="([0-9]*)"', line)
            if match:
                svnversion = match.groups()[0]
    except IOError:
        casaconfigfile = open(os.path.join(casapath, 'bin', 'casa-config'))
        for line in casaconfigfile:
            match = re.match('.*\$revision="([0-9]*)";.*', line)
            if match:
                svnversion = match.groups()[0]

    if svnversion is None:
        raise IOError('Can not find casa version.')

    return svnversion

try:
    # If this was imported from casa we need to ensure that the version
    # of stacker is compatible with the version of casa.
    from taskinit import casa
    _svnversion = _casa_svnversion(casa['dirs']['root'])
    _libpath = os.path.join(clib_path, 'libstacker-r{0}.so'.format(_svnversion))
    in_casapy = True
except ImportError:
    # We are in a pure python session.
    # Not that there is anything wrong with that.
    _libname = 'libstacker.so'
    _libpath = os.path.join(clib_path, _libname)
    in_casapy = False

import os

if in_casapy and _svnversion is not None and os.access(_libpath, os.F_OK):
    print('Loading stacking library for casapy svn revision {0}'.format(
        _svnversion))
elif in_casapy:
    print('warning, no precompiled library compatible with your version of casa exists.')
    print('It is recommended to recompile stacker for your casa version.')
    import re
    stacker_clib_ls = os.listdir(clib_path)
    stackerlibs = [((re.match('libstacker-r([0-9]*).so', f).group(1)), f)
                   for f in stacker_clib_ls
                   if re.match('libstacker-r[0-9]*.so', f)]
    print(stackerlibs)
    vdiff = [(abs(int(_svnversion)-int(v)), lib, v)
             for (v, lib) in stackerlibs]
    vdiff.sort()
    print(vdiff)
    print('Trying to use svn revision {0}.'.format(vdiff[0][2]))
    _libpath = os.path.join(clib_path, vdiff[0][1])


libstacker = cdll.LoadLibrary(_libpath)


class CoordList(list):
    """
        Extended list to contain list of coordinates.

    """
    def __init__(self, imagenames=[], coord_type='physical', unit='rad'):
        """
        Requires an image list in case of pixel coordinates to work properly.
        """

        super(CoordList, self).__init__()

        if isinstance(imagenames, str):
            imagenames = [imagenames]

        self.coords = []
        self.imagenames = imagenames
        self.coord_type = coord_type
        self.unit = unit

    def __getitem__(self, i):
        return self.coords[i]

    def __setitem__(self, i, x):
        self.coords[i] = x

    def append(self, x):
        self.coords.append(x)

    def __len__(self):
        return len(self.coords)

    def __iter__(self):
        for x in self.coords:
            yield x

    def __getslice__(self, i, j):
        new_a = CoordList(self.imagenames, self.coord_type, self.unit)
        new_a.coords = self.coords.__getslice__(i, j)
        return new_a

    def __repr__(self):
        ret = []
        for x in self.coords:
            ret.append('(' + x.__str__() + ')')
        return '\n'.join(ret)

    def __str__(self):
        ret = []
        for x in self.coords:
            ret.append('(' + x.__str__() + ')')
        return '\n'.join(ret)
#         return '{0}, {1}'.format(self.x, self.y)


class Coord:
    """
        Describes a stacking position.

        Class used internally to represent coordinates. May describe a
        physical coordinate or a pixel coordinate.
    """

    def __init__(self, x, y, weight=1., image=0):
        """
            Create a coordinate. A pixel coordinate should always specify
            to which image it belongs. Physical coordinates should be in
            J2000 radians.
        """
        self.x = x
        self.y = y
        self.weight = weight
        self.image = image

    def __str__(self):
        return '{0}, {1}'.format(self.x, self.y)


def readCoords(coordfile, unit='deg'):
    """
        Reads a coordinate file from disk and produces a list.

        coordfile:
            Path to coordinate file. A file in csv format. x and y should
            be in J2000 . A weight may be added in third column
            to weight positions for stacking. If no weight are wanted
            put no third column in coordinate file.
        unit:
            Unit of input coordinates. Allows two values, 'deg' and 'rad'.
    """
    import csv

    coordreader = csv.reader(open(coordfile, 'rb'), delimiter=',')
    coords = CoordList()
    for row in coordreader:
        x = float(row[0])
        y = float(row[1])

        if unit == 'deg':
            x = math.pi/180.*x
            y = math.pi/180.*y

        if len(row) > 2:
            weight = float(row[2])
        else:
            weight = 1.

        if x > 2*math.pi:
            x -= 2*math.pi
        if y > math.pi:
            y -= 2*math.pi

        coords.append(Coord(x, y, weight))

    return coords


def _checkfile(filename, datacolumn):
    import re
    # Currently this supports only ms files
    # As such there is no reason to check filetype.
    # If it cannot be opened as ms it will not be supported.
#     if re.match('^.*[mM][sS]/*$', filename) is not None:
    try:
        from taskinit import ms
        ms.open(filename)
        ms.done()
    except ImportError:
        # This probably means that it was run from a pure python session.
        # We will relegate any checks that it is a valid ms file to the 
        # stacker.
        if not os.access(filename, os.F_OK):
            raise RuntimeError('Could not find data file "{}".'.format(
                filename))

    filename = filename
    filetype = FILE_TYPE_MS
    fileoptions = 0
    if datacolumn == 'data':
        fileoptions = MS_DATACOLUMN_DATA
    elif datacolumn == 'model' or datacolumn == 'model_data':
        fileoptions = MS_MODELCOLUMN_DATA
#     elif re.match('^.*[fF][iI][tT][sS]$', filename) is not None:
#         raise NotImplementedError('FITS format is currently not supported.')
    return filetype, filename, fileoptions


def coordsTocl(name, flux, coords):
    from taskinit import cl, qa

    flux = qa.quantity(flux)
    cl.done()
    cl.rename(name)
    for coord in coords:
        clpars = {}
        clpars['flux'] = -flux['value']
        clpars['fluxunit'] = flux['unit']
        clpars['dir'] = ['J2000', str(coord.x)+'rad', str(coord.y)+'rad']
        clpars['shape'] = 'point'
        cl.addcomponent(**clpars)
    cl.done()


def randomCoords(imagenames, ncoords=10):
    import random
    from taskinit import ia, qa

    xmin, xmax = [], []
    ymin, ymax = [], []
    for image in imagenames:
        ia.open(image)
        print image, ia.boundingbox()
        trc = ia.boundingbox()['trcf'].split(', ')
        blc = ia.boundingbox()['blcf'].split(', ')
        xmin.append(qa.convert(qa.quantity(trc[0]), 'rad')['value'])
        xmax.append(qa.convert(qa.quantity(blc[0]), 'rad')['value'])
        ymin.append(qa.convert(qa.quantity(blc[1]), 'rad')['value'])
        ymax.append(qa.convert(qa.quantity(trc[1]), 'rad')['value'])
        ia.done()

    randomcoords = CoordList(imagenames)
    for i in range(ncoords):
        imageid = random.randint(0, len(imagenames)-1)
        x = random.uniform(xmin[imageid], xmax[imageid])
        y = random.uniform(ymin[imageid], ymax[imageid])
        c = Coord(x, y, 1.0)
        randomcoords.append(c)

    return randomcoords


def randomizeCoords(coords, beam):
    import random
    import math

    randomcoords = CoordList(coords.imagenames, coords.coord_type,
                             unit=coords.unit)

    for coord in coords:
        dr = random.uniform(beam, 5*beam)
        dphi = random.uniform(0, 2*math.pi)
        x = coord.x + dr*math.cos(dphi)
        y = coord.y + dr*math.sin(dphi)
        randomcoords.append(Coord(x, y, coord.weight, coord.image))

    return randomcoords


def _getPixelCoords1ImSimpleProj(coords, imagename):
    from taskinit import ia
    from interval import interval

    ia.open(imagename)
    cs = ia.coordsys()
    imshape = ia.shape()
    ia.done()

    pixcoords = []
    for coord in coords:
        p = cs.convert(coordin=[coord.x, coord.y, 0, 0], absin=[True]*4,
                       unitsin=[coords.unit, coords.unit, 'pix', 'pix'],
                       absout=[True]*4, unitsout=['pix']*4)
        x = p[0]
        y = p[1]

        if x in interval[0, imshape[0]-1] and y in interval[0., imshape[1]-1]:
            c = Coord(x, y)
            try:
                c.index = coord.index
            except AttributeError:
                pass
            pixcoords.append(c)

    return pixcoords


def _getPixelCoords1Im(coords, imagename):
    from interval import interval
    import math

    try:
        from taskinit import ia
        ia.open(imagename)
        cs = ia.coordsys()
        Nx = ia.shape()[0]
        Ny = ia.shape()[1]
        ia.done()
        x0 = cs.referencevalue()['numeric'][0]
        y0 = cs.referencevalue()['numeric'][1]
        x_pix_ref = cs.referencepixel()['numeric'][0]
        y_pix_ref = cs.referencepixel()['numeric'][1]
        x_pix_inc = cs.increment()['numeric'][0]
        y_pix_inc = cs.increment()['numeric'][1]

# If we fail to load ia, we will use pyrap instead.
# This probably means stacker was loaded from outside casapy.
    except ImportError:
        from pyrap.images import image
        im = image(imagename)
        cs = im.coordinates().get_coordinate('direction')
        dir_axis_index = im.coordinates().get_axes().index(cs.get_axes())
        imshape = im.shape()
        try:
            x_axis_index = cs.get_axes().index('Right Ascension')
        except ValueError:
            raise ValueError('Could not find direction coordinate: '\
                              'RightAscension')
        try:
            y_axis_index = cs.get_axes().index('Declination')
        except ValueError:
            raise ValueError('Could not find direction coordinate: '\
                              'Declination')



        Nx = im.shape()[dir_axis_index+x_axis_index]
        Ny = im.shape()[dir_axis_index+y_axis_index]
        x0 = cs.get_referencevalue()[x_axis_index]
        y0 = cs.get_referencevalue()[y_axis_index]
        x_pix_ref = cs.get_referencepixel()[x_axis_index]
        y_pix_ref = cs.get_referencepixel()[y_axis_index]
        x_pix_inc = cs.get_increment()[x_axis_index]
        y_pix_inc = cs.get_increment()[y_axis_index]

    pixcoords = []
    for coord in coords:
        dx = (coord.x - x0)*math.cos(coord.y)
        dy = math.asin(math.sin(coord.y)/math.cos(dx)) - y0
        x = dx/x_pix_inc+x_pix_ref
        y = dy/y_pix_inc+y_pix_ref

        if x in interval[0, Nx-1] and y in interval[0., Ny-1]:
#             pixcoords.append(Coord(x,y, coord.weight))
            c = Coord(x, y, coord.weight)

            try:
                c.index = coord.index
            except AttributeError:
                pass

            pixcoords.append(c)

    return pixcoords


def make_pbfile(vis, pbfile):
    from taskinit import im, ms, ia, qa, tb
    import numpy as np
    from scipy.constants import c

    ms.open(vis)
    fields = ms.range('field_id')['field_id']
    ms.done()
    im.open(vis)
    im.selectvis(field=fields[0])
    ms.open(vis)
    freq = np.mean(ms.range('chan_freq')['chan_freq'])
    phase_dir = ms.range('phase_dir')['phase_dir']['direction']
    ms.done()

    phase_dir = phase_dir[0][0], phase_dir[1][0]
    phase_dir = [qa.formxxx(str(phase_dir[0])+'rad', format='hms'),
                 qa.formxxx(str(phase_dir[1])+'rad', format='dms')]
    phase_dir = 'J2000 '+' '.join(phase_dir)

    tb.open(vis+'/ANTENNA/')
    dishdia = np.min(tb.getcol('DISH_DIAMETER'))
    tb.done()

    # pb of 512 pix cover pb down to 0.001
    # ensure largest pixel to pixel var to .01
    minpb = 0.001
    nx = 512
    cellconv = (nx*np.sqrt(np.log(2)/np.log(1/minpb)))**-1

    beam = c/freq/dishdia
    cell = {}
    cell['value'] = beam*cellconv
    cell['unit'] = 'rad'

#     nx = int(3*3e8/freq/dishdia*1.22*180/
#              math.pi*3600/qa.convert(advise['cell'],
#              'arcsec')['value'])
    # Chosen as to be 3 times fwhm of primary beam,
    # should include up to approximately .01 of peak flux

    im.defineimage(nx=nx, ny=nx, cellx=cell, celly=cell, phasecenter=phase_dir)
    im.setvp(dovp=True)
    im.makeimage(type='pb', image=pbfile)
    im.done()
    ia.open(pbfile)
    cs = ia.coordsys()
    cs.setreferencevalue(type='direction', value=[0., 0.])
    ia.setcoordsys(cs.torecord())
    ia.maskhandler('delete', 'mask0')
    ia.done()


def getPixelCoords(coords, imagenames):
    """
        Creates pixel coordinate list from a physical coordinate
        list and a list of images.
    """

    pixcoords = CoordList(imagenames, 'pixel', unit='pix')

    for (i, imagename) in enumerate(pixcoords.imagenames):
        for coord in _getPixelCoords1Im(coords, imagename):
            coord.image = i
            pixcoords.append(coord)

    return pixcoords

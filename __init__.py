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
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

import math
import os
from ctypes import c_double, POINTER, c_char_p, cdll, c_int

clib_path = os.path.join(os.path.abspath(__path__[0]),
                         'stacker_clib')
libstacker = cdll.LoadLibrary(os.path.join(clib_path, 'libstacker.so'))

"""
    Library to stack interferometric images.
"""

class CoordList(list):
    """
        Extended list to contain list of coordinates.

    """
    def __init__(self, imagenames = [], coord_type = 'physical', unit = 'rad'):
        """ Requires an image list in case of pixel coordinates to work properly. """

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


#     def __getitem__(self, x):
#         new = CoordList(self.imagenames, self.coord_type, self.unit)
# #         new.coords = self.coords#.__getitem__(x)
#         return new


    def __len__(self):
        return len(self.coords)


    def __iter__(self):
        for x in self.coords:
            yield x


    def __getslice__(self, i, j):
        new_a = CoordList(self.imagenames, self.coord_type, self.unit)
        new_a.coords = self.coords.__getslice__(i,j)
        return new_a


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

    def __init__(self, x, y, weight = 1., image = 0):
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


def readCoords(coordfile, unit = 'deg'):
    """
        Reads a coordinate file from disk and produces a list. 
    
        coordfile: 
            Path to coordinate file. A file in csv format. x and y should
            be in J2000 . A weight may be added in third column
            to weight positions for stacking.
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

        if x > math.pi:
            x -= 2*math.pi
        if y > math.pi:
            y -= 2*math.pi

        coords.append(Coord(x,y,weight))

    return coords


def calculateSigma2Weights(coords, imagenames):
    pass


def coordsTocl(name, flux, coords):
    from taskinit import cl,qa

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


def randomCoords(imagenames, ncoords = 10):
    import random
    import math
    import copy
    from taskinit import ia,qa

    
    xmin,xmax = [],[]
    ymin,ymax = [],[]
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
    import copy

    randomcoords = CoordList(coords.imagenames, coords.coord_type, 
                             unit = coords.unit)

    for coord in coords:
        dr = random.uniform(beam, 5*beam)
        dphi = random.uniform(0,2*math.pi)
        x = coord.x + dr*math.cos(dphi)
        y = coord.y + dr*math.sin(dphi)
        randomcoords.append(Coord(x,y, coord.weight, coord.image))


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
        p = cs.convert(coordin=[coords.x, coord.y,0,0], absin = [True]*4, 
                unitsin=[coords.unit,coords.unit, 'pix', 'pix'], absout = [True]*4, unitsout=['pix']*4)
        x = p[0]
        y = p[1]

        if x in interval[0,imshape[0]-1] and y in interval[0.,imshape[1]-1]:
            c = Coord(x,y)
            try:
                c.index = coord.index
            except AttributeError:
                pass
            pixcoords.append(c)

    return pixcoords


def _getPixelCoords1Im(coords, imagename):
    from taskinit import ia
    from interval import interval
    import math

    ia.open(imagename)
    cs = ia.coordsys()
    imshape = ia.shape()
    ia.done()

    pixcoords = []
    for coord in coords:
        tmpx = coord.x
        if tmpx < 0:
            tmpx += 2*math.pi
        x0 = cs.referencevalue()['numeric'][0]
        y0 = cs.referencevalue()['numeric'][1]
        dx = (tmpx - x0)*math.cos(coord.y)
        dy = math.asin(math.sin(coord.y)/math.cos(dx)) - y0
        x_pix_ref = cs.referencepixel()['numeric'][0]
        y_pix_ref = cs.referencepixel()['numeric'][1]
        x_pix_inc = cs.increment()['numeric'][0]
        y_pix_inc = cs.increment()['numeric'][1]
        x = dx/x_pix_inc+x_pix_ref
        y = dy/y_pix_inc+y_pix_ref

#         p = cs.convert(coordin=[coord.x, coord.y,0,0], absin = [True]*4, 
#                 unitsin=[coords.unit,coords.unit, 'pix', 'pix'], absout = [True]*4, unitsout=['pix']*4)
#         x = p[0]
#         y = p[1]

        if x in interval[0,imshape[0]-1] and y in interval[0.,imshape[1]-1]:
#             pixcoords.append(Coord(x,y, coord.weight))
            c = Coord(x,y, coord.weight)

            try:
                c.index = coord.index
            except AttributeError:
                pass

            pixcoords.append(c)

    return pixcoords


def make_pbfile(vis, pbfile):
    from taskinit import im,ms,ia,qa,tb
    import numpy as np

    im.open(vis)
    im.selectvis(field='0')
    advise = im.advise()
    ms.open(vis)
    freq = np.mean(ms.range('chan_freq')['chan_freq'])
    ms.done()

    tb.open(vis+'/ANTENNA/')
    dishdia = np.min(tb.getcol('DISH_DIAMETER'))
    tb.done()
    nx = int(3*3e8/freq/dishdia*1.22*180/math.pi*3600/qa.convert(advise['cell'], 'arcsec')['value'])
    # Chosen as to be 3 times fwhm of primary beam, should include up to approximately .01 of peak flux

    im.defineimage(nx =nx, ny=nx, cellx=advise['cell'], celly=advise['cell'], phasecenter = advise['phasecenter'])
    im.setvp(dovp=True)
    im.makeimage(type='pb', image=pbfile)
    ia.open(pbfile)
    cs = ia.coordsys()
    cs.setreferencevalue(type='direction', value=[0.,0.])
    ia.setcoordsys(cs.torecord())
    ia.maskhandler('delete', 'mask0')
    ia.done()


def calculate_pb_weight(coords, pbfile):
    pass


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




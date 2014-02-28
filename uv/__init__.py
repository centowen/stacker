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
from ctypes import c_double, POINTER, c_char_p, cdll, c_int
import numpy as np
import stacker
import os
from rmtables import rmtables

# lib = cdll.LoadLibrary('stacker/uv/libuvstack.so')
c_stack = stacker.libstacker.stack
c_stack.restype = c_double
c_stack.argtype = [c_char_p, c_char_p, c_char_p,
                   POINTER(c_double), POINTER(c_double), POINTER(c_double), 
                   c_int]
# c_stack.argtype = [c_char_p, c_char_p, c_char_p,
#                    POINTER(c_double), POINTER(c_double), POINTER(c_double), 
#                    c_int, c_int, c_double]


#  WARNING DO NOT USE THIS FUNCTION!!!!!!
# Only left for convience if a new function is to be implemented
# Will not produce proper results in current form (will probably crash)
def ___DANGER___stack_random(vis, npos, imagenames, outvis, 
                 flux = {'unit': 'Jy', 'value': 1.0}, 
                 stampsize = 32, psfmode='point',
                 model = 'mcmodel.cl'):

    from modsub_cl import modsub
    import stacker.image


    coords = stacker.randomCoords(imagenames, npos)
    if os.access(model, os.F_OK): rmtables(model)
    stacker.coordsTocl(model, flux, coords)
    modsub(model, vis, outvis)
    coords = stacker.image.calculate_sigma2_weights(coords, imagenames, stampsize)
    stack(coords, outvis, outvis, stampsize=stampsize)


def stack(coords, vis, outvis='', imagename='', cell = '1arcsec', stampsize = 32):
    import shutil
    import os

    if outvis != '':
        if not os.access(outvis, os.F_OK):
            shutil.copytree(vis, outvis)


# primary beam
    from taskinit import ms,tb,qa
    ms.open(vis)
    freq = int(np.mean(ms.range('chan_freq')['chan_freq'])/1e8)/10.
    ms.done()
    tb.open(vis+'/OBSERVATION')
    telescope = tb.getcol('TELESCOPE_NAME')[0]
    tb.done()
    pbfile = '{0}-{1}GHz.pb'.format(telescope, freq)
    if not os.access(pbfile, os.F_OK):
        stacker.make_pbfile(vis, pbfile)
    

    x = [p.x for p in coords]
    y = [p.y for p in coords]
    weight = [p.weight for p in coords]

    x = (c_double*len(x))(*x)
    y = (c_double*len(y))(*y)
    weight = (c_double*len(weight))(*weight)

    from taskinit import qa

    flux = c_stack(c_char_p(vis), c_char_p(outvis), c_char_p(pbfile), x, y, 
                   weight, c_int(len(coords)))

    if imagename != '':
        import clean
        import clearcal
        clearcal.clearcal(vis=outvis)
        clean.clean(vis=outvis, imagename = imagename, field='0', mode='mfs',
                    cell=cell, imsize=stampsize, weighting='natural')

    return flux


def noise(coords, vis, imagenames, weighting = 'sigma2', nrand = 50, stampsize = 32):
    import stacker
    import stacker.image
    from taskinit import ia, qa

    ia.open(imagenames[0])
    beam = qa.convert(ia.restoringbeam()['major'], 'rad')['value']
    ia.done()

    dist = []
    for i in range(nrand):
        random_coords = stacker.randomizeCoords(coords, beam=beam)
        print 'peti',len(random_coords), len(coords)
        if weighting == 'sigma2':
            random_coords = stacker.image.calculate_sigma2_weights( random_coords, imagenames, stampsize)
        dist.append(stack(random_coords, vis))

    return np.array(dist)

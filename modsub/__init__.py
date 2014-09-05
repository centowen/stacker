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

# lib = cdll.LoadLibrary('libmodsub.so')
lib = stacker.libstacker
c_modsub = lib.modsub
c_modsub.restype = c_int
c_modsub.argtype = [c_char_p, c_char_p, c_char_p,c_char_p]

def modsub(model, vis, outvis='', datacolumn='corrected'):
    import shutil
    import os

    if outvis != '':
        if not os.access(outvis, os.F_OK):
            shutil.copytree(vis, outvis)


# primary beam
    from taskinit import ms,tb,qa
    ms.open(vis)
    freq = int(np.mean(ms.range('chan_freq')['chan_freq'])/1e9*100)/100.
    ms.done()
    tb.open(vis+'/OBSERVATION')
    telescope = tb.getcol('TELESCOPE_NAME')[0]
    tb.done()
    pbfile = '{0}-{1}GHz.pb'.format(telescope, freq)

    if not os.access(pbfile, os.F_OK):
        stacker.make_pbfile(vis, pbfile)

    infiletype, infilename, infileoptions = stacker._checkfile(vis, datacolumn)
    outfiletype, outfilename, outfileoptions = stacker._checkfile(outvis, datacolumn)

    flux = c_modsub(infiletype, c_char_p(infilename), infileoptions,
                    outfiletype, c_char_p(outfilename), outfileoptions,
                    c_char_p(model), c_char_p(pbfile))

    return 0


def cl_from_im(image, clname=None, threshold=None):
    from taskinit import qa,ia, cl

    ia.open(image)
    data = ia.getregion()
    cs = ia.coordsys()
    ia.done()

    data = np.where(np.isnan(data), 0, data)
    if threshold is None:
        datanz = np.nonzero(data)
    else:
        datanz = np.nonzero(np.abs(data) > threshold)

    modelfluxes = data[datanz]
    modelpixels = np.array(datanz)
    modellist = cs.convertmany(modelpixels, 
                               unitsin=['pix', 'pix', 'pix', 'pix'], 
                               unitsout=['rad', 'rad', '', 'Hz'])

    cl.done()

    for i in range(modellist.shape[1]):
        x = qa.formxxx(str(modellist[0, i])+'rad', format='hms', prec=6)
        y = qa.formxxx(str(modellist[1, i])+'rad', format='dms', prec=6)
        pos = ' '.join(['J2000', x, y])
        freq = str(modellist[3, i])+'Hz'
        flux = modelfluxes[i]
        cl.addcomponent(flux=flux, fluxunit='Jy', dir=pos, freq=freq)

    if clname is not None:
        cl.rename(clname)
        cl.done()
    else:
        return cl

        

import math


skymap = []
data = []
oldimagenames = []
stampsize = 0
imagesizes = []


def calculate_sigma2_weights(coords, imagenames=[], stampsize=32, maxmaskradius=None):
    import stacker

    for i, coord in enumerate(coords):
        coord.index = i

    if coords.coord_type == 'physical':
        pixcoords = stacker.getPixelCoords(coords, imagenames)
    else:
        pixcoords = coords

    _allocate_buffers(pixcoords.imagenames, stampsize, len(pixcoords))
    _load_stack(pixcoords)
    pixcoords =  _calculate_sigma2_weights(pixcoords, maxmaskradius)

#     for i in range(len(coords)):
    for coord in coords:
        coord.weight = 0.
        for pixcoord in pixcoords:
            if coord.index == pixcoord.index and pixcoord.weight > coord.weight:
                coord.weight = pixcoord.weight

    return coords


def calculate_flux_weights(coords, imagenames=[]):
    import stacker

    for i, coord in enumerate(coords):
        coord.index = i

    if coords.coord_type == 'physical':
        pixcoords = stacker.getPixelCoords(coords, imagenames)
    else:
        pixcoords = stacker.CoordList(coords.imagenames, 'pixel', unit='pix')

# This works for now. But it is not good if more atrributes exists from somewhere else.
        for coord in coords:
            c = stacker.Coord(coord.x, coord.y, coord.weight)
            c.image = coord.image
            c.index = coord.index
            pixcoords.append(c)

    pixcoords =  _calculate_flux_weights(pixcoords)

    for coord in coords:
        coord.weight = 0.
        for pixcoord in pixcoords:
            if coord.index == pixcoord.index and pixcoord.weight > coord.weight:
                coord.weight = pixcoord.weight

#     for coord in pixcoords:
#         print coord.weight
# 
#     for coord in coords:
#         print coord.weight

    return coords


def stack(coords, outfile, stampsize = 32, imagenames= [], method = 'mean',
        weighting = 'sigma2', maxmaskradius=None, psfmode = 'point'):

    from ..interval import interval
    import stacker
    import os
    import shutil
    import numpy as np
    from taskinit import ia

    global skymap
    global data
#     global stampsize
    global oldimagenames

    if coords.coord_type == 'physical':
        coords = stacker.getPixelCoords(coords, imagenames)

# Important that len(coords) here is for the pixel coordinates, not physical!
    _allocate_buffers(coords.imagenames, stampsize, len(coords))

#	This is infact not neccessary in many case, just leave it like this for simplicity.
    ia.open(coords.imagenames[0])
    cs = ia.coordsys()
    outnchans = ia.boundingbox()['trc'][2]+1
    outnstokes = ia.boundingbox()['trc'][3]+1
    ia.done()


    for imagename in coords.imagenames:
            ia.open(imagename)
            if ia.shape()[2] != outnchans or ia.shape()[3] != outnstokes:
                    print('Channels do not match in all images!')
                    return
            ia.done()

    _load_stack(coords, psfmode)

    if method == 'mean' and weighting == 'sigma2':
        coords = _calculate_sigma2_weights(coords, maxmaskradius)
    elif method == 'mean' and weighting == 'sigma':
        coords = _calculate_sigma_weights(coords, maxmaskradius)

    print 'npos: {0}'.format(len([c.weight for c in coords if c.weight > 1e-6]))

    stacked_im  = _stack_stack(method, coords)
#     print 'peti: {0}'.format(coords[3].weight)


    _write_stacked_image(outfile, stacked_im, coords.imagenames[0])
    return stacked_im[int(stampsize/2+0.5), int(stampsize/2+0.5),0,0]

            
def noise(coords, nrandom = 50, imagenames=[], stampsize=32,
        method = 'mean', weighting = 'simga2', maxmaskradius=None,
        psfmode = 'point'):

    import stacker
    import numpy as np
    from taskinit import ia, qa

    ia.open(imagenames[0])
    beam = qa.convert(ia.restoringbeam()['major'], 'rad')['value']
    ia.done()

#     if coords.coord_type == 'physical':
#         coords = stacker.getPixelCoords(coords, imagenames)

    _allocate_buffers(imagenames, stampsize, len(coords)*len(imagenames))

    dist = []

#     coords = stacker.randomizeCoords(coords, beam=5)
    for i in range(nrandom):
        random_coords = stacker.randomizeCoords(coords, beam=beam)
        random_coords = stacker.getPixelCoords(random_coords, imagenames)
        _load_stack(random_coords, psfmode)

        if method == 'mean' and weighting == 'sigma2':
            random_coords = _calculate_sigma2_weights(random_coords, maxmaskradius)
        elif method == 'mean' and weighting == 'sigma':
            random_coords = _calculate_sigma_weights(random_coords, maxmaskradius)

        stacked_im  = _stack_stack(method, random_coords)

        dist.append(stacked_im[int(stampsize/2+0.5), int(stampsize/2+0.5),0,0])

    return dist


def getFlux(imagename):
    from taskinit import ia,rg
    ia.open(imagename)
    cs = ia.coordsys()
    x = int(cs.referencepixel()['numeric'][0])
    y = int(cs.referencepixel()['numeric'][1])
    ia.done()
    return float(ia.getregion(region=rg.box([x,y], [x,y])))


def _write_stacked_image(imagename, pixels, template_image):
    import os
    import shutil
    import numpy as np
    from taskinit import ia
    global stampsize
    if os.access(imagename, os.F_OK): shutil.rmtree(imagename)

    ia.open(template_image)
    beam = ia.restoringbeam()
    cs = ia.coordsys()
    ia.done()

    csnew = cs.copy()
    csnew.setreferencevalue([0.]*2, 'dir')
    csnew.setreferencepixel([int(stampsize/2+0.5)]*2, 'dir')
    ia.fromarray(imagename, pixels=pixels, csys = csnew.torecord())
    ia.open(imagename)
    ia.setrestoringbeam(beam=beam)
    ia.done()


def _calculate_sigma_weights(coords, maxmaskradius=None):
    import numpy as np
    from taskinit import ia
    global stampsize

    if coords.physical and coords.imagenames[0]:
        ia.open(coords.imagenames[0])
#         beam = qa.convert(ia.restoringbeam()['major'], 'rad')['value']
        if ia.restoringbeam() == {}:
            masksize = 10
        else:
            masksize = 2*np.abs(qa.convert(ia.restoringbeam()['major'],'rad')['value']/ ia.coordsys().increment()['numeric'][0])
        ia.done()

    X = np.arange(0, stampsize)-stampsize/2
    Y = np.arange(0, stampsize)-stampsize/2
    X,Y = np.meshgrid(X,Y)

    for i,coord in enumerate(coords):
        tmpdata = data[i,:,:,:,:]
        for j in range(tmpdata.shape[2]):
            for k in range(tmpdata.shape[3]):
                tmpdata[:,:,j,k]  = (tmpdata[:,:,j,k]*np.double( np.sqrt(X**2+Y**2)>masksize))
        sigma = np.std(tmpdata)
        if sigma == 0:
            coord.weight = 0.
        else:
            coord.weight = 1/sigma
    if maxmaskradius and maxmaskradius < masksize:
        masksize = maxmaskradius

    return coords


def _calculate_flux_weights(coords):
    from taskinit import ia
    import re

    fluxmap = []

    r = re.compile('(.*)\.image/*')
    for imagename in coords.imagenames:
        match = r.match(imagename)

        if match:
            fluximage = match.group(1)+'.flux'

            if ia.open(fluximage):
                fluxmap.append(ia.getregion())

        ia.done()

    for i,coord in enumerate(coords):
        coord.weight = (fluxmap[coord.image][int(coord.x+0.5), int(coord.y+0.5), 0, 0])**2
        
    return coords


def _calculate_sigma2_weights(coords, maxmaskradius=None):
    import numpy as np
    from taskinit import ia,qa
    global stampsize

    masksize=10
    if coords.coord_type == 'pixel' and coords.imagenames[0]:
        ia.open(coords.imagenames[0])
#         beam = qa.convert(ia.restoringbeam()['major'], 'rad')['value']
        if ia.restoringbeam() == {}:
            masksize = 10
        else:
            masksize = int(2*np.abs(qa.convert(ia.restoringbeam()['major'],'rad')['value']/ ia.coordsys().increment()['numeric'][0]))
        ia.done()
    if maxmaskradius and maxmaskradius < masksize:
        masksize = maxmaskradius
    print 'masksize set: {0}'.format(masksize)

    X = np.arange(0, stampsize)-int(stampsize/2)
    Y = np.arange(0, stampsize)-int(stampsize/2)
    X,Y = np.meshgrid(X,Y)

    for i,coord in enumerate(coords):
        tmpdata = np.copy(data[i,:,:,:,:])
        for j in range(tmpdata.shape[2]):
            for k in range(tmpdata.shape[3]):
                tmpdata[:,:,j,k]  = (tmpdata[:,:,j,k]*np.double( np.sqrt(X**2+Y**2)>masksize))
        sigma = np.std(tmpdata[np.nonzero(tmpdata)])
        if sigma == 0:
            coord.weight = 0.
        else:
            coord.weight = 1/sigma**2

    return coords


def _allocate_buffers( imagenames, new_stampsize, nstackpos):
    import numpy as np
    from taskinit import ia

    global skymap
    global data
    global oldimagenames
    global stampsize
    global imagesizes

    ia.open(imagenames[0])
    cs = ia.coordsys()
    outnchans = ia.boundingbox()['trc'][2]+1
    outnstokes = ia.boundingbox()['trc'][3]+1
    ia.done()
    
# To improve performance this module will keep buffers between run.
# This following code resets these buffers if they have grown obsolete.
    if oldimagenames == []:
            oldimagenames = imagenames

    if oldimagenames != imagenames:
            oldimagenames = imagenames
            skymap = []
            data = []

    if stampsize == 0:
            stampsize = new_stampsize

    elif stampsize != new_stampsize:
            stampsize = new_stampsize
            data = []
            skymap = []
    
    if not(data == []) and nstackpos != data.shape[0]:
        data = []

# If there is no data buffer create one.
# The data buffer is used to save the right stacking positions before stacking them.
# During stacking this is where the full stack will actually be saved.
    if data == []:
            data = np.zeros((nstackpos, new_stampsize, new_stampsize, outnstokes, outnchans))
    else:
            data = 0.*data
    
# If there is no skymap buffer create one.
# This is the data that is most important to buffer.
# Reading a skymap from disk is a time consuming task and we don't want to do this too much.
    if skymap == []:
        for imagename in imagenames:
            ia.open(imagename)
            skymap.append(ia.getregion())
            ia.done()

    imagesizes = []
    for imagename in imagenames:
        ia.open(imagename)
        imagesizes.append((ia.shape()[0], ia.shape()[1]))
        ia.done()


def _load_stack(coords, psfmode='point'):
    from ..interval import interval

    global skymap
    global data
    global stampsize
    global imagesizes

    if len(coords) > data.shape[0]:
        _allocate_buffers(coords.imagenames, stampsize, len(coords))

    for (i,coord) in enumerate(coords):
        blcx = int(coord.x - stampsize/2 + 0.5)
        blcy = int(coord.y - stampsize/2 + 0.5)
#         blcx = int(coord.x) - int(stampsize/2 + 0.5)
#         blcy = int(coord.y) - int(stampsize/2 + 0.5)
        trcx = blcx + stampsize
        trcy = blcy + stampsize

# Currently the way that positions on the edge of the skymap are handled is simple.
# They are blanked out, this means in a simple mean stacking 
# they will decrease the flux in the stacked image.
# This could certainly be improved upon.

        if (interval[blcx, trcx] in interval[0, imagesizes[coord.image][0]-1] 
                and interval[blcy, trcy] in interval[0, imagesizes[coord.image][1]-1]):

            data[i,:,:,0,0] = skymap[coord.image][blcx:trcx,blcy:trcy,0,0]

        else:

            data[i,:,:,0,0] = 0

        if psfmode == 'star':

            for x,y in [(1,0), (0,1), (-1,0), (0,-1)]:
                blcx_sub = blcx+x*int(stampsize/2 + 0.5)
                blcy_sub = blcy+y*int(stampsize/2 + 0.5)
                trcx_sub = trcx+x*int(stampsize/2 + 0.5)
                trcy_sub = trcy+y*int(stampsize/2 + 0.5)

                if (interval[blcx_sub, trcx_sub] 
                        in interval[0, imagesizes[coord.image][0]-1] 
                        and interval[blcy_sub, trcy_sub] 
                        in interval[0, imagesizes[coord.image][1]-1]):

                    data[i,:,:,0,0] -= 0.25*skymap[coord.image][blcx_sub:trcx_sub,blcy_sub:trcy_sub,0,0]

                else:
#                     data[i,:,:,0,0] = 0
                    pass



def _stack_stack(method, coords):
    import numpy as np
    """
        Performs the actual stacking on the data in the stack. 
        All data should be loaded in to stack before calling this function.
    """
    pixels = np.zeros(data.shape[1:])

    if method == 'median':
        pixels = np.median(data[0:len(coords),:,:,:,:], 0)
    elif method == 'mean':
        pixels = np.average(data[0:len(coords),:,:,:,:], 0, [coord.weight for coord in coords])

    return pixels



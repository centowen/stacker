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
import csv 
import math
import random
import time
import os 
import shutil
import numpy as np
from taskinit import *

skymap = []
data = []
oldimagenames = []
oldoutimagesize = 0


def readCoords1Image(imagename, coordfile):
#	This part of the code reads in a coordinate list from file and save it as a list of pixel coordinates in the image.
#	Also here a buffer could be used to avoid redoing this unless a new coordinate list is given.
#	Doing this does not save alot of performance and adds to much complexity to the code.

	ia.open(imagename)
	cs = ia.coordsys()
	imshape = ia.shape()
	ia.done()
	coordreader = csv.reader(open(coordfile, 'rb'), delimiter=',');
	coords = []
	for row in coordreader:
		p = cs.convert(coordin=[float(row[0]),float(row[1]),0,0],absin=[True,True,True,True],unitsin=["deg", "deg", "pix", "pix"],absout=[True,True,True,True],unitsout=["pix", "pix", "pix", "pix"])
		if p[0] >= 0 and p[1] >= 0 and p[0] < imshape[0] and p[1] < imshape[1]:
			coords.append((p[0], p[1]))

# 		if(p[0] > outimagesize and p[1] > outimagesize and p[0] < imagesize - outimagesize and p[1] < imagesize - outimagesize) : 
	return coords


def readCoords(imagenames, coordfile):
	coords = []
	if isinstance(imagenames, str):
		imagenames = [imagenames]

	for (i, imagename) in enumerate(imagenames):
		for coord in readCoords1Image(imagename, coordfile):
			coords.append((i, coord[0], coord[1]))
	return coords


def stack(imagenames, coords, outfile, outimagesize, stacking_method = 'median', weighting='none'):
	global skymap
	global data
	global oldimagenames
	global oldoutimagesize
	
	if isinstance(imagenames, str):
		imagenames = [imagenames]

#	To improve performance this module will keep buffers between run.
#	This following code resets these buffers if they have grown obsolete.
	if oldimagenames == []:
		oldimagenames = imagenames

	if oldimagenames != imagenames:
		oldimagenames = imagenames
		skymap = []
		data = []

	if oldoutimagesize == 0:
		oldoutimagesize = outimagesize
	elif oldoutimagesize != outimagesize:
		oldoutimagesize = outimagesize
		data = []
		skymap = []
	
# 	if not(data == []):
# 		print len(coords), data.shape[0], 'lol'
	if not(data == []) and len(coords) != data.shape[0]:
		data = []
		print 'peti reset'

#	This is infact not neccessary in many case, just leave it like this for simplicity.
	ia.open(imagenames[0])
	cs = ia.coordsys()
	outnchans = ia.boundingbox()['trc'][2]+1
	outnstokes = ia.boundingbox()['trc'][3]+1
	ia.done()

	imagesizes = []
	for imagename in imagenames:
		ia.open(imagename)
		imagesizes.append(ia.shape()[0])
		if ia.shape()[2] != outnchans or ia.shape()[3] != outnstokes:
			print('Channels do not match in all images!')
			return
		ia.done()


#	If there is no data buffer create one.
#	The data buffer is used to save the right stacking positions before stacking them.
#	During stacking this is where the full stack will actually be saved.
	if data == []:
		print 'allocating',len(coords)
		data = np.zeros((len(coords), outimagesize, outimagesize, outnstokes, outnchans))
	else:
		data = 0.*data
	
#	If there is no skymap buffer create one.
#	This is the data that is most important to buffer.
#	Reading a skymap from disk is a time consuming task and we don't want to do this to much.
	if skymap == []:
# 		skymap = imval(imagename, box = '{0},{1},{2},{3}'.format(0,0,ia.boundingbox()['trc'][0], ia.boundingbox()['trc'][1]))
		for imagename in imagenames:
			ia.open(imagename)
			skymap.append(ia.getregion())
			ia.done()

	sigma = []

#	Create the stack, also calculates a sigma for each image in the stack.
	for (i,coord) in enumerate(coords):
		blcx = int(coord[1]) - int(outimagesize/2 + 0.5)
		blcy = int(coord[2]) - int(outimagesize/2 + 0.5)
		trcx = blcx + outimagesize
		trcy = blcy + outimagesize

#		Currently the way that positions on the edge of the skymap are handled is simple.
#		They are blanked out, this means in a simple mean stacking they will decrease the flux in the stacked image.
#		This could certainly be improved upon.
		if not (blcx < 0 or blcy < 0 or trcx > imagesizes[coord[0]] or trcy > imagesizes[coord[0]]):
			data[i,:,:,0,0] = skymap[coord[0]][blcx:trcx,blcy:trcy,0,0]
			sigma.append(np.std(data[i,:,:,0,0].reshape(-1)))
		else:
#			It could be useful in stacking to know that this position is not useful.
			sigma.append(float(99999999))

	sigma = np.array(sigma)
        f = open('noise.list', 'w')
        for i in range(len(coords)):
            f.write('{2}\n'.format(coords[i][0], coords[i][1], 1/(sigma[i])**2))
        f.close()


	pixels = np.zeros((outimagesize, outimagesize, outnstokes, outnchans))
	if stacking_method == 'median':
		pixels = np.median(data, 0)
	elif stacking_method == 'mean':
		if weighting == 'sigma':
			norm_sigma = 0.
			for i in range(data.shape[0]):
				norm_sigma = norm_sigma + 1./float(sigma[i]/float(np.mean(sigma)))
				data[i,:,:,:,:] = 1./float(sigma[i]/float(np.mean(sigma)))*data[i,:,:,:,:]
			pixels = np.sum(data, 0)*(1/float(norm_sigma))
		else:
			pixels = np.mean(data, 0)

		
	if os.access(outfile, os.F_OK): shutil.rmtree(outfile)
	ia.open(imagenames[0])
	beam = ia.restoringbeam()
	ia.done()
	csnew = cs.copy()
	csnew.setreferencevalue([0.]*2, 'dir')
	csnew.setreferencepixel([int(outimagesize/2+0.5)]*2, 'dir')
        print csnew.torecord()
	ia.fromarray(outfile, pixels=pixels, csys = csnew.torecord())
	ia.open(outfile)
 	ia.setrestoringbeam(beam=beam)
	ia.done()


def calculateNoise(imagenames, coords, nrandom = 50, stacking_method = 'median', weighting = 'none'):
	random.seed(time.time())
	outimagesize = 3

	ia.open(imagenames[0])
	cs = ia.coordsys()
 	beam = qa.convert(ia.restoringbeam()['major'], 'deg')['value']
	beam = cs.convert(coordin=[float(beam),float(beam),0,0],absin=[False,False,False,False],unitsin=["deg", "deg", "pix", "pix"],absout=[False,False,False,False],unitsout=["pix", "pix", "pix", "pix"])
	beam = abs(beam[0])

	px = []
	for ifg in range(nrandom):
		randomcoords = []
		for coord in coords:
			dr = random.uniform(beam, 5*beam);
			dphi = random.uniform(0,2*math.pi);
			x = coord[1]+dr*math.cos(dphi)
			y = coord[2]+dr*math.sin(dphi)
			randomcoords.append((coord[0],x,y))
		
		stack(imagenames, randomcoords, 'stacked_noise.image', outimagesize, stacking_method, weighting)
		ia.open('stacked_noise.image')
		px.append(getFlux('stacked_noise.image'))
		ia.done()

	px = np.array(px).reshape(-1)
	ia.done()
	return np.std(px)


def getFlux(imagename):
    ia.open(imagename)
    cs = ia.coordsys()
    x = int(cs.referencepixel()['numeric'][0])
    y = int(cs.referencepixel()['numeric'][1])
    return float(ia.getregion(region=rg.box([x+1,y+1], [x+1,y+1])))



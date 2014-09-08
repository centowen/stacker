import stacker
import stacker.image
import stacker.uv
import stacker.modsub
import numpy as np

# Do Monte Carlo noise estimate, can be time consuming
donoise = False


if os.access('output', os.F_OK): shutil.rmtree('output')
os.mkdir('output')

# --- Image the data set to locate bright sources. ---
clean('testdata.ms', 'output/full', 
      imagermode='csclean',
      phasecenter = 'J2000 3h49m10.987 -30d00m00.00',
      cell = '.25arcsec', imsize=1600,
      pbcor = True, niter=1000, minpb=0.05,
      threshold = '0.1mJy', #Only clean bright sources, not target sources.
      psfmode='hogbom')

# --- Create a residual data set where bright sources are removed. ---

# First produce component list of model
stacker.modsub.cl_from_im('output/full.model', 'output/full.cl',
                  threshold = 1e-6)
# and subtract the component list from the data.
stacker.modsub.modsub('output/full.cl', 
        'testdata.ms', 'output/residual.ms',
        primarybeam='constant') # No pb since .model is not pbcorrected.
# Finally image the residual for image stacking and local noise estimation.
clean('output/residual.ms', 'output/residual',
      imagermode='csclean',
      phasecenter = 'J2000 3h49m10.987 -30d00m00.00',
      cell = '.25arcsec', imsize=1600,
      pbcor = True, niter=0, minpb=0.05)

# --- Do the stacking ---

stampsize = 64 # Size of stacked stamp

# Create a coordinate descriptor
coords = stacker.readCoords('coordinates.list')

# Calculate position specific weigths from noise in residual image
coords = stacker.image.calculate_sigma2_weights(coords,
        imagenames = ['output/residual.image'],
        stampsize = stampsize,
        maskradius = 5) # Excludes pixels closer to centre than 5 pixels.

# Actual stack, in both image and uv domain.
flux = {}

flux['image'] = stacker.image.stack(coords, 'output/imstacked.image', 
                    imagenames=['output/residual.image'], 
                    stampsize=stampsize, 
                    method='mean', # As opposed to median.
                    weighting=None) # This ensures that the weights set 
                                    # are used. weighting='sigma2' would
                                    # would recalculate weights.

flux['uv'] = stacker.uv.stack(coords, 'output/residual.ms', 'output/uvstacked.ms')

# Calculate the resulting noises
noise = {}
simplenoise = {}

ms.open('output/uvstacked.ms')
corrected_data = ms.getdata(['corrected_data'])['corrected_data']
ms.done()

simplenoise['image'] = np.sum([c.weight for c in coords])**-.5
simplenoise['uv'] = np.std(np.real(corrected_data))/np.sqrt(np.prod(corrected_data.shape))

# image flux
if donoise:
    noise['image'] = stacker.image.noise(coords, imagenames=['output/residual.image'], stampsize = stampsize, maskradius=5)
else:
    noise['image'] = simplenoise['image']

# uv flux
if donoise:
    noise['uv'] = stacker.uv.noise(coords, 'output/residual.ms', 
            ['output/residual.image'], stampsize = stampsize, maskradius=5)
else:
    noise['uv'] = simplenoise['uv']

print('Stacking results:')
if not donoise:
    print('warning, noise estimate based on image noise')
    print('may not be accurate (simulation indicate typically')
    print('within 20\%)')
print('')
print('image-stacking flux: {0:.1f}+-{1:.1f}'.format(flux['image']*1e6, noise['image']*1e6))
print('uv-stacking flux: {0:.1f}+-{1:.1f}'.format(np.real(flux['uv'])*1e6, noise['uv']*1e6))

# Image the uv-stacked data to produce an image.
clean('output/uvstacked.ms', 'output/uvstacked',
      cell = '.25arcsec', imsize=stampsize,
      mask = [int(stampsize/2)-2, int(stampsize/2)-2,
              int(stampsize/2)+2, int(stampsize/2)+2])


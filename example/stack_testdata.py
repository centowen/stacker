import stacker
import stacker.image
import stacker.uv
import stacker.modsub
import numpy as np

# --- Constants ---
donoise = True # Do Monte Carlo noise estimate, can be time consuming.


# --- First define some useful functions. ---

# Estimate the flux and noise for the source in the phase centre
# of the measurement vis
def estimate_flux_and_noise(vis):
    ms.open(vis)
    # This selection is not needed here, but required for mosaiced data sets
    # (and it does no harm here).
    ms.select({'field': 0}) 
    corrected_data = ms.getdata(['corrected_data'])['corrected_data']
    weight = ms.getdata(['weight'])['weight']
    flag = ms.getdata(['flag'])['flag']
    ms.done()

    # Easier to work with a mask in place of a set of flags.
    mask = np.logical_not(flag)

    # Set all flagged visibilites to zero before averaging.
    # Then sum over all the channels (only 1 weight per channel)>
    masked_data = np.real(np.sum(corrected_data*mask, axis=1))

    # We need to normalize with the mask, however, we want to this in combination
    # with the weight normalization.
    normalized_mask = np.sum(mask, axis=1)

    flux = np.sum(masked_data*weight)/np.sum(normalized_mask*weight)
    # Noise can approximated by the "weighted" standard deviation.
    noise = np.sqrt(np.sum((masked_data-flux*normalized_mask)**2*weight))
    noise = noise / np.sum(normalized_mask*weight)

    return flux, noise


def gaussian_residual(p, data, weight, flag):
    weighted_data = np.sum((data-p[0])*np.logical_not(flag), axis=1)*weight

def model_fit(vis):
    ms.open(vis)
    corrected_data = ms.getdata(['corrected_data'])['corrected_data']
    weight = ms.getdata(['weight'])['weight']
    flag = ms.getdata(['flag'])['flag']
    ms.done()

    mask = np.logical_not(flag)


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
# Note the use of primarybeam='constant', as the model is NOT primary-beam
# corrected this will result in the correct subtraction
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

print("Starting to stack.")
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
                    method='mean',  # As opposed to median.
                    weighting=None) # This ensures that the weights set 
                                    # are used. weighting='sigma2' would
                                    # would recalculate weights.

flux['uv'] = stacker.uv.stack(coords, 'output/residual.ms', 'output/uvstacked.ms')

# Calculate the resulting noises
noise = {}
simplenoise = {}

simplenoise['image'] = np.sum([c.weight for c in coords])**-.5
flux['uv'], simplenoise['uv'] = estimate_flux_and_noise('output/uvstacked.ms')

# Do a uv-model fitting for comparison
# model = 'output/modelfit.cl'
# uvmodelfit('output/uvstacked.ms', comptype='G', 
#            sourcepar=[1., 0., 0., 1., 1., 0.], outfile=model)
# cl.open(model)
# flux['uv_model'] = cl.getcomponent(0)['flux']['value'][0]
# size = qa.convert(cl.getcomponent(0)['shape']['majoraxis'], 'arcsec')['value']
# cl.done()

# image flux
if donoise:
    noise['image'] = stacker.image.noise(coords, imagenames=['output/residual.image'], stampsize = stampsize, maskradius=5)
else:
    noise['image'] = simplenoise['image']

# uv flux
if donoise:
    noise['uv'] = stacker.uv.noise(coords, 'output/residual.ms', 
            imagenames = ['output/residual.image'], stampsize = stampsize, maskradius=5)
else:
    noise['uv'] = simplenoise['uv']

print('Stacking results:')
if not donoise:
    print('warning, noise estimate based on image noise')
    print('may not be accurate (simulation indicate typically')
    print('within 20\%)')

print('')
print('image-stacking flux: {0:.1f}+-{1:.1f} uJy'.format(flux['image']*1e6,
                                                         noise['image']*1e6))
print('uv-stacking flux: {0:.1f}+-{1:.1f} uJy'.format(np.real(flux['uv'])*1e6,
                                                      noise['uv']*1e6))
# print('uv-model fitting flux: {0:.1f} uJy, size: {0:.1f}arcsec'.format(
#     np.real(flux['uv_model'])*1e6, size))

# Image the uv-stacked data to produce an image.
clean('output/uvstacked.ms', 'output/uvstacked',
      cell = '.25arcsec', imsize=stampsize,
      mask = [int(stampsize/2)-2, int(stampsize/2)-2,
              int(stampsize/2)+2, int(stampsize/2)+2])


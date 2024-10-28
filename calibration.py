from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np


# Display image given by the FITS file
# Vary vmin and vmax until image is visible 
# vmin and vmax are scaling factors, vmin > 1, vmax < 1
# (restricting the scale for the linear color map).
def display_img(file, vmin, vmax):
    img = file.data
    plt.figure(figsize=[6,6])

    fig = plt.imshow(img.data,vmin=vmin*np.min(img),vmax=vmax*np.max(img),cmap='plasma')
    plt.colorbar(fig,fraction=0.046,pad=0.04)

# Median combine an array of images as FITS data arrays
def median_combine(imgs):
    stack = np.dstack(tuple(imgs))
    master = np.median(stack,axis=2)

    return master

# input open fits files in primary hdu [0]
# filter as string 'ha' or 'oiii'
def calib_bias(bias_files):
    bias_imgs = [bias_img.data for bias_img in bias_files]
    master_bias = median_combine(bias_imgs)
    hdu = fits.PrimaryHDU(master_bias)
    hdu.writeto('master_bias.fit', overwrite=True)
    return fits.ImageHDU(master_bias, bias_files[0].header)

# Calibrate dark images
# Subtract the master bias from each raw dark image
# then median combine into a master dark. 
def calib_darks(dark_files, master_bias_fname):
    raw_darks = [file.data for file in dark_files]

    with fits.open(master_bias_fname) as mb:
        master_bias = mb[0].data

    bias_subtracted_darks = [dark - master_bias for dark in raw_darks]
    master_dark = median_combine(bias_subtracted_darks)

    hdu = fits.PrimaryHDU(master_dark)
    
    hdu.header = dark_files[0].header
    hdu.writeto('master_dark.fit', overwrite=True)

    return fits.ImageHDU(master_dark, dark_files[0].header)

# Calibrate flat images
# Subtract master bias from each flat image
# Subtract master dark scaled by exposure time
# Median combine and normalize by dividing by median
def calib_flats(flat_files, master_bias_fname, master_dark_fname, filter):
    #bias_subtracted_flats = [flat - master_bias for flat in raw_flats]
    #dark_subtracted_flats = [flat - master_dark*flat.header["EXPTIME"]/dark.header["EXPTIME"]]
    with fits.open(master_bias_fname) as mb:
        master_bias = mb[0].data

    with fits.open(master_dark_fname) as md:
        master_dark = md[0].data
        dark_exptime = md[0].header['EXPTIME']

    calibrated_flats = []

    for raw_flat in flat_files:
        bias_subtracted_flat = raw_flat.data - master_bias
        dark_subtracted_flat = bias_subtracted_flat - master_dark*raw_flat.header["EXPTIME"]/dark_exptime
        calibrated_flats.append(dark_subtracted_flat)

    median_combined_flats = median_combine(calibrated_flats)

    # Normalized master flat
    master_flat = median_combined_flats/np.median(median_combined_flats)

    hdu = fits.PrimaryHDU(master_flat)
    hdu.header = flat_files[0].header
    hdu.writeto('master_flat_' + filter + '.fits', overwrite=True)
    return fits.ImageHDU(master_flat, hdu.header)

# Calibrate science images
def calib_lights(lightHDUs, master_bias, master_dark, master_flat):
    calibratedFiles = []
    for img in lightHDUs:
        calibratedFiles.append(fits.ImageHDU((img.data-master_bias.data-img.header["EXPTIME"]/master_dark.header["EXPTIME"]*master_dark.data)/master_flat.data,
                               img.header))
    return calibratedFiles
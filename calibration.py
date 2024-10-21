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

# Median combine an array of images as FITS files
def median_combine(files):
    imgs = [file.data for file in files]

    stack = np.dstack(tuple(imgs))
    master = np.median(stack,axis=2)

    return master

def calib_bias(bias_files):
    master_bias = median_combine(bias_files)
    hdu = fits.PrimaryHDU(master_bias)
    hdu.writeto('master_bias.fits', overwrite=True)

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
    hdu.writeto('master_dark.fits', overwrite=True)

    md = fits.open("master_dark.fit")[0]
    with fits.open("master_dark.fits") as hdul:
        hdr = hdul[0].header
        hdr["EXPTIME"] = dark_files.header["EXPTIME"]

# Calibrate flat images
# Subtract master bias from each flat image
# Subtract master dark scaled by exposure time
# Median combine and normalize by dividing by median
def calib_flats(flat_files, master_bias_fname, master_dark_fname):
    #bias_subtracted_flats = [flat - master_bias for flat in raw_flats]
    #dark_subtracted_flats = [flat - master_dark*flat.header["EXPTIME"]/dark.header["EXPTIME"]]
    with fits.open(master_bias_fname) as mb:
        master_bias = mb[0].data

    with fits.open(master_dark_fname) as md:
        master_bias = md[0].data

    for raw_flat in flat_files:
        bias_subtracted_flat = raw_flat.data - master_bias
        dark_subtracted_flat = bias_subtracted_flat - master_dark*raw_flat.header["EXPTIME"]/master_dark.header["EXPTIME"]
        # get master dark with FITS header
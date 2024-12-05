from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np

'''
Calibration routines
'''


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
    #hdu = fits.PrimaryHDU(master_bias)
    #hdu.writeto('master_bias.fit', overwrite=True)
    return fits.ImageHDU(master_bias, bias_files[0].header)

# Calibrate dark images
# Subtract the master bias from each raw dark image
# then median combine into a master dark. 
def calib_darks(dark_files, master_bias):
    raw_darks = [file.data for file in dark_files]

    bias_subtracted_darks = [dark - master_bias.data for dark in raw_darks]
    master_dark = median_combine(bias_subtracted_darks)

    return fits.ImageHDU(master_dark, dark_files[0].header)

# Calibrate flat images
# Subtract master bias from each flat image
# Subtract master dark scaled by exposure time
# Median combine and normalize by dividing by median
def calib_flats(flat_files, master_bias, master_dark):
    calibrated_flats = []

    for raw_flat in flat_files:
        bias_subtracted_flat = raw_flat.data - master_bias.data
        dark_subtracted_flat = bias_subtracted_flat - master_dark.data*raw_flat.header["EXPTIME"]/master_dark.header['EXPTIME']
        calibrated_flats.append(dark_subtracted_flat)

    median_combined_flats = median_combine(calibrated_flats)

    # Normalized master flat
    master_flat = median_combined_flats/np.median(median_combined_flats)

    return fits.ImageHDU(master_flat, flat_files[0].header)

# Calibrate science images
def calib_lights(lightHDUs, master_bias, master_dark, master_flat):
    calibratedFiles = []
    for img in lightHDUs:
        calibratedFiles.append(fits.ImageHDU((img.data-master_bias.data-img.header["EXPTIME"]/master_dark.header["EXPTIME"]*master_dark.data)/master_flat.data,
                               img.header))
    return calibratedFiles

def zeroPadLeft(size, index):
    index = str(index)
    while(len(index)<size):
        index = "0"+index
    return index

# Calibrate science images for a data set
def calib_driver(lightPath, biasPath, flatPath, biasSuffix, session):
    hLights = []
    oLights = []
    for i in range(1,11):
        s = zeroPadLeft(3, i)
        hPath = lightPath + s + "-ha.fit"
        oPath = lightPath + s + "-o.fit"
        hLights.append(fits.open(hPath)[0])
        oLights.append(fits.open(oPath)[0])

    biases = []
    darks = []
    hFlats = []
    oFlats = []
    for i in range(1, 8):
        s = zeroPadLeft(3,i)
        biPath = biasPath + s + biasSuffix +".fit"
        dPath = biasPath + s + "-d.fit"
        flatHPath = flatPath + s + "-ha.fit"
        flatOPath = flatPath + s + "-o.fit"
        biases.append(fits.open(biPath)[0])
        darks.append(fits.open(dPath)[0])
        hFlats.append(fits.open(flatHPath)[0])
        oFlats.append(fits.open(flatOPath)[0])

    #calibrations
    master_bias = calib_bias(biases)
    master_dark = calib_darks(darks, master_bias)
    master_flat_ha = calib_flats(hFlats,master_bias, master_dark)
    master_flat_o = calib_flats(oFlats, master_bias, master_dark)
    fits.writeto("master_bias_" + str(session) + ".fit", master_bias.data, master_bias.header, overwrite=True)
    fits.writeto("master_dark_" + str(session) + ".fit", master_dark.data, master_dark.header, overwrite=True)
    fits.writeto("master_flat_ha_" + str(session) + ".fit", master_flat_ha.data, master_flat_ha.header, overwrite=True)

    c_hLights = calib_lights(hLights, master_bias, master_dark, master_flat_ha)
    c_oLights = calib_lights(oLights, master_bias, master_dark, master_flat_o)

    return master_bias, master_dark, master_flat_ha, master_flat_o, c_hLights, c_oLights


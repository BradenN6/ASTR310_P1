# Imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
import aperE
from mpl_toolkits.mplot3d import Axes3D

# EGAIN in FITS Header
egain = 1.2999999523162842
Kccd = 1/egain

# Run Matthew's Aperture photometry routine
def aper_driver(lightHDU, coords, base_aperture, rad1_range, rad2_range, angle):

    apertures = []
    for i in np.linspace(-rad1_range, rad1_range, 0.1):
        for j in np.linspace(-rad2_range, rad2_range, 0.1):
            apertures.append((base_aperture[0]+i, base_aperture[1]+j, base_aperture[2], \
                              base_aperture[3], base_aperture[4], base_aperture[5]))

    img = lightHDU.data

    vals = []
    for aperture in apertures:
        flux, err = aperE.photometry(img, coords[0], coords[1], aperture[0], aperture[1], \
                          aperture[2], aperture[3], aperture[4], aperture[5], egain, angle=angle, verbose=1)
        
        vals.append((flux, err))

    return coords, apertures, vals

def aperE_driver(lightHDU, coords, base_aperture, rad1_range, rad2_range, rad_step):
    '''
    Given lightHDU (FITS[0]), center coords of the object (x, y), and
    a list of test apertures. Uses the aperE function to output a list
    of fluxes and uncertainties using each aperture.

    apertures = [(rad1, rad2, ir1, ir2, or1, or2), ...]
    '''

    apertures = []
    for i in np.arange(-rad1_range, rad1_range+1, rad_step): 
        for j in np.arange(-rad2_range, rad2_range+1, rad_step):
            apertures.append((base_aperture[0]+i, base_aperture[1]+j, base_aperture[2], \
                              base_aperture[3], base_aperture[4], base_aperture[5]))

    img = lightHDU.data

    flux_arr = []
    err_arr = []
    for aperture in apertures:
        flux, err = aperE.aperE(img, coords[0], coords[1], aperture[0], aperture[1], \
                          aperture[2], aperture[3], aperture[4], aperture[5], Kccd)
        
        flux_arr.append(flux)
        err_arr.append(err)

    return coords, apertures, flux_arr, err_arr
        
# Plot flux or SNR as an image or surface
# as function of rad1 and rad2
# type=0 -> imshow, type=1 -> surface
def surface_plot(vals, base_aperture, rad1_range, rad2_range, rad_step, type=0):
    rad1_len = rad1_range*2+1
    rad2_len = rad2_range*2+1

    rad1_min = base_aperture[0]-rad1_range
    rad1_max = base_aperture[0]+rad1_range

    rad2_min = base_aperture[1]-rad2_range
    rad2_max = base_aperture[1]+rad2_range

    extent = [rad1_min, rad1_max, rad2_min, rad2_max]

    y = np.linspace(rad1_min, rad1_max, rad1_len)
    x= np.linspace(rad2_min, rad2_max, rad2_len)

    # cols correspond to rad2, rows to rad1
    im = np.reshape(vals, (rad1_len, rad2_len))

    if type == 0:
        plt.imshow(im, extent=extent, vmin=min(vals), vmax=max(vals))
        plt.xlabel('Rad2 [Pix]')
        plt.ylabel('Rad1 [Pix]')
        plt.colorbar()
    else: 
        X, Y = np.meshgrid(x, y)

        # Create figure and axes
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        #ax = Axes3D(fig)

        # Plot the surface with extent
        surf = ax.plot_surface(X, Y, im, cmap='viridis')

        # Set labels and title
        ax.set_xlabel('Rad2 [Pix]')
        ax.set_ylabel('Rad1 [Pix]')
        ax.set_zlabel('Val') #TODO
        plt.title('Surface Plot')

        # Add a colorbar
        fig.colorbar(surf)

        #fig = plt.figure()
        #ax = Axes3D(fig)
        #surf = ax.plot_trisurf(x, y, z, cmap=cm.jet, linewidth=0.1)
        #fig.colorbar(surf)#, shrink=0.5, aspect=5)


# TODO
# Aperture photometry on nebulae
# testing apertures
# 3d-plot SNR to determine aperture to use
    
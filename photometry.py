# Imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
import aperE

# EGAIN in FITS Header
egain = 1.2999999523162842
Kccd = 1/egain

# Run Matthew's Aperture photometry routine
def aper_driver():
    return None

def aperE_driver(lightHDU, coords, apertures):
    '''
    Given lightHDU (FITS[0]), center coords of the object (x, y), and
    a list of test apertures. Uses the aperE function to output a list
    of fluxes and uncertainties using each aperture.

    apertures = [(rad1, rad2, ir1, ir2, or1, or2), ...]
    '''

    img = lightHDU.data

    vals = []
    for aperture in apertures:
        flux, err = aperE(img, coords[0], coords[1], aperture[0], aperture[1], \
                          aperture[2], aperture[3], aperture[4], aperture[5], Kccd)
        
        vals.append((flux, err))

    return coords, apertures, vals
        


    
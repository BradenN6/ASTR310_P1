from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import align
from align import dispFITS
import calibration as calib

def zeroPadLeft(size, index):
    index = str(index)
    while(len(index)<size):
        index = "0"+index
    return index



lightPath = "data\\20241008_07in_A310_NGC604\\science\\NGC 604-"
hLights = []
oLights = []
for i in range(1,11):
    s = zeroPadLeft(3, i)
    hPath = lightPath + s + "-ha.fit"
    oPath = lightPath + s + "-o.fit"
    hLights.append(fits.open(hPath)[0])
    oLights.append(fits.open(oPath)[0])
    
biasPath = "data\\20241008_07in_A310_NGC604\\calibration\\calib-"
flatPath = "data\\20241008_07in_A310_NGC604\\calibration\\flat-"
biases = []
darks = []
hFlats = []
oFlats = []
for i in range(1, 8):
    s = zeroPadLeft(3,i)
    biPath = biasPath + s + "-bi.fit"
    dPath = biasPath + s + "-d.fit"
    flatHPath = flatPath + s + "-ha.fit"
    flatOPath = flatPath + s + "-o.fit"
    biases.append(fits.open(biPath)[0])
    darks.append(fits.open(dPath)[0])
    hFlats.append(fits.open(flatHPath)[0])
    oFlats.append(fits.open(flatOPath)[0])


#calibrations
master_bias = calib.calib_bias(biases)
master_dark = calib.calib_darks(darks, "master_bias.fit")
master_flat_ha = calib.calib_flats(hFlats,"master_bias.fit", "master_dark.fit", "ha")
master_flat_o = calib.calib_flats(oFlats, "master_bias.fit", "master_dark.fit", "o")

#dispFITS(master_bias, 1, 1)
#dispFITS(master_dark, 1, 1)
#dispFITS(master_flat_ha, 2, 2)
#dispFITS(master_flat_o, 2, 2)

c_hLights = calib.calib_lights(hLights, master_bias, master_dark, master_flat_ha)
c_oLights = calib.calib_lights(oLights, master_bias, master_dark, master_flat_o)

starCenter = (1511, 822)

#dispFITS(c_hLights[0], 1, 2)
#dispFITS(c_oLights[0], 1, 2)

align.alignFrames(c_hLights, starCenter[0], starCenter[1], markup=True)

plt.show()
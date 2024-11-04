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


(master_bias, master_dark, master_flat_ha, master_flat_o, c_hLights, c_oLights) = calib.calib_driver("data\\20241008_07in_A310_NGC604\\science\\NGC 604-", "data\\20241008_07in_A310_NGC604\\calibration\\calib-","data\\20241008_07in_A310_NGC604\\calibration\\flat-")

starCenter = (1511, 822)

#dispFITS(c_hLights[0], 1, 2)
#dispFITS(c_oLights[0], 1, 2)

dispFITS(c_hLights[0], 1, 10)
dispFITS(fits.ImageHDU(align.imshift(c_hLights[0].data, 10, 10), c_hLights[0].header),1,10)

aligned = align.alignFrames(c_hLights, starCenter[0], starCenter[1],20,20)
added = aligned[0].data
for i in range(1, len(aligned),1):
    added += aligned[i].data
dispFITS(fits.ImageHDU(added, aligned[0].header), 1, 1, "Added")

plt.show()
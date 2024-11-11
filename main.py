from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import align
import dataReduce
from align import dispFITS
import calibration as calib

def zeroPadLeft(size, index):
    index = str(index)
    while(len(index)<size):
        index = "0"+index
    return index

SCH1 = (1501, 827)
SCO1 = (1507, 826)
SCH2 = (1240,908)
SCO2 = (1238, 909)
SCH3 = (1227, 890)
SCO3 = (1225, 890)

# (master_bias, master_dark, master_flat_ha, master_flat_o, c_hLights, c_oLights) = \
# calib.calib_driver("data\\20241008_07in_A310_NGC604\\science\\NGC 604-",\
#                     "data\\20241008_07in_A310_NGC604\\calibration\\calib-",\
#                     "data\\20241008_07in_A310_NGC604\\calibration\\flat-",\
#                     "-bi", 1)


# ### First data set
# aligned = align.alignFrames(c_hLights, SCH1[0], SCH1[1],20,20)
# addedFits = dataReduce.sum(aligned)
# dispFITS(addedFits, 1, 1, "Added Session 1 Ha")
# fits.writeto("DS-1_Ha.fit",addedFits.data, addedFits.header, overwrite=True)

# aligned = align.alignFrames(c_oLights, SCO1[0], SCO1[1], 20, 20)
# addedFits = dataReduce.sum(aligned)
# dispFITS(addedFits, 1, 1, "Added Session 1 O[III]")
# fits.writeto("DS-1_OIII.fit",addedFits.data, addedFits.header, overwrite=True)

# ### second data set
# (master_bias, master_dark, master_flat_ha, master_flat_o, c_hLights, c_oLights) = \
# calib.calib_driver("data\\20241017_07in_A310_NGC604\\science\\NGC604-",\
#                     "data\\20241017_07in_A310_NGC604\\callibration\\callib-",\
#                     "data\\20241017_07in_A310_NGC604\\callibration\\flats-",\
#                     "-b", 2)
# aligned = align.alignFrames(c_hLights, SCH2[0], SCH2[1])
# addedFits = dataReduce.sum(aligned)
# dispFITS(addedFits, 1, 1, "Added Session 2 Ha")
# fits.writeto("DS-2_Ha.fit",addedFits.data, addedFits.header, overwrite=True)

# aligned = align.alignFrames(c_oLights, SCO2[0], SCO2[1])
# addedFits = dataReduce.sum(aligned)
# dispFITS(addedFits, 1, 1, "Added Session 2 O[III]")
# fits.writeto("DS-2_OIII.fit",addedFits.data, addedFits.header, overwrite=True)

# ### third data set
# (master_bias, master_dark, master_flat_ha, master_flat_o, c_hLights, c_oLights) = \
# calib.calib_driver("data\\20241027_07in_A310_NGC604\\science\\NGC 604-",\
#                     "data\\20241027_07in_A310_NGC604\\calibration\\calib-",\
#                     "data\\20241027_07in_A310_NGC604\\calibration\\flats-",\
#                     "-bi", 3)
# aligned = align.alignFrames(c_hLights, SCH3[0], SCH3[1])
# addedFits = dataReduce.sum(aligned)
# dispFITS(addedFits, 1, 1, "Added Session 3 Ha")
# fits.writeto("DS-3_Ha.fit",addedFits.data, addedFits.header, overwrite=True)

# aligned = align.alignFrames(c_oLights, SCO3[0], SCO3[1])
# addedFits = dataReduce.sum(aligned)
# dispFITS(addedFits, 1, 1, "Added Session 3 O[III]")
# fits.writeto("DS-3_OIII.fit",addedFits.data, addedFits.header, overwrite=True)

plt.show()
Ha1 = fits.open("DS-1_Ha.fit")
Ha2 = fits.open("DS-2_Ha.fit")
Ha3 = fits.open("DS-3_Ha.fit")
O1 = fits.open("DS-1_OIII.fit")
O2 = fits.open("DS-2_OIII.fit")
O3 = fits.open("DS-3_OIII.fit")
dispFITS()
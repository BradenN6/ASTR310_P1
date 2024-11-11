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


Ha1 = fits.open("DS-1_Ha.fit")[0]
Ha2 = fits.open("DS-2_Ha.fit")[0]
Ha3 = fits.open("DS-3_Ha.fit")[0]
O1 = fits.open("DS-1_OIII.fit")[0]
O2 = fits.open("DS-2_OIII.fit")[0]
O3 = fits.open("DS-3_OIII.fit")[0]


# dispFITS(Ha1, 1, 10)
# dispFITS(Ha2, 1, 10)
# dispFITS(Ha3, 1, 10)

RS1 = (660, 547)
RS2 = (401, 622)
RS3 = (388, 602)

H1x2, H1y2 = align.findMaxPixelCoord(Ha1.data, RS1[0], RS1[1],10,10)
H2x2, H2y2 = align.findMaxPixelCoord(Ha2.data, RS2[0], RS2[1],10,10)
H3x2, H3y2 = align.findMaxPixelCoord(Ha3.data, RS3[0], RS3[1],10,10)

H1x, H1y = align.findMaxPixelCoord(Ha1.data, SCH1[0],SCH1[1],10,10)
H2x, H2y = align.findMaxPixelCoord(Ha2.data, SCH2[0],SCH2[1],10,10)
H3x, H3y = align.findMaxPixelCoord(Ha3.data, SCH3[0],SCH3[1],10,10)

shift2 = (H1x-H2x,H1y-H2y)
shift3 = (H1x-H3x,H1y-H3y)
angle1 = np.atan2(H1y2-H1y, H1x2-H1x)
angle2 = np.atan2(H2y2-H2y, H2x2-H2x)
angle3 = np.atan2(H3y2-H3y, H3x2-H3x)
rot2 = angle2-angle1
rot3 = angle3-angle1

shifted = [Ha1]
shifted.append(fits.ImageHDU(align.imshift(Ha2.data, shift2[1], shift2[0]),Ha2.header))
shifted.append(fits.ImageHDU(align.imshift(Ha3.data, shift3[1], shift3[0]),Ha3.header))
H = dataReduce.sum(shifted)
dispFITS(H,1,1)
plt.show()
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import align
import dataReduce
from align import dispFITS
import calibration as calib
import scipy.ndimage as ndimage
from aperE import photometry
import size
import csv

def zeroPadLeft(size, index):
    index = str(index)
    while(len(index)<size):
        index = "0"+index
    return index

SCH1 = (1501, 827)
SCO1 = (1507, 826)
SCH2 = (1240, 908)
SCO2 = (1238, 909)
SCH3 = (1227, 890)
SCO3 = (1225, 890)

SC = (SCH1, SCO1, SCH2, SCO2, SCH3, SCO3)

def drawRings(X, Y, radX, radY, irX, irY, orX, orY, angle=0):
    #target ellipse
    targetEllipse = patches.Ellipse([X,Y],2*radX, 2*radY,fill=False,color="green")
    targetEllipse.set_angle(angle*180/np.pi)
    plt.gca().add_patch(targetEllipse)
    #inner ellipse for the sky annulus
    innerEllipse = patches.Ellipse([X,Y],2*irX, 2*irY,fill=False,color="blue")
    innerEllipse.set_angle(angle*180/np.pi)
    plt.gca().add_patch(innerEllipse)
    #outter ellipse for the sky annulus
    outerEllipse = patches.Ellipse([X,Y],2*orX, 2*orY,fill=False,color="red")
    outerEllipse.set_angle(angle*180/np.pi)
    plt.gca().add_patch(outerEllipse)

'''
Image Calibration
'''
'''
###########################
# First Data Set 10/08/2024

# Calibrate science images
(master_bias, master_dark, master_flat_ha, master_flat_o, c_hLights, c_oLights) = \
calib.calib_driver("data\\20241008_07in_A310_NGC604\\science\\NGC 604-",\
                    "data\\20241008_07in_A310_NGC604\\calibration\\calib-",\
                    "data\\20241008_07in_A310_NGC604\\calibration\\flat-",\
                    "-bi", 1)

# Align Ha images in the data set, co-add, save
aligned = align.alignFrames(c_hLights, SCH1[0], SCH1[1],20,20)
addedFits = dataReduce.sum(aligned)
dispFITS(addedFits, 1, 1, "Added Session 1 Ha")
fits.writeto("DS-1_Ha.fit", addedFits.data, addedFits.header, overwrite=True)

# Align OIII images, co-add, save
aligned = align.alignFrames(c_oLights, SCO1[0], SCO1[1], 20, 20)
addedFits = dataReduce.sum(aligned)
dispFITS(addedFits, 1, 1, "Added Session 1 O[III]")
fits.writeto("DS-1_OIII.fit",addedFits.data, addedFits.header, overwrite=True)

############################
# Second Data Set 10/17/2024

# Calibrate Science images
(master_bias, master_dark, master_flat_ha, master_flat_o, c_hLights, c_oLights) = \
calib.calib_driver("data\\20241017_07in_A310_NGC604\\science\\NGC604-",\
                    "data\\20241017_07in_A310_NGC604\\callibration\\callib-",\
                    "data\\20241017_07in_A310_NGC604\\callibration\\flats-",\
                    "-b", 2)

# Align, co-add, save Ha and OIII
aligned = align.alignFrames(c_hLights, SCH2[0], SCH2[1])
addedFits = dataReduce.sum(aligned)
dispFITS(addedFits, 1, 1, "Added Session 2 Ha")
fits.writeto("DS-2_Ha.fit",addedFits.data, addedFits.header, overwrite=True)

aligned = align.alignFrames(c_oLights, SCO2[0], SCO2[1])
addedFits = dataReduce.sum(aligned)
dispFITS(addedFits, 1, 1, "Added Session 2 O[III]")
fits.writeto("DS-2_OIII.fit",addedFits.data, addedFits.header, overwrite=True)

###########################
# Third Data Set 10/27/2024

# Calibrate Science Images
(master_bias, master_dark, master_flat_ha, master_flat_o, c_hLights, c_oLights) = \
calib.calib_driver("data\\20241027_07in_A310_NGC604\\science\\NGC 604-",\
                    "data\\20241027_07in_A310_NGC604\\calibration\\calib-",\
                    "data\\20241027_07in_A310_NGC604\\calibration\\flats-",\
                    "-bi", 3)

# Align, co-add, save Ha and OIII
aligned = align.alignFrames(c_hLights, SCH3[0], SCH3[1])
addedFits = dataReduce.sum(aligned)
dispFITS(addedFits, 1, 1, "Added Session 3 Ha")
fits.writeto("DS-3_Ha.fit",addedFits.data, addedFits.header, overwrite=True)

aligned = align.alignFrames(c_oLights, SCO3[0], SCO3[1])
addedFits = dataReduce.sum(aligned)
dispFITS(addedFits, 1, 1, "Added Session 3 O[III]")
fits.writeto("DS-3_OIII.fit",addedFits.data, addedFits.header, overwrite=True)

###########################
# Align/Co-add All Datasets

Ha1 = fits.open("DS-1_Ha.fit")[0]
Ha2 = fits.open("DS-2_Ha.fit")[0]
Ha3 = fits.open("DS-3_Ha.fit")[0]
O1 = fits.open("DS-1_OIII.fit")[0]
O2 = fits.open("DS-2_OIII.fit")[0]
O3 = fits.open("DS-3_OIII.fit")[0]

dispFITS(Ha1, 1, 10)
dispFITS(Ha2, 1, 10)
dispFITS(Ha3, 1, 10)

# Align and co-add all three datasets
H, O = align.align_datasets(Ha1, Ha2, Ha3, O1, O2, O3, SC)

dispFITS(H,1,1, "H-alpha")
dispFITS(O,1,1, "O-III")

fits.writeto("H.fit", H.data, H.header, overwrite=True)
fits.writeto("O.fit", O.data, O.header, overwrite=True)
'''

'''
Brightness Analysis
'''

# Clockwise angle
coordsH = [[1408.,572.,110.,110.,140.,140.,0],
           [509.,198.,130.,120.,160.,170.,np.pi/4],
           [1508.,190.,50.,80.,100.,110.,-np.pi/6],
           [960.,888.,50.,50.,70.,70.,0],
           [786.,1255.,40.,40.,60.,60.,0],
           [1901.,142.,60.,90.,100.,120.,-np.pi/8]]
coordsO = [[1406,570.,52.,50.,70.,65.,0],
           [509.,198.,50.,50.,70.,70., np.pi/4],
           [1506.,201.,50.,50.,70.,80.,0],
           [959.,887.,40.,40.,70.,70.,0],
           [780.,1251.,45.,45.,60.,60.,0],
           [1896.,151.,30.,35.,60.,50.,np.pi/8]
           ]

H = fits.open("H.fit")[0]
O = fits.open("O.fit")[0]

dispFITS(H,1,1,"H")




dispFITS(O,0.2,0.3,"O")


'''
Size Analysis
'''

#rX, rY = size.findSize(O.data,coordsO[0][0],coordsO[0][1],52,50,70,65)
#photometry(O.data, coordsO[0][0],coordsO[0][1],rX,rY,52,50,70,65,1.3, 0,True)
dataO = []
dataH = []
for c in coordsO:
     rX, rY = size.findSize(O.data,c[0], c[1],c[2],c[3],c[4],c[5],c[6],5)
     e, s_e =photometry(O.data,c[0],c[1],rX,rY,c[1],c[2],c[3],c[4],c[5],c[6],True)
     dataO.append([c[0],c[1],rX,rY,e,s_e])
     print(dataO[-1])
for c in coordsH:
    rX, rY = size.findSize(H.data,c[0], c[1],c[2],c[3],c[4],c[5],c[6],5)
    photometry(H.data,c[0],c[1],rX,rY,c[1],c[2],c[3],c[4],c[5],c[6],True)
    dataH.append([c[0],c[1],rX,rY,e,s_e])
    print(dataH[-1])
with open("photometry.csv", "w", newline="") as csvfile:
    writer = csv.writer(csvfile, delimiter=",",quotechar="|",quoting=csv.QUOTE_MINIMAL)
    writer.writerow(["X", "Y", "rX", "rY", "e", "s_e"])
    for i in range(len(dataO)):
        writer.writerow([dataO[i][0], dataO[i][1],dataO[i][2],dataO[i][3],dataO[i][4],dataO[i][5]])
    for i in range(len(dataH)):
        writer.writerow([dataH[i][0], dataH[i][1],dataO[i][2],dataH[i][3],dataH[i][4],dataH[i][5]])
print(dataO)
print(dataH)
plt.show()
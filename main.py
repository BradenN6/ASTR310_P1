from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

def zeroPadLeft(size, index):
    index = str(i)
    while(len(index)<size):
        index = "0"+index
    return index

def dispFITS(hdu, medMinCoeff, medMaxCoeff):
    plt.figure(figsize=[6,6])
    fig = plt.imshow(hdu.data,vmin=medMinCoeff*np.median(hdu.data),vmax=medMaxCoeff*np.median(hdu.data),cmap='plasma')
    plt.colorbar(fig,fraction=0.046,pad=0.04)

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

dispFITS(hLights[0], 0.9, 1.1)
plt.show()
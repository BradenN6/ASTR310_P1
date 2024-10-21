from astropy.io import fits

def zeroPadLeft(size, index):
    index = str(i)
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
flatsH = []
flatsO = []
for i in range(1, 8):
    s = zeroPadLeft(3,i)
    biPath = biasPath + s + "-bi.fit"
    dPath = biasPath + s + "-d.fit"
    flatHPath = flatPath + s + "-ha.fit"
    flatOPath = flatPath + s + "-o.fit"
    biases.append(fits.open(biPath)[0])
    darks.append(fits.open(dPath)[0])
    flatsH.append(fits.open(flatHPath)[0])
    flatsO.append(fits.open(flatOPath)[0])


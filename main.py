from astropy.io import fits

def zeroPadLeft(size, index):
    index = str(i)
    while(len(index)<size):
        index = "0"+index
    return index

startPath = "data\\20241008_07in_A310_NGC604\\science\\NGC 604-"
hLights = []
oLights = []
for i in range(1,11):
    s = zeroPadLeft(3, i)
    hPath = startPath + s + "-ha.fit"
    oPath = startPath + s + "-o.fit"
    hLights.append(fits.open(hPath)[0])
    oLights.append(fits.open(oPath)[0])
    


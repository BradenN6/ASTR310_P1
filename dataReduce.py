from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from align import dispFITS
#function which takes in an array of aligned fits image HDUs and returns a fits image of their added values
def sum(imgs):
    newHeader = imgs[0].header
    expTime = newHeader["EXPTIME"]
    added = imgs[0].data
    for i in range(1, len(imgs),1):
        added += imgs[i].data
        expTime += imgs[i].header["EXPTIME"]
    newHeader["EXPTIME"] = expTime
    return fits.ImageHDU(added, newHeader)

def average(imgs):
    data = imgs[0].data
    num = np.ones(data.shape)
    # dispFITS(fits.ImageHDU(data),1,1,"sum1")
    # dispFITS(imgs[1],1,1,"i2")
    # dispFITS(imgs[3],1,1,"i3")
    # dispFITS(imgs[4],1,1,"i4")
    # dispFITS(imgs[5],1,1,"i5")
    # dispFITS(fits.ImageHDU(num),1,1,"num1")
    for i in range(1,len(imgs)):
        for y in range(imgs[i].data.shape[1]):
            for x in range(imgs[i].data.shape[0]):
                if imgs[i].data[x,y] > 0:
                    num[x,y] += 1
                    data[x,y] += imgs[i].data[x,y]
    #dispFITS(fits.ImageHDU(data),1,1,"sum")
    #dispFITS(fits.ImageHDU(num),1,1,"num")
    data /= num
    # dispFITS(fits.ImageHDU(data),1,1,"avg")
    # plt.show()
    return fits.ImageHDU(data)


def stack(imgs, mul, sHigh=3, sLow=3):
    avg = average(imgs)
    dispFITS(avg, 1,1, "Average")
    sum = np.zeros(imgs[0].data.shape)
    num = np.zeros(imgs[0].data.shape)
    for i in range(len(imgs)):
        for y in range(imgs[i].data.shape[1]):
            for x in range(imgs[i].data.shape[0]):
                if imgs[i].data[x,y] > 0:
                    sum[x,y] += (avg.data[x,y]-imgs[i].data[x,y])**2
                    num[x,y] += 1
    stdDev = np.sqrt(sum/num)
    dispFITS(fits.ImageHDU(stdDev),3,3, "StdDev0")

    sums = np.zeros(imgs[0].data.shape)
    included = np.zeros(imgs[0].data.shape)
    for i in range(len(imgs)):
        for y in range(imgs[i].data.shape[1]):
            for x in range(imgs[i].data.shape[0]):
                if imgs[i].data[x,y] > 0 and imgs[i].data[x,y]-avg.data[x,y]>-sLow*stdDev[x,y] and imgs[i].data[x,y]-avg.data[x,y]<sHigh*stdDev[x,y]:
                    sums[x,y] += imgs[i].data[x,y]
                    included[x,y] += 1
    sums *= np.max(included)/included
    dispFITS(fits.ImageHDU(included),1,3, "Included")
    return fits.ImageHDU(sums), fits.ImageHDU(mul.data * np.max(included)/included)

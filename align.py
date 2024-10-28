# from imshift.py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import main
from astropy.io import fits

def imshift(im,nr,nc):
    """Shifts an image by nr rows and nc columns
    (which can be either positive or negative)"""

    a,b=im.shape
    imr=np.zeros(im.shape)
    ir1 = max(0, -nr)
    ir2 = min(a, a-nr)
    it1 = max(0, -nc)
    it2 = min(b, b-nc)
    r1=max(0,nr)
    r2=min(a,nr+a)
    c1=max(0,nc)
    c2=min(b,nc+b)
    imr[r1:r2, c1:c2] = im[ir1:ir2, it1:it2]
    return imr

#finds the coordinates of a the greatest ADU value in a range of values in a fits image.
def findMaxPixelCoord(data, guessX, guessY, rX, rY):
    maxVal = float("-inf")
    maxCoord = (-1,-1)
    for x in range(int(guessX-rX), int(guessX+rX)):
        for y in range(int(guessY-rY), int(guessY + rY)):
            if data[y,x] > maxVal:
                maxVal = data[y,x]
                maxCoord = (x,y)
    return maxCoord[0], maxCoord[1]

#function which takes in an array of calibrated light frames and returns an array of aligned frames
def alignFrames(hduArray, guessX, guessY, rX=10, rY=10, markup = False):
    shifts = []
    x0, y0 = findMaxPixelCoord(hduArray[0].data, guessX, guessY, rX, rY)
    testX, testY = x0, y0
    for i in range(1, len(hduArray)):
        x, y = findMaxPixelCoord(hduArray[i].data, testX, testY, rX, rY)
        shifts.append([testX-x,testY-y])
        if markup:
            main.dispFITS(hduArray[i],1, 10)
            plt.gca().add_patch(patches.Rectangle((testX-rX, testY-rY),2*rX, 2*rY), color="red", fill=False)
            plt.gca().add_patch(patches.Rectangle((x, y), 1, 1), color="black", fc=(0,0,0,0.3))
        testX += testX-x
        testY += testY-y
    shifted = []
    for img in hduArray:
        shifted.append(fits.ImageHDU(imshift(img.data, shifts[i][0], shifts[i][1]),img.header))
    return shifted
            
# from imshift.py
import numpy as np

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

def alignFrames(hduArray, guessX, guessY, rX=10, rY=10, markup = False):
    shifts = []
    x0, y0 = findMaxPixelCoord(hduArray[0].data, guessX, guessY, rX, rY)
    for i
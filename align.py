# from imshift.py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from astropy.io import fits
import scipy.ndimage as ndimage
import dataReduce

'''
Code for aligning images
'''

figNum = 1
def dispFITS(hdu, medMinCoeff, medMaxCoeff, title=None):
    global figNum
    plt.figure(figsize=[6,6])
    fig = plt.imshow(hdu.data,vmin=np.median(hdu.data)-medMinCoeff*np.std(hdu.data),vmax=np.median(hdu.data)+medMaxCoeff*np.std(hdu.data),cmap='plasma',extent=(0,hdu.header["NAXIS1"],hdu.header["NAXIS2"],0))
    plt.colorbar(fig,fraction=0.046,pad=0.04)
    if title == None:
        plt.title("Figure " + str(figNum))
    else:
        plt.title(title)
    plt.xlabel("Pixel x coordinate")
    plt.ylabel("Pixel y coordinate")
    figNum += 1

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

def rotate_image(image_array, angle):
    """Rotates an image array by a specified angle (in degrees)."""
    return ndimage.rotate(image_array, angle, reshape=False)

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

def trackStar(hduArray, startX, startY, rX=10, rY = 10):
    x0, y0 = findMaxPixelCoord(hduArray[0].data, startX, startY, rX, rY)
    positions = [(x0, y0)]
    testX, testY = x0, y0
    for i in range(1, len(hduArray)):
        x, y = findMaxPixelCoord(hduArray[i].data, testX, testY, rX, rY)
        positions.append((x,y))
        testX = x
        testY = y
    return positions

#function which takes in an array of calibrated light frames and returns an array of aligned frames
def alignFrames(hduArray, guessX, guessY, rX=10, rY=10, markup = False):
    coords = trackStar(hduArray, guessX, guessY, rX, rY)
    shifted = []
    if markup:
        plt.gca().add_patch(patches.Rectangle((guessX, guessY), 1, 1, color="red", fc=(1,0,0,0.3)))
        plt.gca().add_patch(patches.Rectangle((guessX-rX, guessY-rY),2*rX, 2*rY, color="red", fill=False))
        plt.gca().add_patch(patches.Rectangle((coords[0][0], coords[0][0]), 1, 1, color="black", fc=(0,0,0,0.3)))
    for i in range(1, len(hduArray)):
        shifted.append(fits.ImageHDU(imshift(hduArray[i].data, coords[0][1]-coords[i][1], coords[0][0]-coords[i][0]),hduArray[i].header))
        if markup:
            dispFITS(hduArray[i],1, 10)
            plt.gca().add_patch(patches.Rectangle((coords[i-1][0], coords[i-1][1]), 1, 1, color="red", fc=(1,0,0,0.3)))
            plt.gca().add_patch(patches.Rectangle((coords[i-1][0]-rX, coords[i-1][1]-rY),2*rX, 2*rY, color="red", fill=False))
            plt.gca().add_patch(patches.Rectangle((coords[i][0], coords[i][1]), 1, 1, color="black", fc=(0,0,0,0.3)))
  
    return shifted
            
# Manually aligning and rotating the Ha and OIII images
# from each dataset into single combined Ha and OIII images
def align_datasets(Ha1, Ha2, Ha3, O1, O2, O3, SC):
    SCH1, SCO1, SCH2, SCO2, SCH3, SCO3 = SC

    #SCH1 = (1501, 827)
    #SCO1 = (1507, 826)
    #SCH2 = (1240, 908)
    #SCO2 = (1238, 909)
    #SCH3 = (1227, 890)
    #SCO3 = (1225, 890)

    RS1 = (660, 547)
    RS2 = (401, 622)
    RS3 = (388, 602)

    H1x2, H1y2 = findMaxPixelCoord(Ha1.data, RS1[0], RS1[1],10,10)
    H2x2, H2y2 = findMaxPixelCoord(Ha2.data, RS2[0], RS2[1],10,10)
    H3x2, H3y2 = findMaxPixelCoord(Ha3.data, RS3[0], RS3[1],10,10)

    H1x, H1y = findMaxPixelCoord(Ha1.data, SCH1[0],SCH1[1],10,10)
    H2x, H2y = findMaxPixelCoord(Ha2.data, SCH2[0],SCH2[1],10,10)
    H3x, H3y = findMaxPixelCoord(Ha3.data, SCH3[0],SCH3[1],10,10)

    shift2 = (H1x-H2x,H1y-H2y)
    shift3 = (H1x-H3x,H1y-H3y)
    angle1 = np.atan2(H1y2-H1y, H1x2-H1x)
    angle2 = np.atan2(H2y2-H2y, H2x2-H2x)
    angle3 = np.atan2(H3y2-H3y, H3x2-H3x)
    rot2 = angle2-angle1
    rot3 = angle3-angle1

    shifted = [Ha1]
    shifted.append(fits.ImageHDU(rotate_image(imshift(Ha2.data, shift2[1], shift2[0]),rot2*180/np.pi),Ha2.header))
    shifted.append(fits.ImageHDU(rotate_image(imshift(Ha3.data, shift3[1], shift3[0]),rot3*180/np.pi),Ha3.header))
    H = dataReduce.sum(shifted)
    shifted = [O1]
    shifted.append(fits.ImageHDU(rotate_image(imshift(O2.data, shift2[1], shift2[0]),rot2*180/np.pi),O2.header))
    shifted.append(fits.ImageHDU(rotate_image(imshift(O3.data, shift3[1], shift3[0]),rot3*180/np.pi),O3.header))

    O = dataReduce.sum(shifted)
    #dispFITS(H,1,1, "H-alpha")
    #dispFITS(O,1,1, "O-III")

    return H, O
import numpy as np
from aperE import photometry
from align import dispFITS

'''
Size determination code through edge-detection.
'''

#returns the average sky value, the number of pixels, the error in the individual sky value, and the error in the average sky value
#all values are in ADU
def annulus(data, X, Y, irX, irY, orX, orY, angle=0, subsampleRate=10):
    annulusR = np.max([orX, orY])
    annulusTestPts = np.mgrid[int(X-annulusR):int(X+annulusR)+1, int(Y-annulusR):int(Y+annulusR)+1]
    annulusEdge=[]
    for t in np.linspace(0.,2*np.pi,int(628*annulusR)): #two sets of edges for an annulus
        xi = np.floor(irX*np.cos(t)*np.cos(angle)-irY*np.sin(t)*np.sin(angle)+X)
        yi = np.floor(irY*np.sin(t)*np.cos(angle)+irX*np.cos(t)*np.sin(angle)+Y)
        xo = np.floor(orX*np.cos(t)*np.cos(angle)-orY*np.sin(t)*np.sin(angle)+X)
        yo = np.floor(orY*np.sin(t)*np.cos(angle)+orX*np.cos(t)*np.sin(angle)+Y)
        xi2 = np.floor((irX-0.05)*np.cos(t)*np.cos(angle)-(irY-0.05)*np.sin(t)*np.sin(angle)+X)
        yi2 = np.floor((irY-0.05)*np.sin(t)*np.cos(angle)+(irX-0.05)*np.cos(t)*np.sin(angle)+Y)
        xo2 = np.floor((orX+0.05)*np.cos(t)*np.cos(angle)-(orY+0.05)*np.sin(t)*np.sin(angle)+X)
        yo2 = np.floor((orY+0.05)*np.sin(t)*np.cos(angle)+(orX+0.05)*np.cos(t)*np.sin(angle)+Y)
        if [xi, yi] not in annulusEdge:
            annulusEdge.append([xi,yi])
        if [xi2, yi2] not in annulusEdge:
            annulusEdge.append([xi2,yi2])
        if [xo, yo] not in annulusEdge:
            annulusEdge.append([xo,yo])
        if [xo2, yo2] not in annulusEdge:
            annulusEdge.append([xo2,yo2])
    annulusInterior = []
    for x in range(annulusTestPts[0,0,0],annulusTestPts[0,-1,0]+1):
        for y in range(annulusTestPts[1,0,0],annulusTestPts[1,0,-1]+1):
            x_transformed1 = ((x-X)*np.cos(-angle)-(y-Y)*np.sin(-angle))/irX
            y_transformed1 = ((y-Y)*np.cos(-angle)+(x-X)*np.sin(-angle))/irY
            x_transformed2 = ((x-X)*np.cos(-angle)-(y-Y)*np.sin(-angle))/orX
            y_transformed2 = ((y-Y)*np.cos(-angle)+(x-X)*np.sin(-angle))/orY
            #pixels have to be outside the inner radius of the sky annulus, but inside the outer radius of the sky anulus
            if ((x_transformed1)**2 + (y_transformed1)**2) >= 1 \
                and ((x_transformed2)**2 + (y_transformed2)**2) <= 1 \
                and [x,y] not in annulusEdge:
                annulusInterior.append([x,y])
    
    annulusEdgeWeights = []
    for coords in annulusEdge:
        interior = 0
        for sx in range(subsampleRate):
            for sy in range(subsampleRate):
                x_transformed1 = ((coords[0]+sx/subsampleRate-X)*np.cos(-angle)-(coords[1]+sy/subsampleRate-Y)*np.sin(-angle))/irX
                y_transformed1 = ((coords[1]+sy/subsampleRate-Y)*np.cos(-angle)+(coords[0]+sx/subsampleRate-X)*np.sin(-angle))/irY
                x_transformed2 = ((coords[0]+sx/subsampleRate-X)*np.cos(-angle)-(coords[1]+sy/subsampleRate-Y)*np.sin(-angle))/orX
                y_transformed2 = ((coords[1]+sy/subsampleRate-Y)*np.cos(-angle)+(coords[0]+sx/subsampleRate-X)*np.sin(-angle))/orY
                if ((x_transformed1)**2 + (y_transformed1)**2) >= 1 and ((x_transformed2)**2 + (y_transformed2)**2) <= 1:
                    interior += 1
        annulusEdgeWeights.append(interior/subsampleRate**2)

    skySum = 0 #ADU
    skyPixels = 0 
    #Sum up all the whole pixels in the sky annulus and keep track of how many pixels you summed
    for coords in annulusInterior:
        skySum += data[coords[1],coords[0]]
        skyPixels += 1
    #Sum up all the signal in the edge pixels inside the target ellipse and keep track of how many total pixels you summed
    for i in range(len(annulusEdge)):
        skySum += data[int(annulusEdge[i][1]),int(annulusEdge[i][0])]*annulusEdgeWeights[i]
        skyPixels += annulusEdgeWeights[i]
    skyValue = skySum/skyPixels #The average ADU value of the sky annulus
    sumSquareDiff = 0 #variable for finding the square portion of the Root-Mean-square formula for standard deviation
    #find the square portion for the whole pixels
    for coords in annulusInterior:
        sumSquareDiff += (data[coords[1],coords[0]]-skyValue)**2
    #Find the square portion for the partial pixels
    for i in range(len(annulusEdge)):
        sumSquareDiff += (data[int(annulusEdge[i][1]),int(annulusEdge[i][0])]-skyValue)**2 * annulusEdgeWeights[i]
    s_sky = np.sqrt(sumSquareDiff/skyPixels) #The standard deviation = uncertainty of each pixel of the sky background in ADU
    s_skyValue = s_sky/np.sqrt(skyPixels) #The uncertainty of the average sky background in ADU

    return skyValue, skyPixels, s_sky, s_skyValue


def findSize(data, X, Y, sirx, siry, sorx, sory, angle=0, sigma=5, subsamplerate=10):
    sky, skyPixels, s_skyP, s_skyValue = annulus(data, X, Y, sirx, siry, sorx, sory, angle, subsamplerate)
    s_sky = np.sqrt(s_skyP**2 + s_skyValue**2)

    innerX=7
    innerY=7
    switched = False
    axis = 0 #0=X, 1=Y

    print(sky,s_sky)
    while True:
        stepVal, stepPix, stepP, s_stepVal = annulus(data, X, Y,innerX, innerY,innerX+(axis==0), innerY+(axis==1), angle, subsamplerate)
        if stepVal*stepPix-sky*stepPix > sigma * np.sqrt(stepPix)*s_sky:
            if axis == 0:
                innerX += 1
            else:
                innerY += 1
        else:
            if switched:
                break
            switched = True
            axis = 1-axis
            continue
        switched = False

    return innerX, innerY
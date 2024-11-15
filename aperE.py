# aperE.py

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

def aperE(im, col, row, rad1, rad2, ir1, ir2, or1, or2, Kccd, saturation=np.inf):
    """Original code by Professor Alberto Bolatto, edited by Alyssa Pagan, and
    translated to Python by ChongChong He, further edited by Orion Guiffreda.

    Before using aperE.m, rotate your image using imrotate(im,angle) so the
    major axis of your object is perpendicular or parallel to your x or y axis.
    
    APER(im,col,row,rad1,rad1,ir1,ir2,or1,or2,Kccd) Do aperture photometry of image "im"
    for a star, galaxy or nebula centered at the "row,col" coordinates, For an ellipse 
    with a major and minor axis of "rad1,rad2" and an inner sky ellipse with a 
    major and minor axis of (ir1,ir2)and outer sky ellipse of "or1,or2" with CCD
    gain of Kccd ADU/electron. Optionally, a 11th parameter can be passed
    with the saturation value for the CCD.
    """

    a, b = im.shape
    xx, yy = np.meshgrid(range(b), range(a))
    ixsrc = ((xx - col) / rad1) ** 2 + ((yy - row) / rad2) ** 2 <= 1
    ixsky = np.logical_and(
        (((xx - col) / or1) ** 2) + (((yy - row) / or2) ** 2) <= 1,
        (((xx - col) / ir1) ** 2) + (((yy - row) / ir2) ** 2) >= 1
    )
    length = max(ixsky.shape)
    sky = np.median(im[ixsky], axis=0)
    imixsrc = im[ixsrc]
    pix = imixsrc - sky
    sig = np.sqrt(imixsrc / Kccd)
    ssig = np.std(im[ixsky]) / np.sqrt(length) / Kccd
    flx = np.sum(pix) / Kccd
    err = np.sqrt(np.sum(sig) ** 2 + ssig ** 2)
    issat = 0
    if max(imixsrc) > saturation:
        issat = 1
    fw = np.copy(or1)
    ix = np.where(
        np.logical_and(
            np.logical_and(
                np.logical_and(xx >= col - 2 * fw, xx <= col + 2 * fw),
                yy >= row - 2 * fw
            ),
            yy <= row + 2 * fw
        )
    )
    aa = np.sum(np.logical_and(xx[0, :] >= col - 2 * fw,
                               xx[0, :] <= col + 2 * fw))
    bb = np.sum(np.logical_and(yy[:, 0] >= row - 2 * fw,
                               yy[:, 0] <= row + 2 * fw))
    px = np.reshape(xx[ix], (bb, aa))
    py = np.reshape(yy[ix], (bb, aa))
    pz = np.reshape(im[ix], (bb, aa))
    plt.figure()
    plt.imshow(pz, extent=[px[0, 0], px[0, -1], py[0, 0], py[-1, 0]])
    plt.tight_layout()
    # if not np.isempty(imixsrc):
    #     np.caxis(np.concatenate((sky, np.array([max(imixsrc)]))))

    p = np.arange(360) * np.pi / 180
    xc = np.cos(p)
    yc = np.sin(p)
    plt.plot(col+rad1*xc, row+rad2*yc, 'w')
    plt.plot(col+ir1*xc, row+ir2*yc, 'r')
    plt.plot(col+or1*xc, row+or2*yc, 'r')
    if issat:
        plt.text(col, row, 'CHECK SATURATION', ha='center', color='w',
                 va='top', fontweight='bold')
        print('At the peak this source has {:0.0f} counts.'.format(
            max(imixsrc)))
        print('Judging by the number of counts, if this is a single exposure the')
        print('source is likely to be saturated. If this is the coadding of many')
        print('short exposures, check in one of them to see if this message appears.')
        print('If it does, you need to flag the source as bad in this output file.')
    plt.tight_layout()
    plt.savefig("aperE_img.pdf")
    return flx, err




#Matthew's code:

#A function which takes in a 2-dimentional data array of ADU values, data; the column of the center of the photometry annulus;
#the row of the center of the photometry annulus; the x-radius of the target ellipse in pixels; the y-radius of the target ellipse in pixels;
#the x-radius of the inner annulus in pixels; the y-radius of the inner annulus in pixels; the x-radius of the outer annulus in pixels;
#the y-radius of the outer annulus in pixels; the gain of the CCD in e-/ADU;
# and the clockwise angle in radians that the photometry annulus should be rotated by.
#The optional markupImage parameter will draw the photometric ellipses and the ellipse center on the current matplotlib figure
#in the matplotlib state machine if True.
#If set above 0 the optional verbose parameter will draw boxes around the sampled pixels. If set above 1, it will draw each subsample box
#All radii and angle values can be fractional.
#Subsample rate is the number of subsamples used for edge pixel weighting to fit along a single pixel in a single dimention.
#The number of tested subsamples grows as the square of this value
#Function returns the total electron flux interior to the target ellipse, calibrated by the sky background value of pixels exterior to the
#inner radius, but interior to the outer radius along with the uncertainty in the total electron flux. 
def photometry(data, X, Y, radX, radY, irX, irY, orX, orY, gain, angle=0, markupImage=True, verbose=0, subsampleRate=10):
    #draw the ellipses marking out the target, and the inner and outer regions of the sky annulus
    if markupImage:
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

    #create a grid of pixels to test for if they are inside the target ellipse
    innerSquareR = np.max([radX, radY])
    innerTestPts = np.mgrid[int(X-innerSquareR):int(X+innerSquareR), int(Y-innerSquareR):int(Y+innerSquareR)]
    #list of [x,y] pixel coordinates for pixels on the edge of the target ellipse
    if verbose>0 and markupImage:
        plt.gca().add_patch(patches.Rectangle([X-innerSquareR,Y-innerSquareR],2*innerSquareR,2*innerSquareR,color="blue"))

    targetEdge = []
    for t in np.linspace(0,2*np.pi,int(628*innerSquareR)):
        #the x and y coordinates for a bunch of points along the ellipse. 
        #Formulae derived using the parametric form of an ellipse and complex multiplication for rotation
        x = np.floor(radX*np.cos(t)*np.cos(angle)-radY*np.sin(t)*np.sin(angle)+X)
        y = np.floor(radY*np.sin(t)*np.cos(angle)+radX*np.cos(t)*np.sin(angle)+Y)
        #Only add new edge pixels to array
        if [x,y] not in targetEdge:
            targetEdge.append([x,y])
            if verbose>0: #mark edge pixel if requested
                plt.gca().add_patch(patches.Rectangle([x, y],1,1,color="darkgreen", fc=(0,0.5,0,0.3)))
    interiorPts = [] #[x,y] pixel coordinates for piexls inside the target ellipse
    #for each pixel in a grid around the target ellipse, test to see if they are inside the ellipse, and if they are
    #and they are not part of the edge, add them to the interior points array.
    for x in range(innerTestPts[0,0,0],innerTestPts[0,-1,0]):
        for y in range(innerTestPts[1,0,0],innerTestPts[1,0,-1]):
            x_transformed = ((x-X)*np.cos(-angle)-(y-Y)*np.sin(-angle))/radX
            y_transformed = ((y-Y)*np.cos(-angle)+(x-X)*np.sin(-angle))/radY
            if ((x_transformed)**2 + (y_transformed)**2) <= 1 and [x,y] not in targetEdge:
                interiorPts.append([x,y])
                if verbose>0:
                    plt.gca().add_patch(patches.Rectangle([x, y],1,1,color="purple",fc=(0.5, 0, 0.5, 0.3)))
    
    #Same series of steps as for the target ellipse, but for the sky annulus
    annulusR = np.max([orX, orY])
    annulusTestPts = np.mgrid[int(X-annulusR):int(X+annulusR), int(Y-annulusR):int(Y+annulusR)]
    annulusEdge=[]
    for t in np.linspace(0,2*np.pi,628*annulusR): #two sets of edges for an annulus
        xi = np.floor(irX*np.cos(t)*np.cos(angle)-irY*np.sin(t)*np.sin(angle)+X)
        yi = np.floor(irY*np.sin(t)*np.cos(angle)+irX*np.cos(t)*np.sin(angle)+Y)
        xo = np.floor(orX*np.cos(t)*np.cos(angle)-orY*np.sin(t)*np.sin(angle)+X)
        yo = np.floor(orY*np.sin(t)*np.cos(angle)+orX*np.cos(t)*np.sin(angle)+Y)
        if [xi, yi] not in annulusEdge:
            annulusEdge.append([xi,yi])
            if verbose>0:
                plt.gca().add_patch(patches.Rectangle([xi, yi],1,1,color="red",fc=(0.5, 0, 0, 0.3)))
        if [xo, yo] not in annulusEdge:
            annulusEdge.append([xo,yo])
            if verbose>0:
                plt.gca().add_patch(patches.Rectangle([xo, yo],1,1,color="blue",fc=(0, 0, 0.5, 0.3)))
    annulusInterior = []
    for x in range(annulusTestPts[0,0,0],annulusTestPts[0,-1,0]):
        for y in range(annulusTestPts[1,0,0],annulusTestPts[1,0,-1]):
            x_transformed1 = ((x-X)*np.cos(-angle)-(y-Y)*np.sin(-angle))/irX
            y_transformed1 = ((y-Y)*np.cos(-angle)+(x-X)*np.sin(-angle))/irY
            x_transformed2 = ((x-X)*np.cos(-angle)-(y-Y)*np.sin(-angle))/orX
            y_transformed2 = ((y-Y)*np.cos(-angle)+(x-X)*np.sin(-angle))/orY
            #pixels have to be outside the inner radius of the sky annulus, but inside the outer radius of the sky anulus
            if ((x_transformed1)**2 + (y_transformed1)**2) >= 1 \
                and ((x_transformed2)**2 + (y_transformed2)**2) <= 1 \
                and [x,y] not in annulusEdge:
                annulusInterior.append([x,y])
                if verbose>0:
                    plt.gca().add_patch(patches.Rectangle([x, y],1,1,color="cyan",fc=(0, 0.75, 0.75, 0.3)))


    #lists of the fraction of a pixel on the edge of an ellipse is inside a region
    targetEdgeWeights = []
    annulusEdgeWeights = []
    for coords in targetEdge:
        interior = 0 #number of subsamples inside the target ellipse
        for sx in range(subsampleRate):
            for sy in range(subsampleRate):
                #test to see if a subsample is inside the target ellipse, and if so, add it to the counter
                x_transformed = ((coords[0]+sx/subsampleRate-X)*np.cos(-angle)-(coords[1]+sy/subsampleRate-Y)*np.sin(-angle))/radX
                y_transformed = ((coords[1]+sy/subsampleRate-Y)*np.cos(-angle)+(coords[0]+sx/subsampleRate-X)*np.sin(-angle))/radY
                if ((x_transformed)**2 + (y_transformed)**2) <= 1:
                    interior += 1
                    if verbose>1:
                        plt.gca().add_patch(patches.Rectangle([coords[0]+sx/subsampleRate, coords[1]+sy/subsampleRate],1/subsampleRate,1/subsampleRate,color="black", fc=(0,0,0,0.3)))
        targetEdgeWeights.append(interior/subsampleRate**2) #weight is the number of subsamples inside the target ellipse divided by the total number of tested subsamples
    #essentially the same thing as above, but for the sky annulus instead of the target ellipse.
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
                    if verbose>1:
                        plt.gca().add_patch(patches.Rectangle([coords[0]+sx/subsampleRate, coords[1]+sy/subsampleRate],1/subsampleRate,1/subsampleRate,color="black", fc=(0,0,0,0.3)))
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

    targetValue = 0 #ADU
    targetPixels = 0
    #sum up the ADU values of the pixels corected for the sky brightness completely inside the target ellipse and keep track of the number of pixels
    for coords in interiorPts:
        targetValue += data[coords[1],coords[0]]-skyValue
        targetPixels += 1
    #do equivalent work for the partial pixels
    for i in range(len(targetEdge)):
        targetValue += (data[int(targetEdge[i][1]),int(targetEdge[i][0])]-skyValue) * targetEdgeWeights[i]
        targetPixels += targetEdgeWeights[i]
    targetValue *= gain #targetValue is now in e-, not ADU
    s_sky *= gain #s_sky is now in e-, not ADU
    s_skyValue *= gain #s_skyValue is now in e-, not ADU
    uncertainty = np.sqrt(targetValue+targetPixels*s_sky**2+targetPixels*s_skyValue**2) #uncertainty in e-

    return targetValue, uncertainty
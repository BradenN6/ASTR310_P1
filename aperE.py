# aperE.py

import numpy as np
import matplotlib.pyplot as plt

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

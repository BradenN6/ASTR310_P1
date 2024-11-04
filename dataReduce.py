from astropy.io import fits
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
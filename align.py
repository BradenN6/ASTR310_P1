
def findMaxPixelCoord(data, guessX, guessY, rX, rY):
    maxVal = float("-inf")
    maxCoord = (-1,-1)
    for x in range(int(guessX-rX), int(guessX+rX)):
        for y in range(int(guessY-rY), int(guessY + rY)):
            if data[y,x] > maxVal:
                maxVal = data[y,x]
                maxCoord = (x,y)
    return maxCoord[0], maxCoord[1]
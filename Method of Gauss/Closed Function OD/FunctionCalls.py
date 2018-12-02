# Maxwell Holloway
# Method of Gauss
from __future__ import division
import numpy as np
from math import sqrt, sin, cos, acos, degrees, radians, atan2, pi
from FuncFile import *

k = 0.01720209895
mu = 1

#Importing text file as an array
txt = np.loadtxt("MoGInputs.txt")
# Getting values from each item in the array
RAH = np.array([])
RAM = np.array([])
RAS = np.array([])
DecDeg = np.array([])
DecM = np.array([])
DecS = np.array([])
Year = np.array([])
Month = np.array([])
Day = np.array([])
tH = np.array([])
tM = np.array([])
tS = np.array([])
RX = np.array([])
RY = np.array([])
RZ = np.array([])
line = 0
while line < 3:
    RAH = np.append(RAH, txt.item(15 * line))
    RAM = np.append(RAM, txt.item(15 * line + 1))
    RAS = np.append(RAS, txt.item(15 * line + 2))
    DecDeg = np.append(DecDeg, txt.item(15 * line + 3))
    DecM = np.append(DecM, txt.item(15 * line + 4))
    DecS = np.append(DecS, txt.item(15 * line + 5))
    Year = np.append(Year, txt.item(15 * line + 6))
    Month = np.append(Month, txt.item(15 * line + 7))
    Day = np.append(Day, txt.item(15 * line + 8))
    tH = np.append(tH, txt.item(15 * line + 9))
    tM = np.append(tM, txt.item(15 * line + 10))
    tS = np.append(tS, txt.item(15 * line + 11))
    RX = np.append(RX, txt.item(15 * line + 12))
    RY = np.append(RY, txt.item(15 * line + 13))
    RZ = np.append(RZ, txt.item(15 * line + 14))
    line += 1
# Values for input to MoG
RA1Deg = RAHMStoDecDegrees(RAH[0], RAM[0], RAS[0])
Dec1Deg = DecDMStoDecDegrees(DecDeg[0], DecM[0], DecS[0])
RA2Deg = RAHMStoDecDegrees(RAH[1], RAM[1], RAS[2])
Dec2Deg = DecDMStoDecDegrees(DecDeg[1], DecM[1], DecS[1])
RA3Deg = RAHMStoDecDegrees(RAH[2], RAM[2], RAS[2])
Dec3Deg = DecDMStoDecDegrees(DecDeg[2], DecM[2], DecS[2])
RA1Rad = degrees(RA1Deg)
Dec1Rad = degrees(Dec1Deg)
RA2Rad = degrees(RA2Deg)
Dec2Rad = degrees(Dec2Deg)
RA3Rad = degrees(RA3Deg)
Dec3Rad = degrees(Dec3Deg)
RVec1 = np.array([RX[0], RY[0], RZ[0]])
RVec2 = np.array([RX[1], RY[1], RZ[1]])
RVec3 = np.array([RX[2], RY[2], RZ[2]])
# Time values of observations 1, 2, and 3
t1 = CivToJulian(tH[0], tM[0], tS[0], Year[0], Month[0], Day[0])
t2 = CivToJulian(tH[1], tM[1], tS[1], Year[1], Month[1], Day[1])
t3 = CivToJulian(tH[2], tM[2], tS[2], Year[2], Month[2], Day[2])
# Gaussian times using Gaussian constant, k
tau1 = k * (t1 - t2)
tau3 = k * (t3 - t2)
tau = k * (t3 - t1)

# Calling rhoHat function to acquire rhoHat values 1, 2, and 3 respectively
rhoHat1 = rhoHat(RA1Rad, Dec1Rad)
rhoHat2 = rhoHat(RA2Rad, Dec2Rad)
rhoHat3 = rhoHat(RA3Rad, Dec3Rad)


# Getting numbers for the Scalar Equation of Lagrange
A1 = tau3 / tau
A3 = (-tau1) / tau
B1 = (A1 / 6) * ((tau ** 2) - (tau3 ** 2))
B3 = (A3 / 6) * ((tau ** 2) - (tau1 ** 2)) 

# Getting D1j, D2j, and D3j values
D0 = np.dot(rhoHat1, np.cross(rhoHat2, rhoHat3))
D11 = np.dot(np.cross(RVec1, rhoHat2), rhoHat3)
D21 = np.dot(np.cross(rhoHat1, RVec1), rhoHat3)
D31 = np.dot(rhoHat1, np.cross(rhoHat2, RVec1))
D12 = np.dot(np.cross(RVec2, rhoHat2), rhoHat3)
D22 = np.dot(np.cross(rhoHat1, RVec2), rhoHat3)
D32 = np.dot(rhoHat1, np.cross(rhoHat2, RVec2))
D13 = np.dot(np.cross(RVec3, rhoHat2), rhoHat3)
D23 = np.dot(np.cross(rhoHat1, RVec3), rhoHat3)
D33 = np.dot(rhoHat1, np.cross(rhoHat2, RVec3))


# necessary to find coefficients of scalar equation of lagrange
A = (A1 * D21 - D21 + A3 * D23) / (-D0)
B = (B1 * D21 + B3 * D23) / (-D0)
F = np.linalg.norm(RVec2) ** 2
E = (-2) * np.dot(rhoHat2, RVec2)

# coefficients of the scalar equation of lagrange; necessary to find roots
a = -(A ** 2 + A * E + F)
b = (-mu) * (2 * A * B + B * E)
c = -((mu ** 2) * (B ** 2))

# get realRoots values from rootSolve function; necessary for the following for loop
realRoots = rootSolve(a, b, c)
# getting empty arrays to allow for appending at the end of the for loop
r2VecAr = np.array([])
r2DotVec = np.array([])
r2DotMags = np.array([])
rhoMags = np.array([])
f1Ar = np.array([])
f3Ar = np.array([])
g1Ar = np.array([])
g3Ar = np.array([])
i = 0
for realRoot in realRoots:
    # calculating truncated f and g functions
    u = 1 / (realRoot ** 3)
    fT1 = 1 - (1/2) * u * tau1 ** 2
    fT3 = 1 - (1/2) * u * tau3 ** 2
    gT1 = tau1 - (1 / 6) * u * (tau1) ** 3
    gT3 = tau1 - (1 / 6) * u * (tau1) ** 3

    # calculating c values give then f and g functions
    c1 = gT3 / (fT1 * gT3 - gT1 * fT3)
    c2 = -1
    c3 = (-gT1) / (fT1 * gT3 - gT1 * fT3)

    # now calculate values for rho magnitude (scalar equations of range)
    rhoMag1 = (c1 * D11 + c2 * D12 + c3 * D13) / (c1 * D0)
    rhoMag2 = (c1 * D21 + c2 * D22 + c3 * D23) / (c2 * D0)
    rhoMag3 = (c1 * D31 + c2 * D32 + c3 * D33) / (c3 * D0)
    rhoMags = np.append(rhoMags, [rhoMag1, rhoMag2, rhoMag3])

    # get r1 and r3 position vectors in order to get r2 dot vector
    r1Vec = rhoMag1 * rhoHat1 - RVec1
    r3Vec = rhoMag3 * rhoHat3 - RVec3
    r2Vec = c1 * r1Vec + c2 * r3Vec
    r2VecAr = np.append(r2VecAr, r2Vec)

    # get d1 and d3 values
    d1 = (-fT3) / (fT1 * gT3 - fT3 * gT1)
    d3 = fT1 / (fT1 * gT3 - fT3 * gT1)
    
    # get r2 velocity vector and feel like a badass
    r2DotVec = np.append(r2DotVec, d1 * r1Vec + d3 * r3Vec)
    r2DotMags = np.append(r2DotMags, np.linalg.norm(r2DotVec.item(i)))
    r2Mags = rootSolve(a, b, c)
    EA1 = EAFunc(r2Mags, r2DotMags, r2DotVec, tau1)
    EA3 = EAFunc(r2Mags, r2DotMags, r2DotVec, tau3)
    f1, g1 = closedFuncs(r2Mags, EA1, tau1, r2DotVec)
    f3, g3 = closedFuncs(r2Mags, EA3, tau3, r2DotVec)
    i += 1


# initial values for the iterative function
times = np.array([t1, t2, t3])
VLA = np.array([r2VecAr, r2Mags, r2DotVec, r2DotMags, rhoMags, times])
print VLA
print VLA[4]
Rho1 = 2000000000 # large value that will obviously have a large difference from the actual rhoNorm
oRho1 = 5
counter = 0
while abs(Rho1 - oRho1) > 1:
    oVLA = VLA
    oRhoMags = oVLA[4]
    oRho1 = oRhoMags[0]
    VLA = itFunc(oVLA.item(0), oVLA.item(1), oVLA.item(2), oVLA.item(3), oVLA.item(4), oVLA.item(5))
    RhoMags = VLA[4]
    Rho1 = RhoMags[0]
    counter += 1
    print counter
    print "VLA", i, VLA


r2VecArEq = VLA[0][0]
r2DotVecEq = VLA[0][2]
# equatorial to ecliptic transform before calculating orbital elements
r2VecArEc = EqToEc(r2VecArEq)
xPos, yPos, zPos = r2VecArEc[0][0], r2VecArEc[0][1], r2VecArEc[0][2]
r2DotVecEc = EqToEc(r2DotVecEq)
xVel, yVel, zVel = r2DotVecEc[0][0], r2DotVecEc[0][1], r2DotVecEc[0][2]

print "ORBITAL ELEMENTS!!!!!!!!!!!!!!!!!!!!!", OrbElementz(xPos, yPos, zPos, xVel, yVel, zVel)    


        


##if len(realRoots) == 3:
##    r2VecsAll = VLA.item(0)
##    r2Vecs1 = np.array([r2VecsAll.item(0), r2VecsAll.item(3), r2VecsAll.item(6)])
##    r2Vecs2 = np.array([r2VecsAll.item(1), r2VecsAll.item(4), r2VecsAll.item(7)])
##    r2Vecs3 = np.array([r2VecsAll.item(2), r2VecsAll.item(5), r2VecsAll.item(8)])
##    r2Vec11, r2Vec12, r2Vec13 = r2Vecs1[0][0], r2Vecs1[0][1], r2Vecs1[0][2]
##    r2Vec21, r2Vec22, r2Vec23 = r2Vecs2[0][0], r2Vecs2[0][1], r2Vecs2[0][2]
##    r2Vec31, r2Vec32, r2Vec33 = r2Vecs3[0][0], r2Vecs3[0][1], r2Vecs3[0][2]
##    for i in range(3):
##        rhoMagOld = VLA.item(4)
##        rhoMag1 = 20000
##        while abs(rhoMag1Old - rhoMag1) > 10 ** -1:
##            VLA = itFunc(VLA.item(0), VLA.item(1), VLA.item(2), VLA.item(3), VLA.item(4), VLA.item(5), VLA.item(6))
##            rhoMag1 = VLA.item(4)
##        
    
    

    
##elif len(realRoots) == 2:
##elif len(realRoots) == 1:
##else:
##    print "you're shit out of luck; no roots; I'll just guess 1.5 AU"
##    realRoots[0][0] = 1.5

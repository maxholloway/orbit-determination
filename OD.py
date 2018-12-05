# Maxwell Holloway
# Method of Gauss
from __future__ import division
import numpy as np
from math import sqrt, sin, cos, acos, degrees, atan2, pi
from Functions import *

k = 0.01720209895
mu = 1

#Importing text file as an array
data = np.genfromtxt("Inputs.txt", dtype = None)
times, RARadAr, DecRadAr, Year, Month, Day, RX, RY, RZ = readInput(data)

tH = np.array([])
tM = np.array([])
tS = np.array([])
for i in range(len(times)):
    tH = np.append(tH, times[i][0])
    tM = np.append(tM, times[i][1])
    tS = np.append(tS, times[i][2])
RA1Rad = RARadAr[0]
Dec1Rad = DecRadAr[0]
RA2Rad = RARadAr[1]
Dec2Rad = DecRadAr[1]
RA3Rad = RARadAr[2]
Dec3Rad = DecRadAr[2]
RA4Rad = RARadAr[3]
Dec4Rad = DecRadAr[3]
RVec1 = np.array([RX[0], RY[0], RZ[0]]) # good
RVec2 = np.array([RX[1], RY[1], RZ[1]]) # good
RVec3 = np.array([RX[2], RY[2], RZ[2]]) # good
RVec4 = np.array([RX[3], RY[3], RZ[3]]) # good






observation = input("Observation 2 or observation 3 as center? ")


if observation == 2:
    t2 = CivToJulian(tH[1], tM[1], tS[1], Year[1], Month[1], Day[1])
    rhoHat2 = EqToEc(rhoHat(RA2Rad, Dec2Rad))

elif observation == 3:
    t2 = CivToJulian(tH[2], tM[2], tS[2], Year[2], Month[2], Day[2])
    rhoHat2 = EqToEc(rhoHat(RA3Rad, Dec3Rad))
    RVec2 = RVec3

t1 = CivToJulian(tH[0], tM[0], tS[0], Year[0], Month[0], Day[0])
t3 = CivToJulian(tH[3], tM[3], tS[3], Year[3], Month[3], Day[3])
rhoHat1 = EqToEc(rhoHat(RA1Rad, Dec1Rad))
rhoHat3 = EqToEc(rhoHat(RA4Rad, Dec4Rad))

# Gaussian times using Gaussian constant, k
tau1 = k * (t1 - t2)
tau3 = k * (t3 - t2)
tau = k * (t3 - t1)

# Getting numbers for the Scalar Equation of Lagrange
A1 = tau3 / tau
A3 = (-tau1) / tau
B1 = (A1 / 6) * ((tau ** 2) - (tau3 ** 2))
B3 = (A3 / 6) * ((tau ** 2) - (tau1 ** 2))  

# Getting D1j, D2j, and D3j values
D0 = np.dot(rhoHat1, np.cross(rhoHat2, rhoHat3))
D11 = np.dot(np.cross(RVec1, rhoHat2), rhoHat3)
D12 = np.dot(np.cross(RVec2, rhoHat2), rhoHat3)
D13 = np.dot(np.cross(RVec4, rhoHat2), rhoHat3)
D21 = np.dot(np.cross(rhoHat1, RVec1), rhoHat3)
D22 = np.dot(np.cross(rhoHat1, RVec2), rhoHat3)
D23 = np.dot(np.cross(rhoHat1, RVec4), rhoHat3)
D31 = np.dot(rhoHat1, np.cross(rhoHat2, RVec1))
D32 = np.dot(rhoHat1, np.cross(rhoHat2, RVec2))
D33 = np.dot(rhoHat1, np.cross(rhoHat2, RVec4))


# necessary to find coefficients of scalar equation of lagrange
A = (A1 * D21 - D22 + A3 * D23) / (-D0)
B = (B1 * D21 + B3 * D23) / (-D0)
F = np.linalg.norm(RVec2) ** 2
E = (-2) * np.dot(rhoHat2, RVec2)

# coefficients of the scalar equation of lagrange; necessary to find roots
a = -(A ** 2 + A * E + F)
b = (-mu) * (2 * A * B + B * E)
c = -((mu ** 2) * (B ** 2))

# get realRoots values from rootSolve function; necessary for the following for loop
realRoots = rootSolve(a, b, c)
print "realRoots", realRoots
# getting empty arrays to allow for appending at the end of the for loop
r2VecAr = np.array([])
r2DotVec = np.array([])
r2DotMags = np.array([])
rhoMags = np.array([])
f1Ar = np.array([])
f3Ar = np.array([])
g1Ar = np.array([])
g3Ar = np.array([])
realRoot = realRoots

# calculating truncated f and g functions
u = 1 / (realRoot ** 3)
fT1 = 1 - (1/2) * u * tau1 ** 2
fT3 = 1 - (1/2) * u * tau3 ** 2
gT1 = tau1 - (1 / 6) * u * (tau1) ** 3
gT3 = tau3 - (1 / 6) * u * (tau3) ** 3

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
r3Vec = rhoMag3 * rhoHat3 - RVec4
r2Vec = rhoMag2 * rhoHat2 - RVec2

# get d1 and d3 values
d1 = (-fT3) / (fT1 * gT3 - fT3 * gT1)
d3 = fT1 / (fT1 * gT3 - fT3 * gT1)

# get r2 velocity vector and feel like a badass
r2DotVec = d1 * r1Vec + d3 * r3Vec
r2DotMags = np.append(r2DotMags, np.linalg.norm(r2DotVec))
r2Mags = realRoot




times = np.array([t1, t2, t3])
VLA = np.array([r2Vec, r2Mags, r2DotVec, r2DotMags, rhoMags, times])
Rho1 = 2000000000 # large value that will obviously have a large difference from the actual rhoNorm
oRho1 = 5
while abs(oRho1 - Rho1) > 1*10**(-4):
    #print "counter", counter
    oVLA = VLA    
    oRhoMags = oVLA[4]
    oRho1 = oRhoMags[0]
    VLA = itFuncTay(oVLA[0], oVLA[1], oVLA[2], oVLA[3], oVLA[4], oVLA[5])
    RhoMags = VLA[4]
    Rho1 = RhoMags[0]
    
r2VecArEc = VLA[0]
r2DotVecEc = VLA[2]
xPos, yPos, zPos = r2VecArEc[0], r2VecArEc[1], r2VecArEc[2]
xVel, yVel, zVel = r2DotVecEc[0], r2DotVecEc[1], r2DotVecEc[2]

print "\nPosition:", r2VecArEc
print "\nVelocity (AU/Day):", r2DotVecEc * k
print "\nRange:", VLA[4][1]
print "\nORBITAL ELEMENTS!!!!!!!!!!!!!!!!!!!!!", OrbElementz(xPos, yPos, zPos, xVel, yVel, zVel)






# Maxwell Holloway
# Method of Gauss
from __future__ import division
import numpy as np
from math import sqrt, sin, cos, acos, degrees, atan2, pi
from Functions import *


def itFuncTay(r2Vec, r2Mag, r2DotVec, r2DotMags, rhoMags, times):
    # do time correction
    times = [tCorrection(times[i], rhoMags[i]) for i in range(len(times))]

    tau1 = k * (times[0] - times[1])
    tau3 = k * (times[2] - times[1])
    tau = k * (times[2] - times[0])

    f1, g1 = FnGTaylor(tau1, r2Mag, r2Vec, r2DotVec) # calling function to get f1 and g1 closed functions
    f3, g3 = FnGTaylor(tau3, r2Mag, r2Vec, r2DotVec) # calling function to get f3 and g3 closed functions
    
    # finding new c values for iterations
    c1 = g3 / (f1 * g3 - g1 * f3)
    c2 = -1
    c3 = (-g1) / (f1 * g3 - g1 * f3)
 
    #calculating rho magnitude; D values should not change much with time, so not re-doing D values; GET CHECKED
    rhoMag1 = (c1 * D11 + c2 * D12 + c3 * D13) / (c1 * D0)
    rhoMag2 = (c1 * D21 + c2 * D22 + c3 * D23) / (c2 * D0)
    rhoMag3 = (c1 * D31 + c2 * D32 + c3 * D33) / (c3 * D0)
    rhoMags = np.array([rhoMag1, rhoMag2, rhoMag3])
    
    # get r1 and r3 position vectors in order to get r2 dot vector
    rVecs = [rhoMags[i]*rhoHat[i]-RVecs[i]]
    # r1Vec = rhoMag1 * rhoHat1 - RVec1
    # r2Vec = rhoMag2 * rhoHat2 - RVec2
    # r3Vec = rhoMag3 * rhoHat3 - RVec4

    # get d1 and d3 values
    d1 = (-f3) / (f1 * g3 - f3 * g1)
    d3 = f1 / (f1 * g3 - f3 * g1)
    
    # get r2 velocity vector and feel like a badass
    r2Mags = np.linalg.norm(rVecs[1])
    r2DotVec = d1 * rVecs[0] + d3 * rVecs[2]
    r2DotMags = np.linalg.norm(r2DotVec)

    xPos, yPos, zPos = r2Vec[0], r2Vec[1], r2Vec[2]
    xVel, yVel, zVel = r2DotVec[0], r2DotVec[1], r2DotVec[2]
    output = np.array([r2Vec, r2Mags, r2DotVec, r2DotMags, rhoMags, times])
    return output


k = 0.01720209895
mu = 1

#Importing text file as an array
data = np.genfromtxt("Inputs.txt", dtype = None)
times, RARadAr, DecRadAr, Year, Month, Day, RX, RY, RZ = readInput(data)

tH = []
tM = []
tS = []
for i in range(len(times)):
    tH += times[i][0]
    tM += times[i][1]
    tS += times[i][2]

RA1Rad, RA2Rad, RA3Rad, RA4Rad = RARadAr
Dec1Rad, Dec2Rad, Dec3Rad, Dec4Rad = DecRadAr

RVec1 = np.array([RX[0], RY[0], RZ[0]]) 
RVec2 = np.array([RX[1], RY[1], RZ[1]]) 
RVec3 = np.array([RX[2], RY[2], RZ[2]]) 
RVec4 = np.array([RX[3], RY[3], RZ[3]]) 

observation = input("Observation 2 or observation 3 as center? ")
assert(observation in [2, 3]) # raise an error if the observation is invalid
RVecMid = 0 # initialized in the if statements
rhoHat2 = 0 # initialized in the if statements
if observation == 2:
    t2 = CivToJulian(tH[1], tM[1], tS[1], Year[1], Month[1], Day[1])
    rhoHat2 = EqToEc(rhoHat(RA2Rad, Dec2Rad))
    RVecMid = RVec2
else:
    t2 = CivToJulian(tH[2], tM[2], tS[2], Year[2], Month[2], Day[2])
    rhoHat2 = EqToEc(rhoHat(RA3Rad, Dec3Rad))
    RVecMid = RVec3

# Conversions
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
D = [[0 for i in range(3)] for j in range(3)]
RVecs = [RVec1, RVecMid, RVec4]
D[0] = [np.dot(np.cross(RVecs[i], rhoHat2), rhoHat3) for i in range(3)]
D[1] = [np.dot(np.cross(rhoHat1, RVecs[i]), rhoHat3) for i in range(3)]
D[2] = [np.dot(rhoHat1, np.cross(rhoHat2, RVecs[i])) for i in range(3)]

D21 = np.dot(np.cross(rhoHat1, RVec1), rhoHat3)
D22 = np.dot(np.cross(rhoHat1, RVecMid), rhoHat3)
D23 = np.dot(np.cross(rhoHat1, RVec4), rhoHat3)
D31 = np.dot(rhoHat1, np.cross(rhoHat2, RVec1))
D32 = np.dot(rhoHat1, np.cross(rhoHat2, RVecMid))
D33 = np.dot(rhoHat1, np.cross(rhoHat2, RVec4))

# necessary to find coefficients of scalar equation of lagrange
A = (A1 * D[1][0] - D[1][1] + A3 * D[1][2]) / (-D0)
B = (B1 * D[1][0] + B3 * D[1][2]) / (-D0)
F = np.linalg.norm(RVecMid) ** 2
E = (-2) * np.dot(rhoHat2, RVecMid)

# coefficients of the scalar equation of lagrange; necessary to find roots
a = -(A ** 2 + A * E + F)
b = (-mu) * (2 * A * B + B * E)
c = -((mu ** 2) * (B ** 2))

# get realRoots values from rootSolve function; necessary for the following for loop
realRoot = rootSolve(a, b, c)

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
cs = [c1, c2, c3]

# now calculate values for rho magnitude (scalar equations of range)
rhoMags = [sum([cs[j]*D[i][j] for j in range(3)])/(cs[i]*D0) for i in range(3)]
# rhoMag1 = (c1 * D[0][0] + c2 * D[0][1] + c3 * D[0][2]) / (c1 * D0)
# rhoMag2 = (c1 * D[1][0] + c2 * D[1][1] + c3 * D[1][2]) / (c2 * D0)
# rhoMag3 = (c1 * D[2][0] + c2 * D[2][1] + c3 * D[2][2]) / (c3 * D0)
# rhoMags = np.array([rhoMag1, rhoMag2, rhoMag3])

rhoHats = [rhoHat1, rhoHat2, rhoHat3]
# get r1 and r3 position vectors in order to get r2 dot vector
r1Vec, r2Vec, r3Vec = [rhoMags[i]*rhoHats[i] - RVecs[i]]
# r1Vec = rhoMag1 * rhoHat1 - RVec1
# r2Vec = rhoMag2 * rhoHat2 - RVecMid
# r3Vec = rhoMag3 * rhoHat3 - RVec4

# get d1 and d3 values
d1 = (-fT3) / (fT1 * gT3 - fT3 * gT1)
d3 = fT1 / (fT1 * gT3 - fT3 * gT1)

# get r2 velocity vector and feel like a badass
r2DotVec = d1 * r1Vec + d3 * r3Vec
r2DotMags = np.append(r2DotMags, np.linalg.norm(r2DotVec))
r2Mag = realRoot




times = np.array([t1, t2, t3])
Rho1 = 2000000000 # large value that will obviously have a large difference from the actual rhoNorm
oRho1 = 5
while abs(oRho1 - Rho1) > 10**(-4):
    oRho1 = rhoMags[0]
    r2Vec, r2Mag, r2DotVec, r2DotMags, rhoMags, times = itFuncTay(r2Vec, r2Mag, r2DotVec, r2DotMags, rhoMags, times)
    Rho1 = rhoMags[0]

xPos, yPos, zPos = r2Vec[0], r2Vec[1], r2Vec[2]
xVel, yVel, zVel = r2DotVec[0], r2DotVec[1], r2DotVec[2]

print("\nPosition:", r2VecAr)
print("\nVelocity (AU/Day):", r2DotVec * k)
print("\nRange:", VLA[4][1])
print("\nORBITAL ELEMENTS!!!!!!!!!!!!!!!!!!!!!", OrbElementz(xPos, yPos, zPos, xVel, yVel, zVel))






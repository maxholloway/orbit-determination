# Maxwell Holloway
# Method of Gauss
from __future__ import division
import numpy as np
from math import sqrt, sin, cos, acos, degrees, atan2, pi
from Functions import *


k = 0.01720209895
mu = 1

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
 
    #calculating rho magnitude; D values should not change much with time, so not re-doing D values
    rhoMags = [sum([cs[j]*D[i][j] for j in range(3)])/(cs[i]*D0) for i in range(3)]
    
    # get r1 and r3 position vectors in order to get r2 dot vector
    rVecs = [(rhoMags[i]*rhoHats[i]-RVecs[i]) for i in range(3)]

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

def OrbElementz(xPos, yPos, zPos, xVel, yVel, zVel):

    answer = 2
    if answer == 2:
        t2 = CivToJulian(tH[1], tM[1], tS[1], Year[1], Month[1], Day[1])
        rhoHat2 = EqToEc(rhoHat(RA2Rad, Dec2Rad))
    elif answer == 3:
        t2 = CivToJulian(tH[2], tM[2], tS[2], Year[2], Month[2], Day[2])
        rhoHat2 = EqToEc(rhoHat(RA3Rad, Dec3Rad))
    
    posVec = np.array([xPos, yPos, zPos])
    posDot = np.dot(posVec, posVec)
    posMag = np.linalg.norm(posVec)
    velVec = np.array([xVel, yVel, zVel])
    velDot = np.dot(velVec, velVec)

    hVec = np.cross(posVec, velVec)
    hMag = np.linalg.norm(hVec)
    hx = hVec.item(0)
    hy = hVec.item(1)
    hz = hVec.item(2)

    # Calculating orbital element: a
    a = (2 / np.linalg.norm(posVec) - velDot / mu) ** -1

    # Calculating orbital element: e
    e = sqrt((1 - (np.linalg.norm(np.cross(posVec, velVec)) ** 2) / (mu * a)))
    # Calculating orbital element: i
    i = acos(hz / hMag)

    # Calculating orbital element: Omega
    sinOmega = hx / (hMag * sin(i))
    cosOmega = -hy / (hMag * sin(i))
    Omega = atan2(sinOmega, cosOmega)
    if Omega < 0:
        Omega = Omega + 2 * pi

    # Calculating orbital element: omega
    sinfomega = (zPos / (posMag * sin(i)))
    cosfomega = (1 / cos(Omega)) * ((xPos / posMag) + cos(i) * (zPos / (posMag * sin(i))) * sin(Omega))
    fomega = atan2(sinfomega, cosfomega)
    cosf = (1 / e) * (((a * (1 - e ** 2)) / posMag) - 1)
    sinf = (1 / e) * sqrt(a * (1 - e ** 2)) * (np.dot(posVec, velVec) / posMag)
    f = atan2(sinf, cosf)
    if f < 0:
        f = f + 2 * pi
    omega = fomega - f
    omega = omega % (2 * pi)

    # Calculating orbital element: M
    E = acos((1 / e) * (1 - posMag / a))
    M = E - e * sin(E)
    if f < pi and M > pi:
        M = (2 * pi) - M
    elif f > pi and M < pi:
        M = (2 * pi) - M
    t0 = CivToJulian(6, 0, 0, 2017, 7, 22)
    n = sqrt(mu / (a**3))
    M = M + n * k * (t0 - t2)
    M = M % (2*pi)
    return(a, e, degrees(i), degrees(Omega), degrees(omega), degrees(M))

#Importing text file as an array
data = np.genfromtxt("Inputs.txt", dtype = str)

times, RARadAr, DecRadAr, Year, Month, Day, RX, RY, RZ = readInput(data)
tH = []
tM = []
tS = []
for i in range(len(times)):
    tH += [times[i][0]]
    tM += [times[i][1]]
    tS += [times[i][2]]

RA1Rad, RA2Rad, RA3Rad, RA4Rad = RARadAr
Dec1Rad, Dec2Rad, Dec3Rad, Dec4Rad = DecRadAr

RVec1 = np.array([RX[0], RY[0], RZ[0]]) 
RVec2 = np.array([RX[1], RY[1], RZ[1]]) 
RVec3 = np.array([RX[2], RY[2], RZ[2]]) 
RVec4 = np.array([RX[3], RY[3], RZ[3]]) 

observation = input("Observation 2 or observation 3 as center? ")
assert(int(observation) in [2, 3]) # raise an error if the observation is invalid
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

rhoHats = [rhoHat1, rhoHat2, rhoHat3]
# get r1 and r3 position vectors in order to get r2 dot vector
r1Vec, r2Vec, r3Vec = [rhoMags[i]*rhoHats[i] - RVecs[i] for i in range(3)]

# get d1 and d3 values
d1 = (-fT3) / (fT1 * gT3 - fT3 * gT1)
d3 = fT1 / (fT1 * gT3 - fT3 * gT1)

# get r2 velocity vector and feel like a badass
r2DotVec = d1 * r1Vec + d3 * r3Vec
r2DotMags = np.linalg.norm(r2DotVec)
r2Mag = realRoot


times = np.array([t1, t2, t3])
Rho1 = 2000000000 # large value that will obviously have a large difference from the actual rhoNorm
oRho1 = 5
while abs(oRho1 - Rho1) > 10**(-4):
    oRho1 = rhoMags[0]
    r2Vec, r2Mag, r2DotVec, r2DotMags, rhoMags, times = itFuncTay(r2Vec, r2Mag, r2DotVec, r2DotMags, rhoMags, times)
    Rho1 = rhoMags[0]

xPos, yPos, zPos = r2Vec
xVel, yVel, zVel = r2DotVec

print("\nPosition:", r2Vec)
print("\nVelocity (AU/Day):", r2DotVec * k)
print("\nRange:", rhoMags[1])
print("\nORBITAL ELEMENTS!!!!!!!!!!!!!!!!!!!!!", OrbElementz(xPos, yPos, zPos, xVel, yVel, zVel))






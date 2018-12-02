
# Maxwell Holloway
# Method of Gauss
from __future__ import division
import numpy as np
from math import sqrt, sin, cos, acos, degrees, radians, atan2, pi
from FuncFile import RAH

k = 0.01720209895
mu = 1 # get this checked

#Importing text file as an array
txt = np.loadtxt("MoGInputs.txt")
# Getting values from each item in the array
RAH = np.array([])
RAM = np.array([])
RAS = np.array([])
DecDeg = np.array([])
DecM = np.array([])
DecS = np.array([])
tH = np.array([])
tM = np.array([])
tS = np.array([])
RX = np.array([])
RY = np.array([])
RZ = np.array([])
line = 0
while line < 3:
    RAH = np.append(RAH, txt.item(12 * line))
    RAM = np.append(RAM, txt.item(12 * line + 1))
    RAS = np.append(RAS, txt.item(12 * line + 2))
    DecDeg = np.append(DecDeg, txt.item(12 * line + 3))
    DecM = np.append(DecM, txt.item(12 * line + 4))
    DecS = np.append(DecS, txt.item(12 * line + 5))
    tH = np.append(tH, txt.item(12 * line + 6))
    tM = np.append(tM, txt.item(12 * line + 7))
    tS = np.append(tS, txt.item(12 * line + 8))
    RX = np.append(RX, txt.item(12 * line + 9))
    RY = np.append(RY, txt.item(12 * line + 10))
    RZ = np.append(RZ, txt.item(12 * line + 11))
    line += 1
# Julian days function
def CivToJulian(tUTHour, tUTMin, tUTSec, Y, M, D):
    tUTMinDecimal = tUTMin + (tUTSec / 60)
    UTDecimal = tUTHour + (tUTMinDecimal / 60)
    a = 367 * Y
    b = int((M + 9)/12)
    c = 7 * (Y + b)
    d = int(c / 4)
    e = int((275 * M) / 9)
    J0 = a - d + e + D + 1721013.5
    return J0

# Functions to convert to RA and Dec decimal degrees
def DecDMStoDecDegrees(D, M, S):
    if D < 0:
        M = -M
        S = -S
    MDec = M + (S / 60)
    DecDegrees = D + (MDec / 60)
    return DecDegrees

def RAHMStoDecDegrees(H, M, S):
    if H < 0:
        M = -M
        S = -S
    MDec = M + (S / 60)
    HDec = H + (MDec / 60)
    DecDegrees = HDec * 15
    return DecDegrees
def rootSolve(a, b, c):
    coeffs = np.array([1, 0, a, 0, 0, b, 0, 0, c])
    roots = np.roots(coeffs)
    realRootsOld = roots[np.isreal(roots)]
    realRootsAll = np.array([])
    i = 0 # dummy variable to use in the following for loop because its elements are not from 0 to some number
    for element in realRootsOld:
        h = realRootsOld.item(i)
        realRootsAll= np.append(realRootsAll, abs(h))
        i += 1
    realRoots = np.unique(realRootsAll) # used to get rid of negative roots that were converted to positive roots and then had the same value as an already existing root value
    realRoots = np.array([realRoots.item(0)])
    return realRoots
# Function for finding rhoHat
def rhoHat(RARad, decRad):
    rhoHat = np.array([cos(RARad)*cos(decRad), sin(RARad) * cos(decRad), sin(decRad)])
    return rhoHat
def EAFunc(r2Mag, r2DotVec, TAUi):
    for element in range(len(r2Mags)):
        # first must find a in order and get n
        dotProdr2Dot = np.dot(r2DotVec.item(element), r2DotVec.item(element))
        r2Mag = r2Mags.item(element)
        a = abs((2 / r2Mag - dotProdr2Dot) ** (-1))
        n = k * sqrt(1 / (a ** 3))
        
        # getting items from arrays to make calculations cleaner
        r2DotMag = r2DotMags.item(element)

        # starting iterations now!
        M = n * (TAUi) # get this checked; should give elements for time of time at tau3
        EA = M
        absVar = 1 # arbitrary value > 10 ^ -5
        counter = 0
        while absVar > 10 ** -2: #### GET THIS CHECKED; IT'S WRONG; must iterate more
            EOld = EA
            fx = EA - (1 - r2Mag / a) * sin(EA) + (r2Mag * r2DotMag / (n * (a ** 2))) * (1 - cos(EA)) - (n * TAUi) 
            fpx = 1 - (1 - r2Mag / a) * cos(EA) + (r2Mag * r2DotMag / (n * (a ** 2))) * sin(EA)
            EA = EA - fx / fpx
            absVar = abs(EOld - EA)
            counter += 1
            print counter
    return EA
# functions for getting closed function form of f and g
def closedFuncs(a, r2Mags, EA, TAUi):
    k = .01720209894
    a = abs(a) # FIX THIS; 'a' WAS NEGATIVE AND SO I DID THIS AS A QUICK FIX
    n = k * sqrt(1 / (a ** 3))
    f = np.array([])
    g = np.array([])
    i = 0
    for r2Mag in r2Mags:
        r2Mag = r2Mags.item(i)
        f = np.append(f, 1 - (a / r2Mag) * (1 - cos(EA)))
        g = np.append(g, TAUi + (1 / n) * (sin(EA) - EA))
        i += 1
    return f, g
# function for time correction
def tCorrection(JTimei, rhoMagi):
    c = 173 # speed of light in AU per day
    tCorrection = JTimei - rhoMagi / c
    return tCorrection
def itFunc(r2Vec, r2Mags, r2DotVec, r2DotMags, rhoMag1, rhoMag2, rhoMag3):
    # do time correction
    tCorrected1 = tCorrection(t1, rhoMag1)
    tCorrected2 = tCorrection(t2, rhoMag2)
    tCorrected3 = tCorrection(t3, rhoMag3)
    tau1 = k * (tCorrected1 - tCorrected2)
    tau3 = k * (tCorrected3 - tCorrected2)
    tau = k * (tCorrected3 - tCorrected1)
    # jump through some hoops to get f and g functions
    for r2Mag in r2Mags:
        EA1 = EAFunc(r2Mags, r2DotVec, tau1) #calling function to get deltaE for the first interval
        EA3 = EAFunc(r2Mags, r2DotVec, tau3) # calling function to get deltaE for the third interval
        f1CFunc, g1CFunc = closedFuncs(a, r2Mags, EA1, tau1) # calling function to get f1 and g1 closed functions
        f3CFunc, g3CFunc = closedFuncs(a, r2Mags, EA3, tau3) # calling function to get f3 and g3 closed functions
        ##f1Ar, g1Ar = np.append(f1Ar, f1CFunc), np.append(g1Ar, g1CFunc)
        ##f3Ar, g3Ar = np.append(f3Ar, f3CFunc), np.append(g3Ar, g3CFunc)

        # finding new c values for iterations
        c1 = g3CFunc / (f1CFunc * g3CFunc - g1CFunc * f3CFunc)
        c2 = -1
        c3 = (-g1CFunc) / (f1CFunc * g3CFunc - g1CFunc * f3CFunc)

        #calculating rho magnitude; D values should not change much with time, so not re-doing D values; GET CHECKED
        rhoMag1 = (c1 * D11 + c2 * D12 + c3 * D13) / (c1 * D0)
        rhoMag2 = (c1 * D21 + c2 * D22 + c3 * D23) / (c2 * D0)
        rhoMag3 = (c1 * D31 + c2 * D32 + c3 * D33) / (c3 * D0)
        
        # get r1 and r3 position vectors in order to get r2 dot vector
        r1Vec = rhoMag1 * rhoHat1 - RVec1
        r3Vec = rhoMag3 * rhoHat3 - RVec3

        # get r2 vector from r1 and r3 vectors
        r2Vec = c1 * r1Vec + c2 * r3Vec

        # get d1 and d3 values
        d1 = (-fT3) / (fT1 * gT3 - fT3 * gT1)
        d3 = fT1 / (fT1 * gT3 - fT3 * gT1)
        
        # get r2 velocity vector and feel like a badass
        r2DotVec = np.append(r2DotVec, d1 * r1Vec + d3 * r3Vec)
        r2DotMags = np.append(r2DotMags, np.linalg.norm(r2DotVec.item(i)))
        output = np.array([r2Vec, r2Mags, r2DotVec, r2DotMags, rhoMag1, rhoMag2, rhoMag3])
    return output

def EcToEq(EcVec):
    ob = radians(23.43701)
    EqVec = np.dot(np.array([[1, 0, 0], [0, cos(ob), sin(ob)], [0, -sin(ob), cos(ob)]]), EcVec)
    return EqVec

def EqToEc(EqVec):
    ob = radians(23.43701)
    EcVec = np.dot(np.array([[1, 0, 0], [0, cos(ob), sin(ob)], [0, -sin(ob), cos(ob)]]), EqVec)
    return EcVec 

def OrbElementz(xPos, yPos, zPos, xVel, yVel, zVel):
    mu = 1
    
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
    return(e, a, i, Omega, omega, M)































# Values for input to MoG
RA1Deg = RAHMStoDecDegrees(15, 26, 34.47)
Dec1Deg = DecDMStoDecDegrees(-8, 56, 42)
RA2Deg = RAHMStoDecDegrees(15, 19, 55.75)
Dec2Deg = DecDMStoDecDegrees(-5, 28, 47.8)
RA3Deg = RAHMStoDecDegrees(15, 16, 49.82)
Dec3Deg = DecDMStoDecDegrees(-2, 53, 38.2)
RA1Rad = degrees(RA1Deg)
Dec1Rad = degrees(Dec1Deg)
RA2Rad = degrees(RA2Deg)
Dec2Rad = degrees(Dec2Deg)
RA3Rad = degrees(RA3Deg)
Dec3Rad = degrees(Dec3Deg)
RVec1 = np.array([-2.918014086392788E-2, 1.015971690494597, -4.441194241296249E-5])
RVec2 = np.array([-1.640387315073723E-1, 1.003339009095004, -3.884818173140164E-5])
RVec3 = np.array([-2.805159632942489E-1, 9.771715754745710E-1, -3.707484471326609E-5])
# Values to plug into CivToJulian function
Y = 2017
June = 6
July = 7
# Time values of observations 1, 2, and 3
t1 = CivToJulian(9, 43, 0, Y, June, 23)
t2 = CivToJulian(9, 56, 5.5, Y, July, 1)
t3 = CivToJulian(11, 25, 52, Y, July, 8)
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
    fT1 = 1 + (-mu) / (2 * (realRoot ** 3)) * tau1
    fT3 = 1 + (-mu) / (2 * (realRoot ** 3)) * tau3
    gT1 = tau1 + (-mu) / (6 * (realRoot ** 3)) * tau1
    gT3 = tau3 + (-mu) / (6 * (realRoot ** 3)) * tau3

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
    r2Mags = np.array([rootSolve(a, b, c)])
    EA1 = EAFunc(r2Mags, r2DotVec, tau1)
    EA3 = EAFunc(r2Mags, r2DotVec, tau3)
    f1, g1 = closedFuncs(a, r2Mags, EA1, tau1)
    f3, g3 = closedFuncs(a, r2Mags, EA3, tau3)
    i += 1


VLA0 = np.array([r2VecAr, r2Mags, r2DotVec, r2DotMags, rhoMag1, rhoMag2, rhoMag3])
VLA = itFunc(VLA0.item(0), VLA0.item(1), VLA0.item(2), VLA0.item(3), VLA0.item(4), VLA0.item(5), VLA0.item(6))
VLA = itFunc(VLA.item(0), VLA.item(1), VLA.item(2), VLA.item(3), VLA.item(4), VLA.item(5), VLA.item(6))
print "VLA2", VLA
print VLA.item(0)

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






############## CONVERT TO ECLIPTIC BEFORE STARTING BABY OD #########################







from __future__ import division
import numpy as np
from math import sqrt, sin, cos, acos, degrees, radians, atan2, pi

print("WELCOME TO MAXWELL'S OD PROGRAM! RECOMMENDED INPUTS ARE 2, 2, AND 2")

k = 0.01720209895
mu = 1

# Converts an eccentric vector into an equatorial vector
def EcToEq(EcVec):
    ob = radians(23.4347)
    EqVec = np.dot(np.array([[1, 0, 0], [0, cos(ob), sin(ob)], [0, -sin(ob), cos(ob)]]), EcVec)
    return EqVec

# Converts an eccentric vector into an equatorial vector
def EqToEc(EqVec):
    ob = radians(23.4347)
    EcVec = np.dot(np.array([[1, 0, 0], [0, cos(ob), sin(ob)], [0, -sin(ob), cos(ob)]]), EqVec)
    return EcVec


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

def readInput(data):
    year, month, day, time, timeH, timeM, timeS, RAH, RAM, RAS, DecDeg, DecM, DecS, RX, RY, RZ = [[], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []]

    for row in data:
        row = list(row)
        year.append(row[0])
        month.append(row[1])
        day.append(row[2])
        timeHMS = str(row[3])
        timeAr = timeHMS.split(":", 2)
        timeH = float(timeAr[0])
        timeM = float(timeAr[1])
        timeS = float(timeAr[2])
        timeNewAr = np.array([timeH, timeM, timeS])
        time = np.append(time, timeNewAr)
        RAH.append(row[4])
        RAM.append(row[5])
        RAS.append(row[6])
        DecDeg.append(row[7])
        DecM.append(row[8])
        DecS.append(row[9])
        RX.append(row[10])
        RY.append(row[11])
        RZ.append(row[12])

    # time1 = time[:3]
    # time2 = time[3:6]
    # time3 = time[6:9]
    # time4 = time[9:12]
    # finalTimes = np.array([time1, time2, time3, time4])  
    finalTimes = np.array([time[3*i:3*(i+1)] for i in range(4)])  
    # Dec1 = np.array([DecDeg[0], DecM[0], DecS[0]])
    # Dec2 = np.array([DecDeg[1], DecM[1], DecS[1]])
    # Dec3 = np.array([DecDeg[2], DecM[2], DecS[2]])
    # Dec4 = np.array([DecDeg[3], DecM[3], DecS[3]])
    RADegs = np.array([])
    DecDegs = np.array([])
    for i in range(4):
        RADegs = np.append(RADegs, RAHMStoDecDegrees(RAH[i], RAM[i], RAS[i]))
        DecDegs = np.append(DecDegs, DecDMStoDecDegrees(DecDeg[i], DecM[i], DecS[i]))
    # RARad1 = radians(RADegs[0])
    # RARad2 = radians(RADegs[1])
    # RARad3 = radians(RADegs[2])
    # RARad4 = radians(RADegs[3])
    RARadAr = np.array([radians(RADegs[i]) for i in range(4)])
    # DecRad1 = radians(DecDegs[0])
    # DecRad2 = radians(DecDegs[1])
    # DecRad3 = radians(DecDegs[2])
    # DecRad4 = radians(DecDegs[3])
    DecRadAr = np.array([radians(DecDegs[i]) for i in range(4)])
    return finalTimes, RARadAr, DecRadAr, year, month, day, RX, RY, RZ


# Julian days function
def CivToJulian(tUTHour, tUTMin, tUTSec, Y, M, D):
    UTDecimal = tUTHour + (tUTMin / 60) + (tUTSec / 3600)
    a = 367 * Y
    b = int((M + 9)/12)
    c = 7 * (Y + b)
    d = int(c / 4)
    e = int((275 * M) / 9)
    J0 = a - d + e + D + 1721013.5
    JD = J0 + UTDecimal / 24
    return JD

# function to get initial r2 vector value; only uses the first root
def rootSolve(a, b, c):
    coeffs = np.array([1, 0, a, 0, 0, b, 0, 0, c])
    roots = np.roots(coeffs)
    realRootsOld = roots[np.isreal(roots)]
    realRootsAll = np.array([])
    i = 0 # counting variable
    for element in realRootsOld:
        h = realRootsOld.item(i)
        realRootsAll= np.append(realRootsAll, abs(h))
        i += 1
    realRoots = np.unique(realRootsAll) # no more duplicates
    print("The real roots are", realRoots)
    prompt = "Which root would you like to choose; must be be one of the following numbers: "+ str(range(len(realRoots)))
    choice = input(prompt)
    realRoot = np.array([realRoots[choice]])
    return realRoot

# Function for finding rhoHat
def rhoHat(RARad, decRad):
    rhoHat = np.array([cos(RARad)*cos(decRad), sin(RARad) * cos(decRad), sin(decRad)])
    return rhoHat

# function for time correction
def tCorrection(JTimei, rhoMagi):
    cLight = 173.144633 # speed of light in AU per day
    tCorrection = JTimei - rhoMagi / cLight
    return tCorrection


# alternative to closed form functions; used in itFuncTay
def FnGTaylor(TAUi, r2Mag, r2Vec, r2DotVec):
    u = 1 / (r2Mag ** 3)
    z = np.dot(r2Vec, r2DotVec) / (r2Mag ** 2)
    q = (np.dot(r2DotVec, r2DotVec) / (r2Mag ** 2)) - u
    fLast = (1/24)*((3*u*q)-(15*u*(z**2))+(u**2))*(TAUi**4)
    f = 1-((1/2)*u*(TAUi**2)) + (1/2)*u*z*(TAUi**3) +fLast
    g = TAUi-(1/6)*u*(TAUi**3) + (1/4)*u*z*(TAUi**4)
    return f, g

def OrbElementz(xPos, yPos, zPos, xVel, yVel, zVel):

    answer = 2
    if answer == 2:
        t2 = CivToJulian(tH[1], tM[1], tS[1], Year[1], Month[1], Day[1])
        rhoHat2 = EqToEc(rhoHat(RA2Rad, Dec2Rad))
    elif answer == 3:
        t2 = CivToJulian(tH[2], tM[2], tS[2], Year[2], Month[2], Day[2])
        rhoHat2 = EqToEc(rhoHat(RA3Rad, Dec3Rad))



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

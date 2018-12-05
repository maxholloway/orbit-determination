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
        year.append(float(row[0]))
        month.append(float(row[1]))
        day.append(float(row[2]))
        timeHMS = str(row[3])
        timeAr = timeHMS.split(":", 2)
        timeH = float(timeAr[0])
        timeM = float(timeAr[1])
        timeS = float(timeAr[2])
        timeNewAr = [timeH, timeM, timeS]
        time += [timeNewAr]
        RAH.append(float(row[4]))
        RAM.append(float(row[5]))
        RAS.append(float(row[6]))
        DecDeg.append(float(row[7]))
        DecM.append(float(row[8]))
        DecS.append(float(row[9]))
        RX.append(float(row[10]))
        RY.append(float(row[11]))
        RZ.append(float(row[12]))


    finalTimes = time #[time[3*i:3*(i+1)] for i in range(4)]
    RADegs = np.array([])
    DecDegs = np.array([])
    for i in range(4):
        RADegs = np.append(RADegs, RAHMStoDecDegrees(RAH[i], RAM[i], RAS[i]))
        DecDegs = np.append(DecDegs, DecDMStoDecDegrees(DecDeg[i], DecM[i], DecS[i]))

    RARadAr = np.array([radians(RADegs[i]) for i in range(4)])
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
    prompt = "Which root would you like to choose; must be be one of the following numbers: "+ str([i for i in range(len(realRoots))])
    choice = int(input(prompt))
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

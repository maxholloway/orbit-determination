from __future__ import division
import numpy as np
from math import sqrt, sin, cos, acos, degrees, radians, atan2, pi

print "WELCOME TO MAXWELL'S OD PROGRAM! RECOMMENDED INPUTS ARE 2, 2, AND 2"

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
    year = []
    month = []
    day = []
    time = np.array([])
    timeH = []
    timeM = []
    timeS = []
    RAH = []
    RAM = []
    RAS = []
    DecDeg = []
    DecM = []
    DecS = []
    RX = []
    RY = []
    RZ = []
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

    time1 = time[:3]
    time2 = time[3:6]
    time3 = time[6:9]
    time4 = time[9:12]
    finalTimes = np.array([time1, time2, time3, time4])    
    Dec1 = np.array([DecDeg[0], DecM[0], DecS[0]])
    Dec2 = np.array([DecDeg[1], DecM[1], DecS[1]])
    Dec3 = np.array([DecDeg[2], DecM[2], DecS[2]])
    Dec4 = np.array([DecDeg[3], DecM[3], DecS[3]])
    RADegs = np.array([])
    DecDegs = np.array([])
    for i in range(4):
        RADegs = np.append(RADegs, RAHMStoDecDegrees(RAH[i], RAM[i], RAS[i]))
        DecDegs = np.append(DecDegs, DecDMStoDecDegrees(DecDeg[i], DecM[i], DecS[i]))
    RARad1 = radians(RADegs[0])
    RARad2 = radians(RADegs[1])
    RARad3 = radians(RADegs[2])
    RARad4 = radians(RADegs[3])
    RARadAr = np.array([RARad1, RARad2, RARad3, RARad4])
    DecRad1 = radians(DecDegs[0])
    DecRad2 = radians(DecDegs[1])
    DecRad3 = radians(DecDegs[2])
    DecRad4 = radians(DecDegs[3])
    DecRadAr = np.array([DecRad1, DecRad2, DecRad3, DecRad4])
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
    print "The real roots are", realRoots
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

observation = input("Observation 2 or observation 3 as center? ")
if observation == 2:
    t2 = CivToJulian(tH[1], tM[1], tS[1], Year[1], Month[1], Day[1])
    rhoHat2 = EqToEc(rhoHat(RA2Rad, Dec2Rad))
    RVec2 = RVec2 # for symmetry, I include this
elif observation == 3:
    t2 = CivToJulian(tH[2], tM[2], tS[2], Year[2], Month[2], Day[2])
    rhoHat2 = EqToEc(rhoHat(RA3Rad, Dec3Rad))
    RVec2 = RVec3
t1 = CivToJulian(tH[0], tM[0], tS[0], Year[0], Month[0], Day[0])
t3 = CivToJulian(tH[3], tM[3], tS[3], Year[3], Month[3], Day[3])
rhoHat1 = EqToEc(rhoHat(RA1Rad, Dec1Rad))
rhoHat3 = EqToEc(rhoHat(RA4Rad, Dec4Rad))


####--- SECTION GETTING IMPORTANT  VARIABLES TAKEN CARE OF ---####

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
RVec1 = np.array([RX[0], RY[0], RZ[0]])
RVec2 = np.array([RX[1], RY[1], RZ[1]])
RVec3 = np.array([RX[2], RY[2], RZ[2]])
RVec4 = np.array([RX[3], RY[3], RZ[3]])

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

####--- END OF SECTION GETTING IMPORTANT  VARIABLES TAKEN CARE OF ---####

# alternative to closed form functions; used in itFuncTay
def FnGTaylor(TAUi, r2Mags, r2Vec, r2DotVec):
    r2Mag = r2Mags
    u = 1 / (r2Mag ** 3)
    fDong = np.dot(r2Vec, r2DotVec)
    z = fDong / (r2Mag ** 2)
    fDing = np.dot(r2DotVec, r2DotVec)
    q = (fDing / (r2Mag ** 2)) - u
    fLast = (1/24)*((3*u*q)-(15*u*(z**2))+(u**2))*(TAUi**4)
    f = 1-((1/2)*u*(TAUi**2)) + (1/2)*u*z*(TAUi**3) +fLast
    g = TAUi-(1/6)*u*(TAUi**3) + (1/4)*u*z*(TAUi**4)
    return f, g

def itFuncTay(r2Vec, r2Mags, r2DotVec, r2DotMags, rhoMags, times):
    # do time correction
    t1, t2, t3 = times[0], times[1], times[2]
    rhoMag1, rhoMag2, rhoMag3 = rhoMags[0], rhoMags[1], rhoMags[2]
    tCorrected1 = tCorrection(t1, rhoMag1)
    tCorrected2 = tCorrection(t2, rhoMag2)
    tCorrected3 = tCorrection(t3, rhoMag3)
    times = np.array([tCorrected1, tCorrected2, tCorrected3])
    tau1 = k * (times[0] - times[1])
    tau3 = k * (times[2] - times[1])
    tau = k * (times[2] - times[0])
    f1, g1 = FnGTaylor(tau1, r2Mags, r2Vec, r2DotVec) # calling function to get f1 and g1 closed functions
    f3, g3 = FnGTaylor(tau3, r2Mags, r2Vec, r2DotVec) # calling function to get f3 and g3 closed functions
    
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
    r1Vec = rhoMag1 * rhoHat1 - RVec1
    r2Vec = rhoMag2 * rhoHat2 - RVec2
    r3Vec = rhoMag3 * rhoHat3 - RVec4

    # get d1 and d3 values
    d1 = (-f3) / (f1 * g3 - f3 * g1)
    d3 = f1 / (f1 * g3 - f3 * g1)
    
    # get r2 velocity vector and feel like a badass
    r2Mags = np.linalg.norm(r2Vec)
    r2DotVec = d1 * r1Vec + d3 * r3Vec
    r2DotMags = np.linalg.norm(r2DotVec)

    VLA = np.array([r2Vec, r2Mags, r2DotVec, r2DotMags, rhoMags, times])
    r2VecArEc = VLA[0]
    r2DotVecEc = VLA[2]
    xPos, yPos, zPos = r2VecArEc[0], r2VecArEc[1], r2VecArEc[2]
    xVel, yVel, zVel = r2DotVecEc[0], r2DotVecEc[1], r2DotVecEc[2]
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

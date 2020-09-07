import math
import numpy as np

# in km^3/s^2
mew = 398600.5

# in km
rEarth = 6378

h1 = 250
h2 = 200
h3 = 300

# in m/s
isp = 1960

# in years
dt = 2 

# lambda for circular velocity provided a given height (h) in km
vCircular = lambda x : math.sqrt( mew / ( rEarth + x ))

#circular velocity for h1
vcircular1 = vCircular(h1)
print 'vCircular for h1 : ' + str(h1) + " = " + str(vcircular1)

#circular velocity for h2
vcircular2 = vCircular(h2)
print 'vCircular for h2 : ' + str(h2) + " = " + str(vcircular2)

#circular velocity for h3
vcircular3 = vCircular(h3)
print 'vCircular for h3 : ' + str(h3) + " = " + str(vcircular3)

# Matrix A
A = [[1,2,3],[4,5,6],[9,7,6]]

# Matrix B
B = [[9,8,7],[6,5,4],[3,2,1]]

# dot product of matrices
dotProduct = lambda a, b : np.dot(a,b)

# cross product of matrices
crossProduct = lambda a, b : np.cross(a,b)

# Test case for dot product
print 'Here is the dot product of 2 Matrices \n' + str(dotProduct(A,B))

# Test case for cross product
print 'Here is the cross product of 2 Matrices \n' + str(crossProduct(A,B))

# Problem 4 - Page 107
# h1
F = 0.5 * 0.00000000006 * math.pow(vcircular1*1000, 2) *2.2 * 1

print 'Force for h1 - ' + str(F)

s = 63072000

mf = (F * s) / isp

print 'The mass fuel needed for ' + str(h1) + ' is equal to ' + str(mf)
# h1 end

# h2
F = 0.5 * 0.00000000006 * math.pow(vcircular2*1000, 2) *2.2 * 1

print 'Force for h1 - ' + str(F)

s = 63072000

mf = (F * s) / isp

print 'The mass fuel needed for ' + str(h2) + ' is equal to ' + str(mf)
# h2 end

# h3
F = 0.5 * 0.00000000006 * math.pow(vcircular3*1000, 2) *2.2 * 1

print 'Force for h1 - ' + str(F)

s = 63072000

mf = (F * s) / isp

print 'The mass fuel needed for ' + str(h3) + ' is equal to ' + str(mf)
# h3 end

# Calculate the magnitude of a vectors components in 3 space
magnitude = lambda x, y, z : math.sqrt(math.pow(x,2) + math.pow(y, 2) + math.pow(z, 2))

# Test case
print 'The magnitude of a velovity vector with components x = 2, y = 4, z = 7 :\n'
print str(magnitude(2,4,7))

# Calculate the magnitude of a vectors components in 2 space
magnitude = lambda x, y : math.sqrt(math.pow(x,2) + math.pow(y, 2))

# Test case
print 'The magnitude of a velovity vector with components x = 2, y = 4:\n'
print str(magnitude(2,4))

# Convert spherical (Geodedic) coordinates to rectangular
# x = rcos(theta)cos(lambda)
# r = radius of Earth + height above earth (elevation)
elevation1 = 300
r = rEarth + elevation1
thet = 45
lam = 120
xrect = r * math.cos(thet) * math.cos(lam)
print 'The rectangular x coordinate is : ' + str(xrect)

yrect = r * math.cos(thet) * math.sin(lam)
print 'The rectangular y coordinate is : ' + str(yrect)

zrect = r * math.sin(thet)
print 'The rectangular z coordinate is : ' + str(zrect)

convertRectangular = lambda latitude, longitude, height : 
	xrect = r * math.cos(thet) * math.cos(lam)
	yrect = r * math.cos(thet) * math.sin(lam)
	zrect = r * math.sin(thet)
	print 'x = ' + str(xrect)
	print 'y = ' + str(yrect)
	print 'z = ' + str(zrect)

def geoToRect(latitude, longitude, height) :
	xrect = r * math.cos(thet) * math.cos(lam)
	yrect = r * math.cos(thet) * math.sin(lam)
	zrect = r * math.sin(thet)
	return [xrect, yrect, zrect]

# Test case for convertRectangular
print 'Converting lat = 45 degrees'
print 'longitude = 120 degrees'
print 'Elevation = 300km'
print 'Radius of Earth = 6378km'
convertRectangular(thet,lam,elevation1)

# Calculate ECEF coordinates from Geodedic 
# Given
# Lat = 28.5 degrees
# long = -82.0 degrees
# altitude above the elipsoid = 150 ft
# Ellipticity (f) = 0.0033523
# Earths equitorial radius in ft = 20925741.57ft
latitude = 28.5
longitude = -82.0
height = 150
ellipticity = 0.0033523
eradiusft = 20925741.57

xef = ((eradiusft/math.sqrt(math.pow(math.cos(latitude),2) + math.pow((1 - ellipticity), 2) * math.pow(math.sin(latitude), 2))) + height) * math.cos(latitude) * math.cos(longitude)
print 'xECEF = ' + str(xef) + ' ft'
print 'Expected xECEF = 2561350.138 ft'

yef = ((eradiusft/math.sqrt(math.pow(math.cos(latitude),2) + math.pow((1 - ellipticity), 2) * math.pow(math.sin(latitude), 2))) + height) * math.cos(latitude) * math.sin(longitude)
print 'yECEF = ' + str(yef) + ' ft'
print 'Expected yECEF = -18224953.22 ft'


zef = (((math.pow((1 - ellipticity), 2) * eradiusft)/math.sqrt(math.pow(math.cos(latitude),2) + math.pow((1 - ellipticity), 2) * math.pow(math.sin(latitude), 2))) + height) * math.sin(latitude)
print 'zECEF = ' + str(zef) + ' ft'
print 'Expected zECEF = 9925705.88 ft'

convertECEF = lambda latitude, longitude, height, ellipticity, eradiusft :
	xef = ((eradiusft/math.sqrt(math.pow(math.cos(latitude),2) + math.pow((1 - ellipticity), 2) * math.pow(math.sin(latitude), 2))) + height) * math.cos(latitude) * math.cos(longitude)
	yef = ((eradiusft/math.sqrt(math.pow(math.cos(latitude),2) + math.pow((1 - ellipticity), 2) * math.pow(math.sin(latitude), 2))) + height) * math.cos(latitude) * math.sin(longitude)
	zef = (((math.pow((1 - ellipticity), 2) * eradiusft)/math.sqrt(math.pow(math.cos(latitude),2) + math.pow((1 - ellipticity), 2) * math.pow(math.sin(latitude), 2))) + height) * math.sin(latitude)
	print 'xECEF = ' + str(xef) + ' ft'
	print 'yECEF = ' + str(yef) + ' ft'
	print 'zECEF = ' + str(zef) + ' ft'


# CONVERT TO RADIANS, DO NOT USE DEGREES
def convertECEFfeet(latitude, longitude, height, ellipticity, eradiusft) :
	xef = ((eradiusft/math.sqrt(math.pow(math.cos(latitude),2) + math.pow((1 - ellipticity), 2) * math.pow(math.sin(latitude), 2))) + height) * math.cos(latitude) * math.cos(longitude)
	yef = ((eradiusft/math.sqrt(math.pow(math.cos(latitude),2) + math.pow((1 - ellipticity), 2) * math.pow(math.sin(latitude), 2))) + height) * math.cos(latitude) * math.sin(longitude)
	zef = (((math.pow((1 - ellipticity), 2) * eradiusft)/math.sqrt(math.pow(math.cos(latitude),2) + math.pow((1 - ellipticity), 2) * math.pow(math.sin(latitude), 2))) + height) * math.sin(latitude)
	return [xef, yef, zef]

# Test case for convertECEF
print 'Testing converting to ECEF coordinates from Geodedic:'
convertECEF(latitude, longitude, height, ellipticity, eradiusft)

#Linear Tangent Law

def linearTangentLawDegrees(timeInterval, timeEnd, initThrustAngle) :
	return math.degrees(math.atan((1 - (timeInterval/timeEnd)) * math.tan(math.radians(initThrustAngle))))

def linearTangentLawRadians(timeInterval, timeEnd, initThrustAngle) :
	return math.atan((1 - (timeInterval/timeEnd)) * math.tan(math.radians(initThrustAngle)))

def linearTangentLawTan(timeInterval, timeEnd, initT hrustAngle) :
	return (1 - (timeInterval/timeEnd)) * math.tan(math.radians(initThrustAngle))

def peg(rnow, vnow, t, rd, tgo, rgrav, S, Q, vgo, magVgo) :
    print('Generating PEG output...')
    print('----------------------------------------------------------------------')
    # a) Refernce thrust vector of LambdaV
    lamv = np.divide(vgo, magVgo)
    print('lamda V          = ', lamv)

    # b) Compute the reference time tlambda
    L = magVgo
    J = L * tgo - S
    tlam = J / L
    print('t lambda         = ', tlam)

    # c) Compute reference thrust turning rate lambda prime
    rgo = np.subtract(np.subtract(np.subtract(rd, rnow) ,rgrav), np.multiply(vnow, tgo))

    lamPrime = (np.subtract(rgo, np.multiply(S, lamv)))/(Q - (S*tlam))

    print('Lambda prime     = ', lamPrime)

    # d) Calculate the linear Tangent thrust vector and the unit vector
    lamf = np.add(lamv, np.multiply(lamPrime, (t - (tlam + t))))
    uf = np.divide(lamf, (np.sqrt((lamf[0] * lamf[0]) + (lamf[1] * lamf[1]) + (lamf[2] * lamf[2]))))
    print('Lambda f         = ', lamf)
    print('u sub f          = ', uf)

    # e) To Do - Question 26

    # f) vmeco, flight path anlge calculation
    a = (r1+r2)/2
    vmeco = math.sqrt(mu *((2/rmeco) - (1/a)))
    vp = math.sqrt(mu * ((2/r1) - (1/a)))
    magh = np.dot(r1, vp)
    gamma = math.degrees(math.acos(magh/(np.dot(vmeco, rmeco))))
    print('Velocity at MECO = ', vmeco)
    print('Gamma            = ', gamma)
    
    print('----------------------------------------------------------------------')


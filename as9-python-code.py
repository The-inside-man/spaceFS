import numpy as np
import math
import sys


# Here we open the file as append so that we may continuously
# write to the file throughout this program
f = open('Brown_Assignment9_SPC5010GNC_10-21-19.txt', 'a')

# Lets print the header of the homework answer sheet
f.write('\nBrown, Jacob - Assignment 9 - SPC5010 GNC - Calculation Result Sheet')
f.write('\n-----------------------------------------------------------------')

################################
# Question 1
################################
f.write('\n\nQuestion 1:')
f.write('\n\tAnswer: L = 0')

################################
# Question 2
################################
f.write('\n\nQuestion 2:')
f.write('\n\tAnswer: Ballistic and gluiding entry.')

################################
# Question 3
################################
f.write('\n\nQuestion 3:')
f.write('\n\tAnswer: \tGravity Turn Descent Guidance.')
f.write('\n\t\tTerminal State Vector Control.')
f.write('\n\t\tPowered Explicit Guidance.')
f.write('\n\t\tApollo Lunar Powered Descent Guidance.')

################################
# Question 4
################################
f.write('\n\nQuestion 4:')

#Knowns
g = 3.7
h = 1500
v = 10
ht = 10
vt = 1.0


# Compute aT
a = 1
b = -(2 * ((math.pow(v,2) + math.pow(vt, 2)) / (4 * g * (h - ht))))
print('b = ' + str(b))

c = -(1 + (((2 * math.pow(v, 2)) - (2 * math.pow(vt, 2))) / (4 * g * (h - ht))))
print('c = ' + str(c))

quadPos = (-b + math.sqrt(math.pow(b, 2) - 4 * a * c)) / 2 * a
print('Quadratic pos = ' + str(quadPos))

quadNeg = (-b - math.sqrt(math.pow(b, 2) - 4 * a * c)) / 2 * a
print('Quadratic Neg = ' + str(quadNeg))

aT = quadPos * g
print('aT = ' + str(aT))

# Solution
f.write('\n\tAnswer: \t' + str(aT) + ' m/s^2')

################################
# Question 5
################################
f.write('\n\nQuestion 5:')

# Knowns
g = 3.7
r = 1500
v = 10
rf = 10
vf = 1.0
tgo = 1550

# Compute aT

# let x = first part of equation
x = (vf - v) / tgo

# let y = second portion
y = (rf - (r + vf * tgo)) / math.pow(tgo, 2)

aT = 4 * x + 6 * y - g

print('aT = ' + str(aT))

# Solution
f.write('\n\tAnswer: \t' + str(aT) + ' m/s^2')

################################
# Question 6
################################
f.write('\n\nQuestion 6:')
f.write('\n\tAnswer: \t1. Instrument an inertial reference frame using gyros')
f.write('\n\t\t2. Measure specific force using accelerometers in the inertial ref frame')
f.write('\n\t\t3. Have a model of the gravitational field and be able to compute the gravitational acceleration')
f.write('\n\t\t4. Integrate the equations of motion and compute the position and velocity')

################################
# Question 7
################################
f.write('\n\nQuestion 7:')
f.write('\n\tAnswer: \tThe gradient of a potential function specifies the direction in 3D space \n\t\tto move to increase our function. The gradient at any location \n\t\tpoints in the direction of greatest increase of a function.')

################################
# Question 8
################################
f.write('\n\nQuestion 8:')
f.write('\n\tAnswer: latitudinal bands around a planet.')

################################
# Question 9
################################
f.write('\n\nQuestion 9:')
f.write('\n\tAnswer:  A planet’s equatorial bulge')

################################
# Question 10
################################
f.write('\n\nQuestion 10:')
f.write('\n\tAnswer:  No, there is no sensed acceleration as the spacecraft and accelerometer are in free fall.')

################################
# Question 11
################################
f.write('\n\nQuestion 11:')
f.write('\n\tAnswer:  Perform alignment/gyrocompass to obtain the orientation of the ins frame to NED frame.')

################################
# Question 12
################################
f.write('\n\nQuestion 12:')

# Knowns
# Lat North
lat = 29.5
# degrees/hour
Re = 15

ReNorth = Re * math.cos(math.radians(lat))
print(str(ReNorth))

ReUp = Re * math.sin(math.radians(lat))
print(str(ReUp))

f.write('\n\tAnswer: Earth Rate North : ' + str(ReNorth))
f.write('\n\tAnswer: Earth Rate up : ' + str(ReUp))

################################
# Question 13
################################
f.write('\n\nQuestion 13:')
f.write('\n\tAnswer:  \tThe INS is now performing all the functions to determine position and \n\t\tvelocity and attitude in inertial space.')

################################
# Question 14
################################
f.write('\n\nQuestion 14:')
f.write('\n\tAnswer:  mathematically and maintains it using gyros measuring inertial angular rates')

################################
# Question 15
################################
f.write('\n\nQuestion 15:')

# Knowns
g = 9.806
fz = 9.80599
fx = 0.0017114
 
# Solution
tilt = math.degrees(math.atan(math.radians(fx)/math.radians(fz)))
print(str(tilt))


f.write('\n\tAnswer: ' + str(tilt) + ' degrees or 36 ARCS')

################################
# Question 16
################################
f.write('\n\nQuestion 16:')
f.write('\n\tAnswer:  Attitude Reference Unit')

################################
# Question 17
################################
f.write('\n\nQuestion 17:')
f.write('\n\tAnswer:  Inertial Measurement Unit')

################################
# Question 18
################################
f.write('\n\nQuestion 18:')
f.write('\n\tAnswer:  Inertial Navigation System')

################################
# Question 19
################################
f.write('\n\nQuestion 19:')

# Knowns

# in mm
rcoil = 40
# in m
lcoil = 100
# in m
lamdaLaser = 632.8 * math.pow(10, -9)
# Accelerometer Mass in kg
mass = 0.005
# Spring Constant in Nm
Ks = 150

insFrameToBody = [[0.5944400, 0.8041300, -0.0025200], [-0.8041300, 0.5944400, 0.0034700], [0.0042900, -0.0000400, 0.9999990]]

# Phase diff (deg)
x = 0.0025
y = 0.03
z = 0.20010

# Mass deflection
dx = 0.00051
dy = 0.00011
dz = 0.00057

vehicleBodyFrameToJ2000 = [[0.5000000, 0.8660254, 0.0000000], [-0.8660254, 0.5000000, 0.0000000], [0.0000000, 0.0000000, 1.0000000]]

magnitude = lambda x, y, z : math.sqrt(math.pow(x,2) + math.pow(y, 2) + math.pow(z, 2))

def ifog(deltaFi, waveL, R, L) :
	c = 3.0 * math.pow(10, 8)
	return (deltaFi * c * waveL/(4 * math.pi * R * L))

# Step 1 - IFOG equeation
rotx = ifog(x, lamdaLaser, lcoil, (rcoil/1000))
print(str(rotx))
roty = ifog(y, lamdaLaser, lcoil, (rcoil/1000))
print(str(roty))
rotz = ifog(z, lamdaLaser, lcoil, (rcoil/1000))
print(str(rotz))
mag = magnitude(rotx, roty, rotz)
print(str(mag))

# Step 2
rotVect = [rotx, roty, rotz]

vehicleBodyFrame = np.dot(insFrameToBody, rotVect)
print(str(vehicleBodyFrame))
mag = magnitude(vehicleBodyFrame[0], vehicleBodyFrame[1], vehicleBodyFrame[2])
print(str(mag))

# Step 3 - Convert to acceleration
varx = ((Ks / mass) * dx) 
print(str(varx))
vary = ((Ks / mass) * dy) 
print(str(vary))
varz = ((Ks / mass) * dz) 
print(str(varz))

accVect = [varx, vary, varz]

accBodyFrame = np.dot(insFrameToBody, accVect)
print(str(accBodyFrame))

# Step 4 - Convert J2000

bodyRateJ2000 = np.dot(vehicleBodyFrameToJ2000, vehicleBodyFrame)
bodyRateMag = magnitude(bodyRateJ2000[0], bodyRateJ2000[1], bodyRateJ2000[2])
print(str(bodyRateJ2000))

accRateJ2000 = np.dot(vehicleBodyFrameToJ2000, accBodyFrame)
accRateMag = magnitude(accRateJ2000[0], accRateJ2000[1], accRateJ2000[2])
print(str(accRateJ2000))

f.write('\n\tAnswer:')
f.write('\n\tBody Rate J2000 :  '  + str(bodyRateJ2000))
f.write('\n\tMagnitude       :  '  + str(bodyRateMag))
f.write('\n\tBody Rate J2000 :  '  + str(accRateJ2000))
f.write('\n\tMagnitude       :  '  + str(accRateMag))

################################
# Question 20
################################
f.write('\n\nQuestion 20:')
f.write('\n\tAnswer:  Alignment')

################################
# Question 21
################################
f.write('\n\nQuestion 21:')
f.write('\n\tAnswer:  MEMS')
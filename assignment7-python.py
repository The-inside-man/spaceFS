import math
#import numpy as np

# Here is the file that will contain the answers
# to the calculation questions
f = open('brown_assignment7_results.txt', 'a')

# Print the header to the file
f.write('Brown, Jacob - Assignment 7 - Calculations\n')
f.write('-----------------------------------------------------')

#####################################
# Question 3
#####################################
print('-----------------------------------------------------')
print('Question 3:')
# Knowns
di = 5
rearth = 6378
h = 200
r = rearth + h
mu = 398600

# Calculate v
def findV(mu, radius) :
    return math.sqrt(mu/radius)

v = findV(mu, r)
# Print v to the screen for debugging
print('v = ', v)

# Calculate dv - circular
def deltaVCirc(v, di) :
    return (2 * v)* math.sin(math.radians(di)/2)

dv = deltaVCirc(v, di)
# Print dv to screen for debugging
# Note we will multiply the answer by 1000 to switch to meters
print('dv = ', dv * 1000)

# Write the answer to the problem in file
f.write('\nQuestion 3:\n')
f.write('dv = ' + str(dv * 1000) + "\n")
print('-----------------------------------------------------')
#####################################
# Question 4
#####################################
print('Question 4:')
# Knowns
rearth = 6378
h = 200
r = rearth + h
di = 5
mu = 398600
inci = 5
incf = 10
raani = 0
raanf = 25
dRaan = 25

# First we need to calculate alpha
def findAlpha(inci, incf, dRaan) :
    return math.degrees(math.acos(((math.cos(inci) * math.cos(incf)) + (math.sin(inci) * math.sin(incf) * math.cos(dRaan)))))

a = findAlpha(inci, incf, dRaan)
# Print to the screen for debugging
print('Alpha = ', a)

# Calculate dv - circular
def deltaVRAANCirc(v, a) :
    return (2 * v) * math.sin(a/2)

dv = deltaVRAANCirc(v, a)
# Write the answer to the problem in file
f.write('Questions 4:\n')
f.write('dv = ' + str(dv) + "\n")

# Print to screen for debugging
print('dv = ', dv)
print('-----------------------------------------------------')
#####################################
# Question 5
#####################################
print('Question 5:')
ra = rearth + 300
rp = rearth + 200

# Calculate Vf and print both vi and vf
vf = math.sqrt(mu * ((2/rp) - (1/ra)))
print('vf = ', vf)
vi = v
print('vi = ', vi)

def combinedIncandAlt(vi, vf, di) :
    return math.sqrt((math.pow(vi, 2) + math.pow(vf, 2)) - (2 * (vi * vf * (math.cos(di)))))

dv = combinedIncandAlt(vi, vf, di)

# Print answer for debugging
print('dv = ', dv)

# Write answer to problem in file
f.write('Question 5:\n')
f.write('dv = ' + str(dv) + "\n")

print('-----------------------------------------------------')
#####################################
# Question 6
#####################################
print('Question 6:')
# Below are the CW equations

# Note n and nt is in radians/sec and radians
def solveVy0(x0, t, y0, n) :
    num = (((6 * x0) * (n*t - math.sin(n*t)) - y0) * n * math.sin(n*t) - (n * x0 * (4-3*math.cos(n*t))*(1 - math.cos(n*t))))
    den = (4 * math.sin(n*t) - (3 * n*t)) * math.sin(n*t) + 4 * math.pow((1 - math.cos(n*t)), 2)
    return num/den

# Note n and nt is in radians/sec and radians
def solveVx0(vy0, x0, y0, t, n):
    num = -(n * x0) * (4 - 3 * math.cos(n*t)) - 2 * ((1 - math.cos(n*t)) * vy0)
    den = math.sin(n*t)
    return num/den

# Note n and nt is in radians/sec and radians
def solveVxf(n, t, vx0, x0, vy0):
    return (vx0 * math.cos(n*t) + ((2 * vy0 + 3 * n * x0) * math.sin(n*t)))

# Note n and nt is in radians/sec and radians
def solveVyf(vx0, vy0, n, t, x0):
    return (((-2 * vx0) * math.sin(n*t)) + ((4 * vy0 + 6 * n * x0) * math.cos(n*t)) - ((3 * vy0 + 6 * n * x0)))

# Solve for a
a = (357 + 404.5 + 2 * (rearth))/2

# Solve for n
n = math.sqrt(mu/math.pow(a, 3))

#solve for t in seconds
t = 36 * 60

# Set x0 and y0 in meters
x0 = -2.5 * 1000
y0 = -5 * 1000

vy0 = solveVy0(x0, math.radians(t), y0, math.radians(n))
print('vy0 = ', vy0)
# Write the answer to file
f.write('Question 6:\n')
f.write('vy0 = ' +  str(vy0) + "\n")

vx0 = solveVx0(vy0, x0, y0, math.radians(t), math.radians(n))
print('vx0 = ', vx0)
# Write the answer to file
f.write('vx0 = ' +  str(vx0) + "\n")

vxf = solveVxf(math.radians(n), math.radians(t), vx0, x0, vy0)
print('vxf = ', vxf)
# Write the answer to file
f.write('vxf = ' + str(vxf) + "\n")

vyf = solveVyf(vx0, vy0, math.radians(n), math.radians(t), x0)
print('vyf = ', vyf)
# Write the answer to file
f.write('vyf = ' + str(vyf) + "\n")


# Close the output file
f.close()


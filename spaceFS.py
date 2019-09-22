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

# dot product of matrices
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

# h2
F = 0.5 * 0.00000000006 * math.pow(vcircular2*1000, 2) *2.2 * 1

print 'Force for h1 - ' + str(F)

s = 63072000

mf = (F * s) / isp

print 'The mass fuel needed for ' + str(h2) + ' is equal to ' + str(mf)

#h3
F = 0.5 * 0.00000000006 * math.pow(vcircular3*1000, 2) *2.2 * 1

print 'Force for h1 - ' + str(F)

s = 63072000

mf = (F * s) / isp

print 'The mass fuel needed for ' + str(h3) + ' is equal to ' + str(mf)



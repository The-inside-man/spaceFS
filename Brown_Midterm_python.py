import math
import numpy as np

# Here we open the file as append so that we may continuously
# write to the file throughout this midterm.
f = open('Brown_Midterm-SPC5010GNC_10-21-19.txt', 'a')

# Lets print the header of the Midterm Answer sheet
f.write('Brown, Jacob - Midterm for SPC5010 GNC - Calculation Result Sheet\n')
f.write('-----------------------------------------------------------------')

#####################################
# Question 1
#####################################
f.write('\nQuestion 1:\n')

#---------------
# Knowns
#---------------
# Degrees
lat = 40
# Degrees
long = 100
# Provided Value
ellip = 0.0012
# km
height = 0.03
# km
radius = 1737.4

#---------------
# Solution
#---------------

def convertMCMFkm(latitude, longitude, height, ellipticity, radiuskm) :
	xmf = ((radiuskm/math.sqrt(math.pow(math.cos(math.radians(latitude)),2) + math.pow((1 - ellipticity), 2) * math.pow(math.sin(math.radians(latitude)), 2))) + height) * math.cos(math.radians(latitude)) * math.cos(math.radians(longitude))
	ymf = ((radiuskm/math.sqrt(math.pow(math.cos(math.radians(latitude)),2) + math.pow((1 - ellipticity), 2) * math.pow(math.sin(math.radians(latitude)), 2))) + height) * math.cos(math.radians(latitude)) * math.sin(math.radians(longitude))
	zmf = (((math.pow((1 - ellipticity), 2) * radiuskm)/math.sqrt(math.pow(math.cos(math.radians(latitude)),2) + math.pow((1 - ellipticity), 2) * math.pow(math.sin(math.radians(latitude)), 2))) + height) * math.sin(math.radians(latitude))
	return [xmf, ymf, zmf]

ans = convertMCMFkm(lat, long, height, ellip, radius)
print('Answer 1: ' + str(ans))

# Answer is returned in the variable ans printed above
# We will now add to our answer sheet
f.write('Answer: ' + str(ans))
f.write('-----------------------------------------------------------------')

#####################################
# Question 2
#####################################
f.write('\nQuestion 2:\n')

#---------------
# Knowns
#---------------

mcmfTJ2000 = np.array([[0.9636305, 0.2672384, 0.0000000],[-0.2672384, 0.9636305, 0.0000000],[0.0000000, 0.0000000, 1.0000000]])
# Transpose the RNP matrix
mcmfTJ2000 = mcmfTJ2000.transpose()
# Convert to a NUMPY array to allow proper arithmatic 
mcmf = np.array(ans)

#---------------
# Solution
#---------------

# Compute the dot product 
j2000 = mcmfTJ2000.dot(mcmf)

# Compute the Magnitude
magnitude = lambda x, y, z : math.sqrt(math.pow(x,2) + math.pow(y, 2) + math.pow(z, 2))
mag = magnitude(j2000[0], j2000[1], j2000[2])

# Display for debugging
print('Answer 2: ' + str(j2000))
print('Answer 2 Magnitude: ' + str(mag))

# Write the answer to the solution sheet
f.write('Answer: ' + str(j2000))
f.write('Magnitude: ' + str(mag))
f.write('-----------------------------------------------------------------')

#####################################
# Question 3
#####################################
f.write('\nQuestion 3:\n')

#---------------
# Knowns
#---------------

# Radians per second
rate = 0.0000026705589
wmcmf = ([0,0,rate])

#---------------
# Solution
#---------------

vj2000 = np.cross(wmcmf, ans)

vj2000 = mcmfTJ2000.dot(vj2000)
# Display the answer for debugging
mag = magnitude(vj2000[0], vj2000[1], vj2000[2])
print('Answer 3: ' + str(vj2000))
print('Answer 3 Magnitude: ' + str(mag))

# Write the answer to the solution sheet
f.write('Answer: ' + str(vj2000))
f.write('Magnitude: ' + str(mag))
f.write('-----------------------------------------------------------------')

#####################################
# Question 4
#####################################
f.write('\nQuestion 4:\n')

#---------------
# Knowns
#---------------

# Seconds
t0 = 0
t1 = 150
t2 = 300
t3 = 350
t4 = 450
t5 = 480

# Total time in seconds
T = 480

# Degrees
y0 = 45

#---------------
# Solution
#---------------

def linearTangentLawDegrees(timeInterval, timeEnd, initThrustAngle) :
    return math.degrees(math.atan((1 - (timeInterval/timeEnd)) * math.tan(math.radians(initThrustAngle))))

fi0 = linearTangentLawDegrees(t0, T, 45)
fi1 = linearTangentLawDegrees(t1, T, 45)
fi2 = linearTangentLawDegrees(t2, T, 45)
fi3 = linearTangentLawDegrees(t3, T, 45)
fi4 = linearTangentLawDegrees(t4, T, 45)
fi5 = linearTangentLawDegrees(t5, T, 45)

print('Answer 4: t0 = ' + str(fi0))
print('Answer 4: t1 = ' + str(fi1))
print('Answer 4: t2 = ' + str(fi2))
print('Answer 4: t3 = ' + str(fi3))
print('Answer 4: t4 = ' + str(fi4))
print('Answer 4: t5 = ' + str(fi5))


# Write the answer to the solution sheet
f.write('Answer 4: t0 = ' + str(fi0))
f.write('Answer 4: t1 = ' + str(fi1))
f.write('Answer 4: t2 = ' + str(fi2))
f.write('Answer 4: t3 = ' + str(fi3))
f.write('Answer 4: t4 = ' + str(fi4))
f.write('Answer 4: t5 = ' + str(fi5))
f.write('-----------------------------------------------------------------')

#####################################
# Question 5
#####################################
f.write('\nQuestion 5:\n')

#---------------
# Knowns
#---------------

# In km's
r1 = 10
r2 = 60
rdesired = 15
rmoon = 1737.4
rp = rmoon + r1
ra = rmoon + r2
rmeco = rmoon + rdesired
# In km3/sec2
mu = 4902.844

#---------------
# Solution
#---------------

a = ((rp+ra)/2)
vmeco = math.sqrt(mu * ((2/rmeco) - (1/a)) )

print('Asnwer 5: ' + str(vmeco))

# Write the answer to the solution sheet
f.write('Answer: ' + str(vmeco))
f.write('-----------------------------------------------------------------')

#####################################
# Question 6
#####################################
f.write('\nQuestion 6:\n')

#---------------
# Knowns
#---------------

# In km's
r1 = 10
r2 = 60
rdesired = 15
rmoon = 1737.4
rmeco = rmoon + rdesired
rp = rmoon + r1
ra = rmoon + r2
# In km3/sec2
mu = 4902.844


#---------------
# Solution
#---------------
vp = math.sqrt(mu * ((2/rp) - (1/a)))
print('Answer 6: vp = ' + str(vp))
magh = rp * vp
print('Answer 6: |h| = ' + str(magh))
gamma = math.degrees(math.acos(magh/(vmeco * rmeco)))
print('Answer 6: Gamma = ' + str(gamma))

# Write the answer to the solution sheet
f.write('Answer: ' + str(gamma))
f.write('-----------------------------------------------------------------')

#####################################
# Question 7
#####################################
f.write('\nQuestion 7:\n')
f.write('Answer : Super G')
f.write('-----------------------------------------------------------------')

#####################################
# Question 8
#####################################
f.write('\nQuestion 8:\n')
f.write('Answer : Cowells Method')
f.write('-----------------------------------------------------------------')

#####################################
# Question 9
#####################################
f.write('\nQuestion 9:\n')
f.write('Answer : the v infinity vector at the departure planet.')     
f.write('-----------------------------------------------------------------')

#####################################
# Question 10
#####################################
f.write('\nQuestion 10:\n')

#---------------
# Knowns
#---------------

# In km/sec
v1 = ([7.5, 25.6,  1.0])
v2 = ([-10.0, -21.6, 0.1])

# In km/sec
vp1 = ([5.5,  10.0,  0.5])
vp2 = ([-5.0,  -15.0, 0.0]) 

#---------------
# Solution
#---------------

v1inf = np.subtract(v1, vp1)
print('Question 10:' + str(v1inf))

# Write the answer to the solution sheet
f.write('Answer: ' + str(v1inf))
f.write('-----------------------------------------------------------------')

#####################################
# Question 11
#####################################
f.write('\nQuestion 11:\n')

#---------------
# Knowns
#---------------

# In km/sec
v1 = ([7.5, 25.6,  1.0])
v2 = ([-10.0, -21.6, 0.1])

# In km/sec
vp1 = ([5.5,  10.0,  0.5])
vp2 = ([-5.0,  -15.0, 0.0])

#---------------
# Solution
#---------------

v2inf = np.subtract(v2, vp2)
print('Question 11:' + str(v2inf))

# Write the answer to the solution sheet
f.write('Answer: ' + str(v2inf))
f.write('-----------------------------------------------------------------')

#####################################
# Question 12
#####################################
f.write('\nQuestion 12:\n')
f.write('Answer : Hyperbolic')     
f.write('-----------------------------------------------------------------')

#####################################
# Question 13
#####################################
f.write('\nQuestion 13:\n')
f.write('Answer : Vinjection velocity is greater (fig 7.11)')       
f.write('-----------------------------------------------------------------')

#####################################
# Question 14
#####################################
f.write('\nQuestion 14:\n')
f.write('Answer : the semi-major axis of the transfer orbit to achieve the desired flight time')     
f.write('-----------------------------------------------------------------')

#####################################
# Question 15
#####################################
f.write('\nQuestion 15:\n')
f.write('Answer : NSR Maneuver')     
f.write('-----------------------------------------------------------------')

#####################################
# Question 16
#####################################
f.write('\nQuestion 16:\n')
f.write('Answer : NC Maneuver')     
f.write('-----------------------------------------------------------------')

#####################################
# Question 17
#####################################
f.write('\nQuestion 17:\n')
f.write('Answer : Positive Down, Negative Up to the target vehicle.')     
f.write('-----------------------------------------------------------------')

#####################################
# Question 18
#####################################
f.write('\nQuestion 18:\n')

#---------------
# Knowns
#---------------

# Circular Oribit in km
alt = 407
# Chaser distance from ISS at same alt
yo = -5
# no offset in altitude
xo = -2.5
# Minutes till rendevous converted to seconds
t = 36 * 60
# mu is in km^3/sec^2
mu = 398600
# Radius of Earth
r_e = 6378

#---------------
# Solution
#---------------

# In km
a = ((2 * alt) + (2 * r_e))/2
print('Question 18: ' + str(a))

a1 = 6783.75
a2 = 407 + 6378

# In radians per second
n1 = math.sqrt(mu/math.pow(a1, 3))

n2 = math.sqrt(mu/math.pow(a2, 3))
print('Question 18: ' + str(n1))

# In radians
s = np.sin(n1*t)
print('Question 18: ' + str(s))

# In radians
c = np.cos(n1*t)
print('Question 18: ' + str(c))

# Note n and nt is in radians/sec and radians
def solveVy0(x0, t, y0, n) :
    num = (((6 * x0) * (n*t - math.sin(n*t)) - y0) * n * math.sin(n*t) - (n * x0 * (4-3*math.cos(n*t))*(1 - math.cos(n*t))))
    den = (4 * math.sin(n*t) - (3 * n*t)) * math.sin(n*t) + 4 * math.pow((1 - math.cos(n*t)), 2)
    return num/den

vy0 = (solveVy0(xo, math.radians(t), yo, math.radians(n1)))

# Note n and nt is in radians/sec and radians
def solveVx0(vy0, x0, y0, t, n):
    num = -(n * x0) * (4 - 3 * math.cos(n*t)) - 2 * ((1 - math.cos(n*t)) * vy0)
    den = math.sin(n*t)
    return num/den

vx0 = (solveVx0(vy0, xo, yo, math.radians(t), math.radians(n1)))

# Note n and nt is in radians/sec and radians
def solveVxf(n, t, vx0, x0, vy0):
    return (vx0 * math.cos(n*t) + ((2 * vy0 + 3 * n * x0) * math.sin(n*t)))

vxf = (solveVxf(math.radians(n2), math.radians(t), vx0, xo, vy0))

# Note n and nt is in radians/sec and radians
def solveVyf(vx0, vy0, n, t, x0):
    return (((-2 * vx0) * math.sin(n*t)) + ((4 * vy0 + 6 * n * x0) * math.cos(n*t)) - ((3 * vy0 + 6 * n * x0)))

vyf = (solveVyf(vx0, vy0, math.radians(n2), math.radians(t), xo))

print('Question 18: Dvy0 = ' + str(vy0 * 1000))
print('Question 18: Dvx0 = ' + str(vx0 * 1000))
print('Question 18: Dvyf = ' + str(vyf * 1000))
print('Question 18: Dvxf = ' + str(vxf * 1000))

# Write the answer to the solution sheet
f.write('Answer: Dvy0 = ' + str(vy0 * 1000))
f.write('Answer: Dvx0 = ' + str(vx0 * 1000))
f.write('Answer: Dvyf = ' + str(vyf * 1000))
f.write('Answer: Dvxf = ' + str(vxf * 1000))

f.write('-----------------------------------------------------------------')

#####################################
# Question 19
#####################################
f.write('\nQuestion 19:\n')
f.write('Answer : L/D > 0')     
f.write('-----------------------------------------------------------------')

#####################################
# Question 20
#####################################
f.write('\nQuestion 20:\n')
f.write('Answer : to achieve the desired flight path angle at entry as the flight path angle at entry establishes the resulting trajectory')     
f.write('-----------------------------------------------------------------')

#####################################
# Question 21
#####################################
f.write('\nQuestion 21:\n')
f.write('Answer : Ballistic entry trajectory')     
f.write('-----------------------------------------------------------------')

#####################################
# Question 22
#####################################
f.write('\nQuestion 22:\n')
f.write('Answer : Apollo PEG')     
f.write('-----------------------------------------------------------------')










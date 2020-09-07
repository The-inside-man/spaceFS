import math

# initial thrust angle 
initThrustAngle = 75
ita = 75

# time intervals in seconds
t0 = 0
t1 = 100
t2 = 200
t3 = 300
t4 = 400
t5 = 500
t6 = 600

# total burn time in seconds
burnTime = 600

# linear tangent law

print('Answer of Tan Thrust.')
ans1 = (1-(t1/burnTime))* math.tan(math.radians(ita))
print(ans1)
ans1 = (math.degrees(math.atan(ans1)))
ans2 = (1-(t2/burnTime))* math.tan(math.radians(ita))
print(ans2)
ans2 = (math.degrees(math.atan(ans2)))
ans3 = (1-(t3/burnTime))* math.tan(math.radians(ita))
print(ans3)
ans3 = (math.degrees(math.atan(ans3)))
ans4 = (1-(t4/burnTime))* math.tan(math.radians(ita))
print(ans4)
ans4 = (math.degrees(math.atan(ans4)))
ans5 = (1-(t5/burnTime))* math.tan(math.radians(ita))
print(ans5)
ans5 = (math.degrees(math.atan(ans5)))
ans6 = (1-(t6/burnTime))* math.tan(math.radians(ita))
print(ans6)
ans6 = (math.degrees(math.atan(ans6)))

print('Answer of Thrust angle in Degrees.')
print(ans1)
print(ans2)
print(ans3)
print(ans4)
print(ans5)
print(ans6)


# Question 12
# BiLinear Constants
C1 = 0.02
C2 = 800
C3 = 8
C4 = 200

# tan(theta) = (C1 * t + C2) / (C3 * t + C4)
print('This is tan(theta) for Q12')
ans1 = (C1 * t1 + C2) / (C3 * t1 + C4)
print(ans1)
ans2 = (C1 * t2 + C2) / (C3 * t2 + C4)
print(ans2)
ans3 = (C1 * t3 + C2) / (C3 * t3 + C4)
print(ans3)
ans4 = (C1 * t4 + C2) / (C3 * t4 + C4)
print(ans4)
ans5 = (C1 * t5 + C2) / (C3 * t5 + C4)
print(ans5)
ans6 = (C1 * t6 + C2) / (C3 * t6 + C4)
print(ans6)


print('This is Theta rads')
ans1 = math.atan(ans1)
print(ans1)
ans2 = math.atan(ans2)
print(ans2)
ans3 = math.atan(ans3)
print(ans3)
ans4 = math.atan(ans4)
print(ans4)
ans5 = math.atan(ans5)
print(ans5)
ans6 = math.atan(ans6)
print(ans6)

print('This is theta degrees')
print(math.degrees(ans1))
print(math.degrees(ans2))
print(math.degrees(ans3))
print(math.degrees(ans4))
print(math.degrees(ans5))
print(math.degrees(ans6))


def linearTangentLaw(timeInterval, timeEnd, initThrustAngle) :
	return math.degrees(math.atan((1 - (timeInterval/timeEnd)) * math.tan(math.radians(initThrustAngle))))


def convertECEF (latitude, longitude, height, ellipticity, eradiusft) :
	xef = ((eradiusft/math.sqrt(math.pow(math.cos(latitude),2) + math.pow((1 - ellipticity), 2) * math.pow(math.sin(latitude), 2))) + height) * math.cos(latitude) * math.cos(longitude)
	yef = ((eradiusft/math.sqrt(math.pow(math.cos(latitude),2) + math.pow((1 - ellipticity), 2) * math.pow(math.sin(latitude), 2))) + height) * math.cos(latitude) * math.sin(longitude)
	zef = (((math.pow((1 - ellipticity), 2) * eradiusft)/math.sqrt(math.pow(math.cos(latitude),2) + math.pow((1 - ellipticity), 2) * math.pow(math.sin(latitude), 2))) + height) * math.sin(latitude)
	return [xef, yef, zef]

    

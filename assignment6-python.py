import math
import numpy as np


# Provided PEG parameters from text: 5-19, Page 105 - Chapter 5
rnow = [570.0964, -3009.942, 1659.351]
vnow = [0.9286415, -0.1318095, 0.155506]

t=152

rd = [1352.658, -2798.241, 1618.057]

tgo = 366.8735
rgrav = [-72.38909, 293.1735, -163.1685]
S = 516.1065
Q = 77277.45

vgo = [3.553292, 0.1262142, 0.2583443]
magVgo = 3.564906

# nm
Re = 3443.15
# nm^2/sec^2
mu = 62750.2776

r1 = 60 + Re
r2 = 300 + Re
rmeco = 90 + Re

f = open('brown_assignment6_results.txt', 'a')
f.write('Name: Jake Brown')
f.write('\nAssignment 6')
f.write('\n\n')
f.write('----------------------START PROGRAM---------------------')
#Questions 20 + below in PEG fnuction

def peg(rnow, vnow, t, rd, tgo, rgrav, S, Q, vgo, magVgo) :
    print('Generating PEG output...')
    print('----------------------------------------------------------------------')
    # a) Refernce thrust vector of LambdaV
    lamv = np.divide(vgo, magVgo)
    print('lamda V          = ', lamv)
    f.write('\nQuestion 20 (Lambda V) - ' + str(lamv))

    # b) Compute the reference time tlambda
    L = magVgo
    J = L * tgo - S
    tlam = J / L
    print('t lambda         = ', tlam)
    f.write('\nQuestion 21 (t Lambda) - ' + str(tlam))

    # c) Compute reference thrust turning rate lambda prime
    rgo = np.subtract(np.subtract(np.subtract(rd, rnow) ,rgrav), np.multiply(vnow, tgo))

    lamPrime = (np.subtract(rgo, np.multiply(S, lamv)))/(Q - (S*tlam))

    print('Lambda prime     = ', lamPrime)
    f.write('\nQuestion 22 (Lambda Prime) - ' + str(lamPrime))

    # d) Calculate the linear Tangent thrust vector and the unit vector
    lamf = np.add(lamv, np.multiply(lamPrime, (t - (tlam + t))))
    uf = np.divide(lamf, (np.sqrt((lamf[0] * lamf[0]) + (lamf[1] * lamf[1]) + (lamf[2] * lamf[2]))))
    print('Lambda f         = ', lamf)
    print('u sub f          = ', uf)
    f.write('\nQuestion 23 (Lambda f) - ' + str(lamf))
    f.write('\nQuestion 24 (u sub f) - ' + str(uf))


    # e) To Do - Question 26
    ir = [[0.1635154],[-0.8639564],[0.4762794]]
    cosGamma = np.dot(uf, ir)
    print('CosGamma = ', math.degrees(cosGamma))
    f.write('\nQuestion 26 (Cos Gamma - Vertical) - ' + str(math.degrees(cosGamma)))
    f.write('\nQuestion 26 (Cos Gamma - horizontal) - ' + str(90 - math.degrees(cosGamma)))

    # f) vmeco, flight path anlge calculation
    a = (r1+r2)/2
    vmeco = math.sqrt(mu *((2/rmeco) - (1/a)))
    vp = math.sqrt(mu * ((2/r1) - (1/a)))
    magh = np.dot(r1, vp)
    gamma = math.degrees(math.acos(magh/(np.dot(vmeco, rmeco))))
    print('Velocity at MECO = ' + str(vmeco))
    print('Gamma            = ' + str(gamma))
    f.write('\nQuestion 30 (V at MECO) - ' + str(vmeco))
    f.write('\nQuestion 30 (gamma) - ' + str(gamma))
    
    print('----------------------------------------------------------------------')

# assume Earth rotates at 15 degrees per hour
GHAVE = 70

# RAAN is omega
RAAN = 40
inc = 28.5
longitude = 40
latitude = 28.5

# Hours
tepoch = 5
rotrate = 15


def calcTOTime(GHAVE, RAAN, inc, longitude, latitude, tepoch, rotrate):
    t = (1/rotrate)*(math.degrees(math.asin((math.tan(latitude)/math.tan(inc)))) - GHAVE + RAAN - longitude)+tepoch
    return t

peg(rnow, vnow, t, rd, tgo, rgrav, S, Q, vgo, magVgo)

print('Calculation to TO Time: ',calcTOTime(GHAVE, RAAN, inc, longitude, latitude, tepoch, rotrate))
f.write('\nQuestion 31 (Calculation for T/O time) - ' + str(calcTOTime(GHAVE, RAAN, inc, longitude, latitude, tepoch, rotrate)))
f.write('\n-----------------------END PROGRAM----------------------\n')
f.close()

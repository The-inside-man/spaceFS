import numpy as np
import math

# Kalman Filtering
# Iterations

sigmaE = 200
k = 0
Ve = 300100
sigmaM = 5.00

Vm1 = 300006.000
Vm2 = 300003.000
Vm3 = 299997.000
Vm4 = 300007.000
Vm5 = 299995.000

def calcK(sigmaE, sigmaM):
    se2 = math.pow(sigmaE, 2)
    sm2 = math.pow(sigmaM, 2)
    return se2 / (se2 + sm2)

def estVe(Ve, k, Vm):
    return (Ve + (k * (Vm - Ve)))

def calcSigmaE(k, sigmaE):
    return math.sqrt((1 - k) * math.pow(sigmaE, 2))

def printIteration(sigmaE, sigmaM, Vm, Ve):
    k = calcK(sigmaE, sigmaM)
    print('k = ' + str(calcK(sigmaE, sigmaM)))
    print('Meas Velocity = ' + str(Vm))
    eVe = estVe(Ve, calcK(sigmaE, sigmaM), Vm)
    print('Kalman Ve = ' + str(eVe))
    sigE = calcSigmaE(calcK(sigmaE, sigmaM),sigmaE)
    print('new SigmaE = ' + str(sigE))
    return [k, Vm, eVe, sigE]

iter1 = printIteration(sigmaE, sigmaM, Vm1, Ve)
print(iter1)

iter2 = printIteration(iter1[3], sigmaM, Vm2, iter1[2])
print(iter2)

iter3 = printIteration(iter2[3], sigmaM, Vm3, iter2[2])
print(iter3)

iter4 = printIteration(iter3[3], sigmaM, Vm4, iter3[2])
print(iter4)

iter5 = printIteration(iter4[3], sigmaM, Vm5, iter4[2])
print(iter5)


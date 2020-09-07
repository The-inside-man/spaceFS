import math
import numpy as np

def q6(inci, incf, ui, uf):
    return (math.cos(inci) * math.cos(incf) + math.sin(inci) * math.sin(incf) * math.cos(uf -ui))

print(math.degrees(q6(51.8, 51.7, 65.38, 65.48)))

a = math.degrees(q6(51.8, 51.7, 65.38, 65.48))

v = 7.9393

dv = (2 * v) * math.sin(a/2)

vp = math.sqrt(398600/(404.5 + 6378))
vf = math.sqrt(398600/(407 + 6378))
print(vp)
print(vf)

dv = (7.9393 - vf)

print(dv)


T = [[0.5944400, 0.8041300, -0.0025200],[-0.8041300, 0.5944400, 0.0034700],[0.0042900, -0.0000400, 0.9999990]]
ab = [5.0, 15.5, 4.2]


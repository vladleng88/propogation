from lib.atmosphere import Atmosphere
from lib.params import Params
from lib.raytrace import Raytrace
import matplotlib.pyplot as plt


atmosphere = Atmosphere()
params = Params()
atmosphere.setAtmoshere()
raytrace = Raytrace(0, params)
print('nox:', raytrace.getNox())
H = atmosphere.getDefaultHeight()
sH = atmosphere.getSHeight()
T = atmosphere.getDefaultTemperature()
sT = atmosphere.getSTemperature()
t_coeff = atmosphere.getSTemperatureCoeff()
print(list(t_coeff.keys()))
print('----')
print(H)
print(sH)
print(T)
print(sT)
h=5
temp = atmosphere.getTemprature(h)
print('Температура на высоте', h, 'равна', temp, 'градусов Цельсия' )
t = {}
hh = {}
for i in range(1, 20000):
    t[i] = atmosphere.getTemprature(i)
    hh[i] = i

plt.plot(list(t.values()), list(hh.values()), linewidth=2.0)
plt.show()

#for key in sd:
#    print(key, '--', sd[key])
#print('-----------------')
#print('Omega:', params.getOmega())
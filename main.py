from lib.atmosphere import Atmosphere
from lib.params import Params

atmosphere = Atmosphere()
params = Params()
atmosphere.setAtmoshere()

H = atmosphere.getHeight()
d = atmosphere.getDensity()
sd = atmosphere.getSDensity()

print(sd)
print(sd)
for key in sd:
    print(key, '--', sd[key])
print('-----------------')
print('Omega:', params.getOmega())
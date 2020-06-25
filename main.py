from lib.atmosphere import Atmosphere

atmosphere = Atmosphere()
atmosphere.setAtmoshere()

H = atmosphere.getHeight()
d = atmosphere.getDensity()
sd = atmosphere.getSDensity()
print(sd)
print(sd)
for key in sd:
    print(key, '--', sd[key])

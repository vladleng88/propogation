from lib.atmosphere import Atmosphere
from lib.pathway import Pathway
from lib.params import Params
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

atmosphere = Atmosphere()
atmosphere.setAtmoshere()

pathway = Pathway(atmosphere, 5000, 0, 0, 0, 1.2, 1.4, 10, 0, 0)
pathway.setStritelinePathway()

#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ax = Axes3D(fig)
x = pathway.getX()
y = pathway.getY()
z = pathway.getZ()
t = pathway.getTime()
V = pathway.getVelocity()
xx = list(x.values())
yy = list(y.values())
zz = list(z.values())
print(len(xx))
print(len(yy))
print(len(zz))
print(len(t))
print(len(V))
#ax.set_xlabel('x')
#ax.set_ylabel('z')
#ax.set_zlabel('y')
#ax.set_xlim(3000)
#plt.plot(xx, yy, zz, 'gray', zdir='y')
#plt.show()
tau = pathway.getTime()

for i in tau:
    mach = pathway.getVelocity()[i] / pathway.soundSpeed(atmosphere.getTemperature(pathway.getY()[i]))
    y0 = pathway.getY()[i]
    x0 = pathway.getX()[i]
    t0 = tau[i]= pathway.getX()[i]
    z0 = pathway.getZ0()
    o = pathway.getOmega()
    o1 = pathway.getOmega1()
    param = Params(y0, x0, z0, t0, o, o1, mach)
    print('i', i, 'x:',pathway.getX()[i], 'y:',pathway.getY()[i], 'z:',pathway.getZ()[i], 't:', tau[i], 'mach:',mach, 'velocity=',pathway.getVelocity()[i] )


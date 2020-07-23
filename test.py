from lib.atmosphere import Atmosphere
from lib.pathway import Pathway
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

atmosphere = Atmosphere()
atmosphere.setAtmoshere()

pathway = Pathway(atmosphere, 15500, 10, 0, 0, 1.2, 1.4, 10, 10, 10)
pathway.setStritelinePathway()

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
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
ax.set_xlabel('x')
ax.set_ylabel('z')
ax.set_zlabel('y')
ax.set_xlim(3000)
plt.plot(xx, yy, zz, 'gray', zdir='y')
plt.show()
ss = summ(3,5)
print(ss)


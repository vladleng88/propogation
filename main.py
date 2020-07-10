from lib.atmosphere import Atmosphere
from lib.params import Params
from lib.raytrace import Raytrace
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

atmosphere = Atmosphere()
params = Params()
atmosphere.setAtmoshere()
raytrace = Raytrace(12, params, atmosphere)
print('nz',raytrace.getNoz())
xx=[]
yy=[]
zz=[]
y_rev = 0
y=params.getY0()
dt = raytrace.getDt()
t = params.gett0()
x = params.getX0()
z = params.getZ0()
xx.append(x)
yy.append(y)
zz.append(z)
ny = raytrace.getNoy()
y_next = y+raytrace.getDy(y, ny)
Integrals0 = raytrace.getInitialIntagrals()
I1 = Integrals0['I1']
I2 = Integrals0['I2']
I3 = Integrals0['I3']
I4 = Integrals0['I4']
while y >= params.getYtarget() and y_next >= params.getYtarget():
    dy = raytrace.getDy(y, ny)
    y_rev = y_rev + dy
    y = params.getY0() + y_rev
    ny_next = raytrace.n_definition(y, ny)['ny']
    dy_next = raytrace.getDy(y, ny_next)
    y_rev_next = y_rev + dy_next
    y_next = params.getY0() + y_rev_next
    #print('y=', y, 'y_next=', y_next)
    if y_next >= params.getYtarget():
        x = x + raytrace.getDx(y, ny, dy_next)
        xx.append(x)
        z = z + raytrace.getDz(y, ny, dy_next)
        zz.append(z)
        yy.append(y)
        t = t + dt
        print('x==', x, 'y==', y, 'z==', z)
        Integrals = raytrace.getDIntegrals(y, ny, dy)
    else:
        dy = y - dy
        print('dy', dy)
        y = params.getYtarget()
        u_final = atmosphere.getWindY(y)+raytrace.soundSpeed(atmosphere.getTemperature(y))*raytrace.n_definition(y,ny)['ny']
        dt_final = dy/u_final
        #x = x + atmosphere.getWindX(y)+raytrace.soundSpeed(atmosphere.getTemperature(y))*raytrace.n_definition(y, ny)['nx']
        x = x + raytrace.getDx(y, ny, dy)
        #z = z + atmosphere.getWindZ(y)+raytrace.soundSpeed(atmosphere.getTemperature(y))*raytrace.n_definition(y, ny)['nz']
        z = z + raytrace.getDz(y, ny, dy)
        t = t + dt_final
        xx.append(x)
        zz.append(z)
        yy.append(params.getYtarget())
        print('x=', x, 'y=', y, 'z=', z)
        Integrals = raytrace.getDIntegrals(y, ny, dy)
    I1 = I1 + Integrals['dI1']
    I2 = I2 + Integrals['dI2']
    I3 = I3 + Integrals['dI3']
    I = raytrace.getIntegral(I1, I2, I3)
    I4 = I4 + raytrace.getDIntegral4(y, ny, I, dy)
    k = raytrace.getK(I4)
    k1 = raytrace.getK1(params.getYtarget(), I, ny)
    k2 = raytrace.getK2()
    ny = raytrace.n_definition(y, ny)['ny']

    #print('Высота:',y,'; k = ', k, '; k1 = ', k1, '; k2 = ', k2, '; I1 = ', I1, '; I2 = ', I2, '; I3 = ', I3, '; I4 = ', I4)
t = {}
hh = {}
for i in range(1, 20000):
    t[i] = atmosphere.getTemperature(i)
    hh[i] = i

fig = plt.figure()
ax = plt.axes(projection="3d")
ax.plot3D(xx, yy, zz, 'gray', zdir='y')
plt.show()
plt.plot(list(t.values()), list(hh.values()), linewidth=2.0)
#plt.show()

#for key in sd:
#    print(key, '--', sd[key])
#print('-----------------')
#print('Omega:', params.getOmega())
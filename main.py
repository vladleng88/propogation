from lib.atmosphere import Atmosphere
from lib.params import Params
from lib.raytrace import Raytrace
from lib.aerodynamics import Aerodynamics
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from math import *

atmosphere = Atmosphere()

params = Params()
atmosphere.setAtmoshere()
raytrace = Raytrace(0, params, atmosphere)
aerodynamics = Aerodynamics(params, raytrace.getFlightMachNumber())
coeff_normal = 0.5*atmosphere.getDensity(params.getY0())*(raytrace.getFlightMachNumber()**2)*\
    (raytrace.soundSpeed(atmosphere.getTemperature(params.getY0())))**2*sqrt(params.getLength())/\
    (sqrt((raytrace.getFlightMachNumber())**2 - 1)*params.getLift()print('coeff_normal=', coeff_normal)
aerodynamics.setAerodynmamicsData(coeff_normal)
ettaRef = aerodynamics.getEttaRef()
withemRef = aerodynamics.getWithemRef()
potentialRef = aerodynamics.getPotentalRef()
dSdX = aerodynamics.getdSdX()
fig, axs = plt.subplots(3, 2)
axs[0, 0].plot(list(ettaRef.values()), list(dSdX.values()), linewidth=2.0, color='red')
axs[0, 0].set_xlabel('etta')  # Add an x-label to the axes.
axs[0, 0].set_ylabel('dS/dX')  # Add a y-label to the axes.

axs[0, 1].plot(list(ettaRef.values()), list(potentialRef.values()), linewidth=2.0, color='red')
axs[0, 1].set_xlabel('etta')  # Add an x-label to the axes.
axs[0, 1].set_ylabel('Potential')  # Add a y-label to the axes.

axs[2, 0].plot(list(ettaRef.values()), list(withemRef.values()), linewidth=2.0, color='red')
axs[2, 0].set_xlabel('etta')  # Add an x-label to the axes.
axs[2, 0].set_ylabel('Withem')  # Add a y-label to the axes.
#plt.show()
P_BANG = {}
T_BANG = {}
FF_DOP = {}
T_DOP = {}
xx = []
yy = []
zz = []
k = 0
k1 = 0
k2 = 0
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
Integrals = raytrace.getInitialIntagrals()
I1 = Integrals['I1']
I2 = Integrals['I2']
I3 = Integrals['I3']
I4 = Integrals['I4']
print('ny_start', ny)
print('I1_start', I1)
print('I2_start', I2)
print('I3_start', I3)
print('I4_start', I4)
i1_test = I1
M_eff = 1/sqrt(1-ny**2)
while y >= params.getYtarget() and y_next >= params.getYtarget():
    dy = raytrace.getDy(y, ny)
    y_rev = y_rev + dy
    y = params.getY0() + y_rev
    ny_next = raytrace.n_definition(y, ny)['ny']
    dy_next = raytrace.getDy(y, ny_next)
    y_rev_next = y_rev + dy_next
    y_next = params.getY0() + y_rev_next
    #print('y', y, 'dy', dy, 'ynext', y_next)
    if y_next >= params.getYtarget():
        x = x + raytrace.getDx(y, ny, dy_next)
        xx.append(x)
        z = z + raytrace.getDz(y, ny, dy_next)
        zz.append(z)
        yy.append(y)
        t = t + dt
    else:
        #dy = y - dy
        dy = y-params.getYtarget() - dy
        y = params.getYtarget()
        u_final = atmosphere.getWindY(y)+raytrace.soundSpeed(atmosphere.getTemperature(y))*raytrace.n_definition(y,ny)['ny']
        dt_final = dy/u_final
        x = x + raytrace.getDx(y, ny, dy)
        z = z + raytrace.getDz(y, ny, dy)
        t = t + dt_final
        xx.append(x)
        zz.append(z)
        yy.append(params.getYtarget())
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

    #print('h=', y, ' k=', k, ' k1=', k1, ' k2=', k2, 'I1=', I1,'I2=', I2,'I3=', I3,'I4=', I4, 'I', I)
    #i1_test += (raytrace.soundSpeed(atmosphere.getTemperature(y)))*dy/ny
    #print('h=', y, ' k=', k, ' k1=', k1, ' k2=', k2, 'I1=', I1, 'ny', ny, 'a', raytrace.soundSpeed(atmosphere.getTemperature(y)), 'i1_test', i1_test)
print(' k=', k, ' k1=', k1, ' k2=', k2, )
print('sigma', raytrace.getSigma())
print('I1_ref', I1*sqrt(M_eff**2 - 1)/raytrace.soundSpeed(atmosphere.getTemperature(params.getY0()))/params.getY0()/M_eff, 'M_eff', M_eff)
print('I2_ref', I2*(sqrt(M_eff**2 - 1))**3/raytrace.soundSpeed(atmosphere.getTemperature(params.getY0()))/params.getY0()/M_eff**3, 'M_eff', M_eff)
for j in range(0, len(withemRef)):
   #P_BANG[j] = k1 * withemRef[j]
   #T_BANG[j] = k2 * (ettaRef[j] - k * withemRef[j])
   FF_DOP[j] = 2*potentialRef[j]/k - withemRef[j]**2
   T_DOP[j] = ettaRef[j] - k*withemRef[j]
   #print('t_dop', T_DOP[j],'p_dop', F_DOP[j])
#Take the minimum and maximum value of the duration
T_MIN = T_DOP[0]
T_MAX = T_DOP[0]
for j in range(0, len(withemRef)):
    if T_MIN > T_DOP[j]:
        T_MIN = T_DOP[j]
    if T_MAX < T_DOP[j]:
        T_MAX = T_DOP[j]
DT = (T_MAX - T_MIN)/(params.getIspace()-1)
P_BOOM = {}
T_BOOM = {}
P_HALF_BOOM = {}
nPT = 1
T = {}
FF_wave = {}
F_wave = {}
print('T_MIN=', T_MIN, ' T_MAX=', T_MAX)
print('DT = ', DT)
for j in range(0, (params.getIspace())):
    T[j] = T_MIN + (j)*DT
    #FF_wave[j] = 0
F_wave[0] = 0
count = 0
t_start = 0
for i in range(1, len(withemRef)):
    for j in range(0, len(T)):
        if T[j]>=T_DOP[i-1] and T[j]<T_DOP[i]:
            a = (withemRef[i]-withemRef[i-1])/(ettaRef[i] - ettaRef[i-1])
            b = withemRef[i-1] + 2 * (potentialRef[i] - potentialRef[i-1])/(ettaRef[i] - ettaRef[i-1]) - withemRef[i]
            FF = FF_DOP[i-1] + (T[j] - T_DOP[i-1]) * (b - 2*a*k*withemRef[i-1] + a*(T[j] - T_DOP[i-1]))/(k * (1-a*k))
            #print('j', j,'t_wave', T[j], 'ff', FF)
            if t_start == 0:
                t_start = T[j]
            if FF_wave.get(T[j]) or FF_wave.get(T[j]) == 0:
                if FF > FF_wave[T[j]]:
                    FF_wave[T[j]] = FF
                    #F_wave[j] = (FF_wave[j] - FF_wave[j-1])/(T[j] - T_DOP[i-1])
            else:
                FF_wave[T[j]] = FF
TT = sorted(list(FF_wave.keys()))
ttt = {}
FFF = {}
counter = 0
#print(len(TT))
#for tt in TT:
#    if tt<t_start:
#        FF_wave[tt] = 0
#counter = 0
for tt in TT:
    if tt>=t_start and tt<2.1:
        FFF[counter] = FF_wave[tt]
        ttt[counter] = tt
        counter+=1
        #print(tt, '\t', FF_wave[tt])
print('t_start', t_start)

axs[1, 0].plot(list(T_DOP.values()), list(FF_DOP.values()), linewidth=2.0, color='red')
axs[1, 0].plot(list(ttt.values()), list(FFF.values()), linewidth=2.0, color='blue')
axs[1, 0].set_xlabel('t_wave')  # Add an x-label to the axes.
axs[1, 0].set_ylabel('Ф_wave')  # Add a y-label to the axes.
axs[1, 0].set_title("FFF_WAVE")  # Add a title to the axes.

#exit(0)
F_wave[0] = 0
for i in range(1, len(FFF)):
    F_wave[i] = (k/2)*(FFF[i] - FFF[i - 1]) / (ttt[i] - ttt[i - 1])
    #print('i', i, 't', ttt[i], 'F_wave', F_wave[i])
print('len_FFF', len(FFF), 'ttt_finish', ttt[len(FFF)-1])

#axs[0, 1].plot(list(ttt.values()), list(F_wave.values()), linewidth=2.0, color='red')
#axs[0, 1].set_xlabel('t_wave')  # Add an x-label to the axes.
#axs[0, 1].set_ylabel('withem')  # Add a y-label to the axes.
#axs[0, 1].set_title("Withem_transform")  # Add a title to the axes.
#withemRefnew = {}
#ettaRefnew = {}
#WT = {}
#for i in range(1, len(F_wave)):
#    for j in range(0, len(ettaRef)):
#        #print('ettaRef = ', ettaRef[j])
#        if ettaRef[j]>=ttt[i-1] and ettaRef[j]<ttt[i]:
#            a = (F_wave[i] - F_wave[i-1])/(ttt[i]-ttt[i-1])
#            withemRefnew[j] = F_wave[i-1] + a * (ettaRef[j]-ttt[i-1])
#            ettaRefnew[j] = ettaRef[j] + k*withemRefnew[j]
#            #print('j', j, 'ettaRef=', ettaRefnew[j], 'withem_new', withemRefnew[j])

#print('len_withem_new', len(withemRefnew))
#print('len_etta', len(ettaRef))
for j in range(0, len(F_wave)):
    P_BOOM[j] = k1*F_wave[j]*2
    T_BOOM[j] = k2*(ttt[j])
    #print('j',j ,'\t', 'etta=', ettaRef[j], '\t','t=', T_BOOM[j], '\t','p=', P_BOOM[j])
axs[1, 1].plot(list(T_BOOM.values()), list(P_BOOM.values()), linewidth=2.0, color='red')
axs[1, 1].set_xlabel('t_boom')  # Add an x-label to the axes.
axs[1, 1].set_ylabel('p_boom')  # Add a y-label to the axes.
#plt.show()


#plt.plot(list(T_DOP.values()), list(F_DOP.values()), linewidth=2.0, color='red')
#plt.plot(list(T.values()), list(FF_wave.values()), linewidth=2.0, color='blue')
#plt.show()
#fig = plt.figure()
#axs[2, 1] = plt.axes(projection="3d")
#axs[2, 1].add_subplot(xx, yy, zz, 'gray', zdir='y')
plt.show()


#plt.plot(list(t.values()), list(hh.values()), linewidth=2.0)
#plt.show()

#for key in sd:
#    print(key, '--', sd[key])
#print('-----------------')
#print('Omega:', params.getOmega())
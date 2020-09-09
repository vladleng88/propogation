from lib.atmosphere import Atmosphere
from lib.params import Params
from lib.raytrace import Raytrace
from lib.aerodynamics import Aerodynamics
from lib.pathway import Pathway
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from math import *
from lib.praphview import *

atmosphere = Atmosphere()
atmosphere.setAtmoshere()

pathway = Pathway(atmosphere, 15500, 0, 0, 0, 1.2, 1.4, 10, 0, 0)
pathway.setStritelinePathway()
fig, ax = plt.subplots()
tauDict = pathway.getTime()
for key in tauDict:
    print('------------------------------------------start raytrace at t=', tauDict[key],'----------------------------------------------')
    mach = pathway.getVelocity()[key]/pathway.soundSpeed(atmosphere.getTemperature(pathway.getY()[key]))
    print('-------y0=',pathway.getY()[key], 'x0=',pathway.getX()[key], 'z0=',pathway.getZ()[key], 't0=',tauDict[key], 'o0=',pathway.getOmega(), 'o1=',pathway.getOmega1(), 'M0=',mach)
    params = Params(pathway.getY()[key], pathway.getX()[key], pathway.getZ()[key], tauDict[key], pathway.getOmega(), pathway.getOmega1(), mach)
    a = {'ax': pathway.getAx()[key],'ay': pathway.getAy()[key], 'az': pathway.getAz()[key]}
    raytrace = Raytrace(0, params, atmosphere, a, pathway.getDaDh()[key], pathway.getDVDh()[key])
    aerodynamics = Aerodynamics(params, raytrace.getFlightMachNumber())
    coeff_normal = 0.5*atmosphere.getDensity(params.getY0())*(raytrace.getFlightMachNumber()**2)* \
                   (raytrace.soundSpeed(atmosphere.getTemperature(params.getY0())))**2*sqrt(params.getLength())/ \
                   (sqrt((raytrace.getFlightMachNumber())**2 - 1)*params.getLift())
    #print(sqrt((raytrace.getFlightMachNumber())**2 - 1))
    print('coeff_normal=', coeff_normal)
    aerodynamics.setAerodynmamicsData(coeff_normal)
    ettaRef = aerodynamics.getEttaRef()
    withemRef = aerodynamics.getWithemRef()
    potentialRef = aerodynamics.getPotentalRef()
    dSdX = aerodynamics.getdSdX()
    P_BANG = {}; T_BANG = {}
    FF_DOP = {}; T_DOP = {}
    xx = []; yy = []; zz = []; II=[]
    k = 0; k1 = 0; k2 = 0
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
    I1 = Integrals['I1']; I2 = Integrals['I2']; I3 = Integrals['I3']; I4 = Integrals['I4']
    II.append(raytrace.getIntegral(I1, I2, I3))
    while y >= params.getYtarget() and y_next >= params.getYtarget():
        # print('sds')
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
            #print('11111')
        else:
            #dy = y - dy
            dy = y-params.getYtarget() - dy
            y = params.getYtarget()
            u_final = atmosphere.getWindY(y)+raytrace.soundSpeed(atmosphere.getTemperature(y))*raytrace.n_definition(y, ny)['ny']
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
        II.append(I)
        I4 = I4 + raytrace.getDIntegral4(y, ny, I, dy)
        #print('I4=', I4)
        k = raytrace.getK(I4)
        k1 = raytrace.getK1(params.getYtarget(), I, ny)
        #print('k1=', k1)
        k2 = raytrace.getK2()
        ny = raytrace.n_definition(y, ny)['ny']

        #print('h=', y, ' k=', k, ' k1=', k1, ' k2=', k2, 'I1=', I1,'I2=', I2,'I3=', I3,'I4=', I4, 'I', I)
        #print('h=', y, ' k=', k, ' k1=', k1, ' k2=', k2, 'I1=', I1, 'ny', ny, 'a', raytrace.soundSpeed(atmosphere.getTemperature(y)), 'i1_test', i1_test)

    #IntegralPlot(II, yy)
    #print(II)
    ax.plot(xx, yy, linewidth=1.0, color='red')
    print(' k=', k, ' k1=', k1, ' k2=', k2, 'I1=', I1,'I2=', I2,'I3=', I3,'I4=', I4, 'I', I, 'sigma', raytrace.getSigma(), 'gamma', raytrace.getGamma() )

    for j in range(0, len(withemRef)):
        FF_DOP[j] = 2*potentialRef[j]/k - withemRef[j]**2
        T_DOP[j] = ettaRef[j] - k*withemRef[j]
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
    #print('T_MIN=', T_MIN, ' T_MAX=', T_MAX)
    #print('DT = ', DT)
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
                if t_start == 0:
                    t_start = T[j]
                if FF_wave.get(T[j]) or FF_wave.get(T[j]) == 0:
                    if FF > FF_wave[T[j]]:
                        FF_wave[T[j]] = FF
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
    #print('t_start', t_start)

    F_wave[0] = 0
    for i in range(1, len(FFF)):
        F_wave[i] = (k/2)*(FFF[i] - FFF[i - 1]) / (ttt[i] - ttt[i - 1])

    for j in range(0, len(F_wave)):
        P_BOOM[j] = k1*F_wave[j]*2
        T_BOOM[j] = k2*(ttt[j])
        #print('j',j ,'\t', 'etta=', ettaRef[j], '\t','t=', T_BOOM[j], '\t','p=', P_BOOM[j])
    #graphPlot(ettaRef, dSdX, potentialRef, withemRef, T_BOOM, P_BOOM,  T_DOP, FF_DOP, ttt, FFF, tauDict[key])
    #if key == 1:
        #break

plt.show()

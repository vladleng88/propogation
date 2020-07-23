from math import *
from lib.atmosphere import Atmosphere

class Pathway:

    def __init__(self, atmosphere: Atmosphere,
                 y0=15500,#метров
                 x0=0,
                 z0=0,
                 t0=0,
                 startM=1.2,
                 finishM=1.4,
                 acceleration=10, # м/с2
                 omega=0,
                 omega1=0,
                 kappa=1.4,
                 R=287.3,
                 dt=0.05
                 ):
        self.__y0 = y0
        self.__x0 = x0
        self.__z0 = z0
        self.__t0 = t0
        self.__M0 = startM
        self.__M1 = finishM
        self.__omega = omega
        self.__R = R
        self.__kappa = kappa
        self.__omega1 = omega1
        self.__acceleration = acceleration
        self.__atmosphere = atmosphere
        self.__dt = dt
        self.__t = {}
        self.__x = {}
        self.__y = {}
        self.__z = {}
        self.__V = {}

    def setStritelinePathway(self):
        V = {}
        x = {}
        y = {}
        z = {}
        t = {}
        Vx = {}
        Vy = {}
        Vz = {}
        T0 = self.__atmosphere.getTemperature(self.__y0)
        V[0] = self.getFlightMachNumber(self.__M0)*self.soundSpeed(T0)
        Vx[0] = V[0]*cos(radians(self.__omega))*cos(radians(self.__omega1))
        #print('Vx0', Vx[0])
        Vz[0] = V[0]*cos(radians(self.__omega))*sin(radians(self.__omega1))
        #print('Vz0', Vz[0])
        Vy[0] = V[0]*sin(radians(self.__omega))
        #print('Vy0', Vy[0])
        x[0] = self.__x0
        y[0] = self.__y0
        z[0] = self.__z0
        t[0] = self.__t0
        mach = V[0] / self.soundSpeed(T0)
        i=1
        while mach < self.__M1:
            V[i] = V[i-1] + self.__acceleration*self.__dt
            Vx[i] = V[i] * cos(radians(self.__omega)) * cos(radians(self.__omega1))
            #print('Vx=',  Vx[i])
            Vz[i] = V[i] * cos(radians(self.__omega)) * sin(radians(self.__omega1))
            #print('Vz', Vz[i])
            Vy[i] = V[i] * sin(radians(self.__omega))
            ax = (Vx[i] + Vx[i-1])/self.__dt
            ay = (Vy[i] + Vy[i-1])/self.__dt
            az = (Vz[i] + Vz[i-1])/self.__dt
            x[i] = x[i-1] + Vx[i-1]*self.__dt + 0.5*ax*(self.__dt)**2
            y[i] = y[i-1] + Vy[i-1]*self.__dt + 0.5*ay*(self.__dt)**2
            z[i] = z[i-1] + Vz[i-1]*self.__dt + 0.5*az*(self.__dt)**2
            t[i] = t[i-1] + self.__dt
            temperature = self.__atmosphere.getTemperature(y[i])
            mach = V[i] / self.soundSpeed(temperature)
            #print('z', z[i])
            if mach < self.__M1:
                i+=1
            else:

                dMdt = ((V[i]/self.__atmosphere.getTemperature(y[i])) - (V[i-1]/self.__atmosphere.getTemperature(y[i-1])))/(t[i]-t[i-1])
                t_final = t[i-1]+(self.__M1-(V[i-1]/self.soundSpeed(self.__atmosphere.getTemperature(y[i-1]))))/(dMdt)
                print('1', self.__M1)
                print('2', V[i-1]/self.__atmosphere.getTemperature(y[i-1]))
                print('t_final', t_final)
                y[i] = y[i-1] + ((y[i] - y[i - 1])/(t[i] - t[i - 1]))*(t_final - t[i-1])
                x[i] = x[i-1] + ((x[i] - x[i - 1])/(t[i] - t[i - 1]))*(t_final - t[i-1])
                z[i] = z[i-1] + ((z[i] - z[i - 1])/(t[i] - t[i - 1]))*(t_final - t[i-1])
                V[i] = self.__M1 * self.soundSpeed(self.__atmosphere.getTemperature(y[i]))
                mach = self.__M1
        self.__t = t
        self.__x = x
        self.__y = y
        self.__z = z
        self.__V = V

    def getFlightMachNumber(self, M0):
        omega = self.__omega
        omega1 = self.__omega1
        T0 = self.__atmosphere.getTemperature(self.__y0)

        V = M0 * self.soundSpeed(T0)
        Vz = V * cos(radians(omega)) * sin(radians(omega1))
        Vx = V * cos(radians(omega)) * cos(radians(omega1))
        Vy = V * sin(radians(omega))
        windX = self.__atmosphere.getWindX(self.__y0)
        windY = self.__atmosphere.getWindY(self.__y0)
        windZ = self.__atmosphere.getWindZ(self.__y0)
        dVx = Vx - windX
        dVy = Vy - windY
        dVz = Vz - windZ
        V_new = sqrt(dVx ** 2 + dVy ** 2 + dVz ** 2)
        M = V_new / self.soundSpeed(T0)

        return M

    def soundSpeed(self, T):
        return sqrt(self.__kappa * self.__R * T)

    def getX(self):
        return self.__x

    def getX0(self):
        return self.__x0
    def setX0(self, x0):
        self.__x0 = x0

    def getY(self):
        return self.__y

    def getY0(self):
        return self.__y0
    def setY0(self, y0):
        self.__y0 = y0

    def getZ(self):
        return self.__z

    def getZ0(self):
        return self.__z0
    def setZ0(self, z0):
        self.__z0 = z0

    def gett0(self):
        return self.__t0
    def sett0(self, t0):
        self.__t0 = t0

    def getKappa(self):
        return self.__kappa
    def setKappa(self, kappa):
        self.__kappa = kappa

    def getR(self):
        return self.__R
    def setR(self, R):
        self.__R = R

    def getAcceleration(self):
        return self.__acceleration
    def setAcceleration(self, a):
        self.__acceleration = a

    def getDt(self):
        return self.__dt
    def setDt(self, dt):
        self.__dt = dt

    def getTime(self):
        return self.__t

    def getVelocity(self):
        return self.__V


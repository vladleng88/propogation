from math import *
from lib.atmosphere import Atmosphere
from lib.params import Params


class Raytrace:

    def __init__(self, tetta, params: Params, atmosphere: Atmosphere):
        self.__tetta = tetta
        self.__params = params
        self.__atmosphere = atmosphere
        self.__nox = self.__nox_definition()
        self.__noy = self.__noy_definition()
        self.__noz = self.__noz_definition()
        self.__dt = 0.05

    def __nox_definition(self):
        omega = self.__params.getOmega()
        omega1 = self.__params.getOmega1()
        myu = self.getMyu()
        nox = sin(radians(myu)) * cos(radians(omega1)) * cos(radians(omega)) - \
              cos(radians(myu)) * sin(radians(omega1)) * sin(radians(self.__tetta)) + \
              cos(radians(myu)) * cos(radians(omega1)) * sin(radians(omega)) * cos(radians(self.__tetta))
        return nox

    def __noy_definition(self):
        omega = self.__params.getOmega()
        myu = self.getMyu()
        noy = sin(radians(myu)) * sin(radians(omega)) - cos(radians(myu)) * cos(radians(omega)) * cos(
            radians(self.__tetta))
        return noy

    def __noz_definition(self):
        omega = self.__params.getOmega()
        omega1 = self.__params.getOmega1()
        myu = self.getMyu()

        noz = sin(radians(myu)) * sin(radians(omega1)) * cos(radians(omega)) + \
              cos(radians(myu)) * cos(radians(omega1)) * sin(radians(self.__tetta)) + \
              cos(radians(myu)) * sin(radians(omega1)) * sin(radians(omega)) * cos(radians(self.__tetta))
        return noz

    def n_definition(self, y, ny_prev=-1):
        n = {}
        nox = self.getNox()
        noy = self.getNoy()
        noz = self.getNoz()
        T0 = self.__atmosphere.getTemperature(self.__params.getY0())
        dWindX = self.__atmosphere.getWindX(y) - self.__atmosphere.getWindX(self.__params.getY0())
        dWindZ = self.__atmosphere.getWindZ(y) - self.__atmosphere.getWindZ(self.__params.getY0())
        a = self.soundSpeed(self.__atmosphere.getTemperature(y))
        a_star = self.soundSpeed(T0) + self.__atmosphere.getWindY(self.__params.getY0()) * noy - \
                 (dWindX * nox + dWindZ * noz)
        # Вычисление коэффициентов квадратного уравнения для определения ny
        c0 = (1 - noy ** 2) * (self.__atmosphere.getWindY(y)) ** 2 + a_star ** 2
        c1 = 2 * (1 - noy ** 2) * (self.__atmosphere.getWindY(y)) * a
        c2 = (1 - noy ** 2) * (a) ** 2 - a_star ** 2
        #print('y=', y, 'd=', c1 ** 2 - 4 * c0 * c2)
        d = sqrt(c1 ** 2 - 4 * c0 * c2)
        ny1 = (-1 * c1 + d) / (2 * c0)
        ny2 = (-1 * c1 - d) / (2 * c0)
        if ny1*ny_prev > 0:
            n['ny'] = ny1
        elif ny2*ny_prev > 0:
            n['ny'] = ny2
        else:
            exit('Ошибка с определнием ny')
        nx = nox * (a + self.__atmosphere.getWindY(y) * n['ny']) / a_star
        n['nx'] = nx
        nz = noz * (a + self.__atmosphere.getWindZ(y) * n['ny']) / a_star
        n['nz'] = nz
        return n

    def getSigma(self):
        noy = self.getNoy()
        omega = radians(self.__params.getOmega())
        myu = self.getMyu()
        sigma = (cos(myu)*sin(omega)+sin(myu)*cos(omega)*cos(radians(self.__tetta)))/(fabs(cos(myu)*sin(omega)+sin(myu)*cos(omega)*cos(radians(self.__tetta)))) * \
                (noy/sqrt(1-noy**2)) * (cos(omega)*sin(radians(self.__tetta)))/sqrt(1-(cos(omega))**2*(sin(radians(self.__tetta)))**2)
        return sigma
    def getGamma(self):
        return 0

    def getNox(self):
        return self.__nox

    def getNoy(self):
        return self.__noy

    def getNoz(self):
        return self.__noz

    def getDt(self):
        return self.__dt

    def getFlightMachNumber(self):
        y0 = self.__params.getY0()
        omega = self.__params.getOmega()
        omega1 = self.__params.getOmega1()
        T0 = self.__atmosphere.getTemperature(y0)
        M0 = self.__params.getMachNumber0()
        V = M0 * self.soundSpeed(T0)
        Vx = V * cos(radians(omega)) * cos(radians(omega1))
        Vy = V * sin(radians(omega))
        Vz = V * cos(radians(omega)) * sin(radians(omega1))
        windX = self.__atmosphere.getWindX(y0)
        windY = self.__atmosphere.getWindY(y0)
        windZ = self.__atmosphere.getWindZ(y0)
        dVx = Vx - windX
        dVy = Vy - windY
        dVz = Vz - windZ
        V_new = sqrt(dVx ** 2 + dVy ** 2 + dVz ** 2)
        M = V_new / self.soundSpeed(T0)
        return M

    def getMyu(self):
        myu = degrees(asin(1 / self.getFlightMachNumber()))
        return myu

    def getDy(self, y, ny_prev):
        return (self.__atmosphere.getWindY(y) + self.soundSpeed(self.__atmosphere.getTemperature(y))*self.n_definition(y, ny_prev)['ny']) * self.getDt()

    def getDx(self, y, ny_prev, dy_next):
        n = self.n_definition(y, ny_prev)
        T = self.__atmosphere.getTemperature(y)
        Vy = self.__atmosphere.getWindY(y)
        k1 = (self.__atmosphere.getWindX(y)) + self.soundSpeed(T) * n['nx']
        #dy_half_h = (self.soundSpeed(T) * n['ny'] + Vy) * self.getDt()/2
        dy_half_h = dy_next/2
        k2 = (self.__atmosphere.getWindX(fabs(y + dy_half_h))) + self.soundSpeed(self.__atmosphere.getTemperature(fabs(y + dy_half_h))) * n['nx']
        k3 = (self.__atmosphere.getWindX(fabs(y + dy_half_h))) + self.soundSpeed(self.__atmosphere.getTemperature(fabs(y + dy_half_h))) * n['nx']
        #dy_h = (self.soundSpeed(T) * n['ny'] + Vy) * self.getDt()
        dy_h = dy_next
        k4 = (self.__atmosphere.getWindX(fabs(y + dy_h))) + self.soundSpeed(self.__atmosphere.getTemperature(y + dy_h)) * n['nx']
        dx = self.getDt()*(k1+2*k2+2*k3+k4)/6
        return dx

    def getDz(self, y, ny_prev, dy_next):
        n = self.n_definition(y, ny_prev)
        T = self.__atmosphere.getTemperature(y)
        Vy = self.__atmosphere.getWindY(y)
        k1 = (self.__atmosphere.getWindZ(y)) + self.soundSpeed(T) * n['nz']
        #dy_half_h = (self.soundSpeed(T) * n['ny'] + Vy) * self.getDt()/2
        dy_half_h = dy_next / 2
        k2 = (self.__atmosphere.getWindZ(fabs(y + dy_half_h))) + self.soundSpeed(self.__atmosphere.getTemperature(fabs(y + dy_half_h))) * n['nz']
        k3 = (self.__atmosphere.getWindZ(fabs(y + dy_half_h))) + self.soundSpeed(self.__atmosphere.getTemperature(fabs(y + dy_half_h))) * n['nz']
        #dy_h = (self.soundSpeed(T) * n['ny'] + Vy) * self.getDt()
        dy_h = dy_next
        k4 = (self.__atmosphere.getWindZ(fabs(y + dy_h))) + self.soundSpeed(self.__atmosphere.getTemperature(fabs(y + dy_h))) * n['nz']
        dz = self.getDt()*(k1+2*k2+2*k3+k4)/6
        return dz

    def getDIntegrals(self, y, ny_prev, dy):
        a0 = self.soundSpeed(self.__atmosphere.getTemperature(self.__params.getY0()))
        a = self.soundSpeed(self.__atmosphere.getTemperature(y))
        T = self.__atmosphere.getTemperature(y)
        nox = self.getNox()
        noy = self.getNoy()
        noz = self.getNoz()
        dWindX = self.__atmosphere.getWindX(y) - self.__atmosphere.getWindX(self.__params.getY0())
        dWindZ = self.__atmosphere.getWindZ(y) - self.__atmosphere.getWindZ(self.__params.getY0())
        ny = self.n_definition(y, ny_prev)['ny']
        a_star = a0 + self.__atmosphere.getWindY(self.__params.getY0()) * noy - (dWindX * nox + dWindZ * noz)
        l0x = -1*noz
        l0z = nox
        # коэффициент u = Vy + a*ny
        u = self.__atmosphere.getWindY(y) + self.soundSpeed(T)*ny
        # коэффициент K = Un/Un0 = (a+Vy*ny)/a_star
        K = (a + self.__atmosphere.getWindY(y) * ny)/a_star
        #dI1 = (a0*(1+(K**2)*((l0x*dWindX+l0z*dWindZ)**2)/(u**2))*a*K/u)*dy
        dI1 = (a0*a*K/u)*dy
        #print('dI1', (a0*(1+(K**2)*((l0x*dWindX+l0z*dWindZ)**2)/(u**2))*a*K/u))

        dI2 = ((a0**3)*a*(K**3)/(u**3))*dy
        dI3 = ((a0**2)*a*(l0x*dWindX+l0z*dWindZ)*K/(u**3))*dy
        dI = {}
        dI['dI1'] = dI1
        dI['dI2'] = dI2
        dI['dI3'] = dI3
        return dI
    def getIntegral(self, I1, I2, I3):
        T0 = self.__atmosphere.getTemperature(self.__params.getY0())
        a0 = self.soundSpeed(T0)
        noy = self.getNoy()
        # коэффициент u0 = Vy0 + a0*ny0
        u0 = self.__atmosphere.getWindY(self.__params.getY0()) + a0 * noy
        sigma = self.getSigma()
        gamma = self.getGamma()
        I = a0*noy*(1-sigma**2)*I1/u0 + u0*(sigma**2)*I2/(a0*noy) + 2*sigma*sqrt(1-(sigma**2))*I3 - gamma*(I1*I2 - I3**2)/a0
        return I
    def getDIntegral4(self, y, ny_prev, I, dy):
        a0 = self.soundSpeed(self.__atmosphere.getTemperature(self.__params.getY0()))
        a = self.soundSpeed(self.__atmosphere.getTemperature(y))
        density = self.__atmosphere.getDensity(y)
        ny = self.n_definition(y, ny_prev)['ny']
        nox = self.getNox()
        noy = self.getNoy()
        noz = self.getNoz()
        dWindX = self.__atmosphere.getWindX(y) - self.__atmosphere.getWindX(self.__params.getY0())
        dWindZ = self.__atmosphere.getWindZ(y) - self.__atmosphere.getWindZ(self.__params.getY0())
        # коэффициент u = Vy + a*ny
        u = self.__atmosphere.getWindY(y) + self.soundSpeed(self.__atmosphere.getTemperature(y)) * ny
        # коэффициент K = Un/Un0 = (a+Vy*ny)/a_star
        a_star = a0 + self.__atmosphere.getWindY(self.__params.getY0()) * noy - (dWindX * nox + dWindZ * noz)
        K = (a + self.__atmosphere.getWindY(y) * ny) / a_star
        dI4 = sqrt(a/(density*(K**3)*fabs((u**3)*I)))*fabs(dy)
        #dI4 = (a0*sqrt(a0))/(a**2 * sqrt(a*density*fabs(ny*I))) * dy/ny
        return dI4

    def getK(self, I4):
        kappa = self.__params.getKappa()
        M = self.getFlightMachNumber()
        sigma = self.getSigma()
        noy = self.getNoy()
        Y = self.__params.getLift()
        l = self.__params.getLength()
        a0 = self.soundSpeed(self.__atmosphere.getTemperature(self.__params.getY0()))
        density0 = self.__atmosphere.getDensity(self.__params.getY0())
        k = 0.5*(kappa+1) * M**1.5 * I4 * sqrt((sigma**2+(1-sigma**2)*noy**2)/(fabs(noy)*a0*density0))*Y/(pi*sqrt(2*l**5))
        #k = 0.5 * (kappa+1) * I4 * sqrt((sigma**2+(1-sigma**2)*noy**2)/(fabs(noy)*a0*density0))/(pi*sqrt(2))
        return k

    def getK1(self, y00, I, ny_prev):
        a00 = self.soundSpeed(self.__atmosphere.getTemperature(y00))
        density00 = self.__atmosphere.getDensity(y00)
        a0 = self.soundSpeed(self.__atmosphere.getTemperature(self.__params.getY0()))
        density0 = self.__atmosphere.getDensity(self.__params.getY0())
        sigma = self.getSigma()
        noy = self.getNoy()
        M = self.getFlightMachNumber()
        Y = self.__params.getLift()
        l = self.__params.getLength()
        ny = self.n_definition(y00, ny_prev)['ny']
        nox = self.getNox()
        noy = self.getNoy()
        noz = self.getNoz()
        dWindX = self.__atmosphere.getWindX(y00) - self.__atmosphere.getWindX(self.__params.getY0())
        dWindZ = self.__atmosphere.getWindZ(y00) - self.__atmosphere.getWindZ(self.__params.getY0())
        # коэффициент K = Un/Un0 = (a+Vy*ny)/a_star
        a_star = a0 + self.__atmosphere.getWindY(self.__params.getY0()) * noy - (dWindX * nox + dWindZ * noz)
        K = (a00 + self.__atmosphere.getWindY(y00) * ny) / a_star
        # коэффициент u = Vy + a*ny
        u00 = self.__atmosphere.getWindY(y00) + self.soundSpeed(self.__atmosphere.getTemperature(y00)) * ny
        K1 = (1/pi)*sqrt(M)*sqrt(((density00*a00**3)/(2*a0*density0*K))*((sigma**2+(1-sigma**2)*noy**2)/(fabs(noy*u00*I))))*(Y/sqrt(l**3))
        return K1

    def getK2(self):
        a0 = self.soundSpeed(self.__atmosphere.getTemperature(self.__params.getY0()))
        noy = self.getNoy()
        M = self.getFlightMachNumber()
        l = self.__params.getLength()
        a_star = a0 + self.__atmosphere.getWindY(self.__params.getY0()) * noy
        k2 = l/(M*a_star)
        return k2
    def getInitialIntagrals(self):
        a0 = self.soundSpeed(self.__atmosphere.getTemperature(self.__params.getY0()))
        T0 = self.__atmosphere.getTemperature(self.__params.getY0())
        nox = self.getNox()
        noy = self.getNoy()
        sigma = self.getSigma()
        density0 = self.__atmosphere.getDensity(self.__params.getY0())
        # коэффициент u = Vy + a*ny
        u0 = self.__atmosphere.getWindY(self.__params.getY0()) + a0 * noy
        dy = (self.__atmosphere.getWindY(self.__params.getY0()) + a0 * noy)*self.getDt()/10
        I1_0 = (a0**2)/(u0**2) * fabs(dy)
        I2_0 = ((a0**4)/(u0**3))*fabs(dy)
        I3_0 = 0
        I4_0 = 2*sqrt((dy*fabs(noy))/(density0 * u0 * a0**2 * (sigma**2+(1-sigma**2)*noy**2)))
        I0 = {}
        I0['I1'] = I1_0
        I0['I2'] = I2_0
        I0['I3'] = I3_0
        I0['I4'] = I4_0
        return I0






    def soundSpeed(self, T):
        return sqrt(self.__params.getKappa() * self.__params.getR() * T)

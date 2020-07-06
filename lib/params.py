from math import *

class Params:

    def __init__(self, y0=15000, x0=0, z0=0, omega=0, omega1=0, M=1.5):
        self.__y0 = y0
        self.__x0 = x0
        self.__z0 = z0
        self.__omega = omega
        self.__omega1 = omega1
        self.__M = M

    def getMachNumber(self):
        return self.__M

    def setMachNumber(self, M):
        self.__M = M

    def getOmega(self):
        return self.__omega

    def setOmega(self, omega):
        self.__omega = omega

    def getOmega1(self):
        return self.__omega1

    def setOmega1(self, omega1):
        self.__omega1 = omega1


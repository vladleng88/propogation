from math import *

class Raytrace:

    def __init__(self, tetta, params):
        self.__tetta = tetta
        self.__params = params
        self.__nox = self.__no_definition()['nox']
        self.__noy = self.__no_definition()['noy']
        self.__noz = self.__no_definition()['noz']
    def __no_definition(self):
        no = {}
        omega = self.__params.getOmega()
        omega1 = self.__params.getOmega1()
        myu = 1/self.__params.getMachNumber()
        nox = sin(radians(myu))*cos(radians(omega1))*cos(radians(omega)) -\
                   cos(radians(myu))*sin(radians(omega1))*sin(radians(self.__tetta)) + \
                   cos(radians(myu))*cos(radians(omega1))*sin(radians(omega))*cos(radians(self.__tetta))
        noy = sin(radians(myu))*sin(radians(omega)) - cos(radians(myu))*cos(radians(omega))*cos(radians(self.__tetta))
        noz = sin(radians(myu))*sin(radians(omega1))*cos(radians(omega)) + \
              cos(radians(myu))*cos(radians(omega1))*sin(radians(self.__tetta)) + \
              cos(radians(myu))*sin(radians(omega1))*sin(radians(omega))*cos(radians(self.__tetta))
        no['nox'] = nox
        no['noy'] = noy
        no['noz'] = noz
        return no
    def getNox(self):
        return self.__nox
    def getNoy(self):
        return self.__noy
    def getNoz(self):
        return self.__noz


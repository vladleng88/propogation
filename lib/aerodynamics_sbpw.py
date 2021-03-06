from lib.params import Params
from math import *

class AerodynamicsSBPW:

    def __init__(self, params: Params, MachNumber, path=r'in/sbpw_atm/case2_0.txt', Ngeoam = 4048):
        self.__path = path
        self.__withemRef = {}
        self.__dSdX = {}
        self.__ettaRef = {}
        self.__potentialRef = {}
        self.__Ngeom = Ngeoam
        self.__params = params
        self.__MachNumber = MachNumber


    def setAerodynmamicsData(self, k111):
        input = self.getInputData()
        inputResample = self.getResampleInputData(input)
        const = (self.__MachNumber**2*self.__params.getKappa())/(2*pi*self.__params.getLength()*sqrt(2*sqrt(self.__MachNumber**2 - 1)))
        self.__ettaRef = inputResample['ettaRefResample']
        print('const', const)
        for i in range(0, len(inputResample['dSdXResample'])):
            self.__dSdX[i] = inputResample['dSdXResample'][i]

        self.__potentialRef[0] = 0
        self.__withemRef[0] = 0

        betta = sqrt(self.__MachNumber ** 2 - 1)
        r = 82.296

        for i in range(1, len(self.__dSdX)):
            potential_integral = 0
            for j in range(0, i):
                if self.__ettaRef[j]>3:
                    continue
                potential_loc = pi*sqrt(2*betta*r)*0.5*(self.__dSdX[j]+self.__dSdX[j+1])*(self.__ettaRef[j+1] - self.__ettaRef[j])*self.__params.getLength()
                potential_integral = potential_integral + potential_loc
            self.__potentialRef[i] =k111*potential_integral

        for i in range(1, len(self.__dSdX)-1):
            withem_loc = (self.__potentialRef[i+1]-self.__potentialRef[i-1])/(self.__ettaRef[i+1] - self.__ettaRef[i-1])
            #print(withem_loc)
            #count += 1
            if self.__ettaRef[i] < 0:
                self.__withemRef[i] = 0
            else:
                self.__withemRef[i] = withem_loc
        self.__withemRef[len(self.__dSdX)-1] = 0
       # betta = sqrt(self.__MachNumber**2-1)
        #r = 82.296
        #for i in range(0, len(self.__dSdX)):
        #    self.__withemRef[i] = pi*sqrt(2*betta*r) * self.__dSdX[i]


    def getInputData(self):
        file = open(self.__path, 'r', encoding='utf-8')
        dSdX = {}
        ettaRef = {}

        # position = int(file.read().find('@'))
        # file.seek(position)
        for line in file:
            if line[len(line) - 2] == '@':
                break
        count = 1
        for line in file:
            list = line.split('\t')
            #print(list)
            #if float(list[0]) < 0.0:
            #    continue
            if float(list[0]) > 3:
                break
            ettaRef[count] = float(list[0])
            dSdX[count] = float(list[1])
            count += 1
        ettaRef[0] = ettaRef[1] - 0.5
        dSdX[0] = 0
        print('etta_before', len(ettaRef))
        print('dSdX_before', len(dSdX))
        return {'ettaRef': ettaRef, 'dSdX': dSdX}

    def getWithemRef(self):
        return self.__withemRef
    def getdSdX(self):
        return self.__dSdX

    def getEttaRef(self):
        return self.__ettaRef
    def getPotentalRef(self):
        return self.__potentialRef

    def getResampleInputData(self, input):
        ettaRef = input['ettaRef']
        dSdX = input['dSdX']
        Ka = {}
        Kb = {}
        for i in range(0, len(dSdX)-1):
            Ka[i] = (dSdX[i+1] - dSdX[i])/(ettaRef[i+1] - ettaRef[i])
            Kb[i] = (dSdX[i]*ettaRef[i+1] - dSdX[i+1]*ettaRef[i])/(ettaRef[i+1] - ettaRef[i])
        Xgeom = ettaRef[len(ettaRef)-1] - ettaRef[0]
        dx = (Xgeom)/(self.__Ngeom - 1)
        ettaRefResample = {}
        dSdXResample = {}
        #c=0
        for j in range(0, self.__Ngeom):
            ettaRefResample[j] = ettaRef[0] + j*dx
        for i in range(0, len(dSdX)-1):
            for j in range(0, self.__Ngeom):
                if (ettaRefResample[j]>=ettaRef[i] and ettaRefResample[j]<=ettaRef[i+1]):
                    dSdXResample[j] = Ka[i]*ettaRefResample[j] + Kb[i]
                    #print('count', c,'min', ettaRef[i], 'max',ettaRef[i+1], '\tetta=', ettaRefResample[j], '\t', dSdXResample[j])
                    #c+=1
                #if ettaRefResample[j] == ettaRef[i+1]:
                #    dSdXResample[j] = Ka[i] * ettaRefResample[j] + Kb[i]
                #    print('count', c, 'min', ettaRef[i], 'max', ettaRef[i + 1], '\tetta=', ettaRefResample[j], '\t',
                #          dSdXResample[j])
        #print(len(ettaRefResample))
        #print(len(dSdXResample))
        #for i in range(0, len(ettaRefResample)):
            #print('i', i, 'etta', ettaRefResample[i], 'ds_dx =', dSdXResample[i])

        #print('etta_after', len(ettaRefResample))
        #print('dSdX_after', len(dSdXResample))
        #print('etta_final', ettaRef[len(dSdX)-1])
        #print('etta_res_final', ettaRefResample[self.__Ngeom-1])
        return {'ettaRefResample': ettaRefResample, 'dSdXResample': dSdXResample}



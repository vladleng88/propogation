

class Atmosphere:
    __path = ''
    __H = {}; __sH = {}
    __T = {}; __sT = {}
    __Vx = {}; __sVx = {}
    __Vy = {}; __sVy = {}
    __Vz = {}; __sVz = {}
    __d = {}; __sd = {}

    __spline_d = 3.0
    __entries = 0

    # Коэффициенты сплайна
    __coeff_T = {}; __coeff_Vx = {}; __coeff_Vy = {}; __coeff_Vz = {}; __coeff_d = {}

    def __init__(self, path=r'in/air.dat'):
        self.__path = path

    def setAtmoshere(self):
        file = open(self.__path, 'r', encoding='utf-8')
        # position = int(file.read().find('@'))
        # file.seek(position)
        for line in file:
            if line[len(line)-2] == '@':
                break
        count = 1
        for line in file:
            list = line.split()
            self.__H[count] = float(list[0])
            self.__T[count] = float(list[1])
            self.__Vx[count] = float(list[2])
            self.__Vy[count] = float(list[3])
            self.__Vz[count] = float(list[4])
            self.__d[count] = float(list[5])
            count += 1
            self.__entries += 1
        self.__setSTemperature()
        self.__setSDensity()
        self.__setSWindX()
        self.__setSWindY()
        self.__setSWindZ()
        self.__setSHeight()

    def __setSTemperature(self):
        spline_coeff = self.spline_fit(self.__H, self.__T)
        self.__sT = spline_coeff['spline']
        spline_coeff.pop('spline')
        self.__coeff_T = spline_coeff
    def getSTemperature(self):
        return self.__sT
    def getSTemperatureCoeff(self):
        return self.__coeff_T

    def __setSDensity(self):
        spline_coeff = self.spline_fit(self.__H, self.__d)
        self.__sd = spline_coeff['spline']
        spline_coeff.pop('spline')
        self.__coeff_d = spline_coeff
    def getSDensity(self):
        return self.__sd
    def getSDensityCoeff(self):
        return self.__coeff_d

    def __setSWindX(self):
        spline_coeff = self.spline_fit(self.__H, self.__Vx)
        self.__sVx = spline_coeff['spline']
        spline_coeff.pop('spline')
        self.__coeff_Vx = spline_coeff
    def getSWindX(self):
        return self.__sVx
    def getSWindXCoeff(self):
        return self.__coeff_Vx

    def __setSWindY(self):
        spline_coeff = self.spline_fit(self.__H, self.__Vy)
        self.__sVy = spline_coeff['spline']
        spline_coeff.pop('spline')
        self.__coeff_Vy = spline_coeff
    def getSWindY(self):
        return self.__sVy
    def getSWindYCoeff(self):
        return self.__coeff_Vy

    def __setSWindZ(self):
        spline_coeff = self.spline_fit(self.__H, self.__Vz)
        self.__sVz = spline_coeff['spline']
        spline_coeff.pop('spline')
        self.__coeff_Vz = spline_coeff
    def getSWindZ(self):
        return self.__sVz
    def getSWindZCoeff(self):
        return self.__coeff_Vz

    def __setSHeight(self):
        self.__sH = self.__spline_height(self.__H)
    def getSHeight(self):
        return self.__sH

    def __spline_height(self, h):
        j = 1
        spline = {}
        dhl = {}
        dhr = {}
        h_i0 = {}
        h_i1 = {}
        for i in range(2, len(h)):
            dh_l = (h[i]-h[i-1])/self.__spline_d
            dh_r = (h[i+1]-h[i])/self.__spline_d
            if (dh_l < dh_r):
                dh_r = dh_l
            else:
                dh_l = dh_r
            dhl[i] = dh_l
            dhr[i] = dh_r
            dh_il = h[i]-dhl[i]
            dh_ir = h[i]+dhr[i]
            h_i0[i] = dh_il
            h_i1[i] = dh_ir
        # h_i0[self.__entries] = h[self.__entries]
        spline[j] = h[1]
        j += 1
        for i in range(2, len(h)):
            spline[j] = h_i0[i]
            j += 1
            spline[j] = h[i]
            j += 1
            spline[j] = h_i1[i]
            j += 1
        spline[j] = h[self.__entries]
        return spline

    def spline_fit(self, h, y):
        j = 1
        spline = {}
        coeff_a0 = {}
        coeff_a1 = {}
        coeff_a2 = {}
        coeff_a3 = {}
        coeff_b = {}

        dhl = {}
        dhr = {}
        h_i0 = {}
        y_i0 = {}
        h_i1 = {}
        y_i1 = {}
        for i in range(2, len(h)):
            dh_l = (h[i]-h[i-1])/self.__spline_d
            dh_r = (h[i+1]-h[i])/self.__spline_d
            if (dh_l < dh_r):
                dh_r = dh_l
            else:
                dh_l = dh_r
            dhl[i] = dh_l
            dhr[i] = dh_r
            dh_il = h[i]-dhl[i]
            dh_ir = h[i]+dhr[i]
            h_i0[i] = dh_il
            h_i1[i] = dh_ir
            y_i0[i] = y[i]-dhl[i]*(y[i]-y[i-1])/(h[i]-h[i-1])
            y_i1[i] = y[i]+dhr[i]*(y[i+1]-y[i])/(h[i+1]-h[i])

        # h_i0[self.__entries] = h[self.__entries]

        # define coefficient for Rergion 0
        coeff_a0[j] = y[1]
        coeff_a1[j] = (y[2]-y[1])/(h[2]-h[1])
        coeff_a2[j] = 0
        coeff_a3[j] = 0
        coeff_b[j] = h[1]
        spline[j] = y[1]
        # spline.setdefault(h[1], y[1])
        j+=1
        for i in range(2, len(h)):
            # define coefficient for Rergion 1
            # region = 1
            coeff_a0[j] = y_i0[i]
            coeff_a1[j] = (y[i]-y[i-1])/(h[i]-h[i-1])
            coeff_a2[j] = 0
            coeff_a3[j] = (y_i0[i]-2*y[i]+y_i1[i])/(6*(dhl[i])**3)
            coeff_b[j] = h_i0[i]
            spline[j] = y_i0[i]
            #spline[j] = coeff_a0[j] + coeff_a1[j]*(hi-h_i0[i])+coeff_a2[j]*(hi-h_i0[i])**2 + coeff_a3[j]*(hi-h_i0[i])**3
            j += 1
            # h_y1 = self.__increment(i, a0_1[i], a1_1[i], a2_1[i], a3_1[i], b1, region, h_i1[i])
            # for key in h_y1:
            #    spline.setdefault(key, h_y1[key])
            # define coefficient for Rergion 2
            # region = 2
            coeff_a0[j] = (y_i0[i]+4*y[i]+y_i1[i])/6
            coeff_a1[j] = (y_i1[i] - y_i0[i]) / (2*dhr[i])
            coeff_a2[j] = (y_i0[i]-2*y[i]+y_i1[i])/(2*(dhr[i])**2)
            coeff_a3[j] = (-1)*(y_i0[i] - 2*y[i] + y_i1[i]) / (6*(dhr[i])**3)
            coeff_b[j] = h[i]
            spline[j] = y[i]
            j += 1
            # h_y2 = self.__increment(i, a0_2[i], a1_2[i], a2_2[i], a3_2[i], b2, region, h_i1[i])
            # for key in h_y2:
            #    spline.setdefault(key, h_y2[key])

            # define coefficient for Rergion 3
            coeff_a0[j] = y_i1[i]
            coeff_a1[j] = (y[i+1] - y[i]) / (h[i+1] - h[i])
            coeff_a2[j] = 0
            coeff_a3[j] = 0
            coeff_b[j] = h_i1[i]
            spline[j] = y_i1[i]
            j += 1
        # spline.setdefault(h[self.__entries], y[self.__entries])
        spline[j] = y[self.__entries]
        spline_return = {}
        spline_return['spline'] = spline
        spline_return['a0'] = coeff_a0
        spline_return['a1'] = coeff_a1
        spline_return['a2'] = coeff_a2
        spline_return['a3'] = coeff_a3
        spline_return['b'] = coeff_b
        return spline_return

    def getTemperature(self, h):
        if h < 0 or h > self.__getHmax():
            exit('Error: Измерения за гранью температурного профиля')
        for key in self.__sH:
            if h >= self.__sH[key] and h < self.__sH[key+1]:
                b = self.__sH[key]
                a0 = self.__coeff_T['a0'][key]
                a1 = self.__coeff_T['a1'][key]
                a2 = self.__coeff_T['a2'][key]
                a3 = self.__coeff_T['a3'][key]
                yi = a0 + a1*(h - b) + a2*(h - b)**2 + a3*(h - b)**3
        return yi

    def getWindX(self, h):
        if h < 0 or h > self.__getHmax():
            exit('Error: Измерения за гранью температурного профиля')
        for key in self.__sH:
            if h >= self.__sH[key] and h < self.__sH[key+1]:
                b = self.__sH[key]
                a0 = self.__coeff_Vx['a0'][key]
                a1 = self.__coeff_Vx['a1'][key]
                a2 = self.__coeff_Vx['a2'][key]
                a3 = self.__coeff_Vx['a3'][key]
                yi = a0 + a1 * (h - b) + a2 * (h - b) ** 2 + a3 * (h - b) ** 3
        return yi

    def getWindY(self, h):
        if h < 0 or h > self.__getHmax():
            exit('Error: Измерения за гранью температурного профиля')
        for key in self.__sH:
            if h >= self.__sH[key] and h < self.__sH[key+1]:
                b = self.__sH[key]
                a0 = self.__coeff_Vy['a0'][key]
                a1 = self.__coeff_Vy['a1'][key]
                a2 = self.__coeff_Vy['a2'][key]
                a3 = self.__coeff_Vy['a3'][key]
                yi = a0 + a1 * (h - b) + a2 * (h - b) ** 2 + a3 * (h - b) ** 3
        return yi
    def getWindZ(self, h):
        if h < 0 or h > self.__getHmax():
            exit('Error: Измерения за гранью температурного профиля')
        for key in self.__sH:
            if h >= self.__sH[key] and h < self.__sH[key+1]:
                b = self.__sH[key]
                a0 = self.__coeff_Vz['a0'][key]
                a1 = self.__coeff_Vz['a1'][key]
                a2 = self.__coeff_Vz['a2'][key]
                a3 = self.__coeff_Vz['a3'][key]
                yi = a0 + a1 * (h - b) + a2 * (h - b) ** 2 + a3 * (h - b) ** 3
        return yi

    def getDensity(self, h):
        if h < 0 or h > self.__getHmax():
            exit('Error: Измерения за гранью температурного профиля')
        for key in self.__sH:
            if h >= self.__sH[key] and h < self.__sH[key+1]:
                b = self.__sH[key]
                a0 = self.__coeff_d['a0'][key]
                a1 = self.__coeff_d['a1'][key]
                a2 = self.__coeff_d['a2'][key]
                a3 = self.__coeff_d['a3'][key]
                yi = a0 + a1 * (h - b) + a2 * (h - b) ** 2 + a3 * (h - b) ** 3
        return yi

    def __increment(self, i, a0, a1, a2, a3, b, region, h_i1):
        FINENESS = 5.0
        if (region == 1):
            end = self.__H[i]
        else:
            end = h_i1
        inc = (end - b)/FINENESS
        h_y = {}
        for ii in range(0, int(FINENESS)+1, 1):
            hi = b + ii*inc
            yi = a0 + a1*(hi-b)+a2*(hi-b)**2 + a3*(hi-b)**3
            h_y[hi] = yi
        return h_y

    def getDefaultHeight(self):
        return self.__H
    def getDefaultTemperature(self):
        return self.__T
    def getDefaultWindX(self):
        return self.__Vx
    def getDefaultWindY(self):
        return self.__Vy
    def getDefaultWindZ(self):
        return self.__Vz
    def getDefaultDensity(self):
        return self.__d
    def __getHmax(self):
        max = 0
        for key in self.__H:
            if max < self.__H[key]:
                max = self.__H[key]
        return max


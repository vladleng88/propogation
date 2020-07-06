

class Atmosphere:
    __path = ''
    __H = {}
    __T = {}
    __Vx = {}
    __Vy = {}
    __d = {}

    __spline_d = 3.0
    __entries = 0

    __sH = []
    __sT = {}
    __sVx = {}
    __sVy = {}
    __sd = {}

    def __init__(self, path=r'in/air.dat', T={}, H={}, Vx={}, Vy={}, d={}):
        self.__path = path
        self.__T = T
        self.__H = H
        self.__Vx = Vx
        self.__Vy = Vy
        self.__d = d

    def setAtmoshere(self):
        file = open(self.__path, 'r', encoding='utf-8')
        # position = int(file.read().find('@'))
        # print(position)
        # file.seek(position)
        print(self.__path)
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
            self.__d[count] = float(list[4])
            count += 1
            self.__entries += 1
            # print(list)
        # self.__setSTemperature()
        self.__setSDensity()
        # self.__setSWindX()
        # self.__setSWindY()

    def __setSTemperature(self):
        self.__sT = self.spline_fit(self.__H, self.__T)

    def getSTemperature(self):
        return self.__T

    def __setSDensity(self):
        self.__sd = self.spline_fit(self.__H, self.__d)

    def getSDensity(self):
        return self.__sd

    def __setSWindX(self):
        self.__sVx = self.spline_fit(self.__H, self.__Vx)

    def __setSWindY(self):
        self.__sVy = self.spline_fit(self.__H, self.__Vy)

    def setSplineHeight(self):
        return 1

    def spline_fit(self, h, y):
        spline = {}
        dhl = {}
        dhr = {}
        h_i0 = {}
        y_i0 = {}
        h_i1 = {}
        y_i1 = {}
        print('-----')
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

        h_i0[self.__entries] = h[self.__entries]
        a0_1 = {}; a0_2 = {}
        a1_1 = {}; a1_2 = {}
        a2_1 = {}; a2_2 = {}
        a3_1 = {}; a3_2 = {}
        # coefficient for Rergion 0
        a0_0 = y[1]
        a1_0 = (y[2]-y[1])/(h[2]-h[1])
        a2_0 = 0
        a3_0 = 0
        spline.setdefault(h[1], y[1])
        for i in range(2, len(h)):
            # coefficient for Rergion 1
            region = 1
            a0_1[i] = y_i0[i]
            a1_1[i] = (y[i]-y[i-1])/(h[i]-h[i-1])
            a2_1[i] = 0
            a3_1[i] = (y_i0[i]-2*y[i]+y_i1[i])/(6*(dhl[i])**3)
            b1 = h_i0[i]
            h_y1 = self.__increment(i, a0_1[i], a1_1[i], a2_1[i], a3_1[i], b1, region, h_i1[i])
            for key in h_y1:
                spline.setdefault(key, h_y1[key])

            # coefficient for Rergion 2
            region = 2
            a0_2[i] = (y_i0[i]+4*y[i]+y_i1[i])/6
            a1_2[i] = (y_i1[i] - y_i0[i]) / (2*dhr[i])
            a2_2[i] = (y_i0[i]-2*y[i]+y_i1[i])/(2*(dhr[i])**2)
            a3_2[i] = (-1)*(y_i0[i] - 2*y[i] + y_i1[i]) / (6*(dhr[i])**3)
            b2 = h[i]
            h_y2 = self.__increment(i, a0_2[i], a1_2[i], a2_2[i], a3_2[i], b2, region, h_i1[i])
            for key in h_y2:
                spline.setdefault(key, h_y2[key])
            # coefficient for Rergion 3
            a0_3 = y_i1[i]
            a1_3 = (y[2] - y[1]) / (h[2] - h[1])
            a2_3 = 0
            a3_3 = 0
        spline.setdefault(h[self.__entries], y[self.__entries])
        return spline
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

    def getHeight(self):
        return self.__H
    def getWindX(self):
        return self.__Vx
    def getWindY(self):
        return self.__Vy
    def getDensity(self):
        return self.__d


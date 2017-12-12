# Math functions and constants
from math import sqrt, sin, cos, acos, pi, tan
# Interval analysis
from interval import interval, imath
# Operator
import operator
# Min Float
import sys
# GlobalOptimization by using Brute Force
from scipy import optimize

# This script finds the extrema for a box by using interval analysis
# (getPhi function)

class Phiwos(object):
    def __init__(self, rhomin, rhomax, thetamin, thetamax):
        # Initialize the parameters of the planar parallel robot
        self.__rho = [rhomin, rhomax]
        self.__theta_min = thetamin
        self.__theta_max = thetamax

    def getRho(self):
        return self.__rho

    def getTheta_Bounds(self):
        return (self.__theta_min, self.__theta_max)

    def getCi(self):
        return [0,0]

############################################################################
#                               Interval Analysis
############################################################################
    ########################################################################
    #                             The Link Constraints
    ########################################################################
    def gRhoi(self, x):
        return [x[0]**2 + x[1]**2 - self.__rho[1]**2,\
               self.__rho[0]**2 - x[0]**2 - x[1]**2]

    ########################################################################
    #                            The Angle Constraints
    ########################################################################
    @staticmethod
    @interval.function
    def arccos(x):
        if x[1] < -1.0 or x[0] > 1.0:
            print 'Invalid argument in acos function. The argument is out of this interval[-1, 1]'
        l = -1.0 if x[0] < -1.0 else x[0]
        r = 1.0 if x[1] > 1.0 else x[1]
        return interval([acos(r), acos(l)])

    def gThetaPlus(self, invec):
        x = invec[0]
        y = invec[1]
        return [self.arccos(x/(imath.sqrt(x**2+y**2))) - self.__theta_max,\
                self.__theta_min - self.arccos(x/(imath.sqrt(x**2+y**2)))]

    def gThetaMinus(self, invec):
        x = invec[0]
        y = invec[1]
        return [2*pi - self.arccos(x/(imath.sqrt(x**2+y**2))) - self.__theta_max,\
                self.__theta_min - 2*pi + self.arccos(x/(imath.sqrt(x**2+y**2)))]

    def gTheta(self, invec):
        x = invec[0]
        y = invec[1]
        if y[0][0] >= 0.0 and y[0][1] > 0.0:
            return self.gThetaPlus(invec)
        if y[0][0] < 0.0 and y[0][1] <= 0.0:
            return self.gThetaMinus(invec)
        if y[0][0] < 0.0 and y[0][1] > 0.0:
            raise ValueError

    ########################################################################
    #                           Extrema Finding
    ########################################################################

    def getIntlRes(self, iBox):
        vec = iBox.getInterval()
        flist = ['gRhoi', 'gTheta']

        x = vec[0]
        y = vec[1]
        if y[0][0] < 0.0 and y[0][1] > 0.0:
            leq0 = [x, interval([y[0][0], 0.0])]
            geq0 = [x, interval([0.0, y[0][1]])]

            resleq0 = [interval(intls) for f in flist for intls in getattr(self, f)(leq0)]
            resgeq0 = [interval(intls) for f in flist for intls in getattr(self, f)(geq0)]
            return resleq0, resgeq0
        else:
            return [interval(intls) for f in flist for intls in getattr(self, f)(vec)]

    def getPhi(self, iBox):
        res = self.getIntlRes(iBox)

        if len(res) == 2:
            phileq0 = interval[max(x[0][0] for x in res[0]),\
                               max(x[0][1] for x in res[0])]
            phigeq0 = interval[max(x[0][0] for x in res[1]),\
                               max(x[0][1] for x in res[1])]

            return min(phileq0[0][0],phigeq0[0][0]), max(phileq0[0][1],phigeq0[0][1])
        else:
            phi = interval[max(x[0][0] for x in res),\
                           max(x[0][1] for x in res)]

            return phi[0][0], phi[0][1]

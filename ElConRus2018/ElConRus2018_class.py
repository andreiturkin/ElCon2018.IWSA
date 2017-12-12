import numpy as np
from numpy import linalg as la
import datetime
# Pi
from math import pi
# Interval analysis
from interval import interval, imath
#Plotting
import matplotlib.pyplot as plt
import matplotlib.patches as patches
#Saving data to files
import os
import json

# from ElConRus2018_Phi import Phi
from ElConRus2018_Phiwos import Phiwos
from ElConRus2018_Phiws import Phiws

# Additional Dependences
# Utils folder
import sys
from os.path import join as join_path
abspath = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.abspath(join_path(abspath, '..', 'Utils')))
from ElConRus2018Utils_NUC import Box, CoveringTree
from ElConRus2018Utils_Plotting import PlottingTree


class IAwoS(CoveringTree, PlottingTree, Phiwos):
############################################################################################
# Public Methods
############################################################################################
    def __init__(self, idelta, rhomin, rhomax, thetamin, thetamax, ShowCovPrc=False):

        # Define the Initial Rectangle P
        corner = 0.0
        side = 40.0

        iBox = Box((-side/2.0, -side/2.0,), (side, side))
        pBox = Box((-side/2.0, -side/2.0,), (side, side))

        # Call the following initialization routines
        CoveringTree.__init__(self, iBox, idelta, ShowCovPrc, 0.0)
        PlottingTree.__init__(self, pBox, ShowCovPrc)
        Phiwos.__init__(self, rhomin, rhomax, thetamin, thetamax)

    def SaveResults(self, fileName, AddRings):
        super(IAwoS, self).saveResultAsImage(self.getTree()['Tree'], fileName, AddRings=AddRings)

    def SaveSolution(self, fName):
        super(IAwoS, self).SaveSolution(fName)

    def isFileExist(self, fName):
        return super(IAwoS, self).isFileExist(fName)

    def LoadSolution(self, fName):
        super(IAwoS, self).LoadSolution(fName)

############################################################################################
# Private Methods
############################################################################################

############################################################################################
# Abstract Methods
############################################################################################

    ########################################################################################
    # CoveringTree
    ########################################################################################
    def getMinMaxVal(self, iBox):

        phiMin, phiMax = self.getPhi(iBox)
        # return minimum and maximum values
        return phiMin, phiMax

    ########################################################################################
    # PlottingTree
    ########################################################################################
    def AdditionalPlotting(self, ax):
        ci = self.getCi()

        for c in ci:
            for r in self.getRho():
                self.drawCircle(c, r)

        theta_bounds = self.getTheta_Bounds()
        theta_intl = interval[theta_bounds[0], theta_bounds[1]]

        # theta_a = [a & b for a, b in zip(theta_a_intl,theta_b_intl)]
        # theta_a = [theta_a_intli & theta_b_intli for theta_a_intli, theta_b_intli in theta_a_intl, theta_b_intl]

        psi = (theta_intl[0][0], theta_intl[0][1])

        for idx, c in enumerate(ci):
            for r in self.getRho():
                self.drawArc(c, r, psi[idx])

        for idx, c in enumerate(ci):
            self.drawLine(c, self.getRho(), psi[idx])
        return

############################################################################################
# !Abstract Methods
############################################################################################

class IAwS(CoveringTree, PlottingTree, Phiws):
############################################################################################
# Public Methods
############################################################################################
    def __init__(self, idelta, rhomin, rhomax, thetamin, thetamax, ShowCovPrc=False):

        # Define the Initial Rectangle P
        corner = 0.0
        side = 40.0

        iBox = Box((-side/2.0, -side/2.0,), (side, side))
        pBox = Box((-side/2.0, -side/2.0,), (side, side))

        # Call the following initialization routines
        CoveringTree.__init__(self, iBox, idelta, ShowCovPrc, 0.0)
        PlottingTree.__init__(self, pBox, ShowCovPrc)
        Phiws.__init__(self, rhomin, rhomax, thetamin, thetamax)

    def SaveResults(self, fileName, AddRings):
        super(IAwS, self).saveResultAsImage(self.getTree()['Tree'], fileName, AddRings=AddRings)

    def SaveSolution(self, fName):
        super(IAwS, self).SaveSolution(fName)

    def isFileExist(self, fName):
        return super(IAwS, self).isFileExist(fName)

    def LoadSolution(self, fName):
        super(IAwS, self).LoadSolution(fName)

############################################################################################
# Private Methods
############################################################################################

############################################################################################
# Abstract Methods
############################################################################################

    ########################################################################################
    # CoveringTree
    ########################################################################################
    def getMinMaxVal(self, iBox):

        phiMin, phiMax = self.getPhi(iBox)
        # return minimum and maximum values
        return phiMin, phiMax

    ########################################################################################
    # PlottingTree
    ########################################################################################
    def AdditionalPlotting(self, ax):
        ci = self.getCi()

        for c in ci:
            for r in self.getRho():
                self.drawCircle(c, r)

        theta_bounds = self.getTheta_Bounds()
        theta_intl = interval[theta_bounds[0], theta_bounds[1]]

        # theta_a = [a & b for a, b in zip(theta_a_intl,theta_b_intl)]
        # theta_a = [theta_a_intli & theta_b_intli for theta_a_intli, theta_b_intli in theta_a_intl, theta_b_intl]

        psi = (theta_intl[0][0], theta_intl[0][1])

        for idx, c in enumerate(ci):
            for r in self.getRho():
                self.drawArc(c, r, psi[idx])

        for idx, c in enumerate(ci):
            self.drawLine(c, self.getRho(), psi[idx])
        return

############################################################################################
# !Abstract Methods
############################################################################################

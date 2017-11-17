
import timeit
import datetime
import matplotlib
# http://matplotlib.org/faq/usage_faq.html#what-is-a-backend
matplotlib.use('TkAgg')

import numpy as np
import os, glob, gc, sys
from os.path import join as join_path
from math import pi

# A Parallel Robot
abspath = os.path.dirname(os.path.abspath(__file__))
sys.path.append(join_path(abspath, 'ElConRus2018'))
from ElConRus2018_class import SimpleRobot

def SR_GetWorkspace(iDelta, rhomin, rhomax, thetamin, thetamax, fig, ShowRes=False):
    # Define a system of inequalities
    System = SimpleRobot(iDelta, rhomin, rhomax, thetamin, thetamax, ShowCovPrc=False)
    maxLevels = 64
    System.getSolution(maxLevels)

    if ShowRes:
        System.SaveResults('./Images/ElConRus2018_Figure_{}.pdf'.format(fig), AddRings=False)

def CreateOrClear(DirName, bClean):
    if not os.path.isdir(DirName):
        os.makedirs(DirName)
    else:
        if bClean:
            for f in glob.glob(DirName + '*.*'):
                os.remove(f)

if __name__ == '__main__':
    ImageDir = 'Images/'
    # Cleanning up the directories
    CreateOrClear(ImageDir, False)
    print'\n#############################################################################'
    print'                            Parallel Robot                                     '
    print'#############################################################################\n'
    print'                              Figure 1a                                        '
    # Get the workspace
    SR_GetWorkspace(0.05, rhomin = 8.0, rhomax = 12.0, \
                         thetamin = pi/16.0, thetamax = pi - pi/16.0, \
                         fig = '1a', ShowRes=True)
    print'                              Figure 1b                                        '
    # Get the workspace
    SR_GetWorkspace(0.05, rhomin = 8.0, rhomax = 12.0, \
                          thetamin = 3*pi/4.0, thetamax = 2*pi-pi/4.0, \
                          fig = '1b', ShowRes=True)
    print'                              Figure 1c                                        '
    # Get the workspace
    SR_GetWorkspace(0.05, rhomin = 8.0, rhomax = 12.0, \
                          thetamin = pi + pi/16.0, thetamax = 2*pi - pi/16.0, \
                          fig = '1c', ShowRes=True)
    print'                              Figure 1d                                        '
    # Get the workspace
    SR_GetWorkspace(0.05, rhomin = 8.0, rhomax = 12.0, \
                          thetamin = 0.0, thetamax = 2*pi, \
                          fig = '1d', ShowRes=True)
    print'\n#############################################################################'
    print'                                 Done!                                         '
    print'#############################################################################\n'
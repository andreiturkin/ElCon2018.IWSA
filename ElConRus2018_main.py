
import timeit
import datetime
import matplotlib
# http://matplotlib.org/faq/usage_faq.html#what-is-a-backend
matplotlib.use('TkAgg')
import progressbar
from terminaltables import SingleTable

import numpy as np
import os, glob, gc, sys
from os.path import join as join_path
from math import pi

# Systems
abspath = os.path.dirname(os.path.abspath(__file__))
sys.path.append(join_path(abspath, 'ElConRus2018'))
from ElConRus2018_class import IAwS
from ElConRus2018_class import IAwoS

def GetSolution_woS(iDelta, rhomin, rhomax, thetamin, thetamax, fig, ShowRes=False):
    # Define a system of inequalities that is approached by using interval analysis
    # without simplification
    System = IAwoS(iDelta, rhomin, rhomax, thetamin, thetamax, ShowCovPrc=False)
    maxLevels = 64
    System.getSolution(maxLevels)

    if ShowRes:
        System.SaveResults('./Images/ElConRus2018_Figure_{}.pdf'.format(fig), AddRings=False)

    return ['{}'.format(iDelta), \
            '{}'.format(System.getResIterations())]

def GetSolution_wS(iDelta, rhomin, rhomax, thetamin, thetamax, fig, ShowRes=False):
    # Define a system of inequalities that is approached by using interval analysis
    # with simplification
    System = IAwS(iDelta, rhomin, rhomax, thetamin, thetamax, ShowCovPrc=False)
    maxLevels = 64
    System.getSolution(maxLevels)

    if ShowRes:
        System.SaveResults('./Images/ElConRus2018_Figure_{}.pdf'.format(fig), AddRings=False)

    return ['{}'.format(iDelta), \
            '{}'.format(System.getResIterations())]

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
    GetSolution_woS(0.05, rhomin = 8.0, rhomax = 12.0, \
                         thetamin = pi/16.0, thetamax = pi - pi/16.0, \
                         fig = '1a', ShowRes=True)
    print'                              Figure 1b                                        '
    # Get the workspace
    GetSolution_woS(0.05, rhomin = 8.0, rhomax = 12.0, \
                          thetamin = 3*pi/4.0, thetamax = 2*pi-pi/4.0, \
                          fig = '1b', ShowRes=True)
    print'                              Figure 1c                                        '
    # Get the workspace
    GetSolution_woS(0.05, rhomin = 8.0, rhomax = 12.0, \
                          thetamin = pi + pi/16.0, thetamax = 2*pi - pi/16.0, \
                          fig = '1c', ShowRes=True)
    print'                              Figure 1d                                        '
    # Get the workspace
    GetSolution_woS(0.05, rhomin = 8.0, rhomax = 12.0, \
                          thetamin = 0.0, thetamax = 2*pi, \
                          fig = '1d', ShowRes=True)
    print'#############################################################################\n'
    print'                              Figure 3a                                        '
    # Get the workspace
    GetSolution_wS(0.05, rhomin = 8.0, rhomax = 12.0, \
                         thetamin = pi/16.0, thetamax = pi - pi/16.0, \
                         fig = '3a', ShowRes=True)
    print'                              Figure 3b                                        '
    # Get the workspace
    GetSolution_wS(0.05, rhomin = 8.0, rhomax = 12.0, \
                          thetamin = 3*pi/4.0, thetamax = 2*pi-pi/4.0, \
                          fig = '3b', ShowRes=True)
    print'                              Figure 3c                                        '
    # Get the workspace
    GetSolution_wS(0.05, rhomin = 8.0, rhomax = 12.0, \
                          thetamin = pi + pi/16.0, thetamax = 2*pi - pi/16.0, \
                          fig = '3c', ShowRes=True)
    print'                              Figure 3d                                        '
    # Get the workspace
    GetSolution_wS(0.05, rhomin = 8.0, rhomax = 12.0, \
                          thetamin = 0.0, thetamax = 2*pi, \
                          fig = '3d', ShowRes=True)
    print'###############################################################################'
    print'                             Table 1'
    Table1_deltas = [0.5, 0.3, 0.25, 0.2, 0.14, 0.08, 0.07, 0.05, 0.03]
    pbar = progressbar.ProgressBar(maxval=len(Example2_deltas), \
                                  widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
    pbar.start()
    rows = []
    rows.append(['     ', 'An approach without simplification','The proposed approach with simplified system'])
    rows.append(['Delta', 'Number of Iterations'              ,'Number of Iterations'])
    for i, d in enumerate(Table1_deltas):
        wos = GetSolution_woS( d, rhomin = 8.0, rhomax = 12.0, \
                               thetamin = pi/16.0, thetamax = pi - pi/16.0, \
                               fig = None, ShowRes=False)
        ws = GetSolution_wS( d, rhomin = 8.0, rhomax = 12.0, \
                             thetamin = pi/16.0, thetamax = pi - pi/16.0, \
                             fig = None, ShowRes=False)
        rows.append([ wos[0], wos[1], ws[1] ])
        pbar.update(i+1)
    pbar.finish()
    table_instance = SingleTable(rows)
    table_instance.inner_heading_row_border = False
    print(table_instance.table)
    print'\n#############################################################################'
    print'                                 Done!                                         '
    print'#############################################################################\n'
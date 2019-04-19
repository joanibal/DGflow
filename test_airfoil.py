""" main script for proj4 """
import numpy as np
import numpy
import matplotlib.pyplot as plt
from meshes import Mesh

from cfdsolvers import DGSolver
import pickle
import copy


def naca0012(x, y):
    if y > 0:
        newY = 0.6*(0.2969*np.sqrt(x) - 0.1260*x - 0.3516*x**2 + 0.2843*x**3 - 0.1036*x**4)
    #     newY =
    elif y <= 0:
        newY = -0.6*(0.2969*np.sqrt(x) - 0.1260*x - 0.3516*x**2 + 0.2843*x**3 - 0.1036*x**4)

    return newY


airfoil = Mesh('meshes/naca0012.gri', wallGeomFunc=naca0012)
airfoil.refine()
alpha = 5.0
dalpha = 1e-4
CFDSolver = DGSolver(airfoil, order=2, alpha=alpha)
CFDSolver.solve()
cl = CFDSolver.postprocess()
CFDSolver.writeSolution('airfoil')

# dFdX_adjoint = CFDSolver.solveAdjoint()
# print(dFdX_adjoint)

# CFDSolver2 = DGSolver(airfoil, order=0,alpha=alpha+dalpha)
# CFDSolver2.solve()
# cl2 = CFDSolver2.postprocess()
# dFdX_FD = (cl2-cl)/dalpha
# print(dFdX_FD)

""" main script for proj4 """
import numpy as np
import matplotlib.pyplot as plt
from meshes import Mesh

from cfdsolvers import DGSolver
import pickle
import copy

airfoil = Mesh('meshes/naca0012.gri')
alpha = 5.0
dalpha = 1e-4
CFDSolver = DGSolver(airfoil, order=0,alpha=alpha)
CFDSolver.solve()
cl = CFDSolver.postprocess()
dFdX_adjoint = CFDSolver.solveAdjoint()

CFDSolver2 = DGSolver(airfoil, order=0,alpha=alpha+dalpha)
CFDSolver2.solve()
cl2 = CFDSolver2.postprocess()
dFdX_FD = (cl2-cl)/dalpha
print(dFdX_FD)
print(dFdX_adjoint)
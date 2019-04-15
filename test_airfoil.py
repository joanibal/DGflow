""" main script for proj4 """
import numpy as np
import matplotlib.pyplot as plt
from meshes import Mesh

from cfdsolvers import DGSolver
import pickle
import copy

def get_dFdx(CFDSolver):
    dalpha = 1e-4
    cl = CFDSolver.postprocess()
    CFDSolver2 = copy.deepcopy(CFDSolver)
    CFDSolver2.alpha = CFDSolver.alpha + dalpha
    cl2 = CFDSolver2.postprocess()
    dCldAlpha = (cl2-cl)/dalpha
    alpha = CFDSolver.alpha
    CF = CFDSolver.F/(CFDSolver.gamma/2*CFDSolver.P_inf*CFDSolver.mach_Inf**2)
    analytic = -CF[1]*np.sin(alpha)-CF[0]*np.cos(alpha)
    return analytic

def get_dFdU(CFDSolver):
    step = 1e-4
    U = CFDSolver.U
    U_shape = U.shape
    U = U.flatten()
    cl0 = CFDSolver.postprocess()
    dFdU = np.zeros_like(U)
    for iu in xrange(U.size): # TODO: do this sparsely, by perturbing only boundary elements
        Up = U.copy()
        dU = np.maximum(step,step*Up[iu])
        Up[iu] += dU
        CFDSolver.U = Up.reshape(U_shape)
        dFdU[iu] = (CFDSolver.postprocess() - cl0)/dU


airfoil = Mesh('meshes/naca0012.gri')
alpha = 5.0
CFDSolver = DGSolver(airfoil, order=0,alpha=alpha)
CFDSolver.solve()
dFdx = get_dFdx(CFDSolver)
print(dFdx)

dFdU = get_dFdU(CFDSolver)
print(dFdU)
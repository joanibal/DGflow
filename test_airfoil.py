""" main script for proj4 """
import numpy as np
import numpy
import matplotlib.pyplot as plt
from meshes import Mesh

from cfdsolvers import DGSolver
import pickle
import copy

def unpickleFile(fname):
    with open(fname, 'rb') as handle:
        b = pickle.load(handle)
    return b
def pickleFile(fname,obj):
    with open(fname, 'wb') as handle:
        pickle.dump(obj, handle)

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

order = int(0)
CFDSolver = DGSolver(airfoil, order=order, alpha=alpha,mach=1.5)
CFDSolver.solve()
cl = CFDSolver.postprocess()

dFdX_adjoint = CFDSolver.solveAdjoint()
print(dFdX_adjoint)
CFDSolver.writeSolution('airfoil_'+str(order)+'_high_mach')
solution = {
    'U' : CFDSolver.U,
    'psi' : CFDSolver.psi,
}
pickleFile('airfoil_'+str(order)+'_high_mach.pkl',solution)

# CFDSolver2 = DGSolver(airfoil, order=order,alpha=alpha+dalpha)
# CFDSolver2.solve()
# cl2 = CFDSolver2.postprocess()
# dFdX_FD = (cl2-cl)/dalpha
# print(dFdX_FD)
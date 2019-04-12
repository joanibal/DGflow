""" main script for proj4 """
import numpy as np
import matplotlib.pyplot as plt
from meshes import Mesh

from cfdsolvers import DGSolver
import pickle
# task 1


airfoil = Mesh('meshes/naca0012.gri')
# airfoil.plot()
# plt.show()

CFDSolver = DGSolver(airfoil, order=0,alpha=2)
CFDSolver.solve()
CFDSolver.writeSolution('airfoil')
CFDSolver.postprocess()
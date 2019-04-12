""" main script for proj2 """
import numpy as np
import matplotlib.pyplot as plt
from meshes import Mesh

from cfdsolvers import DGSolver
import pickle
# task 1


airfoil = Mesh('meshes/naca0012.gri')
# airfoil.plot()
# plt.show()

CFDSolver = DGSolver(airfoil, order=0)
CFDSolver.solve()
CFDSolver.writeSolution('airfoil')

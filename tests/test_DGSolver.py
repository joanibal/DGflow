from __future__ import print_function, division, absolute_import

import unittest
import numpy as np

from mesh import Mesh

from dg_solver import fluxes
from cfdsolvers import DGSolver


class TestBasic(unittest.TestCase):

    def setUp(self):
        test_mesh = Mesh('meshes/test0.gri')
        self.CFDSolver = DGSolver(test_mesh)

    def test_freestream(self):
        self.CFDSolver.solve(maxIter=1, freeStreamTest=True)
        self.assertAlmostEqual(self.CFDSolver.Rmax[0], 0.0)

    def test_freestream_preservation(self):
        self.CFDSolver.solve(maxIter=1000, freeStreamTest=True)
        self.assertAlmostEqual(self.CFDSolver.Rmax[0], 0.0)


class TestFluxes(unittest.TestCase):
    # for now these just test if it runs without error, not if the output is right

    def setUp(self):
        self.uL = np.array([2.0, 2.2, 0.1, 9.9])
        self.uR = np.array([0.9, 0.1, 2.2, 11.0])
        self.n = np.array([1, 1])*np.sqrt(2)/2
        # n = np.array([1, 0])
        # Ul = np.array([0.5,  0.1, 10.55])

    def test_consistency(self):
        U = np.vstack((self.uL, self.uL)).T
        F, _ = fluxes.roeflux(U, self.n)
        F_analytical = fluxes.analyticFlux(self.uL, self.n)

        print('test_consistency', F, F_analytical)
        np.testing.assert_array_almost_equal(F, F_analytical, decimal=6)

    def test_direction(self):
        U = np.vstack((self.uL, self.uR)).T
        F, _ = fluxes.roeflux(U, self.n)
        U = np.vstack((self.uR, self.uL)).T
        F_flipped, _ = fluxes.roeflux(U, -self.n)

        print('test_direction', F, F_flipped)
        np.testing.assert_array_almost_equal(F, F_flipped, decimal=6)

    def test_supersonic(self):
        uL = np.array([0.3, 2.25, 0.1, 10.55])
        gam = 1.4

        p = (gam - 1)*(uL[3] - 0.5*np.linalg.norm(uL[1:3])**2/uL[1])

        c = np.sqrt(gam*p/uL[1])
        M = np.linalg.norm(uL[1:3])/uL[0]/c
        print('M', M)

        U = np.vstack((self.uR, uL)).T
        F, _ = fluxes.roeflux(U, self.n)

        F_analytical = fluxes.analyticFlux(self.uL, self.n)

        print('test_supersonic', F, F_analytical)

        np.testing.assert_array_almost_equal(F, F_analytical, decimal=6)


class TestGeoModification(unittest.TestCase):

    def setUp(self):
        X = pyFoil._readCoordFile('rae2822.dat')
        self.foil = pyFoil.Airfoil(X)

    def test_reorder(self):
        Xorig = self.foil.spline.X
        newfoil = pyFoil.Airfoil(Xorig[::-1, :])
        Xreordered = newfoil.spline.X
        np.testing.assert_array_almost_equal(Xorig, Xreordered, decimal=6)

    def test_round_TE(self):
        self.foil.roundTE(k=4)
        refTE = np.array([0.990393, 0.0013401])
        newTE = self.foil.TE
        np.testing.assert_array_almost_equal(refTE, newTE, decimal=6)

    def test_blunt_TE(self):
        self.foil.makeBluntTE()
        refTE = np.array([0.97065494, 0.00352594])
        newTE = self.foil.TE
        np.testing.assert_array_almost_equal(refTE, newTE, decimal=8)


if __name__ == '__main__':
    unittest.main()

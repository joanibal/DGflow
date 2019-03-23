from __future__ import print_function, division, absolute_import

import unittest
import numpy as np

from meshes import Mesh

from dg_solver import fluxes
from cfdsolvers import DGSolver


# class TestBasic(unittest.TestCase):

#     def setUp(self):
#         test_mesh = Mesh('meshes/test.gri')
#         self.CFDSolver = DGSolver(test_mesh)

#     def test_freestream(self):
#         self.CFDSolver.solve(maxIter=1, freeStreamTest=True)
#         self.assertAlmostEqual(self.CFDSolver.Rmax[0], 0.0)

#     def test_freestream_preservation(self):
#         self.CFDSolver.solve(maxIter=1000, freeStreamTest=True)
#         self.assertAlmostEqual(self.CFDSolver.Rmax[0], 0.0)


class TestFluxes(unittest.TestCase):
    # for now these just test if it runs without error, not if the output is right

    def setUp(self):
        self.uL = np.array([2.0, 2.2, 0.1, 9.9])
        self.uR = np.array([0.9, 0.1, 2.2, 11.0])
        self.n = np.array([1, 1])*np.sqrt(2)/2

    def test_consistency(self):
        U = np.vstack((self.uL, self.uL)).T
        F, _ = fluxes.roeflux(U, self.n)
        F_analytic = fluxes.analyticflux(self.uL)
        F_analytic = F_analytic[0]*self.n[0] + F_analytic[1]*self.n[1]
        # print('test_consistency', F, F_analytic)
        np.testing.assert_array_almost_equal(F, F_analytic, decimal=6)

    def test_direction(self):
        U = np.vstack((self.uL, self.uR)).T
        F, _ = fluxes.roeflux(U, self.n)
        U = np.vstack((self.uR, self.uL)).T
        F_flipped, _ = fluxes.roeflux(U, -self.n)

        # print('test_direction', F, F_flipped)
        np.testing.assert_array_almost_equal(F, -1*F_flipped, decimal=6)

    def test_supersonic(self):
        uL = np.array([1.63, 2.53, 6.53, 25.88])

        # gam = 1.4
        # p = (gam - 1)*(uL[3] - 0.5*np.linalg.norm(uL[1:3])**2/uL[1])
        # c = np.sqrt(gam*p/uL[1])
        # M = np.linalg.norm(uL[1:3])/uL[0]/c
        # print('M', M)

        U = np.vstack((uL, self.uR)).T
        F, _ = fluxes.roeflux(U, self.n)

        F_analytic = fluxes.analyticflux(uL)
        F_analytic = F_analytic[0]*self.n[0] + F_analytic[1]*self.n[1]

        # print('test_supersonic', F, F_analytic)

        np.testing.assert_array_almost_equal(F, F_analytic, decimal=6)


class TestCurvedElement(unittest.TestCase):

    def setUp(self):
        self.linearMesh = Mesh('meshes/test0_1.gri')
        self.durvedMesh = Mesh('meshes/test0_2.gri')

    def test_linear_jacobian(self):
        invJ, detJ = self.linearMesh.getLinearJacobian()
        np.testing.assert_array_equal(invJ[0], np.eye(2))
        assert(detJ[0], 1)
    # def test_round_TE(self):
    #     self.foil.roundTE(k=4)
    #     refTE = np.array([0.990393, 0.0013401])
    #     newTE = self.foil.TE
    #     np.testing.assert_array_almost_equal(refTE, newTE, decimal=6)

    # def test_blunt_TE(self):
    #     self.foil.makeBluntTE()
    #     refTE = np.array([0.97065494, 0.00352594])
    #     newTE = self.foil.TE
    #     np.testing.assert_array_almost_equal(refTE, newTE, decimal=8)


if __name__ == '__main__':
    unittest.main()

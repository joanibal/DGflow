from __future__ import print_function, division, absolute_import

import unittest
import numpy as np

from meshes import Mesh

# from dg_solver import fluxes
from cfdsolvers import DGSolver
import dg_solver


def bumpShape(x):
    return 0.0625*np.exp(-25*x**2)


if True:
    bump = Mesh('meshes/bump0_kfid.gri', wallGeomFunc=bumpShape)


class TestResiduals(unittest.TestCase):

    def setUp(self):
        self.bump = Mesh('meshes/bump0_kfid.gri', wallGeomFunc=bumpShape)

    def test_freestream_preservation(self):
        for order in range(3):
            CFDSolver = DGSolver(self.bump, order=order)
            CFDSolver.testFreestream()
            self.assertAlmostEqual(CFDSolver.Rmax[0], 0.0)


class TestFluxes(unittest.TestCase):
    # for now these just test if it runs without error, not if the output is right

    def setUp(self):
        self.uL = np.array([2.0, 2.2, 0.1, 9.9])
        self.uR = np.array([0.9, 0.1, 2.2, 11.0])
        self.n = np.array([1, 1])*np.sqrt(2)/2

        # set external conditions for boundary fluxes
        dg_solver.constants.temptot_inf = 2.0
        dg_solver.constants.ptot_inf = 1.255
        dg_solver.constants.alpha = 0.1

        dg_solver.constants.p_inf = 0.66

# Roe flux

    def test_consistency(self):
        # U = np.vstack((self.uL, self.uL)).T
        F, _ = dg_solver.fluxes.roeflux(self.uL, self.uL, self.n)
        F_analytic = dg_solver.fluxes.analyticflux(self.uL).T
        F_analytic = F_analytic[0]*self.n[0] + F_analytic[1]*self.n[1]
        # print('test_consistency', F, F_analytic)
        np.testing.assert_array_almost_equal(F, F_analytic, decimal=6)

    def test_direction(self):
        F, _ = dg_solver.fluxes.roeflux(self.uL, self.uR, self.n)
        F_flipped, _ = dg_solver.fluxes.roeflux(self.uR, self.uL, -self.n)

        # print('test_direction', F, F_flipped)
        np.testing.assert_array_almost_equal(F, -1*F_flipped, decimal=6)

    def test_supersonic(self):
        uL = np.array([1.63, 2.53, 6.53, 25.88])

        # uncomment to know what mach number the state is at
        # gam = 1.4
        # p = (gam - 1)*(uL[3] - 0.5*np.linalg.norm(uL[1:3])**2/uL[1])
        # c = np.sqrt(gam*p/uL[1])
        # M = np.linalg.norm(uL[1:3])/uL[0]/c
        # print('M', M)

        F, _ = dg_solver.fluxes.roeflux(uL, self.uR, self.n)

        F_analytic = dg_solver.fluxes.analyticflux(uL).T
        F_analytic = F_analytic[0]*self.n[0] + F_analytic[1]*self.n[1]

        # print('test_supersonic', F, F_analytic)

        np.testing.assert_array_almost_equal(F, F_analytic, decimal=6)

# boundary fluxes
    # ***Caution*** the boundary condition fluxes are tested against previous values

    def test_wall(self):

        flux, s = dg_solver.fluxes.wallflux(self.uL, self.n)

        # old values
        flux_correct = np.array([0., 2.64422581, 2.64422581, 0.])
        s_correct = 2.37282019441

        np.testing.assert_array_almost_equal(flux, flux_correct, decimal=8)
        np.testing.assert_array_almost_equal(s, s_correct, decimal=8)

    def test_inflow(self):

        flux, s = dg_solver.fluxes.inflowflux(self.uL, self.n)

        # old values
        flux_correct = np.array([0.17487024, 0.92177938, 0.86403723, 1.22409169])
        s_correct = 1.950734432

        np.testing.assert_array_almost_equal(flux, flux_correct, decimal=8)
        np.testing.assert_array_almost_equal(s, s_correct, decimal=8)

    def test_outflow(self):

        flux, s = dg_solver.fluxes.outflowflux(self.uL, self.n)

        # old values
        flux_correct = np.array([1.50232517,  3.8692354,  2.29179397, 10.64555511])
        s_correct = 3.69070036

        np.testing.assert_array_almost_equal(flux, flux_correct, decimal=8)
        np.testing.assert_array_almost_equal(s, s_correct, decimal=8)
#


# class TestCurvedElement(unittest.TestCase):

#     def setUp(self):

#         def flatWall(x):
#             return 0

#         def curvWall(x):
#             return -x*(x-1)*0.2

#         # just a mesh of the reference element
#         self.refElem = Mesh('meshes/refElem.gri', wallGeomFunc=flatWall)

#         # both element types
#         self.curvMesh = Mesh('meshes/twoElem.gri', wallGeomFunc=curvWall)
#         # _, self.curvBasis = quadrules.getTriLagrangeBasis2D(q)

#     def test_linear_jacobian(self):
#         invJ, detJ = self.refElem.getLinearJacobian()
#         np.testing.assert_array_equal(invJ[0], np.eye(2))
#         assert detJ[0], 1

    # def test_curved_jacobian(self):
    #     # get the quadrature points for the test

    #     # order of the geom
    #     q = 3
    #     quadPts2D, quadWts2D = quadrules.getQuadPtsTri(q+1)

    #     self.curvMesh.getCurvedJacobian(0, quadPts2D, curvBasis)


if __name__ == '__main__':
    unittest.main()

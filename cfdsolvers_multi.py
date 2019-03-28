from __future__ import print_function
# from __future__ import division
# from readgri import readgri

import numpy as np
import matplotlib.pyplot as plt
from meshes import Mesh
# import flux_lib
import time
# from solver import fvsolver as solver
# from solver import mesh as solver_mesh
# from solver import constants as solver_constants
# from solver import postprocess as solver_post

from dg_solver import fluxes

import quadrules
# import solver

# TODO
"""

    # convert python to fortran to python for speed
    # do post processing for outputs

    #

"""


class DGSolver(object):
    def __init__(self, mesh, order=1):
        """
            class to solver for the flow in a given mesh
        """

        # self.mesh = mesh
        # self.U = np.zeros((mesh.nElem, 4))
        self.order = order
        self.mesh = mesh
        self.mesh.curvOrder = order+1
        self.nStates = 4

        # set BC data
        self.gamma = 1.4
        self.machInf = 0.5
        self.tempTotInf = 1 + (self.gamma - 1)/2 * self.machInf**2
        self.R = 1.0
        self.rhoInf = 1.
        self.P_inf = 1.
        self.tempInf = 1.

        self.Tt = (1 + 0.5*(self.gamma - 1)*self.machInf**2)*self.tempInf
        self.Ptot_inf = self.Tt**(self.gamma/(self.gamma - 1))

        # set mesh data

        # this method assumes that there is only one of each bc type

        self.curvWall = -1
        self.wall = -1
        self.inlet = -1
        self.outlet = -1
        for idx, bcname in enumerate(mesh.BCs.keys()):
            if 'curv' in bcname and 'wall' in bcname:
                self.curvWall = idx
            elif 'wall' in bcname:
                self.wall = idx
            elif 'inlet' in bcname:
                self.inlet = idx
            elif 'outlet' in bcname:
                self.outlet = idx
        # ====================================================
        #  Basis
        # ====================================================
        # get the basis function that are used to represent the solution
        self.nSolBasis, self.solBasis = quadrules.getTriLagrangeBasis2D(order)

        # number of degrees of freedom
        self.dof = self.mesh.nElem*self.nSolBasis

        # basis fucntion to represent the curved geometric elements
        _, self.curvBasis = quadrules.getTriLagrangeBasis2D(self.mesh.curvOrder)

        self.initFreestream()

        # ====================================================
        #  Qaud rules for the mesh (get seperate rules for the linear and curved elements)
        # ====================================================

        self.mesh.getHighOrderNodes()

        # for maximum python efficency these two nearly identical blocks should be a loop
        # but it works better with fortran this way

        # --------------------- Linear Element Quadrature ------------------------

        # these are 1d points form 0 to 1 (used for edge integration)
        self.linQuadPts1D, self.linQuadWts1D = quadrules.getQuadPts1D(order+1, 0, 1)
        self.nLinQuadPts1D = len(self.linQuadWts1D)

        # these are 2d points of a tri element in reference space
        self.linQuadPts2D, self.linQuadWts2D = quadrules.getQuadPtsTri(order+2)
        self.nLinQuadPts2D = len(self.linQuadWts2D)

        self.linPhi, self.linGPhi = self.getPhiMat(self.linQuadPts2D)
        self.lLinEdgePhi, self.rLinEdgePhi = self.getEdgePhi(self.linQuadPts1D)
        # get jacobian
        # ** note this is done for **all** elements because the edge residuals for the stright edges of the
        # curved element can be calculated using a linear jacobian
        self.linInvJ, self.linDetJ = self.mesh.getLinearJacobian()

        # calculate inverse mass matrix for each element
        self.invM = np.zeros((self.mesh.nElem, self.nSolBasis, self.nSolBasis))
        for k in self.mesh.linElem:
            M = self.linDetJ[k] * np.matmul(np.matmul(self.linPhi.T,
                                                      np.diag(self.linQuadWts2D)), self.linPhi)
            self.invM[k] = np.linalg.inv(M)

        # --------------------- Curved Element Quadrature -------------------------

        # these are 1d points form 0 to 1 (used for edge integration)
        self.curvQuadPts1D, self.curvQuadWts1D = quadrules.getQuadPts1D(self.mesh.curvOrder+1, 0, 1)
        self.nCurvQuadPts1D = len(self.curvQuadWts1D)

        # these are 2d points of a tri element in reference space
        self.curvQuadPts2D, self.curvQuadWts2D = quadrules.getQuadPtsTri(self.mesh.curvOrder+2)
        self.nCurvQuadPts2D = len(self.curvQuadWts2D)

        self.curvPhi, self.curvGPhi = self.getPhiMat(self.curvQuadPts2D)
        self.lCurvEdgePhi, self.rCurvEdgePhi = self.getEdgePhi(self.curvQuadPts1D)

        # get jacobian
        self.mesh.nCurvElem = len(self.mesh.curvElem)
        self.curvJ = np.zeros((self.mesh.nCurvElem, self.nCurvQuadPts2D, 2, 2))
        self.curvInvJ = np.zeros((self.mesh.nCurvElem, self.nCurvQuadPts2D, 2, 2))
        self.curvDetJ = np.zeros((self.mesh.nCurvElem, self.nCurvQuadPts2D))

        for idx_elem, elem in enumerate(self.mesh.curvElem):
            # print elem
            self.curvJ[idx_elem], self.curvInvJ[idx_elem], self.curvDetJ[idx_elem] = self.mesh.getCurvedJacobian(
                idx_elem, self.curvQuadPts2D, self.curvBasis)
        # calculate inverse mass matrix for each element
        for idx_elem, k in enumerate(self.mesh.curvElem):
            M = np.matmul(np.matmul(self.curvPhi.T, np.diag(
                self.curvDetJ[idx_elem]*self.curvQuadWts2D)), self.curvPhi)
            self.invM[k] = np.linalg.inv(M)

        self.curvDetJEdge, self.curvNormal = self.mesh.getEdgeJacobain(
            self.curvQuadPts1D, self.curvBasis)

        # TODO
        # tranfer precomputed values to the fortran layer

    def getPhiMat(self, quadPts2D):
        nQuadPts2D = len(quadPts2D)

        Phi = np.zeros((nQuadPts2D, self.nSolBasis))
        dPhi_dXi = np.zeros((nQuadPts2D, self.nSolBasis, self.mesh.nDim))

        # import ipdb
        # ipdb.set_trace()
        for idx, pt in enumerate(quadPts2D):
            # the basis function value at each of the quad points
            Phi[idx], dPhi_dXi[idx] = self.solBasis(pt)

        return Phi, dPhi_dXi

    # precompute the values of the basis functions are each edge of the reference element
    def getEdgePhi(self, quadPts1D):
        # 3 because there are three faces of a triangle
        nQuadPts1D = len(quadPts1D)

        leftEdgePhi = np.zeros((3, nQuadPts1D, self.nSolBasis))
        rightEdgePhi = np.zeros((3, nQuadPts1D, self.nSolBasis))
        for edge in range(3):
            # map 2D fave coordinates to
            pts = np.zeros((nQuadPts1D, 2))
            if edge == 0:
                pts[:, 0] = 1 - quadPts1D
                pts[:, 1] = quadPts1D
            elif edge == 1:
                pts[:, 0] = 0
                pts[:, 1] = 1 - quadPts1D

            elif edge == 2:
                pts[:, 0] = quadPts1D
                pts[:, 1] = 0

            for idx, pt in enumerate(pts):
                # the basis function value at each of the quad points
                leftEdgePhi[edge, idx], _ = self.solBasis(pt)

            for idx, pt in enumerate(pts[::-1]):
                # the self.solBasis function value at each of the quad points
                rightEdgePhi[edge, idx], _ = self.solBasis(pt)

        return leftEdgePhi, rightEdgePhi

    def initFreestream(self):
        # calculate conversed qualities
        c = np.sqrt(self.gamma*self.R*self.tempInf)

        u = self.machInf*c
        Ub = np.array([self.rhoInf,
                       self.rhoInf*u,
                       0.0,
                       self.P_inf/(self.gamma-1) + 0.5*self.rhoInf*u**2])

        # Ub = np.array([1, 0.5, 0, 2.5])
        self.U = np.zeros((self.mesh.nElem, self.nSolBasis, self.nStates))

        self.U[:, :] = Ub
        self.Ub = Ub
        print(self.U)

    # def initFromSolve(self, coarseMesh):

    #     for idx, _ in enumerate(coarseMesh.nElem):
    #         # set the four element that were created from uniform refinement
    #         self.U[:, 4*(idx-1):4*(idx)] = coarseMesh.U[idx]

    def setResidual(self, U):

        # loop over elements and compute residual contribution from interrior
        self.R = np.zeros(U.shape)
        self.S = np.zeros(self.mesh.nElem)

        self.setInteralResiduals(U)
        # print('internal', self.R)
        # quit()
        self.setEdgeResiduals(U)
        # print('sum', self.R)

    def setInteralResiduals(self, U):

        for idx_elem, elem in enumerate(self.mesh.linElem):

            for idx_basis in range(self.nSolBasis):
                Rtot = 0
                Uq = np.matmul(self.linPhi, U[elem])
                for qPt in range(self.nLinQuadPts2D):
                    flux = fluxes.analyticflux(Uq[qPt])
                    Rtot += np.dot(np.matmul(self.linGPhi[qPt, idx_basis], self.linInvJ[idx_elem]), flux) *\
                        self.linQuadWts2D[qPt]*self.linDetJ[idx_elem]

                self.R[elem, idx_basis] -= Rtot

        for idx_elem, elem in enumerate(self.mesh.curvElem):
            for idx_basis in range(self.nSolBasis):
                Rtot = 0
                Uq = np.matmul(self.curvPhi, U[elem])
                for qPt in range(self.nCurvQuadPts2D):
                    flux = fluxes.analyticflux(Uq[qPt])
                    Rtot += np.dot(np.matmul(self.curvGPhi[qPt, idx_basis], self.curvInvJ[idx_elem, qPt]), flux) *\
                        self.curvQuadWts2D[qPt]*self.curvDetJ[idx_elem, qPt]

                self.R[elem, idx_basis] -= Rtot

        # return R

    def setEdgeResiduals(self, U):
        for idx_edge in range(self.mesh.nInEdge):
            # ! get the elements connected to the edge

            idx_elem_left = self.mesh.inEdge2Elem[idx_edge, 0]
            idx_edge_left = self.mesh.inEdge2Elem[idx_edge, 1]

            idx_elem_right = self.mesh.inEdge2Elem[idx_edge, 2]
            idx_edge_right = self.mesh.inEdge2Elem[idx_edge, 3]

            uL = np.matmul(self.lLinEdgePhi[idx_edge_left], U[idx_elem_left])
            uR = np.matmul(self.rLinEdgePhi[idx_edge_right], U[idx_elem_right])

            for idx_basis in range(self.nSolBasis):
                # import ipdb
                # ipdb.set_trace()
                Rtot_left = np.zeros(self.nStates)
                Rtot_right = np.zeros(self.nStates)
                Stot_left = 0
                Stot_right = 0
                for q in range(self.nLinQuadPts1D):
                    U_edge = np.vstack((uL[q], uR[q])).T

                    flux, s = fluxes.roeflux(U_edge, self.mesh.inNormal[idx_edge])

                    # flux * delta L * wq
                    tmp = flux*self.mesh.inLength[idx_edge] * self.linQuadWts1D[q]

                    # import ipdb
                    # ipdb.set_trace()
                    Rtot_left += self.lLinEdgePhi[idx_edge_left, q, idx_basis] * tmp
                    Rtot_right += self.rLinEdgePhi[idx_edge_right, q, idx_basis] * tmp

                    Stot_left += s*self.mesh.inLength[idx_edge] * self.linQuadWts1D[q]
                    Stot_right += s*self.mesh.inLength[idx_edge] * self.linQuadWts1D[q]

                self.R[idx_elem_left, idx_basis] += Rtot_left
                self.R[idx_elem_right, idx_basis] -= Rtot_right

                self.S[idx_elem_left] += Stot_left
                self.S[idx_elem_right] += Stot_right

        # print('in edge', self.R)
        # import ipdb
        # ipdb.set_trace()
        # print(self.mesh.bcEdge2Elem)
        for idx_edge in range(self.mesh.nBCEdge):
            idx_elem = self.mesh.bcEdge2Elem[idx_edge, 0]

            idx_edge_loc = self.mesh.bcEdge2Elem[idx_edge, 1]
            bc = self.mesh.bcEdge2Elem[idx_edge, 2]
            # print(idx_elem, idx_edge,  bc)

            if bc == self.curvWall:
                edgePhi = self.lCurvEdgePhi[idx_edge_loc]
                nQuadPts = self.nCurvQuadPts1D
                quadWts = self.curvQuadWts1D
                elem = np.where(self.mesh.curvElem == idx_elem)[0][0]
                # print(bc, idx_edge, idx_elem)
                # import ipdb; ipdb.set_trace()
            else:
                edgePhi = self.lLinEdgePhi[idx_edge_loc]
                nQuadPts = self.nLinQuadPts1D
                quadWts = self.linQuadWts1D

            uL = np.matmul(edgePhi, U[idx_elem])

            # it was written this way to avoid checking the bc type in the loop of basis and q

            if bc == self.wall or bc == self.curvWall:
                for idx_basis in range(self.nSolBasis):
                    Rtot_left = np.zeros(self.nStates)
                    Stot_left = 0

                    for q in range(nQuadPts):

                        if bc == self.curvWall:
                            flux, s = fluxes.wallflux(uL[q], self.curvNormal[elem, idx_edge_loc, q])
                            tmp = self.curvDetJEdge[elem, idx_edge_loc, q]*quadWts[q]
                        else:
                            flux, s = fluxes.wallflux(uL[q], self.mesh.bcNormal[idx_edge])
                            tmp = self.mesh.bcLength[idx_edge]*quadWts[q]

                        Rtot_left += edgePhi[q, idx_basis] * flux * tmp
                        Stot_left += s*tmp

                    self.R[idx_elem, idx_basis] += Rtot_left
                    self.S[idx_elem] += Stot_left

                    # print('wall', bc,  Rtot_left)
                    # if bc == self.curvWall:
                    #     print('wall', elem,  bc == self.curvWall, Rtot_left)

            elif bc == self.inlet:
                for idx_basis in range(self.nSolBasis):
                    Rtot_left = np.zeros(self.nStates)
                    Stot_left = 0

                    for q in range(nQuadPts):

                        flux, s = fluxes.inflowflux(
                            self.tempTotInf, self.Ptot_inf, 0.0, uL[q], self.mesh.bcNormal[idx_edge])

                        # flux * delta L * wq
                        tmp = flux*self.mesh.bcLength[idx_edge] * self.linQuadWts1D[q]

                        Rtot_left += edgePhi[q, idx_basis] * tmp
                        Stot_left += s*self.mesh.bcLength[idx_edge]*self.linQuadWts1D[q]

                    self.R[idx_elem, idx_basis] += Rtot_left
                    self.S[idx_elem] += Stot_left
                    # print('inlet', Rtot_left)

            elif bc == self.outlet:
                for idx_basis in range(self.nSolBasis):
                    Rtot_left = np.zeros(self.nStates)
                    Stot_left = 0
                    for q in range(nQuadPts):

                        flux, s = fluxes.outflowflux(
                            self.P_inf, uL[q], self.mesh.bcNormal[idx_edge])

                        # flux * delta L * wq
                        tmp = flux*self.mesh.bcLength[idx_edge] * self.linQuadWts1D[q]

                        Rtot_left += edgePhi[q, idx_basis] * tmp
                        Stot_left += s*self.mesh.bcLength[idx_edge]*self.linQuadWts1D[q]

                    self.R[idx_elem, idx_basis] += Rtot_left
                    self.S[idx_elem] += Stot_left
                    # print('outlet', Rtot_left)

            else:
                print(bc)
                raise NotImplementedError

    def TVDRK2(self, cfl):
        # self.U.reshape()

        # dt = 0.004
        U_FE = np.zeros(self.U.shape)

        # print self.U

        self.setResidual(self.U)
        # print(self.R)
        # quit()
        dt = 2*self.mesh.area*cfl/self.S
        # print dt
        for idx_elem in range(self.mesh.nElem):
            # import ipdb
            # ipdb.set_trace()
            U_FE[idx_elem] = self.U[idx_elem] - dt[idx_elem] * \
                np.matmul(self.invM[idx_elem], self.R[idx_elem])

        self.setResidual(U_FE)
        for idx_elem in range(self.mesh.nElem):
            self.U[idx_elem] = 0.5*(self.U[idx_elem] + U_FE[idx_elem] -
                                    dt[idx_elem] * np.matmul(self.invM[idx_elem], self.R[idx_elem]))

        # print self.U

    def FE(self, cfl):
        # self.U.reshape()

        # dt = 0.004
        # U_FE = np.zeros(self.U.shape)

        # print self.U

        self.setResidual(self.U)
        # print(self.R)
        # quit()

        dt = 2*self.mesh.area*cfl/self.S
        # print dt
        for idx_elem in range(self.mesh.nElem):
            # import ipdb
            # ipdb.set_trace()
            self.U[idx_elem] = self.U[idx_elem] - dt[idx_elem] * \
                np.matmul(self.invM[idx_elem], self.R[idx_elem])

        # print

    def TVDRK3(self, cfl):
        # self.U.reshape()

        # dt = 0.004
        U_1 = np.zeros(self.U.shape)
        U_2 = np.zeros(self.U.shape)

        # print self.U

        self.setResidual(self.U)
        dt = 2*self.mesh.area*cfl/self.S
        # print(self.R)
        # quit()

        # print dt

        for idx_elem in range(self.mesh.nElem):
            U_1[idx_elem] = self.U[idx_elem] - dt[idx_elem] * \
                np.matmul(self.invM[idx_elem], self.R[idx_elem])

        self.setResidual(U_1)
        for idx_elem in range(self.mesh.nElem):
            U_2[idx_elem] = 3.0/4*self.U[idx_elem] + 1.0/4*U_1[idx_elem] - 1.0/4 * dt[idx_elem] * \
                np.matmul(self.invM[idx_elem], self.R[idx_elem])

        self.setResidual(U_2)
        for idx_elem in range(self.mesh.nElem):
            self.U[idx_elem] = 1.0/3*self.U[idx_elem] + 2.0/3*U_2[idx_elem] - 2.0/3 * dt[idx_elem] * \
                np.matmul(self.invM[idx_elem], self.R[idx_elem])

        # print self.U

    def JRK(self, cfl, nStages=3):
        U_stage = np.zeros(self.U.shape)
        # U_2 = np.zeros(self.U.shape)

        # print self.U

        self.setResidual(self.U)
        dt = 2*self.mesh.area*cfl/self.S
        # print(self.R)
        # quit()

        # print dt
        for ii in range(nStages, 1, -1):
            # U_stage = self.U
            for idx_elem in range(self.mesh.nElem):
                U_stage[idx_elem] = self.U[idx_elem] - dt[idx_elem]/ii * \
                    np.matmul(self.invM[idx_elem], self.R[idx_elem])
            self.setResidual(U_stage)

        for idx_elem in range(self.mesh.nElem):
            self.U[idx_elem] = self.U[idx_elem] - dt[idx_elem] * \
                np.matmul(self.invM[idx_elem], self.R[idx_elem])

        # print self.U

    def solve(self, maxIter=10000, tol=1e-7, cfl=0.95, freeStreamTest=False):

        # solver.cfl = cfl
        # solver_constants.mode = freeStreamTest

        # t = time.time()

        # if self.order == 1:
        #     self.Rmax, self.R2 = solver.solve_1storder(maxIter, tol)
        # elif self.order == 2:
        #     solver_constants.recon_bc_flux = self.recon_bc_flux
        # solver.u = np.load('bump0_kfid_2_sol.npz')['q'].T

        #     self.Rmax, self.R2 = solver.solve_2ndorder(maxIter, tol)
        # else:
        #     Exception('order not found')
        for iter in range(maxIter):
            self.TVDRK3(cfl)
            # self.FE(cfl)
            # self.TVDRK2(cfl)
            # self.JRK3(cfl)
            print(iter, np.max(np.max(np.max(np.abs(self.R)))))
        # self.wallTime = time.time() - t

        # self.U = solver.u.T
        # print(self.Rmax[:self.nIter])
        # import ipdb
        # ipdb.set_trace()
        # self.Rmax = self.Rmax[self.Rmax > 0]
        # self.nIter = len(self.Rmax)
        # print('wall time', self.wallTime, 'iters', self.nIter, 'Rmax', self.Rmax[-1])

    def postprocess(self):
        """
            get pressure and M for each cell
            """
        # self.Es = solver_post.getfeildvaribles()

        # get the edges of the bump
        idxs_bumpEdge = np.where(self.mesh.bcEdge2Elem[:, 2] == 1)[0]

        # initalize the force vector
        self.F = np.zeros(2)

        self.cp_wall = np.zeros(idxs_bumpEdge.shape[0])
        self.x_wall = np.zeros(idxs_bumpEdge.shape[0])

        # for each cell along the wall
        for idx, idx_edge in enumerate(idxs_bumpEdge):
            length = self.mesh.bcLength[idx_edge]
            normal = self.mesh.bcNormal[idx_edge]

            nodes = self.mesh.elem2Node[self.mesh.bcEdge2Elem[idx_edge, 0]]
            nodes = np.delete(nodes, [self.mesh.bcEdge2Elem[idx_edge, 1]])

            self.x_wall[idx] = np.mean(self.mesh.node2Pos[nodes], axis=0)[0]

            self.cp_wall[idx] = (self.p_wall[idx] - self.P_inf)
            self.F += -1*(self.P_inf - self.p_wall[idx])*normal*length

        # get state at each edge quad pts on the wall edge of the cell

        # map the quad pts in reffernce space to global (for plotting x vs Cp)

        # for second order use du_dx to find the value at the edge
        if self.order == 1 or not self.recon_p:
            # if True:
            bump_elem = self.mesh.bcEdge2Elem[idxs_bumpEdge, 0]
            # bump_local_edges = self.mesh.bcEdge2Elem[idxs_bumpEdge, 1]
            u_wall = self.U[bump_elem, 1:3]
            r_wall = self.U[bump_elem, 0]
            # import ipdb
            # ipdb.set_trace()

            unL = np.sum(u_wall*self.mesh.bcNormal[idxs_bumpEdge], axis=1)/r_wall
            qL = np.linalg.norm(u_wall, axis=1)/r_wall
            utL = np.sqrt(qL**2 - unL**2)

            self.p_wall = (self.gamma-1)*(self.U[bump_elem, 3] - 0.5*r_wall*utL**2)

            # self.p_wall = solver_post.p[self.mesh.bcEdge2Elem[idxs_bumpEdge, 0]]
        else:
            self.p_wall = solver_post.getwallpressure(idxs_bumpEdge+1)
        for idx, idx_edge in enumerate(idxs_bumpEdge):
            length = self.mesh.bcLength[idx_edge]
            normal = self.mesh.bcNormal[idx_edge]

            nodes = self.mesh.elem2Node[self.mesh.bcEdge2Elem[idx_edge, 0]]
            nodes = np.delete(nodes, [self.mesh.bcEdge2Elem[idx_edge, 1]])

            self.x_wall[idx] = np.mean(self.mesh.node2Pos[nodes], axis=0)[0]

            self.cp_wall[idx] = (self.p_wall[idx] - self.P_inf)
            self.F += -1*(self.P_inf - self.p_wall[idx])*normal*length

        h = 0.0625
        self.cp_wall /= (self.gamma/2*self.P_inf*self.machInf**2)
        self.cd, self.cl = self.F/(self.gamma/2*self.P_inf*self.machInf**2*h)

        print('cd', self.cd, 'cl', self.cl, 'Es', self.Es)
        idxs_sorted = np.argsort(self.x_wall)
        self.x_wall = self.x_wall[idxs_sorted]
        self.p_wall = self.p_wall[idxs_sorted]
        self.cp_wall = self.cp_wall[idxs_sorted]

    def plotResiduals(self):
        plt.semilogy(range(1, self.nIter+1), self.Rmax, label='order: ' + str(self.order))
        plt.xlabel('Iteration', fontsize=16)
        plt.ylabel(r'$|R_{\infty }|$',  fontsize=16)
        plt.title('Residual History',  fontsize=16)

    def plotCP(self):
        plt.plot(self.x_wall, self.cp_wall, '-',  label='order: ' + str(self.order))
        plt.xlabel(r'$X$', fontsize=16)
        plt.ylabel(r'$C_p$', fontsize=16)
        plt.title(r'$C_p$ along the bump')
        plt.gca().invert_yaxis()

    def writeSolution(self, fileName):
        # def writeLine(cmd)

        # so all the array values are printed
        np.set_printoptions(threshold=np.inf)

        # for each set of elements of the same geoemtric order
        with open(fileName + '.dat', 'w') as fid:

            fid.write('TITLE = "bump"\n')
            # fid.write('Variables="X", "Y", "U", "V", "rho", "P", "M", "rhoRes"\n')
            fid.write('Variables="X", "Y", "U", "V", "rho"\n')

            elemOrder = {
                1: self.mesh.linElem,
                self.mesh.curvOrder: self.mesh.curvElem
            }

            for q in elemOrder.keys():

                # get the basis functions for the mapping
                # this is a little bit of extra work for the linear elements, but by not utilizing the constant jacobian
                #  the next loop can be written without conditional statments (if q==1) which is nice and clean

                #
                nBasis, basis = quadrules.getTriLagrangeBasis2D(q)

                # these are the points(in reference space) where the function will be evaluated
                interpOrder = min([self.order+3, q+4])
                Xi = quadrules.getTriLagrangePts2D(interpOrder)

                # create element connectivity
                N = len(Xi)

                nodes = np.arange(1, N+1)

                rows = []
                k = 0
                for ii in range(1, interpOrder+1)[::-1]:
                    rows.append(nodes[k:k+ii+1])
                    k += ii+1

                # conn = np.zeros(interpOrder**2, 3)
                conn = []
                # idx_conn = 0
                for idx, row in enumerate(rows):
                    # do add the top and bottom elements to the connectivity matrix
                    k = len(row)

                    if idx > 0:
                        # add the elements below this row of nodes
                        for i in range(k-1):
                            conn.append([row[i], row[i+1], row[i]-k])

                    # add the elements above the row
                    for i in range(k-1):
                        conn.append([row[i], row[i+1], row[i+1]+k-1])
                conn = np.array(conn)
                # import ipdb
                # ipdb.set_trace()

                geomPhi = np.zeros((len(Xi), nBasis))
                solPhi = np.zeros((len(Xi), self.nSolBasis))

                for idx, pt in enumerate(Xi):
                    # the basis function value at each of the quad points
                    geomPhi[idx], _ = basis(pt)
                    solPhi[idx], _ = self.solBasis(pt)

                # loop over elements now
                for idx_elem, elem in enumerate(elemOrder[q]):
                    # nodesPos = self.mesh.node2Pos[self.mesh.elem2Node[elem]]
                    if q == 1:
                        nodesPos = self.mesh.node2Pos[self.mesh.elem2Node[elem]]
                    else:
                        nodesPos = self.mesh.curvNodes[idx_elem]
                    # import ipdb
                    # ipdb.set_trace()
                    Upts = np.matmul(solPhi, self.U[elem])
                    Xpts = np.matmul(geomPhi, nodesPos)

                    # writeZone()
                    # print(elem)
                    fid.write('ZONE\n')
                    fid.write('T="elem '+str(elem)+'"\n')
                    fid.write('DataPacking=Block\n')
                    fid.write('ZoneType=FETRIANGLE\n')
                    fid.write('N=' + str(N) +
                              ' E=' + str(interpOrder**2)+'\n')
                    # fid.write('VarLocation=([3-9]=CellCentered)\n')

                    fid.write('#XData (Grid Variables must be nodal)\n')
                    fid.write(str(Xpts[:, 0])[1:-1]+'\n')
                    fid.write('#YData (Grid Variables must be nodal)\n')
                    fid.write(str(Xpts[:, 1])[1:-1]+'\n')

                    fid.write('#U Data \n')
                    fid.write(str(Upts[:, 1]/Upts[:, 0])[1:-1]+'\n')

                    fid.write('#V Data \n')
                    fid.write(str(Upts[:, 2]/Upts[:, 0])[1:-1]+'\n')

                    fid.write('#rho Data \n')
                    fid.write(str(Upts[:, 0])[1:-1]+'\n')

                    # fid.write('#P Data \n')
                    # fid.write(str(solver_post.p)[1:-1]+'\n')
                    # fid.write('#M Data \n')
                    # fid.write(str(solver_post.mach)[1:-1]+'\n')
                    # fid.write('#rhoRes Data \n')
                    # fid.write(str(solver.res[0, :])[1:-1]+'\n')

                    # fid.write('#number Data \n')
                    # fid.write(str(np.arange(self.mesh.nElem)+1)[1:-1]+'\n')

                    fid.write('#Connectivity List\n')
                    for idx in range(len(conn)):
                        fid.write(str(conn[idx])[1:-1]+'\n')

        # set np back to normal
        np.set_printoptions(threshold=1000)

    def writeSolution_lite(self, fileName):
        # def writeLine(cmd)

        # so all the array values are printed
        np.set_printoptions(threshold=np.inf)

        with open(fileName + '.dat', 'w') as fid:
            fid.write('TITLE = "bump"\n')
            fid.write(
                'Variables="X", "Y", "U", "V", "rho", "rhoRes", "num"\n')
            fid.write('ZONE\n')
            fid.write('T="Zone Title"\n')
            fid.write('DataPacking=Block\n')
            fid.write('ZoneType=FETRIANGLE\n')
            fid.write('N='+str(self.mesh.nNodes) +
                      ' E=' + str(self.mesh.nElem)+'\n')
            fid.write('VarLocation=([3-7]=CellCentered)\n')

            fid.write('#XData (Grid Variables must be nodal)\n')
            fid.write(str(self.mesh.node2Pos[:, 0])[1:-1]+'\n')
            fid.write('#YData (Grid Variables must be nodal)\n')
            fid.write(str(self.mesh.node2Pos[:, 1])[1:-1]+'\n')

            fid.write('#U Data \n')
            fid.write(str(self.U[:, 0, 1]/self.U[:, 0, 0])[1:-1]+'\n')

            fid.write('#V Data \n')
            fid.write(str(self.U[:, 0, 2]/self.U[:, 0, 0])[1:-1]+'\n')

            fid.write('#rho Data \n')
            fid.write(str(self.U[:, 0, 0])[1:-1]+'\n')

            fid.write('#rhoRes Data \n')
            fid.write(str(self.R[:, 0, 0])[1:-1]+'\n')

            fid.write('#number Data \n')
            fid.write(str(np.arange(self.mesh.nElem)+1)[1:-1]+'\n')

            fid.write('#Connectivity List\n')
            for jj in range(self.mesh.nElem):
                fid.write(str(self.mesh.elem2Node[jj]+1)[1:-1]+'\n')

        # set np back to normal
        np.set_printoptions(threshold=1000)


if __name__ == '__main__':
    # bump = Mesh('meshes/test0_21.gri')

    # def flatWall(x):
    #     return -x*(x-1)*0.2
    #     # return 0
    def bumpShape(x):
        return 0.0625*np.exp(-25*x**2)

    bump = Mesh('meshes/bump0_kfid.gri', wallGeomFunc=bumpShape)
    # test = Mesh('meshes/test0_2.gri', wallGeomFunc=flatWall)
    DGSolver = DGSolver(bump, order=2)
    # DGprint(FVSolver.getResidual())
    DGSolver.solve(maxIter=2000)
    # DGSolver.postprocess()
    DGSolver.writeSolution('test_multi_2')

    # DGFVSolver.
    # DGFVSolver.solve()
    # DGFVSolver.postprocess()
    # DGFVSolver.writeSolution('lv0')
    # DGFVSolver.plotCP()
    # DGFVSolver.plotResiduals()

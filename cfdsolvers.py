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
    # need for freestream test
    calculate interrior residual
    calc edge residual
        - interrroir residual
        - edge residual

    # convert python to fortran to python for speed

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
        temp_wall = []
        temp_inlet = []
        temp_outlet = []
        for idx, bcname in enumerate(mesh.BCs.keys()):
            if 'wall' in bcname:
                temp_wall.append(idx)
            if 'inlet' in bcname:
                temp_inlet.append(idx)
            if 'outlet' in bcname:
                temp_outlet.append(idx)

        self.wall = temp_wall
        self.inlet = temp_inlet
        self.outlet = temp_outlet

        # % quadrature points sufficient to integrate 2*p in 1D and 2D

        self.nQuadPts1D = order+1
        # these are 1d points form 0 to 1 (used for edge integration)
        self.quadPts1D, self.quadWts1D = quadrules.getQuadPts1D(order+1, 0, 1)

        # these are 2d points of a tri element in reference space
        self.quadPts2D, self.quadWts2D = quadrules.getQuadPtsTri(order+1)
        self.nQuadPts2D = len(self.quadWts2D)

        # nQuadPts2D, quadPts2D, quadWts2D = quadrules.quad2d(nQuadPts1D, 0, 1)

        self.nBasis, self.basis = quadrules.getTriLagrangeBasis2D(order)
        # number of degrees of freedom
        # Ndof = nelem*nn

        self.initFreestream()

        # % calculate inverse mass matrix for each element
        self.Phi = np.zeros((self.nQuadPts2D, self.nBasis))
        self.dPhi_dXi = np.zeros((self.nQuadPts2D, self.nBasis, self.mesh.nDim))

        # import ipdb
        # ipdb.set_trace()
        for idx, pt in enumerate(self.quadPts2D):
            # the basis function value at each of the quad points
            self.Phi[idx], self.dPhi_dXi[idx] = self.basis(pt)

        self.invM = np.zeros((self.mesh.nElem, self.nBasis, self.nBasis))

        self.invJ, self.detJ = self.mesh.getLinearJacobian()
        for k in range(self.mesh.nElem):
            M = self.detJ[k] * np.matmul(np.matmul(self.Phi.T, np.diag(self.quadWts2D)), self.Phi)
            self.invM[k] = np.linalg.inv(M)

        # precompute the values of the basis functions are each edge of the reference element

        # 3 because there are three faces of a triangel
        self.leftEdgePhi = np.zeros((3, self.nQuadPts1D, self.nBasis))
        self.rightEdgePhi = np.zeros((3, self.nQuadPts1D, self.nBasis))
        for edge in range(3):
            # map 2D fave coordinates to
            pts = np.zeros((self.nQuadPts1D, 2))
            if edge == 0:
                pts[:, 0] = 1 - self.quadPts1D
                pts[:, 1] = self.quadPts1D
            elif edge == 1:
                pts[:, 0] = 0
                pts[:, 1] = 1 - self.quadPts1D

            elif edge == 2:
                pts[:, 0] = self.quadPts1D
                pts[:, 1] = 0

            for idx, pt in enumerate(pts):
                # the basis function value at each of the quad points
                self.leftEdgePhi[edge, idx], _ = self.basis(pt)

            for idx, pt in enumerate(pts[::-1]):
                # the basis function value at each of the quad points
                self.rightEdgePhi[edge, idx], _ = self.basis(pt)

        # import ipdb
        # ipdb.set_trace()

    def initFreestream(self):
        # calculate conversed qualities
        c = np.sqrt(self.gamma*self.R*self.tempInf)

        u = self.machInf*c
        Ub = np.array([self.rhoInf,
                       self.rhoInf*u,
                       0.0,
                       self.P_inf/(self.gamma-1) + 0.5*self.rhoInf*u**2])

        # Ub = np.array([1, 0.5, 0, 2.5])
        self.U = np.zeros((self.mesh.nElem, self.nBasis, self.nStates))

        self.U[:, :] = Ub
        self.Ub = Ub
        print(self.U)

    # def initFromSolve(self, coarseMesh):

    #     for idx, _ in enumerate(coarseMesh.nElem):
    #         # set the four element that were created from uniform refinement
    #         self.U[:, 4*(idx-1):4*(idx)] = coarseMesh.U[idx]

    def getResidual(self):

        # loop over elements and compute residual contribution from interrior
        self.R = np.zeros(self.U.shape)
        self.S = np.zeros(self.mesh.nElem)

        self.getInteralResiduals()
        # print R
        self.getEdgeResiduals()

        # return R
        # calculate analytical flux

    def getInteralResiduals(self):

        for idx_elem in range(self.mesh.nElem):

            for idx_basis in range(self.nBasis):
                Rtot = 0
                Uq = np.matmul(self.Phi, self.U[idx_elem])
                for q in range(self.nQuadPts2D):
                    flux = fluxes.analyticflux(Uq[q])
                    Rtot += np.dot(np.matmul(self.dPhi_dXi[q, idx_basis], self.invJ[idx_elem]), flux) *\
                        self.quadWts2D[q]*self.detJ[idx_elem]

                self.R[idx_elem, idx_basis] -= Rtot
        # return R

    def getEdgeResiduals(self):
        for idx_edge in range(self.mesh.nInEdge):
            # ! get the elements connected to the edge

            idx_elem_left = self.mesh.inEdge2Elem[idx_edge, 0]
            idx_edge_left = self.mesh.inEdge2Elem[idx_edge, 1]
            idx_elem_right = self.mesh.inEdge2Elem[idx_edge, 2]
            idx_edge_right = self.mesh.inEdge2Elem[idx_edge, 3]

            uL = np.matmul(self.leftEdgePhi[idx_edge_left], self.U[idx_elem_left])
            uR = np.matmul(self.rightEdgePhi[idx_edge_right], self.U[idx_elem_right])

            for idx_basis in range(self.nBasis):
                import ipdb
                ipdb.set_trace()
                Rtot_left = np.zeros(self.nStates)
                Rtot_right = np.zeros(self.nStates)
                Stot_left = 0
                Stot_right = 0
                for q in range(self.nQuadPts1D):
                    U_edge = np.vstack((uL[q], uR[q])).T

                    flux, s = fluxes.roeflux(U_edge, self.mesh.inNormal[idx_edge])

                    # flux * delta L * wq
                    tmp = flux*self.mesh.inLength[idx_edge] * self.quadWts1D[q]

                    # import ipdb
                    # ipdb.set_trace()
                    Rtot_left += self.leftEdgePhi[idx_edge_left, q, idx_basis] * tmp
                    Rtot_right += self.rightEdgePhi[idx_edge_right, q, idx_basis] * tmp

                    Stot_left += s*self.mesh.inLength[idx_edge] * self.quadWts1D[q]
                    Stot_right += s*self.mesh.inLength[idx_edge] * self.quadWts1D[q]

                self.R[idx_elem_left, idx_basis] += Rtot_left
                self.R[idx_elem_right, idx_basis] -= Rtot_right

                self.S[idx_elem_left] += Stot_left
                self.S[idx_elem_right] += Stot_right

        for idx_edge in range(self.mesh.nBCEdge):
            idx_elem_left = self.mesh.bcEdge2Elem[idx_edge, 0]
            idx_edge_left = self.mesh.bcEdge2Elem[idx_edge, 1]
            bc = self.mesh.bcEdge2Elem[idx_edge, 2]

            uL = np.matmul(self.leftEdgePhi[idx_edge_left], self.U[idx_elem_left])

            if any(bc == self.wall):
                for idx_basis in range(self.nBasis):
                    Rtot_left = np.zeros(self.nStates)
                    Stot_left = 0

                    for q in range(self.nQuadPts1D):

                        flux, s = fluxes.wallflux(uL[q], self.mesh.bcNormal[idx_edge])

                        # flux * delta L * wq
                        tmp = flux*self.mesh.bcLength[idx_edge] * self.quadWts1D[q]

                        Rtot_left += self.leftEdgePhi[idx_edge_left, q, idx_basis] * tmp
                        Stot_left += s*self.mesh.bcLength[idx_edge]*self.quadWts1D[q]

                    self.R[idx_elem_left, idx_basis] += Rtot_left
                    self.S[idx_elem_left] += Stot_left

            elif any(bc == self.inlet):
                for idx_basis in range(self.nBasis):
                    Rtot_left = np.zeros(self.nStates)
                    Stot_left = 0

                    for q in range(self.nQuadPts1D):

                        flux, s = fluxes.inflowflux(
                            self.tempTotInf, self.Ptot_inf, 0.0, uL[q], self.mesh.bcNormal[idx_edge])

                        # flux * delta L * wq
                        tmp = flux*self.mesh.bcLength[idx_edge] * self.quadWts1D[q]

                        Rtot_left += self.leftEdgePhi[idx_edge_left, q, idx_basis] * tmp
                        Stot_left += s*self.mesh.bcLength[idx_edge]*self.quadWts1D[q]

                    self.R[idx_elem_left, idx_basis] += Rtot_left
                    self.S[idx_elem_left] += Stot_left

            elif any(bc == self.outlet):
                for idx_basis in range(self.nBasis):
                    Rtot_left = np.zeros(self.nStates)
                    Stot_left = 0
                    for q in range(self.nQuadPts1D):

                        flux, s = fluxes.outflowflux(
                            self.P_inf, uL[q], self.mesh.bcNormal[idx_edge])

                        # flux * delta L * wq
                        tmp = flux*self.mesh.bcLength[idx_edge] * self.quadWts1D[q]

                        Rtot_left += self.leftEdgePhi[idx_edge_left, q, idx_basis] * tmp
                        Stot_left += s*self.mesh.bcLength[idx_edge]*self.quadWts1D[q]

                    self.R[idx_elem_left, idx_basis] += Rtot_left
                    self.S[idx_elem_left] += Stot_left

            else:
                print(bc)
                raise NotImplementedError

            # for idx_basis in range(self.nBasis):
            #     Rtot_left = np.zeros(self.nStates)
            #     Stot_left = 0
            #     for q in range(self.nQuadPts1D):

            #         U_edge = np.vstack((uL[q], self.Ub)).T

            #         flux, s = fluxes.roeflux(U_edge, self.mesh.bcNormal[idx_edge])

            #         # flux * delta L * wq
            #         tmp = flux*self.mesh.bcLength[idx_edge] * self.quadWts1D[q]

            #         Rtot_left += self.leftEdgePhi[idx_edge_left, q, idx_basis] * tmp
            #         Stot_left += s*self.mesh.bcLength[idx_edge]*self.quadWts1D[q]

            #     self.R[idx_elem_left, idx_basis] += Rtot_left
            #     self.S[idx_elem_left] += Stot_left

        # return R

    def RK2(self, cfl):
        # self.U.reshape()

        # dt = 0.004
        U_FE = np.zeros(self.U.shape)

        # print self.U

        self.getResidual()
        dt = 2*self.mesh.area*cfl/self.S
        # print dt
        for idx_elem in range(self.mesh.nElem):
            # import ipdb
            # ipdb.set_trace()
            U_FE[idx_elem] = self.U[idx_elem] - dt[idx_elem] * \
                np.matmul(self.invM[idx_elem], self.R[idx_elem])

        self.getResidual()
        for idx_elem in range(self.mesh.nElem):
            self.U[idx_elem] = 0.5*(self.U[idx_elem] + U_FE[idx_elem] -
                                    dt[idx_elem] * np.matmul(self.invM[idx_elem], self.R[idx_elem]))

        # print self.U

    def solve(self, maxIter=2000, tol=1e-7, cfl=0.9, freeStreamTest=False):

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
            self.RK2(cfl)
            print iter, np.max(np.max(np.max(np.abs(self.R))))
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
        self.Es = solver_post.getfeildvaribles()

        # get the edges of the bump
        idxs_bumpEdge = np.where(self.mesh.bcEdge2Elem[:, 2] == 1)[0]

        # initalize the force vector
        self.F = np.zeros(2)

        self.cp_wall = np.zeros(idxs_bumpEdge.shape[0])
        self.x_wall = np.zeros(idxs_bumpEdge.shape[0])

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

        with open(fileName + '.dat', 'w') as fid:
            fid.write('TITLE = "bump"\n')
            fid.write(
                'Variables="X", "Y", "U", "V", "rho", "P", "M", "rhoRes", "num"\n')
            fid.write('ZONE\n')
            fid.write('T="Zone Title"\n')
            fid.write('DataPacking=Block\n')
            fid.write('ZoneType=FETRIANGLE\n')
            fid.write('N='+str(self.mesh.nNodes) +
                      ' E=' + str(self.mesh.nElem)+'\n')
            fid.write('VarLocation=([3-9]=CellCentered)\n')

            fid.write('#XData (Grid Variables must be nodal)\n')
            fid.write(str(self.mesh.node2Pos[:, 0])[1:-1]+'\n')
            fid.write('#YData (Grid Variables must be nodal)\n')
            fid.write(str(self.mesh.node2Pos[:, 1])[1:-1]+'\n')

            fid.write('#U Data \n')
            fid.write(str(self.U[:, 1]/self.U[:, 0])[1:-1]+'\n')

            fid.write('#V Data \n')
            fid.write(str(self.U[:, 2]/self.U[:, 0])[1:-1]+'\n')

            fid.write('#rho Data \n')
            fid.write(str(self.U[:, 0])[1:-1]+'\n')

            fid.write('#P Data \n')
            fid.write(str(solver_post.p)[1:-1]+'\n')
            fid.write('#M Data \n')
            fid.write(str(solver_post.mach)[1:-1]+'\n')
            fid.write('#rhoRes Data \n')
            fid.write(str(solver.res[0, :])[1:-1]+'\n')

            fid.write('#number Data \n')
            fid.write(str(np.arange(self.mesh.nElem)+1)[1:-1]+'\n')

            fid.write('#Connectivity List\n')
            for jj in range(self.mesh.nElem):
                fid.write(str(self.mesh.elem2Node[jj]+1)[1:-1]+'\n')

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
    bump = Mesh('meshes/bump0.gri')
    # bump = Mesh('meshes/test2.gri')

    DGSolver = DGSolver(bump, order=0)
    # DGprint(FVSolver.getResidual())
    DGSolver.solve()
    # DGSolver.writeSolution_lite('test2')
    # DGFVSolver.
    # DGFVSolver.solve()
    # DGFVSolver.postprocess()
    # DGFVSolver.writeSolution('lv0')
    # DGFVSolver.plotCP()
    # DGFVSolver.plotResiduals()

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
        self.geomOrder = order+1
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

        # ====================================================
        #  Qaud rules for the solution
        # ====================================================
        # % quadrature points sufficient to integrate 2*p in 1D and 2D

        # these are 1d points form 0 to 1 (used for edge integration)
        self.quadPts1D, self.quadWts1D = quadrules.getQuadPts1D(order+1, 0, 1)
        self.nQuadPts1D = len(self.quadWts1D)

        # these are 2d points of a tri element in reference space
        self.quadPts2D, self.quadWts2D = quadrules.getQuadPtsTri(order+1)
        self.nQuadPts2D = len(self.quadWts2D)

        # nQuadPts2D, quadPts2D, quadWts2D = quadrules.quad2d(nQuadPts1D, 0, 1)

        self.nBasis, self.basis = quadrules.getTriLagrangeBasis2D(order)
        # number of degrees of freedom
        # Ndof = nelem*nn

        self.Phi = np.zeros((self.nQuadPts2D, self.nBasis))
        self.dPhi_dXi = np.zeros((self.nQuadPts2D, self.nBasis, self.mesh.nDim))

        # import ipdb
        # ipdb.set_trace()
        for idx, pt in enumerate(self.quadPts2D):
            # the basis function value at each of the quad points
            self.Phi[idx], self.dPhi_dXi[idx] = self.basis(pt)

        self.initFreestream()

        # ====================================================
        #  Qaud rules for the mesh (assuming all elements are high order)
        # ====================================================

        # can't use dictionaries in fortran so unpack
        self.elemOrder = self.mesh.elemOrder.keys()
        self.elems = self.mesh.elemOrder. values()

        self.mesh.getHighOrderNodes()

       # these are 1d points form 0 to 1 (used for edge integration)
        self.quadPts1DElem, self.quadWts1DElem = quadrules.getQuadPts1D(self.geomOrder+1, 0, 1)
        self.nQuadPts1DElem = len(self.quadWts1DElem)

        # these are 2d points of a tri element in reference space
        self.quadPts2DElem, self.quadWts2DElem = quadrules.getQuadPtsTri(self.geomOrder+1)
        self.nQuadPts2DElem = len(self.quadWts2DElem)

        # nQuadPts2DElem, quadPts2DElem, quadWts2DElem = quadrules.quad2d(nQuadPts1DElem, 0, 1)

        self.nBasisElem, self.basisElem = quadrules.getTriLagrangeBasis2D(self.geomOrder)

        # % calculate inverse mass matrix for each element
        self.Phi_Elem = np.zeros((self.nQuadPts2DElem, self.nBasis))
        self.dPhi_dXi_Elem = np.zeros((self.nQuadPts2DElem, self.nBasis, self.mesh.nDim))

        # import ipdb
        # ipdb.set_trace()
        for idx, pt in enumerate(self.quadPts2DElem):
            # the basis function value at each of the quad points
            self.Phi_Elem[idx], self.dPhi_dXi_Elem[idx] = self.basis(pt)

        self.invM = np.zeros((self.mesh.nElem, self.nBasis, self.nBasis))

        # self.invJ, self.detJ = self.mesh.getLinearJacobian()

        # get jacobian
        self.J = np.zeros((self.mesh.nElem, self.nQuadPts2DElem, 2, 2))
        self.invJ = np.zeros((self.mesh.nElem, self.nQuadPts2DElem, 2, 2))
        self.detJ = np.zeros((self.mesh.nElem, self.nQuadPts2DElem))

        for elem in range(self.mesh.nElem):
            # print elem
            self.J[elem], self.invJ[elem], self.detJ[elem] = self.mesh.getCurvedJacobian(self.geomOrder,
                                                                                         elem, self.quadPts2DElem, self.basisElem)[:3]
        for k in range(self.mesh.nElem):
            M = np.matmul(np.matmul(self.Phi_Elem.T, np.diag(
                self.detJ[k]*self.quadWts2DElem)), self.Phi_Elem)
            self.invM[k] = np.linalg.inv(M)

        # precompute the values of the basis functions are each edge of the reference element

        # 3 because there are three faces of a triangle
        self.leftEdgePhi = np.zeros((3, self.nQuadPts1DElem, self.nBasis))
        self.rightEdgePhi = np.zeros((3, self.nQuadPts1DElem, self.nBasis))
        for edge in range(3):
            # map 2D fave coordinates to
            pts = np.zeros((self.nQuadPts1DElem, 2))
            if edge == 0:
                pts[:, 0] = 1 - self.quadPts1DElem
                pts[:, 1] = self.quadPts1DElem
            elif edge == 1:
                pts[:, 0] = 0
                pts[:, 1] = 1 - self.quadPts1DElem

            elif edge == 2:
                pts[:, 0] = self.quadPts1DElem
                pts[:, 1] = 0

            for idx, pt in enumerate(pts):
                # the basis function value at each of the quad points
                self.leftEdgePhi[edge, idx], _ = self.basis(pt)

            for idx, pt in enumerate(pts[::-1]):
                # the basis function value at each of the quad points
                self.rightEdgePhi[edge, idx], _ = self.basis(pt)

        # calculate Jedge
        # for each element, for each edge, for each quad pt
        self.detJEdge = np.zeros((self.mesh.nElem,  3, self.nQuadPts1DElem))
        self.normalEdge = np.zeros((self.mesh.nElem,  3, self.nQuadPts1DElem, 2))
        for elem in range(self.mesh.nElem):
            for edge in range(3):
                pts = np.zeros((self.nQuadPts1DElem, 2))
                if edge == 0:
                    pts[:, 0] = 1 - self.quadPts1DElem
                    pts[:, 1] = self.quadPts1DElem
                    dXi_dX = np.array([-1, 1])

                elif edge == 1:
                    pts[:, 0] = 0
                    pts[:, 1] = 1 - self.quadPts1DElem
                    dXi_dX = np.array([0, -1])

                elif edge == 2:
                    pts[:, 0] = self.quadPts1DElem
                    pts[:, 1] = 0
                    dXi_dX = np.array([1, 0])

                J = self.mesh.getCurvedJacobian(
                    self.geomOrder, elem, pts, self.basisElem)[0]

                for q in range(len(J)):
                    tang_vec = J[q][:, 0]*dXi_dX[0] + J[q][:, 1]*dXi_dX[1]
                    self.normalEdge[elem][edge][q] = np.array([tang_vec[1], -tang_vec[0]])

                    self.detJEdge[elem][edge][q] = np.linalg.norm(self.normalEdge[elem][edge][q])
                    self.normalEdge[elem][edge][q] /= self.detJEdge[elem][edge][q]

                print(edge, self.normalEdge[elem][edge])

        print(self.normalEdge)

        print(self.detJEdge)
        # quit()
        # ==================================
        # intialize curved element stuff
        # ==================================

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
        # print('internal', self.R)
        self.getEdgeResiduals()

        # return R
        # calculate analytical flux

    def getInteralResiduals(self):

        for idx_elem in range(self.mesh.nElem):

            for idx_basis in range(self.nBasis):
                Rtot = 0
                Uq = np.matmul(self.Phi_Elem, self.U[idx_elem])
                for q in range(self.nQuadPts2DElem):
                    flux = fluxes.analyticflux(Uq[q])
                    Rtot += np.dot(np.matmul(self.dPhi_dXi_Elem[q, idx_basis], self.invJ[idx_elem, q]), flux) *\
                        self.quadWts2DElem[q]*self.detJ[idx_elem, q]

                self.R[idx_elem, idx_basis] -= Rtot

        # return R

    def getEdgeResiduals(self):
        for idx_edge in range(self.mesh.nInEdge):
            # ! get the elements connected to the edge

            idx_elem_left = self.mesh.inEdge2Elem[idx_edge, 0]
            idx_edge_left = self.mesh.inEdge2Elem[idx_edge, 1]
            idx_elem_right = self.mesh.inEdge2Elem[idx_edge, 2]
            idx_edge_right = self.mesh.inEdge2Elem[idx_edge, 3]

            import ipdb
            ipdb.set_trace()
            uL = np.matmul(self.leftEdgePhi[idx_edge_left], self.U[idx_elem_left])
            uR = np.matmul(self.rightEdgePhi[idx_edge_right], self.U[idx_elem_right])

            for idx_basis in range(self.nBasis):
                # import ipdb
                # ipdb.set_trace()
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
            idx_elem = self.mesh.bcEdge2Elem[idx_edge, 0]
            idx_edge_left = self.mesh.bcEdge2Elem[idx_edge, 1]
            bc = self.mesh.bcEdge2Elem[idx_edge, 2]

            uL = np.matmul(self.leftEdgePhi[idx_edge_left], self.U[idx_elem])

            for idx_basis in range(self.nBasis):
                Rtot_left = np.zeros(self.nStates)
                Stot_left = 0

                for q in range(self.nQuadPts1DElem):

                    if any(bc == self.wall):
                        flux, s = fluxes.wallflux(uL[q], self.normalEdge[idx_elem, idx_edge][[q]])

                    elif any(bc == self.inlet):
                        flux, s = fluxes.inflowflux(
                            self.tempTotInf, self.Ptot_inf, 0.0, uL[q], self.normalEdge[idx_elem, idx_edge][[q]])
                    elif any(bc == self.outlet):
                        flux, s = fluxes.outflowflux(
                            self.P_inf, uL[q], self.normalEdge[idx_elem, idx_edge][[q]])
                    else:
                        print(bc)
                        raise NotImplementedError

                    # flux * delta L * wq
                    # import ipdb
                    # ipdb.set_trace()
                    tmp = flux*self.detJEdge[idx_elem, idx_edge, q] * self.quadWts1DElem[q]
                    # if any(bc == self.outlet):
                    # print('tmp', tmp)
                    # print('flux', flux)
                    # print('detJ', self.detJEdge[idx_elem, idx_edge, q])

                    Rtot_left += self.leftEdgePhi[idx_edge_left, q, idx_basis] * tmp
                    Stot_left += s*self.detJEdge[idx_elem, idx_edge, q]*self.quadWts1DElem[q]

                # print(bc, 'Rtot_left', Rtot_left)
                # print(self.R)

                self.R[idx_elem, idx_basis] += Rtot_left
                # print(self.R)
                # import ipdb
                # ipdb.set_trace()
                self.S[idx_elem] += Stot_left
            # print(self.R)
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
        # print(self.R)
        # quit()

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

    def solve(self, maxIter=100, tol=1e-7, cfl=0.9, freeStreamTest=False):

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

        # for each set of elements of the same geoemtric order
        with open(fileName + '.dat', 'w') as fid:

            fid.write('TITLE = "bump"\n')
            # fid.write('Variables="X", "Y", "U", "V", "rho", "P", "M", "rhoRes"\n')
            fid.write('Variables="X", "Y", "U", "V", "rho"\n')
            for q in self.mesh.elemOrder.keys():

                # get the basis functions for the mapping
                # this is a little bit of extra work for the linear elements, but it always us to
                # write the next loop without conditional statments (if q==1) which is nice and clean

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
                solPhi = np.zeros((len(Xi), self.nBasis))

                for idx, pt in enumerate(Xi):
                    # the basis function value at each of the quad points
                    geomPhi[idx], _ = basis(pt)
                    solPhi[idx], _ = self.basis(pt)

                # loop over elements now
                for idx_elem in self.mesh.elemOrder[q]:
                    # nodesPos = self.mesh.node2Pos[self.mesh.elem2Node[idx_elem]]
                    nodesPos = self.mesh.elem2HighOrderNode[q][idx_elem]
                    Upts = np.matmul(solPhi, self.U[idx_elem])
                    Xpts = np.matmul(geomPhi, nodesPos)

                    # writeZone()
                    print(idx_elem)
                    fid.write('ZONE\n')
                    fid.write('T="elem '+str(idx_elem)+'"\n')
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
    # bump = Mesh('meshes/bump0.gri')
    # bump = Mesh('meshes/test0_21.gri')

    def flatWall(x):
        return x*(x-1)*0.3
        # return 0

    test = Mesh('meshes/test0_21.gri', wallGeomFunc=flatWall)
    DGSolver = DGSolver(test, order=1)
    # DGprint(FVSolver.getResidual())
    DGSolver.solve()
    DGSolver.writeSolution('test2')
    # DGFVSolver.
    # DGFVSolver.solve()
    # DGFVSolver.postprocess()
    # DGFVSolver.writeSolution('lv0')
    # DGFVSolver.plotCP()
    # DGFVSolver.plotResiduals()

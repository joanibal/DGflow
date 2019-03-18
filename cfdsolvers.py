import numpy as np
import matplotlib.pyplot as plt
from meshes import Mesh
# import flux_lib
import time
# from solver import fvsolver as solver
# from solver import mesh as solver_mesh
# from solver import constants as solver_constants
# from solver import postprocess as solver_post


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

        # % quadrature points sufficient to integrate 2*p in 1D and 2D

        nQuadPts1D = order+1
        # these are 1d points
        quadPts1D, quadWts1D = quadrules.lgwt(nQuadPts1D, 0, 1)
        # these are 2d points
        # nQuadPts2D, quadPts2D, quadWts2D = quadrules.quad2d(nQuadPts1D, 0, 1)

        self.nBasis, self.basis = quadrules.getTriLagrangeBasis2D(order)
        # number of degrees of freedom
        # Ndof = nelem*nn

        self.initFreestream()

        # % calculate inverse mass matrix for one element
        Phi = np.zeros((nQuadPts2D, self.nBasis))
        dPhi_dXi = np.zeros((nQuadPts2D, self.nBasis, self.mesh.nDim))

        import ipdb
        ipdb.set_trace()
        for idx, pt in enumerate(quadPts2D):
            # the basis function value at each of the quad points
            Phi[idx], dPhi_dXi[idx] = self.basis(pt)

        self.invM = np.zeros((self.nBasis*self.nStates, self.nBasis*self.nStates))

        for k in range(self.mesh.nElem):
            M = mesh.area[k] * np.matmul(np.matmul(Phi.T, np.diag(quadWts2D)), Phi)
            self.invM[k] = np.linalg.inv(M)

        # M = h*h*Phi'* diag(wq)*Phi
        # iM = inv(M)

    def initFreestream(self):
        # calculate conversed qualities
        c = np.sqrt(self.gamma*self.R*self.tempInf)

        u = self.machInf*c
        Ub = np.array([self.rhoInf,
                       self.rhoInf*u,
                       0.0,
                       self.P_inf/(self.gamma-1) + 0.5*self.rhoInf*u**2])

        self.U = np.zeros((self.mesh.nElem, self.nBasis, self.nStates))

        self.U[:, :] = Ub
        print(self.U)

    # def initFromSolve(self, coarseMesh):

    #     for idx, _ in enumerate(coarseMesh.nElem):
    #         # set the four element that were created from uniform refinement
    #         self.U[:, 4*(idx-1):4*(idx)] = coarseMesh.U[idx]

    def getResidual(self):

        # loop over elements and compute residual contribution from interrior

        for elem in range(nElem):
            R

    def solve(self, maxIter=100000, tol=1e-7, cfl=1.0, freeStreamTest=False):

        solver.cfl = cfl
        solver_constants.mode = freeStreamTest
        t = time.time()
        if self.order == 1:
            self.Rmax, self.R2 = solver.solve_1storder(maxIter, tol)
        elif self.order == 2:
            solver_constants.recon_bc_flux = self.recon_bc_flux
            # solver.u = np.load('bump0_kfid_2_sol.npz')['q'].T

            self.Rmax, self.R2 = solver.solve_2ndorder(maxIter, tol)
        else:
            Exception('order not found')
        self.wallTime = time.time() - t

        self.U = solver.u.T
        # print(self.Rmax[:self.nIter])
        # import ipdb
        # ipdb.set_trace()
        self.Rmax = self.Rmax[self.Rmax > 0]
        self.nIter = len(self.Rmax)
        print('wall time', self.wallTime, 'iters', self.nIter, 'Rmax', self.Rmax[-1])

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

            fid.write('#U Data (Grid Variables must be nodal)\n')
            fid.write(str(self.U[:, 1]/self.U[:, 0])[1:-1]+'\n')

            fid.write('#V Data (Grid Variables must be nodal)\n')
            fid.write(str(self.U[:, 2]/self.U[:, 0])[1:-1]+'\n')

            fid.write('#rho Data (Grid Variables must be nodal)\n')
            fid.write(str(self.U[:, 0])[1:-1]+'\n')

            fid.write('#P Data (Grid Variables must be nodal)\n')
            fid.write(str(solver_post.p)[1:-1]+'\n')
            fid.write('#M Data (Grid Variables must be nodal)\n')
            fid.write(str(solver_post.mach)[1:-1]+'\n')
            fid.write('#rhoRes Data (Grid Variables must be nodal)\n')
            fid.write(str(solver.res[0, :])[1:-1]+'\n')

            fid.write('#number Data (Grid Variables must be nodal)\n')
            fid.write(str(np.arange(self.mesh.nElem)+1)[1:-1]+'\n')

            fid.write('#Connectivity List\n')
            for jj in range(self.mesh.nElem):
                fid.write(str(self.mesh.elem2Node[jj]+1)[1:-1]+'\n')

        # set np back to normal
        np.set_printoptions(threshold=1000)


if __name__ == '__main__':
    bump = Mesh('meshes/test.gri')

    FVSolver = DGSolver(bump, order=1)
    # FVSolver.solve()
    # FVSolver.postprocess()
    # FVSolver.writeSolution('lv0')
    # FVSolver.plotCP()
    # FVSolver.plotResiduals()

""" Mesh is the class that is used to hold information relavent to the mesh"""

from __future__ import print_function
from __future__ import division
# from readgri import readgri
import numpy as np
import copy
import matplotlib.pyplot as plt

import quadrules
# TODO
# add logger
# add error handeling


def getNormal(nodesEdgeXY, nodeXY):
    dX = nodesEdgeXY[0] - nodesEdgeXY[1]
    n = np.array([-1 * dX[1], dX[0]]) / np.linalg.norm(dX)

    # make sure the norm is pointing outward
    return np.sign(np.dot(n, nodesEdgeXY[0] - nodeXY)) * n


class Mesh(object):

    def __init__(self, gridFile=None, wallGeomFunc=None,  elem2Node=None, node2Pos=None, BCs=None, check=False):
        """
            initializes Mesh object from a gri mesh file


            Parameters
            ----------
            gridFile : str
                path the the mesh file
        """
        self.nDim = 2

        if gridFile:
            self.node2Pos, self.elem2Node, self.BCs = self.readGrid(
                gridFile)
        elif (elem2Node.any() and node2Pos.any() and BCs):
            self.elem2Node = elem2Node
            self.node2Pos = node2Pos
            self.BCs = BCs

        self.nElem = len(self.elem2Node)
        self.nNodes = len(self.node2Pos)
        self.nBCs = len(self.BCs)

        self.wallGeomFunc = wallGeomFunc
        # self.inEdge2Elem, self.inNormal, self.inLength, self.bcEdge2Elem, self.bcNormal, self.bcLength, self.area, self.cellCenter

        # for efficency this was written was one big ugly loop/function
        outputs = self.preprocess()
        self.inEdge2Elem, self.inNormal, self.inLength = outputs[:3]
        self.bcEdge2Elem, self.bcNormal, self.bcLength = outputs[3:6]
        self.area, self.cellCenter, self.elem2dX, self.wallEdges = outputs[6:]
        # self.elem2dX, self.cellCenter = self.getDistance2CellCenter()

        self.nInEdge = len(self.inEdge2Elem)
        self.nBCEdge = len(self.bcEdge2Elem)

        self.curvElem, self.linElem, self.elem2CurvElem = self.getElementOrder()
        self.nLinElem = len(self.linElem)
        self.nCurvElem = len(self.curvElem)

        if check:
            self.checkMesh()
        else:
            print('*mesh not checked*')

# IO stuff
    def readGrid(self, fileName):
        with open(fileName, 'r') as fid:

            curvOrder = 1
            if fileName[-3:] == 'gri':
                nNodes, nElem, nDim = [int(s) for s in fid.readline().split()]

                # read vertices
                V = np.array([[float(s) for s in fid.readline().split()]
                              for n in range(nNodes)])

                # read boundaries
                BCs = {}
                nB = int(fid.readline())
                for _ in range(nB):
                    s = fid.readline().split()
                    Nb = int(s[0])
                    BCname = s[2].lower()
                    Bi = np.array(
                        [[int(s) - 1 for s in fid.readline().split()] for n in range(Nb)])
                    # B.append(Bi[0])
                    BCs[BCname] = Bi[0]

                # read elements
                Ne0 = 0
                E = []
                linElem = np.array([], dtype=int)
                curvElem = np.array([], dtype=int)
                while (Ne0 < nElem):
                    s = fid.readline().split()
                    ne = int(s[0])
                    pe = int(s[1])
                    Ei = np.array(
                        [[int(s) - 1 for s in fid.readline().split()] for n in range(ne)])
                    E = Ei if (Ne0 == 0) else np.concatenate((E, Ei), axis=0)

                    # if pe > 1:
                    #     curvOrder = pe
                    #     curvElem = np.append(curvElem, np.arange(Ne0, Ne0 + ne, dtype=int))
                    # else:
                    #     linElem = np.append(linElem, np.arange(Ne0, Ne0 + ne, dtype=int))

                    Ne0 += ne

            elif fileName[-3:] == 'su2':
                nDim = int(fid.readline().split()[-1])
                nElem = int(fid.readline().split()[-1])

                E = -1*np.ones((nElem, 3), dtype=int)
                for ii in range(nElem):
                    E[ii] = fid.readline().split()[1:-1]

                nNodes = int(fid.readline().split()[-1])

                V = -1*np.ones((nNodes, 2))

                for ii in range(nNodes):
                    V[ii] = fid.readline().split()[:-1]

                # BCs
                nBCs = int(fid.readline().split()[-1])

                BCs = {}
                for ii in range(nBCs):
                    BCname = fid.readline().split()[-1].lower()
                    n = int(fid.readline().split()[-1])

                    b = np.array([fid.readline().split()[1:]], dtype=int)
                    for _ in range(n-1):
                        b = np.append(b, int(fid.readline().split()[-1]))

                    BCs[BCname] = b
            else:
                raise Exception
        return V, E, BCs

    def writeGrid(self, fileName, gridFormat='gri'):

        # if gridFormat == 'gri':
        with open(fileName, 'w') as fid:
            fid.write(str(self.nNodes) + ' ' + str(self.nElem) +
                      ' ' + str(self.nDim) + '\n')

            for ii in range(self.nNodes):
                arrStr = np.array2string(
                    self.node2Pos[ii], precision=5, separator=' ')[1:-1]
                fid.write(arrStr + '\n')

            fid.write(str(self.nBCs) + '\n')
            for BCName in self.BCs.keys():
                fid.write(
                    '1 ' + str(len(self.BCs[BCName])) + ' ' + BCName + '\n')
                fid.write(str(self.BCs[BCName]+1).replace('\n', '')[1:-1]+'\n')

            nElemGroup = 1
            order = [1]
            for ii in range(nElemGroup):
                fid.write(str(self.nElem) + ' ' +
                          str(order[ii]) + ' TriLagrange\n')
                for jj in range(self.nElem):
                    fid.write(str(self.elem2Node[jj]+1)[1:-1]+'\n')

# preprocessing
    def getElementOrder(self):
        """
            determines if an element should be linear or curved depending on if
            is touching a wall deseignated as curved
        """
        # get the element that has a curvWallEdge
        curvElem = self.bcEdge2Elem[self.wallEdges][:, 0]

        # all the reamaining elements are linear
        linElem = np.delete(range(self.nElem), curvElem)

        # create an inversee mapping of elem idx to curElem idx
        elem2CurvElem = np.ones(self.nElem, dtype=int)*-1
        for idx, elem in enumerate(curvElem):
            elem2CurvElem[elem] = idx

        return curvElem, linElem, elem2CurvElem

    def preprocess(self):
        node2Edge = {}
        I2E = []
        B2E = []
        Bn = []
        In = []
        bcLength = []
        inLength = []

        area = np.zeros(self.nElem)
        cellCenter = np.zeros((self.nElem, 2))
        elem2dX = np.zeros((self.nElem, 3, 2))

        wallEdges = []

        for idx_elem in range(self.nElem):
            nodes = self.elem2Node[idx_elem]

            # geometric information about the cell
            pos = self.node2Pos[nodes]
            vecs = pos - pos[0]
            area[idx_elem] = np.abs(np.cross(vecs[1], vecs[2])) / 2.0
            cellCenter[idx_elem] = np.mean(pos, axis=0)

            # loop over the edges of the cell, defined by the node not on the face
            for idx_edge in range(3):
                localNodeVec = np.delete(np.arange(3), idx_edge)
                edgeNodes = nodes[localNodeVec]

                # import pdb
                # pdb.set_trace()
                nodesEdgeXY = self.node2Pos[edgeNodes]
                edgeLength = np.linalg.norm(nodesEdgeXY[0] - nodesEdgeXY[1])
                edgeMidpoint = np.mean(nodesEdgeXY, axis=0)
                elem2dX[idx_elem, idx_edge, :] = edgeMidpoint - cellCenter[idx_elem]

                # check to see if both nodes on the edge are boundary nodes
                isBCNode = [set(edgeNodes).issubset(set(x))
                            for x in self.BCs.values()]

                if any(isBCNode):

                    # add the info to B2E
                    bcGroupIdx = isBCNode.index(True)

                    # info = np.array()
                    B2E.append([idx_elem, idx_edge, bcGroupIdx])

                    bn = getNormal(
                        self.node2Pos[edgeNodes], self.node2Pos[nodes[idx_edge]])
                    Bn.append(bn)

                    bcLength.append(edgeLength)

                    # add to the list of wall edges if it is indeed a wall edge\
                    if 'curv' in self.BCs.keys()[bcGroupIdx]:
                        wallEdges.append(len(B2E) - 1)
                    continue  # if the nodes are on a BC don't add it to the list

                # information about the edge used for the I2E array
                info = np.array([idx_elem, idx_edge])

                if frozenset(edgeNodes) in node2Edge:
                    edgeIdx = node2Edge[frozenset(edgeNodes)]

                    if idx_elem > I2E[edgeIdx][0]:
                        I2E[edgeIdx] = np.hstack((I2E[edgeIdx], info))
                    elif idx_elem < I2E[edgeIdx][0]:
                        I2E[edgeIdx] = np.hstack((info, I2E[edgeIdx]))
                        In[edgeIdx] *= -1

                    else:
                        raise Exception
                else:
                    node2Edge[frozenset(edgeNodes)] = len(I2E)
                    I2E.append(info)
                    i_n = getNormal(
                        self.node2Pos[edgeNodes], self.node2Pos[nodes[idx_edge]])
                    In.append(i_n)
                    inLength.append(np.linalg.norm(
                        nodesEdgeXY[0] - nodesEdgeXY[1]))

        return np.array(I2E), np.array(In), np.array(inLength), np.array(B2E), np.array(Bn), np.array(bcLength), area, cellCenter, elem2dX, np.array(wallEdges)

    def getDistance2CellCenter(self):
        """
            reuturns the distance from the midpoint of each edge to the cell center
            this is used in the second order FVM solver

            returns

                elem2dX: numpy array (nElem, 3, 2 )

        """
        elem2dX = np.zeros((self.nElem, 3, 2))
        cellCenter = np.zeros((self.nElem, 2))
        for idx_elem in range(self.nElem):
            nodes = self.elem2Node[idx_elem]
            pos = self.node2Pos[nodes]
            cellCenter[idx_elem] = np.mean(pos, axis=0)
            # nodes = set(nodes)
            for jj, node in enumerate(nodes):
                edgeNodes = nodes - set([node])
                # import pdb
                # pdb.set_trace()
                nodesEdgeXY = self.node2Pos[list(edgeNodes)]
                edgeMidpoint = np.mean(nodesEdgeXY, axis=0)
                # import ipdb
                # ipdb.set_trace()
                elem2dX[idx_elem, jj, :] = edgeMidpoint - cellCenter[idx_elem]
                # print('distElem', edgeMidpoint, cellCenter, 'edge', idx_face)

        return elem2dX, cellCenter

    def checkMesh(self, tol=1e-8):
        """
            Checks the internal consistence of the mesh by taking the sum of the edge length
            in each direction for each element
        """
        # check that all the normals are unit
        if self.inNormal.size:  # check to make sure we have interrior edges
            np.testing.assert_almost_equal(np.linalg.norm(
                self.inNormal, axis=1), np.ones(len(self.inNormal)))

        np.testing.assert_almost_equal(np.linalg.norm(
            self.bcNormal, axis=1), np.ones(len(self.bcNormal)))

        elemRes = np.zeros((self.elem2Node.shape[0], self.nDim))

        for edgeIdx, IE in enumerate(self.inEdge2Elem):
            length = self.inLength[edgeIdx]
            # import ipdb
            # ipdb.set_trace()
            elemRes[IE[0]] += self.inNormal[edgeIdx]*length
            elemRes[IE[2]] -= self.inNormal[edgeIdx]*length

        for edgeIdx, BE in enumerate(self.bcEdge2Elem):
            length = self.bcLength[edgeIdx]
            elemRes[BE[0]] += self.bcNormal[edgeIdx]*length

        if np.sum(sum(np.abs(elemRes))) > tol:
            print('elemRes')
            print(elemRes)
            print(np.max(np.abs(elemRes)))

            raise Exception
        else:
            print('passed mesh check ', np.max(np.abs(elemRes)))

# modification

    def refineElement(self, lv, idx_elem, geomOrder=1):
        """
            used for ploting the solution
            refine a signle element(don't alter the grid though) and return the node list and connectivity
        """
        nodes = self.elem2Node[idx_elem]
        # nodes = set(nodes)
        newNodesPos = np.zeros((3, 2))
        newNodes = np.arange(3, 6)
        elem2Node = np.zeros((4, 3))

        for _ in range(lv):
            for idx_edge in range(3):
                localNodeVec = np.delete(np.arange(3), idx_edge)
                edgeNodes = nodes[localNodeVec]

                newNodesPos[idx_edge] = np.mean(self.node2Pos[edgeNodes], axis=0)

                # check to see if both nodes on the edge are boundary nodes
                isBCNode = [set(edgeNodes).issubset(set(x))
                            for x in self.BCs.values()]
                if any(isBCNode):

                    # add the info to B2W
                    bcGroupName = self.BCs.keys()[isBCNode.index(True)]

                    if self.wallGeomFunc and 'wall' in bcGroupName:
                        # use the analytic function for the wall geometry to snap the point to the wall
                        newNodesPos[idx_edge][1] = self.wallGeomFunc(newNodesPos[idx_edge][0])

            # now we have calculated all the new nodes, loop over the faces and create the elems
            for idx_edge, n in enumerate(nodes):
                temp = copy.copy(newNodes)

                temp[idx_edge] = n
                elem2Node[idx_edge] = temp[::-1]

            # add another node made from all the new nodes
            elem2Node[3] = newNodes

        print(newNodesPos)
        print(elem2Node)

    def refine(self):
        """
            splits every element into 4 elements
        """
        old2NewNodes = {}
        elem2Node = np.ones((self.nElem*4, 3), dtype=int)

        BCs = self.BCs

        nodeIdx = len(self.node2Pos) - 1

        newNodes = np.empty(3, dtype=int)
        node2Pos = self.node2Pos

        for idx_elem in range(self.nElem):
            nodes = self.elem2Node[idx_elem]
            # nodes = set(nodes)
            for idx_edge in range(3):
                localNodeVec = np.delete(np.arange(3), idx_edge)
                edgeNodes = nodes[localNodeVec]

                # ndary nodes for that boundary group
                if frozenset(edgeNodes) in old2NewNodes:

                    newEdgeNodes = old2NewNodes[frozenset(edgeNodes)]
                    nodeIdx = newEdgeNodes[1]

                    # else:
                    #     raise Exception
                else:
                    nodeIdx = len(node2Pos)
                    newEdgeNodes = np.insert(list(edgeNodes), 1, nodeIdx)
                    old2NewNodes[frozenset(edgeNodes)] = newEdgeNodes

                    node2Pos = np.vstack(
                        (node2Pos, np.mean(node2Pos[newEdgeNodes[[0, 2]]], axis=0)))

                # check to see if both nodes on the edge are boundary nodes
                isBCNode = [set(edgeNodes).issubset(set(x))
                            for x in self.BCs.values()]
                if any(isBCNode):

                    # add the info to B2W
                    bcGroupName = BCs.keys()[isBCNode.index(True)]
                    BCs[bcGroupName] = np.append(BCs[bcGroupName], nodeIdx)

                    if self.wallGeomFunc and 'curv' in bcGroupName:
                        # use the analytic function for the wall geometry to snap the point to the wall
                        node2Pos[nodeIdx][1] = self.wallGeomFunc(node2Pos[nodeIdx][0])

                newNodes[idx_edge] = newEdgeNodes[1]
                # add node numbers to the list of bou

            # now we have calculated all the new nodes, loop over the faces and create the elems
            for idx_edge, n in enumerate(nodes):
                temp = copy.copy(newNodes)

                temp[idx_edge] = n
                elem2Node[idx_elem*4+idx_edge] = temp[:: -1]

            # add another node made from all the new nodes
            elem2Node[idx_elem*4+3] = newNodes

        self.__init__(elem2Node=elem2Node, node2Pos=node2Pos, BCs=BCs,
                      wallGeomFunc=self.wallGeomFunc, check=True)

    def plot(self, fileName=None, geomOrder=1):
        plt.figure(figsize=(12, 4))

        bc_colors = ['-b', '-g', '-m', '-r']

        for IE in self.inEdge2Elem:
            localNodeVec = np.delete(np.arange(3),  IE[1])
            edgeNodes = self.elem2Node[IE[0]][localNodeVec]
            nodes = self.node2Pos[edgeNodes]

            plt.plot(nodes[:, 0], nodes[:, 1], '-k')

        for BE in self.bcEdge2Elem:
            localNodeVec = np.delete(np.arange(3),  BE[1])
            edgeNodes = self.elem2Node[BE[0]][localNodeVec]
            nodes = self.node2Pos[edgeNodes]

            plt.plot(nodes[:, 0], nodes[:, 1], bc_colors[BE[2]])

        plt.axis('equal')

        if fileName:
            plt.axis('off')
            plt.savefig(fileName, bbox_inches='tight')

        # print(self.BCs)
        # print(bc_colors)

    def plotElem(self, elemIdxs):
        for edgeIdx in elemIdxs:

            nodes = self.node2Pos[self.elem2Node[edgeIdx]]
            nodes = np.vstack((nodes, nodes[0]))
            plt.plot(nodes[:, 0], nodes[:, 1], 'k')

    def plotElemOrder(self):

        for idx_elem in self.linElem:

            nodes = self.node2Pos[self.elem2Node[idx_elem]]
            nodes = np.vstack((nodes, nodes[0]))
            plt.plot(nodes[:, 0], nodes[:, 1], 'k')

        for idx_elem in self.curvElem:

            nodes = self.node2Pos[self.elem2Node[idx_elem]]
            nodes = np.vstack((nodes, nodes[0]))
            plt.plot(nodes[:, 0], nodes[:, 1], 'b')

        plt.show()

    def getLinearJacobian(self):
        nElem = len(self.linElem)

        detJ = np.zeros(nElem)
        invJ = np.zeros((nElem, self.nDim, self.nDim))

        # if 1 in self.elemOrder:
        for idx_linElem, idx_elem in enumerate(self.linElem):
            nodes = self.node2Pos[self.elem2Node[idx_elem]]

            # from eqn 4.3.8 in notes
            J = np.array([[nodes[1][0] - nodes[0][0], nodes[2][0] - nodes[0][0]],
                          [nodes[1][1] - nodes[0][1], nodes[2][1] - nodes[0][1]]])
            detJ[idx_linElem] = np.linalg.det(J)
            invJ[idx_linElem] = np.linalg.inv(J)

        return invJ, detJ

    def getHighOrderNodes(self):
        """
        finds the location of the high order nodes on the mesh
        """

        # for q in highOrder:
        nHiOrderElem = len(self.curvElem)
        q = self.curvOrder

        Xi = quadrules.getTriLagrangePts2D(q)
        N = len(Xi)

        self.curvNodes = np.zeros((nHiOrderElem, N, 2))

        self.elemIdx2HiOrderElemIdx = np.ones(self.nElem, dtype=int)*10**10

        for idx_elem, elem in enumerate(self.curvElem):
            nodes = self.node2Pos[self.elem2Node[elem]]

            J = np.array([[nodes[1][0] - nodes[0][0], nodes[2][0] - nodes[0][0]],
                          [nodes[1][1] - nodes[0][1], nodes[2][1] - nodes[0][1]]])

            # map the nodes to physcial space using a linear jacobain
            self.curvNodes[idx_elem] = nodes[0] + np.matmul(J, Xi.T).T
            self.elemIdx2HiOrderElemIdx[elem] = idx_elem

        # print(self.curvNodes[0])

        # loop over wall nodes and snap nodes on wall face to edge and move nodes inside up slightly
        for edge in self.bcEdge2Elem[self.wallEdges]:
            elem = edge[0]
            idx_edge = edge[1]
            idx_hiElem = self.elemIdx2HiOrderElemIdx[elem]
            idx_face = edge[1]

            # plt.plot(self.curvNodes[idx_hiElem][:, 0],
            #          self.curvNodes[idx_hiElem][:, 1], '-o')

            # depending on which face is on the wall different nodes need to be adjusted
            nodes2Move = []
            if idx_face == 0:
                nodeIdxs = [2*q]
                for ii in range(q-1, 1, -1):
                    nodeIdxs.append(nodeIdxs[-1] + ii)

                nodes2Move.append(nodeIdxs)
                for ii in range(q - 2, 0, -1):
                    nodes2Move.append([x - 1 for x in nodes2Move[-1][:-1]])

            elif idx_face == 1:
                nodeIdxs = [1 + q]
                for ii in range(q-1, 1, -1):
                    nodeIdxs.append(nodeIdxs[-1] + ii + 1)

                nodes2Move.append(nodeIdxs)
                for ii in range(q - 2, 0, -1):
                    nodes2Move.append([x + 1 for x in nodes2Move[-1][:-1]])

            # different
            elif idx_face == 2:
                nodes2Move.append(range(1, q))
                for ii in range(q+1, 3, -1):
                    nodes2Move.append([x + ii for x in nodes2Move[-1][:-1]])

            else:
                raise NotImplementedError

            # snap wall elements y position
            for idx_row, fact in enumerate(np.linspace(1, 0, q)[:-1]):
                for idx_node in nodes2Move[idx_row]:

                    # determine the diplacement from a linear
                    nodes = self.elem2Node[elem]
                    # import ipdb
                    # ipdb.set_trace()
                    localNodeVec = np.delete(np.arange(3), idx_edge)
                    edgeNodes = self.node2Pos[nodes[localNodeVec]]

                    x = self.curvNodes[idx_hiElem][idx_node][0]
                    y = (x - edgeNodes[0, 0])*(edgeNodes[1, 1] - edgeNodes[0, 1]) / \
                        (edgeNodes[1, 0] - edgeNodes[0, 0]) + edgeNodes[0, 1]
                    y_snap = self.wallGeomFunc(x)
                    # print(idx_row, idx_node, fact)
                    self.curvNodes[idx_hiElem][idx_node][1] += fact*(y_snap - y)

            # plt.plot(self.curvNodes[idx_hiElem][:, 0],
            #          self.curvNodes[idx_hiElem][:, 1], '-o')
            # plt.show()

        # print(self.curvNodes[0])
        # print(idx_face, nodes2Move)
        # plt.xlabel('X')
        # plt.ylabel('Y')
        # plt.legend(['Unperturbed', 'Perturbed'])
        # quit()
        # print(self.curvNodes[q][idx_hiElem][nodes2Move[1]])
        # print(self.curvNodes[q][idx_hiElem][nodes2Move[2]])

    def getCurvedJacobian(self, elem, quadPts, basis):
        """
        returns the value of the jacobian at each of the quadrature points evaluated using the supplied basis functions
        """
        # nBasis, basis = quadrules.getTriLagrangeBasis2D(q)
        # _, quadWts = quadrules.getQuadPts1D(q+1, 0, 1)

        # quadPts2D = quadrules.getQuadPts2D(q+1)
        # for pt in quadPts2D:

        J = np.zeros((len(quadPts), 2, 2))
        invJ = np.zeros((len(quadPts), 2, 2))
        detJ = np.zeros(len(quadPts))
        # tang_vec = np.zeros((len(quadPts), 2))
        # n = np.zeros((len(quadPts), 2))
        for idx, pt in enumerate(quadPts):
            _, gphi = basis(pt)
            J[idx] = np.matmul(self.curvNodes[elem].T, gphi[0])
            invJ[idx] = np.linalg.inv(J[idx])
            detJ[idx] = np.linalg.det(J[idx])
            # print(J[idx])
            # tang_vec = -J[idx][:, 0] + J[idx][:, 1]
            # print(J[idx], pt)

            # # the result of crossing the trangent vector with k hat
            # n[idx] = np.array([tang_vec[0], -tang_vec[1]])
            # # quit()
            # # print

            # arclength += np.linalg.norm(t_ds_dsig*quadWts[idx])
            # print(np.linalg.norm(t_ds_dsig))

        # print(arclength, np.sum(quadWts))

        return J, invJ, detJ
        # for elem in range(self.elemOrder[q]):

    def getEdgeJacobain(self, quadPts1D, basis):
        nElem = len(self.curvElem)
        nQuadPts1D = len(quadPts1D)

        detJEdge = np.zeros((nElem,  3, nQuadPts1D))
        normalEdge = np.zeros((nElem,  3, nQuadPts1D, 2))

        for idx_elem, elem in enumerate(self.curvElem):
            for edge in range(3):
                pts = np.zeros((nQuadPts1D, 2))
                if edge == 0:
                    pts[:, 0] = 1 - quadPts1D
                    pts[:, 1] = quadPts1D
                    dXi_dX = np.array([-1, 1])

                elif edge == 1:
                    pts[:, 0] = 0
                    pts[:, 1] = 1 - quadPts1D
                    dXi_dX = np.array([0, -1])

                elif edge == 2:
                    pts[:, 0] = quadPts1D
                    pts[:, 1] = 0
                    dXi_dX = np.array([1, 0])

                J = self.getCurvedJacobian(idx_elem, pts, basis)[0]

                for q in range(len(J)):
                    # import ipdb
                    # ipdb.set_trace()
                    tang_vec = J[q][:, 0]*dXi_dX[0] + J[q][:, 1]*dXi_dX[1]
                    normalEdge[idx_elem][edge][q] = np.array([tang_vec[1], -tang_vec[0]])

                    detJEdge[idx_elem][edge][q] = np.linalg.norm(
                        normalEdge[idx_elem][edge][q])
                    normalEdge[idx_elem][edge][q] /= detJEdge[idx_elem][edge][q]

                # print(edge, self.normalEdge[elem][edge])
        return detJEdge, normalEdge

        # def plotBCs(self):
        #     for BCname in self.BCs.keys():
if __name__ == '__main__':

    def flatWall(x):
        return -x*(x-1)*0.2

    # test = Mesh('meshes/test0_2.gri', wallGeomFunc=flatWall)
    test = Mesh('meshes/bump0_kfid.gri')
    test.refine()
    test.refine()
    test.plotElemOrder()
    # test.getHighOrderNodes()
    # test.refineElement(0)
    # pt = np.array([0, 5])

    # test.getCurvedJacobian(2, pt)
    # test.getLinearJacobian()
    # print(test.inEdge2Elem)
    # print(test.bcEdge2Elem)
    # print(test.elem2Node)
    # print(test.node2Pos)
    # # test.refine()
    # # test.refine()
    # test.plot()
    # plt.show()

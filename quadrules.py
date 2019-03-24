""" python implementation of lgwt.m and qaud2d.m"""
import numpy as np


def getQuadPts1D(N, a, b):
    """
    adapted from lwgt.m


    This script is for computing definite integrals using Legendre-Gauss
    Quadrature. Computes the Legendre-Gauss nodes and weights  on an interval
    [a,b] with truncation order N

    Suppose you have a continuous function f(x) which is defined on [a,b]
    which you can evaluate at any x in [a,b]. Simply evaluate it at all of
    the values contained in the x vector to obtain a vector f. Then compute
    the definite integral using sum(f.*w)


    Orignially written in matlab by Greg von Winckel - 02/25/2004


    Parameters
    ----------
    N : int
        number of evaluation points
    a : float
        lower interval bound
    b : float
        upper interval bound
    """
    N = N-1
    N1 = N+1
    N2 = N+2

    xu = np.linspace(-1, 1, N1)

    # Initial guess
    y = np.cos((2*np.arange(N+1) + 1)*np.pi/(2*N+2)) + (0.27 / N1) * np.sin(np.pi * xu * N/N2)

    # Legendre-Gauss Vandermonde Matrix
    L = np.zeros((N1, N2))

    # Compute the zeros of the N+1 Legendre Polynomial
    # using the recursion relation and the Newton-Raphson method
    y0 = 2

    # Iterate until new points are uniformly within epsilon of old points
    while np.max(np.abs(y-y0)) > np.finfo(float).eps:

        L[:, 0] = 1
        L[:, 1] = y

        for k in np.arange(2, N1+1):
            L[:, k] = ((2*k-1)*y*L[:, k-1]-(k-1)*L[:, k-2])/k

        Lp = (N2)*(L[:, N1-1] - y*L[:, N2-1])/(1-y**2)

        y0 = y
        y = y0 - L[:, N2-1]/Lp

    # Linear map from[-1, 1] to[a, b]
    x = (a*(1-y)+b*(1+y))/2
    # Compute the weights
    w = (b-a)/((1-y ** 2)*Lp ** 2)*(float(N2)/N1)**2

    return x[::-1], w[::-1]


def getQuadPtsTri(order):
    """
    from quad2d.c

    Dunavant points generated with .m code written by John Burkard

        http://people.scs.fsu.edu/~burkardt/f_src/dunavant/dunavant.html

    1. David Dunavant,
        High Degree Efficient Symmetrical Gaussian Quadrature Rules for the Triangle,
        International Journal for Numerical Methods in Engineering,
        Volume 21, 1985, pages 1129-1148.
    2. James Lyness, Dennis Jespersen,
        Moderate Degree Symmetric Quadrature Rules for the Triangle,
        Journal of the Institute of Mathematics and its Applications,
        Volume 15, Number 1, February 1975, pages 19-32.
    """

    if order == 0 or order == 1:
        quadPts = [[0.333333333333333, 0.333333333333333]]
        quadWts = [0.500000000000000]

    elif order == 2:
        quadPts = [[0.666666666666667, 0.166666666666667],
                   [0.166666666666667, 0.166666666666667],
                   [0.166666666666667, 0.666666666666667]]
        quadWts = [0.166666666666666,
                   0.166666666666666,
                   0.166666666666666]

    elif order == 3:
        quadPts = [[0.333333333333333, 0.333333333333333],
                   [0.600000000000000, 0.200000000000000],
                   [0.200000000000000, 0.200000000000000],
                   [0.200000000000000, 0.600000000000000]]
        quadWts = [-0.281250000000000,
                   0.260416666666667,
                   0.260416666666667,
                   0.260416666666667]

    elif order == 4:
        quadPts = [[0.108103018168070, 0.445948490915965],
                   [0.445948490915965, 0.445948490915965],
                   [0.445948490915965, 0.108103018168070],
                   [0.816847572980459, 0.091576213509771],
                   [0.091576213509771, 0.091576213509771],
                   [0.091576213509771, 0.816847572980459]]
        quadWts = [0.111690794839005,
                   0.111690794839005,
                   0.111690794839005,
                   0.054975871827661,
                   0.054975871827661,
                   0.054975871827661]

    elif order == 5:
        quadPts = [[0.333333333333333, 0.333333333333333],
                   [0.059715871789770, 0.470142064105115],
                   [0.470142064105115, 0.470142064105115],
                   [0.470142064105115, 0.059715871789770],
                   [0.797426985353087, 0.101286507323456],
                   [0.101286507323456, 0.101286507323456],
                   [0.101286507323456, 0.797426985353087]]
        quadWts = [0.112500000000000,
                   0.066197076394253,
                   0.066197076394253,
                   0.066197076394253,
                   0.062969590272414,
                   0.062969590272414,
                   0.062969590272414]

    else:
        raise NotImplementedError

    return np.array(quadPts), np.array(quadWts)


def getTriLagrangeBasis2D(p):
    # % calculates coeffs for full-order Lagrange basis of order p
    # % reference element is a unit isoceles right triangle

    xi = np.linspace(0, 1, p+1)
    eta = xi

    N = (p+1)*(p+2)/2
    # % number of basis functions
    A = np.zeros((N, N))
    C = A

    i = 0
    # % build A-matrix
    for iy in range(p+1):
        for ix in range(p-iy+1):  # loop over nodes
            k = 0

            for s in range(p+1):
                for r in range(p-s+1):  # % loop over monomials
                    A[i, k] = xi[ix]**r * eta[iy]**s
                    k = k+1

            i = i + 1

    C = np.linalg.inv(A)
    C = C.T

    # compute dPhi_dxi and dPhi_deta
    # hard to explain but iw you work through it by hand this is the pattern emerges
    Cx = np.zeros(C.shape)
    Cy = np.zeros(C.shape)
    kk = 0
    ll = p+1
    for ii in range(1, p+1):
        # print(ii, p+1-ii)
        for jj in range(1, p+2-ii):
            # print(jj, ii, kk, ll)
            Cx[:, kk] = jj*C[:, kk+1]
            Cy[:, kk] = ii*C[:, ll]
            kk += 1
            ll += 1
        # print('zeros', kk)

        Cx[:, kk] = np.zeros(N)
        kk += 1

    def basis(Xi):
        # evaluates the jth basis function at (xi, eta) in reference space
        X = np.zeros(N)
        # Phi = np.zeros(N)
        idx = 0
        for s in range(p+1):
            for r in range(p-s+1):  # % loop over monomials
                X[idx] = Xi[0]**r * Xi[1]**s
                idx += 1

        vals = np.sum(C*X, axis=1)
        gradXi = np.sum(Cx*X, axis=1)
        gradEta = np.sum(Cy*X, axis=1)
        # import ipdb
        # ipdb.set_trace()
        return vals, np.dstack((gradXi, gradEta))

    return N, basis


def getTriLagrangePts2D(p):
    xi = np.linspace(0, 1, p+1)
    eta = xi
    N = (p+1)*(p+2)//2

    Xi = np.zeros((N, 2))
    idx = 0
    for iy in range(p+1):
        for ix in range(p-iy+1):  # loop over nodes
            Xi[idx] = np.array([xi[ix], eta[iy]])
            idx += 1
    return Xi


if __name__ == "__main__":
    b = getTriLagrangeBasis2D(p=2)
    print(getTriLagrangePts2D(p=2))

    # print(b([1, 0]))
    # print(b([0, 0]))
    # print(b([0, 1]))
    # print(b([0.5, 0.5]))
    # print(b([0, 0.5]))
    # print(b([0.5, 0]))

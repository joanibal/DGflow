""" python implementation of lgwt.m"""
import numpy as np


def lgwt(N, a, b):
    """

    lgwt.m

    This script is for computing definite integrals using Legendre-Gauss
    Quadrature. Computes the Legendre-Gauss nodes and weights  on an interval
    [a,b] with truncation order N

    Suppose you have a continuous function f(x) which is defined on [a,b]
    which you can evaluate at any x in [a,b]. Simply evaluate it at all of
    the values contained in the x vector to obtain a vector f. Then compute
    the definite integral using sum(f.*w);

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
        # for ii in range(2):

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

    return x, w


def quad2d(nq1, a, b):
    """
    nq1^2 quadrature points on [a,b]^2

    """

    [sq, wq1] = lgwt(nq1, a, b)   # these are 1d points
    nq = nq1*nq1
    xyq = np.zeros((nq, 2))
    wq = np.zeros((nq, 1))
    k = 0
    for j in range(nq1):
        for i in range(nq1):
            xyq[k, :] = np.array([sq[i], sq[j]])
            wq[k] = wq1[i]*wq1[j]
            k = k+1

    return nq, xyq, wq


print(quad2d(3, 0, 1))

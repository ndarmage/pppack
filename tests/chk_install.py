import numpy as np
import itertools as itt
import pppack as ppk


def cmpgtau(tau, g):
    Ni = [len(t) for t in tau]
    N = np.prod(Ni)
    gtau = np.zeros((N,), order='F')
    for x, i in zip(itt.product(*tau), itt.product(*[range(nx) for nx in Ni])):
        gtau[ppk.getidx(i, Ni)] = g(*x)
    return gtau


def chkapprx(tau, f, g):
    for x in itt.product(*tau):
        fx, gx = f(x), g(*x)
        print("%+13.6e %+13.6e %+7.2f" % (fx, gx, (1. - fx / gx) * 100))


def g3(x, y, z):
    return np.sin(x) / (x + 1) + (z + 1) * y**2 + x * np.exp(z) + .25


def g2(x, y):
    return g3(x, y, 0)


def g1(x):
    return g2(x, 0)


def run_test():
    k, t, xi = 4, np.zeros((17,), order='F'), [0., 1.]
    n = len(t) - k
    t[0:k] = 0.
    t[n:] = 10.
    for i in range(1, 10):
        t[k + i - 1] = float(i)
    S1 = ppk.S(k, t)

    P1 = ppk.Fd(ppk.PP(3, xi), name="quadratic")
    P2 = ppk.PP(2, xi)
    P2.name = "linear"
    V1 = P1 * P2
    V1.name = "poly"
    S2 = ppk.S(k, t)  # like a deepcopy of S1
    V2 = ppk.Fd(S1, S2, name="2D-bspline")
    V3 = ppk.Fd(P1, P2, ppk.PP(2, xi, name="linear2"), name="3D-poly")
    V4 = ppk.Fd(P1, S2, P2, name="3D-mix")

    print(" --- test 1D poly ---")
    f1 = ppk.fd.inFd(P1)
    tau = [0., .5, 1.]  # quadratic
    gtau = cmpgtau(tau, g1)
    f1.cmpcoef(tau, gtau)
    chkapprx(tau, f1, g1)

    print(" --- test 2D poly ---")
    f2 = ppk.fd.inFd(V1)
    tau.append([0., 1.])  # linear
    gtau = cmpgtau(tau, g2)
    f2.cmpcoef(tau, gtau)
    chkapprx(tau, f2, g2)

    print(" --- test 2D Bspline ---")
    h2 = ppk.fd.inFd(V2)
    tau = []
    for V in V2:
        tau.append(V.cmptau())
    gtau = cmpgtau(tau, g2)
    h2.cmpcoef(tau, gtau)
    chkapprx(tau, h2, g2)

    print(" --- test 3D poly ---")
    f3 = ppk.fd.inFd(V3)
    tau = [0., .5, 1.]  # quadratic
    tau.append([0., 1.])  # linear
    tau.append([0., 1.])  # linear
    gtau = cmpgtau(tau, g3)
    f3.cmpcoef(tau, gtau)
    chkapprx(tau, f3, g3)

    print(" --- test 3D mix ---")
    f4 = ppk.fd.inFd(V4)
    tau = [0., .5, 1.]  # quadratic
    tau.append(S2.cmptau())
    tau.append([0., 1.])  # linear
    gtau = cmpgtau(tau, g3)
    f4.cmpcoef(tau, gtau)
    chkapprx(tau, f4, g3)


if __name__ == "__main__":

    run_test()

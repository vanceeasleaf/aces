import numpy as np
import aces.tools as tl
import h5py
from numpy.linalg import norm
import os
from ase import io
from scipy import optimize
from aces.f import read_forces, writefc2, writefc3, disp2atoms


def shrink(y, a):
    return np.sign(y) * np.maximum(np.abs(y) - a, 0.0)


def dot(a, b):
    return np.tensordot(a, b, axes=([1], [0]))


def maxeig(A):
    B = np.zeros_like(A)
    for j in A.shape[1]:
        B[:, j] = A[:, j] / A.sum(axis=0)[j]
    C = B.sum(axis=1)
    W = C / C.sum()
    lmax = (A.dot(W) / W).average()
    return lmax


class runner:

    def __init__(self, NAH=3, split=True, mu=0.1, lam=0.9):
        self.NAH = NAH
        self.split = split
        self.mu = mu
        self.lam = lam
        # self.db=h5py.File('force.hdf5')

    def getForce(self, pos, files):
        print("reading vasprun.xml and POSCAR")
        u = []
        for file in files:
            dir = os.path.dirname(file)
            atoms = io.read(dir + '/POSCAR')
            u.append(atoms.positions - pos)
        forces = []
        for file in files:
            forces.append(read_forces(file))
        return np.array(forces), np.array(u)

    def getsatoms(self):
        filename = 'disp_fc3.yaml'
        if (tl.exists(filename)):
            return disp2atoms(filename)
        filename = 'disp.yaml'
        if (tl.exists(filename)):
            return disp2atoms(filename)
        filename = '3RD.SPOSCAR'
        if (tl.exists(filename)):
            from ase import io
            return io.read(filename, format='vasp')
        filename = 'SPOSCAR'
        if (tl.exists(filename)):
            from ase import io
            return io.read(filename, format='vasp')

    def getSupercell(self, atoms):
        from pyspglib import spglib
        s = spglib.get_symmetry(atoms)
        symmetry = []
        print("building symmetry")
        for i, rot in enumerate(s['rotations'][:100]):
            print("symetry :", i)
            trans = s['translations'][i]
            map0 = self.getMap(atoms, rot, trans)
            symmetry.append([rot, map0])
        return symmetry

    def getMap(self, atoms, rot, trans):
        v = atoms.copy()
        v.positions = v.positions.dot(rot.T)
        v.translate(trans.dot(v.cell))
        import itertools
        from scipy.spatial.distance import cdist
        posi = atoms.positions
        d2s = np.empty((27, len(v), len(v)))
        for j, (ja, jb, jc) in enumerate(
                itertools.product(range(-1, 2), range(-1, 2), range(-1, 2))):
            posj = v.positions + np.dot([ja, jb, jc], v.cell)
            d2s[j, :, :] = cdist(posi, posj, "sqeuclidean")
        d2min = d2s.min(axis=0)
        map0 = np.argmin(d2min, axis=1)
        return map0

    def getTrainSets(self, u):
        assert len(u) > 0
        self.L = len(u)
        n = self.natom = len(u[0])
        row = 0
        rowr = [0]
        for i in range(self.NAH):
            row += (n * 3)**i
            rowr.append(row)
        self.rowr = rowr

    def getMatrix(self, F, u):

        print("getting compressive matrix")
        rowr = self.rowr
        A = np.zeros([self.L, rowr[-1]])
        g = self.mulU
        # shape=F.shape
        n = self.natom
        for j in range(self.L):
            for i in range(self.NAH):
                r = range(rowr[i], rowr[i + 1])
                A[j, r] = -g(u[j].flatten(), i)
        F = F.reshape([self.L, 3 * n])
        c = 3 * n
        F = F.T.flatten().T
        A = np.kron(np.eye(c), A)
        return F, A

    def gauss(a):
        m, n = a.shape
        b = np.zeros(n)
        for i in range(0, n - 1):
            for j in range(i + 1, n):
                imax = np.abs(a[i:n, i]).maxarg()
                if imax != i:
                    a[i], a[imax] = a[imax], a[i]
                if a[j, i] != 0.0 and a[i, i] != 0.0:
                    lam = float(a[j, i]) / a[i, i]
                    a[j] = a[j] - lam * a[i]
        for k in range(n - 1, -1, -1):
            b[k] = (b[k] - np.dot(a[k, (k + 1):], b[(k + 1):])) / a[k, k]
        result = b
        return result

    def mulU(self, x, p):
        if p > 0:
            return np.kron(self.mulU(x, p - 1), x) / p
        else:
            return 1.0

    def getCsMat(self, F, u, symmetry):
        self.getTrainSets(u)

        # keep to be the constrain of the newest variables
        Q = []
        n = u.shape[1]
        v = self.rowr[-1]
        p = 3 * n
        step = 0
        nval = p * v

        # the connection between oldest variables and newest variables
        E = None
        for rot, map0 in symmetry:
            step += 1
            print("step:", step)
            for i in range(n):
                print("atom:", i)
                for j in range(n):
                    ii = map0[i]
                    jj = map0[j]
                    for a in range(3):
                        for b in range(3):
                            t = np.zeros(nval)
                            for r in range(3):
                                id = (ii * 3 + a) * v + 1 + jj * 3 + r
                                id1 = (i * 3 + a) * v + 1 + j * 3 + b
                                if E is None:
                                    t[id] += rot[r, b]
                                    t[id1] -= rot[a, r]
                                else:
                                    t[id] += E[id] * rot[r, b]
                                    t[id1] -= E[id1] * rot[a, r]

                            # phi[ii,jj].dot(rot)=rot.dot(phi[i,j])
                            Q.append(t)
                            if (len(Q) == 50):
                                e, c = np.linalg.eig(Q)
                                if E is None:
                                    E = c
                                else:
                                    E = E.dot(c)
                                nval = E.shape[1]
                                print("nval:", nval)
                                Q = []
        self.R = E
        v = norm(u, axis=2)
        u0 = v.flatten().max()
        F, A = self.getMatrix(F, u / u0)

        return F, A.dot(self.R)

    def run(self):
        atoms = self.getsatoms()
        symmetry = self.getSupercell(atoms)
        files = tl.shell_exec(
            'find dirs/dir_* -name vasprun.xml |sort -n').split('\n')
        if len(files) > 100:
            files = files[:100]
        pos = atoms.positions
        f, u = self.getForce(pos, files)
        F, A = self.getCsMat(f, u, symmetry)
        print("start compressive sensing ")
        B = cs(mu=self.mu, split=self.split, lam=self.lam).run(F, A)
        print("rebuilding IFCs ")
        phi = self.rebuild(B)
        print("writing IFCs ")
        v = norm(u, axis=2)
        u0 = v.flatten().max()
        fc2 = np.einsum(phi[1], [1, 0, 3, 2]) / u0

        writefc2(fc2, 'csfc2')
        if self.NAH >= 3:
            a = h5py.File('fc3p.hdf5')
            if 'fc3' in a:
                del a['fc3']
            a['fc3'] = phi[2] / u0 / u0
            a.close()
            self.fc3()

    def fc3(self):
        print("writing csfc3 ")
        a = h5py.File('fc3p.hdf5')
        fc3 = np.einsum(a['fc3'], [0, 2, 1, 3, 5, 4])
        from ase import io
        atoms = io.read('POSCAR')
        satoms = self.getsatoms()
        writefc3(fc3, atoms, satoms, 'csfc3')

    def rebuild(self, B):
        n = self.natom
        rowr = self.rowr
        B = self.R.dot(B).T.reshape([-1, 3 * n]).T
        phi = []
        for i in range(self.NAH):
            r = range(rowr[i], rowr[i + 1])
            x = B[r].reshape([n, 3] * (i + 1))
            idx = np.array([0, i + 1])
            rdx = []
            for j in range(i):
                rdx.extend(idx + (j + 1))
            rdx.extend(idx)

            x = np.einsum(x, rdx)
            phi.append(x)

        return phi


class cssklearn:

    def __init__(self):
        pass

    def initu(self, f, A):
        dim = list(f.shape)
        dim[0] = A.shape[1]
        # so dim is the shape of u
        return np.ones(dim)

    def run(self, f, A):
        # from sklearn import cross_validation
        from sklearn import linear_model
        reg = linear_model.Lasso(
            alpha=1e-15, fit_intercept=False, max_iter=10000, tol=1e-5)
        print(reg.fit([[0, 0, 2], [1, 1, 2]], [[0, 1], [1, 1]]).coef_.T)
        print(A.shape, f.shape)
        return reg.fit(A, f).coef_.T
        # k_fold = cross_validation.KFold(n=len(f), n_folds=10)
        # [svc.fit(X_digits[train], y_digits[train])\
        # .score(X_digits[test], y_digits[test]) for train, test in kfold]


class csfortran:

    def __init__(self):
        pass

    def initu(self, f, A):
        dim = list(f.shape)
        dim[0] = A.shape[1]
        # so dim is the shape of u
        return np.ones(dim)

    def run(self, f, A):
        u = self.initu(f, A)
        import sys
        sys.path.append("/home/xggong/home1/zhouy/soft/bcs-master/wrap")
        import bcs as p
        ebars = np.zeros(len(A[0]))
        sigma2 = np.std(f) / 100.
        p.bcs.do_wrapped(A, f, sigma2, 1e-8, u, ebars)
        return u


class cs:

    def __init__(self, mu=0.7, split=True, lam=0.9):
        self.mu, self.lam = mu, lam
        self.split = split

    def initu(self, f, A):
        dim = list(f.shape)
        dim[0] = A.shape[1]
        # so dim is the shape of u
        return np.ones(dim)

    def testcs(self):
        f = (np.ones(1) * 20.0).reshape(1, 1)
        A = np.array([7.0, 10.0]).reshape(1, 2)
        print(self.run(f, A))

    def test2(self):
        f = np.array([7.0, 8.0])
        A = np.array([[1.0, 0], [1.0, 0]])
        print(self.run(f, A))

    def test3(self):
        f = np.array([7.0, 8.0])
        A = np.array([[1.0, 0, 0], [1.0, 0, 0]])
        print(self.run(f, A))

    def test4(self):
        f = np.array([7.0, 8.0]).reshape(1, 2)
        A = np.array([[1.0, 0]])
        print(self.run(f, A))

    def run(self, f, A):

        # normalize
        print("normalizing sensing matrix")

        # from scipy.sparse.linalg import eigsh
        """
            aA=eigsh(A.T.dot(A),k=6)[0].max()#the largest eigenvalue
            f/=np.sqrt(aA)
            # print norm(f)
            A/=np.sqrt(aA)
            # print norm(A)
        """
        aA = np.double(A.shape[0]**A.max().max())
        # maxeig(A.T.dot(A))
        f /= np.sqrt(aA)
        A /= np.sqrt(aA)
        """

            v=np.eye(len(A.T))-A.T.dot(A)
            for i in range(20):
                v=v.dot(v)
                print norm(v)
        """
        if self.split:
            return self.split_bregman(f, A)
        else:
            return self.bregman(f, A)

    def split_bregman(self, f, A):
        def cc(u1):
            print("CG error:", norm(u1 - self.bb.flatten()) / norm(self.bb))
            tt1 = g(u1, A, f, lam, d, mu, b)
            print("CG target:", (tt1 - self.tt) / self.tt)
            self.tt = tt1

            self.bb = u1

        def g(u, *args):
            A, f, lam, d, mu, b = args
            u = u.reshape(shape)
            return 1.0 / 2 * norm(np.dot(A, u) - f)**2 + \
                lam / 2.0 * norm(d - mu * u - b)**2

        def dg(u, *args):
            A, f, lam, d, mu, b = args
            u = u.reshape(shape)
            return (A.T.dot(A.dot(u) - f) - lam * mu *
                    (d - b - mu * u)).flatten()

        u = self.initu(f, A)
        shape = u.shape
        d = np.zeros_like(u)
        b = np.zeros_like(u)
        deta = 0.001
        erru = 1.0
        lam = self.lam
        t = 1.0
        tt = 1.0
        self.tt = tt
        mu = self.mu
        scale = 1.0 / np.amax(np.abs(f)) * 1000.0
        print("scale=" + str(scale))
        f0 = np.zeros_like(f)
        self.bb = np.zeros_like(u)
        # f*=scale
        print('dimmensions:', A.shape, u.shape)
        while erru > deta:
            # g=lambda u:1.0/2*norm(dot(A,u.reshape(shape))-f)**2\
            # +lam/2.0*norm(d-mu*u.reshape(shape)-b)**2
            f1 = (f * scale - dot(A, u)) + (f0 + dot(A, u)) / 2
            u1 = optimize.fmin_cg(
                g,
                u,
                args=(A, f1, lam, d, mu, b),
                disp=True,
                fprime=dg,
                callback=cc,
                gtol=deta * 10).reshape(shape)

            d1 = shrink(mu * u1 + b, 1.0 / lam)
            b1 = b + mu * u1 - d1
            erru = norm(u1 - u) / norm(u)
            print('split bregman iteration error:', erru)
            b = b1
            u = u1
            d = d1
            f0 = f1
            t1 = norm(d, 1) + tt
            print('change of target func:', (t1 - t) / t)
            t = t1
        return u / scale

    def bregman(self, f, A):
        u = self.initu(f, A)
        f0 = np.zeros_like(f)
        deta = 0.0001
        erru = 1
        scale = 1000.0
        while erru > deta:
            f1 = f * scale + f0 - dot(A, u)
            u1 = self.FCP(f1, A, u)

            erru = norm(u1 - u) / norm(u)
            print('bregman iteration error:', erru)
            u = u1
            f0 = f1
        return u / scale

    def FCP(self, f, A, u=None):
        if u is None:
            u = self.initu(f, A)

        m, n = A.shape
        if len(f.shape) > 1:
            n *= list(f.shape)[1]
        ta = 1.99999  # min(1.999,max(1.1,-1.665*np.float(m)/n+2.665))

        mu = self.mu
        deta = 0.01
        # errg=1
        erru = 1
        while erru > deta:  # or errg>deta:
            p = np.dot(A, u) - f

            g = np.dot(A.T, p)

            u1 = shrink(u - ta * g, mu * ta)
            # errg=1.0/mu*norm(g,np.inf)-1
            erru = norm(u1 - u) / norm(u)
            print('FCP iteration :', erru)
            u = u1

        return u

# encoding: utf-8
import unittest
from aces.modify import *
import numpy as np
from ase import Atoms


class testModify(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def testget_cluster(self):
        n = 4
        dis = np.ones([n, n])
        dis[0, 1] = 0
        dis[2, 3] = 0
        cluster = get_cluster(n, 0, dis, [])
        self.assertEqual(sorted(cluster), [1])

        cluster = get_cluster(n, 3, dis, [])
        self.assertEqual(sorted(cluster), [])

        n = 5
        dis = np.ones([n, n])
        dis[0, 1] = 0
        dis[1, 3] = 0
        dis[2, 4] = 0
        cluster = get_cluster(n, 0, dis, [])
        self.assertEqual(sorted(cluster), [1, 3])

        cluster = get_cluster(n, 2, dis, [])
        self.assertEqual(sorted(cluster), [4])

    def testget_clusters(self):
        n = 4
        dis = np.ones([n, n])
        dis[0, 1] = 0
        dis[2, 3] = 0
        clusters = get_clusters(n, dis)
        self.assertEqual(clusters, [True, False, True, False])

        n = 5
        dis = np.ones([n, n])
        dis[0, 1] = 0
        dis[1, 3] = 0
        dis[2, 4] = 0
        clusters = get_clusters(n, dis)
        self.assertEqual(clusters, [True, False, True, False, False])

    def testget_unique_atoms(self):
        unit = Atoms(
            "C4", positions=[
                (0, 0, 0), (0, 1, 0), (1, 1, 0), (1, 0, 0)])
        unit1 = unit.copy()
        unit1.translate([1, 0, 0])
        atoms = Atoms()
        atoms.extend(unit)
        atoms.extend(unit1)
        self.assertEqual(len(atoms), 8)
        atoms = get_unique_atoms(atoms)
        self.assertEqual(len(atoms), 6)


def test():
    unittest.main()
if __name__ == '__main__':
    test()

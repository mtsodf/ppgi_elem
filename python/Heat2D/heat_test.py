from heat import *
import unittest


class TestHeat(unittest.TestCase):

    def test_Dirichlet(self):
        diff = run_case(0, 40, 40)

        self.assertLess(diff, 1e-3)

    def test_Dirichlet_non_hom(self):

        diff = run_case(2, 40, 40)

        self.assertLess(diff, 1e-3)


    def test_Neumann(self):

        diff = run_case(1, 60, 60)

        self.assertLess(diff, 1e-3)


    def test_Neumann_non_hom(self):

        diff = run_case(3, 60, 60)

        self.assertLess(diff, 1e-3)


        diff = run_case(4, 90, 90)

        self.assertLess(diff, 1e-3)

if __name__ == '__main__':
    unittest.main()
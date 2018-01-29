from Strain import *
import unittest


class TestHeat(unittest.TestCase):

    def test_DeslX(self):

        print "Test Dirichlet Homogeneo"

        residue = run_case(0, 20, 20, 0.0, False, False)

        self.assertLess(residue, 3e-3)


    def test_DeslY(self):

        print "Test Dirichlet Homogeneo"

        residue = run_case(1, 20, 20, 0.0, False, False)

        self.assertLess(residue, 3e-3)

if __name__ == '__main__':
    unittest.main()

# generated with GPT-4 and hand-fixed
import unittest
import itertools
from src.search_n2 import two_compatible_2tuples, n_compatible_2tuples, wrapper


class TestSearchN2Functions(unittest.TestCase):
    def setUp(self):
        self.dist1 = (12, 18)
        self.dist2 = (10, 15)
        self.nsize = 20

        self.A = [(i, i+6) for i in range(5, 15+self.nsize)]
        self.B = [(i, i+6) for i in range(11, 30+self.nsize)]
        self.C = [(i, i+7) for i in range(25, 41+self.nsize)]
        self.test_lists = [self.A, self.B, self.C]
        self.test_dists = [self.dist1, self.dist2]

    def test_two_compatible_2tuples(self):
        I = [(3, 9), (6, 12), (10, 16)]
        J = [(18, 24), (20, 26), (22, 28)]
        dist = [8, 9]

        expected = [((6, 12), (20, 26)),]
        result = two_compatible_2tuples(I, J, dist)
        self.assertEqual(sorted(result), sorted(expected))

    def test_fixed_dist(self):
        I = [(3,9),]
        J = [(9, 15)]
        dist = [0,1]
        result = two_compatible_2tuples(I, J, dist)
        self.assertEqual(result, [((3,9), (9,15))])


    def test_n_compatible_2tuples(self):
        # expected result from manual calculation or reference implementation
        ref = [
            (a, b, c)
            for a, b, c in itertools.product(self.A, self.B, self.C)
            if b[0] - a[1] in range(*self.dist1) and c[0] - b[1] in range(*self.dist2)
        ]
        result = n_compatible_2tuples(self.test_lists, self.test_dists)
        self.assertEqual(sorted(result), sorted(ref))

    def test_wrapper(self):
        # expected result from manual calculation or reference implementation
        ref = [
            (a, b, c)
            for a, b, c in itertools.product(self.A, self.B, self.C)
            if b[0] - a[1] in range(*self.dist1) and c[0] - b[1] in range(*self.dist2)
        ]
        result = wrapper(self.test_lists, self.test_dists)
        self.assertEqual(sorted(result), sorted(ref))

    def test_wrapper_invalid_width(self):
        # Test for ValueError when widths are not the same
        invalid_list = [(5, 11), (6, 12), (7, 14)]
        with self.assertRaises(ValueError):
            wrapper([invalid_list, self.B, self.C], self.test_dists)

    def test_wrapper_invalid_b_greater_than_a(self):
        # Test for ValueError when b is not greater than a
        invalid_list = [(5, 4), (6, 5), (7, 6)]
        with self.assertRaises(ValueError):
            wrapper([invalid_list, self.B, self.C], self.test_dists)


if __name__ == '__main__':
    unittest.main()

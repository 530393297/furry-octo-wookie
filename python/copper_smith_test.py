import unittest
import copper_smith as cs
import rsa
import random
import decimal as d
import math


class TestCopperSmithMethods(unittest.TestCase):

    def test_inner(self):
        x = [1, 2, 3, 4, 5]
        y = [7, 8, 9, 10, 11]
        ans = sum([7, 16, 27, 40, 55])
        self.assertEqual(cs.inner(x, y), ans)

        x = [0.8, -5, 3.2, 6.1, 7]
        y = [11, 3, -14, 0.22, 11]
        ans = sum([8.8, -15, -44.8, 1.342, 77])
        self.assertAlmostEqual(cs.inner(x, y), ans)

    def test_sub(self):
        x = [1, 2, 3, 4, 5]
        y = [7, 8, 9, 10, 11]
        ans = [-6, -6, -6, -6, -6]
        self.assertEqual(cs.sub(x, y), ans)

        x = [0.8, -5, 3.2, 6.1, 7]
        y = [11, 3, -14, 0.22, 11]
        ans = [-10.2, -8, 17.2, 5.88, -4]
        out = cs.sub(x, y)
        for i in range(5):
            self.assertAlmostEqual(ans[i], out[i])

    def test_add(self):
        x = [1, 2, 3, 4, 5]
        y = [7, 8, 9, 10, 11]
        ans = [8, 10, 12, 14, 16]
        self.assertEqual(cs.add(x, y), ans)

        x = [0.8, -5, 3.2, 6.1, 7]
        y = [11, 3, -14, 0.22, 11]
        ans = [11.8, -2, -10.8, 6.32, 18]
        out = cs.add(x, y)
        for i in range(5):
            self.assertAlmostEqual(ans[i], out[i])

    def test_mult(self):
        x = [1, 2, 3, 4, 5]
        y = 2
        ans = [2, 4, 6, 8, 10]
        self.assertEqual(cs.mult(y, x), ans)

        x = [0.8, -5, 3.2, 6.1, 7]
        y = -14
        ans = [-11.2, 70, -44.8, -85.4, -98]
        out = cs.mult(y, x)
        for i in range(5):
            self.assertAlmostEqual(ans[i], out[i])

    def test_lll_one(self):
        data = [[1, 1, 1], [-1, 0, 2], [3, 5, 6]]
        ans = [[0, 1, 0], [1, 0, 1], [-1, 0, 2]]
        self.assertEqual(cs.lll(data, 0.75), ans)

    def test_lll_two(self):
        data = [[201, 37], [1648, 297]]
        ans = [[1, 32], [40, 1]]
        self.assertEqual(cs.lll(data, 0.75), ans)

    def test_lll_three(self):
        data = [[15, 23, 11], [46, 15, 3], [32, 1, 1]]
        ans = [[1, 9, 9], [13, 5, -7], [6, -9, 15]]
        self.assertEqual(cs.lll(data, 0.75), data)

    def test_poly_deriv(self):
        x = [13, 8, 1]
        ans = [8, 2]
        self.assertEqual(ans, cs.poly_derivative(x))

    def test_generate(self):
        h = 3
        N = 35
        f = [19, 14, 1]
        k = 2
        X = math.ceil(
            (pow(2, - 1 / 2) * pow(h * k, -1 / (h * k - 1))) * pow(N, (h - 1) / (h * k - 1))) - 1
        # RUN
        degree = 2
        ans = [[1225.0, 0, 0, 0, 0, 0],
               [0.0, 2450.0, 0, 0, 0, 0],
               [665.0, 980.0, 140.0, 0, 0, 0],
               [0.0, 1330.0, 1960.0, 280.0, 0, 0],
               [361.0, 1064.0, 936.0, 224.0, 16.0, 0],
               [0.0, 722.0, 2128.0, 1872.0, 448.0, 32.0]]
        self.assertEqual(cs.generate(f, N, h, X, degree), ans)

    def test_generate_and_lll(self):
        h = 3
        N = 35
        f = [19, 14, 1]
        degree = 2
        k = 2
        X = math.ceil(
            (pow(2, - 1 / 2) * pow(h * k, -1 / (h * k - 1))) * pow(N, (h - 1) / (h * k - 1))) - 1
        gen = cs.generate(f, N, h, degree, X)

        ans = [[3, 8 * 2, -24 * 2 * 2, -8 * pow(2, 3), -1 * pow(2, 4), 2 * pow(2, 5)],
               [49, 50 * 2, 0 * 2 * 2, 20 *
                   pow(2, 3), 0 * pow(2, 4), 2 * pow(2, 5)],
               [115, -83 * 2, 4 * 2 * 2, 13 *
                   pow(2, 3), 6 * pow(2, 4), 2 * pow(2, 5)],
               [61, 16 * 2, 37 * 2 * 2, -16 *
                   pow(2, 3), 3 * pow(2, 4), 4 * pow(2, 5)],
               [21, -37 * 2, -14 * 2 * 2, 2 *
                   pow(2, 3), 14 * pow(2, 4), -4 * pow(2, 5)],
               [-201, 4 * 2, 33 * 2 * 2, -4 * pow(2, 3), -3 * pow(2, 4), 1 * pow(2, 5)]]

        self.assertEqual(cs.lll(gen, d.Decimal(0.75)), ans)

    def test_solve(self):
        h = 3
        N = 35
        f = [19, 14, 1]
        k = 2
        print(cs.run_copper_smith(f, N, h, k))

    #def test_coron(self):
    #    N = 2122840968903324034467344329510307845524745715398875789936591447337206598081
    #    h = 3
    #    c = 1792963459690600192400355988468130271248171381827462749870651408943993480816
    #    f = [pow(2, 1500) - c, 3 * pow(2, 1000), 3 * pow(2, 1000), 1]
    #    k = 3
    #    print(cs.run_copper_smith(f, N, h, k))


if __name__ == '__main__':
    unittest.main()
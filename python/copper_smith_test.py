import unittest
import copper_smith as cs
import rsa
import random


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

    def test_two_degree_coppersmith_small(self):
        p = 106859
        q = 106661
        Nn = p * q
        e = 3
        n = Nn.bit_length()

        x0 = random.randint(0, pow(2, n // 3 - 10))
        a = random.randint(0, Nn)
        b = -x0 * x0 - a * x0 % Nn
        A = [1, a, b]

        self.assertAlmostEqual(x0, cs.run_copper_smith(A, Nn), x0)

    def test_three_degree_coppersmith_small(self):
        p = rsa.random_prime(1000000)
        q = rsa.random_prime(1000000)
        Nn = p * q
        e = 3
        n = Nn.bit_length()

        x0 = random.randint(0, round(pow(Nn, 1 / 4) / 6))
        a = random.randint(0, Nn)
        b = random.randint(0, Nn)
        c = (-x0 * x0 * x0 - a * x0 * x0 - b * x0) % Nn
        A = [1, a, b, c]

        self.assertAlmostEqual(x0, cs.run_copper_smith(A, Nn), x0)


if __name__ == '__main__':
    unittest.main()

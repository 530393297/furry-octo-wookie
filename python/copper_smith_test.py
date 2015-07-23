import unittest
import copper_smith as cs


class TestStringMethods(unittest.TestCase):

    def test_lll_wikipedia(self):
        data = [[1, 1, 1], [-1, 0, 2], [3, 5, 6]]
        ans = [[0, 1, 0], [1, 0, 1], [-1, 0, 2]]
        out = cs.lll(data, 0.75)
        self.assertEqual(len(out), len(ans))
        for _ in range(3):
            for _ in range(3):
                self.assertEqual(ans, data)

    def test_lll_youtube(self):
        data = [[201, 37], [1648, 297]]
        ans = [[1, 32], [40, 1]]
        out = cs.lll(data, 0.75)
        self.assertEqual(len(out), len(ans))
        for _ in range(2):
            for _ in range(2):
                self.assertEqual(ans, data)

    def test_lll_youtube_two(self):
        data = [[201, 37], [1648, 297]]
        ans = [[1, 32], [40, 1]]
        out = cs.lll(data, 0.75)
        self.assertEqual(len(out), len(ans))
        for _ in range(2):
            for _ in range(2):
                self.assertEqual(ans, data)


if __name__ == '__main__':
    unittest.main()

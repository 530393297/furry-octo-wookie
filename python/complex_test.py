import unittest
import complex as c
import math

class ComplexTest(unittest.TestCase):

	def test_add(self):
		a = c.Complex(5, 6)
		b = c.Complex(-3, 4)
		d = a + b
		self.assertEqual(d.i, 10)
		self.assertEqual(d.r, 2)

	def test_sub(self):
		a = c.Complex(5, 6)
		b = c.Complex(-3, 4)
		d = a - b
		self.assertEqual(d.i, 2)
		self.assertEqual(d.r, 8)

	def test_mult(self):
		a = c.Complex(5, 6)
		b = c.Complex(-3, 4)
		d = a * b
		self.assertEqual(d.i, 2)
		self.assertEqual(d.r, -39)

	def test_conjugate(self):
		a = c.Complex(5, 6)
		a = a.conjugate()
		self.assertEqual(a.i, -6)
		self.assertEqual(a.r, 5)

	def test_print(self):
		a = c.Complex(5, 6)
		ans = str(a)
		self.assertEqual(ans, "5 + 6i")
		a = c.Complex(5, 0)
		ans = str(a)
		self.assertEqual(ans, "5")
		a = c.Complex(0, 6)
		ans = str(a)
		self.assertEqual(ans, "6i")
		a = c.Complex(5, -6)
		ans = str(a)
		self.assertEqual(ans, "5 - 6i")

	def test_sin(self):
		a = c.Complex(5, 6)
		b = a.sin()
		self.assertAlmostEqual(b.r, -193.4300200, 5)
		self.assertAlmostEqual(b.i, 57.21839505, 5)

	def test_cos(self):
		a = c.Complex(5, 6)
		b = a.cos()
		self.assertAlmostEqual(b.r, 57.219098, 5)
		self.assertAlmostEqual(b.i, 193.427643, 5)

	def test_sqrt(self):
		a = c.Complex(5, 6)
		b = a.sqrt()
		self.assertAlmostEqual(b.r, 2.53083, 5)
		self.assertAlmostEqual(b.i, 1.18538, 5)

		a = c.Complex(5, 0)
		b = a.sqrt()
		self.assertAlmostEqual(b.r, 2.2360679, 5)
		self.assertAlmostEqual(b.i, 0, 5)

		a = c.Complex(-5, 0)
		b = a.sqrt()
		self.assertAlmostEqual(b.i, 2.2360679, 5)
		self.assertAlmostEqual(b.r, 0, 5)

	def test_div(self):
		a = c.Complex(5, 6)
		b = c.Complex(-3, 4)
		d = a / b
		self.assertEqual(d.i, -1.52)
		self.assertEqual(d.r, 0.36)




if __name__ == '__main__':
    unittest.main()
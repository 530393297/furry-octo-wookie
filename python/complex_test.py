import unittest
import decimal as d
import complex as c
import math


class ComplexDecimalTest(unittest.TestCase):

    def decimal_almost_equals(self, lhs, rhs, epsilon=7):
        epsilon = pow(10, -epsilon)
        if abs(d.Decimal(lhs - rhs)) < epsilon:
            return True
        return False

    def check_type(self, result):
        if isinstance(result.real, d.Decimal) and isinstance(result.imag, d.Decimal):
            return True
        return False

    def check_answer(self, comp, real, imag, epsilon=7):
        if self.decimal_almost_equals(comp.real, real, epsilon) is False:
            return False
        if self.decimal_almost_equals(comp.imag, imag, epsilon) is False:
            return False
        return True

    def test_add(self):
        a = c.ComplexDecimal(5, 6)
        b = c.ComplexDecimal(-3, 4)
        result = a + b
        self.assertTrue(self.check_answer(result, d.Decimal(2), d.Decimal(10)))
        self.assertTrue(self.check_type(result))

    def test_sub(self):
        a = c.ComplexDecimal(5, 6)
        b = c.ComplexDecimal(-3, 4)
        result = a - b
        self.assertTrue(self.check_answer(result, d.Decimal(8), d.Decimal(2)))
        self.assertTrue(self.check_type(result))

    def test_mult(self):
        a = c.ComplexDecimal(5, 6)
        b = c.ComplexDecimal(-3, 4)
        result = a * b
        self.assertTrue(
            self.check_answer(result, d.Decimal(-39), d.Decimal(2)))
        self.assertTrue(self.check_type(result))

    def test_div(self):
        a = c.ComplexDecimal(5, 6)
        b = c.ComplexDecimal(-3, 4)
        result = a / b
        self.assertTrue(
            self.check_answer(result, d.Decimal(0.36), d.Decimal(-1.52)))
        self.assertTrue(self.check_type(result))

    def test_conjugate(self):
        a = c.ComplexDecimal(5, 6)
        result = a.conjugate()
        self.assertTrue(self.check_answer(result, d.Decimal(5), d.Decimal(-6)))
        self.assertTrue(self.check_type(result))

    def test_sqrt(self):
        a = c.ComplexDecimal(5, 6)
        result = a.sqrt()
        self.assertTrue(
            self.check_answer(result, d.Decimal(2.5083), d.Decimal(1.18538), 1))

        a = c.ComplexDecimal(5, 0)
        result = a.sqrt()
        self.assertTrue(
            self.check_answer(result, d.Decimal(2.23667), d.Decimal(0), 1))

        a = c.ComplexDecimal(-5, 0)
        result = a.sqrt()
        self.assertTrue(
            self.check_answer(result, d.Decimal(0), d.Decimal(2.23667), 1))
        self.assertTrue(self.check_type(result))

    def test_print(self):
        a = c.ComplexDecimal(5, 6)
        ans = str(a)
        self.assertEqual(ans, "5 + 6i")
        a = c.ComplexDecimal(5, 0)
        ans = str(a)
        self.assertEqual(ans, "5")
        a = c.ComplexDecimal(0, 6)
        ans = str(a)
        self.assertEqual(ans, "6i")
        a = c.ComplexDecimal(5, -6)
        ans = str(a)
        self.assertEqual(ans, "5 - 6i")


if __name__ == '__main__':
    unittest.main()

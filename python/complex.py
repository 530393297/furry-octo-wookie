import math
import decimal as d

class Complex:
	def __init__(self, realpart, imagpart):
		self.r = realpart
		self.i = imagpart

	def __add__(self, x):
		real = self.r + x.r
		imag = self.i + x.i
		return Complex(real, imag)

	def __sub__(self, x):
		real = self.r - x.r
		imag = self.i - x.i
		return Complex(real, imag)

	def __mul__(self, x):
		old = self.r
		real = self.r * x.r - self.i * x.i
		imag = old * x.i + self.i * x.r
		return Complex(real, imag)

	def conjugate(self):
		return Complex(self.r, -self.i)

	def __str__(self):
		if self.i == 0:
			return str(self.r)
		if self.r == 0:
			return str(self.i) + "i"
		if self.i < 0:
			return str(self.r) + " - " + str(-self.i) + "i"
		return str(self.r) + " + " + str(self.i) + "i"

	def reciprocal(self):
		scale = self.r * self.r + self.i * self.i
		return Complex(self.r / scale, -self.i / scale)

	def __truediv__(self, x):
		return self * x.reciprocal()


	def sin(self):
		return Complex(math.sin(self.r) * math.cosh(self.i), math.cos(self.r) * math.sinh(self.i))

	def cos(self):
		return Complex(math.cos(self.r) * math.cosh(self.i), -math.sin(self.r) * math.sinh(self.i))

	def sqrt(self):
		if self.i == 0:
			if self.r > 0:
				return Complex(math.sqrt(self.r), 0)
			else:
				return Complex(0, math.sqrt(-self.r))
		
		if self.r == 0:
			print("I SHOULDNT BE HERE")
			return Complex(0, -math.sqrt(self.i))
		
		real = (1 / d.Decimal(2).sqrt()) * d.Decimal(d.Decimal(self.r * self.r + self.i * self.i).sqrt() + self.r).sqrt()
		signb = 1
		if self.i < 0:
			signb = d.Decimal(-1)
		imag =  (signb / d.Decimal(2).sqrt()) * d.Decimal(d.Decimal(self.r * self.r + self.i * self.i).sqrt() - self.r).sqrt()
		return Complex(real, imag)





if __name__ == "__main__":
    import doctest
    doctest.testmod()
import decimal as d


class ComplexDecimal(object):

    def __init__(self, real, imag):
        self.real = d.Decimal(real)
        self.imag = d.Decimal(imag)

    def __add__(self, rhs):
        return ComplexDecimal(self.real + rhs.real, self.imag + rhs.imag)

    def __sub__(self, rhs):
        return ComplexDecimal(self.real - rhs.real, self.imag - rhs.imag)

    def __mul__(self, rhs):
        return ComplexDecimal(self.real * rhs.real - self.imag * rhs.imag,
                              self.real * rhs.imag + self.imag * rhs.real)

    def __truediv__(self, rhs):
        return self * rhs.reciprocal()

    def __neg__ (self):
    	return ComplexDecimal(-self.real, -self.imag)

    def __str__(self):
        if self.imag == 0:
            return str(self.real)
        if self.real == 0:
            return str(self.imag) + "i"
        if self.imag < 0:
            return str(self.real) + " - " + str(-self.imag) + "i"
        return str(self.real) + " + " + str(self.imag) + "i"

    def conjugate(self):
        return ComplexDecimal(self.real, -self.imag)

    def reciprocal(self):
        scale = self.real * self.real + self.imag * self.imag
        return ComplexDecimal(self.real / scale, -self.imag / scale)

    def size(self):
    	return ComplexDecimal(self.real*self.real+self.imag*self.imag, 0)
    def sqrt(self):
        if self.imag == 0:
            if self.real > 0:
                return ComplexDecimal(self.real.sqrt(), 0)
            else:
                return ComplexDecimal(0, d.Decimal(-self.real).sqrt())

        if self.real == 0:
          	if self.imag > 0:
          		return ComplexDecimal((d.Decimal(2).sqrt() / 2) * self.imag.sqrt(), (d.Decimal(2).sqrt() / 2) * self.imag.sqrt())
          	else:
          		return ComplexDecimal((d.Decimal(2).sqrt() / 2) * d.Decimal(-self.imag).sqrt(), -(d.Decimal(2).sqrt() / 2) *  d.Decimal(-self.imag).sqrt())

        real = (1 / d.Decimal(2).sqrt()) * d.Decimal(
            d.Decimal(self.real * self.real + self.imag * self.imag).sqrt() + self.real).sqrt()
        signb = 1
        if self.imag < 0:
            signb = d.Decimal(-1)
        imag = (signb / d.Decimal(2).sqrt()) * d.Decimal(
            d.Decimal(self.real * self.real + self.imag * self.imag).sqrt() - self.real).sqrt()
        return ComplexDecimal(real, imag)


if __name__ == "__main__":
    import doctest
    doctest.testmod()

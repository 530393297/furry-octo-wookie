import decimal as d
import scipy.optimize as opt
import numpy.polynomial.polynomial as poly
import functools
import operator
import random
import rsa
import math
import complex as c


def inner(vector0, vector1):
    """ Calculate the inner product of two vectors """
    return sum(map(operator.mul, vector0, vector1))


def sub(vector0, vector1):
    """ Subtract each in two vectors, returning a new vector """
    return list(map(operator.sub, vector0, vector1))


def add(vector0, vector1):
    """ Add each value in two vectors together, returning a new vector """
    return list(map(operator.add, vector0, vector1))


def mult(coefficient, vector):
    """ Multiply each element in a vector against a coefficient """
    return list(i * coefficient for i in vector)

def div(coefficient, vector):
    """ divide each element in a vector against a coefficient """
    return list(i / coefficient for i in vector)


def gram_schmidt(basis):
    """ Calculate the orthonormal set from a given basis set of vectors

    Arguments:
    basis -- the basis vectors
    """
    ortho = basis.copy()
    for i in range(len(basis)):
        for j in range(i):
            ortho[i] = sub(
                ortho[i], mult(inner(basis[i], ortho[j]) / inner(ortho[j], ortho[j]), ortho[j]))
    return ortho


def lll(basis, delta=0.75):
    """ Run the Lenstra–Lenstra–Lovasz (LLL) lattice basis reduction algorithm

    Arguments:
    basis -- the basis vectors
    delta -- the delta used in the Lovasz condition (default 0.75)
    """
    ortho = gram_schmidt(basis)

    print("lll")

    for i in range(1, len(basis)):
        for j in reversed(range(i)):
            coefficient = round(
                inner(basis[i], ortho[j]) / inner(ortho[j], ortho[j]))
            basis[i] = sub(basis[i], mult(coefficient, basis[j]))
    for i in range(len(basis) - 1):
        product = inner(ortho[i], ortho[i])
        lovasz_inner = add(
            mult(inner(basis[i + 1], ortho[i]) / product, ortho[i]), ortho[i + 1])
        print(lovasz_inner)
        if delta * product > inner(lovasz_inner, lovasz_inner):
            basis[i], basis[i + 1] = basis[i + 1], basis[i]
            return lll(basis, delta)
    return basis

def poly_pad(f, n):
	g = [d.Decimal(0.0)] * n
	for i in range(len(f)):
		g[i] = f[i]
	return g

def poly_mul(f1, f2):
    res = [d.Decimal(0.0)] * (len(f1)+len(f2)-1)
    for o1, i1 in enumerate(f1):
        for o2, i2 in enumerate(f2):
            res[o1 + o2] += i1 * i2
    return res

def poly_pow(x, n):
    if n == 0:
        return [1]
    y = [1]
    while n > 1:
      if n % 2 == 0:
        x = poly_mul(x, x)
        n = n >> 1 
      else:
        y = poly_mul(x, y)
        x = poly_mul(x, x)
        n = (n - 1) / 2
    return poly_mul(x, y)

def poly_derivative(f):
    g = [d.Decimal(0.0)]*(len(f)-1)
    for i in range(1, len(f)):
        g[i - 1] = i * f[i]
    return g

def newton(f, x0):
    f_prime = poly_derivative(f)
    tolerance = pow(10, -7)
    eps = pow(10, -14)

    iteration = 500

    for i in range(iteration):

        y = func(f, x0)
        y_prime = func(f_prime, x0)

        if abs(y_prime) < eps:
            break

        x1 = x0 - y / y_prime

        if abs(x1 - x0) / abs(x1) < tolerance:
            return x1

        x0 = x1
    return -1

def normalize(poly):
    while poly and poly[-1] == 0:
        poly.pop()
    if poly == []:
        poly.append(0)


def poly_divmod(num, den):
    #Create normalized copies of the args
    num = num[:]
    normalize(num)
    den = den[:]
    normalize(den)

    if len(num) >= len(den):
        #Shift den towards right so it's the same degree as num
        shiftlen = len(num) - len(den)
        den = [0] * shiftlen + den
    else:
        return [0], num

    quot = []
    divisor = c.ComplexDecimal(den[-1], 0)
    for i in range(shiftlen + 1):
        #Get the next coefficient of the quotient.
        mult = num[-1] / divisor
        quot = [mult] + quot

        #Subtract mult * den from num, but don't bother if mult == 0
        #Note that when i==0, mult!=0; so quot is automatically normalized.
        if mult != 0:
            dd = [mult * u for u in den]
            num = [u - v for u, v in zip(num, dd)]

        num.pop()
        den.pop(0)

    normalize(num)
    return quot, num

def generate(p, N, h, k, X):
	poly_matrix = [[d.Decimal(0.0) for _ in range(h * k)] for _ in range(h * k)]
	n = [N]
	x = [0, 1]

	for i in range(1, (h * k) + 1):
		v = math.floor((i - 1) / k)
		v = max(v, 0)
		u = (i - 1) - (k * v)
		q = poly_mul(poly_mul(poly_pow(n, h - 1 - v), poly_pow(x, u)), poly_pow(p, v))
		q = poly_pad(q, h * k)
		for j in range(h * k):
			poly_matrix[i - 1][j] = d.Decimal(q[j] * pow(X, j)) 

	return poly_matrix

def power_iteration(basis):
    
    b = [c.ComplexDecimal(random.random(), random.random()) for _ in range(len(basis))]
    norm_sq = c.ComplexDecimal(0, 0)

    for k in range(len(b)):
        norm_sq += b[k].size()
    norm = norm_sq.sqrt()

    b = div(norm, b)

    out = [c.ComplexDecimal(0, 0)] * len(basis)
    previous = c.ComplexDecimal(0, 0)
    new = c.ComplexDecimal(1, 0)

    for _ in range(100000):
        for i in range(len(basis)):
            out[i] = c.ComplexDecimal(0, 0)
            for j in range(len(basis)):
                out[i] = out[i] + basis[i][j].conjugate() * b[j]

        norm_sq = c.ComplexDecimal(0, 0)
        for k in range(len(basis)):
            norm_sq += out[k] * out[k].conjugate()
        norm = norm_sq.sqrt()

        previous = new
        new = out[0] / b[0]
        #print(1 / new)
        b = div(norm, out)

    return new, out


def func(f, x):
    total = 0
    y = d.Decimal(x)
    for i in range(len(f)):
        total = total + pow(y, i) * f[i]
    return total

def solve(p, X, N):
    p = [p[i] / pow(X,i) for i in range(len(p))]
    
    # This is the only line that needs fixing
    x_init = d.Decimal(pow(N,1/3)/pow(2,10))
    return newton(p, x_init)

def eigen_matrix(f):
    out = [[c.ComplexDecimal(d.Decimal(0), d.Decimal(0)) for _ in range(len(f) - 1)] for _ in range(len(f) - 1)]
    for i in range(len(f) - 2):
        out[i + 1][i] = c.ComplexDecimal(1, 0)
    for i in range(len(f) - 1):
        out[0][i] = -f[i + 1] / f[0]
    return out


def find_n_roots(f, n):
    mat = eigen_matrix(f)
    ans = [0 for _ in range (n)]

    ans[0] = power_iteration(mat, d.Decimal(0))

    for i in range(1, n):
        ans[i] = power_iteration(mat, d.Decimal(ans[i - 1]))
        ans[i] = ans[i] + ans[i - 1]

    return [1 / i for i in ans]

def mu_transpose(A,u, mu):
    v = [[c.ComplexDecimal(0,0) for _ in range(len(u))]for _ in range(len(u))]
    for i in range(len(u)):
        for j in range(len(u)):
            v[i][j] = A[i][j]-mu*u[i]*u[j].conjugate()
    return v

def run_copper_smith(f, N, h, k):
    X = math.ceil((pow(2, - 1 / 2) * pow(h * k, -1 / (h * k - 1))) * pow(N, (h - 1) / (h * k - 1))) - 1
    gen = lll(generate(f, N, h, k, X),  d.Decimal(0.75))
    final_poly = [c.ComplexDecimal(gen[0][i] / pow(X, i), 0) for i in range(len(gen[0]))]
    

    out = [0] * (len(final_poly) - 1)
    a = eigen_matrix(final_poly)
    lamb, u = power_iteration(a)
    print(c.ComplexDecimal(1,0)/lamb)
    a = mu_transpose(a,u,c.ComplexDecimal(1,0)/lamb)
    lamb, u = power_iteration(a)
    print(c.ComplexDecimal(1,0)/lamb)
    a = mu_transpose(a,u,c.ComplexDecimal(1,0)/lamb)
    lamb, u = power_iteration(a)
    print(c.ComplexDecimal(1,0)/lamb)
    a = mu_transpose(a,u,c.ComplexDecimal(1,0)/lamb)
    lamb, u = power_iteration(a)
    print(c.ComplexDecimal(1,0)/lamb)
    a = mu_transpose(a,u,c.ComplexDecimal(1,0)/lamb)
    lamb, u = power_iteration(a)
    print(c.ComplexDecimal(1,0)/lamb)

#    for i in range(len(final_poly) - 1):
#        out[i] = c.ComplexDecimal(1, 0)/ power_iteration(eigen_matrix(final_poly)) 
#        print(out[i])
#        final_poly = poly_divmod(final_poly, [-out[i].real, 1])[0]


    return out

#f = [6, 5, 1]

#print(find_n_roots(f, 2))

inp = [[1,1,1], [-1,0,2], [3,5,6]]
print(lll(inp))
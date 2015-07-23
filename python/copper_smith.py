import decimal as d
import scipy.optimize as opt
import functools
import operator


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


def gram_schmidt(basis):
    """ Calculate the orthonormal set from a given basis set of vectors

    Keyword arguments:
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

    Keyword arguments:
    basis -- the basis vectors
    delta -- the delta used in the Lovasz condition (default 0.75)
    """
    ortho = gram_schmidt(basis)

    for i in range(1, len(basis)):
        for j in reversed(range(i)):
            coefficient = round(
                inner(basis[i], ortho[j]) / inner(ortho[j], ortho[j]))
            basis[i] = sub(basis[i], mult(coefficient, basis[j]))
    for i in range(len(basis) - 1):
        product = inner(ortho[i], ortho[i])
        lovasz = add(
            mult(inner(basis[i + 1], ortho[i]) / product, ortho[i]), ortho[i + 1])
        if delta * product > inner(lovasz, lovasz):
            basis[i], basis[i + 1] = basis[i + 1], basis[i]
            return lll(basis, delta)
    return basis


def copper_smith(poly, N):
    """ Perform a Copper-Smith attack

    Keyword arguments:
    poly -- the basis polynomial
    N -- p * q used in the encryption
    """
    degree = len(poly)
    X = d.Decimal(pow(N, 1 / degree) / ((degree - 1) * 2))
    poly_matrix = [[d.Decimal(0.0) for _ in range(degree)] for _ in range(degree)]

    for i in range(degree):
        poly_matrix[i][i] = pow(X, degree - i) * N
        poly_matrix[0][i] = pow(X, degree - i) * poly[i]

    reduced = lll(poly_matrix, d.Decimal(0.75))
    return [reduced[0][i] / pow(X, (degree - 1 - i)) for i in range(degree)]


def poly_func(f, x):
    total = 0
    for i in range(len(f)):
        total = total + pow(x, (len(f) - 1) - i) * round(f[i])
    return total


def run_copper_smith(a, N):
    n = N.bit_length()

    p = copper_smith(a, N)
    f = functools.partial(poly_func, p)

    # This is hard coded for degree two
    return opt.brentq(f, 0, pow(2, n // 3 - 10))

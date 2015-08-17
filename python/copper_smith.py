import operator
import sympy as s


def inner(vector0, vector1):
    return sum(map(operator.mul, vector0, vector1))


def sub(vector0, vector1):
    return list(map(operator.sub, vector0, vector1))


def add(vector0, vector1):
    return list(map(operator.add, vector0, vector1))


def mult(coefficient, vector):
    return list(i * coefficient for i in vector)


def div(coefficient, vector):
    return list(i / coefficient for i in vector)


def gram_schmidt(basis):
    ortho = basis.copy()
    for i in range(len(basis)):
        for j in range(i):
            ortho[i] = sub(
                ortho[i], mult(inner(basis[i], ortho[j]) /
                               inner(ortho[j], ortho[j]), ortho[j]))
    return ortho


def lll(basis, delta=s.S(0.75)):
    ortho = gram_schmidt(basis)

    for i in range(1, len(basis)):
        for j in reversed(range(i)):
            coefficient = round(
                inner(basis[i], ortho[j]) / inner(ortho[j], ortho[j]))
            basis[i] = sub(basis[i], mult(coefficient, basis[j]))
    for i in range(len(basis) - 1):
        product = inner(ortho[i], ortho[i])
        lovasz_inner = add(
            mult(inner(basis[i + 1], ortho[i]) /
                 product, ortho[i]), ortho[i + 1])
        if delta * product > inner(lovasz_inner, lovasz_inner):
            basis[i], basis[i + 1] = basis[i + 1], basis[i]
            return lll(basis, delta)
    return basis


def poly_pad(f, n):
    g = [s.S(0.0)] * n
    for i in range(len(f)):
        g[i] = f[i]
    return g


def poly_mul(f1, f2):
    res = [s.S(0.0)] * (len(f1) + len(f2) - 1)
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
            n = n / 2
        else:
            y = poly_mul(x, y)
            x = poly_mul(x, x)
            n = (n - 1) / 2
    return poly_mul(x, y)


def generate(p, N, h, k, X):
    poly_matrix = [[s.S(0) for _ in range(h * k)] for _ in range(h * k)]
    n = [N]
    x = [s.S(0), s.S(1)]

    for i in range(1, (h * k) + 1):
        v = s.floor((i - 1) / k)
        v = max(v, 0)
        u = (i - 1) - (k * v)
        q = poly_mul(
            poly_mul(poly_pow(n, h - 1 - v), poly_pow(x, u)), poly_pow(p, v))
        q = poly_pad(q, h * k)
        for j in range(h * k):
            poly_matrix[i - 1][j] = s.S(q[j] * s.Pow(X, j))

    return poly_matrix


def run_copper_smith(f, N, h, k):
    X = s.S(s.ceiling((s.Pow(2, - 1 / 2) * s.Pow(h * k, -1 / (h * k - 1))) *
                      s.Pow(N, (h - 1) / (h * k - 1))) - 1)
    gen = lll(generate(f, N, h, k, X), s.S(0.75))
    final_poly = [gen[0][i] / s.Pow(X, i) for i in range(len(gen[0]))][::-1]

    roots = s.roots(final_poly)
    for r in roots:
        print(r)

    return roots

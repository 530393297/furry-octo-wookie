import numpy as np
import numpy.polynomial.polynomial as poly

def matrix_mult(a, b):
	c = [0] * len(a[0])

	for i in range(len(a)):
		for j in range(len(a[0])):
			c[j] = c[j] + (a[i][j] * b[i])

	return c

def vector_sub(a, b):
	return [a[i] - b[i] for i in range(len(a))]

def vector_add(a, b):
	return [a[i] + b[i] for i in range(len(a))]

def scalar_mult(a, b):
	return [i * a for i in b]

def inner_product(a, b):
	c = 0

	for i in range(len(a)):
		c = c + a[i] * b[i]

	return c

def generate(N,C,e,m):
	f = [0]*(e+1)
	f[e] = 1
	f[0] = -C
	
	x = [0,1]
	g = []
	for u in range(e):
		temp = []
		for v in range(m+1):
			temp.append(poly.polymul([pow(N,m-v)], poly.polymul(poly.polypow(f,v), poly.polypow(x,u))).tolist())
		g.append(temp)
	return g

def polypartial(f,B):
	for i in range(len(f)):
		f[i] = f[i]*pow(B,i)
	return f

def unpolypartial(f,B):
	for i in range(len(f)):
		f[i] = f[i]/pow(B,i)
	return f

def poly_pad(f, n):
	g = [0]*n
	for i in range(len(f)):
		g[i] = f[i]
	return g

def generate_basis(g, d, m):
	n = d * (m + 1)
	B = pow(pow(n, 1 / d), m / (m + 1))

	out = []

	for u in range(d):
		for v in range(m + 1):
			out.append(poly_pad(polypartial(g[u][v], B), n)) 

	return out


def gram_schmidt(B):
	Bs = list(B)
	for i in range(len(B)):
		for j in range(i):
			Bs[i] = vector_sub(Bs[i], scalar_mult(inner_product(B[i], Bs[j]) / inner_product(B[j], Bs[j]), Bs[j]))
	
	return Bs

def LLL(B, delta):
	Bs = gram_schmidt(B)

	for i in range(1, len(B)):
		for j in reversed(range (i)):
			c = round(inner_product(B[i], Bs[j]) / inner_product(Bs[j], Bs[j]))
			B[i] = vector_sub(B[i], scalar_mult(c, B[j]))
	for i in range(len(B) - 1):
		if delta * inner_product(Bs[i], Bs[i]) > inner_product(vector_add(scalar_mult(inner_product(B[i + 1], Bs[i]) /  inner_product(Bs[i], Bs[i]), Bs[i]), Bs[i+1]), vector_add(scalar_mult(inner_product(B[i + 1], Bs[i]) /  inner_product(Bs[i], Bs[i]), Bs[i]), Bs[i+1])):
			temp = B[i]
			B[i] = B[i + 1]
			B[i + 1] = temp
			return LLL(B, delta)
	return B	


c = [[1,1,1], [-1,0,2], [3,5,6]]

cc = np.array([ (1,1,1), (-1,0,2), (3,5,6) ], dtype='f')

N = 33
e = 7
m = 4
M = 2
C = 29
d=e
l = LLL(generate_basis(generate(N,C,d,m), d, m), 0.75)
n = d * (m + 1)
B = pow(pow(n, 1 / d), m / (m + 1))
ll = unpolypartial(l[0],B)
print(l)
print(poly.polyval(M,ll[0])%N)
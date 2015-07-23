import decimal as d
import scipy.optimize as scipy
import functools as functools
import random
import operator
from rsa import random_prime

def inner(x, y):
	return sum(map(operator.mul, x, y))

def sub(x, y):
	return list(map(operator.sub, x, y))

def add(x, y):
	return list(map(operator.add, x, y))

def mult(x, y):
	return list(i * x for i in y)

def gram_schmidt(b):
    B = b.copy()
    for i in range(len(b)):
            for j in range(i):
                    ip = inner(b[i], B[j]) / inner(B[j], B[j])
                    B[i] = sub(B[i], mult(ip, B[j]))
    return B

def lll(b, delt):    
    B = gram_schmidt(b)

    for i in range(1, len(b)):
            for j in reversed(range (i)):
                    c = round(inner(b[i], B[j]) / inner(B[j], B[j]))
                    b[i] = sub(b[i], mult(c, b[j]))
    for i in range(len(b) - 1):
            p = inner(B[i], B[i])
            y = add(mult(inner(b[i + 1], B[i]) / p, B[i]), B[i + 1])
            if(delt * p > inner(y, y)):
                    b[i], b[i + 1] = b[i + 1], b[i]
                    return lll(b, delt)
    return b

def copper_smith(a, N):
	n = len(a)
	X = d.Decimal(pow(N, 1 / n) / ((n - 1) * 2))
	M = [[d.Decimal(0.0) for x in range(len(a))] for x in range(len(a))] 
	
	for i in range(len(a)):
		M[i][i] = pow(X,len(a)-i)*N
	for i in range(len(a)):
		M[0][i] = pow(X, len(a) - i) * a[i]
	
	V = lll(M, d.Decimal(0.75))
	v = V[0]
	return [v[i] / pow(X, (len(a)-1 - i)) for i in range(len(a))]

def func(f, x):
	total = 0
	for i in range(len(f)):
		total = total + pow(x, (len(f) - 1) - i) * round(f[i])
	return total

#N = 1147
#a = 1000
#b = 1000
#print("poly test 2:")
#print(CopDeg2(a,b,N))
#print("poly real answers:")
#print("322,507,787,972")

#N = 1147
#a = 1000
#b = 1000
#c = 1000
#print("poly test 3:")
#print(CopDeg3(a,b,c,N))
#print("poly real answers:")
#print("443,598,660")

#N = 1147
#a = 1000
#b = 1000
#c = 259
#f = CopDeg3(a,b,c,N)
#print(f)
#funcf = functools.partial(func3,f)
#X = pow(N, 1 / 4.0) / 6
#print(funcf(10))
#print(scipy.brentq(funcf,0,X))


#N = 2122840968903324034467344329510307845524745715398875789936591447337206598081
#a = 3*pow(2,500)
#b = 3*pow(2,500)*pow(2,500)
#c = pow(2,500)*pow(2,500)*pow(2,500)-1792963459690600192400355988468130271248171381827462749870651408943993480816

#A = [1,a,b,c]
#print("RSA test:")
#print(copper_smith(A,N))


#p = random_prime(1000000)
#q = random_prime(1000000)
#Nn = p*q
#e = 3
#n = 1
#TN = Nn >> 1
#while TN != 0:
#    n = n + 1
#    TN = TN >> 1
#print(round(pow(Nn, 1 / 4) / 6))
#x0 = random.randint(0,round(pow(Nn, 1 / 4) / 6))
#a = random.randint(0,Nn)
#b = random.randint(0,Nn)
#c = (-x0*x0*x0-a*x0*x0-b*x0)%Nn
#A = [1,a,b,c]
#print(x0)
#v = copper_smith(A,Nn)
#funcv = functools.partial(func,v)
#print(scipy.brentq(funcv,0,x0+100))
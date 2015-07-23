import numpy as np
import decimal
import scipy.optimize as scipy
import functools as functools
import random

def gram_schmidt(b):
        B = b.copy()
        for i in range(len(b)):
                for j in range(i):
                        B[i] = B[i] - ((np.inner(b[i], B[j]) / np.inner(B[j], B[j])) * B[j])
        return B

def lll(b, delta):
        B = gram_schmidt(b)

        for i in range(1, len(b)):
                for j in reversed(range (i)):
                        c = round(np.inner(b[i], B[j]) / np.inner(B[j], B[j]))
                        b[i] = b[i] - (b[j] * c)
        for i in range(len(b) - 1):
                y = (np.inner(b[i + 1], B[i]) /  np.inner(B[i], B[i]) * B[i]) + B[i + 1]
                if(delta * np.inner(B[i], B[i]) > np.inner(y, y)):
                        b[[i, i + 1]] = b[[i + 1, i]]
                        return lll(b, delta)
        return b

def CopDeg2(a,b,N):
        n = 1
        TN = N >> 1
        while TN != 0:
                n = n + 1
                TN = TN >> 1
        X=pow(2, (n // 3 - 5))
        M = np.array([(X * X, a * X, b), (0, N * X, 0), (0 ,0 ,N)])

        V = lll(M, 0.75)
        v = V[0]
        return [v[i] / pow(X, (2 - i)) for i in range(3)]

def CopDeg3(a,b,c,N):
        X = pow(N, 1 / 4.0) / 6
        M = np.array([(X * X * X, a * X * X, b * X, c), (0, N * X * X, 0, 0), (0 ,0 ,N * X, 0), (0, 0, 0, N)])

        V = lll(M, 0.75)
        v = V[0]
        return [v[i] / pow(X, (3 - i)) for i in range(4)]


def func3(f,x):
	return x*x*x*f[0]+x*x*f[1]+f[2]*x+f[3]

def func2(f,x):
	return x*x*f[0]+x*f[1]+f[2]

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
#print("RSA test:")
#print(CopDeg3(a,b,c,N))


p = 106859
q = 106661
Nn = p*q
e = 3
n = 1
TN = Nn >> 1
while TN != 0:
    n = n + 1
    TN = TN >> 1
print(pow(2,n//3-10))
x0 = random.randint(0,pow(2,n//3-10))
a = random.randint(0,Nn)
b = -x0*x0-a*x0%Nn
print(x0)
v = CopDeg2(a,b,Nn)
funcv = functools.partial(func2,v)
print(scipy.brentq(funcv,0,pow(2,n//3-10)))
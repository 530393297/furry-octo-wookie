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
                if(delta * np.inner(B[i], B[i]) > np.inner(y,y)):
                        b[[i, i + 1]] = b[[i + 1, i]]
                        return lll(b, delta)
        return b

def copper_smith(a, N):
	X = pow(N, 1 / len(a)) / (len(a) * 2)

	M = [[]]*len(a)
	for i in range(len(a)):
		M[i] = [0]*len(a)
		M[i][i] = pow(X,len(a)-i)*N
	for i in range(len(a)):
		M[0][i] = pow(X,len(a)-i)*a[i]
	print(M)
	V = lll(np.array(M), 0.75)
	v = V[0]

	return [v[i] / pow(X, (len(a)-1 - i)) for i in range(len(a))]

def func(f, x):
	total = 0
	for i in range(len(f)):
		total = total + pow(x, (len(f) - 1) - i) * round(f[i])
	return total

cc = np.array([(1,1,1), (-1,0,2), (3,5,6)],dtype=np.float32)

print(lll(cc, 0.75))

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


#p = 106859
#q = 106661
#Nn = p*q
#e = 3
#n = 1
#TN = Nn >> 1
#while TN != 0:
#    n = n + 1
#    TN = TN >> 1
#print(pow(2,n//3-10))
#x0 = random.randint(0,pow(2,n//3-10))
#a = random.randint(0,Nn)
#b = -x0*x0-a*x0%Nn
#A = [1,a,b]
#print(x0)
#v = copper_smith(A,Nn)
#funcv = functools.partial(func,v)
#print(scipy.brentq(funcv,0,pow(2,n//3-10)))

#print(func([1,1,1],1))
#fff = functools.partial(func,[1,1,1])
#print(fff(1))
#print(copper_smith([1,2,3],5))
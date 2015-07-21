def numpy_gram_schmidt(B):
	print("GS")
	Bs = B.copy()
	for i in range(len(B)):
		for j in range(i):
			Bs[i] = np.float64(Bs[i] - np.float64(((np.float64(np.inner(B[i], Bs[j]) / np.inner(B[j], Bs[j]))) * Bs[j])))

	print(Bs)
	return Bs

def numpy_lll(B, delta):
	print ("LLL - B")
	print(B)
	Bs = numpy_gram_schmidt(B)

	for i in range(1, len(B)):
		for j in reversed(range (i)):
			c = round(np.inner(B[i], Bs[j]) / np.inner(Bs[j], Bs[j]))
			xx = np.float64(np.inner(B[i], Bs[j]) / np.inner(Bs[j], Bs[j]))
			print(xx)
			B[i] = np.float64(B[i] - np.float64((B[j] * c)))
			print(B[i])
	for i in range(len(B) - 1):
		x = delta * np.inner(Bs[i], Bs[i])
		y = (np.inner(B[i + 1], Bs[i]) /  np.inner(Bs[i], Bs[i]) * Bs[i]) + Bs[i+1]
		if(x > np.inner(y,y)):
			print(i)
			B[[i,i+1]] = B[[i+1,i]]
			print(B)
			return numpy_lll(B, delta)
	return B
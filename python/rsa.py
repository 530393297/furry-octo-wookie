

def is_prime(x):
	for i in range(3, len(small_primes)):
		p = small_primes[i]
		q = x // p
		if q < p:
			return True
		if x == q * p:
			return False
	i = 31
	while True:
		q = x // i
		if q < i:
			return True
		if x == q * i:
			return False
		i = i + 6
		if q < i:
			return True
		if x == q * i:
			return False
		i = i + 4
		if q < i:
			return True
		if x == q * i:
			return False
		i = i + 2
		if q < i:
			return True
		if x == q * i:
			return False
		i = i + 4
		if q < i:
			return True
		if x == q * i:
			return False
		i = i + 2
		if q < i:
			return True
		if x == q * i:
			return False
		i = i + 4
		if q < i:
			return True
		if x == q * i:
			return False
		i = i + 6
		if q < i:
			return True
		if x == q * i:
			return False
		i = i + 2
	
def next_prime(n):
	L = 30
	if n <= small_primes[len(small_primes) - 1]:
		return small_primes[bisect.bisect_left(small_primes, n)]
	M = len(indices)
	k0 = n // L
	ind = bisect.bisect_left(indices, n - k0 * L) 
	n = L * k0 + indices[ind]
	while not is_prime(n):
		ind = ind + 1
		if ind == M:
			k0 = k0 + 1
			ind = 0
		n = L * k0 + indices[ind]
	return n
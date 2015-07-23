import random
from fractions import gcd

def composite(a, d, n, s):
    if pow(a, d, n) == 1:
        return False
    for i in range(s):
        if pow(a, 2**i * d, n) == n - 1:
            return False
    return True


def is_prime(n, precis=16):
    if n in PRIMES or n in (0, 1):
        return True
    if any((n % p) == 0 for p in PRIMES):
        return False
    d, s = n - 1, 0
    while not d % 2:
        d, s = d >> 1, s + 1

    if n < 1373653:
        return not any(composite(a, d, n, s) for a in (2, 3))
    if n < 25326001:
        return not any(composite(a, d, n, s) for a in (2, 3, 5))
    if n < 118670087467:
        if n == 3215031751:
            return False
        return not any(composite(a, d, n, s) for a in (2, 3, 5, 7))
    if n < 2152302898747:
        return not any(composite(a, d, n, s) for a in (2, 3, 5, 7, 11))
    if n < 3474749660383:
        return not any(composite(a, d, n, s) for a in (2, 3, 5, 7, 11, 13))
    if n < 341550071728321:
        return not any(composite(a, d, n, s) for a in (2, 3, 5, 7, 11, 13, 17))
    return not any(composite(a, d, n, s)
                   for a in PRIMES[:precis])

PRIMES = [2, 3]
PRIMES += [y for y in range(5, 1000, 2) if is_prime(y)]

def random_prime(n):
    x = 4
    while is_prime(x) == False:
        x = random.randint(n / 2, n)
    return x


def keygen(n=256):
    while True:
        p = random_prime(pow(2, (n // 2)))
        q = random_prime(pow(2, (n // 2)))
        e = 3
        if gcd(e, (p - 1) * (q - 1)) == 1:
            break
    Nn = p * q
    return {'p': p, 'q': q, 'Nn': Nn, 'e': e}

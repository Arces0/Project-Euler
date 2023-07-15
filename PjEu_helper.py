import math
import numpy as np
from random import randint
from sys import stdin
IN = lambda: stdin.readline()

# helper functions

isInt = lambda r: r == int(r)

def value_inc(dct, key, initial = 0, increment = 1):						# assert type(dct) = {key: int}
	'''
	Input - dct:dict, key:Any
	Output - N/A

	If key is not in dct, set dct[key] = initial
	Else, dct[key] += increment
	'''
	if key not in dct: dct[key] = initial
	dct[key] += increment

def polyeval(N, coeffL, x):
	'''
	Input - N:int; coeffL:list of numeric, length N+1; x:numeric
	Output - evaluation of polynomial of degree N at x, with coefficients from coeffL
	Raise assertionError if len(coeffL) != N+1

	Assume p(x) = L[0]x^N + L[1]x^(N-1) + ... + L[N-1]x + L[N]
	Applied Horner's method
	'''
	assert(len(coeffL) == N+1)
	rtn = coeffL[0]
	for i in range(1, N+1):
		rtn *= x; rtn += coeffL[i]
	return rtn

def next_word(w):
	'''
	Input - w:str without whitespace
	Output - Next word of w in length and alphabetical order

	Ex) a -> b, z -> aa, buzz -> bvaa
	'''
	zstrip = w.rstrip('z')
	return zstrip[:-1] + chr(ord(zstrip[-1])+1) + 'a'*(len(w)-len(zstrip)) if zstrip else 'a'*(len(w)+1)

def sumdigit(n):
	'''
	Input - n:int
	Output - sum of digits of n
	'''
	s = str(n)
	return sum(int(d) * s.count(d) for d in "123456789")

def isPrime(n):
	'''
	Input - n:int
	Output - True if n is prime, False otherwise
	'''
	if n == 2 or n == 3: return True
	if n%2 == 0 or n<2: return False
	for _ in range(3, int(n**0.5)+1, 2):
		if n%_ == 0: return False
	return True

def primes(n):
	'''
	Input - n: Int
	Output - list of primes no greater than n
	'''
	if n < 2: return []
	# odd-sieve: index i -> integer 2*i + 1 except 0 -> 2
	sieve = [True] * ((n+1)//2)
	for i in range(1, int(math.sqrt(n)/2)+1):
		if sieve[i]:
			for j in range(2*i*(i+1), (n+1)//2, 2*i+1):
				sieve[j] = False
	return [2] + [2*i+1 for i in range(1, (n+1)//2) if sieve[i]]

def gcd(a, b):									# gcd for positive integers
	if b == 0: return a
	else:
		if a < 0: a = -a
		if b < 0: b = -b
		return gcd(b, a%b)

def Factorial(n):								# List of k! from 0 to n
	res = [1]
	for k in range(1, n+1):
		res.append(res[-1]*k)
	return res

def multinomial(L):
	factL = Factorial(sum(L))
	val = factL[-1]
	for k in L: val //= factL[k]
	return val

def Power(p, e):
	if e == 0: return 1
	elif e == 1: return p
	else:
		temp = Power(p, e//2)
		if e % 2 == 0: res = temp * temp
		else: res = temp * temp * p
		return res

def isPrimeprob(n, k=10): # miller-rabin
	if n < 2: return False
	if n < 800000: return isPrime(n)
	for p in [2,3,5,7,11,13,17,19,23,29]:
		if n % p == 0: return n == p
	s, d = 0, n-1
	while d % 2 == 0:
		s, d = s+1, d//2
	for i in range(k):
		x = pow(randint(2, n-1), d, n)
		if x == 1 or x == n-1: continue
		for r in range(1, s):
			x = (x * x) % n
			if x == 1: return False
			if x == n-1: break
		else: return False
	return True

# def factors(n, b2=-1, b1=10000): # 2,3,5-wheel, then rho
# 	def gcd(a,b): # euclid's algorithm
# 		if b == 0: return a
# 		return gcd(b, a%b)
# 	if -1 <= n <= 1: return dict()
# 	elif n < -1: return factors(-n)
# 	wheel = [1,2,2,4,2,4,2,4,6,2,6]
# 	w, f, fs = 0, 2, {}
# 	while f*f <= n and f < b1:
# 		while n % f == 0:
# 			value_inc(fs, f)
# 			n //= f
# 		f, w = f + wheel[w], w+1
# 		if w == 11: w = 3
# 	if n == 1: return fs
# 	h, t, g, c = 1, 1, 1, 1
# 	while not isPrimeprob(n):
# 		while b2 != 0 and g == 1:
# 			h = (h*h+c)%n # the hare runs
# 			h = (h*h+c)%n # twice as fast
# 			t = (t*t+c)%n # as the tortoise
# 			g = gcd(t-h, n); b2 -= 1
# 		if b2 == 0: return fs
# 		if isPrimeprob(g):
# 			while n % g == 0:
# 				value_inc(fs, g)
# 				n //= g
# 		h, t, g, c = 1, 1, 1, c+1
# 	value_inc(fs, n)
# 	return fs

def factors(N):
	res = {}
	while N % 2 == 0:
		value_inc(res, 2)
		N //= 2
	while N % 3 == 0:
		value_inc(res, 3)
		N //= 3
	while N % 5 == 0:
		value_inc(res, 5)
		N //= 5
	
	# wheel factorization
	inc = [4,2,4,2,4,6,2,6]
	p = 7; idx = 0
	while p * p <= N:
		if N % p == 0:
			value_inc(res, p)
			N //= p
		else:
			p += inc[idx]
			idx = (idx+1) % 8

	if N > 1: value_inc(res, N)
	return res

def divisorsum(n):								# sum of divisior including itself
	pfactorD = factors(n)
	return math.prod((p**(pfactorD[p]+1) -1)//(p-1) for p in pfactorD.keys())

def divisorsum2(n):
	res = 0
	for d in range(1, int(n**0.5)+1):
		if n%d == 0:
			res += d + n // d
			if d == n // d: res -= d
	return res

def Euphi(n, valonly = True):
	res = n; factorD = factors(n)
	for q in factorD:
		res //= q
		res *= q-1
	if valonly: return res
	else: return res, factorD.keys()

def contfracsqr(n):
	if isInt(math.sqrt(n)): return [int(math.sqrt(n)), tuple()]
	b,c = 0,1; k0 = int(math.sqrt(n)); contfracL = []	#the coeff. of sqrt(n) is always 1 in cont.frac. (deriven from last #)
	while 1:
		k = (k0+b)//c
		if not contfracL: contfracL.append(k); contfracL.append(tuple())
		else: contfracL[1] += (k,)
		if k == 2*int(k0): break						#cont.frac. ends iff k == 2k0: sqrt(n)-k0 = ... = k0+1/...+1/(sqrt(n)+k0) (<-2*k0 + sqrt(n)-k0)
		m = b-c*k; r = (n-m*m)//c						#proven that c | (n-m*m)
		b, c = -m, r
	return contfracL
'''
Input - A:np.array, n:int
Output - power of A to n
'''
mtxpow = lambda A, n: np.linalg.matrix_power(A, n)
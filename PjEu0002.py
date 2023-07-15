# Copyright 2023. Arces0. All rights reserved.
# https://github.com/Arces0/Project-Euler

from PjEu_helper import *

def PjEu002(N: int) -> int:
	'''
	Input:	integer N
	Output:	sum of even-valued Fibonacci numbers less than N
		if no such numbers exist, return 0
	'''
	rtn = 0
	# EF_k: k-th even-valued Fibonacci numbers
	# Then EF_k+1 = 4*EF_k + EF_k-1 where EF_1 = 2, EF_0 = 0
	m, n = 0, 2
	while n < N:
		m, n = n, 4*n+m
		rtn += m
	return rtn

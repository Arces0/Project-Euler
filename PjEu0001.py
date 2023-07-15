# Copyright 2023. Arces0. All rights reserved.
# https://github.com/Arces0/Project-Euler

from PjEu_helper import *

def PjEu001(N: int) -> int:
	'''
	Input:	nonnegtive integer N
	Output:	sum of every multiples of 3 or 5 less than N
	'''
	assert(N >= 0)
	def summultiple(K, m):
		'''
		Input:	integer K, integer m < K
		Output:	sum of every multiples of m less than K
		'''
		k = (K-1)//m
		return m * k*(k+1)//2

	return summultiple(N, 3) + summultiple(N, 5) - summultiple(N, 15)

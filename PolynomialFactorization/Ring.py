import math
class Ring:

	def __init__(self, ring_info):

		if 


	def _isPrime(n):
		if n==2:
			return True
		for i in range(math.sqrt(n)):
			if n%i == 0:
				return False

		return True
import math
from random import randint
from copy import copy, deepcopy

def isPrime(n):
	if n == 2:
		return True
	for i in range(2, int(math.sqrt(n))):
		if n%i==0:
			return False
	return True



class Polynomial:

	#def __init__(self, polynomial_ring, coeffs):
	#	self.polynomial_ring = polynomial_ring
	#	self.coefficients = coeffs

	def __init__(self, coeffs, characteristic, indeterminate_symbol):
		assert characteristic >= 0
		assert isPrime(characteristic)
		self.characteristic = characteristic
		self.coefficients = coeffs

		self.coefficients = [coeff % characteristic for coeff in self.coefficients]

		while self.coefficients[-1]==0 and len(self.coefficients)>1:
			self.coefficients = self.coefficients[:-1]


		self.indeterminate_symbol = indeterminate_symbol
		if len(self.coefficients) == 1:
			if self.coefficients[0] == 0:
				self.degree = -1
			else:
				self.degree = 0
		else:
			self.degree = len(self.coefficients)-1

	def __add__(self, other):
		if isinstance(other, int):
			coeffs = self.coefficients.copy()
			coeffs[0] += other
			return Polynomial(coeffs, self.characteristic, self.indeterminate_symbol)
		assert self._compatible(other)
		smaller = []
		bigger = []

		if self.degree > other.degree:
			smaller = other.coefficients
			bigger = self.coefficients
		else:
			smaller = self.coefficients
			bigger = other.coefficients

		sum_coeffs = [smaller[i] + bigger[i] for i in range(len(smaller))]
		if len(bigger)>len(smaller):
			sum_coeffs += bigger[len(smaller):]

		return Polynomial(sum_coeffs, self.characteristic, self.indeterminate_symbol)

	def __rmul__(self, scalar):
		return Polynomial([scalar*coeff for coeff in self.coefficients], self.characteristic, self.indeterminate_symbol)


	def __sub__(self, other):
		return self + ((-1)*other)

	def __mul__(self, other):
		assert self._compatible(other)
		if self.degree == -1 or other.degree == -1:
			return Polynomial([0], self.characteristic, self.indeterminate_symbol)

		bigger = []
		smaller = []
		if self.degree > other.degree:
			bigger = self.coefficients
			smaller = other.coefficients + ([0]*(self.degree - other.degree))
		else:
			bigger = other.coefficients
			smaller = self.coefficients + ([0]*(other.degree - self.degree))

		prod_coeffs = [0 for _ in range(self.degree + other.degree+1)]
		for i in range(len(self.coefficients)):
			for j in range(len(other.coefficients)):
				prod_coeffs[i+j] += self.coefficients[i]*other.coefficients[j]

		return Polynomial(prod_coeffs, self.characteristic, self.indeterminate_symbol)

	def __truediv__(self, other):
		return self.division_with_remainder(other)[0]

	def __mod__(self, other):
		return self.division_with_remainder(other)[1]

	def division_with_remainder(self, other):
		assert self._compatible(other)
		p = self.characteristic
		x = Polynomial.gen(p, self.indeterminate_symbol)
		q = Polynomial.zero(p, self.indeterminate_symbol)

		if self.degree < other.degree:
			return (Polynomial.zero(p, self.indeterminate_symbol), copy(self))
		else:
			r = deepcopy(self)
			inv = modinv(other.coefficients[other.degree], p)
			while r.degree >= other.degree:
				# print('g = ' + str(other))
				# print('r = ' + str(r))
				# print('q = ' + str(q))
				# print('r.degree = ' + str(r.degree))
				# print('other.degree = ' + str(other.degree))
				# if r.degree - other.degree > 1000:
				# 	print('fucking why')
				# 	return
				qterm = (r.coefficients[r.degree])*inv * (x**(r.degree - other.degree))
				q += qterm
				r -= qterm * other
			return (q, r)

	#TODO: square n add?
	# powmod
	# def __pow__(self, exponent):
	# 	if exponent == 0:
	# 		return Polynomial([1], self.characteristic, self.indeterminate_symbol)
	# 	elif exponent == 1:
	# 		return Polynomial(self.coefficients, self.characteristic, self.indeterminate_symbol)
	# 	else:
	# 		return self * (self**(exponent-1))


	def __pow__(self, exponent):
		if exponent < 0:
			return Polynomial.zero(self.characteristic, self.indeterminate_symbol)
		if exponent == 0:
			return Polynomial.one(self.characteristic, self.indeterminate_symbol)
		if exponent == 1:
			return copy(self)

		# is it x?
		if self.coefficients == [0,1]:
			return Polynomial([0]*exponent + [1], self.characteristic, self.indeterminate_symbol)

		# is it a power of x?
		if self.coefficients[:self.degree] == ([0]*self.degree) and self.coefficients[-1] == 1:
			return Polynomial([0]*(self.degree * exponent) + [1], self.characteristic, self.indeterminate_symbol)
		ret = copy(self)
		for _ in range(exponent-1):
			ret *= self
		return ret




	def powmod(self, exponent, modulus):
		assert self._compatible(modulus)
		if modulus.degree == 0:
			return Polynomial.zero(self.characteristic, self.indeterminate_symbol)
		ret = Polynomial.one(self.characteristic, self.indeterminate_symbol)
		base = self % modulus
		while exponent > 0:
			if exponent % 2 == 1:
				ret = (ret * base) % modulus
			exponent = exponent >> 1
			base = (base*base) % modulus
		return ret


	@staticmethod
	def zero(characteristic, indeterminate_symbol):
		return Polynomial([0], characteristic, indeterminate_symbol)


	@staticmethod
	def gen(characteristic, indeterminate_symbol):
		return Polynomial([0, 1], characteristic, indeterminate_symbol)

	@staticmethod
	def one(characteristic, indeterminate_symbol):
		return Polynomial([1], characteristic, indeterminate_symbol)


	def __str__(self):
		ret = ''
		if self.degree == 0 or self.degree == -1:
			return str(self.coefficients[0])
		for i in range(self.degree + 1):
			coeff = self.coefficients[self.degree-i]
			powerstr = ''
			if self.degree-i== 1:
				powerstr += self.indeterminate_symbol
			elif self.degree-i >1:
				powerstr += self.indeterminate_symbol + '^' + str(self.degree-i)
			if coeff != 0:
				if coeff == 1:
					if i == 0:
						ret += powerstr
					elif i == self.degree:
						ret += '+1'
					else:
						ret += '+' + powerstr
				else:
					if i==0:
						ret += str(coeff) + powerstr
					else:
						ret += '+' + str(coeff)+powerstr
		return ret

	def derivative(self):
		coeffs = [(n+1)*coeff for (n, coeff) in enumerate(self.coefficients[1:])]
		return Polynomial(coeffs, self.characteristic, self.indeterminate_symbol)

	def isSquarefree(self):
		return gcd(self, self.derivative()).coefficients == [1]
				
	def _compatible(self, other):
		return self.characteristic == other.characteristic and self.indeterminate_symbol == other.indeterminate_symbol


def egcd(a,b):
	if a==0:
		return (b, 0, 1)
	else:
		g, x, y = egcd(b%a, a)
		return (g, y-(b//a)*x, x)

def modinv(b, n):
	g, x, _ = egcd(b,n)
	assert g==1
	return x%n

#(quotient, remainder)
def division(f,g):
	assert f.characteristic == g.characteristic
	assert f.indeterminate_symbol == g.indeterminate_symbol
	x = Polynomial([0,1], f.characteristic, f.indeterminate_symbol)
	q = Polynomial([0], f.characteristic, f.indeterminate_symbol)
	p = f.characteristic

	if f.degree < g.degree:
		return (Polynomial([0], p, f.indeterminate_symbol), f)
	else:
		r = Polynomial(f.coefficients.copy(), p, f.indeterminate_symbol)
		inv = modinv(g.coefficients[-1], p)
		while r.degree >= g.degree:
			qterm = (r.coefficients[-1])*inv * (x**(r.degree - g.degree))
			q += qterm
			r -= qterm * g
		return (q, r)

#polynomial gcd
def gcd(f,g):
	assert f.characteristic == g.characteristic
	p = f.characteristic
	if f.degree >= g.degree:
		r0 = modinv(f.coefficients[-1], p) * f
		r1 = modinv(g.coefficients[-1], p)*g
	else:
		r0 = modinv(g.coefficients[-1], p)*g
		#print('f = ' + str(f))
		r1 = modinv(f.coefficients[-1], p)*f

	r = r1
	rm = r0
	while r.degree is not -1:
		(qp, rp) = division(rm, r)
		if rp.degree is not -1:
			rp = modinv(rp.coefficients[-1], p)*rp
		rm = r
		r = rp

	return rm

#returns (h, s, t) with
# sf+tg=h
# h a gcd of f and g not necessarily monic
def poly_egcd(f,g):
	rm = f
	r = g
	sm = Polynomial([1], f.characteristic, f.indeterminate_symbol)
	tm = Polynomial([0], f.characteristic, f.indeterminate_symbol)
	s = Polynomial([0], f.characteristic, f.indeterminate_symbol)
	t = Polynomial([1], f.characteristic, f.indeterminate_symbol)

	while r.degree is not -1:
		(q, rp) = gcd(rm, r)
		sp = sm - (q*s)
		tp = tm - (q*t)
		rm = r
		sm = s
		tm = t
		r = rp

	return (rm, sm, tm)


def DDD(f):
	assert f.isSquarefree()

	fd = f
	x = Polynomial([0,1], f.characteristic, f.indeterminate_symbol)
	gd = Polynomial([0, 1], f.characteristic, f.indeterminate_symbol)
	ret = []
	p = f.characteristic

	d = 0
	while fd.degree > 0:
		d += 1
		gd = gd**p
		#gd = Polynomial([0]*(p**d) + [1], f.characteristic, f.indeterminate_symbol)
		hd = gcd(fd, gd - x)
		ret.append((d, hd))
		fd = division(fd, hd)[0]

	return ret


def TrialSplit(h, degree):
	#print('TrialSplit(' + str(h) + ', ' + str(degree) + ')')
	ret = Polynomial([0], h.characteristic, h.indeterminate_symbol)
	if h.degree == 1:
		return ret
	else:
		g1 = Polynomial([0], h.characteristic, h.indeterminate_symbol)
		while g1.degree < 1:
			g1 = randomPoly(h.degree - 1, h.characteristic, h.indeterminate_symbol)
		g2 = gcd(g1, h)
		if g2.degree is not 0:
			ret = g2
		else:
			e = ((h.characteristic**degree) - 1) // 2
			#g3 = division(g1**e, h)[1]
			g3 = g1.powmod(e, h)
			#print('g1 = ' + str(g1))
			#print('g3 = ' + str(g3))
			if g3.coefficients == [1]:
				g4 = gcd(g3+1, h)
			else:
				g4 = gcd(g3-1, h)
			if 0<g4.degree < h.degree:
				ret = g4
			else:
				ret = Polynomial([0], h.characteristic, h.indeterminate_symbol)
	return ret


def Split(h, degree, s):
	ret = Polynomial([0], h.characteristic, h.indeterminate_symbol)
	k = 0
	while ret.degree == -1 and k<s:
		ret = TrialSplit(h, degree)
		k += 1
	return ret

def randomPoly(degree, characteristic, indeterminate_symbol):
	coeffs = [randint(0, characteristic - 1) for _ in range(degree + 1)]
	return Polynomial(coeffs, characteristic, indeterminate_symbol)

def EDD(h, degree, s):
	factorlist = []

	def EDD_helper(h1, degree1, s1):
		#print('EDD_h(' + str(h1) + ', ' + str(degree1) + ', ' + str(s1) + ')')
		g1 = Split(h1, degree1, s1)
		if g1.degree <1:
			factorlist.append(h1)
		else:
			EDD_helper(g1, degree1, s1)
			EDD_helper(h1/g1, degree1, s1)

	EDD_helper(h, degree, s)

	return factorlist


def Factor(f, s):
	x = Polynomial.gen(f.characteristic, f.indeterminate_symbol)
	complete = []
	d = 0
	fd = deepcopy(f)
	gd = Polynomial.gen(f.characteristic, f.indeterminate_symbol)
	while fd.coefficients != [1]:
		#print('d = ' + str(d))
		d += 1
		gd = gd**f.characteristic
		hd = gcd(fd, gd-x)
		#print('hd = ' + str(hd))
		if hd.coefficients != [1]:
			factorlist = EDD(hd, d, s)
			for k in factorlist:
				m = 0
				while (fd%k).degree == -1:
					m += 1
					fd /= k
				complete.append((k, m))
	return complete




####################################################

x = Polynomial([0,1], 3, 'x')

h = (x**12) + (x**11) + 3*(x**10) + 4*(x**8) + 2*(x**7) + 4*(x**5) + 2*(x**4) + (x**3) + 2*(x**2) + (4*x) + 4
#f = (x**12) + (x**8) + (x**6) + (x**4) + (x**2) + 1
f = (x**10) + (x**8) + (x**7) + 2*(x**6) + 2*(x**5) + (x**4) + (x**3) + (x**2) + (2*x)

# modulus = (x**25)-x

# print(h.powmod(50, modulus))

for pair in Factor(f, 6):
	print('(' + str(pair[0]) + ')^' + str(pair[1]))
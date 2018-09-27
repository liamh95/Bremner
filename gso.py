import numpy as np
from fraction import Fraction
# X = MY
# X is original matrix. row span vectors
# M is gso coeffs
# Y is the gso vectors
#
# matrices and fractions should have elements of type fraction.Fraction

X = [[3, -5, 9, 4], [1, -3, 5, 4], [2, -1, -1, 6], [3, 2, 2, 5]]
X = [[Fraction(n) for n in v] for v in X]
M = [[] for _ in X]
Y = []

def dot(u, v):
    res = Fraction(0)
    for i in range(len(u)):
        res += u[i]*v[i]
    return res


for i in range(len(X)):
    yi = np.array(X[i])
    if i>0:
        for j in range(i):
            mij = np.dot(X[i], Y[j])/np.dot(Y[j], Y[j])
            print np.dot(X[i], Y[j]) + "/"+ np.dot(Y[j], Y[j])
            M[i].append(mij)
            yi -= mij*Y[j]
    Y.append(yi)

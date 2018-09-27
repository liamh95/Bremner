from fractions import Fraction



class QMatrix:

    def __init__(self, matrix):
        good = True
        n = len(matrix[0])
        for row in matrix:
            if len(row) is not n:
                good = False
                break
        assert good
        self.shape = (len(matrix), n)
        self.matrix = [[Fraction(entry) for entry in row] for row in matrix]

    def __getitem__(self, index):
        return self.matrix[index]

    def __add__(self, other):
        assert self.shape == other.shape

        mat = [[self.matrix[i][j] + other.matrix[i][j] for j in range(self.shape[1])] for i in range(self.shape[0])]

        return QMatrix(mat)


    #scalar multiplication
    def __rmul__(self, other):

        return QMatrix([[other*entry for entry in row] for row in self.matrix])

    def __sub__(self, other):
        return self + ((-1)*other)

    def __mul__(self, other):
        assert self.shape[1] == other.shape[0]
        mat = []
        for i in range(self.shape[0]):
            row = []
            for j in range(other.shape[1]):
                entry = 0
                for k in range(self.shape[1]):
                    entry += self.matrix[i][k]*other.matrix[k][j]
                row.append(entry)
            mat.append(row)
        return QMatrix(mat)

    def __setitem__(self, key, item):
        self.matrix[key] = item


    def __str__(self):
        ret = '['

        #find longest string in each column
        col_lengths = []
        for j in range(self.shape[1]):
            longest = 1
            for i in range(self.shape[0]):
                length = len(str(self.matrix[i][j]))
                if length > longest:
                    longest = length
            col_lengths.append(longest)

        for i in range(self.shape[0]):
            if i is not 0:
                ret += ' '
            ret += '['
            for j in range(self.shape[1]):
                ret += ' '*( col_lengths[j] - len(str(self.matrix[i][j]))  ) + str(self.matrix[i][j])
                if j is not (self.shape[1]-1):
                    ret += ', '

            ret += ']'
            if i is not(self.shape[0]-1):
                ret += '\n'
        ret += ']'
        return ret

    @staticmethod
    def dot(u, v):
        assert len(u)==len(v)
        return sum(u[i]*v[i] for i in range(len(u)))
    
    def gso(self):
        pass
        ymat = []
        mmat = []
        for i in range(self.shape[0]):
            xi = QMatrix([self.matrix[i]])
            yi = xi
            if i>0:
                mi = []
                for j in range(i):
                    mij = self.dot(xi.matrix[0], ymat[j]) / self.dot(ymat[j], ymat[j])
                    mi.append(mij)
                    yi = yi - (mij*QMatrix([ymat[j]])) 
                mi.append(1)
                mi += [0]*(self.shape[1]-len(mi))
                mmat.append(mi)

            else:
                mmat.append([1] + ([0]*(self.shape[1]-1)))

            ymat.append(yi.matrix[0])
        #print(mmat)
        #print(ymat)
        return (QMatrix(mmat), QMatrix(ymat))

    
    def reduce(Y, M, k, l):
        if abs(M[k][l])>Fraction(1,2):
            yk = QMatrix([Y.matrix[k]])
            yk = yk - (round(M[k][l]) * QMatrix([Y.matrix[l]]))
            Y.matrix[k] = yk.matrix[0]

            for j in range(l):
                M[k][j] = M[k][j] - round(M[k][l])*M[l][j]
            M[k][l] = M[k][l] - round(M[k][l])

    def exchange(Y, Ygso, gamma, M, k):
        z = Y[k-1]
        v = M[k][k-1]
        Y[k-1] = Y[k]
        Y[k] = z

        delta = gamma[k] + (v*v*gamma[k-1])
        M[k][k-1] = v*gamma[k-1] / delta
        gamma[k] = gamma[k] * gamma[k-1] / delta
        gamma[k-1] = delta

        for j in range(k-1):
            t = M[k-1][j]
            M[k-1][j] = M[k][j]
            M[k][j] = t

        for i in range(k+1, Y.shape[1]):
            #print('({}, {}'.format(i, k))
            xi = M[i][k]
            M[i][k] = M[i][k-1]-v*M[i][k]
            M[i][k-1] = M[k][k-1]*M[i][k]+xi

    def LLL(self, alpha):
        Y = QMatrix(self.matrix)
        (M, Ystar) = Y.gso()
        gamma = [QMatrix.dot(yi, yi) for yi in Ystar.matrix]

        k = 1
        while k<= (self.shape[1]-1):
            #print(k)
            QMatrix.reduce(Y, M, k, k-1)
            if gamma[k] >= (alpha - M[k][k-1]**2)*gamma[k-1]:
                for l in range(k-2, -1, -1):
                    QMatrix.reduce(Y, M, k, l)
                k = k+1
            else:
                QMatrix.exchange(Y, Ystar, gamma, M, k)
                if k>1:
                    k = k-1

        return Y



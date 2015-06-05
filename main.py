from math import sqrt
from itertools import count, islice, zip_longest, repeat
import sys

from pprint import pprint

def isPrime(n):
    if n < 2: return False
    return all(n%i for i in islice(count(2), int(sqrt(n)-1)))

def xgcd(a,b):
    """Extended euclidian algorithm"""
    prevx, x = 1, 0
    prevy, y = 0, 1
    while b:
        q = a//b
        x, prevx = prevx - q*x, x
        y, prevy = prevy - q*y, y
        a, b = b, a % b
    return a, prevx, prevy

class FieldElement:
    def get_inverse(self):
        raise NotImplementedError()

    def get_negative(self):
        raise NotImplementedError()

    def __add__(self, other):
        if type(self) is type(other):
            raise NotImplementedError()
        return NotImplemented

    def __mul__(self, other):
        """Multiplies with another element or an
        integer by summing up.
        """
        if type(self) is type(other):
            raise NotImplementedError()
        elif isinstance(other, int):
            # self + self + ... + self
            return sum(repeat(self, other))
        return NotImplemented

    def __sub__(self, other):
        if type(self) is type(other):
            return self + other.get_negative()
        return NotImplemented

    def __truediv__(self, other):
        if type(self) is type(other):
            return self * other.get_inverse()
        return NotImplemented

    def __radd__(self, other):
        """Fixes builtin sum"""
        if other == 0:
            return self
        return NotImplemented



class ResidueFieldElement(FieldElement):
    """Represents an element of the residue field with a given prime number"""

    def __init__(self, value, prime):
        if not isPrime(prime):
            raise ValueError("Expected prime, got {}".format(prime))
        self.prime = prime
        self.value = value % prime

    def get_inverse(self):
        if self.value == 0:
            return None
        else:
            gcd, a, b = xgcd(self.value, self.prime)
            return ResidueFieldElement(a, self.prime)

    def get_negative(self):
        return ResidueFieldElement(-self.value, self.prime)

    def __add__(self, other):
        if isinstance(other, ResidueFieldElement):
            if other.prime != self.prime:
                raise ArithmeticError("Tried to add elements of different residue fields.")
            return ResidueFieldElement(self.value + other.value, self.prime)
        return NotImplemented

    def __mul__(self, other):
        if isinstance(other, ResidueFieldElement):
            if other.prime != self.prime:
                raise ArithmeticError("Tried to multiply elements of different residue fields.")
            return ResidueFieldElement(self.value * other.value, self.prime)
        return NotImplemented


    def __eq__(self, other):
        """check for equality"""
        if isinstance(other, ResidueFieldElement):
            if self.prime == other.prime:
                return self.value == other.value
            else:
                return False
        return NotImplemented

    def __ne__(self, other):
        """check for non-equality"""
        result = self.__eq__(other)
        if result is NotImplemented:
            return result
        return not result

    def __lt__(self,other):
        """less than"""
        if isinstance(other, ResidueFieldElement) and self.prime == other.prime:
                return self.value < other.value
        return NotImplemented

    def __le__(self,other):
        """less or equal"""
        if isinstance(other, ResidueFieldElement) and self.prime == other.prime:
                return self.value <= other.value
        return NotImplemented

    def __gt__(self,other):
        """greater than"""
        if isinstance(other, ResidueFieldElement) and self.prime == other.prime:
                return self.value > other.value
        return NotImplemented

    def __ge__(self,other):
        """greater or equal"""
        if isinstance(other, ResidueFieldElement) and self.prime == other.prime:
                return self.value >= other.value
        return NotImplemented




class Transversal:
    def __init__(self, t, prime):
        """Uses the list t of numbers to represent all residue classes"""
        self.prime = prime
        self.t_dict = dict()

        if not isPrime(prime):
            raise ValueError("Expected prime, got {}".format(prime))

        needed = list(range(prime))
        for i in t:
            residue = i % prime
            if residue in self.t_dict:
               raise MalformedTransversalError("Doubled representative.")

            self.t_dict[residue] = i
            needed.remove(residue)

        if len(needed) != 0:
            raise MalformedTransversalError("Representative missing.")

    def get_repr(self, n):
        """Returns the representative of number in this transversal"""
        return self.t_dict[n.value]

class MalformedTransversalError(Exception):
    pass

class StdTransversal(Transversal):
    """The transversal 0,...,prime - 1
    For example for prime = 5: 0, 1, 2, 3, 4
    """
    def __init__(self, prime):
        self.prime = prime

    def get_repr(self, n):
        return n.value % self.prime

class StdNegTransversal(Transversal):
    """The transversal -prime / 2,...,0,...,prime / 2
    For example for prime = 5: -2, -1, 0, 1, 2
    """
    def __init__(self, prime):
        self.prime = prime

    def get_repr(self, n):
        v = n.value % prime
        if v > prime // 2:
            return v - prime
        else:
            return v


class MalformedMatrixError(Exception):
    pass

class Matrix:
    def __init__(self, m):
        self.m = m
        self.width, self.height = self.__process_new_matrix()

    def create_similar_matrix(self, new_m):
        """Subclasses can use this to let matrix operations
        create matrices of their own type.
        """
        return Matrix(new_m)

    def process_new_element(self, ele):
        """Subclasses can use this to process new entries
        before they are saved in the matrix.
        """
        return ele

    def transpose(self):
        # zip the rows
        transposed = zip(*self.m)
        transposed = [list(c) for c in transposed]
        return self.create_similar_matrix(transposed)

    def __process_new_matrix(self):
        width = None
        for i in range(len(self.m)):
            l = len(self.m[i])
            if width == None:
                width = l
            elif width != l:
                raise MalformedMatrixError("The matrix has multiple row lengths.")

            for j in range(l):
                self.m[i][j] = self.process_new_element(self.m[i][j])

        return width, len(self.m)

    def __mul__(self, other):
        if isinstance(other, int):
            for r in range(self.height):
                for c in range(self.width):
                    self.m[r][c] = other * self.m[r][c]
        elif isinstance(other, Matrix):
            if self.width != other.height:
                error_msg = "Tried to multiplicate matrices with unfitting dimensions: {} and {}"
                raise ArithmeticError(error_msg.format(self.width, other.height))

            tr_other = other.transpose()

            # the columns of other are the rows of tr_other
            product = [[sum(ele_a * ele_b for ele_a, ele_b in zip(row_a, col_b)) 
                        for col_b in tr_other.m] 
                       for row_a in self.m]
            return self.create_similar_matrix(product)
            
        return NotImplemented

    def __add__(self, other):
        if isinstance(other, Matrix):
            if self.width != other.width or self.height != other.height:
                error_msg = "Tried to add matrices with unfitting dimensions: {}x{} and {}x{}"
                raise ArithmeticError(error_msg.format(self.height, self.width, other.height, other.width))

            summed = [[ele_a + ele_b for ele_a, ele_b in zip(row_a, row_b)]
                        for row_a, row_b in zip(other.m, self.m)]

            return self.create_similar_matrix(summed)
        return NotImplemented

    def __sub__(self, other):
        if isinstance(other, Matrix):
            if self.width != other.width or self.height != other.height:
                error_msg = "Tried to subtract matrices with unfitting dimensions: {}x{} and {}x{}"
                raise ArithmeticError(error_msg.format(self.height, self.width, other.height, other.width))

            diff = [[ele_a - ele_b for ele_a, ele_b in zip(row_a, row_b)]
                        for row_a, row_b in zip(other.m, self.m)]

            return self.create_similar_matrix(diff)
        return NotImplemented

    def __or__(self, other):
        """Concatenates matrices to each other like the following
        /1 1 1\ | /2 2\   /1 1 1 2 2\ 
        |1 1 1| | |2 2| = |1 1 1 2 2|
        \1 1 1/ | \2 2/   \1 1 1 2 2/
        """
        if type(self) is type(other):
            if self.height != other.height:
                error_msg = "Tried to concatenate matrices with different heights: {} and {}"
                raise ArithmeticError(error_msg.format(self.height, other.height))
            new_m = [row_a + row_b for row_a, row_b in zip(self.m, other.m)]
            return self.create_similar_matrix(new_m)
        return NotImplemented

    def __radd__(self, other):
        """Fixes builtin sum"""
        if other == 0:
            return self
        return NotImplemented

    def __getitem__(self, pos):
        """Get the element at position (pos_h, pos_w) in the matrix or
        the row/column as new matrix, if one of them is -1
        """
        if len(pos) == 2:
            (pos_h, pos_w) = pos
            if pos_h >= self.height or pos_w >= self.width or max(pos_h, pos_w) < 0:
                raise IndexError("pos is outside of matrix bounds")
            elif pos_h >= 0 and pos_w >= 0:
                return self.m[pos_h][pos_w]
            elif pos_w < 0:
                # return row
                row = [self.m[pos_h]]
                return self.create_similar_matrix(row)
            else:
                # return column
                col = [[row[pos_w]] for row in self.m]
                return self.create_similar_matrix(col)

        return NotImplemented

    def __setitem__(self, pos, value):
        """Set the element at position (pos_h, pos_w) in the matrix or
        replace the specified row/column.
        """
        if len(pos) == 2:
            (pos_h, pos_w) = pos
            if pos_h >= self.height or pos_w >= self.width or max(pos_h, pos_w) < 0:
                raise IndexError("pos is outside of matrix bounds")
            elif pos_h >= 0 and pos_w >= 0:
                self.m[pos_h][pos_w] = self.process_new_element(value)
            elif pos_w < 0:
                # return row
                row = [self.m[pos_h]]
                return self.create_similar_matrix(row)
            else:
                # return column
                col = [[row[pos_w]] for row in self.m]
                return self.create_similar_matrix(col)

        return NotImplemented

    def __str__(self):
        # convert all elements to strings
        str_matrix = [[str(ele) for ele in row] for row in self.m]
        # transpose string matrix
        tr_str_matrix = list(zip(*str_matrix))
        # find max column width for each column
        max_width = [max(len(ele) for ele in col) for col in tr_str_matrix]
        # join rows to strings, using the individual column widths
        lines = ["[" + ", ".join(ele.rjust(col_width) for ele, col_width in zip(row, max_width)) + "]"
                 for row in str_matrix]

        return "[" + ",\n ".join(lines) + "]"

class ResFieldMatrix(Matrix):
    """A matrix over a residue field. All cells are 
    instances of ResidueFieldElement.
    """
    def __init__(self, m, prime, transversal=None):
        """Creates a matrix over a residue field. m can be
        a normal int matrix, which will be converted internally.
        """
        self.prime = prime
        if transversal == None:
            transversal = StdTransversal(prime)
        elif transversal.prime != prime:
            raise ValueError("Transversal does not match the specified residue field.")
        self.transversal = transversal

        super(ResFieldMatrix, self).__init__(m)

                
    def process_new_element(self, ele):
        """Convert normal integers to ResidueFieldElement"""
        if isinstance(ele, int):
            return ResidueFieldElement(ele, self.prime)
        elif isinstance(ele, ResidueFieldElement):
            return ele

        raise ValueError("Unsupported type for matrix element.")

    def create_similar_matrix(self, new_m):
        """Creates a matrix with the same prime and transversal"""
        return ResFieldMatrix(new_m, self.prime, transversal=self.transversal)

    def __str__(self):
        # convert all elements to strings
        str_matrix = [[str(self.transversal.get_repr(ele)) for ele in row] for row in self.m]
        # transpose string matrix
        tr_str_matrix = list(zip(*str_matrix))
        # find max column width for each column
        max_width = [max(len(ele) for ele in col) for col in tr_str_matrix]
        # join rows to strings, using the individual column widths
        lines = ["[" + ", ".join(ele.rjust(col_width) for ele, col_width in zip(row, max_width)) + "]"
                 for row in str_matrix]

        return "[" + ",\n ".join(lines) + "]"


def gauss(A):
    n = len(A)

    for i in range(0, n):
        # Search for maximum in this column
        maxEl = abs(A[i,i])
        maxRow = i
        for k in range(i+1, n):
            if abs(A[k,i]) > maxEl:
                maxEl = abs(A[k,i])
                maxRow = k

        # Swap maximum row with current row (column by column)
        for k in range(i, n+1):
            tmp = A[maxRow,k]
            A[maxRow,k] = A[i,k]
            A[i,k] = tmp

        # Make all rows below this one 0 in current column
        for k in range(i+1, n):
            c = -A[k,i]/A[i,i]
            for j in range(i, n+1):
                if i == j:
                    A[k,j] = 0
                else:
                    A[k,j] += c * A[i,j]

    # Solve equation Ax=b for an upper triangular matrix A
    x = [0 for i in range(n)]
    for i in range(n-1, -1, -1):
        x[i] = A[i][n]/A[i][i]
        for k in range(i-1, -1, -1):
            A[k][n] -= A[k][i] * x[i]
    return x

def input_matrix():
    """Creates a python matrix from command line input"""
    print("Type matrix. Press CTRL-D to end input:")
    text = sys.stdin.read()
    lines = text.split("\n")
    max_width = []
    elements = []
    for i, l in zip(range(len(lines)), lines):
        l = l.split(" ")
        l = [e for e in l if e != ""]
        if len(l) != 0:
            elements.append(l)

            for j, e, old_max in zip_longest(range(max(len(elements[-1]), len(max_width))), 
                                             elements[-1], max_width, fillvalue=0):
                if e != 0:
                    if j >= len(max_width):
                        max_width.append(len(e))
                    else:
                        max_width[j] = max(len(e), old_max)

        
    for i in range(len(elements)):
        elements[i] = [e.rjust(col_width) for e, col_width in zip(elements[i], max_width)]
        elements[i] = "[" + ", ".join(elements[i]) + "]"
    

    return "[" + ",\n ".join(elements) + "]"


#   
#   Examples
#

prime = 7

trans_pr = StdTransversal(prime)
neg_trans_pr = StdNegTransversal(prime)


value = 3
a = ResidueFieldElement(value, prime)
a_i = a.get_inverse()
msg = "The inverse of {} in F_{} is: {}"
print(msg.format(value, prime, trans_pr.get_repr(a_i)))

A = [[ 1, -9,  3,  5],
     [ 7, -3, -6,  8],
     [ 5, -2, -4,  6],
     [-5,  1,  2, -8]]

C = [[ 1, 2, -3,  5],
     [-1, 0,  5, -7],
     [ 0, 2,  4, -4]]

D = [[ 8, 7, 4, -3],
     [-4, 2, 3,  1],
     [ 7, 2, 0, -2]]

E = [[ 1, -1],
     [-1,  2],
     [ 3, -1]]

F = [[2, 7],
     [5, 4],
     [1, 8]]

G = [[0, 1, 1, 1],
     [1, 0, 1, 1],
     [1, 1, 0, 1]]



A = ResFieldMatrix(A, prime, transversal=neg_trans_pr)
C = ResFieldMatrix(C, prime, transversal=neg_trans_pr)
D = ResFieldMatrix(D, prime, transversal=neg_trans_pr)
E = ResFieldMatrix(E, prime, transversal=neg_trans_pr)
F = ResFieldMatrix(F, prime, transversal=neg_trans_pr)

G = ResFieldMatrix(G, prime, transversal=neg_trans_pr)

# select the cell at A_(3,4) (indices start at 0)
print("A_(3,4) = {}".format(A[2,3].value))
print()
# select the 3rd column, A_(-1,3), ie. at index 2. Analog for rows. Result is a new matrix
print("A_(-1,3) = \n{}".format(A[-1,2]))
print()

# append F to E and multiply the result with A
#R = (E | F) * A
R = (E | F) * A
print("(E | F) * A = \n{}".format(R))
print()

print(input_matrix())

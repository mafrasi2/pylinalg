from math import sqrt
from fractions import Fraction
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

class Field:
    def from_representant(self, repr):
        raise NotImplementedError()

    def element_to_str(self, ele):
        raise NotImplementedError()

    def get_zero(self):
        raise NotImplementedError()

    def get_one(self):
        raise NotImplementedError()

    def get_inverse(self, ele):
        raise NotImplementedError()

    def get_negative(self, ele):
        raise NotImplementedError()

    def add(self, a, b):
        raise NotImplementedError()

    def mul(self, a, b):
        raise NotImplementedError()

    def eq(self, a, b):
        """Checks a and b for equality"""
        raise NotImplementedError()

    def sub(self, a, b):
        return a + self.get_negative(b)

    def div(self, a, b):
        return a * self.get_inverse(b)


class FieldElement:
    def __init__(self, value, field):
        self.value = value
        self.field = field

    def __repr__(self):
        return self.field.element_to_str(self)

    def __str__(self):
        return self.__repr__()

    def __add__(self, other):
        if isinstance(other, FieldElement) and other.field == self.field:
            return self.field.add(self, other)
        else:
            return self.field.add(self, self.field.from_representant(other))

    def __mul__(self, other):
        if isinstance(other, FieldElement) and other.field == self.field:
            return self.field.mul(self, other)
        else:
            return self.field.mul(self, self.field.from_representant(other))

    def __sub__(self, other):
        if isinstance(other, FieldElement) and other.field == self.field:
            return self.field.sub(self, other)
        else:
            return self.field.sub(self, self.field.from_representant(other))

    def __truediv__(self, other):
        if isinstance(other, FieldElement) and other.field == self.field:
            return self.field.div(self, other)
        else:
            return self.field.div(self, self.field.from_representant(other))

    def __neg__(self):
        return self.field.get_negative(self)

    def __eq__(self, other):
        """Check for equality"""
        if isinstance(other, FieldElement) and other.field == self.field:
            return self.field.eq(self, other)
        else:
            return self.field.eq(self, self.field.from_representant(other))

    def __rmul__(self, other):
        return self.field.mul(self.field.from_representant(other), self)

    def __radd__(self, other):
        if isinstance(other, int) and other == 0:
            # Fixes builtin sum
            return self
        else:
            return self.field.add(self.field.from_representant(other), self)

    def __rsub__(self, other):
        return self.field.sub(self.field.from_representant(other), self)

    def __ne__(self, other):
        """Check for non-equality"""
        result = self.__eq__(other)
        if result is NotImplemented:
            return result
        return not result

class ResidueField(Field):
    def __init__(self, prime, transversal=None):
        if not isPrime(prime):
            raise ValueError("Expected prime, got {}".format(prime))
        self.prime = prime

        if transversal == None:
            transversal = StdTransversal(prime)
        elif transversal.prime != prime:
            raise ValueError("Transversal does not match the specified residue field.")
        self.transversal = transversal

    def from_representant(self, rep):
        if isinstance(rep, int):
            return FieldElement(rep, self)
        elif isinstance(rep, FieldElement) and rep.field == self:
            return rep

    def element_to_str(self, ele):
        return str(self.transversal.get_repr(ele))

    def get_zero(self):
        return FieldElement(0, self)

    def get_one(self):
        return FieldElement(1, self)

    def get_inverse(self, ele):
        if ele.value % self.prime == 0:
            raise ArithmeticError("Inverse of 0 is not defined.")
        else:
            gcd, a, b = xgcd(ele.value % self.prime, self.prime)
            return FieldElement(a, self)

    def get_negative(self, ele):
        return FieldElement((-ele.value) % self.prime, self)

    def add(self, a, b):
        return FieldElement((b.value + a.value) % self.prime, self)

    def mul(self, a, b):
        return FieldElement((b.value * a.value) % self.prime, self)

    def eq(self, a, b):
        return a.field == b.field and a.value % self.prime == b.value % self.prime

class RationalNumbersField(Field):
    def __init__(self, output_as_fraction=True):
        self.output_as_fraction = output_as_fraction

    def from_representant(self, rep):
        """Please don't use floats to prevent any misbehaviour"""
        if isinstance(rep, tuple):
            return FieldElement(Fraction(*rep), self)
        elif isinstance(rep, int):
            return FieldElement(Fraction(rep), self)
        elif isinstance(rep, Fraction):
            return FieldElement(rep, self)
        elif isinstance(rep, FieldElement) and rep.field == self:
            return rep
        return NotImplemented

    def element_to_str(self, ele):
        if self.output_as_fraction:
            return str(ele.value)
        else:
            return str(float(ele.value))

    def get_zero(self):
        return FieldElement(Fraction(0,1), self)

    def get_one(self):
        return FieldElement(Fraction(1,1), self)

    def get_inverse(self, ele):
        return FieldElement(1 / ele.value, self)

    def get_negative(self, ele):
        return FieldElement(-ele.value, self)

    def add(self, a, b):
        return FieldElement(b.value + a.value, self)

    def mul(self, a, b):
        return FieldElement(b.value * a.value, self)

    def eq(self, a, b):
        return a.value == b.value

class MalformedTransversalError(Exception):
    pass

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
        v = n.value % self.prime
        if v > self.prime // 2:
            return v - self.prime
        else:
            return v


class MalformedMatrixError(Exception):
    pass

class Matrix:
    def __init__(self, m, field):
        self.m = [[cell for cell in row] for row in m]
        self.field = field
        self.width, self.height = self.__process_new_matrix()

    def transpose(self):
        # zip the rows
        transposed = zip(*self.m)
        transposed = [list(c) for c in transposed]
        return Matrix(transposed, self.field)

    def __process_new_matrix(self):
        """Checks consistent width and converts all representatives to field elements"""

        width = None
        for i in range(len(self.m)):
            l = len(self.m[i])
            if width == None:
                width = l
            elif width != l:
                raise MalformedMatrixError("The matrix has multiple row lengths.")

            for j in range(l):
                self.m[i][j] = self.field.from_representant(self.m[i][j])

        return width, len(self.m)

    def copy(self):
        """Returns a new matrix with the same elements"""
        return Matrix([list(row) for row in self.m], self.field)

    def __mul__(self, other):
        """Scalar and matrix multiplication"""
        if isinstance(other, Matrix):
            if self.width != other.height:
                error_msg = "Tried to multiplicate matrices with unfitting dimensions: {} and {}"
                raise ArithmeticError(error_msg.format(self.width, other.height))

            tr_other = other.transpose()

            # the columns of other are the rows of tr_other
            product = [[sum(ele_a * ele_b for ele_a, ele_b in zip(row_a, col_b)) 
                        for col_b in tr_other.m] 
                       for row_a in self.m]
            return Matrix(product, self.field)
        elif not isinstance(other, FieldElement):
            other = self.field.from_representant(other)
        elif not other.field == self.field:
            return NotImplemented

        m = [[cell * other for cell in row] for row in self.m]
        return Matrix(m, self.field)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __pow__(self, other):
        if isinstance(other, int):
            if other == 0:
                return identity_matrix(self.field, self.height, self.width)
            elif other > 0:
                factor = self
            elif other < 0:
                factor = self.get_inverse()
                other = -other

            result = factor
            for _ in range(other-1):
                result *= factor

            return result
        return NotImplemented

    def __add__(self, other):
        if isinstance(other, Matrix):
            if self.width != other.width or self.height != other.height:
                error_msg = "Tried to add matrices with unfitting dimensions: {}x{} and {}x{}"
                raise ArithmeticError(error_msg.format(self.height, self.width, other.height, other.width))

            summed = [[ele_a + ele_b for ele_a, ele_b in zip(row_a, row_b)]
                        for row_a, row_b in zip(other.m, self.m)]

            return Matrix(summed, self.field)
        return NotImplemented

    def __sub__(self, other):
        if isinstance(other, Matrix):
            if self.width != other.width or self.height != other.height:
                error_msg = "Tried to subtract matrices with unfitting dimensions: {}x{} and {}x{}"
                raise ArithmeticError(error_msg.format(self.height, self.width, other.height, other.width))

            diff = [[ele_a - ele_b for ele_a, ele_b in zip(row_a, row_b)]
                        for row_a, row_b in zip(other.m, self.m)]

            return Matrix(diff, self.field)
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
            return Matrix(new_m, self.field)
        return NotImplemented

    def __rmul__(self, other):
        return self.__mul__(other)

    def __radd__(self, other):
        """Fixes builtin sum"""
        if other == 0:
            return self
        return NotImplemented

    def __eq__(self, other):
        if isinstance(other, Matrix):
            return (self.field == other.field and
                   self.height == other.height and
                   self.width == other.width and
                   all(a == b for row_a, row_b in zip(self.m, other.m) 
                              for a, b in zip(row_a, row_b) )
                   )
        return NotImplemented

    def __getitem__(self, pos):
        """Get the element at position (pos_h, pos_w) in the matrix or
        the row/column as new matrix, if one of them is -1
        """
        if isinstance(pos, tuple) and len(pos) == 2:
            (pos_h, pos_w) = pos
            if pos_h >= self.height or pos_w >= self.width or max(pos_h, pos_w) < 0:
                msg = "position ({},{}) is outside of matrix bounds"
                raise IndexError(msg.format(pos_h, pos_w))
            elif pos_h >= 0 and pos_w >= 0:
                return self.m[pos_h][pos_w]
            elif pos_w < 0:
                # return row
                row = [self.m[pos_h]]
                return Matrix(row, self.field)
            else:
                # return column
                col = [[row[pos_w]] for row in self.m]
                return Matrix(col, self.field)
        return NotImplemented

    def __setitem__(self, pos, value):
        """Replace the element at position (pos_h, pos_w) in the matrix or
        replace the specified row/column.
        """
        if len(pos) == 2:
            (pos_h, pos_w) = pos
            if pos_h >= self.height or pos_w >= self.width or max(pos_h, pos_w) < 0:
                raise IndexError("pos is outside of matrix bounds")
            elif pos_h >= 0 and pos_w >= 0:
                self.m[pos_h][pos_w] = self.field.from_representant(value)
            elif pos_w < 0:
                # set row
                if value.width != self.width:
                    raise ValueError("The row doesn't have the same width as the matrix.")
                self.m[pos_h] = [self.field.from_representant(ele) for ele in value.m[0]]
            else:
                # set column
                if value.height != self.height:
                    raise ValueError("The column doesn't have the same height as the matrix.")

                for i in range(self.height):
                    self.m[i][pos_w] = self.field.from_representant(value.m[i][0]) 

        return NotImplemented

    def __str__(self):
        # convert all elements to strings
        str_matrix = [[str(ele) for ele in row] for row in self.m]
        # transpose string matrix
        tr_str_matrix = zip(*str_matrix)
        # find max column width for each column
        max_width = [max(len(ele) for ele in col) for col in tr_str_matrix]
        # join rows to strings, using the individual column widths
        lines = ["[" + ", ".join(ele.rjust(col_width) for ele, col_width in zip(row, max_width)) + "]"
                 for row in str_matrix]

        return "[" + ",\n ".join(lines) + "]"

    def trace(self):
        """Get the product of the diagonal entries"""
        if self.height != self.width:
            msg = "tr is undefined for a {}x{} matrix"
            raise ArithmeticError(msg.format(self.height, self.width))

        tr = self.field.get_one()
        for i in range(self.width):
            tr *= self[i,i]
        return tr

    def to_upper_triangular_matrix(self, normalize=False):
        """Transforms into an upper triangular matrix"""
        zero = self.field.get_zero()
        n = min(self.height, self.width)

        swaps = 0
        for i in range(n):
            # Search for first non-zero in this column
            max_row = i
            for max_row in range(i, n):
                if self[max_row, i] != zero:
                    break
            if max_row != i:
                # swap row i and max_row
                self[i,-1], self[max_row,-1] = self[max_row,-1], self[i,-1]
                swaps += 1

            if self[i,i] != zero:
                # Make all rows below this one 0 in current column
                for k in range(i + 1, self.height):
                    c = -self[k,i]/self[i,i]
                    self[k,i] = 0
                    for j in range(i + 1, self.width):
                        self[k,j] += c * self[i,j]

                if normalize:
                    inv = self.field.get_inverse(self[i,i])
                    self[i,i] = 1
                    for k in range(i + 1, self.width):
                        self[i,k] *= inv

        return swaps

    def to_diagonal_matrix(self, normalize=False):
        """Transforms into a diagonal matrix"""
        self.to_upper_triangular_matrix(normalize=normalize)

        zero = self.field.get_zero()
        n = min(self.height, self.width)

        for i in range(n-1,-1,-1):
            # find first non-zero cell in the current row
            first_nz = i
            for first_nz in range(i, self.width):
                if self[i, first_nz] != zero:
                    break

            if self[i, first_nz] != zero:
                for j in range(i-1,-1,-1):
                    c = -self[j,first_nz]/self[i,first_nz]
                    self[j,first_nz] = 0
                    for k in range(first_nz + 1, self.width):
                        self[j,k] += c * self[i,k]

    def solve(self, x):
        """Solves the equation Ab=x for b"""
        if self.width > self.height:
            raise ArithmeticError("Can't calculate solution for a matrix with width > height")

        tmp = self | x
        tmp.to_upper_triangular_matrix(normalize=False)

        n = tmp.width
        # Solve equation Ax=b for an upper triangular matrix tmp
        for i in range(n-2, -1, -1):
            tmp[i, n-1] = tmp[i,n-1]/tmp[i,i]
            for k in range(i-2, -1, -1):
                tmp[k,n-1] -= tmp[k,i] * tmp[i, n-1]
        return tmp[-1, n-1]

    def get_rank(self):
        """Calculates the rank of this matrix"""
        tmp = self.copy()
        tmp.to_upper_triangular_matrix()

        # Count non-zero rows
        zero = tmp.field.get_zero()
        count = 0
        for row in tmp.m:
            for e in row:
                if e != zero:
                    count += 1
                    break
        return count

    def get_determinant(self):
        """Calculates the determinant of this matrix"""
        if self.width != self.height:
            msg = "The determinant of a {}x{} matrix is undefined"
            raise ArithmeticError(msg.format(self.width, self.height))

        tmp = self.copy()
        swaps = tmp.to_upper_triangular_matrix()

        det = tmp.trace()

        if swaps % 2 != 0:
            return -det
        else:
            return det

    def get_inverse(self):
        """Calculates the inverse of this matrix, if possible"""
        if self.width != self.height:
            msg = "There is no inverse of a {}x{} matrix"
            raise ArithmeticError(msg.format(self.width, self.height))

        # create identity matrix
        id_matrix = identity_matrix(self.field, self.height, self.width)
        tmp = self | id_matrix
        # get the identity matrix to the left
        tmp.to_diagonal_matrix(normalize=True)
        # concatenate all columns on the right
        inv = tmp[-1,self.width]
        for i in range(self.width - 1):
            inv = inv | tmp[-1,i + self.width + 1]

        return inv

def identity_matrix(field, n, m):
    """Creates a identity matrix with height n and width m"""
    one = field.get_one()
    zero = field.get_zero()
    id_matrix = [[one if j == i else zero for j in range(m)] for i in range(n)]
    return Matrix(id_matrix, field)


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

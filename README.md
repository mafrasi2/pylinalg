#   Examples

Import the module:
```python
from pylinalg import *
```

## Fields
We will to use 𝔽<sub>5</sub> for our first examples.

Create the new residue field 𝔽<sub>5</sub> and its transversal {-2,...,2}:

```python
neg_trans_5 = StdNegTransversal(5)
res_field_5 = ResidueField(5, transversal=neg_trans_5)
```

Find the inverse of 3 in 𝔽<sub>5</sub>:

```python
value = 3
a = res_field_5.from_representant(value)
a_i = res_field_5.get_inverse(a)
msg = "The inverse of {} in F_{} is: {}"
print(msg.format(value, 5, a_i))
```

Swap the transversal of 𝔽<sub>5</sub> to {0,...,4}

```python
print(a)  # -2
print("Switching to standard transversal:")
trans_5 = StdTransversal(5)
res_field_5.transversal = trans_5
print(a)  #  3
```

## Basic matrix operations
Create some matrices:

```python
A = [[ 1, -9,  3,  5],
     [ 7, -3, -6,  8],
     [ 5, -2, -4,  6],
     [-5,  1,  2, -8]]

B = [[ 1, -1],
     [-1,  2],
     [ 3, -1]]

C = [[2, 7],
     [5, 4],
     [1, 8]]

A = Matrix(A, res_field_5)
B = Matrix(B, res_field_5)
C = Matrix(C, res_field_5)
```

Indices start at 0 in all functions calls!

Select the cell at A<sub>3,4</sub>:

```python
cell = A[2,3]
print("A_(3,4) = {}".format(cell))
```

Select the 3rd column, A<sub>-1,3</sub>, ie. at index 2. Rows work analog. The result is a new matrix.

```python
column = A[-1,2]
print("A_(-1,3) = \n{}".format(column))
```

Replace the 3rd row of A with ones.

```python
print(A)
print(" -> ")
row = Matrix([[1,1,1,1]], res_field_5)
A[2,-1] = row
print(A)
```

## Matrix arithmetic

Append C to B and multiply the result with A.

```python
R = (B | C) * A
print("(B | C) * A = \n{}".format(R))
```

## Solving systems of linear equations
Transform B | C into an upper triangular matrix. Use normalize=False, if you
don't want to normalize the diagonal.

```python
D = B | C
D.to_upper_triangular_matrix(normalize=True)
```

And into a diagonal matrix (the step above is not necessary to do this).
```python
D.to_diagonal_matrix(normalize=True)
```

Calculate the inverse of A
```python
print(A * A.get_inverse())
```

Solve the equation Eb=x for b in 𝔽<sub>2</sub>.

```python
trans_2 = StdTransversal(2)
res_field_2 = ResidueField(2, transversal=trans_2)

E = [[0, 1, 0],
     [1, 0, 1],
     [1, 1, 0]]

x = [[0],
     [0],
     [1]]
E = Matrix(E, res_field_2)
x = Matrix(x, res_field_2)

print("Extended coefficient matrix for Sol(E,x):")
print(E | x)
R = E.solve(x)
print("Eb=x for b = \n{}".format(R))
```

## Rank and determinant
Calculate the rank and determinant of A

```python
print("rk(A) = {}".format(E.get_rank()))
print("det(A) = {}".format(E.get_determinant()))
```

## Misc
Get a nicely formatted python matrix from console input.

```python
print(input_matrix())
```

The expected syntax is as following:

```
1 31 1 3
3 1 1 -1
12 1 41 0
```

which results in:

```python
[[ 1, 31,  1,  3],
 [ 3,  1,  1, -1],
 [12,  1, 41,  0]]
```

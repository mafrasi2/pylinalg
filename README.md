Examples
--------------

We will to use 𝔽<sub>5</sub> for our first examples.

     from pylinalg import *
     prime = 5

Create the new residue field 𝔽<sub>5</sub> and its transversal (-2,...,2):

     neg_trans_5 = StdNegTransversal(prime)
     res_field_5 = ResidueField(prime, transversal=neg_trans_5)

Find the inverse of 3 in 𝔽<sub>5</sub>:

     value = 3
     a = res_field_5.from_representant(value)
     a_i = res_field_5.get_inverse(a)
     msg = "The inverse of {} in F_{} is: {}"
     print(msg.format(value, prime, a_i))

Create some matrices:

     A = [[ 1, -9,  3,  5],
          [ 7, -3, -6,  8],
          [ 5, -2, -4,  6],
          [-5,  1,  2, -8]]

     E = [[ 1, -1],
          [-1,  2],
          [ 3, -1]]

     F = [[2, 7],
          [5, 4],
          [1, 8]]

     A = Matrix(A, res_field_5)
     E = Matrix(E, res_field_5)
     F = Matrix(F, res_field_5)


Indices start at 0 in all functions calls!

Select the cell at A<sub>3,4</sub>:

     cell = A[2,3]
     print("A_(3,4) = {}".format(cell))

Select the 3rd column, A<sub>-1,3</sub>, ie. at index 2. Rows work analog. The result is a new matrix.

     column = A[-1,2]
     print("A_(-1,3) = \n{}".format(column))

Replace the 3rd row with zeros.

     print(A)
     print(" -> ")
     row = Matrix([[0,0,0,0]], res_field_5)
     A[2,-1] = row
     print(A)

Switch the transversale for a field.

     print(F)
     print("Switching to standard transversal:")
     trans_5 = StdTransversal(prime)
     res_field_5.transversal = trans_5
     print(F)

Append F to E and multiply the result with A.

     R = (E | F) * A
     print("(E | F) * A = \n{}".format(R))


Solve the equation Gb=x for b in 𝔽<sub>2</sub>.

     prime = 2

     trans_2 = StdTransversal(prime)
     res_field_2 = ResidueField(prime, transversal=trans_2)

     G = [[0, 1, 0],
          [1, 0, 1],
          [1, 1, 0]]

     x = [[0],
          [0],
          [1]]
     G = Matrix(G, res_field_2)
     x = Matrix(x, res_field_2)

     print("Extended coefficient matrix for Sol(G,x):")
     print(G | x)
     R = G.solve(x)
     print("Gb=x for b = \n{}".format(R))


Get a nicely formatted python matrix from console input.

     print(input_matrix())

The expected syntax is as following:
    
     1 31 1 3
     3 1 1 -1
     12 1 41 0

which results in:

     [[ 1, 31,  1,  3],
      [ 3,  1,  1, -1],
      [12,  1, 41,  0]]
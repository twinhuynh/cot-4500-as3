#Question 1:
def euler_method(f, a, b, n, initial):
  h = (b - a) / n
  t = a
  y = initial
  for i in range(n):
    y = y + h * f(t, y)
    t = t + h
  return y


def func(t, y):
  return t - y**2


result = euler_method(func, 0, 2, 10, 1)
print(result)
print()


#Question 2:
def runge_kutta(func, t0, y0, h, n):
  t = t0
  y = y0

  for i in range(n):
    k1 = h * func(t, y)
    k2 = h * func(t + h / 2, y + k1 / 2)
    k3 = h * func(t + h / 2, y + k2 / 2)
    k4 = h * func(t + h, y + k3)

    y = y + (k1 + 2 * k2 + 2 * k3 + k4) / 6
    t = t + h

  print(y)

t0 = 0
y0 = 1
h = 0.2
n = 10

runge_kutta(func, t0, y0, h, n)

print()


#Question 3:
from fractions import Fraction

def print_matrix(matrix):
  for row in matrix:
    print(row)

def gaussian_elimination(A):
  n = len(A)

  # Divide the first row by its first element to make it 1
  divisor = A[0][0]
  for j in range(n + 1):
    A[0][j] = Fraction(A[0][j], divisor)

  # Subtract multiples of the first row from the last two rows
  for i in range(1, n):
    factor = A[i][0] / A[0][0]
    for j in range(n + 1):
      A[i][j] -= factor * A[0][j]

  # Divide the second row by its second element to make it 1
  divisor = A[1][1]
  for j in range(1, n + 1):
    A[1][j] = Fraction(A[1][j], divisor)

  # Subtract multiples of the second row from the third row
  factor = A[2][1] / A[1][1]
  for j in range(1, n + 1):
    A[2][j] -= factor * A[1][j]

  # Divide the third row by its third element to make it 1
  divisor = A[2][2]
  for j in range(2, n + 1):
    A[2][j] = Fraction(A[2][j], divisor)

  # Backward elimination
  for i in range(n - 1, 0, -1):
    for j in range(i - 1, -1, -1):
      factor = A[j][i] / A[i][i]
      for k in range(i, n + 1):
        A[j][k] -= factor * A[i][k]

  # Extract solutions
  solutions = [row[-1] for row in A]

  return solutions

A = [[2, -1, 1, 6], [1, 3, 1, 0], [-1, 5, 4, -3]]

solutions = gaussian_elimination(A)

# Print the solutions
print("[", end=" ")
for i, solution in enumerate(solutions):
  if i != len(solutions) - 1:
    print(solution, end=" ")
  else:
    print(solution, end=" ")
print("]")

print()

    
#Question 4:
A = [[1, 1, 0, 3],
     [2, 1, -1, 1],
     [3, -1, -1, 2],
     [-1, 2, 3, -1]]

# calculate the det of a matrix
def determinant(matrix):
    # if matrix is 2x2
    if len(matrix) == 2:
        return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]

    det = 0
    for i in range(len(matrix)):
        submatrix = [row[:i] + row[i+1:] for row in matrix[1:]]
        det += matrix[0][i] * ((-1) ** i) * determinant(submatrix)
    return det

L = [[0] * 4 for _ in range(4)]
for i in range(4):
    L[i][i] = 1

U = [[0] * 4 for _ in range(4)]

for i in range(4):
    # Upper triangular matrix
    for k in range(i, 4):
        U[i][k] = A[i][k] - sum(L[i][j] * U[j][k] for j in range(i))

    # Lower triangular matrix
    for j in range(i+1, 4):
        L[j][i] = (A[j][i] - sum(L[j][k] * U[k][i] for k in range(i))) / U[i][i]

det_A = 1
for i in range(4):
    det_A *= U[i][i]

def print_matrix(matrix):
    for row in matrix:
        print("[", end="")
        for element in row:
            print("{:4}".format(int(element)), end=" ")
        print("]")

# Print out determinant
print("{:.0f}".format(det_A))

# Print out L matrix
print()
print_matrix(L)

# Print out U matrix
print()
print_matrix(U)
print()


#Question 5:
A = [[9, 0, 5, 2, 1], [3, 9, 1, 2, 1], [0, 1, 7, 2, 3], [4, 2, 3, 12, 2],
          [3, 2, 4, 0, 8]]

def is_diagonally_dominant(A):
  for i in range(len(A)):
    diagonal_element = abs(A[i][i])
    off_diagonal_sum = sum(
        abs(A[i][j]) for j in range(len(A)) if j != i)

    if diagonal_element <= off_diagonal_sum:
      return False
  return True


if is_diagonally_dominant(A):
  print("True")
else:
  print("False")

print()


#Question 6:
A = [[2, 2, 1], [2, 3, 0], [1, 0, 2]]

# compute the determinant of a 2x2 matrix
def det2x2(matrix):
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]

# compute the determinant of a 3x3 matrix
def det3x3(matrix):
    return (matrix[0][0] * det2x2([[matrix[1][1], matrix[1][2]], [matrix[2][1], matrix[2][2]]])
            - matrix[0][1] * det2x2([[matrix[1][0], matrix[1][2]], [matrix[2][0], matrix[2][2]]])
            + matrix[0][2] * det2x2([[matrix[1][0], matrix[1][1]], [matrix[2][0], matrix[2][1]]]))

# coefficients of polynomial
a = 1
b = -(A[0][0] + A[1][1] + A[2][2])
c = det3x3(A)

# discriminant
disc = b ** 2 - 4 * a * c

# Check if non-negative
if disc >= 0:
    # Compute the eigenvalues using the quadratic formula
    eigenvalue1 = (-b + disc ** 0.5) / (2 * a)
    eigenvalue2 = (-b - disc ** 0.5) / (2 * a)

    # Check if both eigenvalues are positive
    if eigenvalue1 > 0 and eigenvalue2 > 0:
        print("True")
    else:
        print("False")
else:
    print("False")

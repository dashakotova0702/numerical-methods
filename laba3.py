import sys
import numpy as np
A = np.array([[2.31, 31.49, 1.52],
[4.21, 22.42, 3.85],
[3.49, 4.85, 28.72]])
B = np.array([40.95, 30.24, 42.81])
lengthA = len(A)
lengthB = len(B)
x = np.zeros(lengthA)
A_norm = np.array(A)
B_norm = np.array(B)
if np.linalg.det(A) == 0:
    sys.exit(1)
if lengthA != lengthB:
    sys.exit(1)
for i in range(0, lengthA, 1):
    for j in range(0, lengthA, 1):
        if j - i <= 0:
            continue
        a = A[i][i]
        b = A[j][i] 
        for k in range(0, lengthA, 1):
            A[j][k] = a * A[j][k] - b * A[i][k]
            B[j] = a * B[j] - b * B[i]
for i in range(lengthA-1, -1, -1):
    x[i] = B[i]
    for j in range(0, lengthA, 1):
        if j != i:
            x[i] = x[i] - A[i][j] * x[j]
    x[i] = x[i] / A[i][i]
F = np.zeros(lengthA)
for i in range(0, lengthA, 1):
    for j in range (0, lengthA, 1):
        F[i] = F[i] + A_norm[i][j] * x[j]
    F[i] = F[i] - B_norm[i]
norm = max(F)
print("x:", x)
print("F: ", F)
print("Norm: ", norm)
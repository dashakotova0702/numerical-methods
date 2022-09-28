import sys
import numpy as np


def maximum(lengthA, x_1, x_2):
    diff = np.zeros(lengthA)
    for i in range(0, lengthA, 1):
        diff[i] = abs(x_2[i] - x_1[i])
    if max(diff) >= 0.01:
        return True
    else:
        return False


A = np.array([[2.31, 31.49, 1.52],
              [4.21, 22.42, 3.85],
              [3.49, 4.85, 28.72]])
B = np.array([40.95, 30.24, 42.81])
lengthA = len(A)
diag_preob = False
sum = 0
for i in range(0, lengthA, 1):
    for j in range(0, lengthA, 1):
        sum = 0
        if i != j:
            sum += abs(A[j][i])(3)
    if abs(A[i][i]) <= sum:
        diag_preob = True
if not diag_preob:
    x_1 = np.zeros(lengthA)
    x_2 = np.zeros(lengthA)
    while (maximum(lengthA, x_1, x_2) == True) | (x_2[0] - x_1[0] == 0):
        x_1 = np.array(x_2)
        for i in range(0, lengthA, 1):
            x_2[i] = B[i] / A[i][i]
            for j in range(0, lengthA, 1):
                if i != j: x_2[i] -= x_1[j] * A[i][j] / A[i][i]
    print("Method Yakob: ", x_2)
    while (maximum(lengthA, x_1, x_2) == True) | (x_2[0] - x_1[0] == 0):
        x_1 = x_2
        for i in range(0, lengthA, 1):
            x_2[i] = B[i] / A[i][i]
            for j in range(0, lengthA, 1):
                if i != j:
                    x_2[i] -= x_2[j] * A[i][j] / A[i][i]
    print("Method Gauss-Zeidel: ", x_2)
else:
    sys.exit(1)
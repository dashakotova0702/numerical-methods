import numpy as np 
from numpy import linalg as LA 

	 
e = 0.01 
lyambda = 0 
lyambda_next = 1000 
length = int(input('Matrix size: '))
A = np.zeros((length, length))
for i in range(0, length, 1): 
	for j in range(0, length, 1): 
		A[i][j] = input('') 
u = np.zeros(length) 
u_0 = np.zeros(length) 
print('U: ') 
for i in range(0, length, 1): 
	u_0[i] = input('') 
u = u_0 
u_next = np.zeros(length) 
while abs(lyambda_next - lyambda) > e: 
	lyambda = lyambda_next 
	u_next = A.dot(u) 
	numerator = u_next.dot(u) 
	denominator = u.dot(u) 
	lyambda_next = numerator / denominator                                         
	u = u_next 
	u_next = 0 
print('Lyambda_max = ', lyambda_next) 
A = LA.inv(A) 
lyambda = 0 
lyambda_next = 1000 
u = u_0 
u_next = 0 
while abs(lyambda_next - lyambda) > e: 
	lyambda = lyambda_next 
	u_next = A.dot(u) 
	numerator = u_next.dot(u) 
	denominator = u.dot(u) 
	lyambda_next = numerator / denominator 
	u = u_next 
	u_next = 0 
print('Lyambda_min = ', lyambda_next) 

import numpy as np 
from numpy import linalg as LA 
import math 
	 
e = 0.001 
Apprxmtn_1 = np.zeros(2) 
Apprxmtn_2 = np.zeros(2) 
Apprxmtn_1_next = np.zeros(2)
Apprxmtn_2_next = np.zeros(2) 
discrepancy_vector = np.zeros(2) 
temp = np.zeros(2)
d2_v = np.zeros(2) 
W = np.zeros((2, 2)) 
print('Initial approximation 1: ')
for i in range(0, 2, 1): 
	Apprxmtn_1[i] = float(input('')) 17 print('Initial approximation 2: ')
for i in range(0, 2, 1): 
	Apprxmtn_2[i] = float(input('')) 
NIT = int(input('NIT: ')) 
print('No            d1             d2') 
d1 = 10.00 
d2 = 10.00 
i = 0 
while d1 >= e and d2 >= e: 
	i += 1 
	discrepancy_vector[0] = math.cos(0.40 * Apprxmtn_1[1] + pow(Apprxmtn_1[0], 2)) + pow(Apprxmtn_1[1], 2) + pow(Apprxmtn_1[0], 2) - 1.6 
	discrepancy_vector[1] = 1.5 * pow(Apprxmtn_1[0], 2) - pow(Apprxmtn_1[1], 2) / 0.36 - 1 
	W[0][0] = -2 * Apprxmtn_1[0] * math.sin(0.4 * Apprxmtn_1[1] + pow(Apprxmtn_1[0], 2)) + 2 * Apprxmtn_1[0] 
	W[1][0] = 3 * Apprxmtn_1[0] 
	W[0][1] = -0.4 * math.sin(0.4 * Apprxmtn_1[1] + pow(Apprxmtn_1[0], 2)) + 2 * Apprxmtn_1[1] 
	W[1][1] = -2 * Apprxmtn_1[1] / 0.36 
	W = LA.inv(W) 
	Apprxmtn_1_next = Apprxmtn_1 - W.dot(discrepancy_vector)
	d1 = max(abs(discrepancy_vector[0]), abs(discrepancy_vector[1]))
	for j in range(0, 2, 1): 
		if Apprxmtn_1_next[j] < 1: 
			d2_v[j] = Apprxmtn_1_next[j] - Apprxmtn_1[j]
		else: 
			d2_v[j] = (Apprxmtn_1_next[j] - Apprxmtn_1[j]) / Apprxmtn_1_next[j] 
	d2 = max(abs(d2_v[0]), abs(d2_v[1])) 
	print(i, '        ', d1, '          ', d2) 
	Apprxmtn_1[0] = Apprxmtn_1_next[0] 
	Apprxmtn_1[1] = Apprxmtn_1_next[1]
	if i > NIT: 
		print("IER = 2") 
		break 
print('x1: ', Apprxmtn_1[0], '        x2: ', Apprxmtn_1[1]) 
d1 = 10 
d2 = 10 
i = 0 
while d1 >= e and d2 >= e: 
	i += 1 
	discrepancy_vector[0] = math.cos(0.40 * Apprxmtn_2[1] + pow(Apprxmtn_2[0], 2)) + pow(Apprxmtn_2[1], 2) + pow(Apprxmtn_2[0], 2) - 1.6 
	discrepancy_vector[1] = 1.5 * pow(Apprxmtn_2[0], 2) - pow(Apprxmtn_2[1], 2) / 0.36 - 1 
	W[0][0] = -2 * Apprxmtn_2[0] * math.sin(0.4 * Apprxmtn_2[1] + pow(Apprxmtn_2[0], 2)) + 2 * Apprxmtn_2[0] 
	W[1][0] = 3 * Apprxmtn_2[0] 
	W[0][1] = -0.4 * math.sin(0.4 * Apprxmtn_2[1] + pow(Apprxmtn_2[0], 2)) + 2 * Apprxmtn_2[1] 
	W[1][1] = -2 * Apprxmtn_2[1] / 0.36 
	W = LA.inv(W) 
	Apprxmtn_2_next = Apprxmtn_2 - W.dot(discrepancy_vector)
	d1 = max(abs(discrepancy_vector[0]), abs(discrepancy_vector[1]))
	for j in range (0, 2, 1): 
		if Apprxmtn_2_next[j] < 1: 
			d2_v[j] = Apprxmtn_2_next[j] - Apprxmtn_2[j]
		else: 
			d2_v[j] = (Apprxmtn_2_next[j] - Apprxmtn_2[j]) / Apprxmtn_2_next[j] 
	d2 = max(abs(d2_v[0]), abs(d2_v[1])) 
	print(i, '        ', d1, '          ', d2) 
	Apprxmtn_2[0] = Apprxmtn_2_next[0] 
	Apprxmtn_2[1] = Apprxmtn_2_next[1] 
	if i > NIT: 
	print("IER = 2") 
	break 
print('x1: ', Apprxmtn_2[0], '        x2: ', Apprxmtn_2[1]) 

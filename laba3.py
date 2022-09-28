import math 
import matplotlib.pyplot as plt 
import numpy as np 
import sympy 
from sympy import symbols, diff, sin, exp 
from numpy import linalg as LA 
 	 
def init_u(number): 
	u_start = { 
		0: [0, -0.412], 
		1: [1, 0], 
		2: [1, 1, 1] 
	} 
	u = np.array(u_start[number]) 
	return u 
  
def vect_F(count, number, y, t):
	if count != 3: 
		vector_F = { 
			0: [-y[0] * y[1] + (math.sin(t)) / t, -pow(y[1], 2) + (3.5 * t) / (1 + pow(t, 2))], 
			1: [y[1] - (2.25 * y[0] + y[1]) * y[0], math.exp(y[0]) - (y[0] + 2.25 * y[1]) * y[0]] 
		} 
		F = np.array(vector_F[number])
	else: 
		F = np.array([(2.25 - 1.25) / 1.25 * y[1] * y[2], (1.25 + 2.25) / 2.25 * y[0] * y[2], (1.25 - 2.25) / 1.25 * y[0] * y[1]]) 
	return F 
  
def explicit_method(count, u, tau_max, T, t, eps, number): 
	y = u 
	while t < T: 
		F = vect_F(count, number, y, t)
		tau = np.zeros(count) 
		for i in range(0, count, 1): 
			tau[i] = eps/(abs(F[i])+eps/tau_max) 
		step_tau = min(tau) 
		y = y + F * step_tau 
		plt.plot(t, y[0], ',', color='b') 
		plt.plot(t, y[1], ',', color='r')
		if number == 2: 
			plt.plot(t, y[2], ',', color='y') 
		t = t + step_tau 
  
def Newton(count, x0, F, t_now): 
	u1, u2, u3, t_symb = symbols("u1, u2, u3, t_symb") 
	c = 0 
	if count != 3: 
		x_next = np.array([1000, 1000]) 
		while abs(x_next[0]-x0[0]) > 0.01 and abs(x_next[1]-x0[1]) > 0.01: 
			c = c + 1
			if c > 5:
				return x_next
			else: 
				x_next = x0 
				x0 = x_next 
				f0 = np.array([float(((F[0].subs(u1, x0[0])).subs(u2, x0[1])).subs(t_symb, t_now)), float(((F[1].subs(u1, x0[0])).subs(u2, x0[1])).subs(t_symb, t_now))]) 
				W = np.array([[float(((diff(F[0], u1).subs(u1, x0[0])).subs(u2, x0[1])).subs(t_symb, t_now)), float(((diff(F[0], u2).subs(u1, x0[0])).subs(u2, x0[1])).subs(t_symb, t_now))], [float(((diff(F[1], u1).subs(u1, x0[0])).subs(u2, x0[1])).subs(t_symb, t_now)), float(((diff(F[1], u2).subs(u1, x0[0])).subs(u2, x0[1])).subs(t_symb, t_now))]]) 
				W = LA.inv(W) 
				v = W.dot(f0) 
				x_next = x0 - v
		return x_next
	else: 
		x_next = np.array([1000, 1000, 1000]) 
		while abs(x_next[0] - x0[0]) > 0.01 and abs(x_next[1] - x0[1]) > 0.01 and abs(x_next[2] - x0[2]) > 0.01: 
			c = c + 1
			if c > 5:
				return x_next
			else: 
				x_next = x0 
				x0 = x_next 
				f0 = np.array([float(((F[0].subs(u1, x0[0])).subs(u2, x0[1])).subs(u3, x0[2])), float(((F[1].subs(u1, x0[0])).subs(u2, x0[1])).subs(u3, x0[2])), float(((F[2].subs(u1, x0[0])).subs(u2, x0[1])).subs(u3, x0[2]))]) 
				W = np.array([[float(((diff(F[0], u1).subs(u1, x0[0])).subs(u2, x0[1])).subs(u3, x0[2])), float(((diff(F[0], u2).subs(u1, x0[0])).subs(u2, x0[1])).subs(u3, x0[2])), float(((diff(F[0], u3).subs(u1, x0[0])).subs(u2, x0[1])).subs(u3, x0[2]))], [float(((diff(F[1], u1).subs(u1, x0[0])).subs(u2, x0[1])).subs(u3, x0[2])), float(((diff(F[1], u2).subs(u1, x0[0])).subs(u2, x0[1])).subs(u3, x0[2])), float(((diff(F[1], u3).subs(u1, x0[0])).subs(u2, x0[1])).subs(u3, x0[2]))], [float(((diff(F[2], u1).subs(u1, x0[0])).subs(u2, x0[1])).subs(u3, x0[2])), float(((diff(F[2], u2).subs(u1, x0[0])).subs(u2, x0[1])).subs(u3, x0[2])), float(((diff(F[2], u3).subs(u1, x0[0])).subs(u2, x0[1])).subs(u3, x0[2]))]]) 
				W = LA.inv(W) 
				v = W.dot(f0) 
				x_next = x0 - v 
		return x_next 
 
def implicit_method(count, u, tau_max, tau_min, T, eps, number): 
	t_now = tau_max 
	t_next = t_now 
	y_now = u 
	y_back = u 
	y_next = u 
	tau_back = tau_min 
	tau_now = tau_min 
	u1, u2, u3, t_symb = symbols("u1, u2, u3, t_symb") 
	vctr_F = { 
		0: [-u1 * u2 + (sin(t_symb)) / t_symb, -pow(u2, 2) + (3.5 * t_symb) / (1 + pow(t_symb, 2))], 
		1: [u2 - (2.25 * u1 + u2) * u1, exp(u1) - (u1 + 2.25 * u2) * u1], 
		2: [(2.25 - 1.25) / 1.25 * u2 * u3, (1.25 + 2.25) / 2.25 * u1 * u3, (1.25 - 2.25) / 1.25 * u1 * u2] 
	} 
	F = np.array(vctr_F[number])
	while t_now < T: 
		t_next = t_now + tau_now 
		vector_F = np.array([u1 - y_now[0] - F[0]*tau_now, u2 - y_now[1] - F[1]*tau_now]) 
		if count == 3: 
			vector_F = np.array([u1 - y_now[0] - F[0]*tau_now, u2 - y_now[1] - F[1]*tau_now, u3 - y_now[2] - F[2]*tau_now]) 
		x0 = y_now.transpose() 
		y_next = Newton(count, x0, vector_F, t_now).transpose() 
		eps_k = -(tau_now/(tau_now + tau_back))*(y_next - y_now - (tau_now/tau_back)*(y_now-y_back))
		tau_next_vector = np.zeros(count)
		for i in range (0, count, 1): 
			if abs(eps_k[i]) > eps: 
				tau_next_vector[i] = tau_now/2 
			else: 
				if eps/4 < abs(eps_k[i]) and abs(eps_k[i]) <= eps:
					tau_next_vector[i] = tau_now
				else: 
					tau_next_vector[i] = 2*tau_now
		tau_next = min(tau_next_vector)
		if tau_next > tau_max: 
			tau_next = tau_max 
		plt.plot(t_now, y_next[0], ',', color='b') 
		plt.plot(t_now, y_next[1], ',', color='r')
		if number == 2: 
			plt.plot(t_now, y_next[2], ',', color='y') 
		y_back = y_now 
		y_now = y_next 
		tau_back = tau_now 
		tau_now = tau_next 
		t_now = t_next 
  
def shihman_method(count, u, tau_max, tau, T, eps, number): 
	alpha_0 = 1 
	alpha_1 = 0 
	beta_0 = tau 
	t_now = tau 
	y_now = u 
	y_back = u 
	y_next = u 
	tau_back = tau 
	tau_now = tau 
	u1, u2, u3, t_symb = symbols("u1, u2, u3, t_symb") 
	vctr_F = { 
		0: [-u1 * u2 + (sin(t_symb)) / t_symb, -pow(u2, 2) + (3.5 * t_symb) / (1 + 150 pow(t_symb, 2))], 
		1: [u2 - (2.25 * u1 + u2) * u1, exp(u1) - (u1 + 2.25 * u2) * u1], 
		2: [(2.25 - 1.25) / 1.25 * u2 * u3, (1.25 + 2.25) / 2.25 * u1 * u3, (1.25 - 2.25) / 1.25 * u1 * u2] 
	} 
	F = np.array(vctr_F[number]) 
	c = 0 
	while t_now < T: 
		c = c = 1 
		if c > 2: 
			alpha_0 = pow((tau_now + tau_back), 2) / (tau_back * (2 * tau_now + tau_back)) 
			alpha_1 = -pow(tau_now, 2) / (tau_back * (2 * tau_now + tau_back)) 
			beta_0 = tau_now * (tau_now + tau_back) / (2 * tau_now + tau_back) 
		t_next = t_now + tau_now
		vector_F = np.array([u1 - alpha_1 * y_back[0] - alpha_0 * y_now[0] - beta_0 * F[0], u2 - alpha_1 * y_back[1] - alpha_0 * y_now[1] - beta_0 * F[1]]) 
		if count == 3: 
			vector_F = np.array( 
				[u1 - alpha_1 * y_back[0] - alpha_0 * y_now[0] - beta_0 * F[0], 
				u2 - alpha_1 * y_back[1] - alpha_0 * y_now[1] - beta_0 * F[1], 
				u3 - alpha_1 * y_back[2] - alpha_0 * y_now[2] - beta_0 * F[2]]) 
		x0 = y_now.transpose() 
		y_next = Newton(count, x0, vector_F, t_now).transpose()
		eps_k = -(tau_now/(tau_now + tau_back))*(y_next - y_now - (tau_now/tau_back)*(y_now-y_back))
		tau_next_vector = np.zeros(count)
		for i in range(0, count, 1):
			if abs(eps_k[i]) > eps: 
				tau_next_vector[i] = tau_now / 2
			else: 
				if eps / 4 < abs(eps_k[i]) and abs(eps_k[i]) <= eps:
					tau_next_vector[i] = tau_now
				else: 
					tau_next_vector[i] = 2 * tau_now
		tau_next = min(tau_next_vector)
		if tau_next > tau_max: 
			tau_next = tau_max 
		plt.plot(t_now, y_next[0], ',', color='b')
		plt.plot(t_now, y_next[1], ',', color='r')
		if number == 2: 
			plt.plot(t_now, y_next[2], ',', color='y') 
		y_back = y_now 
		y_now = y_next 
	tau_back = tau_now 
	tau_now = tau_next 
	t_now = t_next 
  
def main(): 
	tau_max = 0.001 
	tau_min = 0.001 
	eps = 0.001 
	T = 1 
	t = 0.01 
	explicit_method(2, init_u(0), tau_max, T, t, eps, 0) 
	plt.title("Явный метод Эйлера") 
	plt.xlabel("t") 
	plt.ylabel("u1(t), u2(t)") 
	plt.figure() 
	implicit_method(2, init_u(0), tau_max, tau_min, T, eps, 0) 
	plt.title("Неявный метод Эйлера") 
	plt.xlabel("t") 
	plt.ylabel("u1(t), u2(t)") 
	plt.figure() 
	shihman_method(2, init_u(0), tau_max, tau_min, T, eps, 0) 
	plt.title("Метод Шихмана") 
	plt.xlabel("t") 
	plt.ylabel("u1(t), u2(t)") 
	plt.figure() 
	explicit_method(2, init_u(1), tau_max, T, t, eps, 1) 
	plt.title("Явный метод Эйлера") 
	plt.xlabel("t") 
	plt.ylabel("u1(t), u2(t)") 
	plt.figure() 
	implicit_method(2, init_u(1), tau_max, tau_min, T, eps, 1) 
	plt.title("Неявный метод Эйлера") 
	plt.xlabel("t") 
	plt.ylabel("u1(t), u2(t)") 
	plt.figure() 
	shihman_method(2, init_u(1), tau_max, tau_min, T, eps, 1) 
	plt.title("Метод Шихмана") 
	plt.xlabel("t") 
	plt.ylabel("u1(t), u2(t)") 
	plt.figure() 
	explicit_method(3, init_u(2), tau_max, T, t, eps, 2) 
	plt.title("Явный метод Эйлера") 
	plt.xlabel("t") 
	plt.ylabel("u1(t), u2(t), u3(t)") 
	plt.figure() 
	implicit_method(3, init_u(2), tau_max, tau_min, T, eps, 2) 
	plt.title("Неявный метод Эйлера") 
	plt.xlabel("t") 
	plt.ylabel("u1(t), u2(t), u3(t)") 
	plt.figure() 
	shihman_method(3, init_u(2), tau_max, tau_min, T, eps, 2) 
	plt.title("Метод Шихмана") 
	plt.xlabel("t") 
	plt.ylabel("u1(t), u2(t), u3(t)") 
	plt.show() 
 
if (__name__ == "__main__"): 
	main() 

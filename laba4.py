import numpy as np
from numpy import linalg as LA
import math
import matplotlib.pyplot as plt


def explicit_method(Lx, Ly, T, c):
    delta_x = 0.1
    delta_y = 0.1
    delta_t = (pow(delta_x, 2)/c)*0.1
    length_x = int(Lx/delta_x + 1)
    length_y = int(Ly/delta_y + 1)
    x = np.linspace(0, Lx, length_x)
    y = np.linspace(0, Ly, length_y)
    u = np.zeros((length_x, length_y))
    for i in range(0, length_x, 1):
        for j in range(0, length_y, 1):
            u[i][j] = 30
    X, Y = np.meshgrid(x, y)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    plt.ion()
    ax.plot_surface(X, Y, u, cmap='jet')
    t = 0
    u_x = np.zeros((length_x, length_y))
    u_y = np.zeros((length_x, length_y))
    while t <= T:
        for i in range (0, length_x, 1):
            for j in range (0, length_y, 1):
                if i == 0:
                    u_x[i][j] = (u[i+1][j] - 2 * u[i][j] + 30) / pow(delta_x, 2)
                else:
                    if i == length_x - 1:
                        u_x[i][j] = (100 - 2 * u[i][j] + u[i-1][j]) / pow(delta_x, 2)
                    else:
                        u_x[i][j] = (u[i+1][j] - 2 * u[i][j] + u[i-1][j]) / pow(delta_x, 2)
                if j == 0:
                    u_y[i][j] = (u[i][j+1] - u[i][j]) / pow(delta_y, 2)
                else:
                    if j == length_y - 1:
                        u_y[i][j] = (-u[i][j] + u[i][j-1]) / pow(delta_y, 2)
                    else:
                        u_y[i][j] = (u[i][j+1] - 2 * u[i][j] + u[i][j-1]) / pow(delta_y, 2)
        u = u + delta_t * (c * (u_x + u_y) + 100*math.sin(10*t))
        ax.clear()
        ax.plot_surface(X, Y, u, cmap='jet')
        plt.draw()
        plt.pause(delta_t/3)
        print(t)
        t = t + delta_t


def implicit_method(Lx, Ly, T, c):
    delta_x = 0.1
    delta_y = 0.1
    delta_t = (pow(delta_x, 2)/c)*0.1
    length_x = int(Lx/delta_x + 1)
    length_y = int(Ly/delta_y + 1)
    x = np.linspace(0, Lx, length_x)
    y = np.linspace(0, Ly, length_y)
    u = np.zeros((length_x, length_y))
    for i in range(0, length_x, 1):
        for j in range(0, length_y, 1):
            u[i][j] = 30
    X, Y = np.meshgrid(x, y)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    plt.ion()
    ax.plot_surface(X, Y, u, cmap='jet')
    t = 0
    while t <= T:
        u_coeff = np.zeros((length_x, length_x))
        u_x_prev = np.zeros(length_x)
        for i in range (0, length_x, 1):
            if i == 0:
                u_coeff[i][i] = 1 + 2 * c * delta_t/pow(delta_x, 2)
                u_coeff[i][i+1] = -c * (delta_t/pow(delta_x, 2))
                u_x_prev[i] = u[i][1] + delta_t * (c * 30 / pow(delta_x, 2) - 100 * math.sin(10 * (t + delta_t)))
            else:
                if i == length_x - 1:
                    u_coeff[i][i-1] = -c * (delta_t/pow(delta_x, 2))
                    u_coeff[i][i] = 1 + 2 * c * delta_t / pow(delta_x, 2)
                    u_x_prev[i] = u[i][1] + delta_t * (c * 100 / pow(delta_x, 2) - 100 * math.sin(10 * (t + delta_t)))
                else:
                    u_coeff[i][i - 1] = -c * (delta_t / pow(delta_x, 2))
                    u_coeff[i][i] = 1 + 2 * c * delta_t / pow(delta_x, 2)
                    u_coeff[i][i + 1] = -c * (delta_t / pow(delta_x, 2))
                    u_x_prev[i] = u[i][1] + delta_t * 100 * math.sin(10 * (t + delta_t))
        u_x_next = LA.inv(u_coeff).dot(u_x_prev.transpose())
        for j in range (0, length_y, 1):
            for i in range (0, length_y, 1):
                u[i][j] = u_x_next[i]
        ax.clear()
        ax.plot_surface(X, Y, u, cmap='jet')
        plt.draw()
        plt.pause(delta_t/3)
        print(t)
        t = t + delta_t


def main():
    Lx = float(input("Lx: "))
    Ly = float(input("Ly: "))
    T = float(input("T: "))
    c = float(input("c: "))
    #explicit_method(Lx, Ly, T, c)
    implicit_method(Lx, Ly, T, c)
    plt.ioff()
    plt.show()


if __name__ == '__main__':
    main()

#include <iostream>
#include <math.h>
using namespace std;

void inv(double** A, int N) {
    double temp;
    double** E = new double* [N];
    for (int i = 0; i < N; i++)
        E[i] = new double[N];
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++) {
            E[i][j] = 0.0;
            if (i == j)
                E[i][j] = 1.0;
        }
    for (int k = 0; k < N; k++) {
        temp = A[k][k];
        for (int j = 0; j < N; j++) {
            A[k][j] /= temp;
            E[k][j] /= temp;
        }
        for (int i = k + 1; i < N; i++) {
            temp = A[i][k];
            for (int j = 0; j < N; j++) {
                A[i][j] -= A[k][j] * temp;
                E[i][j] -= E[k][j] * temp;
            }
        }
    }
    for (int k = N - 1; k > 0; k--) {
        for (int i = k - 1; i >= 0; i--) {
            temp = A[i][k];
            for (int j = 0; j < N; j++) {
                A[i][j] -= A[k][j] * temp;
                E[i][j] -= E[k][j] * temp;
            }
        }
    }
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            A[i][j] = E[i][j];
    for (int i = 0; i < N; i++)
        delete[] E[i];
    delete[] E;
}

void explicit_method(double Lx, double Ly, double T, double c) {
    double delta_x = 0.1, delta_y = 0.1, t = 0;
    double delta_t = 0.1 * pow(delta_x, 2) / c;
    int length_x = Lx / delta_x + 1, length_y = Ly / delta_y + 1;
    double* x = new double [length_x];
    double* y = new double [length_y];
    double** u = new double* [length_x];
    double** u_x = new double* [length_x];
    double** u_y = new double* [length_x];
    double** u_partial_derivativies = new double* [length_x];
    for (int i = 0; i < length_x; i++) {
        u[i] = new double [length_y];
        u_x[i] = new double [length_y];
        u_y[i] = new double [length_y];
        u_partial_derivativies[i] = new double [length_y];
        x[i] = i * delta_x;
        y[i] = i * delta_y;
        for (int j = 0; j < length_y; j++) {
            u[i][j] = 30;
            u_x[i][j] = 0;
            u_y[i][j] = 0;
            cout << u[i][j] << " ";
        }
        cout << endl;
    }
    while (t <= T) {
        system("cls");
        for (int i = 0; i < length_x; i++) {
            for (int j = 0; j < length_y; j++) {
                if (i == 0)
                    u_x[i][j] = (u[i + 1][j] - 2 * u[i][j] + 30) / pow(delta_x, 2);
                else if (i == length_x - 1)
                    u_x[i][j] = (100 - 2 * u[i][j] + u[i - 1][j]) / pow(delta_x, 2);
                else
                    u_x[i][j] = (u[i + 1][j] - 2 * u[i][j] + u[i - 1][j]) / pow(delta_x, 2);
                if (j == 0)
                    u_y[i][j] = (u[i][j + 1] - u[i][j]) / pow(delta_y, 2);
                else if (j == length_y - 1)
                    u_y[i][j] = (-u[i][j] + u[i][j - 1]) / pow(delta_y, 2);
                else
                    u_y[i][j] = (u[i][j + 1] - 2 * u[i][j] + u[i][j - 1]) / pow(delta_y, 2);
                u_partial_derivativies[i][j] = u_x[i][j] + u_y[i][j];
                u[i][j] = u[i][j] + delta_t * (c + u_partial_derivativies[i][j] + 100 * sin(10 * t));
                cout << u[i][j] << " ";
            }
            cout << endl;
        }
    }
}

void implicit_method(double Lx, double Ly, double T, double c) {
    double delta_x = 0.1, delta_y = 0.1, t = 0;
    double delta_t = 0.1 * pow(delta_x, 2) / c;
    int length_x = Lx / delta_x + 1, length_y = Ly / delta_y + 1;
    double* x = new double[length_x];
    double* y = new double[length_y];
    double* u_x_prev = new double[length_x];
    double* u_x_next = new double[length_x];
    double** u = new double* [length_x];
    double** u_coeff = new double* [length_x];
    for (int i = 0; i < length_x; i++) {
        u[i] = new double[length_y];
        u_coeff[i] = new double[length_x];
        x[i] = i * delta_x;
        y[i] = i * delta_y;
        u_x_prev[i] = 0;
        u_x_next[i] = 0;
        for (int j = 0; j < length_y; j++) {
            u[i][j] = 30;
            u_coeff[i][j] = 0;
            cout << u[i][j] << " ";
        }
        cout << endl;
    }
    while (t <= T) {
        for (int i = 0; i < length_x; i++) {
            u_x_next[i] = 0;
            if (i == 0) {
                u_coeff[i][i] = 1 + 2 * c * delta_t / pow(delta_x, 2);
                u_coeff[i][i + 1] = -c * (delta_t / pow(delta_x, 2));
                u_x_prev[i] = u[i][1] + delta_t * (c * 30 / pow(delta_x, 2) - 100 * sin(10 * (t + delta_t)));
            }
            else if (i == length_x - 1) {
                u_coeff[i][i - 1] = -c * (delta_t / pow(delta_x, 2));
                u_coeff[i][i] = 1 + 2 * c * delta_t / pow(delta_x, 2);
                u_x_prev[i] = u[i][1] + delta_t * (c * 100 / pow(delta_x, 2) - 100 * sin(10 * (t + delta_t)));
            }
            else {
                u_coeff[i][i - 1] = -c * (delta_t / pow(delta_x, 2));
                u_coeff[i][i] = 1 + 2 * c * delta_t / pow(delta_x, 2);
                u_coeff[i][i + 1] = -c * (delta_t / pow(delta_x, 2));
                u_x_prev[i] = u[i][1] + delta_t * 100 * sin(10 * (t + delta_t));
            }
        }
        inv(u_coeff, length_x);
        for (int i = 0; i < length_x; i++) {
            for (int j = 0; j < length_y; j++) {
                u_x_next[i] = u_x_next[i] + u_coeff[i][j] * u_x_prev[j];
            }
        }
        system("cls");
        for (int i = 0; i < length_x; i++) {
            for (int j = 0; j < length_y; j++) {
                u[i][j] = u_x_next[i];
                cout << u[i][j] << " ";
            }
            cout << endl;
        }
    }
}

int main() {
    double Lx, Ly, T, c;
    cout << "Lx: " << endl;
    cin >> Lx;
    cout << "Ly: " << endl;
    cin >> Ly;
    cout << "T: " << endl;
    cin >> T;
    cout << "c: " << endl;
    cin >> c;
    explicit_method(Lx, Ly, T, c);
    implicit_method(Lx, Ly, T, c);
}

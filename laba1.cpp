#include <iostream> 
using namespace std;

void inversion(long double** A, int N) {
	long double temp;
	long double** E = new long double* [N];
	for (int i = 0; i < N; i++)
		E[i] = new long double[N];
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
45
int main() {
	double e = 0.00001, lyambda = 0, lyambda_next = 1000; 48     int len = 0;
	cout << "Matrix size: " << endl;
	cin >> len;
	long double** A = new long double* [len];
	for (int i = 0; i < len; i++) {
		A[i] = new long double[len];
	}
	cout << "Matrix: " << endl;
	for (int i = 0; i < len; i++) {
		for (int j = 0; j < len; j++) {
			A[i][j] = 0;
			cin >> A[i][j];
		}
	}
	cout << "U: " << endl;
	long double* u_0 = new long double[len];
	for (int i = 0; i < len; i++) {
		u_0[i] = 0;
		cin >> u_0[i];
	}
	long double* u = new long double[len];
	for (int i = 0; i < len; i++) {
		u[i] = u_0[i];
	}
	long double* u_next = new long double[len];
	for (int i = 0; i < len; i++) {
		u_next[i] = 0;
	}
	while (abs(lyambda_next - lyambda) > e) {
		lyambda = lyambda_next;
		for (int i = 0; i < len; i++) {
			for (int j = 0; j < len; j++) {
				u_next[i] += A[i][j] * u[j];
			}
		}
		long double numerator = 0, denominator = 0;
		for (int i = 0; i < len; i++) {
			denominator += u[i] * u[i];
			numerator += u_next[i] * u[i];
		}
		lyambda_next = numerator / denominator;
		for (int i = 0; i < len; i++) {
			u[i] = u_next[i];
			u_next[i] = 0;
		}
	}
	cout << "Lyambda_max = " << lyambda_next << endl;
	inversion(A, len);
	lyambda = 1000, lyambda_next = 0;
	for (int i = 0; i < len; i++) {
		u[i] = u_0[i];
		u_next[i] = 0;
	}
	while (abs(lyambda_next - lyambda) > e) {
		lyambda = lyambda_next;
		for (int i = 0; i < len; i++) {
			for (int j = 0; j < len; j++) {
				u_next[i] += A[i][j] * u[j];
			}
		}
		long double numerator = 0, denominator = 0;
		for (int i = 0; i < len; i++) {
			denominator += u[i] * u[i];
			numerator += u_next[i] * u[i];
		}
		lyambda_next = numerator / denominator;
		for (int i = 0; i < len; i++) {
			u[i] = u_next[i];
			u_next[i] = 0;
		}
	}
	cout << "Lyambda_min = " << lyambda_next << endl;
}

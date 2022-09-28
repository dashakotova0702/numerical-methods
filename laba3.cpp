#include <iostream> 
using namespace std;
int main() {
	long double A[3][3] = { {2.31, 31.49, 1.52}, {4.21, 22.42, 3.85}, {3.49, 4.85, 28.72} };
	long double B[3] = { 40.95, 30.24, 42.81 };
	int len = sizeof(A) / sizeof(long double);
	len = sqrt(len);
	if (len != sizeof(B) / sizeof(long double)) {
		return 1;
	}
	long double** A_norm = new long double* [len];
	for (int i = 0; i < len; i++) {
		A_norm[i] = new long double[len];
	}
	long double* B_norm = new long double[len];
	for (int i = 0; i < len; i++) {
		for (int j = 0; j < len; j++) {
			A_norm[i][j] = A[i][j];
		}
		B_norm[i] = B[i];
	}
	long double* x = new long double[len];
	for (int i = 0; i < len; i++) {
		x[i] = 0;
	}
	for (int i = 0; i < len; i++) {
		for (int j = 0; j < len; j++) {
			if (j - i <= 0) {
				continue;
			}
			long double a = A[i][i];
			long double b = A[j][i];
			for (int k = 0; k < len; k++) {
				A[j][k] = a * A[j][k] - b * A[i][k];  	 	 	 	 	
			}
			B[j] = a * B[j] - b * B[i];
		}
	}
	for (int i = len - 1; i >= 0; i--) {
		x[i] = B[i];
		for (int j = 1; j < len; j++) {
			if (j != i) {
				x[i] = x[i] - A[i][j] * x[j];
			}
		}
		x[i] = x[i] / A[i][i];
	}
	long double* F = new long double[len];
	for (int i = 0; i < len; i++) {
		F[i] = 0;
	}
	for (int i = 0; i < len; i++) {
		for (int j = 0; j < len; j++) {
			F[i] = F[i] + A_norm[i][j] * x[j];
		}
		F[i] = F[i] - B_norm[i];
	}
	long double norm = 0;
	for (int i = 0; i < len; i++) {
		if (F[i] > norm) {
			norm = F[i];
		}
	}
	for (int i = 0; i < len; i++) {
		cout << "x[" << i << "]: " << x[i] << " " << endl;
	}
	for (int i = 0; i < len; i++) {
		cout << "F[" << i << "]: " << F[i] << " " << endl;
	}
	cout << "Norm: " << norm << endl;
	return 0;
}

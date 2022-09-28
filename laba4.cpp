#include <iostream> 
using namespace std;

bool max(int len, long double* x_1, long double* x_2) {
	long double max = 0;
	for (int i = 0; i < len; i++) {
		if (abs(x_1[i] - x_2[i]) > max)
			max = abs(x_1[i] - x_2[i]);
	}
	if (max >= 0.001)
		return true;
	else
		return false;
}
15
int main() {
	long double A[3][3] = { {2.31, 31.49, 1.52}, {4.21, 22.42, 3.85}, {3.49, 4.85,
	28.72} };
	long double B[3] = { 40.95, 30.24, 42.81 };
	int len = sizeof(A) / sizeof(long double);
	len = sqrt(len);
	bool diag_preob = false;
	long double sum = 0;
	for (int i = 0; i < len; i++) {
		for (int j = 0; j < len; j++) {
			sum = 0;
			if (i != j)
				sum += abs(A[j][i]);
		}
		if (abs(A[i][i]) <= sum)
			diag_preob = true;
	}
	if (!diag_preob) {
		long double* x_1 = new long double[len];
		long double* x_2 = new long double[len];                                      
			for (int i = 0; i < len; i++) {
				x_1[i] = 0;
				x_2[i] = 0;
			}
		while (max(len, x_1, x_2) || x_2[0] - x_1[0] == 0) {
			for (int i = 0; i < len; i++)
				x_1[i] = x_2[i];
			for (int i = 0; i < len; i++) {
				x_2[i] = B[i] / A[i][i];
				for (int j = 0; j < len; j++) {
					if (i != j)
						x_2[i] = x_2[i] - x_1[j] * A[i][j] / A[i][i];
				}
			}
		}
		cout << "Method Yakob: ";
		for (int i = 0; i < len; i++)
			cout << x_2[i] << " ";
		while (max(len, x_1, x_2) || x_2[0] - x_1[0] == 0) {
			for (int i = 0; i < len; i++)
				x_1[i] = x_2[i];
			for (int i = 0; i < len; i++) {
				x_2[i] = B[i] / A[i][i];
				for (int j = 0; j < len; j++) {
					if (i != j)
						x_2[i] = x_2[i] - x_2[j] * A[i][j] / A[i][i];
				}
			}
		}
		cout << endl << "Method Gauss-Zeidel: ";
		for (int i = 0; i < len; i++)
			cout << x_2[i] << " ";
	}
	else {
		return 1;
	}
}

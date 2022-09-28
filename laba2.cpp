#include <iostream> 
#include "math.h" 
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

int main() {
	double e = 0.000000001, d1 = 10.00, d2 = 10.00;
	int NIT = 0, i = 0;
	long double* Apprxmtn_1 = new long double[2];
	long double* Apprxmtn_2 = new long double[2];
	long double* Apprxmtn_1_next = new long double[2];
	long double* Apprxmtn_2_next = new long double[2];
	long double* temp = new long double[2];
	long double* d2_v = new long double[2];
	long double** W = new long double* [2];
	for (int i = 0; i < 2; i++) {
		W[i] = new long double[2];
	}
	long double* discrepancy_vector = new long double[2];
	cout << "Initial approximation 1: " << endl;
	for (int i = 0; i < 2; i++)
		cin >> Apprxmtn_1[i];
	cout << "Initial approximation 2: " << endl;
	for (int i = 0; i < 2; i++)
		cin >> Apprxmtn_2[i];
	cout << "NIT: " << endl;
	cin >> NIT;
	cout << "No         d1          d2" << endl;
	while (d1 >= e && d2 >= e) {
		i += 1;
		discrepancy_vector[0] = cos(0.40 * Apprxmtn_1[1] + pow(Apprxmtn_1[0], 2)) + pow(Apprxmtn_1[1], 2) + pow(Apprxmtn_1[0], 2) - 1.6;
		discrepancy_vector[1] = 1.5 * pow(Apprxmtn_1[0], 2) - pow(Apprxmtn_1[1], 2) / 0.36 - 1;
		W[0][0] = -2 * Apprxmtn_1[0] * sin(0.4 * Apprxmtn_1[1] + pow(Apprxmtn_1[0], 2)) + 2 * Apprxmtn_1[0];
		W[1][0] = 3 * Apprxmtn_1[0];
		W[0][1] = -0.4 * sin(0.4 * Apprxmtn_1[1] + pow(Apprxmtn_1[0], 2)) + 2 * Apprxmtn_1[1];
		W[1][1] = -2 * Apprxmtn_1[1] / 0.36;
		inversion(W, 2);
		temp[0] = W[0][0] * discrepancy_vector[0] + W[0][1] * discrepancy_vector[1];
		temp[1] = W[1][0] * discrepancy_vector[0] + W[1][1] * discrepancy_vector[1];
		Apprxmtn_1_next[0] = Apprxmtn_1[0] - temp[0];
		Apprxmtn_1_next[1] = Apprxmtn_1[1] - temp[1];
		d1 = max(abs(discrepancy_vector[0]), abs(discrepancy_vector[1]));
		for (int j = 0; j < 2; j++) {
			if (Apprxmtn_1_next[j] < 1)
				d2_v[j] = Apprxmtn_1_next[j] - Apprxmtn_1[j];
			else
				d2_v[j] = (Apprxmtn_1_next[j] - Apprxmtn_1[j]) / Apprxmtn_1_next[j];
		}
		d2 = max(abs(d2_v[0]), abs(d2_v[1]));
		cout << i << "        " << d1 << "          " << d2 << endl;
		Apprxmtn_1[0] = Apprxmtn_1_next[0];
		Apprxmtn_1[1] = Apprxmtn_1_next[1];
		if (i > NIT) {
			cout << "IER = 2" << endl;
			break;
		}
	}
	cout << "x1: " << Apprxmtn_1[0] << "    x2: " << Apprxmtn_1[1] << endl;
	d1 = 10.00, d2 = 10.00, i = 0;
	while (d1 >= e && d2 >= e) {
		i += 1;
		discrepancy_vector[0] = cos(0.40 * Apprxmtn_2[1] + pow(Apprxmtn_2[0], 2)) + pow(Apprxmtn_2[1], 2) + pow(Apprxmtn_2[0], 2) - 1.6;
		discrepancy_vector[1] = 1.5 * pow(Apprxmtn_2[0], 2) - pow(Apprxmtn_2[1], 2) / 0.36 - 1;
		W[0][0] = -2 * Apprxmtn_2[0] * sin(0.4 * Apprxmtn_2[1] + pow(Apprxmtn_2[0], 2)) + 2 * Apprxmtn_2[0];
		W[1][0] = 3 * Apprxmtn_2[0];
		W[0][1] = -0.4 * sin(0.4 * Apprxmtn_2[1] + pow(Apprxmtn_2[0], 2)) + 2 * Apprxmtn_2[1];
		W[1][1] = -2 * Apprxmtn_2[1] / 0.36;
		inversion(W, 2);
		temp[0] = W[0][0] * discrepancy_vector[0] + W[0][1] * discrepancy_vector[1];
		temp[1] = W[1][0] * discrepancy_vector[0] + W[1][1] * discrepancy_vector[1];
		Apprxmtn_2_next[0] = Apprxmtn_2[0] - temp[0];
		Apprxmtn_2_next[1] = Apprxmtn_2[1] - temp[1];
		d1 = max(abs(discrepancy_vector[0]), abs(discrepancy_vector[1]));
		for (int j = 0; j < 2; j++) {
			if (Apprxmtn_2_next[j] < 1)
				d2_v[j] = Apprxmtn_2_next[j] - Apprxmtn_2[j];
			else
				d2_v[j] = (Apprxmtn_2_next[j] - Apprxmtn_2[j]) / Apprxmtn_2_next[j];
		}
		d2 = max(abs(d2_v[0]), abs(d2_v[1]));
		cout << i << "        " << d1 << "          " << d2 << endl;
		Apprxmtn_2[0] = Apprxmtn_2_next[0];
		Apprxmtn_2[1] = Apprxmtn_2_next[1];
		if (i > NIT) {
			cout << "IER = 2" << endl;
			break;
		}
	}
	cout << "x1: " << Apprxmtn_2[0] << "    x2: " << Apprxmtn_2[1] << endl;
}

#include <stdio.h>
#include <math.h>
#include <memory.h>
#include <Windows.h>

void rdft(int n, int isgn, double *a, int *ip, double *w);

const double M_PI = 3.14159265358979323846;

int main()
{
	const int n = 8192; // data length (int), n >= 2, n = power of 2
	int ip[2 + 64]; // work area for bit reversal(int *), length of ip >= 2 + sqrt(n / 2)
	double w[n / 2 - 1];

	ip[0] = 0; // first time only

	double a[n]; // signal

	for (int i = 0; i < n; i++) {
		a[i] = sin(2.0 * M_PI * i / n * 10.0) + cos(2.0 * M_PI * i / n * 3.0);
	}

	double a_tmp[n];
	memcpy(a_tmp, a, sizeof(a));

	rdft(n, 1, a_tmp, ip, w);

	// M‰ñ‘ª’è
	const int M = 10000;

	DWORD time1 = GetTickCount();

	for (int i = 0; i < M; i++) {
		memcpy(a_tmp, a, sizeof(a));
		rdft(n, 1, a_tmp, ip, w);
	}

	DWORD time2 = GetTickCount();

	printf("time = %lu ms\n", time2 - time1);

	return 0;
}
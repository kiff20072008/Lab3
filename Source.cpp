#include "rkf45.hpp"

using namespace std;

int func(int n, double t1, double *y1, double *dy1)// Функция для метода rkf45 
{
	dy1[0] = -44 * y1[0] - 160 * y1[1] + cos(t1 + 1);
	dy1[1] = y1[0] + atan(1 + t1 * t1);
	return 0;
}

double adamsf0(double x, double y0, double y1)//первое уравнение для метода Адамса 
{
	return -44 * y0 - 160 * y1 + cos(x + 1);
}

double adamsf1(double x, double y0, double y1)//второе уравнение для метода Адамса 
{
	return y0 + atan(1 + x * x);
}

void rkf(double y0[20],double y1[20])
{
	double h, reller, abserr, x1, x2, y[2], yp[2],dh;
	int n, flag, nfe, maxfe, fail;

	//initiation 
	n = 2;
	flag = 1;
	maxfe = 5000;
	reller = 1.0e-4;
	abserr = 1.0e-4;
	dh = 0.08;
	rkfinit(n, &fail);

	if (!fail)
	{
		y[0] = 2;
		y[1] = 0.5;
		for (int i = 1; i < 21; ++i)
		{
			x2 = dh * i;
			x1 = dh * (i - 1);
			rkf45(func, n, y, yp, &x1, x2, &reller, abserr, &h, &nfe, maxfe, &flag);
			y0[i - 1] = y[0];
			y1[i - 1] = y[1];
			if (flag != 2)
			{
				cout << "Something wrong" << endl;
				break;
			}
		}
		rkfend();
	}
}

void adams1(double y0[20], double y1[20])
{
	double h, reller, abserr, x1, x2, dh, y[2], yp[2];
	int n, flag, nfe, maxfe, fail;

	//initiation 
	n = 2;
	flag = 1;
	maxfe = 5000;
	reller = 1.0e-4;
	abserr = 1.0e-4;
	rkfinit(n, &fail);

	if (!fail)
	{
		y[0] = 2;
		y[1] = 0.5;
		dh = 0.008;
		for (int i = 1; i < 5; ++i)
		{
			x2 = i * dh;
			x1 = (i - 1)*dh;
			rkf45(func, n, y, yp, &x1, x2, &reller, abserr, &h, &nfe, maxfe, &flag);
			y0[i - 1] = y[0];
			y1[i - 1] = y[1];
		}
		for (int i = 5; i <= 210; ++i)
		{
			x2 = i * dh;
			x1 = (i - 1)*dh;
			y[0] = y[0] + dh * (55 * adamsf0(x2, y[0], y[1]) - 59 * adamsf0(x2 - dh, y[0], y[1]) + 37 * adamsf0(x2 - 2 * dh, y[0], y[1]) - 9 * adamsf0(x2 - 3 * dh, y[0], y[1])) / 24;
			y[1] = y[1] + dh * (55 * adamsf1(x2, y[0], y[1]) - 59 * adamsf1(x2 - dh, y[0], y[1]) + 37 * adamsf1(x2 - 2 * dh, y[0], y[1]) - 9 * adamsf1(x2 - 3 * dh, y[0], y[1])) / 24;
			if ((i % 10 - 1) == 0)
			{
				y0[i/10 - 1] = y[0];
				y1[i/10 - 1] = y[1];
			}
		}
		rkfend();
	}
}

void adams2(double y0[20], double y1[20])
{
	double h, reller, abserr, x1, x2, dh, y[2], yp[2];
	int n, flag, nfe, maxfe, fail;

	//initiation 
	n = 2;
	flag = 1;
	maxfe = 5000;
	reller = 1.0e-4;
	abserr = 1.0e-4;
	rkfinit(n, &fail);

	if (!fail)
	{
		y[0] = 2;
		y[1] = 0.5;
		dh = 0.004;
		for (int i = 1; i < 5; ++i)
		{
			x2 = i * dh;
			x1 = (i - 1)*dh;
			rkf45(func, n, y, yp, &x1, x2, &reller, abserr, &h, &nfe, maxfe, &flag);
			y0[i - 1] = y[0];
			y1[i - 1] = y[1];
		}
		for (int i = 5; i <= 420; ++i)
		{
			x2 = i * dh;
			x1 = (i - 1)*dh;
			y[0] = y[0] + dh * (55 * adamsf0(x2, y[0], y[1]) - 59 * adamsf0(x2 - dh, y[0], y[1]) + 37 * adamsf0(x2 - 2 * dh, y[0], y[1]) - 9 * adamsf0(x2 - 3 * dh, y[0], y[1])) / 24;
			y[1] = y[1] + dh * (55 * adamsf1(x2, y[0], y[1]) - 59 * adamsf1(x2 - dh, y[0], y[1]) + 37 * adamsf1(x2 - 2 * dh, y[0], y[1]) - 9 * adamsf1(x2 - 3 * dh, y[0], y[1])) / 24;
			if ((i % 20 - 1) == 0)
			{
				y0[i / 20- 1] = y[0];
				y1[i / 20 - 1] = y[1];
			}
		}
		rkfend();
	}
}

void adams3(double y0[20], double y1[20])
{
	double h, reller, abserr, x1, x2, dh, y[2], yp[2];
	int n, flag, nfe, maxfe, fail;

	//initiation 
	n = 2;
	flag = 1;
	maxfe = 5000;
	reller = 1.0e-4;
	abserr = 1.0e-4;
	rkfinit(n, &fail);

	if (!fail)
	{
		y[0] = 2;
		y[1] = 0.5;
		dh = 0.08;
		for (int i = 1; i < 5; ++i)
		{
			x2 = i * dh;
			x1 = (i - 1)*dh;
			rkf45(func, n, y, yp, &x1, x2, &reller, abserr, &h, &nfe, maxfe, &flag);
			y0[i - 1] = y[0];
			y1[i - 1] = y[1];
		}
		for (int i = 5; i < 21; ++i)
		{
			x2 = i * dh;
			x1 = (i - 1)*dh;
			y[0] = y[0] + dh * (55 * adamsf0(x2, y[0], y[1]) - 59 * adamsf0(x2 - dh, y[0], y[1]) + 37 * adamsf0(x2 - 2 * dh, y[0], y[1]) - 9 * adamsf0(x2 - 3 * dh, y[0], y[1])) / 24;
			y[1] = y[1] + dh * (55 * adamsf1(x2, y[0], y[1]) - 59 * adamsf1(x2 - dh, y[0], y[1]) + 37 * adamsf1(x2 - 2 * dh, y[0], y[1]) - 9 * adamsf1(x2 - 3 * dh, y[0], y[1])) / 24;
			y0[i  - 1] = y[0];
		    y1[i  - 1] = y[1];
			
		}
		rkfend();
	}
}

int main()
{
	double y0rkf[20], y1rkf[20], y0ad1[20], y1ad1[20], y0ad2[20], y1ad2[20], y0ad3[20], y1ad3[20];
	rkf(y0rkf,y1rkf);
	adams1(y0ad1,y1ad1);
	adams2(y0ad2, y1ad2);
	adams3(y0ad3, y1ad3);
	double dh = 0.08;
	cout << "x          y0rkf     y0ad(0.008)  rkf-ad       y0ad(0.004)  rkf-ad       y0ad(0.08)      rkf-ad        y1rkf       y1ad(0.008)   rkf-ad       y1ad(0.004)   rkf-ad       y1ad(0.08)    rkf-ad" << endl;
	
	for(int i = 0; i < 20; ++i)
	{
		cout << setw(10) << left << (i+1)*dh <<setw(9)<<left << y0rkf[i] << "   " << setw(9) << left << y0ad1[i] << "   " << setw(11) << left <<abs(y0rkf[i]-y0ad1[i]) << "   " << setw(9) << left << y0ad2[i] << "   " << setw(11) << left << abs(y0rkf[i] - y0ad2[i]) << "   " << setw(12) << left << y0ad3[i] << "   " << setw(11) << left << abs(y0rkf[i] - y0ad3[i]) << "   " << setw(10) << left << y1rkf[i] << "   " << setw(10) << left << y1ad1[i] << "   " << setw(11) << left << abs(y1rkf[i] - y1ad1[i]) << "   " << setw(10) << left << y1ad2[i] << "   " << setw(11) << left << abs(y1rkf[i] - y1ad2[i]) << "   " << setw(10) << left << y1ad3[i] << "   " << setw(10) << left << abs(y1rkf[i] - y1ad3[i]) << endl;
	}
	system("pause");
	return 0;
}
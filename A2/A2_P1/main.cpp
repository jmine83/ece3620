// 11/21/2014 - ECE 3620 - Meine, Joel
// Assignment #2

#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <fstream>
using namespace std;

vector<double> f1 = {0,1,2,3,2,1};
vector<double> f2 = {-2,-2,-2,-2,-2,-2,-2};
vector<double> f3 = {1,-1,1,-1};
vector<double> f4 = {0,0,0,-3,-3};

vector<double> a = {15,12,5}; // Constants
vector<double> yi = {-3,2,1}; // Values - Initial
int O = a.size(); // Order of Differential Equation
double T = 0.001; // Time - Step Size
double ti = 0; // Time - Initial
double tf = 10; // Time - Final
double Dp = 10; // Data Points, No. of

// NOTE: The number of items in a and yi "must" be equal.

double pi = 3.14159265358979323846264338327950288419716939937510;

vector<double> ft(vector<double> t)
{
	vector<double> FT;
	double f;
	for (int i = 0; i < t.size(); i++)
	{
		f = sin(3*pi*t[i]);
		FT.push_back(f);
		f = 0;
	}
	return(FT);
}

vector<double> ht(vector<double> t)
{
	double LAMBDA1 = -1.1980924212362;
	double LAMBDA2 = -2.6038151575276;
	double OMEGA = 2.079748094628;
	double Ch1 = 0.33386426881482;
	double Ch2 = 0.25516523268429;
	double Ch3 = -0.33386426881443;
	
	vector<double> HT;
	double h;
	for (int i = 0; i < t.size(); i++)
	{
		h = (Ch1*cos(OMEGA*t[i]) + Ch2*sin(OMEGA*t[i]))*exp(LAMBDA1*t[i]) + Ch3*exp(LAMBDA2*t[i]);
		HT.push_back(h);
		h = 0;
	}
	return(HT);
}

vector<double> conv(vector<double> f, vector<double> h)
{
	if (f.size() < h.size())
		f.resize(h.size());
	else if (f.size() > h.size())
		h.resize(f.size());
	double sum;
	int n,m;
	int nU = 2*f.size()-2;
	int mU = f.size();
	vector<double> y;
	for (n = 0; n <= nU; n++)
	{
		sum = 0;
		for (m = 0; m < mU; m++)
		{
			if (n-m < 0 || n-m >= h.size())
				sum += 0;
			else
				sum += f[m]*h[n-m];
		}
		y.push_back(sum);
	}
	while (y.back() == 0)
		y.pop_back();
	return(y);
}

vector<double> matrixOP_vecmul(vector<vector<double>> M, vector<double> u)
{ // Calculates v(n) = M(nxm) * u(m)
	vector<double> v;
	for (int k = 0; k < O; k++)
		v.push_back(0);
	int i; int j;
	for (i = 0; i < O; i++){
		v[i] = 0;
		for (j = 0; j < O; j++){
			v[i] = v[i] + M[i][j] * u[j];
		}
	}
	return(v);
}

void print_vec(vector<double> v)
{
	cout << "{";
	for (int i = 0; i < v.size(); i++)
	{
		if (i != v.size() - 1)
			cout << v[i] << ",";
		else
			cout << v[i] << "}" << endl;
	}
}

vector<vector<double>> ODE_m(vector<double> t)
{
	vector<double> y0;
	vector<vector<double>> Y0;
	y0 = yi;
	Y0.push_back(y0);

	// Z = w + A*dt
	vector<double> z;
	int r, s;
	for (r = 0; r < O; r++) // Column Size, Initialize
		z.push_back(0);
	vector<vector<double>> Z, w, A;
	for (s = 0; s < O; s++) // Row Size, Initialize
	{
		Z.push_back(z);
		w.push_back(z);
		A.push_back(z);
	}
	int i, j;
	for (i = 0; i < O; i++) // w
		for (j = 0; j < O; j++)
			if (i == j)
				w[i][j] = 1;
			else
				w[i][j] = 0;
	int m, n;
	for (m = 0; m < O; m++) // A
		for (n = 0; n < O; n++)
			if (n == m + 1 && m != O - 1)
				A[m][n] = 1;
			else if (m == O - 1)
				A[m][n] = -a[n];
			else
				A[m][n] = 0;
			int p, q;
			for (p = 0; p < O; p++) // Z
				for (q = 0; q < O; q++)
					Z[p][q] = w[p][q] + A[p][q] * T;
			// y = Z*y
			for (int r = 0; r < t.size(); r++)
			{
				y0 = matrixOP_vecmul(Z, y0);
				Y0.push_back(y0);
			}
			Y0.pop_back();
			return(Y0);
}

int main()
{
	std::cout << "Assignment #2 - Problem 3" << std::endl;
	std::cout << "==================================================================" << std::endl;
	std::cout << "For the given functions..." << std::endl;
	std::cout << "f1 = "; print_vec(f1);
	std::cout << "f2 = "; print_vec(f2);
	std::cout << "f3 = "; print_vec(f3);
	std::cout << "f4 = "; print_vec(f4);
	
	vector<double> A, B, C, D, E;
	A = conv(f1, f1);
	B = conv(f1, f2);
	C = conv(f1, f3);
	D = conv(f2, f3);
	E = conv(f1, f4);
	
	std::cout << " " << std::endl;
	std::cout << "the resulting convolutions are..." << std::endl;
	std::cout << "(a) conv(f1,f1) = "; print_vec(A);
	std::cout << "(b) conv(f1,f2) = "; print_vec(B);
	std::cout << "(c) conv(f1,f3) = "; print_vec(C);
	std::cout << "(d) conv(f2,f3) = "; print_vec(D);
	std::cout << "(e) conv(f1,f4) = "; print_vec(E);

	std::cout << " " << std::endl;
	std::cout << "==================================================================" << std::endl;
	std::cout << "For the differential equation..." << std::endl;
	std::cout << "(D^3 + 5D^2 + 12D + 15)y0(t) = (D + 0.5)*f(t)" << std::endl;
	std::cout << "with initial conditions y0(0) = -3, Dy0(0) = 2, D^2y0(0) = 1 and" << std::endl;
	std::cout << "f(t) = sin(3*PI*t)*u(t)" << std::endl;
	std::cout << " " << std::endl;
	std::cout << "h(t) = [(0.334*cos(2.080*t)+0.255*sin(2.080*t))*exp(-1.198*t)..." << std::endl;
	std::cout << "       ...-0.334*exp(-2.604*t)]*u(t)" << std::endl;
	std::cout << " " << std::endl;
	std::cout << "Time - Step Size, dt = " << T << std::endl;
	std::cout << "Time - Initial, ti = " << ti << std::endl;
	std::cout << "Time - Final, tf = " << tf << std::endl;
	std::cout << "------------------------------------------------------------------" << std::endl;
	std::cout << "    t              y0                y             y[tot]         " << std::endl;
	std::cout << "------------------------------------------------------------------" << std::endl;

	double I;
	I = (tf / T) / Dp; // Increment Length, I

	vector<double> t, sT; // Time List, t; Time List (Shortened), sT
	for (int j = 0; j <= tf / T; j++)
	{
		t.push_back(j*T+ti);
	}
	for (int k = 0; k <= Dp; k++)
	{
		sT.push_back(t[k*I]);
	}

	vector<vector<double>> Q;
	vector<double> P; // Results of ODE with only values of y0(t) and not t.
	Q = ODE_m(t); // Results of ODE including t and y0(t).
	int f, g;
	for (f = 0; f < t.size(); f++)
		for (g = 0; g < t.size(); g++)
			if (g == 0)
				P.push_back(Q[f][g]); // Retrieves Only the Results of y0(t).

	vector<double> hC, fC;
	hC = ht(t); // Results of function h(t) at time t.
	fC = ft(t); // Results of function f(t) at time t.
	C = conv(hC,fC); // Results of convolution of h(t) and f(t).

	vector<double> S;
	double a;
	for (int z = 0; z < t.size(); z++)
	{
		a = P[z] + C[z];
		S.push_back(a);
		a = 0;
	} // Returns the total solution.

	vector<double> sP1, sP2, sP3; // Print Lists (Shortened)
	for (int m = 0; m <= Dp; m++)
	{
		sP1.push_back(P[m*I]); // y0
		sP2.push_back(C[m*I]); // y
		sP3.push_back(S[m*I]); // y[tot]
	}
	for (int n = 0; n <= Dp; n++) // Prints results to window.
	{
		if (fmod(sT[0],1) && fmod(sT[1],1) != 0) // When sT list does not contain whole number values.
			printf(" %7.3f     %14.11f   %14.11f   %14.11f \n", sT[n], sP1[n], sP2[n], sP3[n]);
		else
			printf(" %7.0f     %14.11f   %14.11f   %14.11f \n", sT[n], sP1[n], sP2[n], sP3[n]);
	}
	ofstream outfile_y0("t_y0.txt");
	for (int o = 0; o < t.size(); o++)
	{
		outfile_y0 << t[o] << "  " << P[o] << endl; // Writes t and y0(t) values to file.
	}
	ofstream outfile_y("t_y.txt");
	for (int p = 0; p < t.size(); p++)
	{
		outfile_y << t[p] << "  " << C[p] << endl; // Writes t and y(t) values to file.
	}
	ofstream outfile_yTOT("t_yTOT.txt");
	for (int q = 0; q < t.size(); q++)
	{
		outfile_yTOT << t[q] << "  " << S[q] << endl; // Writes t and y[tot](t) values to file.
	}
	ofstream outfile_h("t_h.txt");
	for (int r = 0; r < t.size(); r++)
	{
		outfile_h << t[r] << "  " << hC[r] << endl; // Writes t and h values to file.
	}
	ofstream outfile_f("t_f.txt");
	for (int s = 0; s < t.size(); s++)
	{
		outfile_f << t[s] << "  " << fC[s] << endl; // Writes t and f values to file.
	}

	std::cout << "------------------------------------------------------------------" << std::endl;
	system("pause");
	return 0;
}
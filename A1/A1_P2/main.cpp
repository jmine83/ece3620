// 10/10/2014 - ECE 3620 - Meine, Joel
// Assignment #1 - Problem 3

#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <fstream>
using namespace std;

// vector<double> a = {2.5063,25.1125,0.6}; // Constants
// vector<double> yi = {1.5,2,-1}; // Values - Initial
vector<double> a = {869.5652173913,195.65217391304}; // Constants (P3)
vector<double> yi = {5,847.82608695652}; // Values - Initial (P3)
int O = a.size(); // Order of Differential Equation
double dt = 0.001; // Time - Step Size
double tf = 10; // Time - Final
double Dp = 10; // Data Points, No. of

// NOTE: The number of items in a and yi "must" be equal.

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

vector<vector<double>> ODE_m(vector<double> t)
{
	vector<double> y;
	vector<vector<double>> Y;
	y = yi;
	Y.push_back(y);

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
			Z[p][q] = w[p][q] + A[p][q]*dt;
	// y = Z*y
	for (int r = 0; r < t.size(); r++)
	{
		y = matrixOP_vecmul(Z, y);
		Y.push_back(y);
	}
	Y.pop_back();
	return(Y);
}

int main()
{
	std::cout << "Assignment #1 - Problem 3" << std::endl;
	std::cout << "===================================================================" << std::endl;
	std::cout << "For the differential equation..." << std::endl;
	std::cout << "(D^2 + 195.652D + 869.565)yo(t) = 869.565*f(t)" << std::endl;
	std::cout << "with initial conditions yo(0) = 5, Dyo(0) = 847.826:" << std::endl;
	std::cout << "Constant, a0 = " << a[0] << std::endl;
	std::cout << "Constant, a1 = " << a[1] << std::endl;
	std::cout << "Value - Initial, yo(0) = " << yi[0] << std::endl;
	std::cout << "Value - Initial, Dyo(0) = " << yi[1] << std::endl;
	std::cout << "Time - Step Size, dt = " << dt << std::endl;
	std::cout << "Time - Max, tf = " << tf << std::endl;
	std::cout << "-------------------------------------------------------------------" << std::endl;
	std::cout << " t     y" << std::endl;
	std::cout << "-------------------------------------------------------------------" << std::endl;
	
	double I;
	I = (tf/dt)/Dp; // Increment Length, I

	vector<double> T, sT; // Time List, T; Time List (Shortened), sT
	for (int j = 0; j <= tf / dt; j++)
	{
		T.push_back(j*dt);
	}
	for (int k = 0; k <= Dp; k++)
	{
		sT.push_back(T[k*I]);
	}
	
	vector<double> P, sP; // Print List, P; Print List (Shortened), sP
	
	vector<vector<double>> Q;
	Q = ODE_m(T);
	int f, g;
	for (f = 0; f < T.size(); f++)
		for (g = 0; g < T.size(); g++)
			if (g == 0)
				P.push_back(Q[f][g]); // Retrieves Only the Results of yo(t).
	
	for (int m = 0; m <= Dp; m++)
	{
		sP.push_back(P[m*I]);
	}
	for (int n = 0; n < sP.size(); n++)
	{
		printf(" %2.1f   %12.15f \n", sT[n], sP[n]);
	}
	ofstream outfile("A1_P2.txt");
	for (int o = 0; o < P.size(); o++)
	{
		outfile << T[o] << "  " << P[o] << endl;
	}

	std::cout << "-------------------------------------------------------------------" << std::endl;
	system("pause");
	return 0;
}
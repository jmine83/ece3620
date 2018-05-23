// 10/10/2014 - ECE 3620 - Meine, Joel
// Assignment #1 - Problem 1

#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <fstream>
using namespace std;

double a = -2.5; // Constant
double yi = 3; // Value - Initial
double dt = 0.1; // Time - Step Size
double tf = 10; // Time - Final
double Dp = 10; // Data Points, No. of

vector<double> ODE(vector<double> t)
{
	double y;
	vector<double> Y;
	y = yi;
	Y.push_back(y);
	for (int i = 0; i < t.size(); i++)
	{
		y = (1 + a*dt)*y;
		Y.push_back(y);
	}
	Y.pop_back();
	return(Y);
}

int main()
{
	std::cout << "Assignment #1 - Problem 1" << std::endl;
	std::cout << "============================" << std::endl;
	std::cout << "For the differential equation..." << std::endl;
	std::cout << "(D + 2.5)yo(t) = 0" << std::endl;
	std::cout << "with initial condition yo(0) = 3:" << std::endl;
	std::cout << "Constant, a = " << a << std::endl;
	std::cout << "Value - Initial, yo(0) = " << yi << std::endl;
	std::cout << "Time - Step Size, dt = " << dt << std::endl;
	std::cout << "Time - Max, tf = " << tf << std::endl;
	std::cout << "----------------------------" << std::endl;
	std::cout << " t     y                    " << std::endl;
	std::cout << "----------------------------" << std::endl;

	double I;
	I = (tf/dt)/Dp; // Increment Length, I

	vector<double> T,sT; // Time List, T; Time List (Shortened), sT
	for (int j = 0; j <= tf/dt; j++)
	{
		T.push_back(j*dt);
	}
	for (int k = 0; k <= I; k++)
	{
		sT.push_back(T[k*I]);
	}

	vector<double> P,sP; // Print List, P; Print List (Shortened), sP
	P = ODE(T);
	for (int m = 0; m <= I; m++)
	{
		sP.push_back(P[m*I]);
	}
	for (int n = 0; n < sP.size(); n++)
	{
		printf(" %2.1f   %20.18f \n", sT[n], sP[n]);
	}
	ofstream outfile("A1_P1.txt");
	for (int o = 0; o < P.size(); o++)
	{
		outfile << T[o] << "  " << P[o] << endl;
	}
	
	std::cout << "----------------------------" << std::endl;
	system("pause");
	return 0;
}
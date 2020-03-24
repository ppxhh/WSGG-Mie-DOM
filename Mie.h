/********************************************************
Input:  1£¬Complex refractive index: m = n - ki; 
		2, Refraction index: n(1);
        3, Absorption index: k(1);
		4, Particle radius:  D(¦Ìm);
		5, Incident wavelength: ¦Ë(¦Ìm);
********************************************************
Output: 1, Attenuation factor£ºQext(m-1);
        2, Scattering factor£ºQsca(m-1);
        3, Absorption factor£ºQabs(m-1);
*******************************************************/

#include <iostream>
#include <complex>
#include <cmath>
#include <iomanip>

using namespace std;

#define PI 3.141593
complex<double> I(0, 1);

void Mie( complex<double> m, double D, double lamuda, double *Result){
	double X = 2 * PI * D / lamuda;			// Scale parameter;
	
	int Nstop;								// nstop;
	if (X >= 0.02 && X <= 8){
		Nstop = (int) X + 4 * pow(X, 1 / 3) + 1;
	} 
	else if (X >= 8 && X <= 4200){
		Nstop = (int) X + 4.05 * pow(X, 1 / 3) + 1;
	}
	else{
		Nstop = (int) X + 4 * pow(X, 1 / 3) + 2;
	}
	
	complex<double> z = m * X;				// m * X;
	
	Nstop = Nstop + 2;				// -1 and 0 Not included£» 
	complex<double> Psi_x[Nstop], Psi_x_[Nstop];	// Psi(x) and first derivative;;
	complex<double> Psi_z[Nstop], Psi_z_[Nstop];	// Psi(z) and first derivative;;
	complex<double> Xi_x[Nstop], Xi_x_[Nstop];		// Xi(x) and first derivative;;
	
	Psi_x[0] = cos(X);		// Psi-1(x);
	Psi_x[1] = sin(X);		// Psi0(x);
	Psi_z[0] = cos(z);		// Psi-1(z);
	Psi_z[1] = sin(z);		// Psi0(z);
	Xi_x[0]  = cos(X) - I * sin(X);	//Xi-1(x);
	Xi_x[1]  = sin(X) + I * cos(X);	//Xi-1(x);
	

	for (int i = 2; i <= Nstop; ++ i){
		complex<double> absl1(2 * i - 3, 0);
		complex<double> absl2(- (i - 1), 0);
		Psi_x[i] = absl1 / X * Psi_x[i - 1] - Psi_x[i - 2];
		Psi_x_[i] = absl2 / X * Psi_x[i] + Psi_x[i - 1];
		Psi_z[i] = absl1 / z * Psi_z[i - 1] - Psi_z[i - 2];
		Psi_z_[i] = -absl2 / z * Psi_z[i] + Psi_z[i - 1];
		Xi_x[i] = absl1 / X * Xi_x[i - 1] - Xi_x[i - 2];
		Xi_x_[i] = absl2 / X * Xi_x[i] + Xi_x[i - 1];
	}
	
	complex<double> an[Nstop], bn[Nstop];
	
	for (int i = 2; i <= Nstop; ++ i){
		an[i] = (Psi_z_[i] * Psi_x[i] - m * Psi_z[i] * Xi_x_[i]) / (Psi_z_[i] * Xi_x[i] - m * Psi_z[i] * Xi_x_[i]);
		bn[i] = (m * Psi_z_[i] * Psi_x[i] - Psi_z[i] * Xi_x_[i]) / (m * Psi_z_[i] * Xi_x[i] - Psi_z[i] * Xi_x_[i]);
	}
	
	double Qext = 0.,  Qsca = 0.,  Qabs = 0.;
	
	for (int i = 2; i <= Nstop; ++ i){
		double absy1 = (2 * i - 1) * 2 / (X * X);
		Qext += absy1 * (an[i].real() + bn[i].real());
		Qsca += absy1 * (an[i].real() * an[i].real() + an[i].imag() * an[i].imag() + bn[i].real() * bn[i].real() + bn[i].imag() * bn[i].imag());
	}

	Qabs = fabs(Qext - Qsca);
	
	Result[0] = Qext;
	Result[1] = Qsca;
	Result[2] = Qabs;
	
	cout << Nstop << X << Qext << Qsca << Qabs << endl;
}


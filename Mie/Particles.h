#include <iostream>
#include <math.h>
#include <iomanip>

#include "Complex.h"

using namespace std;

#define PI 3.141593
#define N 200

class Particle
{
	private:
		double D;
		double X;
		Complex x;
		Complex m;
		Complex z;
		
		int nmx;
		int nmax;
		
		Complex Be0;
		Complex Be_0;
		Complex Be[N];
		Complex Be_[N];
		
		Complex Ha0;
		Complex Ha_0;
		Complex Ha[N];
		Complex Ha_[0];
						
	public:
		Complex Ln[N];
		Complex an[N];
		Complex bn[N]; 
		
		double Qext;
		double Qsca;
		double Qabs;
		
		Particle();
		Particle(double X, Complex b){
			x = Complex(X, 0);
			m = b;
		}
		
		void Particle_Init(const double X, Complex z);
		void Ln_Compute(Complex *Ln);
		void Bessel_Hankel();
		void Mie(Complex *Ln, const double x);
};

void Particle::Particle_Init(const double X, Complex z){
	nmx = (int)(2 + X + 4 * pow(X, 0.3333));
	if (nmx > z.Modu()){
		nmax = nmx + 16;
	}
	else{
		nmax = z.Modu() + 16;
	}
}

void Particle::Ln_Compute(Complex *Ln){
	Complex Z;
	Z = m * x;
	Ln[nmax] = Complex(0, 0);
	for (int j = nmax; j > 0; -- j){
		Complex n(j, 0);
		Complex one(1, 0);
		Ln[j - 1] = n / Z - one / (n / Z + Ln[j]);
	}
} 

void Particle::Bessel_Hankel(){
	Complex Z;
	Z = x;
	
	Complex i(0, 1);
	Complex one(1, 0);
	Complex zero(0, 0);
	
	Be0    = Z.Cos();
	Be[0]  = Z.Sin();
	Be[1]  = one / Z * Be[0] - Be0;
	Be_[0] = Be0;
	Be_[1] = zero - one / Z * Be[1] + Be[0];
	
	Ha0    = Z.Cos() + i * Z.Sin();
	Ha[0]  = Z.Sin() - i * Z.Cos();
	Ha[1]  = one / Z * Ha[0] - Ha0;
	Ha_[0] = Ha0;
	Ha_[1] = zero - one / Z * Ha[1] + Ha[0];
	
	for (int j = 2; j <= nmax; ++ j){
		Complex n(j, 0);
		Complex n2_1(2 * j - 1, 0); 
		Be[j]  = n2_1 / Z * Be[j - 1] - Be[j - 2];
		Be_[j] = zero - n / Z * Be[j] + Be[j - 1];
		
		Ha[j]  = n2_1 / Z * Ha[j - 1] - Ha[j - 2];
		Ha_[j] = zero - n / Z * Ha[j] + Ha[j - 1]; 
	}
	
	//for (int j = 0; j <= nmax; ++ j){
	//	cout << j << ":\t" << Be[j] << endl;
	//	cout << j << ":\t" << Ha[j] << endl;
	//}
}

void Particle::Mie(Complex *Ln, const double x){
	
	
	Particle::Ln_Compute(Ln);
	Particle::Bessel_Hankel();

    Complex an_up;
	Complex an_dn;
	Complex bn_up;
	Complex bn_dn;
	Complex M = m;
	Complex X = x;
	
	double qe = 0.0;
	double qs = 0.0;
	double qa = 0.0;
	
	for (int j = 1; j <= nmx + 1; ++ j){		
		Complex n(j, 0);
		an_up = (Ln[j] / M + n / X) * Be[j] - Be[j - 1];
		an_dn = (Ln[j] / M + n / X) * Ha[j] - Ha[j - 1];
		an[j] = an_up / an_dn;

		bn_up = (M * Ln[j] + n / X) * Be[j] - Be[j - 1];
		bn_dn = (M * Ln[j] + n / X) * Ha[j] - Ha[j - 1];
		bn[j] = bn_up / bn_dn;
		
		//cout << j << ":\t " << an[j] << "\t" << bn[j] << endl;
	}
	
	for (int j = 1; j <= nmx + 1; ++ j){
		Complex aan;
		Complex bbn;
		aan = an[j];
		bbn = bn[j];
		
		qe += 2 / x / x * (2 * j + 1) * (aan.Real() + bbn.Real());
		qs += 2 / x / x * (2 * j + 1) * (aan.Modu() * aan.Modu() + bbn.Modu() * bbn.Modu());
	}
	cout << qe << "\t" << qs << "\t" << qe - qs << endl;	
}

void Mie_AB(Complex *Lnmx, Complex *Lnx, Complex *an, Complex *bn, Complex M, Complex X){
	Complex i(0, 1);
	Complex one(1, 0);
	Complex zero(0, 0);
	
	Complex A[N];
	Complex B[N];
	Complex TA[N];
	Complex TB[N];
	
	Complex aan, bbn;
	double Qext = 0.0;
	double Qsca = 0.0;
	
	B[0] = zero - i;
	A[1] = one / (one + i * (X.Cos() + X * X.Sin()) / (X.Sin() - X * X.Cos()));
	B[1] = zero - one / X + one / (one / X - B[0]);
	TA[1] = (Lnmx[1] / M - Lnx[1]) / (Lnmx[1] / M - B[1]);
	TB[1] = (M * Lnmx[1] - Lnx[1]) / (M * Lnmx[1] - B[1]);
	
	for (int j = 2; j <= 38; ++ j){
		Complex n(j, 0);
		B[j] = zero - n / X + one / (n / X - B[j - 1]);
		A[j] = A[j - 1] * (B[j] + n / X) / (Lnx[j] + n / X);
		
		TA[j] = (Lnmx[j] / M - Lnx[j]) / (Lnmx[j] / M - B[j]);
		TB[j] = (M * Lnmx[j] - Lnx[j]) / (M * Lnmx[j] - B[j]);
	}
	for (int j = 1; j <= 38; ++ j){
		an[j] = A[j] * TA[j];
		bn[j] = A[j] * TB[j];
		
		cout << j << ":\t " << an[j] << "\t" << bn[j] << endl;
		
		aan = an[j];
		bbn = bn[j];
		Qext += 2 / 100.0 * (2 * j + 1) * (aan.Real() + bbn.Real());
		Qsca += 2 / 100.0 * (2 * j + 1) * (aan.Modu() * aan.Modu() + bbn.Modu() * bbn.Modu());
	}
	
	cout << Qext << "\t" << Qsca << endl;
}

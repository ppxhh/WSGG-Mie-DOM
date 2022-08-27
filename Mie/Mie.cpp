#include <iostream>
#include <cmath>
#include <iomanip>

#include "Particles.h"
 
using namespace std;
 
int main()
{
 	Complex D[N]; 
	Complex m(1.5, -0.014);
	double x = 3.095;
	Complex X(x, 0);
    Particle obj(x, m);
    Complex Z = m * X;
    obj.Particle_Init(x, Z);
	obj.Mie(D, x);
    return 0;
}

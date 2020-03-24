#include <math.h>

void WSGG_P(double XH, double XC, double L, double P, double T, 
			double *C2ij[4], double *C1ij[4], double *C0ij[4],
			double *K2i, double *K1i, double *K0i,
			double Result[2])
{
	double Tref = 2000;    // Reference Temperature;
	int I = 4;             // Number of grey gases;
	int J = 5;             // Temperature p-level of weight coefficients;
	
	double cij[I][J];      // WSGG model parameters;
	double a[I];           // 4 weight coefficients;
	double ki[I];          // 4 absorption coefficients;
	double Total_emissivity = 0;  // Total emissivity;
	
	for (int i = 0; i < I; ++ i){
		a[i] = 0;
		for (int j = 0; j < J; ++ j){
			cij[i][j] = C2ij[i][j] * P * P + C1ij[i][j]* P + C0ij[i][j];
			a[i] = a[i] + cij[i][j] * pow(Tref, I - j);
		}
		
		ki[i] = K2i[i] * P * P + K1i[i] * P + K0i[i];
		Total_emissivity = Total_emissivity + a[i] * (1 - exp(- ki[i] * (XH / P + XC / P) * L));
	}
	
	double Kaver;               // Average absorption coefficient;
	Kaver = - log(1 - Total_emissivity) / L;
	
	Result[1] = Total_emissivity;
	Result[2] = Kaver;
} 

/********************************************************
Input:  1£¬Wall temperature£ºTw1, Tw2(K); 
		2, Wall emissivity: e1, e2;
        3, Medium temperature: Tg(K);
		4, Medium absorption coefficient: ka(m-1);
		5, Medium scattering coefficient: ks(m-1);
		6, Mean beam length: L(m);
		7, Number of grids: Nx;
		S-8 Solution and Ladder format (fx = 1)
********************************************************
Output: 1, Incident radiation£ºG(Wm-2);
        2, Radiation heat flux£ºq(Wm-2);
        3, Radiation source term£º(Wm-3);
*******************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define Sigma 5.67e-8				// Stefan Boltzmann constant; 
#define M 8							// S-8 Solution;
#define PI 3.141593

double Xi[] = {0.1422555, 0.5773503, 0.8040087, 0.9795543, -0.1422555, -0.5773503, -0.8040087, -0.9795543};
double Omega[] = {2.1637144, 2.6406988, 0.7938272, 0.6849436, 2.1637144, 2.6406988, 0.7938272, 0.6849436};

void DOM_1D(double Tw1, double Tw2, double e1, double e2, double Tg, double ka, double ks, double L, int Nx, double *G, double *q, double *dq){
	double Deltx = L / Nx;		// Mesh size;
	
	int na = Nx + 2;
	double Im[na][M];
	
	// System of linear equations AI = B;
	double I[na] = {0};
	double A[na][na];
	double AA[na * na];
	double B[na] = {0};
	
	double G0[na], G1[na];		// 	Incident radiation; Iteration parameter;
	for (int i = 0; i < na; ++ i){
		G0[i] = 1.0;
	}
	
	int li, liw, lie;
	double aw, ae;
	for (int it = 1; it <= 1; ++ it){ 		// Iteration times;
		for (int i = 0; i < na; ++ i){
			G1[i] = G0[i];
		}
		
		for (int id = 0; id < M; ++ id){	// Every direction;
			for (int i = 0; i < na; ++ i){
				li = i;
				if (i == 0){
					A[li][li] = 1.0;
					B[li] = e1 * Sigma * pow(Tw1, 4) / PI;		// Wall-1;
					continue;
				}
				if (i == (na - 1)){
					A[li][li] = 1.0;
					B[li] = e2 * Sigma * pow(Tw2, 4) / PI;		// Wall-2;
					continue;
				}
				
				liw = li - 1;
				lie = li + 1;
				aw = (-Xi[id] * (-1.0) > 0.0) ? -Xi[id] * (-1.0) : 0.0;
				ae = (-Xi[id] * (1.0) > 0.0) ? -Xi[id] * (1.0) : 0.0;
				
				A[li][li] = aw + ae + (ka + ks) * Deltx;
				A[li][liw] = -aw;
				A[li][lie] = -ae;
				
				B[li] = ka * Sigma * pow(Tg, 4) / PI + ks / 4 / PI * G1[li] * Deltx;	
			}
			
			for (int i = 0; i < na; ++ i){
				for (int j = 0; j < na; ++ j){
					AA[i * na + j] = A[i][j];
				}	
			}
				
			BICGSTAB(AA, I, B);		//
			
			for (int i = 0; i < na; ++ i){
				Im[i][id] = I[i];
			}
			
			for (int i = 0; i < na; ++ i){
				double Sum_G = 0;
				for (int m = 0; m < M; ++ m){
					Sum_G = Sum_G + Im[i][m] * Omega[m];
				}
				G0[i] = Sum_G;
			}		
		}
		
		double Relative_error[na];
		long double Max_error = 1e30; 
		for (int i = 0; i < na; ++ i){
			Relative_error[i] = fabs(G0[i] - G1[i]) / (G1[i] + 1.0e-10);
			if (Relative_error[i] < Max_error){
				Max_error = Relative_error[i];
			}
		}
			
		if (Max_error < 1.0e-9){
			break;
		}
	}
	
	for (int i = 0; i < Nx; ++ i){
		G[i] = G0[i + 1];
		dq[i] = ka * (4 * PI * Sigma * pow(Tg, 4) / PI - G[i]);
	}
	
	double qw[Nx + 2] = {0.}, qe[Nx + 2] = {0.};
	for (int i = 0; i < Nx + 2; ++ i){
		for (int m = 0; m < M; ++ m){
			if (Xi[m] > 0){
				qw[i] = qw[i] + Im[i][m]* Omega[m] * Xi[m];
			}
			else{
				qe[i] = qe[i] + Im[i][m]* Omega[m] * Xi[m];
			} 
		}
	}
	
	for (int i = 0; i < Nx + 1; ++ i){
		q[i] = qw[i] + 	qe[i + 1];
	}
}



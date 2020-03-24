/********************************************************
This code can solve 3D cube radiation heat transfer; 
Input:  1£¬Wall temperature£ºTw1, Tw2, Tw3, Tw4, Tw5, Tw6(K); 
		2, Wall emissivity: e1, e2, e3, e4, e5, e6;
        3, Medium temperature: *Tg(K);
		4, Medium absorption coefficient: *ka(m-1);
		5, Medium scattering coefficient: *ks(m-1);
		6, Cube's Length, Wide, Height: Lx(m), Ly(m), Lz(m);
		7, Number of grids: Nx, Ny, Nz;
		S-6 Solution and Ladder format (fx = 1)
********************************************************
Output: 1, Incident radiation£º*G(Wm-2);
        2, Radiation heat flux£º*q(Wm-2);
        3, Radiation source term£º(Wm-3);
*******************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define Sigma 5.67e-8				// Stefan Boltzmann constant; 
#define M 6							// S-6 Solution;
#define N (M * (M + 2))				// Number of total direction 
#define PI 3.141593

void BICGSTAB(double *A, double *I, double *B){
	
}

double xi[] = {	0.1838670, 0.1838670, 0.1838670, 0.6950514, 0.6950514, 0.9656013,
				0.1838670, 0.1838670, 0.1838670, 0.6950514, 0.6950514, 0.9656013,
				0.1838670, 0.1838670, 0.1838670, 0.6950514, 0.6950514, 0.9656013,
				0.1838670, 0.1838670, 0.1838670, 0.6950514, 0.6950514, 0.9656013, 
				-0.1838670, -0.1838670, -0.1838670, -0.6950514, -0.6950514, -0.9656013,
				-0.1838670, -0.1838670, -0.1838670, -0.6950514, -0.6950514, -0.9656013,
				-0.1838670, -0.1838670, -0.1838670, -0.6950514, -0.6950514, -0.9656013,
				-0.1838670, -0.1838670, -0.1838670, -0.6950514, -0.6950514, -0.9656013};

double eta[] = {0.1838670, 0.6950514, 0.9656013, 0.183867, 0.69505140, 0.1838670,
				0.1838670, 0.6950514, 0.9656013, 0.183867, 0.69505140, 0.1838670,
				-0.1838670, -0.6950514, -0.9656013, -0.183867, -0.69505140, -0.1838670,
				-0.1838670, -0.6950514, -0.9656013, -0.183867, -0.69505140, -0.1838670,
				0.1838670, 0.6950514, 0.9656013, 0.183867, 0.69505140, 0.1838670,
				0.1838670, 0.6950514, 0.9656013, 0.183867, 0.69505140, 0.1838670,
				-0.1838670, -0.6950514, -0.9656013, -0.183867, -0.69505140, -0.1838670,
				-0.1838670, -0.6950514, -0.9656013, -0.183867, -0.69505140, -0.1838670};
				
double miu[] = {0.9656013, 0.6950514, 0.1838670, 0.6950514, 0.1838670, 0.1838670,
				-0.9656013, -0.6950514, -0.1838670, -0.6950514, -0.1838670, -0.1838670,
				0.9656013, 0.6950514, 0.1838670, 0.6950514, 0.1838670, 0.1838670,
				-0.9656013, -0.6950514, -0.1838670, -0.6950514, -0.1838670, -0.1838670,
				0.9656013, 0.6950514, 0.1838670, 0.6950514, 0.1838670, 0.1838670,
				-0.9656013, -0.6950514, -0.1838670, -0.6950514, -0.1838670, -0.1838670,
				0.9656013, 0.6950514, 0.1838670, 0.6950514, 0.1838670, 0.1838670,
				-0.9656013, -0.6950514, -0.1838670, -0.6950514, -0.1838670, -0.1838670};
				
double omega[] = {	0.1609517, 0.3626469, 0.1609517, 0.3626469, 0.3626469, 0.1609517,
					0.1609517, 0.3626469, 0.1609517, 0.3626469, 0.3626469, 0.1609517,
					0.1609517, 0.3626469, 0.1609517, 0.3626469, 0.3626469, 0.1609517,
					0.1609517, 0.3626469, 0.1609517, 0.3626469, 0.3626469, 0.1609517,
					0.1609517, 0.3626469, 0.1609517, 0.3626469, 0.3626469, 0.1609517,
					0.1609517, 0.3626469, 0.1609517, 0.3626469, 0.3626469, 0.1609517,
					0.1609517, 0.3626469, 0.1609517, 0.3626469, 0.3626469, 0.1609517,
					0.1609517, 0.3626469, 0.1609517, 0.3626469, 0.3626469, 0.1609517};
					
#define Nx 17
#define Ny 17
#define Nz 34

void DOM_3D(double *Tw, double *ew, double *Tg, double *ka, 
			double *ks, double Lx, double Ly, double Lz, 
			double *G, double *q, double *dq){
	double dx, dy, dz;
	int nax, nay, naz, na;
	
	dx = Lx / Nx;
	dy = Ly / Ny;
	dz = Lz / Nz;
	
	nax = Nx + 2;
	nay = Ny + 2;
	naz = Nz + 2;
	na  = nax * nay * naz;
	
	double *Im;
	Im = (double*)malloc(na * N * sizeof(double));
	
	double *I, *A, *B;
	I  = (double*)malloc(na * sizeof(double));
	A  = (double*)malloc(na * na * sizeof(double));
	B  = (double*)malloc(na * sizeof(double));
	
	double *G0, *G1;
	G0 = (double*)malloc(na * sizeof(double));
	G1 = (double*)malloc(na * sizeof(double));
	memset(G0, 0., na * sizeof(double));
	memset(G1, 0., na * sizeof(double));
	
	long li, liw, lie, lis, lin, lib, lit;
	double aw, ae, as, an, ab, at, ap;
	for (int it = 0; it <= 100; ++ it){
		for (int i = 0; i < na; ++ i){
			G1[i] = G0[i];
		}
		
		for (int id = 0; id < N; ++ id){
			for (int i = 0; i < nax; ++ i){
				for (int j = 0; j < nay; ++ j){
					for (int k = 0; k < naz; ++ k){
						li = i + j * nax + k * nay * nax;
						if (i == 0){
							A[li * na + li] = 1.0;
							B[li] = ew[0] * Sigma * pow(Tw[0], 4) / PI;		// x-Wall-1;
							continue;
						}
						if (i == nax - 1){
							A[li * na + li] = 1.0;
							B[li] = ew[1] * Sigma * pow(Tw[1], 4) / PI;		// x-Wall-2;
							continue;
						}
						if (j == 0){
							A[li * na + li] = 1.0;
							B[li] = ew[2] * Sigma * pow(Tw[2], 4) / PI;		// y-Wall-1;
							continue;
						}
						if (j == nay - 1){
							A[li * na + li] = 1.0;
							B[li] = ew[3] * Sigma * pow(Tw[3], 4) / PI;		// y-Wall-2;
							continue;
						}
						if (k == 0){
							A[li * na + li] = 1.0;
							B[li] = ew[4] * Sigma * pow(Tw[4], 4) / PI;		// z-Wall-1;
							continue;
						}
						if (k == naz - 1){
							A[li * na + li] = 1.0;
							B[li] = ew[5] * Sigma * pow(Tw[5], 4) / PI;		// z-Wall-2;
							continue;
						}
						
						liw = li - 1;
						lis = li - nax;
						lib = li - nax * nay;
						lie = li + 1;
						lin = li + nax;
						lit = li + nax * nay;
						
						aw  = (-xi[id] * (-1.0) > 0.0) ? -xi[id] * (-1.0) : 0.0;
						ae 	= (-xi[id] * (1.0) > 0.0) ? -xi[id] * (1.0) : 0.0;
						as  = (-eta[id] * (-1.0) > 0.0) ? -eta[id] * (-1.0) : 0.0;
						an 	= (-eta[id] * (1.0) > 0.0) ? -eta[id] * (1.0) : 0.0;
						ab  = (-miu[id] * (-1.0) > 0.0) ? -miu[id] * (-1.0) : 0.0;
						at 	= (-miu[id] * (1.0) > 0.0) ? -miu[id] * (1.0) : 0.0;
						
						A[li * na + li] = aw * dy * dz + ae * dy * dz + as * dx * dz + an * dx * dz + ab * dx * dy + at * dx * dy + (ka[li] + ks[li]) * dx * dy * dz;
						A[li * na + liw] = -aw * dy * dz;
						A[li * na + lie] = -ae * dy * dz;
						A[li * na + lis] = -as * dx * dz;
						A[li * na + lin] = -an * dx * dz;
						A[li * na + lib] = -ab * dx * dy;
						A[li * na + lit] = -at * dx * dy;
						
						B[li] = ka[li] * (Sigma * pow(Tg[li], 4) / PI) * dx * dy * dz + ks[li] / 4 / PI * G1[li] * dx *dy * dz;	
					}
				}
			}
			
			BICGSTAB(A, I, B);
			
			for (int i = 0; i < nax; ++ i){
				for (int j = 0; j < nay; ++ j){
					for (int k = 0; k < naz; ++ k){
						li = i + j * nax + k * nay * nax;
						Im[id * na + li] = I[li];
					}
				}
			}		
		}
		
		for (int i = 0; i < nax; ++ i){
			for (int j = 0; j < nay; ++ j){
				for (int k = 0; k < naz; ++ k){
					double G_Sum = 0.0;
					li = i + j * nax + k * nay * nax;
					for (int m = 0; m < M; ++ m) {
						G_Sum = G_Sum + Im[m * na + li] * omega[m];
					}
					G0[li] = G_Sum;	
				}
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
		
		for (int i = 0; i < na; ++ i){
			G[i] = G0[i + 1];
			dq[i] = ka[i] * (4 * PI * Sigma * pow(Tg[i], 4) / PI - G[i]);
		}
	}  
}

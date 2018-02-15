// NIT_inspiral - code to rapidly compute extreme mass-ratio inspirals using self-force results
// Copyright (C) 2017  Niels Warburton
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

#include "lib_Sch_GSF.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_integration.h>
#include "lib_Sch_geodesic.h"

#define LOAD_FAIL 0
#define LOAD_SUCCESS 1

int j_bar = 4;
int k_bar = 9;

int n_max = 7;

double beta = 0;

double ***a_n_jk;
double ***b_n_jk;
double ***c_n_jk;
double ***d_n_jk;

int load_matrix(double ***output, char *file_location);
double*** alloc_3D_double_array(int ni, int nj, int nk);
void free_alloced_3D_double_array(double*** array, int ni,  int nj);


/** Loads the coefficients of the mode. This function must be called before any others.*/
int lib_Sch_GSF_load_model()
{
	int model_status = MODEL_LOAD_SUCCESS;

	a_n_jk = alloc_3D_double_array(n_max+1, j_bar+1, k_bar+1);
	int a_n_jk_status = load_matrix(a_n_jk, "GSF_model_data/a_n_jk");

	b_n_jk = alloc_3D_double_array(n_max+1, j_bar+1, k_bar+1);
	int b_n_jk_status = load_matrix(b_n_jk, "GSF_model_data/b_n_jk");

	c_n_jk = alloc_3D_double_array(n_max+1, j_bar+1, k_bar+1);
	int c_n_jk_status = load_matrix(c_n_jk, "GSF_model_data/c_n_jk");

	d_n_jk = alloc_3D_double_array(n_max+1, j_bar+1, k_bar+1);
	int d_n_jk_status = load_matrix(d_n_jk, "GSF_model_data/d_n_jk");

	if(a_n_jk_status == LOAD_FAIL) {printf("Loading a_n_jk matrix failed\n"); model_status = MODEL_LOAD_FAIL;}
	if(b_n_jk_status == LOAD_FAIL) {printf("Loading b_n_jk matrix failed\n"); model_status = MODEL_LOAD_FAIL;}
	if(c_n_jk_status == LOAD_FAIL) {printf("Loading c_n_jk matrix failed\n"); model_status = MODEL_LOAD_FAIL;}
	if(d_n_jk_status == LOAD_FAIL) {printf("Loading d_n_jk matrix failed\n"); model_status = MODEL_LOAD_FAIL;}
	
	return model_status;
}

/** Frees the memory used by the model*/
void lib_Sch_GSF_free_model()
{
	free_alloced_3D_double_array(a_n_jk, n_max+1,  j_bar+1);
	free_alloced_3D_double_array(b_n_jk, n_max+1,  j_bar+1);
	free_alloced_3D_double_array(c_n_jk, n_max+1,  j_bar+1);
	free_alloced_3D_double_array(d_n_jk, n_max+1,  j_bar+1);
}

/* Returns the conservative component of the F^r for a given eccentricity e, semi-latus rectum p and orbital phase chi*/
double lib_Sch_GSF_Fr_cons(double e, double p, double chi)
{
	double A_ns_from_fit[n_max+1];

	int n,j,k;
	for(n = 0; n <= n_max; n++){
		A_ns_from_fit[n] = 0;
		for(j=0; j <= j_bar; j++){
			for(k=0; k <= k_bar; k++){
				A_ns_from_fit[n] += a_n_jk[n][j][k] * pow(e, n + 2.0*j)*pow(p,- 2.0 -k);
			}
		}
		A_ns_from_fit[n] *= pow(p, -2.0);
	}

	double Fr_cons = 0;
	Fr_cons += A_ns_from_fit[0]/2.0;
	for(n=1; n < 7; n++){
		 Fr_cons 	+= A_ns_from_fit[n]*cos(n * chi);
	}
	return Fr_cons;

}

/* Returns the dissipative component of the F^r for a given eccentricity e, semi-latus rectum p and orbital phase chi*/
double lib_Sch_GSF_Fr_diss(double e, double p, double chi)
{
	double B_ns_from_fit[n_max+1];

	int n,j,k;
	for(n = 0; n <= n_max; n++){
		B_ns_from_fit[n] = 0;
		for(j=0; j <= j_bar; j++){
			for(k=0; k <= k_bar; k++){
				B_ns_from_fit[n] += b_n_jk[n][j][k] * pow(e, n + 2.0*j)*pow(p,-4.5 -k);
			}
		}
		B_ns_from_fit[n] *= pow(p, -4.5);
	}

	double Fr_diss = 0;
	for(n=1; n < 7; n++){
		 Fr_diss 	+= B_ns_from_fit[n]*sin(n * chi);
	}
	return Fr_diss;

}

/* Returns the conservative component of the F^\phi for a given eccentricity e, semi-latus rectum p and orbital phase chi*/
double lib_Sch_GSF_Fphi_cons(double e, double p, double chi)
{
	double C_ns_from_fit[n_max+1];

	int n,j,k;
	for(n = 0; n <= n_max; n++){
		C_ns_from_fit[n] = 0;
		for(j=0; j <= j_bar; j++){
			for(k=0; k <= k_bar; k++){
				C_ns_from_fit[n] += c_n_jk[n][j][k] * pow(e, n + 2.0*j)*pow(p,-3.0 -k);
			}
		}
		C_ns_from_fit[n] *= pow(p, -3.0);
	}

	double Fphi_cons = 0;
	for(n=1; n < 7; n++){
		 Fphi_cons 	+= C_ns_from_fit[n]*sin(n * chi);
	}
	return Fphi_cons;
}

/* Returns the dissipative component of the F^\phi for a given eccentricity e, semi-latus rectum p and orbital phase chi*/
double lib_Sch_GSF_Fphi_diss(double e, double p, double chi)
{
	double D_ns_from_fit[n_max+1];

	int n,j,k;
	for(n = 0; n <= n_max; n++){
		D_ns_from_fit[n] = 0;
		for(j=0; j <= j_bar; j++){
			for(k=0; k <= k_bar; k++){
				D_ns_from_fit[n] += d_n_jk[n][j][k] * pow(e, n + 2.0*j)*pow(p,-5.5-k);
			}
		}
		D_ns_from_fit[n] *= pow(p, -5.5);
	}

	double Fphi_diss = 0;
	Fphi_diss += D_ns_from_fit[0]/2.0;
	for(n=1; n < 7; n++){
		 Fphi_diss 	+= D_ns_from_fit[n]*cos(n * chi);
	}
	return Fphi_diss;
}

/** Allocate memeory for a real 3D array
*
* @param ni the number of rows
* @param nj the number of columns
* @param nk the depth of the array block
* @return the pointer to the start of the memory block for the 2D array
*/
double*** alloc_3D_double_array(int ni, int nj, int nk)
{
	double*** array = (double***)calloc( ni, sizeof (double**));
	int i,j;
	for(i = 0; i < ni; i++){
		array[i] = (double**)calloc(nj, sizeof(double*));
		for(j = 0; j < nj; j++){
			array[i][j] = (double*)calloc(nk , sizeof(double));
		}
	}

	return array;
}

/** Free the memory associated with a real 3D array
*
* @param array the array to be freed
* @param ni the number of rows in the array
* @param nj the number of columns in the array
*/
void free_alloced_3D_double_array(double*** array, int ni,  int nj)
{
	int i, j;
	for(i = 0; i < ni; i++){
		for(j = 0; j < nj; j++){
			free(array[i][j]);
		}
		free(array[i]);
	}

}

/** Loads the matrices of model coefficients */
int load_matrix(double ***output, char *file_location)
{
	int n = 0;
	int j = 0;
	int k = 0;
	char line[6000];
	FILE* a_njk_file;

	a_njk_file = fopen(file_location,"r");
	if(a_njk_file == NULL){
		return LOAD_FAIL;
	}else{
		// Count the number of lines of data in the file
		while ( fgets(line, 6000, a_njk_file) != NULL){ 
			if(line[0] == '#') continue;
			strtok(line, " | ");
			j=0;
			k = 0;
			while(1){
				char *a_jk = strtok(NULL, " | ");
				if(a_jk == NULL) break;
				if( j > j_bar) break;
				output[n][j][k] = (double)strtod(a_jk, NULL);
				k++;
				if(k == k_bar+1){k=0; j++;}
			}
			n++;
		}
		fclose(a_njk_file);
	}
	return LOAD_SUCCESS;
}


/*
	The code below was written by Sarp Akcay.

	This file uses {p,e} as input to compute the unknown coefficient b2 so we can obtain
	
		alpha = Energy - 1.5*b2

	Since we are only interested in determining alpha here, I will omit the details. But briefly:

        b2 comes from using the method of extended homogeneous solutions.
	This involves inverting a matrix constructed from the homogeneous solutions. 	
	This inverse matrix is called Phi_inv.   
  	 
	b2 is given by the product of Phi_inv with a column vector of sources J= (0,0,J^1,J^3)^T. 
	Because of the zeros, we only need to know the right hand half (last 2 columns) of the matrix. 

	b2 itself is given by the product of the (4,3) and (4,4) elements of the matrix with J^1, J^3.

	CODE details:
	- Currently I use 2 of my existing libraries for things like E_p, L_p, dt_dchi etc., you will need to change these to your own libraries.
	- Since I used a bit of CForm stuff from mathematica, I added some new definitions.
	- As it is the code simply prints out alpha. You will need to modify it however you like to make it
	  a function that spits out alpha. I would imagine this is almost trivial.
	- I compile with standard gcc -Wall -c
*/


#define Power(x, y)	(pow((double)(x), (double)(y)))
#define Sqrt(x)		(sqrt((double)(x)))

#define Abs(x)		(fabs((double)(x)))

#define Exp(x)		(exp((double)(x)))
#define Log(x)		(log((double)(x)))
#define ln(x)		(log((double)(x)))
#define Cos(x)		(cos((double)(x)))

struct integrand_params{
	double M;
	double mu;
	double E_p;
	double L_p;
	double T_r;
};

// Integrand for b2
double integ1( double chi, void *new_params)
{
	struct integrand_params *params = (struct integrand_params*)new_params;

	double M 	= params->M;
	double mu 	= params->mu;
	double E_p	= params->E_p;
	double L_p	= params->L_p;
	double T_r	= params->T_r;

	double r = r_of_chi(chi);
	double Phi_inv_43 = (192*mu*(2.*M - 1.*r)*(4.*M - 1.*r)*Power(r,7))/
 (M*Sqrt(M_PI)*(5.12e-14*Power(M,8) + 2.56e-14*Power(M,7)*r - 
       5.12e-15*Power(M,6)*Power(r,2) + 
       4.e-15*Power(M,5)*Power(r,3) - 
       4.e-15*Power(M,4)*Power(r,4) - 
       4.e-15*Power(M,3)*Power(r,5) - 
       18432.*Power(M,2)*Power(r,6) + 18432.*M*Power(r,7) - 
       4608.*Power(r,8)));

	double Phi_inv_44 = (-192*mu*(2.*M - 1.*r)*Power(r,6)*
     (4.*Power(M,2) - 4.*M*r + Power(r,2)))/
   (M*Sqrt(M_PI)*(5.12e-14*Power(M,8) + 2.56e-14*Power(M,7)*r - 
       5.12e-15*Power(M,6)*Power(r,2) + 
       4.e-15*Power(M,5)*Power(r,3) - 
       4.e-15*Power(M,4)*Power(r,4) - 
       4.e-15*Power(M,3)*Power(r,5) - 
       18432.*Power(M,2)*Power(r,6) + 18432.*M*Power(r,7) - 
       4608.*Power(r,8)));

	double J1, J3;  // Sources for fields h^1 and h^3 for l=0,m=0.
	J1 = -16.0*sqrt(M_PI)/(E_p*pow(r,3.0)*T_r)*(2.0*E_p*E_p*r*r-(1-2/r)*r*r-(1-2/r)*L_p*L_p);
	J3= -16.0*sqrt(M_PI)/(E_p*pow(r,3.0)*T_r)*(r*r+L_p*L_p);

	
	return (Phi_inv_43*J1 + Phi_inv_44*J3)*dt_dchi(chi);
}




double lib_Sch_GSF_alpha(double e, double p)
{
	double M = 1.0;	   // The blackhole mass
	double mu = 1.0;   // small particle mass.

	set_primary_orbital_params(p, e);

	double b2, alpha;  // the unknown coefficient and time scaling factor alpha.

	double E_p, L_p, T_r; // Energy, ang. momentum and radial period.

	
	
	// Computing E_p, L_p, T_r. These are computed in the library lib_Sch_ecc.c
	E_p = orbit_energy();
	L_p = orbit_ang_mom();
	T_r = t_of_chi(2.0*M_PI);



	// Writing in the integral used to compute b2:
	// I used Maple to analytically invert the matrix Phi so that I could retain the functional dependence on r.

	// For sake of clarity I pull the factor of dt_dchi out of J1, J3 and multiply both terms by it inside the integral.


//************* doing the 4 integrals all back-to-back.*************
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (100000);
	
       
     double error;
     
     gsl_function F1;
	 F1.function = &integ1;
	
	struct integrand_params params = {M,mu,E_p,L_p,T_r};
	
	F1.params = &params;
    
	gsl_integration_qag (&F1, 0.0, M_PI, 0, 1e-10, 100000, GSL_INTEG_GAUSS61, w, &b2, &error);
    	 
	gsl_integration_workspace_free (w); 
	
// ****** alpha ************

	alpha = E_p - 1.5*b2;

	return alpha;


}



// Code for the derivative of the self-force w.r.t. p and e

/////////////////////////////////////
// Below are the derivatives w.r.t. p

/* Returns the conservative component of the F^r for a given eccentricity e, semi-latus rectum p and orbital phase chi*/
double lib_Sch_GSF_dFr_cons_dp(double e, double p, double chi)
{
	double A_ns_from_fit[n_max+1];

	int n,j,k;
	for(n = 0; n <= n_max; n++){
		A_ns_from_fit[n] = 0;
		for(j=0; j <= j_bar; j++){
			for(k=0; k <= k_bar; k++){
				A_ns_from_fit[n] += a_n_jk[n][j][k] * pow(e, n + 2.0*j)*pow(p, - 3.0 -k)*(-4.0 - k);
			}
		}
		A_ns_from_fit[n] *= pow(p, -2.0);
	}

	double Fr_cons = 0;
	Fr_cons += A_ns_from_fit[0]/2.0;
	for(n=1; n < 7; n++){
		 Fr_cons 	+= A_ns_from_fit[n]*cos(n * chi);
	}
	return Fr_cons;

}

/* Returns the dissipative component of the F^r for a given eccentricity e, semi-latus rectum p and orbital phase chi*/
double lib_Sch_GSF_dFr_diss_dp(double e, double p, double chi)
{
	double B_ns_from_fit[n_max+1];

	int n,j,k;
	for(n = 0; n <= n_max; n++){
		B_ns_from_fit[n] = 0;
		for(j=0; j <= j_bar; j++){
			for(k=0; k <= k_bar; k++){
				B_ns_from_fit[n] += b_n_jk[n][j][k] * pow(e, n + 2.0*j)*pow(p,-5.5 -k)*(-9.0 - k);
			}
		}
		B_ns_from_fit[n] *= pow(p, -4.5);
	}

	double Fr_diss = 0;
	for(n=1; n < 7; n++){
		 Fr_diss 	+= B_ns_from_fit[n]*sin(n * chi);
	}
	return Fr_diss;

}

/* Returns the conservative component of the F^\phi for a given eccentricity e, semi-latus rectum p and orbital phase chi*/
double lib_Sch_GSF_dFphi_cons_dp(double e, double p, double chi)
{
	double C_ns_from_fit[n_max+1];

	int n,j,k;
	for(n = 0; n <= n_max; n++){
		C_ns_from_fit[n] = 0;
		for(j=0; j <= j_bar; j++){
			for(k=0; k <= k_bar; k++){
				C_ns_from_fit[n] += c_n_jk[n][j][k] * pow(e, n + 2.0*j)*pow(p,-4.0 -k)*(-6.0 - k);
			}
		}
		C_ns_from_fit[n] *= pow(p, -3.0);
	}

	double Fphi_cons = 0;
	for(n=1; n < 7; n++){
		 Fphi_cons 	+= C_ns_from_fit[n]*sin(n * chi);
	}
	return Fphi_cons;
}

/* Returns the dissipative component of the F^\phi for a given eccentricity e, semi-latus rectum p and orbital phase chi*/
double lib_Sch_GSF_dFphi_diss_dp(double e, double p, double chi)
{
	double D_ns_from_fit[n_max+1];

	int n,j,k;
	for(n = 0; n <= n_max; n++){
		D_ns_from_fit[n] = 0;
		for(j=0; j <= j_bar; j++){
			for(k=0; k <= k_bar; k++){
				D_ns_from_fit[n] += d_n_jk[n][j][k] * pow(e, n + 2.0*j)*pow(p,-6.5-k)*(-11.0 - k);
			}
		}
		D_ns_from_fit[n] *= pow(p, -5.5);
	}

	double Fphi_diss = 0;
	Fphi_diss += D_ns_from_fit[0]/2.0;
	for(n=1; n < 7; n++){
		 Fphi_diss 	+= D_ns_from_fit[n]*cos(n * chi);
	}
	return Fphi_diss;
}

/////////////////////////////////////
// Below are the derivatives w.r.t. e

/* Returns the conservative component of the F^r for a given eccentricity e, semi-latus rectum p and orbital phase chi*/
double lib_Sch_GSF_dFr_cons_de(double e, double p, double chi)
{
	double A_ns_from_fit[n_max+1];

	int n,j,k;
	for(n = 0; n <= n_max; n++){
		A_ns_from_fit[n] = 0;
		for(j=0; j <= j_bar; j++){
			for(k=0; k <= k_bar; k++){
				A_ns_from_fit[n] += a_n_jk[n][j][k] * pow(e, n + 2.0*j - 1.0)*pow(p,- 2.0 -k)*(n + 2.0*j);
			}
		}
		A_ns_from_fit[n] *= pow(p, -2.0);
	}

	double Fr_cons = 0;
	Fr_cons += A_ns_from_fit[0]/2.0;
	for(n=1; n < 7; n++){
		 Fr_cons 	+= A_ns_from_fit[n]*cos(n * chi);
	}
	return Fr_cons;

}

/* Returns the dissipative component of the F^r for a given eccentricity e, semi-latus rectum p and orbital phase chi*/
double lib_Sch_GSF_dFr_diss_de(double e, double p, double chi)
{
	double B_ns_from_fit[n_max+1];

	int n,j,k;
	for(n = 0; n <= n_max; n++){
		B_ns_from_fit[n] = 0;
		for(j=0; j <= j_bar; j++){
			for(k=0; k <= k_bar; k++){
				B_ns_from_fit[n] += b_n_jk[n][j][k] * pow(e, n + 2.0*j - 1.0)*pow(p,-4.5 -k)*(n + 2.0*j);
			}
		}
		B_ns_from_fit[n] *= pow(p, -4.5);
	}

	double Fr_diss = 0;
	for(n=1; n < 7; n++){
		 Fr_diss 	+= B_ns_from_fit[n]*sin(n * chi);
	}
	return Fr_diss;

}

/* Returns the conservative component of the F^\phi for a given eccentricity e, semi-latus rectum p and orbital phase chi*/
double lib_Sch_GSF_dFphi_cons_de(double e, double p, double chi)
{
	double C_ns_from_fit[n_max+1];

	int n,j,k;
	for(n = 0; n <= n_max; n++){
		C_ns_from_fit[n] = 0;
		for(j=0; j <= j_bar; j++){
			for(k=0; k <= k_bar; k++){
				C_ns_from_fit[n] += c_n_jk[n][j][k] * pow(e, n + 2.0*j - 1.0)*pow(p,-3.0 -k)*(n + 2.0*j);
			}
		}
		C_ns_from_fit[n] *= pow(p, -3.0);
	}

	double Fphi_cons = 0;
	for(n=1; n < 7; n++){
		 Fphi_cons 	+= C_ns_from_fit[n]*sin(n * chi);
	}
	return Fphi_cons;
}

/* Returns the dissipative component of the F^\phi for a given eccentricity e, semi-latus rectum p and orbital phase chi*/
double lib_Sch_GSF_dFphi_diss_de(double e, double p, double chi)
{
	double D_ns_from_fit[n_max+1];

	int n,j,k;
	for(n = 0; n <= n_max; n++){
		D_ns_from_fit[n] = 0;
		for(j=0; j <= j_bar; j++){
			for(k=0; k <= k_bar; k++){
				D_ns_from_fit[n] += d_n_jk[n][j][k] * pow(e, n + 2.0*j - 1.0)*pow(p,-5.5-k)*(n + 2.0*j);
			}
		}
		D_ns_from_fit[n] *= pow(p, -5.5);
	}

	double Fphi_diss = 0;
	Fphi_diss += D_ns_from_fit[0]/2.0;
	for(n=1; n < 7; n++){
		 Fphi_diss 	+= D_ns_from_fit[n]*cos(n * chi);
	}
	return Fphi_diss;
}


////////////////////////////////////////
// Below are the derivatives w.r.t. v

/* Returns the conservative component of the F^r for a given eccentricity e, semi-latus rectum p and orbital phase chi*/
double lib_Sch_GSF_dFr_cons_dv(double e, double p, double chi)
{
	double A_ns_from_fit[n_max+1];

	int n,j,k;
	for(n = 0; n <= n_max; n++){
		A_ns_from_fit[n] = 0;
		for(j=0; j <= j_bar; j++){
			for(k=0; k <= k_bar; k++){
				A_ns_from_fit[n] += a_n_jk[n][j][k] * pow(e, n + 2.0*j)*pow(p,- 2.0 -k);
			}
		}
		A_ns_from_fit[n] *= pow(p, -2.0);
	}

	double Fr_cons = 0;
	//Fr_cons += A_ns_from_fit[0]/2.0;
	for(n=1; n < 7; n++){
		 Fr_cons 	+= -A_ns_from_fit[n]*n*sin(n * chi);
	}
	return Fr_cons;

}

/* Returns the dissipative component of the F^r for a given eccentricity e, semi-latus rectum p and orbital phase chi*/
double lib_Sch_GSF_dFr_diss_dv(double e, double p, double chi)
{
	double B_ns_from_fit[n_max+1];

	int n,j,k;
	for(n = 0; n <= n_max; n++){
		B_ns_from_fit[n] = 0;
		for(j=0; j <= j_bar; j++){
			for(k=0; k <= k_bar; k++){
				B_ns_from_fit[n] += b_n_jk[n][j][k] * pow(e, n + 2.0*j)*pow(p,-4.5 -k);
			}
		}
		B_ns_from_fit[n] *= pow(p, -4.5);
	}

	double Fr_diss = 0;
	for(n=1; n < 7; n++){
		 Fr_diss 	+= B_ns_from_fit[n]*n*cos(n * chi);
	}
	return Fr_diss;

}

/* Returns the conservative component of the F^\phi for a given eccentricity e, semi-latus rectum p and orbital phase chi*/
double lib_Sch_GSF_dFphi_cons_dv(double e, double p, double chi)
{
	double C_ns_from_fit[n_max+1];

	int n,j,k;
	for(n = 0; n <= n_max; n++){
		C_ns_from_fit[n] = 0;
		for(j=0; j <= j_bar; j++){
			for(k=0; k <= k_bar; k++){
				C_ns_from_fit[n] += c_n_jk[n][j][k] * pow(e, n + 2.0*j)*pow(p,-3.0 -k);
			}
		}
		C_ns_from_fit[n] *= pow(p, -3.0);
	}

	double Fphi_cons = 0;
	for(n=1; n < 7; n++){
		 Fphi_cons 	+= C_ns_from_fit[n]*n*cos(n * chi);
	}
	return Fphi_cons;
}

/* Returns the dissipative component of the F^\phi for a given eccentricity e, semi-latus rectum p and orbital phase chi*/
double lib_Sch_GSF_dFphi_diss_dv(double e, double p, double chi)
{
	double D_ns_from_fit[n_max+1];

	int n,j,k;
	for(n = 0; n <= n_max; n++){
		D_ns_from_fit[n] = 0;
		for(j=0; j <= j_bar; j++){
			for(k=0; k <= k_bar; k++){
				D_ns_from_fit[n] += d_n_jk[n][j][k] * pow(e, n + 2.0*j)*pow(p,-5.5-k);
			}
		}
		D_ns_from_fit[n] *= pow(p, -5.5);
	}

	double Fphi_diss = 0;
	//Fphi_diss += D_ns_from_fit[0]/2.0;
	for(n=1; n < 7; n++){
		 Fphi_diss 	+= -D_ns_from_fit[n]*n*sin(n * chi);
	}
	return Fphi_diss;
}

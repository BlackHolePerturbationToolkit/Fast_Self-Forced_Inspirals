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

#define MODEL_LOAD_FAIL 0
#define MODEL_LOAD_SUCCESS 1

/** Loads the coefficients of the model from the data files. This function must be called first*/
int lib_Sch_GSF_load_model();

/** Frees the memory used by the model */
void lib_Sch_GSF_free_model();

/* Returns the conservative component of the F^r for a given eccentricity e, semi-latus rectum p and orbital phase chi */
double lib_Sch_GSF_Fr_cons(double e, double p, double chi);

/* Returns the dissipative component of the F^r for a given eccentricity e, semi-latus rectum p and orbital phase chi */
double lib_Sch_GSF_Fr_diss(double e, double p, double chi);

/* Returns the conservative component of the F^\phi for a given eccentricity e, semi-latus rectum p and orbital phase chi */
double lib_Sch_GSF_Fphi_cons(double e, double p, double chi);

/* Returns the dissipative component of the F^\phi for a given eccentricity e, semi-latus rectum p and orbital phase chi */
double lib_Sch_GSF_Fphi_diss(double e, double p, double chi);

/* Computes the alpha parameter needed to convert from Lorenz gauge time to time as measured at infinity*/
double lib_Sch_GSF_alpha(double e, double p);

double lib_Sch_GSF_dFr_cons_dp(double e, double p, double chi);
double lib_Sch_GSF_dFr_diss_dp(double e, double p, double chi);

double lib_Sch_GSF_dFphi_cons_dp(double e, double p, double chi);
double lib_Sch_GSF_dFphi_diss_dp(double e, double p, double chi);


double lib_Sch_GSF_dFr_cons_de(double e, double p, double chi);
double lib_Sch_GSF_dFr_diss_de(double e, double p, double chi);

double lib_Sch_GSF_dFphi_cons_de(double e, double p, double chi);
double lib_Sch_GSF_dFphi_diss_de(double e, double p, double chi);


double lib_Sch_GSF_dFr_cons_dv(double e, double p, double chi);
double lib_Sch_GSF_dFr_diss_dv(double e, double p, double chi);

double lib_Sch_GSF_dFphi_cons_dv(double e, double p, double chi);
double lib_Sch_GSF_dFphi_diss_dv(double e, double p, double chi);

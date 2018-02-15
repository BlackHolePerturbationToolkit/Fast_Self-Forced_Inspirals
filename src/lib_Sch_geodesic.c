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

#include <math.h>
#include <gsl/gsl_integration.h>

double e;		//<! The orbital eccentricity
double p;		//<! The orbital semi-latus rectum

/*Sets the primary orbital parameters. This must be set before using any of the below functions*/
void set_primary_orbital_params(double p_new, double e_new)
{
	p = p_new;
	e = e_new;
}

/** The value of r for a given \f$\chi\f$
*
* @param chi \f$\chi\f$
* @return the radius for the supplied value of \f$\chi\f$
*/
double r_of_chi(double chi)
{
	return p/(1.0+e*cos(chi));
}


/** Calculates the orbital energy \f$ = \left( \sqrt{\frac{(p-2-2e)(p-2+2e)}{p(p-3-e^2)}}\right) \f$
*
* @return the orbital energy
*/
double orbit_energy()
{
	return sqrt( (p-2.0-2.0*e)*(p-2.0+2.0*e)/( p*(p-3.0-e*e) ) );
}


/** Calculates the orbital angular momentum \f$ = \left( \frac{p}{\sqrt{p-3-e^2}} \right) \f$
*
* @return the orbital angular momentum
*/
double orbit_ang_mom()
{
	return p/sqrt( (p-3.0-e*e) );
}


/** Calculates \f$d\chi/dt\f$ for eccentric orbit about Schwarzschild black hole
 *
 * @param chi the chi value at which to determine \f$d\chi/dt\f$
 * @return \f$d\chi/dt\f$
 */
double dt_dchi(double chi)
{
	 return p*p/( pow((1.0+e*cos(chi)), 2.0)*(p-2.0-2.0*e*cos(chi)) )*sqrt( (pow(p - 2.0, 2.0) - 4.0*e*e)/ (p - 6.0 - 2.0*e*cos(chi)) );
}


static double t_integrand(double chi, void *params)
{
	return dt_dchi(chi);
}

double t_of_chi(double chi)
{

	double result, error;

	gsl_integration_workspace * w = gsl_integration_workspace_alloc (100000);

	gsl_function F;
	F.function = &t_integrand;
	F.params = NULL;

	gsl_integration_qag (&F, 0.0, chi, 0, 1e-13, 100000, GSL_INTEG_GAUSS61, w, &result, &error); 

	gsl_integration_workspace_free (w);

	return result;
}



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

#include <gsl/gsl_sf_ellint.h>
#include <math.h>

// Definitions for using CForm'ed output from Mathematica
#define Sin(x)          (sin((double)(x)))
#define Cos(x)          (cos((double)(x)))
#define Power(x, y)     (pow((double)(x), (double)(y)))
#define Sqrt(x)         (sqrt((double)(x)))
#define Pi				M_PI

// Osculating element equations below

double dp_dchi(double p, double e, double v, double Fphi, double Fr){
	return 2.0*pow(p, 3.5)*(p-3.0-e*e)*pow(p-6.0-2.0*e*cos(v), 0.5)*(p-3.0-e*e*pow(cos(v), 2.0))/((p-6.0+2.0*e)*(p-6.0-2.0*e)*pow(1.0+e*cos(v),4.0))*Fphi - 2.0*p*p*p*e*(p-3.0-e*e)*sin(v)/((p-6.0+2.0*e)*(p-6.0-2.0*e)*pow(1.0+e*cos(v), 2.0))*Fr;								
}

double de_dchi(double p, double e, double v, double Fphi, double Fr){
	return pow(p,2.5)*(p-3.0-e*e)*((p-6.0-2.0*e*e)*((p-6.0-2.0*e*cos(v))*e*cos(v)+2.0*(p-3.0))*cos(v)+e*(p*p-10.0*p+12.0+4.0*e*e))/((p-6.0+2.0*e)*(p-6.0-2.0*e)*pow(p-6.0-2.0*e*cos(v),0.5)*pow(1.0+e*cos(v),4.0)) * Fphi + p*p*(p-3.0-e*e)*(p-6.0-2.0*e*e)*sin(v)/((p-6.0+2.0*e)*(p-6.0-2.0*e)*pow(1.0+e*cos(v), 2.0)) * Fr;								
}

double dw_dchi(double p, double e, double v, double Fphi, double Fr){
	return pow(p,2.5)*(p-3.0-e*e)*((p-6.0)*((p-6.0-2.0*e*cos(v))*e*cos(v) + 2.0*(p-3.0)) - 4.0*e*e*e*cos(v))*sin(v)/(e*(p-6.0+2.0*e)*(p-6.0-2.0*e)*pow(p-6.0-2.0*e*cos(v),0.5)*pow(1.0+e*cos(v),4.0)) * Fphi - p*p*(p-3.0-e*e)*((p-6.0)*cos(v)+2.0*e)/(e*(p-6.0+2.0*e)*(p-6.0-2.0*e)*pow(1.0+e*cos(v),2.0)) * Fr;									
}

double dphi_dchi(double p, double e, double v){
	return sqrt(p/(p-6-2*e*cos(v)));
}

double dt_dchi(double p, double e, double v){
	return (pow(p,2)*sqrt(((-2 - 2*e + p)*(-2 + 2*e + p))/(-6 + p - 2*e*cos(v))))/((-2 + p - 2*e*cos(v))*pow(1 + e*cos(v),2));
}


// Derivatives of the osculating equations w.r.t. the orbital elements and phases

double dF1p_dp(double p, double e, double v, double Fphi, double Fr, double dFphi_dp, double dFr_dp){
	
	double gpr 		= (-2*e*(3 + Power(e,2) - p)*Power(p,3)*Sin(v))/((6 + 2*e - p)*(-6 + 2*e + p)*Power(1 + e*Cos(v),2));
	double dgpr_dp	= (2*e*Power(p,2)*(324 - 12*Power(e,4) - 216*p + 39*Power(p,2) - 2*Power(p,3) + Power(e,2)*(72 - 8*p + Power(p,2)))*Sin(v))/
   (Power(-6 - 2*e + p,2)*Power(-6 + 2*e + p,2)*Power(1 + e*Cos(v),2));
	
	double gpphi	 = (2*Power(p,3.5)*(-3 - Power(e,2) + p)*Sqrt(-6 + p - 2*e*Cos(v))*(-3 + p - Power(e,2)*Power(Cos(v),2)))/((-6 - 2*e + p)*(-6 + 2*e + p)*Power(1 + e*Cos(v),4));
	double dgpphi_dp = (Power(p,2.5)*(2*p*(-6 - 2*e + p)*(-6 + 2*e + p)*(-3 - Power(e,2) + p)*(-6 + p - 2*e*Cos(v)) + 
       p*(-6 - 2*e + p)*(-6 + 2*e + p)*(-3 - Power(e,2) + p)*(-3 + p - Power(e,2)*Power(Cos(v),2)) + 
       2*p*(-6 - 2*e + p)*(-6 + 2*e + p)*(-6 + p - 2*e*Cos(v))*(-3 + p - Power(e,2)*Power(Cos(v),2)) - 
       2*p*(-6 - 2*e + p)*(-3 - Power(e,2) + p)*(-6 + p - 2*e*Cos(v))*(-3 + p - Power(e,2)*Power(Cos(v),2)) - 
       2*p*(-6 + 2*e + p)*(-3 - Power(e,2) + p)*(-6 + p - 2*e*Cos(v))*(-3 + p - Power(e,2)*Power(Cos(v),2)) + 
       7*(-6 - 2*e + p)*(-6 + 2*e + p)*(-3 - Power(e,2) + p)*(-6 + p - 2*e*Cos(v))*(-3 + p - Power(e,2)*Power(Cos(v),2))))/
   (Power(-6 - 2*e + p,2)*Power(-6 + 2*e + p,2)*Sqrt(-6 + p - 2*e*Cos(v))*Power(1 + e*Cos(v),4));
	
	return (gpr*dFr_dp + dgpr_dp*Fr) + (gpphi*dFphi_dp + dgpphi_dp*Fphi);
}

double dF1p_de(double p, double e, double v, double Fphi, double Fr, double dFphi_de, double dFr_de){
	
	double gpr 		= (-2*e*(3 + Power(e,2) - p)*Power(p,3)*Sin(v))/((6 + 2*e - p)*(-6 + 2*e + p)*Power(1 + e*Cos(v),2));
	
	double dgpr_de = (2*Power(p,3)*(-4*Power(e,4) - Power(-6 + p,2)*(-3 + p) + Power(e,2)*(120 - 40*p + 3*Power(p,2)) + 
       e*(4*Power(e,4) + Power(-6 + p,2)*(-3 + p) + Power(e,2)*(72 - 24*p + Power(p,2)))*Cos(v))*Sin(v))/
   (Power(-6 - 2*e + p,2)*Power(-6 + 2*e + p,2)*Power(1 + e*Cos(v),3));
	
	double gpphi	 = (2*Power(p,3.5)*(-3 - Power(e,2) + p)*Sqrt(-6 + p - 2*e*Cos(v))*(-3 + p - Power(e,2)*Power(Cos(v),2)))/((-6 - 2*e + p)*(-6 + 2*e + p)*Power(1 + e*Cos(v),4));
	
	double dgpphi_de = (-2*Power(p,3.5)*(2*e*(864 - 720*p + 210*Power(p,2) - 25*Power(p,3) + Power(p,4)) + 
       (-3 + p)*(4*Power(e,4)*(-23 + 4*p) + Power(-6 + p,2)*(69 - 35*p + 4*Power(p,2)) + Power(e,2)*(-216 + 72*p + 7*Power(p,2) - 2*Power(p,3)))*Cos(v) - 
       e*(4*Power(e,4)*(-9 + 5*p) + Power(-6 + p,2)*(27 - 24*p + 5*Power(p,2)) + Power(e,2)*(-936 + 480*p - 71*Power(p,2) + Power(p,3)))*Power(Cos(v),2) - 
       Power(e,2)*(4*Power(e,4)*(-7 + 2*p) + Power(e,2)*(-600 + 248*p - 25*Power(p,2)) + Power(-6 + p,2)*(21 - 13*p + 2*Power(p,2)))*Power(Cos(v),3) + 
       Power(e,3)*(12*Power(e,4) + 3*Power(-6 + p,2)*(-3 + p) + Power(e,2)*(120 - 40*p + Power(p,2)))*Power(Cos(v),4)))/
   (Power(-6 - 2*e + p,2)*Power(-6 + 2*e + p,2)*Sqrt(-6 + p - 2*e*Cos(v))*Power(1 + e*Cos(v),5));
	
	return (gpr*dFr_de + dgpr_de*Fr) + (gpphi*dFphi_de + dgpphi_de*Fphi);
	
}

double dF1p_dv(double p, double e, double v, double Fphi, double Fr, double dFphi_dv, double dFr_dv){
	
	double gpr 		= (-2*e*(3 + Power(e,2) - p)*Power(p,3)*Sin(v))/((6 + 2*e - p)*(-6 + 2*e + p)*Power(1 + e*Cos(v),2));
	double dgpr_dv 	= (e*(3 + Power(e,2) - p)*Power(p,3)*(-2*Cos(v) + e*(-3 + Cos(2*v))))/((6 + 2*e - p)*(-6 + 2*e + p)*Power(1 + e*Cos(v),3));
	
	double gpphi  	= (2*Power(p,3.5)*(-3 - Power(e,2) + p)*Sqrt(-6 + p - 2*e*Cos(v))*(-3 + p - Power(e,2)*Power(Cos(v),2)))/((-6 - 2*e + p)*(-6 + 2*e + p)*Power(1 + e*Cos(v),4));
	double dgpphi_dv = (2*e*Power(p,3.5)*(-3 - Power(e,2) + p)*(69 - 35*p + 4*Power(p,2) - e*(-9 + 5*p)*Cos(v) - Power(e,2)*(-7 + 2*p)*Power(Cos(v),2) + 3*Power(e,3)*Power(Cos(v),3))*Sin(v))/
   ((-6 - 2*e + p)*(-6 + 2*e + p)*Sqrt(-6 + p - 2*e*Cos(v))*Power(1 + e*Cos(v),5));
	
	return (gpr*dFr_dv + dgpr_dv*Fr) + (gpphi*dFphi_dv + dgpphi_dv*Fphi);
}

double dF1e_dp(double p, double e, double v, double Fphi, double Fr, double dFphi_dp, double dFr_dp){
	
	double ger = (Power(p,2)*(-6 - 2*Power(e,2) + p)*(-3 - Power(e,2) + p)*Sin(v))/((-6 - 2*e + p)*(-6 + 2*e + p)*Power(1 + e*Cos(v),2));
	
	double dger_dp = (p*(-16*Power(e,6) + 12*Power(e,4)*(4 + p) + Power(-6 + p,2)*(36 - 21*p + 2*Power(p,2)) + Power(e,2)*(720 - 360*p + 56*Power(p,2) - 3*Power(p,3)))*Sin(v))/
   (Power(-6 - 2*e + p,2)*Power(-6 + 2*e + p,2)*Power(1 + e*Cos(v),2));
	
	double gephi = (Power(p,2.5)*(-3 - Power(e,2) + p)*(e*(12 + 4*Power(e,2) - 10*p + Power(p,2)) + (-6 - 2*Power(e,2) + p)*Cos(v)*(2*(-3 + p) + e*Cos(v)*(-6 + p - 2*e*Cos(v)))))/
   ((-6 - 2*e + p)*(-6 + 2*e + p)*Sqrt(-6 + p - 2*e*Cos(v))*Power(1 + e*Cos(v),4));
	
	double dgephi_dp = (Power(p,1.5)*(2*p*(-6 - 2*e + p)*(-6 + 2*e + p)*(-3 - Power(e,2) + p)*(-6 + p - 2*e*Cos(v))*
        (2*e*(-5 + p) - 2*(9 + 2*Power(e,2) - 2*p)*Cos(v) - 2*e*(6 + Power(e,2) - p)*Power(Cos(v),2) - 2*Power(e,2)*Power(Cos(v),3)) - 
       p*(-6 - 2*e + p)*(-6 + 2*e + p)*(-3 - Power(e,2) + p)*(e*(12 + 4*Power(e,2) - 10*p + Power(p,2)) + 
          (-6 - 2*Power(e,2) + p)*Cos(v)*(2*(-3 + p) + e*Cos(v)*(-6 + p - 2*e*Cos(v)))) + 
       2*p*(-6 - 2*e + p)*(-6 + 2*e + p)*(-6 + p - 2*e*Cos(v))*(e*(12 + 4*Power(e,2) - 10*p + Power(p,2)) + 
          (-6 - 2*Power(e,2) + p)*Cos(v)*(2*(-3 + p) + e*Cos(v)*(-6 + p - 2*e*Cos(v)))) - 
       2*p*(-6 - 2*e + p)*(-3 - Power(e,2) + p)*(-6 + p - 2*e*Cos(v))*
        (e*(12 + 4*Power(e,2) - 10*p + Power(p,2)) + (-6 - 2*Power(e,2) + p)*Cos(v)*(2*(-3 + p) + e*Cos(v)*(-6 + p - 2*e*Cos(v)))) - 
       2*p*(-6 + 2*e + p)*(-3 - Power(e,2) + p)*(-6 + p - 2*e*Cos(v))*
        (e*(12 + 4*Power(e,2) - 10*p + Power(p,2)) + (-6 - 2*Power(e,2) + p)*Cos(v)*(2*(-3 + p) + e*Cos(v)*(-6 + p - 2*e*Cos(v)))) + 
       5*(-6 - 2*e + p)*(-6 + 2*e + p)*(-3 - Power(e,2) + p)*(-6 + p - 2*e*Cos(v))*
        (e*(12 + 4*Power(e,2) - 10*p + Power(p,2)) + (-6 - 2*Power(e,2) + p)*Cos(v)*(2*(-3 + p) + e*Cos(v)*(-6 + p - 2*e*Cos(v))))))/
   (2.*Power(-6 - 2*e + p,2)*Power(-6 + 2*e + p,2)*Power(-6 + p - 2*e*Cos(v),1.5)*Power(1 + e*Cos(v),4));
	
	return (ger*dFr_dp + dger_dp*Fr) + (gephi*dFphi_dp + dgephi_dp*Fphi);
}

double dF1e_de(double p, double e, double v, double Fphi, double Fr, double dFphi_de, double dFr_de){
	
	double ger = (Power(p,2)*(-6 - 2*Power(e,2) + p)*(-3 - Power(e,2) + p)*Sin(v))/((-6 - 2*e + p)*(-6 + 2*e + p)*Power(1 + e*Cos(v),2));
	
	double dger_de = (-2*Power(p,2)*(e*(-504 + 8*Power(e,4) - 4*Power(e,2)*Power(-6 + p,2) + 288*p - 52*Power(p,2) + 3*Power(p,3)) + 
       (Power(-6 + p,3)*(-3 + p) - 2*Power(e,4)*(60 - 18*p + Power(p,2)) - 8*Power(e,2)*(18 - 9*p + Power(p,2)))*Cos(v))*Sin(v))/
   (Power(-6 - 2*e + p,2)*Power(-6 + 2*e + p,2)*Power(1 + e*Cos(v),3));
	
	double gephi = (Power(p,2.5)*(-3 - Power(e,2) + p)*(e*(12 + 4*Power(e,2) - 10*p + Power(p,2)) + (-6 - 2*Power(e,2) + p)*Cos(v)*(2*(-3 + p) + e*Cos(v)*(-6 + p - 2*e*Cos(v)))))/
   ((-6 - 2*e + p)*(-6 + 2*e + p)*Sqrt(-6 + p - 2*e*Cos(v))*Power(1 + e*Cos(v),4));
	
	double dgephi_de = (Power(p,2.5)*((-6 - 2*e + p)*(-6 + 2*e + p)*(-3 - Power(e,2) + p)*(-6 + p - 2*e*Cos(v))*(1 + e*Cos(v))*
        (12 + 12*Power(e,2) - 10*p + Power(p,2) - 8*e*(-3 + p)*Cos(v) - (6 + 6*Power(e,2) - p)*(-6 + p)*Power(Cos(v),2) + 4*e*(6 + 4*Power(e,2) - p)*Power(Cos(v),3)) - 
       4*(-6 - 2*e + p)*(-6 + 2*e + p)*(-3 - Power(e,2) + p)*Cos(v)*(-6 + p - 2*e*Cos(v))*
        (e*(12 + 4*Power(e,2) - 10*p + Power(p,2)) + (-6 - 2*Power(e,2) + p)*Cos(v)*(2*(-3 + p) + e*Cos(v)*(-6 + p - 2*e*Cos(v)))) + 
       (-6 - 2*e + p)*(-6 + 2*e + p)*(-3 - Power(e,2) + p)*Cos(v)*(1 + e*Cos(v))*
        (e*(12 + 4*Power(e,2) - 10*p + Power(p,2)) + (-6 - 2*Power(e,2) + p)*Cos(v)*(2*(-3 + p) + e*Cos(v)*(-6 + p - 2*e*Cos(v)))) - 
       2*e*(-6 - 2*e + p)*(-6 + 2*e + p)*(-6 + p - 2*e*Cos(v))*(1 + e*Cos(v))*
        (e*(12 + 4*Power(e,2) - 10*p + Power(p,2)) + (-6 - 2*Power(e,2) + p)*Cos(v)*(2*(-3 + p) + e*Cos(v)*(-6 + p - 2*e*Cos(v)))) - 
       2*(-6 - 2*e + p)*(-3 - Power(e,2) + p)*(-6 + p - 2*e*Cos(v))*(1 + e*Cos(v))*
        (e*(12 + 4*Power(e,2) - 10*p + Power(p,2)) + (-6 - 2*Power(e,2) + p)*Cos(v)*(2*(-3 + p) + e*Cos(v)*(-6 + p - 2*e*Cos(v)))) + 
       2*(-6 + 2*e + p)*(-3 - Power(e,2) + p)*(-6 + p - 2*e*Cos(v))*(1 + e*Cos(v))*
        (e*(12 + 4*Power(e,2) - 10*p + Power(p,2)) + (-6 - 2*Power(e,2) + p)*Cos(v)*(2*(-3 + p) + e*Cos(v)*(-6 + p - 2*e*Cos(v))))))/
   (Power(-6 - 2*e + p,2)*Power(-6 + 2*e + p,2)*Power(-6 + p - 2*e*Cos(v),1.5)*Power(1 + e*Cos(v),5));
	
	return (ger*dFr_de + dger_de*Fr) + (gephi*dFphi_de + dgephi_de*Fphi);	
}

double dF1e_dv(double p, double e, double v, double Fphi, double Fr, double dFphi_dv, double dFr_dv){
	
	double ger 		= (Power(p,2)*(-6 - 2*Power(e,2) + p)*(-3 - Power(e,2) + p)*Sin(v))/((-6 - 2*e + p)*(-6 + 2*e + p)*Power(1 + e*Cos(v),2));
	double dger_dv 	= (Power(p,2)*(-6 - 2*Power(e,2) + p)*(-3 - Power(e,2) + p)*(Cos(v) + e*Power(Cos(v),2) + 2*e*Power(Sin(v),2)))/((-6 - 2*e + p)*(-6 + 2*e + p)*Power(1 + e*Cos(v),3));
	
	double gephi = (Power(p,2.5)*(-3 - Power(e,2) + p)*(e*(12 + 4*Power(e,2) - 10*p + Power(p,2)) + (-6 - 2*Power(e,2) + p)*Cos(v)*(2*(-3 + p) + e*Cos(v)*(-6 + p - 2*e*Cos(v)))))/
   ((-6 - 2*e + p)*(-6 + 2*e + p)*Sqrt(-6 + p - 2*e*Cos(v))*Power(1 + e*Cos(v),4));
	double dgephi_dv = (Power(p,2.5)*(-3 - Power(e,2) + p)*(864 - 1632*Power(e,2) - 694*Power(e,4) - 18*Power(e,6) - 576*p + 1516*Power(e,2)*p + 189*Power(e,4)*p + 120*Power(p,2) - 
       326*Power(e,2)*Power(p,2) - 8*Power(e,4)*Power(p,2) - 8*Power(p,3) + 20*Power(e,2)*Power(p,3) + 
       e*(42*Power(e,4)*(-8 + p) + Power(e,2)*(-1248 + 806*p - 89*Power(p,2)) + 8*(-90 + 99*p - 26*Power(p,2) + 2*Power(p,3)))*Cos(v) - 
       2*Power(e,2)*(360 + 12*Power(e,4) - 234*p + 41*Power(p,2) - 2*Power(p,3) + 4*Power(e,2)*(39 - 16*p + Power(p,2)))*Cos(2*v) - 192*Power(e,3)*Cos(3*v) - 
       64*Power(e,5)*Cos(3*v) + 74*Power(e,3)*p*Cos(3*v) + 14*Power(e,5)*p*Cos(3*v) - 7*Power(e,3)*Power(p,2)*Cos(3*v) - 18*Power(e,4)*Cos(4*v) - 
       6*Power(e,6)*Cos(4*v) + 3*Power(e,4)*p*Cos(4*v))*Sin(v))/(4.*(-6 - 2*e + p)*(-6 + 2*e + p)*Power(-6 + p - 2*e*Cos(v),1.5)*Power(1 + e*Cos(v),5));
	
	return (ger*dFr_dv + dger_dv*Fr) + (gephi*dFphi_dv + dgephi_dv*Fphi);
}

double df1v_dp(double p, double e, double v, double Fphi, double Fr, double dFphi_dp, double dFr_dp){
	
	double gwr = -((Power(p,2)*(-3 - Power(e,2) + p)*(2*e + (-6 + p)*Cos(v)))/(e*(-6 - 2*e + p)*(-6 + 2*e + p)*Power(1 + e*Cos(v),2)));
	double gwphi = (Power(p,2.5)*(-3 - Power(e,2) + p)*(-4*Power(e,3)*Cos(v) + (-6 + p)*(2*(-3 + p) + e*Cos(v)*(-6 + p - 2*e*Cos(v))))*Sin(v))/
   (e*(-6 - 2*e + p)*(-6 + 2*e + p)*Sqrt(-6 + p - 2*e*Cos(v))*Power(1 + e*Cos(v),4));
	
	double dgwr_dp = (p*(-2*e*(-216 - 48*Power(e,2) + 8*Power(e,4) + 144*p - 24*Power(p,2) + Power(p,3)) + 
       (-12*Power(e,4)*(-4 + p) - Power(-6 + p,2)*(36 - 21*p + 2*Power(p,2)) + Power(e,2)*(-288 + 72*p - 8*Power(p,2) + Power(p,3)))*Cos(v)))/
   (e*Power(-6 - 2*e + p,2)*Power(-6 + 2*e + p,2)*Power(1 + e*Cos(v),2));
	
	double dgwphi_dp = (Power(p,1.5)*(2*(-6 + p)*(4*Power(e,4)*(45 - 30*p + 4*Power(p,2)) + 3*Power(-6 + p,2)*(-45 + 45*p - 13*Power(p,2) + Power(p,3)) + 
          Power(e,2)*(-1080 + 792*p - 201*Power(p,2) + 28*Power(p,3) - 2*Power(p,4))) - 
       e*(16*Power(e,6)*(-15 + 2*p) - 8*Power(e,4)*(-540 + 261*p - 42*Power(p,2) + 2*Power(p,3)) - 
          Power(-6 + p,2)*(2160 - 2070*p + 645*Power(p,2) - 77*Power(p,3) + 3*Power(p,4)) + 
          Power(e,2)*(-10800 + 6768*p - 2052*Power(p,2) + 422*Power(p,3) - 49*Power(p,4) + 2*Power(p,5)))*Cos(v) + 
       Power(e,2)*(80*Power(e,6) - 8*Power(e,4)*(240 - 79*p + 8*Power(p,2)) - Power(-6 + p,3)*(180 - 105*p + 11*Power(p,2)) + 
          Power(e,2)*(6480 - 2736*p + 612*Power(p,2) - 98*Power(p,3) + 7*Power(p,4)))*Power(Cos(v),2) + 
       2*Power(e,3)*(4*Power(e,4)*(-30 + 7*p) + Power(-6 + p,2)*(90 - 51*p + 5*Power(p,2)) - 3*Power(e,2)*(-240 + 72*p - 10*Power(p,2) + Power(p,3)))*Power(Cos(v),3))*
     Sin(v))/(e*Power(-6 - 2*e + p,2)*Power(-6 + 2*e + p,2)*Power(-6 + p - 2*e*Cos(v),1.5)*Power(1 + e*Cos(v),4));
	
	return (gwr*dFr_dp + dgwr_dp*Fr) + (gwphi*dFphi_dp + dgwphi_dp*Fphi);
};

double df1v_de(double p, double e, double v, double Fphi, double Fr, double dFphi_de, double dFr_de){
	
	double gwr = -((Power(p,2)*(-3 - Power(e,2) + p)*(2*e + (-6 + p)*Cos(v)))/(e*(-6 - 2*e + p)*(-6 + 2*e + p)*Power(1 + e*Cos(v),2)));
	double gwphi = (Power(p,2.5)*(-3 - Power(e,2) + p)*(-4*Power(e,3)*Cos(v) + (-6 + p)*(2*(-3 + p) + e*Cos(v)*(-6 + p - 2*e*Cos(v))))*Sin(v))/
   (e*(-6 - 2*e + p)*(-6 + 2*e + p)*Sqrt(-6 + p - 2*e*Cos(v))*Power(1 + e*Cos(v),4));
	
	double dgwr_de = (Power(p,2)*(4*Power(e,3)*(48 - 16*p + Power(p,2)) + (16*Power(e,6) + Power(e,4)*(72 - 28*p) + Power(-6 + p,3)*(-3 + p) + 
          Power(e,2)*(-864 + 504*p - 90*Power(p,2) + 5*Power(p,3)))*Cos(v) + 
       e*(-6 + p)*(12*Power(e,4) + 3*Power(-6 + p,2)*(-3 + p) - Power(e,2)*(-24 + 8*p + Power(p,2)))*Power(Cos(v),2)))/
   (Power(e,2)*Power(-6 - 2*e + p,2)*Power(-6 + 2*e + p,2)*Power(1 + e*Cos(v),3));
	
	double dgwphi_de = -((Power(p,2.5)*(2*Power(-6 + p,2)*(-3 + p)*(4*Power(e,4) + Power(-6 + p,2)*(-3 + p) + Power(e,2)*(72 - 24*p + Power(p,2))) + 
         2*e*(-6 + p)*(16*Power(e,6) + (-33 + 5*p)*Power(18 - 9*p + Power(p,2),2) + 12*Power(e,4)*(9 - 8*p + Power(p,2)) + 
            Power(e,2)*(1224 - 744*p + 109*Power(p,2) + 10*Power(p,3) - 2*Power(p,4)))*Cos(v) - 
         Power(e,2)*(16*Power(e,6)*(-9 + 2*p) - 8*Power(e,4)*(-336 + 193*p - 35*Power(p,2) + 2*Power(p,3)) - 
            Power(-6 + p,3)*(-612 + 411*p - 81*Power(p,2) + 4*Power(p,3)) + Power(e,2)*(15984 - 11088*p + 2220*Power(p,2) - 10*Power(p,3) - 33*Power(p,4) + 2*Power(p,5))
            )*Power(Cos(v),2) + Power(e,3)*(80*Power(e,6) - 8*Power(e,4)*(228 - 77*p + 8*Power(p,2)) - Power(-6 + p,3)*(264 - 133*p + 15*Power(p,2)) + 
            Power(e,2)*(-5616 + 3312*p - 372*Power(p,2) - 46*Power(p,3) + 7*Power(p,4)))*Power(Cos(v),3) + 
         2*Power(e,4)*(-6 + p)*(28*Power(e,4) + 7*Power(-6 + p,2)*(-3 + p) + Power(e,2)*(24 - 8*p - 3*Power(p,2)))*Power(Cos(v),4))*Sin(v))/
     (Power(e,2)*Power(-6 - 2*e + p,2)*Power(-6 + 2*e + p,2)*Power(-6 + p - 2*e*Cos(v),1.5)*Power(1 + e*Cos(v),5)));
	
	return (gwr*dFr_de + dgwr_de*Fr) + (gwphi*dFphi_de + dgwphi_de*Fphi);
};

double df1v_dv(double p, double e, double v, double Fphi, double Fr, double dFphi_dv, double dFr_dv){
	
	double gwr = -((Power(p,2)*(-3 - Power(e,2) + p)*(2*e + (-6 + p)*Cos(v)))/(e*(-6 - 2*e + p)*(-6 + 2*e + p)*Power(1 + e*Cos(v),2)));
	double gwphi = (Power(p,2.5)*(-3 - Power(e,2) + p)*(-4*Power(e,3)*Cos(v) + (-6 + p)*(2*(-3 + p) + e*Cos(v)*(-6 + p - 2*e*Cos(v))))*Sin(v))/
   (e*(-6 - 2*e + p)*(-6 + 2*e + p)*Sqrt(-6 + p - 2*e*Cos(v))*Power(1 + e*Cos(v),4));
	
	double dgwr_dv = -(((3 + Power(e,2) - p)*Power(p,2)*(6 + 4*Power(e,2) - p + e*(-6 + p)*Cos(v))*Sin(v))/(e*(6 + 2*e - p)*(-6 + 2*e + p)*Power(1 + e*Cos(v),3)));
	double dgwphi_dv = (Power(p,2.5)*(-3 - Power(e,2) + p)*((Cos(v)*(-6 + p - 2*e*Cos(v))*(1 + e*Cos(v))*(-4*Power(e,3)*Cos(v) + (-6 + p)*(2*(-3 + p) + e*Cos(v)*(-6 + p - 2*e*Cos(v)))))/e - 
       (1 + e*Cos(v))*(6 - p + 2*e*Cos(v))*(4*Power(e,2) - Power(-6 + p,2) + 4*e*(-6 + p)*Cos(v))*Power(Sin(v),2) + 
       4*(-6 + p - 2*e*Cos(v))*(-4*Power(e,3)*Cos(v) + (-6 + p)*(2*(-3 + p) + e*Cos(v)*(-6 + p - 2*e*Cos(v))))*Power(Sin(v),2) - 
       (1 + e*Cos(v))*(-4*Power(e,3)*Cos(v) + (-6 + p)*(2*(-3 + p) + e*Cos(v)*(-6 + p - 2*e*Cos(v))))*Power(Sin(v),2)))/
   ((-6 - 2*e + p)*(-6 + 2*e + p)*Power(-6 + p - 2*e*Cos(v),1.5)*Power(1 + e*Cos(v),5));
	
	return (gwr*dFr_dv + dgwr_dv*Fr) + (gwphi*dFphi_dv + dgwphi_dv*Fphi);
};

// U0 and V0 equations, elliptic integrals use Mathematica's definitions

double EllipticK(double k){
	return gsl_sf_ellint_Kcomp(sqrt(k), GSL_PREC_DOUBLE);
}

double EllipticF(double phi, double k){
	return gsl_sf_ellint_F(phi, sqrt(k), GSL_PREC_DOUBLE) ;
}

double EllipticE(double k){
	return gsl_sf_ellint_Ecomp(sqrt(k), GSL_PREC_DOUBLE);
}

double EllipticEIncomp(double phi, double k){
	return gsl_sf_ellint_E(phi, sqrt(k), GSL_PREC_DOUBLE) ;
}

double EllipticPi(double n, double k){
	return gsl_sf_ellint_Pcomp(sqrt(k), -n, GSL_PREC_DOUBLE);
}

double EllipticPiIncomp(double n, double phi, double k){
	return gsl_sf_ellint_P(phi, sqrt(k), -n, GSL_PREC_DOUBLE);
}

double V0(double p, double e, double v){
	return (2*Sqrt(p/(-6 + 2*e + p))*(Pi*EllipticF((Pi - v)/2.,(4*e)/(-6 + 2*e + p)) + (-Pi + v)*EllipticK((4*e)/(-6 + 2*e + p))))/Pi;
}

double dV0_dp(double p, double e, double v){
	return ((p*Power(-6 + 2*e + p,2)*EllipticE((4*e)/(-6 + 2*e + p)))/(-6 - 2*e + p) + (p*Power(-6 + 2*e + p,2)*v*EllipticE((4*e)/(-6 + 2*e + p)))/((6 + 2*e - p)*Pi) - 
     (p*Power(-6 + 2*e + p,2)*EllipticEIncomp((Pi - v)/2.,(4*e)/(-6 + 2*e + p)))/(-6 - 2*e + p) + Power(-6 + 2*e + p,2)*EllipticF((Pi - v)/2.,(4*e)/(-6 + 2*e + p)) - 
     Power(-6 + 2*e + p,2)*EllipticK((4*e)/(-6 + 2*e + p)) + (Power(-6 + 2*e + p,2)*v*EllipticK((4*e)/(-6 + 2*e + p)))/Pi + 
     (2*e*p*Power(-6 + 2*e + p,1.5)*Sin(v))/((-6 - 2*e + p)*Sqrt(-6 + p - 2*e*Cos(v))))/(Sqrt(p)*Power(-6 + 2*e + p,2.5));
}

double dV0_de(double p, double e, double v){
	return (Sqrt(p)*(-2*(-6 + 2*e + p)*EllipticF((Pi - v)/2.,(4*e)/(-6 + 2*e + p)) + (2*(-6 + 2*e + p)*(Pi - v)*EllipticK((4*e)/(-6 + 2*e + p)))/Pi + 
       ((-6 + p)*(-6 + 2*e + p)*(Pi - v)*((-6 + 2*e + p)*EllipticE((4*e)/(-6 + 2*e + p)) + (6 + 2*e - p)*EllipticK((4*e)/(-6 + 2*e + p))))/(e*(6 + 2*e - p)*Pi) + 
       (-6 + p)*(-((Power(-6 + 2*e + p,2)*EllipticEIncomp((Pi - v)/2.,(4*e)/(-6 + 2*e + p)))/(e*(6 + 2*e - p))) - 
          ((-6 + 2*e + p)*EllipticF((Pi - v)/2.,(4*e)/(-6 + 2*e + p)))/e + (2*Power(-6 + 2*e + p,1.5)*Sin(v))/((6 + 2*e - p)*Sqrt(-6 + p - 2*e*Cos(v))))))/
   Power(-6 + 2*e + p,2.5);
}

double dV0_dv(double p, double e, double v){
	return -(Sqrt(p)/Sqrt(-6 + p - 2*e*Cos(v))) + (2*Sqrt(p/(-6 + 2*e + p))*EllipticK((4*e)/(-6 + 2*e + p)))/Pi;
}

double U0(double p, double e, double v){
	return -((p*Sqrt((-4*Power(e,2) + Power(-2 + p,2))*(-6 + 2*e + p))*EllipticEIncomp((Pi - v)/2.,(4*e)/(-6 + 2*e + p)))/((-1 + Power(e,2))*(-4 + p))) +
   (p*Sqrt((-4*Power(e,2) + Power(-2 + p,2))/(-6 + 2*e + p))*EllipticF((Pi - v)/2.,(4*e)/(-6 + 2*e + p)))/(-1 + Power(e,2)) +
   ((-Pi + v)*((-2*p*Sqrt((-4*Power(e,2) + Power(-2 + p,2))*(-6 + 2*e + p))*EllipticE((4*e)/(-6 + 2*e + p)))/((-1 + Power(e,2))*(-4 + p)) +
        (2*p*Sqrt((-4*Power(e,2) + Power(-2 + p,2))/(-6 + 2*e + p))*EllipticK((4*e)/(-6 + 2*e + p)))/(-1 + Power(e,2)) -
        (4*Sqrt((-4*Power(e,2) + Power(-2 + p,2))/(-6 + 2*e + p))*(8 + p - Power(p,2) + Power(e,2)*(-8 + 3*p))*EllipticPi((2*e)/(-1 + e),(4*e)/(-6 + 2*e + p)))/
         (Power(-1 + e,2)*(1 + e)*(-4 + p)) + (16*Sqrt((-4*Power(e,2) + Power(-2 + p,2))/(-6 + 2*e + p))*EllipticPi((4*e)/(-2 + 2*e + p),(4*e)/(-6 + 2*e + p)))/
         (-2 + 2*e + p)))/(2.*Pi) + (2*Sqrt((-4*Power(e,2) + Power(-2 + p,2))/(-6 + 2*e + p))*(8 + p - Power(p,2) + Power(e,2)*(-8 + 3*p))*
      EllipticPiIncomp((2*e)/(-1 + e),(-Pi + v)/2.,(4*e)/(-6 + 2*e + p)))/(Power(-1 + e,2)*(1 + e)*(-4 + p)) -
   (8*Sqrt((-4*Power(e,2) + Power(-2 + p,2))/(-6 + 2*e + p))*EllipticPiIncomp((4*e)/(-2 + 2*e + p),(-Pi + v)/2.,(4*e)/(-6 + 2*e + p)))/(-2 + 2*e + p) -
   (e*p*Sqrt((-4*Power(e,2) + Power(-2 + p,2))*(-6 + p - 2*e*Cos(v)))*Sin(v))/((-1 + Power(e,2))*(-4 + p)*(1 + e*Cos(v)));
}

double dU0_dp(double p, double e, double v){
	return (Sqrt((-4*Power(e,2) + Power(-2 + p,2))*(-6 + 2*e + p))*EllipticEIncomp(Pi/2. - v/2.,(4*e)/(-6 + 2*e + p)))/((1 - Power(e,2))*(-4 + p)) - 
   (p*Sqrt((-4*Power(e,2) + Power(-2 + p,2))*(-6 + 2*e + p))*EllipticEIncomp(Pi/2. - v/2.,(4*e)/(-6 + 2*e + p)))/((1 - Power(e,2))*Power(-4 + p,2)) + 
   (p*(-4*Power(e,2) + Power(-2 + p,2) + 2*(-2 + p)*(-6 + 2*e + p))*EllipticEIncomp(Pi/2. - v/2.,(4*e)/(-6 + 2*e + p)))/
    (2.*(1 - Power(e,2))*(-4 + p)*Sqrt((-4*Power(e,2) + Power(-2 + p,2))*(-6 + 2*e + p))) - 
   (p*Sqrt((-4*Power(e,2) + Power(-2 + p,2))*(-6 + 2*e + p))*(EllipticEIncomp(Pi/2. - v/2.,(4*e)/(-6 + 2*e + p)) - EllipticF(Pi/2. - v/2.,(4*e)/(-6 + 2*e + p))))/
    (2.*(1 - Power(e,2))*(-4 + p)*(-6 + 2*e + p)) + (Sqrt(-4*Power(e,2) + Power(-2 + p,2))*p*EllipticF(Pi/2. - v/2.,(4*e)/(-6 + 2*e + p)))/
    (2.*(1 - Power(e,2))*Power(-6 + 2*e + p,1.5)) - (Sqrt(-4*Power(e,2) + Power(-2 + p,2))*EllipticF(Pi/2. - v/2.,(4*e)/(-6 + 2*e + p)))/
    ((1 - Power(e,2))*Sqrt(-6 + 2*e + p)) - ((-2 + p)*p*EllipticF(Pi/2. - v/2.,(4*e)/(-6 + 2*e + p)))/
    ((1 - Power(e,2))*Sqrt(-4*Power(e,2) + Power(-2 + p,2))*Sqrt(-6 + 2*e + p)) + 
   ((-Pi + v)*((2*Sqrt((-4*Power(e,2) + Power(-2 + p,2))*(-6 + 2*e + p))*EllipticE((4*e)/(-6 + 2*e + p)))/((1 - Power(e,2))*(-4 + p)) - 
        (2*p*Sqrt((-4*Power(e,2) + Power(-2 + p,2))*(-6 + 2*e + p))*EllipticE((4*e)/(-6 + 2*e + p)))/((1 - Power(e,2))*Power(-4 + p,2)) + 
        (p*(-4*Power(e,2) + Power(-2 + p,2) + 2*(-2 + p)*(-6 + 2*e + p))*EllipticE((4*e)/(-6 + 2*e + p)))/
         ((1 - Power(e,2))*(-4 + p)*Sqrt((-4*Power(e,2) + Power(-2 + p,2))*(-6 + 2*e + p))) - 
        (p*Sqrt((-4*Power(e,2) + Power(-2 + p,2))*(-6 + 2*e + p))*(EllipticE((4*e)/(-6 + 2*e + p)) - EllipticK((4*e)/(-6 + 2*e + p))))/
         ((1 - Power(e,2))*(-4 + p)*(-6 + 2*e + p)) + (Sqrt(-4*Power(e,2) + Power(-2 + p,2))*p*EllipticK((4*e)/(-6 + 2*e + p)))/
         ((1 - Power(e,2))*Power(-6 + 2*e + p,1.5)) - (2*Sqrt(-4*Power(e,2) + Power(-2 + p,2))*EllipticK((4*e)/(-6 + 2*e + p)))/((1 - Power(e,2))*Sqrt(-6 + 2*e + p)) - 
        (2*(-2 + p)*p*EllipticK((4*e)/(-6 + 2*e + p)))/((1 - Power(e,2))*Sqrt(-4*Power(e,2) + Power(-2 + p,2))*Sqrt(-6 + 2*e + p)) + 
        (Sqrt(-4*Power(e,2) + Power(-2 + p,2))*p*(EllipticE((4*e)/(-6 + 2*e + p)) - (1 - (4*e)/(-6 + 2*e + p))*EllipticK((4*e)/(-6 + 2*e + p))))/
         ((1 - Power(e,2))*Power(-6 + 2*e + p,1.5)*(1 - (4*e)/(-6 + 2*e + p))) - 
        (4*Sqrt(-4*Power(e,2) + Power(-2 + p,2))*(1 + 3*Power(e,2) - 2*p)*EllipticPi((-2*e)/(1 - e),(4*e)/(-6 + 2*e + p)))/
         ((1 - e)*(1 - Power(e,2))*(-4 + p)*Sqrt(-6 + 2*e + p)) + (2*Sqrt(-4*Power(e,2) + Power(-2 + p,2))*(8*(1 - Power(e,2)) + (1 + 3*Power(e,2) - p)*p)*
           EllipticPi((-2*e)/(1 - e),(4*e)/(-6 + 2*e + p)))/((1 - e)*(1 - Power(e,2))*(-4 + p)*Power(-6 + 2*e + p,1.5)) + 
        (4*Sqrt(-4*Power(e,2) + Power(-2 + p,2))*(8*(1 - Power(e,2)) + (1 + 3*Power(e,2) - p)*p)*EllipticPi((-2*e)/(1 - e),(4*e)/(-6 + 2*e + p)))/
         ((1 - e)*(1 - Power(e,2))*Power(-4 + p,2)*Sqrt(-6 + 2*e + p)) - 
        (4*(-2 + p)*(8*(1 - Power(e,2)) + (1 + 3*Power(e,2) - p)*p)*EllipticPi((-2*e)/(1 - e),(4*e)/(-6 + 2*e + p)))/
         ((1 - e)*(1 - Power(e,2))*Sqrt(-4*Power(e,2) + Power(-2 + p,2))*(-4 + p)*Sqrt(-6 + 2*e + p)) + 
        (8*e*Sqrt(-4*Power(e,2) + Power(-2 + p,2))*(8*(1 - Power(e,2)) + (1 + 3*Power(e,2) - p)*p)*
           (EllipticE((4*e)/(-6 + 2*e + p))/(-1 + (4*e)/(-6 + 2*e + p)) + EllipticPi((-2*e)/(1 - e),(4*e)/(-6 + 2*e + p))))/
         ((1 - e)*(1 - Power(e,2))*(-4 + p)*Power(-6 + 2*e + p,2.5)*((-2*e)/(1 - e) - (4*e)/(-6 + 2*e + p))) - 
        (16*Sqrt(-4*Power(e,2) + Power(-2 + p,2))*EllipticPi((4*e)/(-2 + 2*e + p),(4*e)/(-6 + 2*e + p)))/(Sqrt(-6 + 2*e + p)*Power(-2 + 2*e + p,2)) - 
        (8*Sqrt(-4*Power(e,2) + Power(-2 + p,2))*EllipticPi((4*e)/(-2 + 2*e + p),(4*e)/(-6 + 2*e + p)))/(Power(-6 + 2*e + p,1.5)*(-2 + 2*e + p)) + 
        (16*(-2 + p)*EllipticPi((4*e)/(-2 + 2*e + p),(4*e)/(-6 + 2*e + p)))/(Sqrt(-4*Power(e,2) + Power(-2 + p,2))*Sqrt(-6 + 2*e + p)*(-2 + 2*e + p)) + 
        (16*Sqrt(-4*Power(e,2) + Power(-2 + p,2))*((-2*e*(EllipticE((4*e)/(-6 + 2*e + p))/(-1 + (4*e)/(-6 + 2*e + p)) + 
                  EllipticPi((4*e)/(-2 + 2*e + p),(4*e)/(-6 + 2*e + p))))/(Power(-6 + 2*e + p,2)*((-4*e)/(-6 + 2*e + p) + (4*e)/(-2 + 2*e + p))) - 
             (2*e*(EllipticE((4*e)/(-6 + 2*e + p)) + ((-2 + 2*e + p)*((4*e)/(-6 + 2*e + p) - (4*e)/(-2 + 2*e + p))*EllipticK((4*e)/(-6 + 2*e + p)))/(4.*e) + 
                  ((-2 + 2*e + p)*((-4*e)/(-6 + 2*e + p) + (16*Power(e,2))/Power(-2 + 2*e + p,2))*EllipticPi((4*e)/(-2 + 2*e + p),(4*e)/(-6 + 2*e + p)))/(4.*e)))/
              (Power(-2 + 2*e + p,2)*((4*e)/(-6 + 2*e + p) - (4*e)/(-2 + 2*e + p))*(-1 + (4*e)/(-2 + 2*e + p)))))/(Sqrt(-6 + 2*e + p)*(-2 + 2*e + p))))/(2.*Pi) + 
   (2*Sqrt(-4*Power(e,2) + Power(-2 + p,2))*(1 + 3*Power(e,2) - 2*p)*EllipticPiIncomp((-2*e)/(1 - e),-Pi/2. + v/2.,(4*e)/(-6 + 2*e + p)))/
    ((1 - e)*(1 - Power(e,2))*(-4 + p)*Sqrt(-6 + 2*e + p)) - (Sqrt(-4*Power(e,2) + Power(-2 + p,2))*(8*(1 - Power(e,2)) + (1 + 3*Power(e,2) - p)*p)*
      EllipticPiIncomp((-2*e)/(1 - e),-Pi/2. + v/2.,(4*e)/(-6 + 2*e + p)))/((1 - e)*(1 - Power(e,2))*(-4 + p)*Power(-6 + 2*e + p,1.5)) - 
   (2*Sqrt(-4*Power(e,2) + Power(-2 + p,2))*(8*(1 - Power(e,2)) + (1 + 3*Power(e,2) - p)*p)*EllipticPiIncomp((-2*e)/(1 - e),-Pi/2. + v/2.,(4*e)/(-6 + 2*e + p)))/
    ((1 - e)*(1 - Power(e,2))*Power(-4 + p,2)*Sqrt(-6 + 2*e + p)) + 
   (2*(-2 + p)*(8*(1 - Power(e,2)) + (1 + 3*Power(e,2) - p)*p)*EllipticPiIncomp((-2*e)/(1 - e),-Pi/2. + v/2.,(4*e)/(-6 + 2*e + p)))/
    ((1 - e)*(1 - Power(e,2))*Sqrt(-4*Power(e,2) + Power(-2 + p,2))*(-4 + p)*Sqrt(-6 + 2*e + p)) + 
   (8*Sqrt(-4*Power(e,2) + Power(-2 + p,2))*EllipticPiIncomp((4*e)/(-2 + 2*e + p),-Pi/2. + v/2.,(4*e)/(-6 + 2*e + p)))/(Sqrt(-6 + 2*e + p)*Power(-2 + 2*e + p,2)) + 
   (4*Sqrt(-4*Power(e,2) + Power(-2 + p,2))*EllipticPiIncomp((4*e)/(-2 + 2*e + p),-Pi/2. + v/2.,(4*e)/(-6 + 2*e + p)))/(Power(-6 + 2*e + p,1.5)*(-2 + 2*e + p)) - 
   (8*(-2 + p)*EllipticPiIncomp((4*e)/(-2 + 2*e + p),-Pi/2. + v/2.,(4*e)/(-6 + 2*e + p)))/(Sqrt(-4*Power(e,2) + Power(-2 + p,2))*Sqrt(-6 + 2*e + p)*(-2 + 2*e + p)) + 
   (4*e*Sqrt(-4*Power(e,2) + Power(-2 + p,2))*p*(-((-6 + 2*e + p)*EllipticEIncomp(Pi/2. - v/2.,(4*e)/(-6 + 2*e + p)))/(8.*e*(-1 + (4*e)/(-6 + 2*e + p))) - 
        ((-6 + 2*e + p)*EllipticF(Pi/2. - v/2.,(4*e)/(-6 + 2*e + p)))/(8.*e) + 
        Sin(2*(Pi/2. - v/2.))/(4.*(-1 + (4*e)/(-6 + 2*e + p))*Sqrt(1 - (4*e*Power(Cos(v/2.),2))/(-6 + 2*e + p)))))/((1 - Power(e,2))*Power(-6 + 2*e + p,2.5)) - 
   (4*e*Sqrt(-4*Power(e,2) + Power(-2 + p,2))*(8*(1 - Power(e,2)) + (1 + 3*Power(e,2) - p)*p)*
      (-(EllipticEIncomp(Pi/2. - v/2.,(4*e)/(-6 + 2*e + p))/(-1 + (4*e)/(-6 + 2*e + p))) + EllipticPiIncomp((-2*e)/(1 - e),-Pi/2. + v/2.,(4*e)/(-6 + 2*e + p)) - 
        (2*e*Sin(2*(-Pi/2. + v/2.)))/((-6 + 2*e + p)*(-1 + (4*e)/(-6 + 2*e + p))*Sqrt(1 - (4*e*Power(Cos(v/2.),2))/(-6 + 2*e + p)))))/
    ((1 - e)*(1 - Power(e,2))*(-4 + p)*Power(-6 + 2*e + p,2.5)*((-2*e)/(1 - e) - (4*e)/(-6 + 2*e + p))) - 
   (8*Sqrt(-4*Power(e,2) + Power(-2 + p,2))*((-2*e*(-(EllipticEIncomp(Pi/2. - v/2.,(4*e)/(-6 + 2*e + p))/(-1 + (4*e)/(-6 + 2*e + p))) + 
             EllipticPiIncomp((4*e)/(-2 + 2*e + p),-Pi/2. + v/2.,(4*e)/(-6 + 2*e + p)) - 
             (2*e*Sin(2*(-Pi/2. + v/2.)))/((-6 + 2*e + p)*(-1 + (4*e)/(-6 + 2*e + p))*Sqrt(1 - (4*e*Power(Cos(v/2.),2))/(-6 + 2*e + p)))))/
         (Power(-6 + 2*e + p,2)*((-4*e)/(-6 + 2*e + p) + (4*e)/(-2 + 2*e + p))) - 
        (2*e*(-EllipticEIncomp(Pi/2. - v/2.,(4*e)/(-6 + 2*e + p)) - 
             ((-2 + 2*e + p)*((4*e)/(-6 + 2*e + p) - (4*e)/(-2 + 2*e + p))*EllipticF(Pi/2. - v/2.,(4*e)/(-6 + 2*e + p)))/(4.*e) + 
             ((-2 + 2*e + p)*((-4*e)/(-6 + 2*e + p) + (16*Power(e,2))/Power(-2 + 2*e + p,2))*EllipticPiIncomp((4*e)/(-2 + 2*e + p),-Pi/2. + v/2.,(4*e)/(-6 + 2*e + p)))/
              (4.*e) - (2*e*Sqrt(1 - (4*e*Power(Cos(v/2.),2))/(-6 + 2*e + p))*Sin(2*(-Pi/2. + v/2.)))/((-2 + 2*e + p)*(1 - (4*e*Power(Cos(v/2.),2))/(-2 + 2*e + p)))))/
         (Power(-2 + 2*e + p,2)*((4*e)/(-6 + 2*e + p) - (4*e)/(-2 + 2*e + p))*(-1 + (4*e)/(-2 + 2*e + p)))))/(Sqrt(-6 + 2*e + p)*(-2 + 2*e + p)) + 
   (e*Sqrt((-4*Power(e,2) + Power(-2 + p,2))*(-6 + p - 2*e*Cos(v)))*Sin(v))/((1 - Power(e,2))*(-4 + p)*(1 + e*Cos(v))) - 
   (e*p*Sqrt((-4*Power(e,2) + Power(-2 + p,2))*(-6 + p - 2*e*Cos(v)))*Sin(v))/((1 - Power(e,2))*Power(-4 + p,2)*(1 + e*Cos(v))) + 
   (e*p*(-4*Power(e,2) + Power(-2 + p,2) + 2*(-2 + p)*(-6 + p - 2*e*Cos(v)))*Sin(v))/
    (2.*(1 - Power(e,2))*(-4 + p)*Sqrt((-4*Power(e,2) + Power(-2 + p,2))*(-6 + p - 2*e*Cos(v)))*(1 + e*Cos(v)));
}

double dU0_de(double p, double e, double v){
	return (2*e*p*Sqrt((-4*Power(e,2) + Power(-2 + p,2))*(-6 + 2*e + p))*EllipticEIncomp(Pi/2. - v/2.,(4*e)/(-6 + 2*e + p)))/(Power(1 - Power(e,2),2)*(-4 + p)) + 
   (p*(2*(-4*Power(e,2) + Power(-2 + p,2)) - 8*e*(-6 + 2*e + p))*EllipticEIncomp(Pi/2. - v/2.,(4*e)/(-6 + 2*e + p)))/
    (2.*(1 - Power(e,2))*(-4 + p)*Sqrt((-4*Power(e,2) + Power(-2 + p,2))*(-6 + 2*e + p))) + 
   (p*(-6 + 2*e + p)*Sqrt((-4*Power(e,2) + Power(-2 + p,2))*(-6 + 2*e + p))*((-8*e)/Power(-6 + 2*e + p,2) + 4/(-6 + 2*e + p))*
      (EllipticEIncomp(Pi/2. - v/2.,(4*e)/(-6 + 2*e + p)) - EllipticF(Pi/2. - v/2.,(4*e)/(-6 + 2*e + p))))/(8.*e*(1 - Power(e,2))*(-4 + p)) + 
   (Sqrt(-4*Power(e,2) + Power(-2 + p,2))*p*EllipticF(Pi/2. - v/2.,(4*e)/(-6 + 2*e + p)))/((1 - Power(e,2))*Power(-6 + 2*e + p,1.5)) + 
   (4*e*p*EllipticF(Pi/2. - v/2.,(4*e)/(-6 + 2*e + p)))/((1 - Power(e,2))*Sqrt(-4*Power(e,2) + Power(-2 + p,2))*Sqrt(-6 + 2*e + p)) - 
   (2*e*Sqrt(-4*Power(e,2) + Power(-2 + p,2))*p*EllipticF(Pi/2. - v/2.,(4*e)/(-6 + 2*e + p)))/(Power(1 - Power(e,2),2)*Sqrt(-6 + 2*e + p)) + 
   ((-Pi + v)*((4*e*p*Sqrt((-4*Power(e,2) + Power(-2 + p,2))*(-6 + 2*e + p))*EllipticE((4*e)/(-6 + 2*e + p)))/(Power(1 - Power(e,2),2)*(-4 + p)) + 
        (p*(2*(-4*Power(e,2) + Power(-2 + p,2)) - 8*e*(-6 + 2*e + p))*EllipticE((4*e)/(-6 + 2*e + p)))/
         ((1 - Power(e,2))*(-4 + p)*Sqrt((-4*Power(e,2) + Power(-2 + p,2))*(-6 + 2*e + p))) + 
        (p*(-6 + 2*e + p)*Sqrt((-4*Power(e,2) + Power(-2 + p,2))*(-6 + 2*e + p))*((-8*e)/Power(-6 + 2*e + p,2) + 4/(-6 + 2*e + p))*
           (EllipticE((4*e)/(-6 + 2*e + p)) - EllipticK((4*e)/(-6 + 2*e + p))))/(4.*e*(1 - Power(e,2))*(-4 + p)) + 
        (2*Sqrt(-4*Power(e,2) + Power(-2 + p,2))*p*EllipticK((4*e)/(-6 + 2*e + p)))/((1 - Power(e,2))*Power(-6 + 2*e + p,1.5)) + 
        (8*e*p*EllipticK((4*e)/(-6 + 2*e + p)))/((1 - Power(e,2))*Sqrt(-4*Power(e,2) + Power(-2 + p,2))*Sqrt(-6 + 2*e + p)) - 
        (4*e*Sqrt(-4*Power(e,2) + Power(-2 + p,2))*p*EllipticK((4*e)/(-6 + 2*e + p)))/(Power(1 - Power(e,2),2)*Sqrt(-6 + 2*e + p)) - 
        (Sqrt(-4*Power(e,2) + Power(-2 + p,2))*p*Sqrt(-6 + 2*e + p)*((-8*e)/Power(-6 + 2*e + p,2) + 4/(-6 + 2*e + p))*
           (EllipticE((4*e)/(-6 + 2*e + p)) - (1 - (4*e)/(-6 + 2*e + p))*EllipticK((4*e)/(-6 + 2*e + p))))/(4.*e*(1 - Power(e,2))*(1 - (4*e)/(-6 + 2*e + p))) - 
        (4*Sqrt(-4*Power(e,2) + Power(-2 + p,2))*(-16*e + 6*e*p)*EllipticPi((-2*e)/(1 - e),(4*e)/(-6 + 2*e + p)))/
         ((1 - e)*(1 - Power(e,2))*(-4 + p)*Sqrt(-6 + 2*e + p)) + (4*Sqrt(-4*Power(e,2) + Power(-2 + p,2))*(8*(1 - Power(e,2)) + (1 + 3*Power(e,2) - p)*p)*
           EllipticPi((-2*e)/(1 - e),(4*e)/(-6 + 2*e + p)))/((1 - e)*(1 - Power(e,2))*(-4 + p)*Power(-6 + 2*e + p,1.5)) + 
        (16*e*(8*(1 - Power(e,2)) + (1 + 3*Power(e,2) - p)*p)*EllipticPi((-2*e)/(1 - e),(4*e)/(-6 + 2*e + p)))/
         ((1 - e)*(1 - Power(e,2))*Sqrt(-4*Power(e,2) + Power(-2 + p,2))*(-4 + p)*Sqrt(-6 + 2*e + p)) - 
        (8*e*Sqrt(-4*Power(e,2) + Power(-2 + p,2))*(8*(1 - Power(e,2)) + (1 + 3*Power(e,2) - p)*p)*EllipticPi((-2*e)/(1 - e),(4*e)/(-6 + 2*e + p)))/
         ((1 - e)*Power(1 - Power(e,2),2)*(-4 + p)*Sqrt(-6 + 2*e + p)) - 
        (4*Sqrt(-4*Power(e,2) + Power(-2 + p,2))*(8*(1 - Power(e,2)) + (1 + 3*Power(e,2) - p)*p)*EllipticPi((-2*e)/(1 - e),(4*e)/(-6 + 2*e + p)))/
         (Power(1 - e,2)*(1 - Power(e,2))*(-4 + p)*Sqrt(-6 + 2*e + p)) - 
        (4*Sqrt(-4*Power(e,2) + Power(-2 + p,2))*(8*(1 - Power(e,2)) + (1 + 3*Power(e,2) - p)*p)*
           ((((-8*e)/Power(-6 + 2*e + p,2) + 4/(-6 + 2*e + p))*(EllipticE((4*e)/(-6 + 2*e + p))/(-1 + (4*e)/(-6 + 2*e + p)) + 
                  EllipticPi((-2*e)/(1 - e),(4*e)/(-6 + 2*e + p))))/(2.*((-2*e)/(1 - e) - (4*e)/(-6 + 2*e + p))) + 
             ((-2/(1 - e) - (2*e)/Power(1 - e,2))*(EllipticE((4*e)/(-6 + 2*e + p)) - 
                  ((1 - e)*((2*e)/(1 - e) + (4*e)/(-6 + 2*e + p))*EllipticK((4*e)/(-6 + 2*e + p)))/(2.*e) - 
                  ((1 - e)*((4*Power(e,2))/Power(1 - e,2) - (4*e)/(-6 + 2*e + p))*EllipticPi((-2*e)/(1 - e),(4*e)/(-6 + 2*e + p)))/(2.*e)))/
              (2.*(-1 - (2*e)/(1 - e))*((2*e)/(1 - e) + (4*e)/(-6 + 2*e + p)))))/((1 - e)*(1 - Power(e,2))*(-4 + p)*Sqrt(-6 + 2*e + p)) - 
        (32*Sqrt(-4*Power(e,2) + Power(-2 + p,2))*EllipticPi((4*e)/(-2 + 2*e + p),(4*e)/(-6 + 2*e + p)))/(Sqrt(-6 + 2*e + p)*Power(-2 + 2*e + p,2)) - 
        (16*Sqrt(-4*Power(e,2) + Power(-2 + p,2))*EllipticPi((4*e)/(-2 + 2*e + p),(4*e)/(-6 + 2*e + p)))/(Power(-6 + 2*e + p,1.5)*(-2 + 2*e + p)) - 
        (64*e*EllipticPi((4*e)/(-2 + 2*e + p),(4*e)/(-6 + 2*e + p)))/(Sqrt(-4*Power(e,2) + Power(-2 + p,2))*Sqrt(-6 + 2*e + p)*(-2 + 2*e + p)) + 
        (16*Sqrt(-4*Power(e,2) + Power(-2 + p,2))*((((-8*e)/Power(-6 + 2*e + p,2) + 4/(-6 + 2*e + p))*
                (EllipticE((4*e)/(-6 + 2*e + p))/(-1 + (4*e)/(-6 + 2*e + p)) + EllipticPi((4*e)/(-2 + 2*e + p),(4*e)/(-6 + 2*e + p))))/
              (2.*((-4*e)/(-6 + 2*e + p) + (4*e)/(-2 + 2*e + p))) + 
             (((-8*e)/Power(-2 + 2*e + p,2) + 4/(-2 + 2*e + p))*(EllipticE((4*e)/(-6 + 2*e + p)) + 
                  ((-2 + 2*e + p)*((4*e)/(-6 + 2*e + p) - (4*e)/(-2 + 2*e + p))*EllipticK((4*e)/(-6 + 2*e + p)))/(4.*e) + 
                  ((-2 + 2*e + p)*((-4*e)/(-6 + 2*e + p) + (16*Power(e,2))/Power(-2 + 2*e + p,2))*EllipticPi((4*e)/(-2 + 2*e + p),(4*e)/(-6 + 2*e + p)))/(4.*e)))/
              (2.*((4*e)/(-6 + 2*e + p) - (4*e)/(-2 + 2*e + p))*(-1 + (4*e)/(-2 + 2*e + p)))))/(Sqrt(-6 + 2*e + p)*(-2 + 2*e + p))))/(2.*Pi) + 
   (2*Sqrt(-4*Power(e,2) + Power(-2 + p,2))*(-16*e + 6*e*p)*EllipticPiIncomp((-2*e)/(1 - e),-Pi/2. + v/2.,(4*e)/(-6 + 2*e + p)))/
    ((1 - e)*(1 - Power(e,2))*(-4 + p)*Sqrt(-6 + 2*e + p)) - (2*Sqrt(-4*Power(e,2) + Power(-2 + p,2))*(8*(1 - Power(e,2)) + (1 + 3*Power(e,2) - p)*p)*
      EllipticPiIncomp((-2*e)/(1 - e),-Pi/2. + v/2.,(4*e)/(-6 + 2*e + p)))/((1 - e)*(1 - Power(e,2))*(-4 + p)*Power(-6 + 2*e + p,1.5)) - 
   (8*e*(8*(1 - Power(e,2)) + (1 + 3*Power(e,2) - p)*p)*EllipticPiIncomp((-2*e)/(1 - e),-Pi/2. + v/2.,(4*e)/(-6 + 2*e + p)))/
    ((1 - e)*(1 - Power(e,2))*Sqrt(-4*Power(e,2) + Power(-2 + p,2))*(-4 + p)*Sqrt(-6 + 2*e + p)) + 
   (4*e*Sqrt(-4*Power(e,2) + Power(-2 + p,2))*(8*(1 - Power(e,2)) + (1 + 3*Power(e,2) - p)*p)*EllipticPiIncomp((-2*e)/(1 - e),-Pi/2. + v/2.,(4*e)/(-6 + 2*e + p)))/
    ((1 - e)*Power(1 - Power(e,2),2)*(-4 + p)*Sqrt(-6 + 2*e + p)) + 
   (2*Sqrt(-4*Power(e,2) + Power(-2 + p,2))*(8*(1 - Power(e,2)) + (1 + 3*Power(e,2) - p)*p)*EllipticPiIncomp((-2*e)/(1 - e),-Pi/2. + v/2.,(4*e)/(-6 + 2*e + p)))/
    (Power(1 - e,2)*(1 - Power(e,2))*(-4 + p)*Sqrt(-6 + 2*e + p)) + 
   (16*Sqrt(-4*Power(e,2) + Power(-2 + p,2))*EllipticPiIncomp((4*e)/(-2 + 2*e + p),-Pi/2. + v/2.,(4*e)/(-6 + 2*e + p)))/(Sqrt(-6 + 2*e + p)*Power(-2 + 2*e + p,2)) + 
   (8*Sqrt(-4*Power(e,2) + Power(-2 + p,2))*EllipticPiIncomp((4*e)/(-2 + 2*e + p),-Pi/2. + v/2.,(4*e)/(-6 + 2*e + p)))/(Power(-6 + 2*e + p,1.5)*(-2 + 2*e + p)) + 
   (32*e*EllipticPiIncomp((4*e)/(-2 + 2*e + p),-Pi/2. + v/2.,(4*e)/(-6 + 2*e + p)))/(Sqrt(-4*Power(e,2) + Power(-2 + p,2))*Sqrt(-6 + 2*e + p)*(-2 + 2*e + p)) - 
   (Sqrt(-4*Power(e,2) + Power(-2 + p,2))*p*((-8*e)/Power(-6 + 2*e + p,2) + 4/(-6 + 2*e + p))*
      (-((-6 + 2*e + p)*EllipticEIncomp(Pi/2. - v/2.,(4*e)/(-6 + 2*e + p)))/(8.*e*(-1 + (4*e)/(-6 + 2*e + p))) - 
        ((-6 + 2*e + p)*EllipticF(Pi/2. - v/2.,(4*e)/(-6 + 2*e + p)))/(8.*e) + 
        Sin(2*(Pi/2. - v/2.))/(4.*(-1 + (4*e)/(-6 + 2*e + p))*Sqrt(1 - (4*e*Power(Cos(v/2.),2))/(-6 + 2*e + p)))))/((1 - Power(e,2))*Sqrt(-6 + 2*e + p)) + 
   (2*Sqrt(-4*Power(e,2) + Power(-2 + p,2))*(8*(1 - Power(e,2)) + (1 + 3*Power(e,2) - p)*p)*
      ((((-8*e)/Power(-6 + 2*e + p,2) + 4/(-6 + 2*e + p))*(-(EllipticEIncomp(Pi/2. - v/2.,(4*e)/(-6 + 2*e + p))/(-1 + (4*e)/(-6 + 2*e + p))) + 
             EllipticPiIncomp((-2*e)/(1 - e),-Pi/2. + v/2.,(4*e)/(-6 + 2*e + p)) - 
             (2*e*Sin(2*(-Pi/2. + v/2.)))/((-6 + 2*e + p)*(-1 + (4*e)/(-6 + 2*e + p))*Sqrt(1 - (4*e*Power(Cos(v/2.),2))/(-6 + 2*e + p)))))/
         (2.*((-2*e)/(1 - e) - (4*e)/(-6 + 2*e + p))) + ((-2/(1 - e) - (2*e)/Power(1 - e,2))*
           (-EllipticEIncomp(Pi/2. - v/2.,(4*e)/(-6 + 2*e + p)) + ((1 - e)*((2*e)/(1 - e) + (4*e)/(-6 + 2*e + p))*EllipticF(Pi/2. - v/2.,(4*e)/(-6 + 2*e + p)))/
              (2.*e) - ((1 - e)*((4*Power(e,2))/Power(1 - e,2) - (4*e)/(-6 + 2*e + p))*EllipticPiIncomp((-2*e)/(1 - e),-Pi/2. + v/2.,(4*e)/(-6 + 2*e + p)))/(2.*e) + 
             (e*Sqrt(1 - (4*e*Power(Cos(v/2.),2))/(-6 + 2*e + p))*Sin(2*(-Pi/2. + v/2.)))/((1 - e)*(1 + (2*e*Power(Cos(v/2.),2))/(1 - e)))))/
         (2.*(-1 - (2*e)/(1 - e))*((2*e)/(1 - e) + (4*e)/(-6 + 2*e + p)))))/((1 - e)*(1 - Power(e,2))*(-4 + p)*Sqrt(-6 + 2*e + p)) - 
   (8*Sqrt(-4*Power(e,2) + Power(-2 + p,2))*((((-8*e)/Power(-6 + 2*e + p,2) + 4/(-6 + 2*e + p))*
           (-(EllipticEIncomp(Pi/2. - v/2.,(4*e)/(-6 + 2*e + p))/(-1 + (4*e)/(-6 + 2*e + p))) + 
             EllipticPiIncomp((4*e)/(-2 + 2*e + p),-Pi/2. + v/2.,(4*e)/(-6 + 2*e + p)) - 
             (2*e*Sin(2*(-Pi/2. + v/2.)))/((-6 + 2*e + p)*(-1 + (4*e)/(-6 + 2*e + p))*Sqrt(1 - (4*e*Power(Cos(v/2.),2))/(-6 + 2*e + p)))))/
         (2.*((-4*e)/(-6 + 2*e + p) + (4*e)/(-2 + 2*e + p))) + (((-8*e)/Power(-2 + 2*e + p,2) + 4/(-2 + 2*e + p))*
           (-EllipticEIncomp(Pi/2. - v/2.,(4*e)/(-6 + 2*e + p)) - ((-2 + 2*e + p)*((4*e)/(-6 + 2*e + p) - (4*e)/(-2 + 2*e + p))*
                EllipticF(Pi/2. - v/2.,(4*e)/(-6 + 2*e + p)))/(4.*e) + 
             ((-2 + 2*e + p)*((-4*e)/(-6 + 2*e + p) + (16*Power(e,2))/Power(-2 + 2*e + p,2))*EllipticPiIncomp((4*e)/(-2 + 2*e + p),-Pi/2. + v/2.,(4*e)/(-6 + 2*e + p)))/
              (4.*e) - (2*e*Sqrt(1 - (4*e*Power(Cos(v/2.),2))/(-6 + 2*e + p))*Sin(2*(-Pi/2. + v/2.)))/((-2 + 2*e + p)*(1 - (4*e*Power(Cos(v/2.),2))/(-2 + 2*e + p)))))/
         (2.*((4*e)/(-6 + 2*e + p) - (4*e)/(-2 + 2*e + p))*(-1 + (4*e)/(-2 + 2*e + p)))))/(Sqrt(-6 + 2*e + p)*(-2 + 2*e + p)) - 
   (e*p*Cos(v)*Sqrt((-4*Power(e,2) + Power(-2 + p,2))*(-6 + p - 2*e*Cos(v)))*Sin(v))/((1 - Power(e,2))*(-4 + p)*Power(1 + e*Cos(v),2)) + 
   (2*Power(e,2)*p*Sqrt((-4*Power(e,2) + Power(-2 + p,2))*(-6 + p - 2*e*Cos(v)))*Sin(v))/(Power(1 - Power(e,2),2)*(-4 + p)*(1 + e*Cos(v))) + 
   (p*Sqrt((-4*Power(e,2) + Power(-2 + p,2))*(-6 + p - 2*e*Cos(v)))*Sin(v))/((1 - Power(e,2))*(-4 + p)*(1 + e*Cos(v))) + 
   (e*p*(-2*(-4*Power(e,2) + Power(-2 + p,2))*Cos(v) - 8*e*(-6 + p - 2*e*Cos(v)))*Sin(v))/
    (2.*(1 - Power(e,2))*(-4 + p)*Sqrt((-4*Power(e,2) + Power(-2 + p,2))*(-6 + p - 2*e*Cos(v)))*(1 + e*Cos(v)));
}

double dU0_dv(double p, double e, double v){
	return (Sqrt(-4*Power(e,2) + Power(-2 + p,2))*p)/(2.*(1 - Power(e,2))*Sqrt(-6 + 2*e + p)*Sqrt(1 - (4*e*Power(Cos(v/2.),2))/(-6 + 2*e + p))) + 
   (Sqrt(-4*Power(e,2) + Power(-2 + p,2))*(8*(1 - Power(e,2)) + (1 + 3*Power(e,2) - p)*p))/
    ((1 - e)*(1 - Power(e,2))*(-4 + p)*Sqrt(-6 + 2*e + p)*(1 + (2*e*Power(Cos(v/2.),2))/(1 - e))*Sqrt(1 - (4*e*Power(Cos(v/2.),2))/(-6 + 2*e + p))) - 
   (p*Sqrt((-4*Power(e,2) + Power(-2 + p,2))*(-6 + 2*e + p))*Sqrt(1 - (4*e*Power(Cos(v/2.),2))/(-6 + 2*e + p)))/(2.*(1 - Power(e,2))*(-4 + p)) - 
   (4*Sqrt(-4*Power(e,2) + Power(-2 + p,2)))/
    (Sqrt(-6 + 2*e + p)*(-2 + 2*e + p)*Sqrt(1 - (4*e*Power(Cos(v/2.),2))/(-6 + 2*e + p))*(1 - (4*e*Power(Cos(v/2.),2))/(-2 + 2*e + p))) + 
   (e*p*Cos(v)*Sqrt((-4*Power(e,2) + Power(-2 + p,2))*(-6 + p - 2*e*Cos(v))))/((1 - Power(e,2))*(-4 + p)*(1 + e*Cos(v))) + 
   ((2*p*Sqrt((-4*Power(e,2) + Power(-2 + p,2))*(-6 + 2*e + p))*EllipticE((4*e)/(-6 + 2*e + p)))/((1 - Power(e,2))*(-4 + p)) - 
      (2*Sqrt(-4*Power(e,2) + Power(-2 + p,2))*p*EllipticK((4*e)/(-6 + 2*e + p)))/((1 - Power(e,2))*Sqrt(-6 + 2*e + p)) - 
      (4*Sqrt(-4*Power(e,2) + Power(-2 + p,2))*(8*(1 - Power(e,2)) + (1 + 3*Power(e,2) - p)*p)*EllipticPi((-2*e)/(1 - e),(4*e)/(-6 + 2*e + p)))/
       ((1 - e)*(1 - Power(e,2))*(-4 + p)*Sqrt(-6 + 2*e + p)) + (16*Sqrt(-4*Power(e,2) + Power(-2 + p,2))*EllipticPi((4*e)/(-2 + 2*e + p),(4*e)/(-6 + 2*e + p)))/
       (Sqrt(-6 + 2*e + p)*(-2 + 2*e + p)))/(2.*Pi) + (Power(e,2)*p*Sqrt((-4*Power(e,2) + Power(-2 + p,2))*(-6 + p - 2*e*Cos(v)))*Power(Sin(v),2))/
    ((1 - Power(e,2))*(-4 + p)*Power(1 + e*Cos(v),2)) + (Power(e,2)*(-4*Power(e,2) + Power(-2 + p,2))*p*Power(Sin(v),2))/
    ((1 - Power(e,2))*(-4 + p)*Sqrt((-4*Power(e,2) + Power(-2 + p,2))*(-6 + p - 2*e*Cos(v)))*(1 + e*Cos(v)));
}

double Phi(double p, double e){
	return 4*Sqrt(p/(-6 + 2*e + p))*EllipticK((4*e)/(-6 + 2*e + p));
}

double T_r(double p, double e){
	return (2*p*Sqrt((-4*Power(e,2) + Power(-2 + p,2))*(-6 + 2*e + p))*EllipticE((4*e)/(-6 + 2*e + p)))/((1 - Power(e,2))*(-4 + p)) - (2*Sqrt(-4*Power(e,2) + Power(-2 + p,2))*p*EllipticK((4*e)/(-6 + 2*e + p)))/((1 - Power(e,2))*Sqrt(-6 + 2*e + p)) - 
   (4*Sqrt(-4*Power(e,2) + Power(-2 + p,2))*(8*(1 - Power(e,2)) + (1 + 3*Power(e,2) - p)*p)*EllipticPi((-2*e)/(1 - e),(4*e)/(-6 + 2*e + p)))/((1 - e)*(1 - Power(e,2))*(-4 + p)*Sqrt(-6 + 2*e + p)) + 
   (16*Sqrt(-4*Power(e,2) + Power(-2 + p,2))*EllipticPi((4*e)/(-2 + 2*e + p),(4*e)/(-6 + 2*e + p)))/(Sqrt(-6 + 2*e + p)*(-2 + 2*e + p));
}
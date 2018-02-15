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

double dp_dchi(double p, double e, double v, double Fphi, double Fr);			// F1p
double de_dchi(double p, double e, double v, double Fphi, double Fr);			// F1e
double dw_dchi(double p, double e, double v, double Fphi, double Fr);			// f1v
double dt_dchi(double p, double e, double v);
double dphi_dchi(double p, double e, double v);

double dF1p_dp(double p, double e, double v, double Fphi, double Fr, double dFphi_dp, double dFr_dp);
double dF1p_de(double p, double e, double v, double Fphi, double Fr, double dFphi_de, double dFr_de);
double dF1p_dv(double p, double e, double v, double Fphi, double Fr, double dFphi_dv, double dFr_dv);

double dF1e_dp(double p, double e, double v, double Fphi, double Fr, double dFphi_dp, double dFr_dp);
double dF1e_de(double p, double e, double v, double Fphi, double Fr, double dFphi_de, double dFr_de);
double dF1e_dv(double p, double e, double v, double Fphi, double Fr, double dFphi_dv, double dFr_dv);

double df1v_dp(double p, double e, double v, double Fphi, double Fr, double dFphi_dp, double dFr_dp);
double df1v_de(double p, double e, double v, double Fphi, double Fr, double dFphi_de, double dFr_de);
double df1v_dv(double p, double e, double v, double Fphi, double Fr, double dFphi_dv, double dFr_dv);

double V0(double p, double e, double v);
double dV0_dp(double p, double e, double v);
double dV0_de(double p, double e, double v);
double dV0_dv(double p, double e, double v);

double U0(double p, double e, double v);
double dU0_dp(double p, double e, double v);
double dU0_de(double p, double e, double v);
double dU0_dv(double p, double e, double v);

double Phi(double p, double e);
double T_r(double p, double e);
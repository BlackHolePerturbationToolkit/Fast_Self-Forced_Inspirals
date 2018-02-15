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

#include <iostream>
#include <fstream>
#include <sstream>
#include <complex>
#include <vector>
#include <algorithm>
#include <chrono>
#include <math.h>
#include <cstring>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_matrix.h>

#include <fftw3.h>

#include <Interpolant.h>
#include <EoM.h>

using namespace std;
using namespace std::literals;
using namespace std::chrono;

typedef complex<double> Complex;

extern "C" {
	#include "lib_Sch_GSF.h"
}

// Used to set the mode the code runs in
#define FULL_INSPIRAL  0
#define NIT_INSPIRAL   1
#define DECOMPOSE      2
#define CONSTRUCT_Fs   3
#define WAVEFORM_FULL  4
#define WAVEFORM_NIT   5

// Functions for computing the RHS of the NIT EoM
void FFT_self_force_over_parameter_space();
void FFT_self_force(double p, double e, fftw_plan *plan, int N, fftw_complex *in, fftw_complex *out, ofstream *fv_file, ofstream *Fp_file, ofstream *Fe_file, ofstream *dFp_dp_file, ofstream *dFp_de_file, ofstream *dFe_dp_file, ofstream *dFe_de_file, ofstream *dV0_dp_file, ofstream *dV0_de_file, ofstream *dV0_dv_file, ofstream *dU0_dp_file, ofstream *dU0_de_file, ofstream *dU0_dv_file);
void construct_tilde_Fs();

// Fuctions for integrating the full EoM
void integrate_osc_eqs(double p0, double e0);
void integrate_osc_eqs_implicit(double p0, double e0);

// Function for integrating the NIT EoM
void interpolate_Fs_and_integrate_NIT_EoM(double p0, double e0);
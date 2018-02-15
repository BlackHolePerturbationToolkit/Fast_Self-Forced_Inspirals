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

#include <NIT_inspiral.h>
#include <libconfig.h++>

using namespace libconfig;

void compute_waveform(string insp_filename, string out_filename);
Complex waveform_h(double p, double e, double v, double phi, double Z, double Phi);

int mode;				// Selects the behavior of the code. Set to either FUll_INSPIRAL, NIT_INSPIRAL, DECOMPOSE or CONSTRUCT_Fs
double q;				// The mass ratio q = m_1/m_2
ofstream fout;			// The output file
Config cfg;

int main(int argc, char* argv[]){	
	double p0, e0;			// The initial (p,e) values for the evolutions
	
	if(argc  > 1){	
		if( !strcmp(argv[1], "-f") ){
			mode = FULL_INSPIRAL;
			cout << "# Mode: inspiral using full EoM" << endl;
		} else if( !strcmp(argv[1], "-n") ){
			mode = NIT_INSPIRAL;
			cout << "# Mode: inspiral using NIT EoM" << endl;
		} else if( !strcmp(argv[1], "-d") ){
			mode = DECOMPOSE;
			cout << "# Mode: decompose the self-force data to compute the F's and f's in NIT EoM" << endl;
		}else if( !strcmp(argv[1], "-c") ){
			mode = CONSTRUCT_Fs;
			cout << "# Mode: Compute the F's and f's in the NIT EoM" << endl;
		}else if( !strcmp(argv[1], "-w") ){
			if(argc == 6){
				p0 = atof(argv[2]);
				e0 = atof(argv[3]);	
				q  = atof(argv[4]);
				if( !strcmp(argv[5], "-n") ){
					mode = WAVEFORM_NIT; 
					cout << "Compute the NIT waveform" << endl;
				}else if( !strcmp(argv[5], "-f") ){
					 mode = WAVEFORM_FULL;
					 cout << "Compute the Full waveform" << endl;
				}else{
					cout << "Unrecognized waveform flag" << endl;
					exit(0);
				}
			}else{
				cout << "For the waveform mode please enter initial p, e and q values and a '-f' or '-n' flag for Full or NIT." << endl;
				exit(0);
			}
		}else{
			cout << "Unrecognized flag. Run with no arguments for instructions." << endl;
		}
		if(mode == FULL_INSPIRAL || mode == NIT_INSPIRAL){
			if(argc == 5){
				p0 = atof(argv[2]);
				e0 = atof(argv[3]);	
				q  = atof(argv[4]);
			}else{
				cout << "For inspiral modes please enter initial p, e value and a q value." << endl;
				exit(0);
			}
		}
	}else{
		cout << "Necessary parameters:" << endl;
		cout << "\t1. flag      '-f', '-n', 'w', '-d' or '-c' " << endl;
		cout << "\t   '-f'      Full inspiral" << endl;
		cout << "\t   '-n'      NIT inspiral" << endl;
		cout << "\t   '-w'      Compute the waveform" << endl;
		cout << "\t   '-d'      Decompose the self-force data into Fourier modes" << endl;
		cout << "\t   '-c'      Compute F's and f's in NIT EoM" << endl;
		cout << "\t2. p         Initial semi-latus rectum for '-f' or '-n' inspiral options" << endl;
		cout << "\t3. e         Initial eccentricity for '-f' or '-n' inspiral options" << endl;
		cout << "\t4. q         Mass ratio for '-f' or '-n' inspiral options" << endl;
		cout << "\t5. flag      Waveform only, '-f' for full, '-n' for NIT" << endl;
		exit(0);
	}	
	
	
	// Read in the parameters from the configuration file
	cout << "Reading configuration file" << endl;
	
	// Check the config file exists and is well formatted
    try{
		 cfg.readFile("config/parameters.cfg");
	} catch(const FileIOException &fioex){
	     std::cerr << "Config file could not be found at 'config/parameters.cfg'" << std::endl;
		 exit(0);
    } catch(const ParseException &pex){
    	std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine() << " - " << pex.getError() << std::endl;
  	    exit(0);
    }
	
	
	// Load the model and check that it loaded correctly	
	int model_status = lib_Sch_GSF_load_model();
	if(model_status == MODEL_LOAD_FAIL) exit(0);
	
	// Disable the GSL error handler (only for production version of code)
	// gsl_set_error_handler_off();
	
	// Output all the data at 10 digits of precision using scientific notation
	fout.precision(10);
	cout.precision(10);
	fout << scientific;
	cout << scientific;
	
	// Perform the calculation depending upon which mode is selected
	if(mode == FULL_INSPIRAL){
		
		ostringstream filename;
		filename << "output/Inspiral_Full_p" << p0 << "_e" << e0 << "_q" << q <<".dat";
		cout << "Outputting inspiral trajectory to " << filename.str() << endl;	
		
		fout.open(filename.str());
		fout << "# Full Inspiral" << endl;
		integrate_osc_eqs(p0, e0);	
			
	}else if(mode == NIT_INSPIRAL){
		
		ostringstream filename;
		filename << "output/Inspiral_NIT_p" << p0 << "_e" << e0 << "_q" << q <<".dat";
		cout << "Outputting inspiral trajectory to " << filename.str() << endl;
	
		fout.open(filename.str());
		fout << "# NIT Inspiral" << endl;
		interpolate_Fs_and_integrate_NIT_EoM(p0, e0);
		
	}else if(mode == WAVEFORM_NIT){
		
		ostringstream out_filename, insp_filename;
		out_filename << "output/Waveform_NIT_p" << p0 << "_e" << e0 << "_q" << q <<".dat";
		insp_filename << "output/Inspiral_NIT_p" << p0 << "_e" << e0 << "_q" << q <<".dat";
		compute_waveform(insp_filename.str(), out_filename.str());
			
	}else if(mode == WAVEFORM_FULL){
		
		ostringstream out_filename, insp_filename;
		out_filename << "output/Waveform_Full_p" << p0 << "_e" << e0 << "_q" << q <<".dat";
		insp_filename << "output/Inspiral_Full_p" << p0 << "_e" << e0 << "_q" << q <<".dat";
		compute_waveform(insp_filename.str(), out_filename.str());
			
	}else if(mode == DECOMPOSE){
		FFT_self_force_over_parameter_space();
	}else if(mode == CONSTRUCT_Fs){
		construct_tilde_Fs();
	}
	
	
	// Close the output file
	fout.close();
}

// Functions to integrate the NIT EoM

// Used to pass the interpolants to the NIT_EoM function
struct interp_params{
	Interpolant *Fp1;
	Interpolant *Fe1;
	Interpolant *fv1;
	Interpolant *Fp2;
	Interpolant *Fe2;
	Interpolant *U1;
	Interpolant *V1;
};

int NIT_EoM (double v, const double y[], double f[], void *params){
	(void)(v); /* avoid unused parameter warning */
	
	struct interp_params *interps = (struct interp_params *)params;

	double p = y[0];
	double e = y[1];

	f[0] = q*interps->Fp1->eval(p-2*e, e) + q*q*interps->Fp2->eval(p-2*e, e);
	f[1] = q*interps->Fe1->eval(p-2*e, e) + q*q*interps->Fe2->eval(p-2*e, e);
	f[2] = 1.0 + q*interps->fv1->eval(p-2*e, e);
	f[3] = T_r(p, e)/(2.0*M_PI) + q*interps->U1->eval(p-2*e, e);
	f[4] = Phi(p, e)/(2.0*M_PI) + q*interps->V1->eval(p-2*e, e);
	
	if(p-6-2*e > 0.1) return GSL_SUCCESS;
	else return GSL_SUCCESS + 1;
}

void interpolate_Fs_and_integrate_NIT_EoM(double p0, double e0){
	ifstream Ftilde_file("data/Ftildes.dat");
	
	// Load the data for the F/f/U/V on the RHS of the NIT EoM
	string Ftilde_string;
	vector<double> ys, es, Fp1s, Fe1s, fv1s, Fp2s, Fe2s, U1s, V1s, Xv1s, Yp1s, Ye1s;
	double y, e, Fp1, Fe1, fv1, Fp2, Fe2, U1, V1, Xv1, Yp1, Ye1;
	while(getline(Ftilde_file, Ftilde_string)){
		
		stringstream Ftilde_ss(Ftilde_string);
		
		Ftilde_ss >> y >> e >> Fp1 >> Fe1 >> fv1 >> Fp2 >> Fe2 >> U1 >> V1 >> Xv1 >> Yp1 >> Ye1;
		
		ys.push_back(y);
		es.push_back(e);
		Fp1s.push_back(Fp1);
		Fe1s.push_back(Fe1);
		fv1s.push_back(fv1);
		Fp2s.push_back(Fp2);
		Fe2s.push_back(Fe2);
		U1s.push_back(U1);
		V1s.push_back(V1);
		Xv1s.push_back(Xv1);
		Yp1s.push_back(Yp1);
		Ye1s.push_back(Ye1);
	}	
	
	// Interpolate the data
	Interpolant Fp1_interp(ys, es, Fp1s);
	Interpolant Fe1_interp(ys, es, Fe1s);
	Interpolant fv1_interp(ys, es, fv1s);
	Interpolant Fp2_interp(ys, es, Fp2s);
	Interpolant Fe2_interp(ys, es, Fe2s);
	Interpolant U1_interp(ys, es, U1s);
	Interpolant V1_interp(ys, es, V1s);
	Interpolant Xv1_interp(ys, es, Xv1s);
	Interpolant Yp1_interp(ys, es, Yp1s);
	Interpolant Ye1_interp(ys, es, Ye1s);
	
	// Numerically integrate the NIT inspiral below
	
	struct interp_params interps = {&Fp1_interp, &Fe1_interp, &fv1_interp, &Fp2_interp, &Fe2_interp, &U1_interp, &V1_interp};
		
	double t0 = 0;
	double phi0 = 0;
	
	double chi = 0;

	// Initial conditions
	double ptilde0 = p0 + q*Yp1_interp.eval(p0-2*e0, e0);
	double etilde0 = e0 + q*Ye1_interp.eval(p0-2*e0, e0);
	double vtilde0 = 0  + q*Xv1_interp.eval(p0-2*e0, e0);
	double y1[5] = {ptilde0, etilde0, vtilde0, t0, phi0};
	
	// Output the initial parameters
	fout << "# Format: chi    p    e    xi    t    phi" << endl;
	fout << chi << " " << y1[0] << " " << y1[1] << " " << y1[2] << " " << y1[3] << " " << y1[4] << " " << endl;
	
	int dense_output;
	try{
		dense_output 		= cfg.lookup("Dense_output");		
 	}catch(const SettingNotFoundException &nfex){
   		cerr << "'Dense_output' setting missing from configuration file." << endl; exit(0);
 	}
	
	if(dense_output == 1){				// Dense output
		
		double n_per_orbit;
		try{
			n_per_orbit 		= cfg.lookup("n_per_orbit");		
	 	}catch(const SettingNotFoundException &nfex){
	   		cerr << "'n_per_orbit' setting missing from configuration file." << endl; exit(0);
	 	}

		double chi_i = 0;
		double delta_chi = 2.0*M_PI/n_per_orbit;
	
	    gsl_odeiv2_system sys 	= {NIT_EoM, NULL, 5, &interps};	
	    gsl_odeiv2_driver *d 	= gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd, 1e-10, 1e-10, 0.0);
	
	    high_resolution_clock::time_point t1 = high_resolution_clock::now();
	
		while(1){
			chi_i += delta_chi;
			int status = gsl_odeiv2_driver_apply (d, &chi, chi_i, y1);

			if (status != GSL_SUCCESS){
				//printf ("error, return value=%d\n", status);
				break;
			}
			fout << chi << " " << y1[0] << " " << y1[1] << " " << y1[2] << " " << y1[3] << " " << y1[4] << " " << endl;
		}
		
	  high_resolution_clock::time_point t2 = high_resolution_clock::now();
	  
	  duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	  fout << "# Computing the inspiral took: " << time_span.count() << " seconds." << endl;

		gsl_odeiv2_driver_free (d);
		
	}else{		// Sparse output
	    const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rk8pd;

	    gsl_odeiv2_step *s 		= gsl_odeiv2_step_alloc (T, 5);
	    gsl_odeiv2_control *c 	= gsl_odeiv2_control_y_new (1e-10, 1e-10);
	    gsl_odeiv2_evolve *e 	= gsl_odeiv2_evolve_alloc (5);

	    gsl_odeiv2_system sys = {NIT_EoM, NULL, 5, &interps};

	    double chi_max = 1000000000.0;	//FIXME make this something like 100/q
	    double h = 1e-1;

	    high_resolution_clock::time_point t1 = high_resolution_clock::now();
		
		double chi_prev = 1.;
	    while (chi < chi_max) {
	        int status = gsl_odeiv2_evolve_apply (e, c, s, &sys, &chi, chi_max, &h, y1);
			if (status != GSL_SUCCESS){
				//printf ("error, return value=%d\n", status);
				break;
			}

			fout << chi << " " << y1[0] << " " << y1[1] << " " << y1[2] << " " << y1[3] << " " << y1[4] << " " << endl;
			
			// Stop output of lots of data near the separatrix when the time step gets very small
			if(fabs(chi_prev/chi - 1.0) < 1e-9) break;
			
		    chi_prev = chi;	
	      }
		  high_resolution_clock::time_point t2 = high_resolution_clock::now();
		  
		  duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
		  fout << "# Computing the inspiral took: " << time_span.count() << " seconds." << endl;
		  
		  
		  
		  gsl_odeiv2_evolve_free (e);
		  gsl_odeiv2_control_free (c);
		  gsl_odeiv2_step_free (s);
	}
	
	
}


// Functions for computing the RHS of the NIT EoM
void FFT_self_force_over_parameter_space(){
	double p, e;
	int N = 50;
    fftw_complex *in, *out;
    fftw_plan plan;

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    plan = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	
	// Open file to output the Fourier data
    ofstream fv_file, Fp_file, Fe_file;
	ofstream dFp_dp_file, dFp_de_file, dFe_dp_file, dFe_de_file;
	ofstream dV0_dp_file, dV0_de_file, dV0_dv_file;
	ofstream dU0_dp_file, dU0_de_file, dU0_dv_file;
	
	fv_file.open("data/Fourier_fv.dat");
	Fp_file.open("data/Fourier_Fp.dat");
	Fe_file.open("data/Fourier_Fe.dat");
	
	dFp_dp_file.open("data/Fourier_dFp_dp.dat");	
	dFp_de_file.open("data/Fourier_dFp_de.dat");
	dFe_dp_file.open("data/Fourier_dFe_dp.dat");
	dFe_de_file.open("data/Fourier_dFe_de.dat");
	
	dV0_dp_file.open("data/Fourier_dV0_dp.dat");
	dV0_de_file.open("data/Fourier_dV0_de.dat");
	dV0_dv_file.open("data/Fourier_dV0_dv.dat");
	
	dU0_dp_file.open("data/Fourier_dU0_dp.dat");
	dU0_de_file.open("data/Fourier_dU0_de.dat");
	dU0_dv_file.open("data/Fourier_dU0_dv.dat");

	double de = 0.002;
	double dp = 0.05;
	e = 0;
	for(int i=1; i <= 100; i++){
		e += de;
		for(int j=0; j < 100; j++){
			p = 6 + 2*e + 0.05  + j*dp;
			FFT_self_force(p, e, &plan, N, in, out, &fv_file, &Fp_file, &Fe_file, &dFp_dp_file, &dFp_de_file, &dFe_dp_file, &dFe_de_file, &dV0_dp_file, &dV0_de_file, &dV0_dv_file, &dU0_dp_file, &dU0_de_file, &dU0_dv_file);
		}
	}
	
	// Close all the files after writing the data out    
	fv_file.close();
	Fp_file.close();
	Fe_file.close();	
	
	dFp_dp_file.close();	
	dFp_de_file.close();
	dFe_dp_file.close();
	dFe_de_file.close();
	
	dV0_dp_file.close();
	dV0_de_file.close();
	dV0_dv_file.close();
	
	dU0_dp_file.close();
	dU0_de_file.close();
	dU0_dv_file.close();
		
    fftw_destroy_plan(plan);
    fftw_free(in); fftw_free(out);
	
}

void FFT_self_force(double p, double e, fftw_plan *plan, int N, fftw_complex *in, fftw_complex *out, ofstream *fv_file, ofstream *Fp_file, ofstream *Fe_file, ofstream *dFp_dp_file, ofstream *dFp_de_file, ofstream *dFe_dp_file, ofstream *dFe_de_file,ofstream *dV0_dp_file, ofstream *dV0_de_file, ofstream *dV0_dv_file, ofstream *dU0_dp_file, ofstream *dU0_de_file, ofstream *dU0_dv_file){
		
	double *Fr = new double[N];
	double *Fphi = new double[N];
	
	double *dFr_dp = new double[N];
	double *dFr_de = new double[N];
	double *dFphi_dp = new double[N];
	double *dFphi_de = new double[N];
	
	double *v = new double[N];
	
	// The number of Fourier modes to output
	int N_out = 10;
	
	// Output the coordinates in phase space, y = p-2e, e 
	*Fp_file << p-2*e << " " << e << " ";
	*Fe_file << p-2*e << " " << e << " ";
	*fv_file << p-2*e << " " << e << " ";
	
	*dFp_dp_file << p-2*e << " " << e << " ";
	*dFp_de_file << p-2*e << " " << e << " ";
	*dFe_dp_file << p-2*e << " " << e << " ";
	*dFe_de_file << p-2*e << " " << e << " ";
	
	*dV0_dp_file << p-2*e << " " << e << " ";
	*dV0_de_file << p-2*e << " " << e << " ";
	*dV0_dv_file << p-2*e << " " << e << " ";
	
	*dU0_dp_file << p-2*e << " " << e << " ";
	*dU0_de_file << p-2*e << " " << e << " ";
	*dU0_dv_file << p-2*e << " " << e << " ";
	
	
	// Cache SF data and compute fv1
	for(int i = 0; i < N; i++){
		v[i] = i*2.0*M_PI/N;
		Fr[i] 	= lib_Sch_GSF_Fr_diss(e, p, v[i]) + lib_Sch_GSF_Fr_cons(e, p, v[i]);
		Fphi[i] = lib_Sch_GSF_Fphi_diss(e, p, v[i]) + lib_Sch_GSF_Fphi_cons(e, p, v[i]);
		
		in[i][0] = dw_dchi(p, e, v[i], Fphi[i], Fr[i]);
	}

    fftw_execute(*plan); /* repeat as needed */
	for(int i = 0; i < N_out; i++){
   	 	*fv_file << out[i][0]/N << " " << out[i][1]/N << " ";
	}
	*fv_file << endl;
	
	// Compute F1p		
	for(int i = 0; i < N; i++){
		in[i][0] = dp_dchi(p, e, v[i], Fphi[i], Fr[i]);
	}
    fftw_execute(*plan); /* repeat as needed */
	for(int i = 0; i < N_out; i++){
   	 	*Fp_file << out[i][0]/N << " " << out[i][1]/N << " ";
	}
	*Fp_file << endl;
	
	// Compute F1e
	for(int i = 0; i < N; i++){
		in[i][0] = de_dchi(p, e, v[i], Fphi[i], Fr[i]);
	}
    fftw_execute(*plan); /* repeat as needed */
	for(int i = 0; i < N_out; i++){
   	 	*Fe_file << out[i][0]/N << " " << out[i][1]/N << " ";
	}
	*Fe_file << endl;
	
	// The derivatives of F and U w.r.t. p, e, and v are below
	
	// Compute dF1p_dp
	for(int i = 0; i < N; i++){
		dFr_dp[i] = lib_Sch_GSF_dFr_cons_dp(e, p, v[i]) + lib_Sch_GSF_dFr_diss_dp(e, p, v[i]);
		dFr_de[i] = lib_Sch_GSF_dFr_cons_de(e, p, v[i]) + lib_Sch_GSF_dFr_diss_de(e, p, v[i]);
		
		dFphi_dp[i] = lib_Sch_GSF_dFphi_cons_dp(e, p, v[i]) + lib_Sch_GSF_dFphi_diss_dp(e, p, v[i]);
		dFphi_de[i] = lib_Sch_GSF_dFphi_cons_de(e, p, v[i]) + lib_Sch_GSF_dFphi_diss_de(e, p, v[i]);
		
		in[i][0] = dF1p_dp(p, e, v[i], Fphi[i], Fr[i], dFphi_dp[i], dFr_dp[i]);
	}
	fftw_execute(*plan); /* repeat as needed */
	for(int i = 0; i < N_out; i++){
   	 	*dFp_dp_file << out[i][0]/N << " " << out[i][1]/N << " ";
	}
	*dFp_dp_file << endl;
	
	//Compute dF1p_de
	for(int i = 0; i < N; i++){
		in[i][0] = dF1p_de(p, e, v[i], Fphi[i], Fr[i], dFphi_de[i], dFr_de[i]);
	}
	fftw_execute(*plan); /* repeat as needed */
	for(int i = 0; i < N_out; i++){
   	 	*dFp_de_file << out[i][0]/N << " " << out[i][1]/N << " ";
	}
	*dFp_de_file << endl;
	
	//Compute dF1e_dp
	for(int i = 0; i < N; i++){
		in[i][0] = dF1e_dp(p, e, v[i], Fphi[i], Fr[i], dFphi_dp[i], dFr_dp[i]);
	}
	fftw_execute(*plan); /* repeat as needed */
	for(int i = 0; i < N_out; i++){
   	 	*dFe_dp_file << out[i][0]/N << " " << out[i][1]/N << " ";
	}
	*dFe_dp_file << endl;
	
	//Compute dF1e_de
	for(int i = 0; i < N; i++){
		in[i][0] = dF1e_de(p, e, v[i], Fphi[i], Fr[i], dFphi_de[i], dFr_de[i]);
	}
	fftw_execute(*plan); /* repeat as needed */
	for(int i = 0; i < N_out; i++){
   	 	*dFe_de_file << out[i][0]/N << " " << out[i][1]/N << " ";
	}
	*dFe_de_file << endl;
	
	
	
	
	//Compute dV0_dp
	for(int i = 0; i < N; i++){
		in[i][0] = dV0_dp(p, e, v[i]);
	}
	fftw_execute(*plan); /* repeat as needed */
	for(int i = 0; i < N_out; i++){
   	 	*dV0_dp_file << out[i][0]/N << " " << out[i][1]/N << " ";
	}
	*dV0_dp_file << endl;
	
	//Compute dV0_de
	for(int i = 0; i < N; i++){
		in[i][0] = dV0_de(p, e, v[i]);
	}
	fftw_execute(*plan); /* repeat as needed */
	for(int i = 0; i < N_out; i++){
   	 	*dV0_de_file << out[i][0]/N << " " << out[i][1]/N << " ";
	}
	*dV0_de_file << endl;
	
	//Compute dV0_dv
	for(int i = 0; i < N; i++){
		in[i][0] = dV0_dv(p, e, v[i]);
	}
	fftw_execute(*plan); /* repeat as needed */
	for(int i = 0; i < N_out; i++){
   	 	*dV0_dv_file << out[i][0]/N << " " << out[i][1]/N << " ";
	}
	*dV0_dv_file << endl;
	
	
	
	
	
	//Compute dU0_dp
	for(int i = 0; i < N; i++){
		in[i][0] = dU0_dp(p, e, v[i]);
	}
	fftw_execute(*plan); /* repeat as needed */
	for(int i = 0; i < N_out; i++){
   	 	*dU0_dp_file << out[i][0]/N << " " << out[i][1]/N << " ";
	}
	*dU0_dp_file << endl;
	
	//Compute dU0_de
	for(int i = 0; i < N; i++){
		in[i][0] = dU0_de(p, e, v[i]);
	}
	fftw_execute(*plan); /* repeat as needed */
	for(int i = 0; i < N_out; i++){
   	 	*dU0_de_file << out[i][0]/N << " " << out[i][1]/N << " ";
	}
	*dU0_de_file << endl;
	
	//Compute dU0_dv
	for(int i = 0; i < N; i++){
		in[i][0] = dU0_dv(p, e, v[i]);
	}
	fftw_execute(*plan); /* repeat as needed */
	for(int i = 0; i < N_out; i++){
   	 	*dU0_dv_file << out[i][0]/N << " " << out[i][1]/N << " ";
	}
	*dU0_dv_file << endl;
	
}

void construct_tilde_Fs(){
	// Open the Fourier data files
	ifstream fv_file("data/Fourier_fv.dat");
	ifstream Fp_file("data/Fourier_Fp.dat");
	ifstream Fe_file("data/Fourier_Fe.dat");
	
	ifstream dFpdp_file("data/Fourier_dFp_dp.dat");
	ifstream dFpde_file("data/Fourier_dFp_de.dat");
	ifstream dFedp_file("data/Fourier_dFe_dp.dat");
	ifstream dFede_file("data/Fourier_dFe_de.dat");
	
	ifstream dU0dp_file("data/Fourier_dU0_dp.dat");
	ifstream dU0de_file("data/Fourier_dU0_de.dat");
	ifstream dU0dv_file("data/Fourier_dU0_dv.dat");
	ifstream dV0dp_file("data/Fourier_dV0_dp.dat");
	ifstream dV0de_file("data/Fourier_dV0_de.dat");
	ifstream dV0dv_file("data/Fourier_dV0_dv.dat");
	
	// The output file for the Ftilde data
	ofstream Ftilde_file("data/Ftildes.dat");
	Ftilde_file.precision(10);
	
	string fv_string, Fp_string, Fe_string, dFpdp_string, dFpde_string, dFedp_string, dFede_string;
	string dU0dp_string, dU0de_string, dU0dv_string, dV0dp_string, dV0de_string, dV0dv_string;
	double p, e, fv1, Fp1, Fe1, null;
	while(getline(fv_file, fv_string)){			// We assume all the Fourier files are on identical grids
		getline(Fp_file, Fp_string);
		getline(Fe_file, Fe_string);
		getline(dFpdp_file, dFpdp_string);
		getline(dFpde_file, dFpde_string);
		getline(dFedp_file, dFedp_string);
		getline(dFede_file, dFede_string);	
		getline(dU0dp_file, dU0dp_string);
		getline(dU0de_file, dU0de_string);
		getline(dU0dv_file, dU0dv_string);
		getline(dV0dp_file, dV0dp_string);
		getline(dV0de_file, dV0de_string);
		getline(dV0dv_file, dV0dv_string);	
		
		stringstream fv_ss(fv_string);
		stringstream Fp_ss(Fp_string);
		stringstream Fe_ss(Fe_string);
		stringstream dFpdp_ss(dFpdp_string);
		stringstream dFpde_ss(dFpde_string);
		stringstream dFedp_ss(dFedp_string);
		stringstream dFede_ss(dFede_string);
		stringstream dU0dp_ss(dU0dp_string);
		stringstream dU0de_ss(dU0de_string);
		stringstream dU0dv_ss(dU0dv_string);
		stringstream dV0dp_ss(dV0dp_string);
		stringstream dV0de_ss(dV0de_string);
		stringstream dV0dv_ss(dV0dv_string);
		
		Fp_ss >> p >> e >> Fp1 >> null;
		Fe_ss >> p >> e >> Fe1 >> null;
		fv_ss >> p >> e >> fv1 >> null;
		fv1 = -fv1;			// Note the sign change as xi = chi - chi0
		
		dFpdp_ss >> p >> e >> null >> null;
		dFpde_ss >> p >> e >> null >> null;
		dFedp_ss >> p >> e >> null >> null;
		dFede_ss >> p >> e >> null >> null;
		
		dU0dp_ss >> p >> e >> null >> null;
		dU0de_ss >> p >> e >> null >> null;
		dU0dv_ss >> p >> e >> null >> null;
		dV0dp_ss >> p >> e >> null >> null;
		dV0de_ss >> p >> e >> null >> null;
		dV0dv_ss >> p >> e >> null >> null;
		
		double re, im;
		double Fp2 = 0, Fe2 = 0, U1 = 0, V1 = 0, Xv1 = 0, Yp1 = 0, Ye1 = 0;
		Complex Fpk, Fek, fvk, dFpdpk, dFpdek, dFedpk, dFedek;
		Complex dU0dpk, dU0dek, dU0dvk, dV0dpk, dV0dek, dV0dvk;
		for(int k = 1; k < 10; k++){
			Fp_ss >> re >> im;
			Fpk = re + 1i*im;
			
			Fe_ss >> re >> im;
			Fek = re + 1i*im;
			
			fv_ss >> re >> im;
			fvk = -re - 1i*im;			// Note the sign change as xi = chi - chi0
			
			dFpdp_ss >> re >> im;
			dFpdpk = re + 1i*im;
			
			dFpde_ss >> re >> im;
			dFpdek = re + 1i*im;
			
			dFedp_ss >> re >> im;
			dFedpk = re + 1i*im;
			
			dFede_ss >> re >> im;
			dFedek = re + 1i*im;
			
			dU0dp_ss >> re >> im;
			dU0dpk = re + 1i*im;
			
			dU0de_ss >> re >> im;
			dU0dek = re + 1i*im;
			
			dU0dv_ss >> re >> im;
			dU0dvk = re + 1i*im;
			
			dV0dp_ss >> re >> im;
			dV0dpk = re + 1i*im;
			
			dV0de_ss >> re >> im;
			dV0dek = re + 1i*im;
			
			dV0dv_ss >> re >> im;
			dV0dvk = re + 1i*im;
			
			Fp2 += 2.*( 1.i/((double)k)*(dFpdpk*conj(Fpk) + dFpdek*conj(Fek)) - Fpk*conj(fvk) ).real();
			Fe2 += 2.*( 1.i/((double)k)*(dFedpk*conj(Fpk) + dFedek*conj(Fek)) - Fek*conj(fvk) ).real();
			
			U1 += 2.*(dU0dpk*conj(Fpk) + dU0dek*conj(Fek) + dU0dvk*conj(fvk)).real();
			V1 += 2.*(dV0dpk*conj(Fpk) + dV0dek*conj(Fek) + dV0dvk*conj(fvk)).real();
			
			// The below is for the initial condition matching and assumes that v0 = 0
			Xv1 += 2.*(1.i/((double)k)*fvk).real();
			Yp1 += 2.*(1.i/((double)k)*Fpk).real();
			Ye1 += 2.*(1.i/((double)k)*Fek).real();			
		}
		Ftilde_file  << p << " " << e << " " << Fp1 << " " << Fe1 << " " << fv1 << " " << Fp2 << " " << Fe2 << " " << U1 << " " << V1 << " " << Xv1 << " " << Yp1 << " " << Ye1 << endl;
	};
	
	
	fv_file.close();
	Fp_file.close();
	Fe_file.close();
	Ftilde_file.close();
}


// Code for explicit integration of the {p, e, v, t, phi} equations using the full self-force

int osc_eqs (double chi, const double y[], double f[], void *params){
	double p = y[0];
	double e = y[1];
	double chi0 = y[2];
	double v = chi - chi0;
	
	double Fr 	= q*(lib_Sch_GSF_Fr_diss(e, p, v) + lib_Sch_GSF_Fr_cons(e, p, v));
	double Fphi = q*(lib_Sch_GSF_Fphi_diss(e, p, v) + lib_Sch_GSF_Fphi_cons(e, p, v));

	f[0] = dp_dchi(p, e, v, Fphi, Fr);
	f[1] = de_dchi(p, e, v, Fphi, Fr);
	f[2] = dw_dchi(p, e, v, Fphi, Fr);
	f[3] = dt_dchi(p, e, v);
	f[4] = dphi_dchi(p, e, v);	
	
	if(p-6-2*e > 0.1) return GSL_SUCCESS;
	else return GSL_SUCCESS + 1;
}

void integrate_osc_eqs(double p0, double e0){
	
	double chi = 0;
	
	// y[0] = p, y[1] = e, y[2] = chi0, y[3] = t, y[4] = phi
	double y[5] = {p0, e0, 0, 0, 0};
	
	// Output the initial parameters
	fout << "# Format: chi    p    e    chi0    t    phi" << endl;
	fout << chi << " " << y[0] << " " << y[1] << " " << y[2] << " " << y[3] << " " << y[4] << " " << endl;


	int dense_output;
	try{
		dense_output 		= cfg.lookup("Dense_output");		
 	}catch(const SettingNotFoundException &nfex){
   		cerr << "'Dense_output' setting missing from configuration file." << endl; exit(0);
 	}


	if(dense_output == 1){
		
		double n_per_orbit;
		try{
			n_per_orbit 		= cfg.lookup("n_per_orbit");		
	 	}catch(const SettingNotFoundException &nfex){
	   		cerr << "'n_per_orbit' setting missing from configuration file." << endl; exit(0);
	 	}
		
	    gsl_odeiv2_system sys 	= {osc_eqs, NULL, 5, NULL};	
	    gsl_odeiv2_driver *d 	= gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd, 1e-10, 1e-10, 0.0);
		
		double chi_i = 0;
		double delta_chi = 2.0*M_PI/n_per_orbit;
		
		while(1){
			chi_i += delta_chi;
			int status = gsl_odeiv2_driver_apply (d, &chi, chi_i, y);

			if (status != GSL_SUCCESS){
				printf ("error, return value=%d\n", status);
				break;
			}

			fout << chi << " " << y[0] << " " << y[1] << " " << y[2] << " " << y[3] << " " << y[4] << " " << endl;
			
		}

		gsl_odeiv2_driver_free (d);
	}else{		// Sparse output
	    const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rk8pd;

	    gsl_odeiv2_step *s 		= gsl_odeiv2_step_alloc (T, 5);
	    gsl_odeiv2_control *c 	= gsl_odeiv2_control_y_new (1e-10, 1e-10);
	    gsl_odeiv2_evolve *e 	= gsl_odeiv2_evolve_alloc (5);

	    gsl_odeiv2_system sys = {osc_eqs, NULL, 5, NULL};

	    double chi_max = 1000000000.0;	//FIXME make this something like 100/q
	    double h = 1e-6;

	    high_resolution_clock::time_point t1 = high_resolution_clock::now();

	    while (chi < chi_max) {
	        int status = gsl_odeiv2_evolve_apply (e, c, s, &sys, &chi, chi_max, &h, y);
			if (status != GSL_SUCCESS){
				//printf ("error, return value=%d\n", status);
				break;
			}

			fout << chi << " " << y[0] << " " << y[1] << " " << y[2] << " " << y[3] << " " << y[4] << " " << endl;
			
	      }
		  high_resolution_clock::time_point t2 = high_resolution_clock::now();
		  
		  duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
		  fout << "# Computing the inspiral took: " << time_span.count() << " seconds." << endl;
		  
		  
		  
		  gsl_odeiv2_evolve_free (e);
		  gsl_odeiv2_control_free (c);
		  gsl_odeiv2_step_free (s);
	}
	
}


// Waveform generation functions below

void compute_waveform(string insp_filename, string out_filename){

	// Check if the assocaited inspiral trajectory file exists
	ifstream insp(insp_filename);
	if(!insp){
		cout << "Inspiral file: " << insp_filename << " does not exist." << endl;
		exit(0);
	}

	// Load and interpolate the inspiral trajectory
	cout << "Loading and interpolating inspiral trajectory data"  << endl;
	
	// FIXME only load data up to t_max (rather than the entire phase space trajectory as is currently done)	
	string insp_string;
	vector<double> chis, ps, es, vs, ts, phis;
	double chi, p, e, v, t, phi;
	while(getline(insp, insp_string)){
				
		if(insp_string.at(0) == '#') continue;
		
		stringstream insp_ss(insp_string);
		
		insp_ss >> chi >> p >> e >> v >> t >> phi;
				
		chis.push_back(chi);
		ps.push_back(p);
		es.push_back(e);
		vs.push_back(v);
		ts.push_back(t);
		phis.push_back(phi);
	}
	
	Interpolant p_interp(chis, ps);
	Interpolant e_interp(chis, es);
	Interpolant v_interp(chis, vs);
	Interpolant t_interp(chis, ts);
	Interpolant phi_interp(chis, phis);
	
	// Read in the parameters and give an error if they do not exist
	double M_solar, Deltat_sec;
	int i_max;
	try{
		M_solar 		= cfg.lookup("M_solar");		// Mass of the primary in solar masses}
 	}catch(const SettingNotFoundException &nfex){
   		cerr << "'M_solar' setting missing from configuration file." << endl; exit(0);
 	}
	
	try{
		Deltat_sec 	= cfg.lookup("Deltat_sec");		// Time step in seconds
	}catch(const SettingNotFoundException &nfex){
   		cerr << "'Delta_sec' setting missing from configuration file." << endl; exit(0);
 	}
	
	try{
		i_max 	= cfg.lookup("i_max");				// The number of time steps to take
	}catch(const SettingNotFoundException &nfex){
   		cerr << "'i_max' setting missing from configuration file." << endl; exit(0);
 	}
			 	
	double Msolar_sec	= 4.9253908380897e-6;	// Solar mass in seconds				
	double M 			= M_solar * Msolar_sec;	// Mass of the primary in seconds
	double Deltat		= Deltat_sec/M;			// Time step in units of M	
	
	cout << "Computing the waveform" << endl;
	fout.open(out_filename);
	
	fout << "# M     = " << M_solar << " [Solar Masses]" << endl;
	fout << "# q     = " << q << endl;
	fout << "# Δt    = " << Deltat_sec << " [Seconds]" << endl;
	fout << "# t_end = " << Deltat_sec*i_max/60./60./24./30. << " [30 day months]" << endl;
	
	
	// Resample t more densely
	// FIXME only resample up to t_max
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
	
	double t_max = i_max*Deltat;
	vector<double> t_dense, chi_dense;
	chi = 0;
	while(chi < chis[chis.size() - 1]){
	
		p = p_interp.eval(chi);
		e = e_interp.eval(chi);
		v = v_interp.eval(chi);
	
		chi_dense.push_back(chi);
		if(mode == WAVEFORM_NIT) t = t_interp.eval(chi) - U0(p,e,v);
		else t = t_interp.eval(chi);
		
		t_dense.push_back(t);
		
		if(t > t_max) break;	
		
		chi += 2.0*M_PI/10.;
	
	}
	Interpolant chi_interp(t_dense, chi_dense);
	high_resolution_clock::time_point t2 = high_resolution_clock::now();

	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	fout << "# Resampling t took: " << time_span.count() << " seconds." << endl;
	fout << "# Format: t h+ h×" << endl;

    t1 = high_resolution_clock::now();
	
	Complex h;
	vector<double> hplus, hcross;
	for(int i = 0; i <= i_max; i ++){
		t 		= i*Deltat;
		chi 	= chi_interp.eval(t);
		
		p = p_interp.eval(chi);
		e = e_interp.eval(chi);
		phi = phi_interp.eval(chi);
		
		if(mode == WAVEFORM_NIT){
			v = v_interp.eval(chi);
			// Note: The below is really part of reconstructing the physical trajectory, but we don't need it until now. 
			// It's quick to evaluate (in comparison to U0) so it has only a small effect on the timing
			phi -= V0(p,e,v);						
		}else{
			v = chi - v_interp.eval(chi);
		}
				
		h = waveform_h(p, e, v, phi, 0., 0.);
		
		
		hplus.push_back(h.real());
		hcross.push_back(h.imag());
	}
	
	t2 = high_resolution_clock::now();

	time_span = duration_cast<duration<double>>(t2 - t1);
	fout << "# Computing the waveform took: " << time_span.count() << " seconds." << endl;
	
	// Output the waveform data
	cout << "Outputting waveform to " << out_filename << endl;
	double t_sec;
	for(int i = 0; i <= i_max; i ++){
		
		t_sec 	= i*Deltat_sec;
		
		fout << t_sec << " " << hplus[i] << " " << hcross[i] << endl;
	}
	
}

// Definitions for using CForm'ed output from Mathematica
#define Sin(x)          (sin((double)(x)))
#define Cos(x)          (cos((double)(x)))
#define Power(x, y)     (pow((double)(x), (double)(y)))
#define Sqrt(x)         (sqrt((double)(x)))
#define Conjugate(x)	(conj(x))
#define Pi				M_PI

// Returns the quadrupolar waveform, Z = Cos Theta and Theta and Phi are the angles between the detector and the source
Complex waveform_h(double p, double e, double v, double phi, double Z, double Phi){
	
	Complex XX = (5*Power(e,4) + Power(e,2)*(52 - 10*p) - 8*(-2 + p)*p)*Sqrt(-6 + p - 2*e*Cos(v)) + 
   e*(48 - Power(e,2)*(-56 + p) + 4*(4 - 3*p)*p)*Cos(v)*Sqrt(-6 + p - 2*e*Cos(v)) + 
   2*Power(e,2)*(30 + 4*Power(e,2) + (5 - 2*p)*p)*Sqrt(-6 + p - 2*e*Cos(v))*Cos(2*v) + Power(e,3)*(24 + p)*Sqrt(-6 + p - 2*e*Cos(v))*Cos(3*v) + 
   3*Power(e,4)*Sqrt(-6 + p - 2*e*Cos(v))*Cos(4*v) - Complex(0,8)*e*Sqrt(p)*(-6 + p - 2*e*Cos(v))*(-1 + p - e*Cos(v))*(1 + e*Cos(v))*Sin(v);
	
	return ((exp(2.*Complex(0,1.)*(-phi + Phi))*XX*Power(1 - Z,2) + exp(2.*Complex(0,1.)*(phi - Phi))*Power(1 + Z,2)*Conjugate(XX))*(-2 + p - 2*e*Cos(v)))/
    (4.*(-4*Power(e,2) + Power(-2 + p,2))*Power(p,1.5)*Sqrt(-6 + p - 2*e*Cos(v))) - 
   (e*(-1 + Power(Z,2))*(2 - p + 2*e*Cos(v))*((Power(e,2)*(56 - 13*p) + 4*(12 - 8*p + Power(p,2)))*Cos(v) + 
        e*(52 + 5*Power(e,2) - 34*p + 4*Power(p,2) + 2*(30 + 4*Power(e,2) - 7*p)*Cos(2*v) - 3*e*(-8 + p)*Cos(3*v) + 3*Power(e,2)*Cos(4*v))))/
    (2.*(-4*Power(e,2) + Power(-2 + p,2))*Power(p,2));
}
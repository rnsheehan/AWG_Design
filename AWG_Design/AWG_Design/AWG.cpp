#ifndef ATTACH_H
#include "Attach.h"
#endif

// Definitions for the AWG_Params class
// R. Sheehan 16 - 11 - 2018

AWG_Params::AWG_Params()
{
	// Default constructor for the params class
	m = N_aw = N_ch = N_in = 0;
	params_defined = false;
	dev_name = empty_string;
	lambda_c = delta_lambda_ch = delta_lambda_fsr = n_c = n_s = n_g = d = delta_x = W = L_min = 0.0;
	delta_L = L_f = L_f_2 = grat_disp_fac = awg_disp_fac = two_theta = one_theta = aperture_arc_length = image_plane_arc_length = 0.0;
	theta_n = theta_a = theta_j = d_theta_aw = d_theta_io = S_n = R_min = slab_sep = half_slab_sep = delta_L_2 = 0.0;
	We = Wgff = Lo = Lp = Lu = BW = 0.0;
}

AWG_Params::AWG_Params(double wl_centre, double wl_ch_spac, double wl_fsr, double wg_eff_indx, double wg_grp_indx, double slab_indx, double grating_pitch, double wg_width)
{
	// Primary Constructor for the params class
	// Inputs
	// wl_centre = central wavelength of AWG, units of um
	// wl_ch_spac = wavelength channel spacing, units of um
	// wl_fsr = wavelength free spectral range, units of um
	// wg_eff_indx = effective index of waveguide in array, this should be wavelength dependent
	// wg_grp_indx = group index of waveguide in array, this should be wavelength dependent
	// slab_indx = effective index of slab region, this should be wavelength dependent
	// grating_pitch = waveguide separation in array, should be large to prevent coupling between arrayed waveguides, units of um
	// output waveguide spacing = grating pitch, adopted custom

	set_params(wl_centre, wl_ch_spac, wl_fsr, wg_eff_indx, wg_grp_indx, slab_indx, grating_pitch, wg_width);
}

AWG_Params::AWG_Params(double wl_centre, double wl_ch_spac, double wl_fsr, double wg_eff_indx, double wg_grp_indx, double slab_indx, double grating_pitch, double output_pitch, double wg_width)
{
	// Primary Constructor for the params class
	// Inputs
	// wl_centre = central wavelength of AWG, units of um
	// wl_ch_spac = wavelength channel spacing, units of um
	// wl_fsr = wavelength free spectral range, units of um
	// wg_eff_indx = effective index of waveguide in array, this should be wavelength dependent
	// wg_grp_indx = group index of waveguide in array, this should be wavelength dependent
	// slab_indx = effective index of slab region, this should be wavelength dependent
	// grating_pitch = waveguide separation in array, should be large to prevent coupling between arrayed waveguides, units of um
	// output waveguide spacing = grating pitch, adopted custom

	set_params(wl_centre, wl_ch_spac, wl_fsr, wg_eff_indx, wg_grp_indx, slab_indx, grating_pitch, output_pitch, wg_width);
}

void AWG_Params::set_params(double wl_centre, double wl_ch_spac, double wl_fsr, double wg_eff_indx, double wg_grp_indx, double slab_indx, double grating_pitch, double wg_width)
{
	// method for setting all the parameters for the class
	// this will be the only set method since all the parameters are inter-related
	// a change to one parameter usually means a change to all
	// Inputs
	// wl_centre = central wavelength of AWG, units of um
	// wl_ch_spac = wavelength channel spacing, units of um
	// wl_fsr = wavelength free spectral range, units of um
	// wg_eff_indx = effective index of waveguide in array, this should be wavelength dependent
	// wg_grp_indx = group index of waveguide in array, this should be wavelength dependent
	// slab_indx = effective index of slab region, this should be wavelength dependent
	// grating_pitch = waveguide separation in array, should be large to prevent coupling between arrayed waveguides, units of um
	// output waveguide spacing = grating pitch, adopted custom

	// test inputs 
	bool c1 = (wl_centre > 0.0 && wl_centre < 10 ? true : false); // wl_centre in units of um?
	bool c2 = (wl_ch_spac > 0.0 ? true : false);
	bool c3 = (wl_fsr > 0.0 && wl_fsr > wl_ch_spac ? true : false);
	bool c4 = (wg_eff_indx > 0.0 ? true : false);
	bool c5 = (wg_grp_indx > 0.0 ? true : false);
	bool c6 = (slab_indx > 0.0 ? true : false);
	bool c7 = (grating_pitch > 1.0 ? true : false); // enforce minimum waveguide separation
	bool c8 = (wg_width > 0.0 ? true : false);
	bool c9 = (c1 && c2 && c3 && c4 && c5 && c6 && c7 && c8 ? true : false);

	// use exception handling to catch bad inputs
	try {
		if (c9) {
			// input parameters
			lambda_c = wl_centre;
			delta_lambda_ch = wl_ch_spac;
			delta_lambda_fsr = wl_fsr;
			n_c = wg_eff_indx;
			n_g = wg_grp_indx;
			n_s = slab_indx;
			d = grating_pitch;
			delta_x = grating_pitch;
			W = wg_width;

			N_in = 7; // specify a fixed number of input waveguides

					  // derived parameters
			set_disp_fac();
			set_m();
			set_Nchnnls();
			set_OPD();
			set_Lslab();
			set_diff_angle();
			set_Naw();
			set_arclengths();

			params_defined = true;
		}
		else {
			params_defined = false;
			std::string reason;
			reason = "Error: AWG_Params::set_params(double wl_centre, double wl_ch_spac, double wl_fsr, double wg_eff_indx, double wg_grp_indx, double slab_indx, double grating_pitch)\n";
			reason += "Input parameters not defined correctly\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		std::cerr << e.what();
	}
}

void AWG_Params::set_params(double wl_centre, double wl_ch_spac, double wl_fsr, double wg_eff_indx, double wg_grp_indx, double slab_indx, double grating_pitch, double output_pitch, double wg_width)
{
	// method for setting all the parameters for the class
	// this will be the only set method since all the parameters are inter-related
	// a change to one parameter usually means a change to all
	// Inputs
	// wl_centre = central wavelength of AWG, units of um
	// wl_ch_spac = wavelength channel spacing, units of um
	// wl_fsr = wavelength free spectral range, units of um
	// wg_eff_indx = effective index of waveguide in array, this should be wavelength dependent
	// wg_grp_indx = group index of waveguide in array, this should be wavelength dependent
	// slab_indx = effective index of slab region, this should be wavelength dependent
	// grating_pitch = waveguide separation in array, should be large to prevent coupling between arrayed waveguides, units of um
	// output_pitch = waveguide separation on image plane of AWG
	// output waveguide spacing = grating pitch, adopted custom

	// test inputs 
	bool c1 = (wl_centre > 0.0 && wl_centre < 10 ? true : false); // wl_centre in units of um?
	bool c2 = (wl_ch_spac > 0.0 ? true : false);
	bool c3 = (wl_fsr > 0.0 && wl_fsr > wl_ch_spac ? true : false);
	bool c4 = (wg_eff_indx > 0.0 ? true : false);
	bool c5 = (wg_grp_indx > 0.0 ? true : false);
	bool c6 = (slab_indx > 0.0 ? true : false);
	bool c7 = (grating_pitch > 1.0 ? true : false); // enforce minimum waveguide separation
	bool c8 = (wg_width > 0.0 ? true : false);
	bool c10 = (output_pitch > 0.0 ? true : false);
	bool c9 = (c1 && c2 && c3 && c4 && c5 && c6 && c7 && c8 && c10 ? true : false);

	// use exception handling to catch bad inputs
	try {
		if (c9) {
			// input parameters
			lambda_c = wl_centre;
			delta_lambda_ch = wl_ch_spac;
			delta_lambda_fsr = wl_fsr;
			n_c = wg_eff_indx;
			n_g = wg_grp_indx;
			n_s = slab_indx;
			d = grating_pitch;
			delta_x = output_pitch;
			W = wg_width;

			N_in = 7; // specify a fixed number of input waveguides

					  // derived parameters
			set_disp_fac();
			set_m();
			set_Nchnnls();
			set_OPD();
			set_Lslab();
			set_diff_angle();
			set_Naw();
			set_arclengths();

			params_defined = true;
		}
		else {
			params_defined = false;
			std::string reason;
			reason = "Error: AWG_Params::set_params(double wl_centre, double wl_ch_spac, double wl_fsr, double wg_eff_indx, double wg_grp_indx, double slab_indx, double grating_pitch)\n";
			reason += "Input parameters not defined correctly\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		std::cerr << e.what();
	}
}

void AWG_Params::set_slab_sep(double Lstr_min, double Rmin, bool loud)
{
	// Compute the separation between the slab sections based on knowledge of 	
	// Lstr_min = length of initial straight waveguide section that connects slab region to waveguide array, units of um
	// R_min = minimum bend radius of first waveguide in array, units of um
	// the angle theta is usually specified while drawing the AWG itself, however, theta = PI_4 is usually a good approximation
	// R. Sheehan 23 - 11 - 2018

	try {
		if (fabs(theta_n) > 0.0 && L_f > 0.0 && Lstr_min >= 0.0 && Rmin > 100) {
			R_min = Rmin; // this value must be updated if it does not satisfy the constraints
			S_n = Lstr_min; // choose S_n so that Tan(theta_j) \approx 1 => S_n = R_min - L_f
			//S_n = R_min - L_f; // this works in but not all cases? Why? 
			double ll = L_f + S_n;
			theta_a = atan(R_min / ll) - theta_n;
			theta_j = theta_a + theta_n;
			L_min = S_n + R_min * theta_j; // length of shortest arrayed waveguide
			slab_sep = 2.0*(R_min*sin(theta_j) + ll * cos(theta_j)); // distance between slab sections
			half_slab_sep = 0.5*slab_sep;

			if (loud) {
				std::cout << "Initial Layout Report\n";
				std::cout << "S_n: " << S_n << " um\n";
				std::cout << "R_n: " << R_min << " um\n";
				std::cout << "theta_a: " << theta_a << " rad, " << theta_a * (180.0 / PI) << " deg\n";
				std::cout << "theta_a + theta_n: " << theta_a + theta_n << " rad, " << (theta_a + theta_n) * (180.0 / PI) << " deg\n";
				std::cout << "L_min: " << L_min << " um\n";
				std::cout << "L_s: " << slab_sep << " um\n\n";
			}

			if (!test_Ls_valid()) {
				params_defined = false;
				std::string reason;
				reason = "Error: AWG_Params::set_slab_sep()\n";
				reason += "Input parameters not defined correctly\n";
				reason += "Need to increase the value of Rmin\n\n";
				throw std::invalid_argument(reason);
			}

			if (!test_constraint_satisfaction()) {
				params_defined = false;
				std::string reason;
				reason = "Error: AWG_Params::set_slab_sep()\n";
				reason += "Input parameters not defined correctly\n";
				reason += "Length contraint not satisfied\n\n";
				throw std::invalid_argument(reason);
			}
		}
		else {
			std::string reason;
			reason = "Error: AWG_Params::set_slab_sep()\n";
			reason += "Input parameters not defined correctly\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		std::cerr << e.what();
	}
}

double AWG_Params::get_L_j(int j)
{
	// Compute the half-length of the jth arrayed waveguide
	// R. Sheehan 29 - 11 - 2018
	try {
		if (j > -1 && L_min > 0.0 && delta_L_2 > 0.0) {
			return j > 0 ? L_min + j * delta_L_2 : L_min;
		}
		else {
			std::string reason;
			reason = "Error: AWG_Params::get_L_j()\n";
			reason += "Input parameters not defined correctly\n";
			throw std::invalid_argument(reason);

			return 0.0;
		}
	}
	catch (std::invalid_argument &e) {
		std::cerr << e.what();

		return 0.0;
	}
}

double AWG_Params::get_T_j(int j)
{
	// Compute the angle of jth arrayed waveguide along the arrayed waveguide arc
	// R. Sheehan 29 - 11 - 2018
	try {
		if (d_theta_aw > 0.0 && theta_j > 0.0) {
			return j > 0 ? theta_j + j * d_theta_aw : theta_j;
		}
		else {
			std::string reason;
			reason = "Error: AWG_Params::get_T_j()\n";
			reason += "Input parameters not defined correctly\n";
			throw std::invalid_argument(reason);

			return 0.0;
		}
	}
	catch (std::invalid_argument &e) {
		std::cerr << e.what();

		return 0.0;
	}
}

double AWG_Params::get_R_j(int j)
{
	// Compute the bend radius of the jth arrayed waveguide
	// R. Sheehan 29 - 11 - 2018
	try {
		if (j > -1 && slab_sep > 0.0 && L_f > 0.0) {
			double angle, cangle, sangle, numer, denom;

			angle = get_T_j(j);
			cangle = cos(angle);
			sangle = sin(angle);
			numer = half_slab_sep - (L_f + get_L_j(j))*cangle;
			denom = sangle - angle * cangle;

			return denom > 0.0 ? numer / denom : 0.0;
		}
		else {
			std::string reason;
			reason = "Error: AWG_Params::get_R_j()\n";
			reason += "Input parameters not defined correctly\n";
			throw std::invalid_argument(reason);

			return 0.0;
		}
	}
	catch (std::invalid_argument &e) {
		std::cerr << e.what();

		return 0.0;
	}
}

double AWG_Params::get_S_j(int j)
{
	// Compute the straight section of the jth arrayed waveguide
	// R. Sheehan 29 - 11 - 2018
	try {
		if (j > -1 && slab_sep > 0.0) {
			double S_j = get_L_j(j) - get_R_j(j) * get_T_j(j);

			return S_j;
		}
		else {
			std::string reason;
			reason = "Error: AWG_Params::get_S_j()\n";
			reason += "Input parameters not defined correctly\n";
			throw std::invalid_argument(reason);

			return 0.0;
		}
	}
	catch (std::invalid_argument &e) {
		std::cerr << e.what();

		return 0.0;
	}
}

double AWG_Params::get_F_j(int j)
{
	// Compute the angle of jth waveguide along the object plane / image plane
	// R. Sheehan 29 - 11 - 2018
	try {
		if (d_theta_aw > 0.0 && theta_j > 0.0) {
			return j > 0 ? theta_j + j * d_theta_io : theta_j;
		}
		else {
			std::string reason;
			reason = "Error: AWG_Params::get_T_j()\n";
			reason += "Input parameters not defined correctly\n";
			throw std::invalid_argument(reason);

			return 0.0;
		}
	}
	catch (std::invalid_argument &e) {
		std::cerr << e.what();

		return 0.0;
	}
}

bool AWG_Params::test_Ls_valid()
{
	// test if a computed Ls satisfies the minimum path-length condition
	// R. Sheehan 29 - 11 - 2018

	try {
		if (slab_sep > 0.0 && N_aw > 0) {
			int j = 0;
			bool valid = true;
			while (j < N_aw) {
				if (slab_sep < 2.0*(L_f + get_L_j(j)) * cos(get_T_j(j))) {
					valid = false;
					break;
				}
				j++;
			}
			return valid;
		}
		else {
			std::string reason;
			reason = "Error: AWG_Params::test_Ls_valid()\n";
			reason += "Input parameters not defined correctly\n";
			throw std::invalid_argument(reason);

			return false;
		}
	}
	catch (std::invalid_argument &e) {
		std::cerr << e.what();

		return false;
	}
}

bool AWG_Params::test_constraint_satisfaction()
{
	// Check all computed Sj, Rj to see if they satisfy the AWG length constraints
	// R. Sheehan 29 - 11 - 2018

	try {
		if (params_defined) {
			int j = 0;
			bool valid = true;
			double angle, ss, rr, ll, t1, t2, err;
			while (j < N_aw) {
				angle = get_T_j(j);
				ss = get_S_j(j);
				rr = get_R_j(j);
				ll = get_L_j(j);
				// First Constraint: Cannot have Sj < 0 or Sj > Lj or Rj < 0
				if (ss < 0.0 || ss > ll || rr < 0.0) {
					std::cout << "First Constraint Not Satisfied\n";
					std::cout << "j: " << j << " , theta_j: " << angle << " rad, S_j: " << ss << " um, R_j: " << rr << " um, L_j: " << ll << " um\n";
					valid = false;
					break;
				}
				t1 = (L_f + ss) * cos(angle);
				t2 = rr * sin(angle);
				err = fabs((t1 + t2) - half_slab_sep);
				// Second Constraint: Must have (Lf+Sj) CosTj + Rj Sin Tj = Ls/2
				if (err > 1.0e-3) {
					std::cout << "Second Constraint Not Satisfied\n";
					std::cout << "j: " << j << " , theta_j: " << angle << " rad, S_j: " << ss << " um, R_j: " << rr << " um, err: " << 1000 * err << " nm\n";
					valid = false;
					break;
				}
				// Third Constraint: Optical path length difference must be constant
				if (j > 0) {
					t1 = ss + rr * angle;
					t2 = get_S_j(j - 1) + get_R_j(j - 1)*get_T_j(j - 1);
					err = fabs(fabs(t1 - t2) - delta_L_2);
					if (err > 1.0e-3) {
						std::cout << "Third Constraint Not Satisfied\n";
						std::cout << "j: " << j << " , L_j: " << t1 << " um, L_j-1: " << t2 << " um, OPD: " << delta_L_2 << " um, err: " << 1000 * err << " nm\n";
						valid = false;
						break;
					}
				}
				j++;
			}
			return valid;
		}
		else {
			std::string reason;
			reason = "Error: AWG_Params::test_constraint_satisfaction()\n";
			reason += "Input parameters not defined correctly\n";
			throw std::invalid_argument(reason);

			return false;
		}
	}
	catch (std::invalid_argument &e) {
		std::cerr << e.what();

		return false;
	}
}

void AWG_Params::set_disp_fac()
{
	// compute the dispersion factors

	try {
		if (n_c > 0.0 && n_g > 0.0 && n_s > 0.0) {
			grat_disp_fac = (n_c / n_g);
			awg_disp_fac = n_s * grat_disp_fac;
		}
		else {
			std::string reason;
			reason = "Error: AWG_Params::set_disp_fac()\n";
			reason += "Input parameters not defined correctly\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		std::cerr << e.what();
	}
}

void AWG_Params::set_m()
{
	// compute the grating order

	try {
		if (grat_disp_fac > 0.0 && lambda_c > 0.0 && delta_lambda_fsr > 0.0) {
			m = static_cast<int>((lambda_c*grat_disp_fac) / delta_lambda_fsr);
		}
		else {
			std::string reason;
			reason = "Error: AWG_Params::set_disp_fac()\n";
			reason += "Input parameters not defined correctly\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		std::cerr << e.what();
	}
}

void AWG_Params::set_Nchnnls()
{
	// compute no. output channels
	// for the sake of simplicity in drawing the AWG I'm going to force N_ch to be odd, similarly for N_aw
	// this can be backed up by a physical argument
	// you want the central input / output waveguide to align with the central waveguide of the array, 
	// easiest way to ensure this is to have an odd number of wg

	try {
		if (delta_lambda_fsr > delta_lambda_ch && delta_lambda_ch > 0.0) {
			N_ch = static_cast<int>(delta_lambda_fsr / delta_lambda_ch);
			if (N_ch % 2 == 0) {
				N_ch += 1;
			}
		}
		else {
			std::string reason;
			reason = "Error: AWG_Params::set_Nchnnls()\n";
			reason += "Input parameters not defined correctly\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		std::cerr << e.what();
	}
}

void AWG_Params::set_OPD()
{
	// compute grating optical pathlength difference

	try {
		if (m > 0 && n_c > 0.0 && lambda_c > 0.0) {
			delta_L = (m*lambda_c) / n_c;
			delta_L_2 = 0.5*delta_L;
		}
		else {
			std::string reason;
			reason = "Error: AWG_Params::set_OPD()\n";
			reason += "Input parameters not defined correctly\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		std::cerr << e.what();
	}
}

void AWG_Params::set_Lslab()
{
	// compute slab focal length

	try {
		if (m > 0 && d > 0.0 && delta_x > 0.0 && delta_lambda_ch > 0.0 && awg_disp_fac > 0.0) {
			L_f = ((d*delta_x*awg_disp_fac) / (m*delta_lambda_ch));
			L_f_2 = 0.5 * L_f;
			d_theta_aw = d / L_f;
			d_theta_io = delta_x / L_f_2;
		}
		else {
			std::string reason;
			reason = "Error: AWG_Params::set_Lslab()\n";
			reason += "Input parameters not defined correctly\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		std::cerr << e.what();
	}
}

void AWG_Params::set_diff_angle()
{
	// compute slab diffraction angle

	try {
		if (n_s > 0.0 && lambda_c > 0.0 && d > 0.0) {
			two_theta = (4.0*lambda_c) / (n_s*d);
			one_theta = 0.5*two_theta;
		}
		else {
			std::string reason;
			reason = "Error: AWG_Params::set_diff_angle()\n";
			reason += "Input parameters not defined correctly\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		std::cerr << e.what();
	}
}

void AWG_Params::set_Naw()
{
	// set no. arrayed waveguides
	// for the sake of simplicity in drawing the AWG I'm going to force N_aw to be odd, similarly for N_ch
	// this can be backed up by a physical argument
	// you want the central input / output waveguide to align with the central waveguide of the array, 
	// easiest way to ensure this is to have an odd number of wg

	try {
		if (two_theta > 0.0 && W > 0.0 && d > 0.0 && L_f > 0.0 && d_theta_aw > 0.0) {
			N_aw = static_cast<int>((two_theta * L_f) / (W + d));
			if (N_aw % 2 == 0) {
				N_aw += 1;
			}
			theta_n = -0.5*(N_aw - 1)*d_theta_aw;
		}
		else {
			std::string reason;
			reason = "Error: AWG_Params::set_Naw()\n";
			reason += "Input parameters not defined correctly\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		std::cerr << e.what();
	}
}

void AWG_Params::set_arclengths()
{
	// set length of input/output to arrayed waveguide section, units of um
	// and length of output coupling section, units of um
	// aperture_arc_length should be consistent with value of two_theta*L_f

	try {
		if (N_aw > 0.0 && W > 0.0 && d > 0.0) {
			double t = (W + d);
			aperture_arc_length = N_aw * t;
			image_plane_arc_length = N_ch * t;
		}
		else {
			std::string reason;
			reason = "Error: AWG_Params::set_arclengths()\n";
			reason += "Input parameters not defined correctly\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		std::cerr << e.what();
	}
}

void AWG_Params::estimate_loss_BW(double ncore, double nsub, double X, bool loud)
{
	// estimate the nonuniformity loss, insertion loss and X dB bandwidth
	// R. Sheehan 3 - 12 - 2018

	try {
		if (ncore > nsub && nsub > 0.0) {
			set_Weff(ncore, nsub);
			set_Wgff();
			set_Lu();
			set_Lo();
			set_BW(X);

			if (loud) {
				std::cout << "Loss and Bandwidth Report\n";
				std::cout << "Recall high numerical aperture gives poor performance\n";
				std::cout << "Mode Effective Width: " << We << " um\n";
				std::cout << "Gaussian Far Field Width: " << Wgff << " rad\n";
				std::cout << "Nonuniformity Loss: " << Lu << " dB\n";
				std::cout << "Insertion Loss: " << Lo << " dB\n";
				std::cout << "BW: " << cnvrt_dNU_GHz_dWL_nm(BW, lambda_c) << " nm\n\n";
				//std::cout << "BW: " << cnvrt_dNU_GHz_dWL_nm(100, 1.55) << " nm\n\n";
			}
		}
		else {
			std::string reason;
			reason = "Error: AWG_Params::estimate_loss_BW()\n";
			reason += "Input parameters not defined correctly\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		std::cerr << e.what();
	}
}

void AWG_Params::set_Weff(double ncore, double nsub)
{
	// compute waveguide mode effective width
	// actually the mode half-width-half-max

	try {
		if (W > 0.0 && lambda_c > 0.0 && ncore > nsub && nsub > 0.0) {
			double k0 = (2.0*PI) / lambda_c;
			double V = 0.5 * k0 * W * sqrt((ncore*ncore) - (nsub*nsub));
			V -= (3.0 / 5.0);
			We = V > 0.0 ? W * (0.5 + (1.0 / V)) : 0.0;
		}
		else {
			std::string reason;
			reason = "Error: AWG_Params::set_Weff()\n";
			reason += "Input parameters not defined correctly\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		std::cerr << e.what();
	}
}

void AWG_Params::set_Wgff()
{
	// compute width of Gaussian far field

	try {
		if (We > 0.0 && lambda_c > 0.0 && n_s > 0.0) {
			Wgff = (lambda_c / n_s)*(1.0 / (We * sqrt(2.0*PI)));
		}
		else {
			std::string reason;
			reason = "Error: AWG_Params::set_Wgff()\n";
			reason += "Input parameters not defined correctly\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		std::cerr << e.what();
	}
}

void AWG_Params::set_Lu()
{
	// estimate the nonuniformity loss, units of dB

	try {
		if (one_theta > 0.0) {
			double ratio = one_theta / Wgff;
			ratio *= ratio;
			Lu = 8.7*ratio;
		}
		else {
			std::string reason;
			reason = "Error: AWG_Params::set_Lu()\n";
			reason += "Input parameters not defined correctly\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		std::cerr << e.what();
	}
}

void AWG_Params::set_Lo()
{
	// estimate the insertion loss, units of dB

	try {
		if (We > 0.0 && d > 0.0) {
			double ratio = We / d;
			ratio *= ratio;
			Lp = 2; // total propagation loss due to absorption and scattering
			Lo = Lp + 17.0*exp(-4.0*PI*ratio);
		}
		else {
			std::string reason;
			reason = "Error: AWG_Params::set_Lo()\n";
			reason += "Input parameters not defined correctly\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		std::cerr << e.what();
	}
}

void AWG_Params::set_BW(double X)
{
	// estimate the insertion loss, units of dB

	try {
		if (X > 0.0 && We > 0.0 && delta_x > 0.0 && delta_lambda_ch > 0.0) {
			double delta_nu = cnvrt_dWL_nm_dNU_GHz(delta_lambda_ch / 1000.0, lambda_c);
			BW = X == 1.0 ? 0.77 * (We / delta_x) * delta_nu : 0.77 * (We / delta_x) * delta_nu * sqrt(X);
		}
		else {
			std::string reason;
			reason = "Error: AWG_Params::set_Lo()\n";
			reason += "Input parameters not defined correctly\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		std::cerr << e.what();
	}
}

double AWG_Params::cnvrt_WL_um_Nu_THz(double lambda)
{
	// convert wavelength in units of um to Frequency in units of THz
	// recall 1.55 um = 193.548 THz

	try {
		if (lambda > 0.0) {
			return speedLight / lambda;
		}
		else {
			std::string reason;
			reason = "Error: AWG_Params::cnvrt_WL_um_Nu_THz()\n";
			reason += "Input parameters not defined correctly\n";
			throw std::invalid_argument(reason);

			return 0.0;
		}
	}
	catch (std::invalid_argument &e) {
		std::cerr << e.what();

		return 0.0;
	}
}

double AWG_Params::cnvrt_dNU_GHz_dWL_nm(double delta_nu_GHz, double lambda_um)
{
	// convert delta frequency Ghz to delta wavelength um
	// recall 100 GHz = 0.8 nm

	try {
		if (lambda_um > 0.0) {
			double nu_THz = cnvrt_WL_um_Nu_THz(lambda_um);
			double delta_lambda_nm = lambda_um * (delta_nu_GHz / nu_THz);
			//double delta_lambda_um = delta_lambda_nm / 1000.0; 

			return delta_lambda_nm;
		}
		else {
			std::string reason;
			reason = "Error: AWG_Params::cnvrt_dNU_GHz_dWL_nm()\n";
			reason += "Input parameters not defined correctly\n";
			throw std::invalid_argument(reason);

			return 0.0;
		}
	}
	catch (std::invalid_argument &e) {
		std::cerr << e.what();

		return 0.0;
	}
}

double AWG_Params::cnvrt_dWL_nm_dNU_GHz(double delta_WL_nm, double lambda_um)
{
	// convert delta wavelength um to delta frequency GHz
	// recall 0.8 nm = 100 GHz

	try {
		if (lambda_um > 0.0 && delta_WL_nm > 0.0) {
			double nu_THz = cnvrt_WL_um_Nu_THz(lambda_um);
			double delta_Nu_GHz = nu_THz * (delta_WL_nm / lambda_um);

			return delta_Nu_GHz;
		}
		else {
			std::string reason;
			reason = "Error: AWG_Params::cnvrt_dWL_nm_dNU_GHz()\n";
			reason += "Input parameters not defined correctly\n";
			throw std::invalid_argument(reason);

			return 0.0;
		}
	}
	catch (std::invalid_argument &e) {
		std::cerr << e.what();

		return 0.0;
	}
}

void AWG_Params::report()
{
	// print the class parameters values to the screen
	std::cout << "Input Parameters\n";
	std::cout << "Central Wavelength: " << lambda_c << " um\n";
	std::cout << "Channel Spacing: " << delta_lambda_ch << " um\n";
	std::cout << "AWG FSR: " << delta_lambda_fsr << "\n";
	std::cout << "WG neff:" << n_c << "\n";
	std::cout << "WG ngrp: " << n_g << "\n";
	std::cout << "Slab neff: " << n_s << "\n";
	std::cout << "Grating pitch: " << d << " um\n";
	std::cout << "WG width: " << W << " um\n";
	std::cout << "\nDerived Parameters\n";
	std::cout << "Grating order: " << m << "\n";
	std::cout << "No. output channels: " << N_ch << "\n";
	std::cout << "No. arrayed waveguides: " << N_aw << "\n";
	std::cout << "Optical pathlength difference: " << delta_L << " um \n";
	std::cout << "Slab focal length: " << L_f << "um \n";
	std::cout << "AWG dispersion factor: " << awg_disp_fac << "\n";
	std::cout << "Input / Output Angular Separation: " << d_theta_io << " rad\n";
	std::cout << "Arrayed Waveguide Angular Separation: " << d_theta_aw << " rad \n";
	std::cout << "Image Plane Arc Length: " << image_plane_arc_length << " um\n";
	std::cout << "Aperture Arc Length: " << aperture_arc_length << " um\n";
	std::cout << "WG angle relative to horizontal theta_n: " << theta_n << " rad = " << theta_n * (180.0 / PI) << " deg\n";
	std::cout << "\n";
}

void AWG_Params::set_name(std::string theName)
{
	try {
		if (theName != empty_string) {
			dev_name = theName;
		}
		else {
			std::string reason;
			reason = "Error: AWG_Params::set_name()\n";
			reason += "Input parameters not defined correctly\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		std::cerr << e.what();
	}
}
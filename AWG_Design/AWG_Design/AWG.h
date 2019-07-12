#ifndef AWG_H
#define AWG_H

// Need the ability to draw AWG using existing PICDraw structure
// Going to define a class that controls all the AWG params
// Will attempt to follow the AWG layout described in Okamoto
// For terminology consult the paper Smit and van Dam, ``PHASAR-Based WDM Devices: Principles, Design and Applications'', STQE, 2(2), 1996
// R. Sheehan 16 - 11 - 2018

class AWG_Params {
public:
	AWG_Params();
	AWG_Params(double wl_centre, double wl_ch_spac, double wl_fsr, double wg_eff_indx, double wg_grp_indx, double slab_indx, double grating_pitch, double wg_width);
	AWG_Params(double wl_centre, double wl_ch_spac, double wl_fsr, double wg_eff_indx, double wg_grp_indx, double slab_indx, double grating_pitch, double output_pitch, double wg_width);

	// Setters
	void set_params(double wl_centre, double wl_ch_spac, double wl_fsr, double wg_eff_indx, double wg_grp_indx, double slab_indx, double grating_pitch, double wg_width);
	void set_params(double wl_centre, double wl_ch_spac, double wl_fsr, double wg_eff_indx, double wg_grp_indx, double slab_indx, double grating_pitch, double output_pitch, double wg_width);
	void set_slab_sep(double Lstr_min, double R_min, bool loud = false);
	//void set_smin_slab_sep(); 
	void estimate_loss_BW(double ncore, double nsub, double X = 1.0, bool loud = false);
	void set_name(std::string theName);
	void report();

	// Getters
	inline int get_Nchnnls() { return N_ch; }
	inline int get_Naw() { return N_aw; }
	inline int get_Nin() { return N_in; }

	inline bool status() { return params_defined; }

	inline double get_OPD() { return delta_L; }
	inline double get_half_OPD() { return delta_L_2; }
	inline double get_Lslab() { return L_f; }
	inline double get_half_Lslab() { return L_f_2; }
	inline double get_pitch() { return d; }
	inline double get_W() { return W; }
	inline double get_2T() { return two_theta; }
	inline double get_T() { return one_theta; }
	inline double get_ss() { return slab_sep; }
	inline double get_half_ss() { return half_slab_sep; }
	inline double get_d_theta_aw() { return d_theta_aw; }
	inline double get_d_theta_io() { return d_theta_io; }
	inline double get_theta_n() { return theta_n; }
	inline double get_theta_a() { return theta_a; }

	inline std::string get_name() { return dev_name != empty_string ? dev_name : "NULL"; }

	double get_L_j(int j);
	double get_T_j(int j);
	double get_R_j(int j);
	double get_S_j(int j);
	double get_F_j(int j);

private:
	bool test_Ls_valid(); // check computed Ls to see if it si valid
	bool test_constraint_satisfaction(); // check all computed Sj, Rj to see if they are valid

	void set_disp_fac(); // compute the dispersion factors
	void set_m(); // compute the grating order
	void set_Nchnnls(); // compute no. output channels
	void set_OPD(); // compute grating optical pathlength difference
	void set_Lslab(); // compute slab focal length
	void set_diff_angle(); // compute slab diffraction angle
	void set_Naw(); // set no. arrayed waveguides
	void set_arclengths(); // compute aperture and image plane arc lengths
	void set_Weff(double ncore, double nsub); // compute waveguide mode effective width
	void set_Wgff(); // compute width of Gaussian far field
	void set_Lu(); // estimate the nonuniformity loss
	void set_Lo(); // estimate the insertion loss
	void set_BW(double X = 1.0); // estimate the X dB Bw of the AWG

	double cnvrt_WL_um_Nu_THz(double lambda_um); // convert wavelength in units of um to Frequency in units of THz
	double cnvrt_dNU_GHz_dWL_nm(double delta_nu_GHz, double lambda_um); // convert delta frequency Ghz to delta wavelength nm
	double cnvrt_dWL_nm_dNU_GHz(double delta_WL_nm, double lambda_um); // convert delta wavelength nm to delta frequency GHz
private:
	// Input parameters
	double lambda_c; // central wavelength of AWG, units of um
	double delta_lambda_ch; // wavelength channel spacing, units of um
	double delta_lambda_fsr; // wavelength free spectral range, units of um
	double n_c; // effective index of waveguide in array, this should be wavelength dependent
	double n_g; // group index of waveguide in array, this should be wavelength dependent
	double n_s; // effective index of slab region, this should be wavelength dependent
	double d; // grating pitch or waveguide separation in array, should be large to prevent coupling between arrayed waveguides, units of um
	double delta_x; // output waveguide spacing = grating pitch, adopted custom
	double W; // ridge width of wg in array
	double S_n; // length of initial straight waveguide section that connects slab region to waveguide array, units of um
	double R_min; // minimum bend radius of first waveguide in array, units of um

	std::string dev_name; // name of the device being designed

	// Derived parameters
	int m; // grating order
	int N_in; // no. of input waveguides
	int N_ch; // no. of output channels
	int N_aw; // no. of arrayed waveguides

	bool params_defined; // boolean that will be true / false depending on whether or not parameters have been defined

	double delta_L; // optical pathlength difference for waveguides in array, units of um
	double delta_L_2; // half optical pathlength difference for waveguides in array, units of um
	double L_f; // focal length of slab region, units of um
	double L_f_2; // half focal length, units of um
	double grat_disp_fac; // combined effective index and group index values for grating, this should be wavelength dependent
	double awg_disp_fac; // combined slab index and grat_disp_fac which accounts for dispersion in AWG, this should be wavelength dependent
	double two_theta; // diffraction angle in slab region, units of rad
	double one_theta; // half diffraction angle in slab region, units of rad
	double aperture_arc_length; // length of input/output to arrayed waveguide section, units of um
	double image_plane_arc_length; // length of output coupling section, units of um
	double d_theta_aw; // angular separation between waveguides on arrayed waveguide section
	double d_theta_io; // angular separation between waveguides on input / output sections
	double We; // waveguide mode effective width, units of um, more like the mode half-width-half-max
	double Wgff; // width of Gaussian far-field in slab region, units of um
	double BW; // AWG bandwidth, units of dB
	double Lu; // nonuniformity loss, units of dB
	double Lo; // insertion loss, units of dB
	double Lp; // waveguide propagation loss, units of dB

	// Layout Parameters, following Okamoto, sect. 9.3.3
	double theta_n; // angle that first waveguide in array makes with horizontal, units of rad
	double theta_a; // angle through the slab section is rotated
	double theta_j; // angle that jth waveguide in array makes with horizontal, includes effect of rotation of slab
	double slab_sep; // distance that separates the slab sections, determined by choice of Sn, Rmin
	double half_slab_sep; // half the distance that separates the slab sections
	double L_min; // length of shortest arrayed waveguide
};

#endif
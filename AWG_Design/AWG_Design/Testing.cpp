#ifndef ATTACH_H
#include "Attach.h"
#endif

void testing::compute_AWG_Params()
{
	// check the AWG param calculator to see if it works
	// R. Sheehan 16 - 11 - 2018

	int Nch = 11;
	double wl_centre, wl_ch_spac, wl_fsr_spac, wg_eff_indx, wg_grp_indx, slab_indx, grating_pitch, wg_width, smin, rmin;

	// SiN AWG Parameters
	// Want 5 channels minimum
	// Channel spacing should be 100 GHz, 200 GHz, 300 GHz, 400 GHz
	wl_centre = 1.535; // units of um

	//wl_ch_spac = (0.785408 / 1000.0); // 100 GHz in units of um
	//wl_fsr_spac = ((7.85408) / 1000.0); // in units of um

	//wl_ch_spac = (1.57082 / 1000.0); // 200 GHz in units of um
	//wl_fsr_spac = ((15.7082) / 1000.0); // in units of um

	//wl_ch_spac = (2.35622 / 1000.0); // 300 GHz in units of um
	//wl_fsr_spac = ((23.5622) / 1000.0); // in units of um

	wl_ch_spac = (3.14163 / 1000.0); // 400 GHz in units of um
	wl_fsr_spac = ((31.4163) / 1000.0); // in units of um

	wg_eff_indx = 1.377;
	wg_grp_indx = 1.718;
	slab_indx = 1.6841;
	grating_pitch = 2.0;
	wg_width = 1.0;

	AWG_Params theAWG;

	theAWG.set_params(wl_centre, wl_ch_spac, wl_fsr_spac, wg_eff_indx, wg_grp_indx, slab_indx, grating_pitch, wg_width);

	theAWG.report();

	smin = 120;
	rmin = 200;
	theAWG.set_slab_sep(smin, rmin, true);

	std::cout << "\nNaw = " << theAWG.get_Naw() << " , Naw / 2 = " << theAWG.get_Naw() / 2 << "\n";

	double ncore = 3.5405, nsub = 1.7176, X = 3.0;
	theAWG.estimate_loss_BW(ncore, nsub, X, true);
}

void testing::REDFINCH_AWG_Params(int set_number)
{
	// Compute the known REDFINCH AWG Params
	// four different flavours were required
	// R. Sheehan 30 - 11 - 2018

	// Define the AWG input parameters and compute all other AWG parameters
	double wl_centre, wl_ch_spac, wl_fsr_spac, wg_eff_indx, wg_grp_indx, slab_indx, grating_pitch, wg_width, smin, rmin;
	std::string awg_name, dev_name;
	std::string home = "C:\\Users\\Robert\\Research\\CAPPA\\MASKS\\REDFINCH_AWG_Dec_2018\\";

	if (set_number == 1) {
		// Set 1, TE
		awg_name = "Set_1_TE"; // name the AWG according to the usual names for the parameter sets
		dev_name = "Set 1 TE";
		wl_centre = 6.45; // units of um
		wl_ch_spac = (5.0 / 1000.0); // units of um
		wl_fsr_spac = ((71.0 + 10.0) / 1000.0); // units of um
		wg_eff_indx = 3.35;
		wg_grp_indx = 3.7;
		slab_indx = 3.41;
		grating_pitch = 10.0; // units of um
		wg_width = 4.0;	// units of um
	}
	else if (set_number == 2) {
		awg_name = "Set_2_TE"; // name the AWG according to the usual names for the parameter sets
		dev_name = "Set 2 TE";
		wl_centre = 6.1;
		wl_ch_spac = (5.0 / 1000.0);
		wl_fsr_spac = ((115.0 + 10.0) / 1000.0);
		wg_eff_indx = 3.37;
		wg_grp_indx = 3.68;
		slab_indx = 3.42;
		grating_pitch = 10.0;
		wg_width = 4.0;	
	}
	else if (set_number == 3) {
		awg_name = "Set_1_TM"; // name the AWG according to the usual names for the parameter sets
		dev_name = "Set 1 TM";
		wl_centre = 6.45;
		wl_ch_spac = (5.0 / 1000.0);
		wl_fsr_spac = ((71.0 + 10.0) / 1000.0);
		wg_eff_indx = 3.30;
		wg_grp_indx = 3.80;
		slab_indx = 3.36;
		grating_pitch = 10.0;
		wg_width = 4.0;
	}
	else if (set_number == 4) {
		awg_name = "Set_2_TM"; // name the AWG according to the usual names for the parameter sets
		dev_name = "Set 2 TM";
		wl_centre = 6.1;
		wl_ch_spac = (5.0 / 1000.0);
		wl_fsr_spac = ((115.0 + 10.0) / 1000.0);
		wg_eff_indx = 3.32;
		wg_grp_indx = 3.76;
		slab_indx = 3.38;
		grating_pitch = 10.0;
		wg_width = 4.0;
	}
	else {
		// Set 1, TE
		awg_name = "Set_1_TE"; // name the AWG according to the usual names for the parameter sets
		dev_name = "Set 1 TE";
		wl_centre = 6.45;
		wl_ch_spac = (5.0 / 1000.0);
		wl_fsr_spac = ((71.0 + 10.0) / 1000.0);
		wg_eff_indx = 3.35;
		wg_grp_indx = 3.7;
		slab_indx = 3.41;
		grating_pitch = 10.0;
		wg_width = 4.0;
	}	

	AWG_Params theAWG(wl_centre, wl_ch_spac, wl_fsr_spac, wg_eff_indx, wg_grp_indx, slab_indx, grating_pitch, wg_width);

	theAWG.set_name(dev_name);

	theAWG.report(); // print computed AWG parameters to the screen

	smin = 20;
	rmin = 1500;
	theAWG.set_slab_sep(smin, rmin, false);
}
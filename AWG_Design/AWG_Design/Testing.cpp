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
	grating_pitch = 2.5;
	wg_width = 1.0;

	AWG_Params theAWG;

	theAWG.set_params(wl_centre, wl_ch_spac, wl_fsr_spac, wg_eff_indx, wg_grp_indx, slab_indx, grating_pitch, wg_width);

	theAWG.report();

	smin = 20;
	rmin = 500;
	theAWG.set_slab_sep(smin, rmin, true);

	std::cout << "\nNaw = " << theAWG.get_Naw() << " , Naw / 2 = " << theAWG.get_Naw() / 2 << "\n";

	double ncore = 3.5405, nsub = 1.7176, X = 3.0;
	theAWG.estimate_loss_BW(ncore, nsub, X, true);
}
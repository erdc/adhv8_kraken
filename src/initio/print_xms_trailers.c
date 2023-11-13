/* ADH Version 2.0.0 6-04 */
/* writes the header to a data set file */

#include "global_header.h"

// this should really be in shallow water structs ...

void print_xms_trailers(SMODEL * mod)
{

    if (mod->flag.SW_FLOW) {
        if (mod->flag.WIND && mod->file_output.wind == ON) {
            print_trailer(mod->io->fout_sw_winds.fp);
        }
        if (mod->flag.WAVE && mod->file_output.wave == ON) {
            print_trailer(mod->io->fout_sw_waves.fp);
        }
    }

    if (mod->flag.SW2_FLOW) {
        print_trailer(mod->io->fout_sw2_head.fp);
        print_trailer(mod->io->fout_sw2_old_head.fp);
        print_trailer(mod->io->fout_sw2_vel.fp);
        print_trailer(mod->io->fout_sw2_old_vel.fp);
        print_trailer(mod->io->fout_sw2_error.fp);
        print_trailer(mod->io->fout_sw2_error_hydro.fp);
    }

    if (mod->flag.SW3_FLOW) {
        // 3d grid files
        print_trailer(mod->io->fout_sw3_displacement.fp);
        print_trailer(mod->io->fout_sw3_vel.fp);
        print_trailer(mod->io->fout_sw3_error.fp);
        print_trailer(mod->io->fout_sw3_error_hydro.fp);
        if (mod->file_output.grid_speed == ON) print_trailer(mod->io->fout_sw3_grid_speed.fp);
        if (mod->file_output.pressure == ON) print_trailer(mod->io->fout_sw3_pressure.fp);
		if (mod->file_output.hyd_vis == ON) print_trailer(mod->io->fout_sw3_hyd_viscosity.fp); //GSAVANT
		if (mod->file_output.trn_dif == ON) print_trailer(mod->io->fout_sw3_trn_diffusivity.fp); //GSAVANT
        
        // 2d grid files
        print_trailer(mod->io->fout_sw3_depth.fp);
        if (mod->file_output.depth_avg_velocity == ON) print_trailer(mod->io->fout_sw3_depth_avg_vel.fp);
        if (mod->file_output.surface_velocity == ON) print_trailer(mod->io->fout_sw3_surface_vel.fp);
        if (mod->file_output.bed_velocity == ON) print_trailer(mod->io->fout_sw3_bottom_vel.fp);
    }

    if (mod->flag.TRANSPORT) {
        int itrns = 0;
        for (itrns=0; itrns<mod->ntransport; itrns++) {
             print_trailer(mod->io->fout_con[itrns].fp);
             print_trailer(mod->io->fout_error_con[itrns].fp);
        }
    }

#ifdef _SEDIMENT
    if (mod->flag.SEDIMENT) {
        print_trailer(mod->io->fout_bed.fp);
        print_trailer(mod->io->fout_bed_flux.fp);
        print_trailer(mod->io->fout_active_layer.fp);
        int ilayer = 0;
        for (ilayer=0; ilayer<mod->sed->nlayers; ilayer++) {
            print_trailer(mod->io->fout_bed_layer[ilayer].fp);
        }
        int ised = 0;
        for (ised=0; ised<mod->sed->nsed; ised++) {
            print_trailer(mod->io->fout_bl_grain[ised].fp);
            print_trailer(mod->io->fout_sl_grain[ised].fp);
        }
    }
#endif

#ifdef _ADH_GROUNDWATER
    if (mod->flag.GW_FLOW) {
        print_trailer(mod->io->fout_gw_phead.fp);
        print_trailer(mod->io->fout_gw_thead.fp);
        print_trailer(mod->io->fout_gw_density.fp);
        print_trailer(mod->io->fout_gw_sat.fp);
        print_trailer(mod->io->fout_gw_flx.fp);
        print_trailer(mod->io->fout_gw_vel.fp);
        print_trailer(mod->io->fout_gw_error.fp);
 
    }
#endif
}

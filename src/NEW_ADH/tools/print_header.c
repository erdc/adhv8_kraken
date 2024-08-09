//  print_header.c
//  Created by Corey Trahan on 10/17/13.
#include "adh.h"

/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

void pdata(char *proj_name, char *run_name, FILE *fp_out, int o_flag, char *title, char *initials, char *begtype, char *meshdim, int lnodes, int lnelems)
{
    
    fprintf(fp_out, "OBJTYPE \"%s\"\n",meshdim);
    fprintf(fp_out, "%s\n",begtype);
    fprintf(fp_out, "ND %d\n", lnodes);
    fprintf(fp_out, "NC %d\n", lnelems);
    
    if (o_flag == 0) {
        fprintf(fp_out, "NAME \"%s%s\"\n", title, run_name);
    }
    else {
        fprintf(fp_out, "NAME \"%s_%s\"\n", proj_name, initials);
    }
}

/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

void pdata_index(SIO *io, FILE *fp_out, int o_flag, char *title, char *initials, char *begtype, char *meshdim, int lnodes, int lnelems, int index)
{

    fprintf(fp_out, "OBJTYPE \"%s\"\n",meshdim);
    fprintf(fp_out, "%s\n",begtype);
    fprintf(fp_out, "ND %d\n", lnodes);
    fprintf(fp_out, "NC %d\n", lnelems);

    if (o_flag == 0) {
        fprintf(fp_out, "NAME \"%s %d%s\"\n", title, index+1, io->run_name);
    }
    else {
        fprintf(fp_out, "NAME \"%s_%s %d\"\n", io->proj_name, initials, index+1);
    }
}

/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

void phot(FILE *fp_out, char *title, char *begtype, char *meshdim, int lnodes, int lelems)
{
    fprintf(fp_out, "OBJTYPE \"%s\"\n",meshdim);
    fprintf(fp_out, "%s\n",begtype);
    fprintf(fp_out, "ND %d\n", lnodes);
    fprintf(fp_out, "NC %d\n", lelems);
    fprintf(fp_out, "NAME \"%s\"\n", title);
}

/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

void phot_index(FILE *fp_out, char *title, char *begtype, char *meshdim, int lnodes, int lelems, int index)
{
    fprintf(fp_out, "OBJTYPE \"%s\"\n",meshdim);
    fprintf(fp_out, "%s\n",begtype);
    fprintf(fp_out, "ND %d\n", lnodes);
    fprintf(fp_out, "NC %d\n", lelems);
    fprintf(fp_out, "NAME \"%s %d\"\n", title, index+1);
}

/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

void print_header_sw2d(
        SMODEL_SUPER *mod,
        FILE * fp_out,   /* the output file */
        int ps_flag      /* the data set flag (PS_FLAG_xxx) */
)
{
   switch (ps_flag) {
           
           /* xms data files */
       case PS_FLAG_DEPTH:
           pdata(mod->io->proj_name, mod->io->run_name, fp_out, mod->o_flag, "Depth", "DEP", "BEGSCL", "mesh2d", mod->grid->macro_nnodes, mod->grid->macro_nelems2d);
           break;
       case PS_FLAG_PDEPTH:
           pdata(mod->io->proj_name, mod->io->run_name, fp_out, mod->o_flag, "Old Depth", "PDEP", "BEGSCL", "mesh2d", mod->grid->macro_nnodes, mod->grid->macro_nelems2d);
           break;
       case PS_FLAG_OLVEL:
           pdata(mod->io->proj_name, mod->io->run_name, fp_out, mod->o_flag, "Depth-Averaged Velocity", "OLVEL", "BEGVEC", "mesh2d", mod->grid->macro_nnodes, mod->grid->macro_nelems2d);
           break;
       case PS_FLAG_POLVEL:
           pdata(mod->io->proj_name, mod->io->run_name, fp_out, mod->o_flag, "Old Depth-Averaged Velocity", "POLVEL", "BEGVEC", "mesh2d", mod->grid->macro_nnodes, mod->grid->macro_nelems2d);
           break;
       case PS_FLAG_ERR:
           pdata(mod->io->proj_name, mod->io->run_name, fp_out, mod->o_flag, "Error", "ERR", "BEGSCL", "mesh2d", mod->grid->macro_nnodes, mod->grid->macro_nelems2d);
           break;
       case PS_FLAG_ERR_HYDRO:
           pdata(mod->io->proj_name, mod->io->run_name, fp_out, mod->o_flag, "Hydro Error", "HYDRO_ERR", "BEGSCL", "mesh2d", mod->grid->macro_nnodes, mod->grid->macro_nelems2d);
           break;
       case PS_FLAG_SWAVE_HEIGHT:
           pdata(mod->io->proj_name, mod->io->run_name, fp_out, mod->o_flag, "Wind-Wave Heights", "SWAVE_HEIGHT", "BEGSCL", "mesh2d", mod->grid->macro_nnodes, mod->grid->macro_nelems2d);
           break;
       case PS_FLAG_WAVES:
           pdata(mod->io->proj_name, mod->io->run_name, fp_out, mod->o_flag, "Wave Stress", "WAVE_STRESS", "BEGVEC", "mesh2d", mod->grid->macro_nnodes, mod->grid->macro_nelems2d);
           break;
       case PS_FLAG_WINDS:
           pdata(mod->io->proj_name, mod->io->run_name, fp_out, mod->o_flag, "Wind Stress", "WIND_STRESS", "BEGVEC", "mesh2d", mod->grid->macro_nnodes, mod->grid->macro_nelems2d);
           break;
       case PS_FLAG_BED_ELV:
           pdata(mod->io->proj_name, mod->io->run_name, fp_out, mod->o_flag, "Bed Elevation", "BED_ELV", "BEGSCL", "mesh2d", mod->grid->macro_nnodes, mod->grid->macro_nelems2d);
           break;
       case PS_FLAG_CON:
           pdata_index(mod->io, fp_out, mod->o_flag, "Concentration", "CON", "BEGSCL", "mesh2d", mod->grid->macro_nnodes, mod->grid->macro_nelems2d, mod->itrns);
           break;
       case PS_FLAG_PCON:
           pdata_index(mod->io, fp_out, mod->o_flag, "Old Concentration", "PCON", "BEGSCL", "mesh2d", mod->grid->macro_nnodes, mod->grid->macro_nelems2d, mod->itrns);
           break;
	   case PS_FLAG_ERR_CON:
		   pdata_index(mod->io, fp_out, mod->o_flag, "Concentration Error", "Transport Error", "BEGSCL", "mesh2d", mod->grid->initial_nnodes, mod->grid->initial_nelems, mod->itrns);
           break;


           /* hotstart files */
       case PS_FLAG_DEP_HOT:
           phot(fp_out, "ioh", "BEGSCL", "mesh2d", mod->grid->macro_nnodes, mod->grid->macro_nelems2d);
           break;
       case PS_FLAG_PDP_HOT:
           phot(fp_out, "iph", "BEGSCL", "mesh2d", mod->grid->macro_nnodes, mod->grid->macro_nelems2d);
           break;
       case PS_FLAG_OVL_HOT:
           phot(fp_out, "iov", "BEGVEC", "mesh2d", mod->grid->macro_nnodes, mod->grid->macro_nelems2d);
           break;
       case PS_FLAG_POV_HOT:
           phot(fp_out, "ipv", "BEGVEC", "mesh2d", mod->grid->macro_nnodes, mod->grid->macro_nelems2d);
           break;
       case PS_FLAG_CON_HOT:
           phot_index(fp_out, "icon", "BEGSCL", "mesh2d", mod->grid->macro_nnodes, mod->grid->macro_nelems2d, mod->itrns);
           break;
       case PS_FLAG_PCON_HOT:
           phot_index(fp_out, "ipc", "BEGSCL", "mesh2d", mod->grid->macro_nnodes, mod->grid->macro_nelems2d, mod->itrns);
           break;
 
   }
}

/**************************************************************************/
/**************************************************************************/
/**************************************************************************/
/* This print assumes there is a 2d grid file associated with the 3d grid */
void print_header_sw3d(
        SMODEL_SUPER *mod,
        FILE * fp_out,   /* the output file */
        int ps_flag      /* the data set flag (PS_FLAG_xxx) */
)
{
    int nnodes = UNSET_INT;
    int nelems = UNSET_INT;
    
    
    /* 3d xms data files */
    nnodes = mod->grid->macro_nnodes;
    nelems = mod->grid->macro_nelems3d;
    switch (ps_flag) {
            
	    case PS_FLAG_HYV:
			pdata(mod->io->proj_name, mod->io->run_name, fp_out, mod->o_flag, "Eddy Viscosity", "DPL", "BEGSCL", "mesh3d", nnodes, nelems);
            break;
		case PS_FLAG_TRD:
			pdata(mod->io->proj_name, mod->io->run_name, fp_out, mod->o_flag, "Eddy Diffusivity", "DPL", "BEGSCL", "mesh3d", nnodes, nelems);
            break;
        case PS_FLAG_DPL:
            pdata(mod->io->proj_name, mod->io->run_name, fp_out, mod->o_flag, "Displacement", "DPL", "BEGSCL", "mesh3d", nnodes, nelems);
            break;
        case PS_FLAG_PDPL:
            pdata(mod->io->proj_name, mod->io->run_name, fp_out, mod->o_flag, "Old Displacement", "PDPL", "BEGSCL", "mesh3d", nnodes, nelems);
            break;
        case PS_FLAG_PRS:
            pdata(mod->io->proj_name, mod->io->run_name, fp_out, mod->o_flag, "Pressure", "PRS", "BEGSCL", "mesh3d", nnodes, nelems);
            break;
        case PS_FLAG_GSP:
            pdata(mod->io->proj_name, mod->io->run_name, fp_out, mod->o_flag, "Grid Speed", "GSP", "BEGSCL", "mesh3d", nnodes, nelems);
            break;
        case PS_FLAG_ERR:
            pdata(mod->io->proj_name, mod->io->run_name, fp_out, mod->o_flag, "Error", "ERR", "BEGSCL", "mesh3d", nnodes, nelems);
            break;
        case PS_FLAG_ERR_HYDRO:
            pdata(mod->io->proj_name, mod->io->run_name, fp_out, mod->o_flag, "Hydro Error", "HYDRO_ERR", "BEGSCL", "mesh3d", nnodes, nelems);
            break;
        case PS_FLAG_FLX:
            pdata(mod->io->proj_name, mod->io->run_name, fp_out, mod->o_flag, "Nodal Flux", "FLX", "BEGSCL", "mesh3d", nnodes, nelems);
            break;
        case PS_FLAG_VEL:
            pdata(mod->io->proj_name, mod->io->run_name, fp_out, mod->o_flag, "Velocity", "VEL", "BEGVEC", "mesh3d", nnodes, nelems);
            break;
        case PS_FLAG_PVEL:
            pdata(mod->io->proj_name, mod->io->run_name, fp_out, mod->o_flag, "Old Velocity", "PVEL", "BEGVEC", "mesh3d", nnodes, nelems);
            break;
        case PS_FLAG_CON:
            pdata_index(mod->io, fp_out, mod->o_flag, "Concentration", "CON", "BEGSCL", "mesh3d", nnodes, nelems, mod->itrns);
            break;
        case PS_FLAG_PCON:
            pdata_index(mod->io, fp_out, mod->o_flag, "Old Concentration", "PCON", "BEGSCL", "mesh3d", nnodes, nelems, mod->itrns);
            break;
		case PS_FLAG_ERR_CON:
            pdata_index(mod->io, fp_out, mod->o_flag, "Concentration Error ", "TRAN_ERR", "BEGSCL", "mesh3d", nnodes, nelems, mod->itrns);
            break;

        /* hotstart files */
        case PS_FLAG_PRS_HOT:
            phot(fp_out, "p", "BEGSCL", "mesh3d", nnodes, nelems);
            break;
        case PS_FLAG_GSP_HOT:
            phot(fp_out, "igs", "BEGSCL", "mesh3d", nnodes, nelems);
            break;
        case PS_FLAG_VEL_HOT:
            phot(fp_out, "iv", "BEGVEC", "mesh3d", nnodes, nelems);
            break;
            
        case PS_FLAG_CON_HOT:
            phot_index(fp_out, "icon", "BEGSCL", "mesh3d", nnodes, nelems, mod->itrns);
            break;
        case PS_FLAG_PCON_HOT:
            phot_index(fp_out, "ipc", "BEGSCL", "mesh3d", nnodes, nelems, mod->itrns);
            break;
            
        default:
            break;
            
    }
    
    /* 2d xms data files */
    nnodes = mod->grid->macro_nnodes_bed;
    nelems = mod->grid->macro_nelems2d_bed;
    switch (ps_flag) {
            
        case PS_FLAG_DEPTH:
            pdata(mod->io->proj_name, mod->io->run_name, fp_out, mod->o_flag, "Depth", "DEP", "BEGSCL", "mesh2d", nnodes, nelems);
            break;
        case PS_FLAG_SWAVE_HEIGHT:
            pdata(mod->io->proj_name, mod->io->run_name, fp_out, mod->o_flag, "Wind-Wave Heights", "SWAVE_HEIGHT", "BEGSCL", "mesh2d", nnodes, nelems);
            break;
        case PS_FLAG_WAVES:
            pdata(mod->io->proj_name, mod->io->run_name, fp_out, mod->o_flag, "Wind-Wave Forces", "SWAVE_FORCES", "BEGVEC", "mesh2d", nnodes, nelems);
            break;
        case PS_FLAG_WINDS:
            pdata(mod->io->proj_name, mod->io->run_name, fp_out, mod->o_flag, "Wind Forces", "SWIND_FORCES", "BEGVEC", "mesh2d", nnodes, nelems);
            break;
        case PS_FLAG_DAVG_VEL:
            pdata(mod->io->proj_name, mod->io->run_name, fp_out, mod->o_flag, "Depth-Averaged Velocities", "DA_VEL", "BEGVEC", "mesh2d", nnodes, nelems);
            break;
        case PS_FLAG_SURF_VEL:
            pdata(mod->io->proj_name, mod->io->run_name, fp_out, mod->o_flag, "Surface Velocities", "SURF_VEL", "BEGVEC", "mesh2d", nnodes, nelems);
            break;
        case PS_FLAG_BOTT_VEL:
            pdata(mod->io->proj_name, mod->io->run_name, fp_out, mod->o_flag, "BOTTOM Velocities", "BOTT_VEL", "BEGVEC", "mesh2d", nnodes, nelems);
            break;
    }
    
}

/**************************************************************************/
/**************************************************************************/
/**************************************************************************/
/* This print assumes there is a 2d grid file associated with the 3d grid */
void print_header_sediment(
        SMODEL_SUPER *mod,
        FILE *fp_out,   /* the output file */
        int ps_flag,     /* the data set flag (PS_FLAG_xxx) */
        int index
)
{
    int nnodes = 0;
    int nelems = 0;
    char *mesh;

    fprintf(fp_out, "DATASET\n");
    
    nnodes = mod->grid->macro_nnodes_bed;
    nelems = mod->grid->initial_nelems_bed;
    switch (ps_flag) {

        case PS_FLAG_BED:
            pdata(mod->io->proj_name, mod->io->run_name, fp_out, mod->o_flag, "Bed Variables", "BED", "BEGSCL", "mesh2d", nnodes, nelems);
            break;
        case PS_FLAG_BED_FLUX:
            pdata(mod->io->proj_name, mod->io->run_name, fp_out, mod->o_flag, "Bed Load Flux", "BED_FLUX", "BEGVEC", "mesh2d", nnodes, nelems);
            break;
        case PS_FLAG_ALAYER:
            pdata(mod->io->proj_name, mod->io->run_name, fp_out, mod->o_flag, "Active Layer Variables", "ALAYER", "BEGSCL", "mesh2d", nnodes, nelems);
            break;
        case PS_FLAG_BED_LAYER:
            pdata(mod->io->proj_name, mod->io->run_name, fp_out, mod->o_flag, "Bed Layer Variables", "BED_LAYER", "BEGSCL", "mesh2d", nnodes, nelems);
            break;
        case PS_FLAG_BL_GRAIN:
            pdata(mod->io->proj_name, mod->io->run_name, fp_out, mod->o_flag, "Bed Load Grain Variables", "BL_GRAIN", "BEGSCL", "mesh2d", nnodes, nelems);
            break;

    }

    nnodes = mod->grid->macro_nnodes;
    nelems = mod->grid->initial_nelems;
    switch (ps_flag) {
        case PS_FLAG_SL_GRAIN:
            pdata(mod->io->proj_name, mod->io->run_name, fp_out, mod->o_flag, "Suspended Load Grain Variables", "SL_GRAIN", "BEGSCL", "mesh2d", nnodes, nelems);
            break;
    }
    
    
    tc_timeunits(fp_out, mod->series_out->outfact);
}

/**************************************************************************/
/**************************************************************************/
/**************************************************************************/
#ifdef _ADH_GROUNDWATER
void print_header_gw(
        SMODEL_SUPER *mod,
        FILE * fp_out,   /* the output file */
        int ps_flag      /* the data set flag (PS_FLAG_xxx) */
)
{
    int nnodes = UNSET_INT;
    int nelems = UNSET_INT;
    
    
    /* 3d xms data files */
    nnodes = mod->grid->macro_nnodes;
    nelems = mod->grid->macro_nelems3d;
    switch (ps_flag) {
            
	    case PS_FLAG_PHEAD:
			pdata(mod->io->proj_name, mod->io->run_name, fp_out, mod->o_flag, "Pressure Head", "PHEAD", "BEGSCL", "mesh3d", nnodes, nelems);
            break;
        case PS_FLAG_THEAD:
            pdata(mod->io->proj_name, mod->io->run_name, fp_out, mod->o_flag, "Total Head", "THEAD", "BEGSCL", "mesh3d", nnodes, nelems);
            break;
        case PS_FLAG_DENS:
            pdata(mod->io->proj_name, mod->io->run_name, fp_out, mod->o_flag, "Density", "DENS", "BEGSCL", "mesh3d", nnodes, nelems);
            break;
        case PS_FLAG_SAT:
            pdata(mod->io->proj_name, mod->io->run_name, fp_out, mod->o_flag, "Saturation", "SAT", "BEGSCL", "mesh3d", nnodes, nelems);
            break;
        case PS_FLAG_ERR:
            pdata(mod->io->proj_name, mod->io->run_name, fp_out, mod->o_flag, "Error", "ERR", "BEGSCL", "mesh3d", nnodes, nelems);
            break;
        case PS_FLAG_FLX:
            pdata(mod->io->proj_name, mod->io->run_name, fp_out, mod->o_flag, "Nodal Flux", "FLX", "BEGSCL", "mesh3d", nnodes, nelems);
            break;
        case PS_FLAG_VEL:
            pdata(mod->io->proj_name, mod->io->run_name, fp_out, mod->o_flag, "Velocity", "VEL", "BEGVEC", "mesh3d", nnodes, nelems);
            break;
        case PS_FLAG_CON:
            pdata_index(mod->io, fp_out, mod->o_flag, "Concentration", "CON", "BEGSCL", "mesh3d", nnodes, nelems, mod->itrns);
            break;
        case PS_FLAG_PCON:
            pdata_index(mod->io, fp_out, mod->o_flag, "Old Concentration", "PCON", "BEGSCL", "mesh3d", nnodes, nelems, mod->itrns);
            break;
		case PS_FLAG_ERR_CON:
            pdata_index(mod->io, fp_out, mod->o_flag, "Concentration Error ", "TRAN_ERR", "BEGSCL", "mesh3d", nnodes, nelems, mod->itrns);
            break;

        /* hotstart files */
        case PS_FLAG_PRS_HOT:
            phot(fp_out, "p", "BEGSCL", "mesh3d", nnodes, nelems);
            break;
        case PS_FLAG_GSP_HOT:
            phot(fp_out, "igs", "BEGSCL", "mesh3d", nnodes, nelems);
            break;
        case PS_FLAG_VEL_HOT:
            phot(fp_out, "iv", "BEGVEC", "mesh3d", nnodes, nelems);
            break;
            
        case PS_FLAG_CON_HOT:
            phot_index(fp_out, "icon", "BEGSCL", "mesh3d", nnodes, nelems, mod->itrns);
            break;
        case PS_FLAG_PCON_HOT:
            phot_index(fp_out, "ipc", "BEGSCL", "mesh3d", nnodes, nelems, mod->itrns);
            break;
            
        default:
            break;
            
    }
    
    /* 2d xms data files */
    nnodes = mod->grid->macro_nnodes_bed;
    nelems = mod->grid->macro_nelems2d_bed;
    switch (ps_flag) {
        default:
            break;
    }    
}
#endif

/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

void print_header( SMODEL_SUPER *mod, FILE * fp_out, int ps_flag) {
 
    fprintf(fp_out, "DATASET\n");
    if (mod->flag.SW2_FLOW == 1) {
        print_header_sw2d(mod, fp_out, ps_flag);
        tc_timeunits(fp_out, mod->series_out->outfact);
    }
    if (mod->flag.SW3_FLOW == 1) {
        print_header_sw3d(mod, fp_out, ps_flag);
        tc_timeunits(fp_out, mod->series_out->outfact);
    }
#ifdef _ADH_GROUNDWATER
    if (mod->flag.GW_FLOW == 1) {
        print_header_gw(mod, fp_out, ps_flag);
        tc_timeunits(fp_out, mod->series_out->outfact);
    }
#endif

}

/**************************************************************************/
/**************************************************************************/
/**************************************************************************/


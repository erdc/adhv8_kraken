#include "global_header.h"

#ifdef _ADH_GROUNDWATER

void read_bc_GW_MP(SMODEL *mod, char *data) {
  char line[MAXLINE];           /* the input line */
  char *subdata = NULL;         /* the data after the second card   is read */
  char *subsubdata = NULL;
    
  int imat = 0, itrn = 0;
  SIO info = *(mod->io);    // alias
  int nmat = mod->nmat;                 // alias
  int ntrn = mod->ntransport; // alias
  SMAT mat;
  /** should only be called if GW detected **/
  assert(mod->flag.GW_FLOW);
  
  switch (parse_card(data, &subdata)) {
    case CARD_ML:
      imat = get_material_id(info, &subdata, nmat);
      mat = mod->mat[imat];
      mat.gw->max_lev = read_int_field(info, &subdata);
      if (mat.gw->max_lev > 0) mod->flag.GRID_ADAPTION = ON;
      break;
    case CARD_FRT:
      imat = get_material_id(info, &subdata, nmat); mat = mod->mat[imat];
      mat.gw->refine_tolerance = read_dbl_field(info, &subdata);
      if (mat.gw->refine_tolerance < SMALL)
	tl_error("Refinement tolerance is too small for the precision of the machine.");
      if (mat.gw->refine_tolerance == 0.0)
	tl_error("Refinement tolerance must be greater than zero.");
      break;
    case CARD_SS:
      imat = get_material_id(info, &subdata, nmat); mat = mod->mat[imat];
      mat.gw->s_s = read_dbl_field(info,&subdata);
      if (mat.gw->s_s < 0.0)
	tl_error("Specific storage must be non-negative.");
      break;
    case CARD_POR:
      imat = get_material_id(info, &subdata, nmat); mat = mod->mat[imat];
      mat.gw->porosity = read_dbl_field(info,&subdata);
      if (mat.gw->porosity < 0.0)
	tl_error("Porosity must be non-negative.");
      break;
    case CARD_FVS:
      imat = get_material_id(info, &subdata, nmat); mat = mod->mat[imat];
      mat.gw->water_vol = read_dbl_field(info,&subdata);
      if (mat.gw->water_vol < 0.0)
	tl_error("Water volume must be non-negative.");
      break;
    case CARD_K:
      imat = get_material_id(info, &subdata, nmat); mat = mod->mat[imat];
      mat.gw->k.xx = read_dbl_field(info,&subdata);
      mat.gw->k.yy = read_dbl_field(info,&subdata);
      mat.gw->k.zz = read_dbl_field(info,&subdata);
      mat.gw->k.xy = read_dbl_field(info,&subdata);
      mat.gw->k.xz = read_dbl_field(info,&subdata);
      mat.gw->k.yz = read_dbl_field(info,&subdata);
      if (mat.gw->k.xx + mat.gw->k.yy + mat.gw->k.zz < NOT_QUITE_SMALL) {
	tl_error("Hydraulic conductivity cannot be negative");
      }
      break;
    case CARD_RSD:
      imat = get_material_id(info, &subdata, nmat); mat = mod->mat[imat];
      mat.gw->residual_sat = read_dbl_field(info,&subdata);
      if (mat.gw->residual_sat < 0.0)
	tl_error("Residual sat must be non-negative.");
      break;
    case CARD_VGA:
      imat = get_material_id(info, &subdata, nmat); mat = mod->mat[imat];
      mat.gw->vangen_alpha = read_dbl_field(info,&subdata);
      if (mat.gw->vangen_alpha <= NOT_QUITE_SMALL)
	tl_error("Van Genuchten alpha too small");
      break;
    case CARD_VGP:
      imat = get_material_id(info, &subdata, nmat); mat = mod->mat[imat];
      mat.gw->vangen_max_cp = read_dbl_field(info,&subdata);
      break;
    case CARD_VGN:
      imat = get_material_id(info, &subdata, nmat); mat = mod->mat[imat];
      mat.gw->vangen_n = read_dbl_field(info,&subdata);
      if (mat.gw->vangen_n <= NOT_QUITE_SMALL)
	tl_error("Van Genuchten n too small");
      break;
    case CARD_VGX:
      imat = get_material_id(info, &subdata, nmat); mat = mod->mat[imat];
      mat.gw->vangen_num_xy = read_int_field(info,&subdata);
      if (mat.gw->vangen_n < 0)
	tl_error("Van Genuchten number of entries in xy series should be non-negative");
      break;
    case CARD_BCL:
      imat = get_material_id(info, &subdata, nmat); mat = mod->mat[imat];
      mat.gw->brooks_lambda = read_dbl_field(info,&subdata);
      if (mat.gw->brooks_lambda < 0)
	tl_error("Brooks Corey Lambda should be non-negative");
      break;
    case CARD_BCE:
      imat = get_material_id(info, &subdata, nmat); mat = mod->mat[imat];
      mat.gw->brooks_pd = read_dbl_field(info,&subdata);
      break;
    case CARD_BCP:
      imat = get_material_id(info, &subdata, nmat); mat = mod->mat[imat];
      mat.gw->brooks_max_cp = read_dbl_field(info,&subdata);
      break;
    case CARD_BCX:
      imat = get_material_id(info, &subdata, nmat); mat = mod->mat[imat];
      mat.gw->brooks_num_xy = read_int_field(info,&subdata);
      if (mat.gw->brooks_num_xy < 0)
	tl_error("Brooks Corey number of entries in series should be non-negative");
      break;
    case CARD_KR:
      imat = get_material_id(info, &subdata, nmat); mat = mod->mat[imat];
      mat.gw->ikr = read_int_field(info,&subdata);
      if (mat.gw->ikr < 0)
	tl_error("relative permeability series index should be non-negative");
      break;
    case CARD_SAT:
      imat = get_material_id(info, &subdata, nmat); mat = mod->mat[imat];
      mat.gw->isat = read_int_field(info,&subdata);
      if (mat.gw->isat < 0)
	tl_error("capillary pressure series index should be non-negative");
      break;
    case CARD_TOR:
      imat = get_material_id(info, &subdata, nmat); mat = mod->mat[imat];
      mat.gw->tortuosity = read_dbl_field(info,&subdata);
      if (mat.gw->tortuosity < 0)
	tl_error("Tortuosity should be non-negative");
      break;
    case CARD_DPL:
      imat = get_material_id(info, &subdata, nmat); mat = mod->mat[imat];
      mat.gw->d_l = read_dbl_field(info,&subdata);
      if (mat.gw->d_l < 0)
	tl_error("Longitudinal dispersivity should be non-negative");
      break;
    case CARD_DPT:
      imat = get_material_id(info, &subdata, nmat); mat = mod->mat[imat];
      mat.gw->d_t = read_dbl_field(info,&subdata);
      if (mat.gw->d_t < 0)
	tl_error("Transverse dispersivity should be non-negative");
      break;
    case CARD_SSA:
      imat = get_material_id(info, &subdata, nmat); mat = mod->mat[imat];
      mat.gw->ss_area = read_dbl_field(info,&subdata);
      if (mat.gw->ss_area < 0)
	tl_error("Specific surface area should be non-negative");
      break;
    case CARD_BUL:
      imat = get_material_id(info, &subdata, nmat); mat = mod->mat[imat];
      mat.gw->bulk_density = read_dbl_field(info,&subdata);
      if (mat.gw->bulk_density < 0)
	tl_error("Bulk density should be non-negative");
      break;
    default:
      break;
  }
}
#endif

#include "adh.h"

#ifdef _ADH_BREACH
#include "breach_type.h"
#endif

/********************************************************/
/********************************************************/

void initialize_input_data_type(INPUT_DATA * input)
{
    assert(input);    /* should exist */
    input->bc_flag = UNSET_INT;
    input->src_flag = UNSET_INT;
    input->isigma = UNSET_INT;
    input->iu_0 = UNSET_INT;
    input->iq = UNSET_INT;
    input->ici = UNSET_INT;
    input->ivx = UNSET_INT;
    input->ivy = UNSET_INT;
    input->ivz = UNSET_INT;
    input->ceq_0 = UNSET_FLT;
}

/********************************************************/
/********************************************************/

void initialize_friction_type(FRICTION * frc)
{
    
    assert(frc);  /* should exist */
    frc->manningsn = UNSET_FLT;
    frc->skfricoef = UNSET_FLT;
    frc->dragcoef = UNSET_FLT;
    frc->eqrheight = UNSET_FLT;
    frc->rheight = UNSET_FLT;
    frc->hghtstem = UNSET_FLT;
    frc->diamstem = UNSET_FLT;
    frc->densstem = UNSET_FLT;
    frc->icethick = UNSET_FLT;
    frc->icedense = UNSET_FLT;
    frc->bedrhght = UNSET_FLT;
    frc->icerhght = UNSET_FLT;
    frc->mng_flag = NO;
    frc->erh_flag = NO;
    frc->sav_flag = NO;
    frc->urv_flag = NO;
    frc->icemors = UNSET_INT;
}

/********************************************************/
/********************************************************/
#ifdef _ADH_BREACH
void initialize_breach_info_type(BREACH_INFO * brc)
{
    assert(brc);  /* should exist */
    brc->form = UNSET_BREACH_FORM;
    brc->dpl_brch = UNSET_INT;
    brc->type = UNSET_BREACH_TYPE;
    brc->erodibility = UNSET_ERODE_TYPE;
    brc->first_node = UNSET_INT;
    brc->second_node = UNSET_INT;
    brc->fail_time = UNSET_FLT;
    brc->dam_ele = -DBL_MAX;
    brc->width = UNSET_FLT;
    brc->min_ele = -DBL_MAX;
    brc->exp = 1.;
    brc->max_depth = UNSET_FLT;
    brc->res_vol = UNSET_FLT;
}
#endif
/********************************************************/
/********************************************************/

void sstr_value_init(STR_VALUE **string, int nstring, int ntransport, int nsed) {
    
    
    int istring;
    (*string) = (STR_VALUE *) tl_alloc(sizeof(STR_VALUE), nstring);
    STR_VALUE *str = *(string);   // alias
    
    for (istring = 0; istring < nstring; istring++) {
        initialize_input_data_type(&(str[istring].flow));
        initialize_input_data_type(&(str[istring].heat));
        initialize_input_data_type(&(str[istring].ol_flow));
        initialize_input_data_type(&(str[istring].ch_flow));
        initialize_input_data_type(&(str[istring].pressure));
        initialize_input_data_type(&(str[istring].displacement));
        initialize_input_data_type(&(str[istring].bed));
        initialize_input_data_type(&(str[istring].Qs));
        initialize_input_data_type(&(str[istring].Td));
        initialize_friction_type(&(str[istring].fterms));
        str[istring].roughness = 0.;
        str[istring].conveyance = 0.;
        str[istring].height = UNSET_FLT;
        str[istring].ref_height = UNSET_FLT;
        str[istring].ref_pressure = DBL_MIN;
        str[istring].ps_flag = UNSET_INT;
        str[istring].string_type = UNSET_INT;
        str[istring].ice_string = UNSET_INT;
        str[istring].phys_flag = UNSET_INT;
        str[istring].flux_flag = FALSE;
        str[istring].sed_div_flag = FALSE;
        str[istring].link = UNSET_INT;
        str[istring].weir_num = UNSET_INT;
        str[istring].flap_num = UNSET_INT;
        str[istring].total_area = 0.0;
#ifdef _ADH_BREACH
        str[istring].breach = NULL;
#endif
        
        // transport
        str[istring].trans = NULL;
        if (ntransport > 0) {
            int itrns;
            str[istring].trans = (INPUT_DATA *) tl_alloc(sizeof(INPUT_DATA), ntransport);
            for (itrns=0; itrns<ntransport; itrns++) {
                initialize_input_data_type(&(str[istring].trans[itrns]));
            }
        }
        
        // sediment
        str[istring].sed = NULL;
        if (nsed > 0) {
            int ised;
            str[istring].sed = (INPUT_DATA *) tl_alloc(sizeof(INPUT_DATA), nsed);
            for (ised=0; ised<nsed; ised++) {
                initialize_input_data_type(&(str[istring].sed[ised]));
            }
        }
    }
}

/********************************************************/
/********************************************************/

//void sstr_value_print(STR_VALUE string) {

void sstr_value_free(STR_VALUE *str, int nstring, int ntransport, int nsed) {
    
    // transport
    int istring;
    for (istring = 0; istring < nstring; istring++) {
        if (ntransport > 0) {
            str[istring].trans = (INPUT_DATA *) tl_free(sizeof(INPUT_DATA), ntransport, str[istring].trans);
        }
    }
    
    // sediment
    for (istring = 0; istring < nstring; istring++) {
        if (nsed > 0) {
            str[istring].sed = (INPUT_DATA *) tl_free(sizeof(INPUT_DATA), nsed, str[istring].sed);
        }
    }
    
    str = (STR_VALUE *) tl_free(sizeof(STR_VALUE), nstring, str);
}


/********************************************************/
/********************************************************/

void sstr_value_check(STR_VALUE str, int istring)
{
    
    /* check node strings */
    if (str.string_type == STR_NODE) {
        
        if (str.phys_flag == OL_FLAG) {
            printf("\nERROR: STRING %d\n", istring + 1);
            printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
            tl_error(">> assigned overland flow to a node string.");
        }
        
        if (str.phys_flag == CH_FLAG) {
            printf("\nERROR: STRING %d\n", istring + 1);
            printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
            tl_error(">> assigned channel flow to a node string.");
        }
        
        if (str.flow.bc_flag == BCT_ROB ||
            str.flow.bc_flag == BCT_NEU ||
            str.flow.bc_flag == BCT_OUTFLOW ||
            str.flow.bc_flag == BCT_FLUX ||
            str.flow.bc_flag == BCT_NO_FLUX ||
            str.flow.bc_flag == BCT_PRS_NEU ||
            str.flow.bc_flag == BCT_VEL_NEU ||
            str.ol_flow.bc_flag == BCT_DIS_NEU || str.ol_flow.bc_flag == BCT_ROB || str.ol_flow.bc_flag == BCT_NEU || str.ol_flow.bc_flag == BCT_OUTFLOW || str.ol_flow.bc_flag == BCT_FLUX || str.ol_flow.bc_flag == BCT_NO_FLUX || str.ol_flow.bc_flag == BCT_PRS_NEU || str.ol_flow.bc_flag == BCT_VEL_NEU) {
            printf("\nERROR: STRING %d\n", istring + 1);
            printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
            tl_error(">> assigned flux boundary conditions for the flow on a node string.");
        }
        
    } else if (str.string_type == STR_EDGE) {
        
        if (str.flow.bc_flag == BCT_ROB || str.flow.bc_flag == BCT_NEU || str.flow.bc_flag == BCT_OUTFLOW || str.flow.bc_flag == BCT_FLUX || str.flow.bc_flag == BCT_NO_FLUX || str.flow.bc_flag == BCT_PRS_NEU || str.flow.bc_flag == BCT_VEL_NEU) {
            printf("\nERROR: STRING %d\n", istring + 1);
            printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
            tl_error(">> assigned flux boundary conditions for the flow on an edge string.");
        }
        
        if (str.ch_flow.bc_flag == BCT_ROB || str.ch_flow.bc_flag == BCT_NEU || str.ch_flow.bc_flag == BCT_OUTFLOW || str.ch_flow.bc_flag == BCT_FLUX || str.ch_flow.bc_flag == BCT_NO_FLUX || str.ch_flow.bc_flag == BCT_PRS_NEU || str.ch_flow.bc_flag == BCT_VEL_NEU) {
            printf("\nERROR: STRING %d\n", istring + 1);
            printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
            tl_error(">> assigned flux boundary conditions for channel flow on an edge string.");
        }
        
        if (str.flow.bc_flag == BCT_DIR ||
            str.flow.bc_flag == BCT_PRS_DIR || str.flow.bc_flag == BCT_VEL_DIR || str.ol_flow.bc_flag == BCT_DIR || str.ol_flow.bc_flag == BCT_PRS_DIR || str.ol_flow.bc_flag == BCT_VEL_DIR || str.ch_flow.bc_flag == BCT_DIR || str.ch_flow.bc_flag == BCT_PRS_DIR || str.ch_flow.bc_flag == BCT_VEL_DIR) {
            printf("\nERROR: STRING %d\n", istring + 1);
            printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
            tl_error(">> assigned Dirichlet boundary conditions on an edge string.");
        }
        
        if (str.ch_flow.bc_flag == BCT_ROB || str.ch_flow.bc_flag == BCT_NEU || str.ch_flow.bc_flag == BCT_OUTFLOW || str.ch_flow.bc_flag == BCT_FLUX || str.ch_flow.bc_flag == BCT_NO_FLUX || str.ch_flow.bc_flag == BCT_PRS_NEU || str.ch_flow.bc_flag == BCT_VEL_NEU) {
            printf("\nERROR: STRING %d\n", istring + 1);
            printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
            tl_error(">> assigned 1d flux bc to an edge string.\n");
        }
        
    } else if (str.string_type == STR_FACE) {
        
        if (str.phys_flag == CH_FLAG) {
            printf("\nERROR: STRING %d\n", istring + 1);
            printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
            tl_error(">> assigned channel flow to a face string.");
        }
        
        if (str.flow.bc_flag == BCT_DIR ||
            str.flow.bc_flag == BCT_PRS_DIR || str.flow.bc_flag == BCT_VEL_DIR || str.ol_flow.bc_flag == BCT_DIR || str.ol_flow.bc_flag == BCT_PRS_DIR || str.ol_flow.bc_flag == BCT_VEL_DIR || str.ch_flow.bc_flag == BCT_DIR || str.ch_flow.bc_flag == BCT_PRS_DIR || str.ch_flow.bc_flag == BCT_VEL_DIR) {
            printf("\nERROR: STRING %d\n", istring + 1);
            printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
            tl_error(">> assigned Dirichlet boundary conditions for the flow on a face string.");
        }
        if (str.ch_flow.bc_flag == BCT_ROB || str.ch_flow.bc_flag == BCT_NEU || str.ch_flow.bc_flag == BCT_OUTFLOW || str.ch_flow.bc_flag == BCT_FLUX || str.ch_flow.bc_flag == BCT_NO_FLUX || str.ch_flow.bc_flag == BCT_PRS_NEU || str.ch_flow.bc_flag == BCT_VEL_NEU) {
            printf("\nERROR: STRING %d\n", istring + 1);
            printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
            tl_error(">> assigned 1d flux bc to a face string.\n");
        }
        
        if (str.flow.src_flag == EXTRACTION_WELL || str.flow.src_flag == INJECTION_WELL) {
            printf("\nERROR: STRING %d\n", istring + 1);
            printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
            tl_error("Assigned a well to a face string.");
        }
        
    }  else {
        
        fprintf(stderr, "  String %d is either a face in a 3D run or a bad string type.\n", istring + 1);
        
    }
}

/********************************************************/
/********************************************************/

void sstr_value_checkall(STR_VALUE *str, int nstring) {
    int istring;
    for (istring = 0; istring < nstring; istring++) {
        sstr_value_check(str[istring], istring);
    }
}

/********************************************************/
/********************************************************/

void sstr_value_set_physics(STR_VALUE *str, int nstring, SFLAGS flag) {
    int istring = UNSET_INT;
    for (istring = 0; istring < nstring; istring++) {
        if (str[istring].phys_flag != OFF && str[istring].string_type != STR_MID) {
            if (flag.SW2_FLOW == ON){
                str[istring].phys_flag = SW2_FLAG;
                if (flag.DIFFUSIVE_WAVE == ON){
                    str[istring].phys_flag = OL_FLAG;
                }
            } else if (flag.SW3_FLOW == ON){
                str[istring].phys_flag = SW3_FLAG;
            } else if (flag.GW_FLOW == ON){
                str[istring].phys_flag = GW_FLAG;
            } else if (flag.DIFFUSIVE_WAVE == ON){
                str[istring].phys_flag = OL_FLAG;
            } else if (flag.NS_FLOW == ON){
                str[istring].phys_flag = NS_FLAG;
            }
        }
    }
}

/********************************************************/
/********************************************************/

void sstr_value_set_roughness(STR_VALUE *str, int nstring, int SW2_FLOW, int SW3_FLOW, double g, double muc) {
    int istring = UNSET_INT;
    for (istring = 0; istring < nstring; istring++) {
        if (SW2_FLOW || SW3_FLOW) {
            if (str[istring].fterms.mng_flag == YES && str[istring].fterms.erh_flag == NO) {
                str[istring].fterms.eqrheight = fr_manningsn_to_rheight(str[istring].fterms.manningsn, muc, g);
            } else if (str[istring].fterms.mng_flag == NO && str[istring].fterms.erh_flag == YES) {
                str[istring].fterms.manningsn = fr_rheight_to_manningsn(str[istring].fterms.eqrheight, g);
            }
        }
    }
}

/********************************************************/
/********************************************************/
// cjt :: should I add displacement to z here???  The extrusion code right now adds displacement to the grid and sets dpl = 0.;
// cjt :: Yes, but for now, since initial dpls are only given in test case preps, and those cases have zero inflow when dpl
//          is used, it's ok.  (djac_fixed is only used for inflows)

void sstr_value_set_total_area(STR_VALUE *str, SGRID *grid) {
    int ie = UNSET_INT;
    
    assert(grid->ndim == 3);
    
    for (ie = 0; ie < grid->nelems2d; ie++) {
        if (grid->elem2d[ie].string > NORMAL) {
            if (grid->elem2d[ie].nnodes == 3) {
                str[grid->elem2d[ie].string].total_area += grid->elem2d[ie].djac3d_fixed;
            } else {
                // use cross product to get area of quadrilateral (cannot use this for djac)
                int idof;
                SVECT nd[grid->elem2d[ie].nnodes];
                for (idof=0; idof<grid->elem2d[ie].nnodes; idof++) {
                    nd[idof].x = grid->node[grid->elem2d[ie].nodes[idof]].x;
                    nd[idof].y = grid->node[grid->elem2d[ie].nodes[idof]].y;
                    nd[idof].z = grid->node[grid->elem2d[ie].nodes[idof]].z;
                }
                SVECT v1; svect_subtract_array2(&v1, &nd[0], &nd[1], 1);
                SVECT v2; svect_subtract_array2(&v2, &nd[0], &nd[3], 1);
                SVECT w1; svect_subtract_array2(&w1, &nd[2], &nd[1], 1);
                SVECT w2; svect_subtract_array2(&w2, &nd[2], &nd[3], 1);
                str[grid->elem2d[ie].string].total_area += one_2 * (svect_mag(svect_cross(v1,v2)) + svect_mag(svect_cross(w1,w2)));
            }
        }
    }
}

/********************************************************/
/********************************************************/



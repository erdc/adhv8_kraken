#include "global_header.h"

/***********************************************************/
/***********************************************************/

/* ielem1d :: input :: counting total of 1d elements */
/* ielem2d :: input :: counting total of 2d elements */

void read_bc_STRINGS(SMODEL *mod, CARD card, char *data, int *ielem1d, int *ielem2d) {
    
    char line[MAXLINE];           /* the input line */
    char *subdata = NULL;         /* the data after the second card   is read */
    char *subsubdata = NULL;
    int i, j, k;
    int ind1=0, ind2=0, nnode=0, imat=0, istring=0;
    int local_index1, local_index2, local_elem;
    SIO info = *(mod->io);   //alias
    int nstring = mod->nstring;           // alias
    int nmat;                             // alias
    
    SGRID *grid = mod->grid;             // alias
    nmat = mod->nmat;
    int nnodes = grid->macro_nnodes; // alias
#ifdef _MESSG
    int myid = grid->smpi->myid;
    int npes =  grid->smpi->npes;
#else
    int myid = 0, npes = 1;
#endif
    int gid_1d=0;
    
    switch (card) {
            
        case CARD_NDS:
            ind1 = get_node_id(info, &data, nnodes);
            istring = get_string_id(info, &data, nstring);
            local_index1 = is_my_node(ind1, grid);
            if (mod->str_values[istring].string_type != UNSET_INT) {
                if (mod->str_values[istring].string_type != STR_NODE) {
                    printf("\n WARNING: NDString number is already used for string type %d\n",mod->str_values[istring].string_type);
                }
            }
            mod->str_values[istring].string_type = STR_NODE;
            if (local_index1 < 0){
                break;
            }
            else {
                ind1 = local_index1;
            }
            grid->node[ind1].string = istring;
            
            //nsnodes++;
            break;
            
        case CARD_EGS:
            gid_1d++;
            ind1 = get_node_id(info, &data, nnodes);
            ind2 = get_node_id(info, &data, nnodes);
            istring = get_string_id(info, &data, nstring);
            local_index1 = is_my_node(ind1, grid);
            local_index2 = is_my_node(ind2, grid);
            if (mod->str_values[istring].string_type != UNSET_INT) {
                if (mod->str_values[istring].string_type != STR_EDGE) {
                    printf("\n WARNING...String number is already used for a different string type.");
                }
            }
            mod->str_values[istring].string_type = STR_EDGE;
            if ((local_index1 < 0) || (local_index2 < 0) ){
                break;
            } else {  
                ind1 = local_index1;
                ind2 = local_index2;
            }
            local_elem = is_my_elem(ind1,ind2, grid);
            if(local_elem<0) break;
            if(*ielem1d>=grid->nelems1d) *ielem1d = elem1d_new(grid, mod->nalloc_inc);
            grid->elem1d[*ielem1d].id = *ielem1d;
            grid->elem1d[*ielem1d].string = istring;
            grid->elem1d[*ielem1d].mat = UNSET_INT;
            grid->elem1d[*ielem1d].nodes[0] = ind1;
            grid->elem1d[*ielem1d].nodes[1] = ind2;
            grid->node[ind1].edge_string = istring;
            grid->node[ind2].edge_string = istring;
            grid->elem1d[*ielem1d].gid=gid_1d;
            get_elem1d_linear_djac_gradPhi(grid, &grid->elem1d[*ielem1d]);
            //printf("myid %d ielem %d elem1d.id %d\n", grid->smpi->myid, *ielem1d, grid->elem1d[*ielem1d].id);
            (*ielem1d)++;
            break;
            
        case CARD_MDS:
            gid_1d++;
            ind1 = get_node_id(info, &data, nnodes);
            ind2 = get_node_id(info, &data, nnodes);
            istring = get_string_id(info, &data, nstring);
            local_index1 = is_my_node(ind1, grid);
            local_index2 = is_my_node(ind2,grid);
            if (mod->str_values[istring].string_type != UNSET_INT) {
                if (mod->str_values[istring].string_type != STR_EDGE) {
                    printf("\n WARNING...String number is already used for a different string type.");
                }
            }
            mod->str_values[istring].string_type = STR_MID;
            if ((local_index1 < 0) || (local_index2 < 0) ){
                break;
            }
            else {
                
                ind1 = local_index1;
                ind2 = local_index2;
            }
            local_elem = is_my_elem(ind1,ind2, grid);
            if(local_elem<0) break;
            if(*ielem1d>=grid->nelems1d) *ielem1d = elem1d_new(grid, mod->nalloc_inc);
            grid->elem1d[*ielem1d].id = *ielem1d;
            grid->elem1d[*ielem1d].string = istring;
            grid->elem1d[*ielem1d].mat = UNSET_INT;
            grid->elem1d[*ielem1d].nodes[0] = ind1;
            grid->elem1d[*ielem1d].nodes[1] = ind2;
            grid->node[ind1].edge_string = istring;
            grid->node[ind2].edge_string = istring;
            grid->elem1d[*ielem1d].gid=gid_1d;
            get_elem1d_linear_djac_gradPhi(grid, &grid->elem1d[*ielem1d]);
            (*ielem1d)++;
            break;
            
            // faces for 3d code now have their own files
            //case CARD_FCS:
            //  break;
            
        case CARD_MTS:
            imat = get_material_id(info, &data, nmat);
            istring = get_string_id(info, &data, nstring);
            if (grid->nelems2d <= 0) {
                tl_error("MTS is to assign material to 2d elements, but no 2d elements.");
            }
            if (mod->str_values[istring].string_type != UNSET_INT) {
                if (mod->str_values[istring].string_type != STR_FACE) {
                    printf("\n WARNING...String number is already used for a different string type.");
                }
            }
            mod->str_values[istring].string_type = STR_FACE;
            for (i = 0; i < grid->nelems2d; i++) {
                if (grid->elem2d[i].mat == imat) {
                    grid->elem2d[i].string = istring;
                }
            }
            break;
        default:
            break;
    }
    
}

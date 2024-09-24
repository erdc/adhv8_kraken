/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  read_bc_STRINGS.c  This file collects methods to read an AdH SuperModel boundary condition input file       */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
static int DEBUG = OFF;
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Reads a SuperModel boundary condition file for reading strings
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] *sm (SUPER_MODEL *)  a double pointer to an array of AdH supermodels
 * @param[in] **token (CHAR) a BC file line string token
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void read_bc_STRINGS(SMODEL_SUPER *mod, char **token) {
    
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
    
    if (strcmp(*token, "NDS") == 0) {
            ind1 = get_id(token,nnodes,"NDS NODE Error\n");
            istring = get_id(token,nstring,"NDS String Error\n");
            local_index1 = is_my_node(ind1, grid);
            if (mod->str_values[istring].string_type != UNSET_INT) {
                if (mod->str_values[istring].string_type != STR_NODE) {
                    printf("\n WARNING: NDString number is already used for string type %d\n",mod->str_values[istring].string_type);
                }
            }
            mod->str_values[istring].string_type = STR_NODE;
            if (local_index1 < 0){break;}
            else {ind1 = local_index1;}
            grid->node[ind1].string = istring;
            
    } else if (strcmp(*token, "EGS") == 0) {
            gid_1d++;
            ind1 = get_id(token,nnodes,"EGS NODE Error\n");
            ind2 = get_id(token,nnodes,"EGS NODE Error\n");
            istring = get_id(token,nstring,"MNG String Error\n");
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
            
    } else if (strcmp(*token, "MDS") == 0) {
            gid_1d++;
            ind1 = get_id(token,nnodes,"MDS NODE Error\n");
            ind2 = get_id(token,nnodes,"MDS NODE Error\n");
            istring = get_id(token,nstring,"MNG String Error\n");
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

    } else if (strcmp(*token, "MTS") == 0) {
            imat = get_id(token,nmat,"MNG MAT Error\n");
            istring = get_id(token,nstring,"MNG STRING Error\n");
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
    }
    
}

#include "global_header.h"

// note :: this file creates edge strings when DB OVL boundary conditions are used.
//      :: As far as I can ascertain, these edge strings are to enforce a weak 1d elem bc on the mass conservation
//      :: Without it, rhs.c_eq is 0 on the element edge.
// note :: adds 1d elements! (cjt)
// note :: all newly created edge strings (even for multiple NDS) will only need one new string.
//      :: This is because we are not going to apply a time-series to this string, just call sw_1d_convection

void create_dbovl_elem1d(SMODEL *mod) {

    int i, j, inode, ie;
    int istring = -1, jstring = -1;
    int nnodes = -1;
    int nd1 = -1, nd2 = -1, nd3 = -1;

    // this is just for SW2 applications
    if (mod->flag.SW2_FLOW != ON) return;

    // alias
    SGRID *grid = mod->grid;
    nnodes = grid->nnodes;

    //first runtime string number (for DB OVL)
    jstring = mod->nstring - 1;

    // allocated maximum array
    int *db_vel_flag = (int *) tl_alloc(sizeof(int), nnodes);
    for (inode = 0; inode < nnodes; inode++) {
        db_vel_flag[inode] = 0;
    }

    //find nodes associated with DB OVL
    int db_vel_count = 0;
    for (istring = 0; istring < mod->nstring; istring++) {
        if (mod->str_values[istring].ol_flow.bc_flag == BCT_VEL_DIR) {
            for (inode = 0; inode < nnodes; inode++) {
                if (grid->node[inode].string == istring) {
                  //db_vel_flag[db_vel_count] = inode;
                  db_vel_flag[inode] = 1;
                  db_vel_count++;
                }
            }
        }
    }

    // for every BCT_VEL_DIR string, create an additional string
    mod->str_values[jstring].ol_flow.bc_flag = BCT_DB_VEL;
    mod->str_values[jstring].string_type = STR_EDGE;
    if (mod->str_values[jstring].phys_flag != SW2_FLAG) {
        mod->str_values[jstring].phys_flag = SW2_FLAG;
        mod->str_values[jstring].roughness = 0.0;
    }

    // loop over 2d elems and find edge
    int found = 0;
    for (ie = 0; ie < grid->nelems2d; ie++) {

        nd1 = grid->elem2d[ie].nodes[0];
        nd2 = grid->elem2d[ie].nodes[1];
        nd3 = grid->elem2d[ie].nodes[2];

        found = 0;
        if (db_vel_flag[nd1] + db_vel_flag[nd2] + db_vel_flag[nd3] > 1) {

            // we found and edge (assuming db ovl are given as domain edges, this should not duplicate)
            grid->nelems1d++;
            grid->elem1d = (SELEM_1D *) tl_realloc(sizeof(SELEM_1D), grid->nelems1d, grid->nelems1d-1, grid->elem1d);
            selem1d_init(&grid->elem1d[grid->nelems1d-1]);
            
            if (db_vel_flag[nd1] == 1) {
                grid->elem1d[grid->nelems1d-1].nodes[found] = nd1;
                found++;
            }
            if (db_vel_flag[nd2] == 1) {
                grid->elem1d[grid->nelems1d-1].nodes[found] = nd2;
                found++;
            } 
            if (db_vel_flag[nd3] == 1) {
                grid->elem1d[grid->nelems1d-1].nodes[found] = nd3;
                found++;
            }

            nd1 = grid->elem1d[grid->nelems1d-1].nodes[0];
            nd2 = grid->elem1d[grid->nelems1d-1].nodes[1];
            get_elem1d_linear_djac_gradPhi(grid, &grid->elem1d[grid->nelems1d-1]);
            grid->elem1d[grid->nelems1d-1].id = grid->nelems1d-1;
            grid->elem1d[grid->nelems1d-1].string = jstring;
            grid->elem1d[grid->nelems1d-1].mat = UNSET_INT;
            grid->elem1d[grid->nelems1d-1].elem2d = ie;
        }
    }

    grid->max_nelems1d = grid->nelems1d;

    /*
    printf("db_vel_count: %d \n",db_vel_count);
    for (ie = 0; ie < grid->nelems1d; ie++) {
        selem1d_printScreen(grid->elem1d[ie]);
    }
    exit(-1);
    */

    db_vel_flag = (int *) tl_free(sizeof(int), nnodes, db_vel_flag);
}

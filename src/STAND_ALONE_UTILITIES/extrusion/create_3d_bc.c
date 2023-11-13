#include "extrusion.h"

void create_3d_bc_adh(SMODEL *mod, SGRID *geo3d, int *node_count, int *node_start, int mesh_type, int *node_flags2d, int *node_flags3d, int *new_node_number2d, int *new_node_number3d) {
    
    int i,j,k,kk,nds[2],istring,ie,ie3d,node[NDONQUAD],ordered_node[NDONQUAD],face_flag,imat;
    
    int *face_nodes = (int *) tl_alloc(sizeof(int), geo3d->nnodes);
    for (i=0; i<geo3d->nnodes; i++) face_nodes[i] = 0;
    
    open_output_file(&(mod->io->fout_bc), "3D bc file", FALSE);
    open_output_file(&(mod->io->fout_faces), "3D faces file", FALSE);
    
    fprintf(mod->io->fout_bc.fp, "OP SW3\n");
    if (mod->ntransport > 0) fprintf(mod->io->fout_bc.fp, "OP TRN %d\n",mod->ntransport + 1);
    fprintf(mod->io->fout_bc.fp,"\n");
    fprintf(mod->io->fout_bc.fp, "IP NTL 1e-6\n");
    fprintf(mod->io->fout_bc.fp, "IP ITL 1e-6\n");
    fprintf(mod->io->fout_bc.fp, "IP NIT 10\n");
    fprintf(mod->io->fout_bc.fp, "IP MIT 100\n");
    fprintf(mod->io->fout_bc.fp,"\n");
    fprintf(mod->io->fout_bc.fp, "MP MUC 1.0\n");
    fprintf(mod->io->fout_bc.fp, "MP MU 0.0000001\n");
    fprintf(mod->io->fout_bc.fp, "MP RHO 1000\n");
    fprintf(mod->io->fout_bc.fp, "MP G 9.8\n");
    fprintf(mod->io->fout_bc.fp,"\n");
    if (mod->ntransport > 0) {
        for (j=0; j<mod->ntransport; j++) {
            if (mod->con[j].type == CON) fprintf(mod->io->fout_bc.fp, "CN CON %d 0.0\n",j+1);
            if (mod->con[j].type == SAL) fprintf(mod->io->fout_bc.fp, "CN SAL %d 0.0\n",j+1);
            if (mod->con[j].type == TMP) fprintf(mod->io->fout_bc.fp, "CN TMP %d 0.0\n",j+1);
        }
    }
    
    for (i=0; i<mod->nmat; i++) {
        fprintf(mod->io->fout_bc.fp, "MP ML %d 0\n",i+1);
        fprintf(mod->io->fout_bc.fp, "MP SRT %d 100\n\n",i+1);
        if (mod->ntransport > 0) {
            for (j=0; j<mod->ntransport; j++) {
                fprintf(mod->io->fout_bc.fp, "MP DF  %d  %d  0.0\n",i+1,j+1);
                fprintf(mod->io->fout_bc.fp, "MP TRT %d  %d 100\n",i+1,j+1);
            }
        }
    }
    fprintf(mod->io->fout_bc.fp,"\n");
    fprintf(mod->io->fout_bc.fp,"# EVS card :: {xx, yy, xy, zz, xz, yz} \n");
    for (i=0; i<mod->nmat; i++) {
        fprintf(mod->io->fout_bc.fp, "MP EVS %d 0.0 0.0 0.0 0.0 0.0 0.0\n",i+1);
    }
    fprintf(mod->io->fout_bc.fp,"\n");
    fprintf(mod->io->fout_bc.fp,"TC T0 %10.5f 0\n",mod->t_init);
    fprintf(mod->io->fout_bc.fp,"TC TF %10.5f 0\n",mod->t_final);
    fprintf(mod->io->fout_bc.fp,"\n");
    
    fprintf(mod->io->fout_bc.fp,"#***********************************************\n");
    fprintf(mod->io->fout_bc.fp,"# SERIES ***************************************\n");
    fprintf(mod->io->fout_bc.fp,"SERIES DT  1  2  0  0  0  0\n");
    fprintf(mod->io->fout_bc.fp,"%10.5f  300.0\n",mod->t_init);
    fprintf(mod->io->fout_bc.fp,"%10.5f  300.0\n",mod->t_final);
    fprintf(mod->io->fout_bc.fp,"\n");
    fprintf(mod->io->fout_bc.fp,"SERIES AWRITE  2  1  0  0  0  0\n");
    fprintf(mod->io->fout_bc.fp,"%10.5f %10.5f 1800 0\n",mod->t_init,mod->t_final);
    fprintf(mod->io->fout_bc.fp,"\n");
    fprintf(mod->io->fout_bc.fp,"SERIES BC 3 2 0 0 0 0\n");
    fprintf(mod->io->fout_bc.fp,"%10.5f 0\n",mod->t_init);
    fprintf(mod->io->fout_bc.fp,"%10.5f 0\n",mod->t_final);
    
    fprintf(mod->io->fout_bc.fp,"#***********************************************\n");
    fprintf(mod->io->fout_bc.fp,"# SIDEWALL BOUNDARY CONDITIONS *****************\n");
    k=0;
    for (i=0; i<mod->nstring; i++) {
        if (mod->str_values[i].ol_flow.bc_flag == BCT_VEL_NEU) {
            k++;
            fprintf(mod->io->fout_bc.fp,"NB VEL %d 3\n",k);
        }
        if (mod->str_values[i].ol_flow.bc_flag == BCT_OUTFLOW) {
            k++;
            fprintf(mod->io->fout_bc.fp,"NB OF %d\n",k);
        }
        if (mod->str_values[i].ol_flow.bc_flag == BCT_DIS_NEU) {
            k++;
            fprintf(mod->io->fout_bc.fp,"NB DIS %d 3\n",k);
        }
        if (mod->str_values[i].ol_flow.bc_flag == BCT_PRS_NEU) {
            k++;
            fprintf(mod->io->fout_bc.fp,"NB OTW %d 3\n",k);
        }
    }
    
    fprintf(mod->io->fout_bc.fp,"#***********************************************\n");
    fprintf(mod->io->fout_bc.fp,"# BED AND SURFACE BOUNDARY CONDITIONS **********\n");
    for (i=0; i<mod->nmat; i++) {
        k++;
        fprintf(mod->io->fout_bc.fp,"NB BED %d 3\n",k);
    }
    for (i=0; i<mod->nmat; i++) {
        k++;
        fprintf(mod->io->fout_bc.fp,"NB FRS %d 3\n",k);
    }
    
    // CJT :: if NS, add node strings for surface pressure condition
    fprintf(mod->io->fout_bc.fp,"\n");
    k++;
    for (i=0; i<mod->grid->nnodes; i++) {
        fprintf(mod->io->fout_bc.fp,"NDS %d %d\n",node_start[i]+1,k);
    }
    fprintf(mod->io->fout_bc.fp,"\n");
    
    fprintf(mod->io->fout_bc.fp, "END\n");
    
    
    //-------------------------------------------------------------------------//
    //-------------------------------------------------------------------------//
    //-------------------------------------------------------------------------//
    face_flag = 2;
    fprintf(mod->io->fout_faces.fp, "! LISTING OF 3D FACES PRODUCTED BY EXTRUSION CODE \n");
    fprintf(mod->io->fout_faces.fp, "! FORMAT :: 3D ELEMENT ID, NODE1, NODE2, ..., BOUNDARY_FLAG, MATERIAL\n\n");
    fprintf(mod->io->fout_faces.fp,"\n!***********************************************\n");
    fprintf(mod->io->fout_faces.fp,  "!************ WALL FACE LISTINGS ***************\n");
    for (istring=0; istring<mod->nstring; istring++) {
        for (ie=0; ie<mod->grid->nelems1d; ie++) {
            if (mod->grid->elem1d[ie].string == istring) {
                nds[0] = mod->grid->elem1d[ie].nodes[0];
                nds[1] = mod->grid->elem1d[ie].nodes[1];
                
                // set the flag for nodes that lie beneath the given nodes
                for (i=0; i<geo3d->nnodes; i++) face_nodes[i] = 0;
                for (i=0; i<2; i++) {
                    for (j=node_start[nds[i]]; j<node_start[nds[i]] + node_count[nds[i]]; j++) {
                        face_nodes[j] = 1;
                        //printf("face node under this edge: %d\n",j+1);
                    }
                }
                
                for (ie3d=0; ie3d<geo3d->nelems3d; ie3d++) {
                    assert(geo3d->elem3d[ie3d].nnodes == 4 || geo3d->elem3d[ie3d].nnodes == 6);
                    /* ncnt is the number of nodes of the element that lie beneath the two surface edge node.
                     * So, if ncnt = 3, we have a face corresponding to the surface edge */
                    int ncnt = 0;
                    for (i=0; i<geo3d->elem3d[ie3d].nnodes; i++) {
                        ncnt += face_nodes[geo3d->elem3d[ie3d].nodes[i]];
                    }
                    //printf("# of nodes that lie beneath two surface edge nodes: %d\n",ncnt);
                    
                    if (ncnt == 3 || ncnt == 4) {
                        if (geo3d->elem3d[ie3d].nnodes == NDONTET)   assert(ncnt==3); // if 3d element is tet, side wall should be triangle
                        if (geo3d->elem3d[ie3d].nnodes == NDONPRISM) assert(ncnt==4); // if 3d element is prism, side wall should be quadrilateral
                        
                        // find 2d face nodes
                        kk=0;
                        for (i=0; i<geo3d->elem3d[ie3d].nnodes; i++) {
                            if (face_nodes[geo3d->elem3d[ie3d].nodes[i]] == 1) {
                                node[kk] = geo3d->elem3d[ie3d].nodes[i];
                                kk++;
                            }
                        }
                        
                        // make sure numbering is correct
                        SVECT triCenter = get_triangle_centroid(geo3d->node[geo3d->elem3d[ie3d].nodes[0]],geo3d->node[geo3d->elem3d[ie3d].nodes[1]],geo3d->node[geo3d->elem3d[ie3d].nodes[2]]);
                        SVECT ref_vector;
                        ref_vector.x = triCenter.x - geo3d->node[node[0]].x;
                        ref_vector.y = triCenter.y - geo3d->node[node[0]].y;
                        ref_vector.z = triCenter.z - geo3d->node[node[0]].z;
                        
                        
                        if (geo3d->elem3d[ie3d].nnodes == NDONPRISM) {  // get numbering right for quads
                            ordered_node[0] = node[1];
                            ordered_node[1] = node[0];
                            ordered_node[2] = node[2];
                            ordered_node[3] = node[3];
                            double check = get_triangle_orientation(geo3d->node[ordered_node[0]], geo3d->node[ordered_node[1]], geo3d->node[ordered_node[3]], ref_vector);
                            if (check > 0) {
                                ordered_node[0] = node[0];
                                ordered_node[1] = node[1];
                                ordered_node[2] = node[3];
                                ordered_node[3] = node[2];
                                check = get_triangle_orientation(geo3d->node[ordered_node[0]], geo3d->node[ordered_node[1]], geo3d->node[ordered_node[3]], ref_vector);
                                if (check > 0) {
                                    printf("ERROR in QUADRILATERAL FACE node arrangement!!!!\n");
                                    exit(-1);
                                }
                            }
                            fprintf(mod->io->fout_faces.fp, "QUAD %7d %7d %7d %7d %7d %7d %7d ## EGS %d %d %d\n", ie3d+1, ordered_node[0]+1, ordered_node[1]+1, ordered_node[2]+1, ordered_node[3]+1,
                                    face_flag, mod->grid->elem1d[ie].string, nds[0]+1, nds[1]+1, istring+1);
                        } else {
                            ordered_node[0] = node[0];
                            ordered_node[1] = node[1];
                            ordered_node[2] = node[2];
                            double check = get_triangle_orientation(geo3d->node[ordered_node[0]], geo3d->node[ordered_node[1]], geo3d->node[ordered_node[2]], ref_vector);
                            if (check > 0) {
                                ordered_node[0] = node[0];
                                ordered_node[1] = node[2];
                                ordered_node[2] = node[1];
                                check = get_triangle_orientation(geo3d->node[ordered_node[0]], geo3d->node[ordered_node[1]], geo3d->node[ordered_node[2]], ref_vector);
                                if (check > 0) {
                                    printf("ERROR in TRIANGULAR FACE node arrangement!!!!\n");
                                    exit(-1);
                                }
                            }
                            
                            // check a few things
                            SVECT nodes[NDONTRI];
                            nodes[0].x = geo3d->node[ordered_node[0]].x; nodes[0].y = geo3d->node[ordered_node[0]].y; nodes[0].z = geo3d->node[ordered_node[0]].z;
                            nodes[1].x = geo3d->node[ordered_node[1]].x; nodes[1].y = geo3d->node[ordered_node[1]].y; nodes[1].z = geo3d->node[ordered_node[1]].z;
                            nodes[2].x = geo3d->node[ordered_node[2]].x; nodes[2].y = geo3d->node[ordered_node[2]].y; nodes[2].z = geo3d->node[ordered_node[2]].z;
                            SVECT nrml = get_elem2d_normals(nodes);
                            if (fabs(nrml.z) > SMALL && face_flag == 2) {
                                printf("sidewall nrml.z: %20.10f \n",nrml.z);
                                 tl_error("ERROR :: Sidewall z-normal should be 0");
                            }
                            
                            
                            fprintf(mod->io->fout_faces.fp, "TRI %7d %7d %7d %7d %7d %7d ## EGS %d %d %d\n", ie3d+1, ordered_node[0]+1, ordered_node[1]+1, ordered_node[2]+1,
                                    face_flag, mod->grid->elem1d[ie].string, nds[0]+1, nds[1]+1, istring+1);

                        }
                    }
                }
            }
        }
    }
    int string_count = mod->nstring;
    
    // cjt :: we must find the 3D element each boundary face belongs too ... there is probably a more efficient way to do this,
    // but I'm busy ... 
    int count_bed = 0, count_sur = 0, sur_node = UNSET_INT, bed_node = UNSET_INT;
    int map_2d_to_3d_sur[mod->grid->nelems2d];
    int map_2d_to_3d_bed[mod->grid->nelems2d];
    for (ie=0; ie<mod->grid->nelems2d; ie++) {
        map_2d_to_3d_sur[ie] = UNSET_INT;
        for (ie3d=0; ie3d<geo3d->nelems3d; ie3d++) {
            count_sur = 0;
            for (i=0; i<mod->grid->elem2d[ie].nnodes; i++) {
                for (j=0; j<geo3d->elem3d[ie3d].nnodes; j++) {
                    if (node_start[mod->grid->elem2d[ie].nodes[i]] == geo3d->elem3d[ie3d].nodes[j] ) {
                        count_sur++;
                        break;
                    }
                }
            }
            if (count_sur == mod->grid->elem2d[ie].nnodes) {
                map_2d_to_3d_sur[ie] = ie3d;
                break;
            }
        }
    }
    for (ie=0; ie<mod->grid->nelems2d; ie++) {
        map_2d_to_3d_bed[ie] = UNSET_INT;
        for (ie3d=0; ie3d<geo3d->nelems3d; ie3d++) {
            count_bed = 0;
            for (i=0; i<mod->grid->elem2d[ie].nnodes; i++) {
                bed_node = node_start[mod->grid->elem2d[ie].nodes[i]] + node_count[mod->grid->elem2d[ie].nodes[i]] - 1;
                for (j=0; j<geo3d->elem3d[ie3d].nnodes; j++) {
                    if (bed_node == geo3d->elem3d[ie3d].nodes[j] ) {
                        count_bed++;
                        break;
                    }
                }
            }
            if (count_bed == mod->grid->elem2d[ie].nnodes) {
                map_2d_to_3d_bed[ie] = ie3d;
                break;
            }

        }
    }
    //exit(-1);

        
    // FIX LATER :: this only works for one string for all surface/bed materials.
    // TRIANGULAR SURFACE FACES :: relies on 2d grid proper element node numbering
    face_flag = 0;
    fprintf(mod->io->fout_faces.fp, "\n!***********************************************\n");
    fprintf(mod->io->fout_faces.fp,   "!**************** SURFACE FACES ****************\n");
    for (imat=0; imat<mod->nmat; imat++) {
        fprintf(mod->io->fout_faces.fp, "!------------------MATERIAL # %d ----------------\n",imat+1);
        for (ie=0; ie<mod->grid->nelems2d; ie++) {
            if (mod->grid->elem2d[ie].mat == imat) {
                //fprintf(mod->io->fout_faces.fp, "TRI %7d %7d %7d %7d %7d %7d\n", ie + 1,
                fprintf(mod->io->fout_faces.fp, "TRI %7d %7d %7d %7d %7d %7d\n", map_2d_to_3d_sur[ie] + 1,
                        node_start[ mod->grid->elem2d[ie].nodes[0] ] + 1,
                        node_start[ mod->grid->elem2d[ie].nodes[1] ] + 1,
                        node_start[ mod->grid->elem2d[ie].nodes[2] ] + 1, face_flag, string_count);
            }
        }
    }
    string_count++;
    
    
    // TRIANGULAR BED FACES :: relies on 2d grid proper element node numbering :: CJT :: reverse so nz is -!
    face_flag = 1;
    fprintf(mod->io->fout_faces.fp, "\n!***********************************************\n");
    fprintf(mod->io->fout_faces.fp,   "!****************** BED FACES ******************\n");
    for (imat=0; imat<mod->nmat; imat++) {
        fprintf(mod->io->fout_faces.fp, "!------------------MATERIAL # %d ----------------\n",imat+1);
        for (ie=0; ie<mod->grid->nelems2d; ie++) {
                            ordered_node[0] = mod->grid->elem2d[ie].nodes[0];
                            ordered_node[1] = mod->grid->elem2d[ie].nodes[2];
                            ordered_node[2] = mod->grid->elem2d[ie].nodes[1];
            if (mod->grid->elem2d[ie].mat == imat) {
                fprintf(mod->io->fout_faces.fp, "TRI %7d %7d %7d %7d %7d %7d\n", map_2d_to_3d_bed[ie] + 1,
                        node_start[ordered_node[0]] + node_count[ordered_node[0]],
                        node_start[ordered_node[1]] + node_count[ordered_node[1]],
                        node_start[ordered_node[2]] + node_count[ordered_node[2]], 
						face_flag, string_count);
            }
        }
    }     
    
    fclose(mod->io->fout_faces.fp);
    fclose(mod->io->fout_bc.fp);
}

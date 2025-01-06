#include "adh.h"
//can expose to frontend later, but hard coding strategy for now
static char reorder_strat[] = "g{pass=10}";
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Renumbers nodes to minimize bandwidth using SCOTCH library, by default
 * 			   uses Gibbs-Poole algorithm.
 *  		   Should only be called for serial runs. Upon completion, node ids will
 * 	           be permuted and an inverse permutation array will be filled out 
 * 			   for convenient output.
 * 			   Designed for only serial and initial mesh, will not update be valid if grid is
 *             partitioned or adapted.
 *  \author    Corey Trahan, Ph.D.
 *  \author    Mark Loveand, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[in,out] pgrid (SGRID *)  pointer to an AdH grid
 * @returns int - 0 if succesful
 *
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

int sgrid_reorder(SGRID *grid){
#ifdef _MPI
	return 0;
#else
	//maybe debug only
	assert(sizeof(SCOTCH_Num) == sizeof(int));
//	if (sizeof (SCOTCH_Num) > sizeof (int)) {
//    	printf("PT-Scotch users, beware: here be dragons!!!\n\n\n");
//    }else{
//    	printf("SCOTCH NUM size = int \n");
//    }
	//serial runs
	int i,j,l,err;
	int nelems3d= grid->nelems3d;
	int nelems2d= grid->nelems2d;
	int nelems1d= grid->nelems1d;
	int nelems = nelems3d + nelems2d + nelems1d;
	SCOTCH_Num k;
	//translate to SGRID to SCOTCH_Mesh object
	//node counts always start at base 0
	SCOTCH_Num velmbas = 0;
	//always start node numbering at nelem
	SCOTCH_Num vnodbas = nelems;
	//number of elements in mesh
	SCOTCH_Num velmnbr = nelems;
	//number of nodes in mesh
	SCOTCH_Num vnodnbr = grid->nnodes;
	//now need number of arcs in mesh
	//this by construction will be 2*sum(nnodes_on_elem) for each element
	//since each element on this mesh has same number of nodes, it is easy to find
	SCOTCH_Num edgenbr;
	//create vertex table, this will be the node connectivity appended by node to node connectivity
	//unfortunately must store in new array
	SCOTCH_Num *verttab;
	//verttab = (SCOTCH_Num *) malloc(sizeof(SCOTCH_Num)*(vnodnbr+velmnbr+1));
	verttab = (SCOTCH_Num *) tl_alloc(sizeof(SCOTCH_Num),vnodnbr+velmnbr+1);
	//short routine to assign appropriate values, the first velemnbr will be easy
	j=0;k=velmbas;
	for(i=0;i<nelems3d;i++){
		verttab[j] = k; 
		k += grid->elem3d[i].nnodes;
		j+=1; 
	}
	for(i=0;i<nelems2d;i++){ 
		verttab[j] = k;
		k += grid->elem2d[i].nnodes;
		j+=1;
	}
	for(i=0;i<nelems1d;i++){
		verttab[j] = k; 
		k += grid->elem1d[i].nnodes;
		j+=1;
 	}
 	verttab[j] = k;
 	edgenbr = 2*k;
 	//printf("EDGENBR = %d\n",edgenbr);

 	//the nodes will not be as trivial, need to know how many elements each node is connected to
	//loop through elements and add up
	//need nodal array to store this
	int *temp_nodal_connections;
	//temp_nodal_connections = (int *) malloc(sizeof(int)*vnodnbr);
	(temp_nodal_connections) = (int *) tl_alloc(sizeof(int), vnodnbr);
	//int temp_nodal_connections[nnodes];
	//initialize to 0
	sarray_init_int(temp_nodal_connections,vnodnbr);
	//simple way to get the nodal connections
	//loop through each element and add
	int nnodes_on_elem,nd1;
	j=0;
	for(i=0;i<nelems3d;i++){
		nnodes_on_elem = grid->elem3d[i].nnodes;
		for(l=0;l<nnodes_on_elem;l++){
			nd1 = grid->elem3d[i].nodes[l];
			temp_nodal_connections[nd1]+=1;
		}
	}
	for(i=0;i<nelems2d;i++){
		nnodes_on_elem = grid->elem2d[i].nnodes;
		for(l=0;l<nnodes_on_elem;l++){
			nd1 = grid->elem2d[i].nodes[l];
			temp_nodal_connections[nd1]+=1;
		}
	}
	for(i=0;i<nelems1d;i++){
		nnodes_on_elem = grid->elem1d[i].nnodes;
		for(l=0;l<nnodes_on_elem;l++){
			nd1 = grid->elem1d[i].nodes[l];
			temp_nodal_connections[nd1]+=1;
		}
	}
	//dump this nodal connections into verttab
	int idx = verttab[velmnbr];
	for(i=0;i<vnodnbr;i++){
		idx += temp_nodal_connections[i];
		verttab[velmnbr+i+1] = idx;
	}
	//destroy temp_nodal_connections
	temp_nodal_connections = (int *) tl_free(sizeof(int), vnodnbr, temp_nodal_connections);
	//now make the edge table
	//dynamically allocate
	SCOTCH_Num *edgetab;
	//edgetab = (SCOTCH_Num *) malloc(sizeof(SCOTCH_Num)*(edgenbr));
	(edgetab) = (SCOTCH_Num *) tl_alloc(sizeof(SCOTCH_Num), edgenbr);
	//printf("Vertex edge table allocated\n");
	//fill in edge table
	//first part is just element to node connections
	//second part is node to element connections
	idx = 0;
	for(j=0;j<nelems3d;j++){
		nnodes_on_elem = grid->elem3d[j].nnodes;
		for(i=0;i<nnodes_on_elem;i++){
			//offset by nelems
			edgetab[idx] = grid->elem3d[j].nodes[i]+nelems;
			idx+=1;
		}
	}
	for(j=0;j<nelems2d;j++){
		nnodes_on_elem = grid->elem2d[j].nnodes;
		for(i=0;i<nnodes_on_elem;i++){
			//offset by nelems
			edgetab[idx] = grid->elem2d[j].nodes[i]+nelems;
			idx+=1;
		}
	}
	for(j=0;j<nelems1d;j++){
		nnodes_on_elem = grid->elem1d[j].nnodes;
		for(i=0;i<nnodes_on_elem;i++){
			//offset by nelems
			edgetab[idx] = grid->elem1d[j].nodes[i]+nelems;
			idx+=1;
		}
	}
	//second part is more tricky, not sure what is quickest
	//initialize to 0
	int *nodal_offsets;
	//nodal_offsets = (int*) malloc(sizeof(int)*vnodnbr);
	(nodal_offsets) = (int *) tl_alloc(sizeof(int), vnodnbr);
	sarray_init_int(nodal_offsets,vnodnbr);
	int offset, temp_val;
	for(j=0;j<nelems3d;j++){
		nnodes_on_elem = grid->elem3d[j].nnodes;
		for(i=0;i<nnodes_on_elem;i++){
			//temp_val=0;
			nd1 = grid->elem3d[j].nodes[i];
			//use nodal offsets to get proper index
			offset = nodal_offsets[nd1];
			//use verttab to get proper temp_val
			temp_val = verttab[nelems+nd1];
			idx = temp_val+offset;
			edgetab[idx] = j;
			nodal_offsets[nd1]+=1;
		}
	}
	for(j=0;j<nelems2d;j++){
		nnodes_on_elem = grid->elem2d[j].nnodes;
		for(i=0;i<nnodes_on_elem;i++){
			//temp_val=0;
			nd1 = grid->elem2d[j].nodes[i];
			//use nodal offsets to get proper index
			offset = nodal_offsets[nd1];
			//use verttab to get proper temp_val
			temp_val = verttab[nelems+nd1];
			idx = temp_val+offset;
			edgetab[idx] = j;
			nodal_offsets[nd1]+=1;
		}
	}
	for(j=0;j<nelems1d;j++){
		nnodes_on_elem = grid->elem1d[j].nnodes;
		for(i=0;i<nnodes_on_elem;i++){
			//temp_val=0;
			nd1 = grid->elem1d[j].nodes[i];
			//use nodal offsets to get proper index
			offset = nodal_offsets[nd1];
			//use verttab to get proper temp_val
			temp_val = verttab[nelems+nd1];
			idx = temp_val+offset;
			edgetab[idx] = j;
			nodal_offsets[nd1]+=1;
		}
	}
	//clear nodal_offsets
	nodal_offsets = (int *) tl_free(sizeof(int), vnodnbr, nodal_offsets);
	//create pointer for vendtab, not necessary NULL can be passed
	//SCOTCH_Num *vendtab;
    //vendtab = verttab + 1;
    SCOTCH_Mesh meshdat;
    err = SCOTCH_meshInit(&meshdat);
    SCOTCH_Strat  stratdat;
    err = SCOTCH_stratInit (&stratdat);
    //either routine seems to work? Idk what is good
    //g is Gibbs-Poole (default)
    //d Block Halo Approximate Minimum Degree method (not recommended)
    //f Block Halo Approximate Minimum Fill method (not recommended)
    err = SCOTCH_stratMeshOrder(&stratdat,reorder_strat); 
    //this all seems to not work great, breaks above 10x10 mesh ?
    //using defaults in source code
    //SCOTCH_Num mesh_strat = SCOTCH_STRATQUALITY;
    //other options
    //also doesnt work??
    //mesh_strat = SCOTCH_STRATSPEED;
    //default breaks?
    //mesh_strat = SCOTCH_STRATDEFAULT;
    //also doesnt work?
    //mesh_strat = SCOTCH_STRATBALANCE;
    //err = SCOTCH_stratMeshOrderBuild(&stratdat,SCOTCH_STRATQUALITY,0.1);
    err = SCOTCH_meshBuild (&meshdat, velmbas, vnodbas, velmnbr, vnodnbr, verttab, NULL,
	NULL, NULL, NULL, edgenbr, edgetab);
#ifdef _DEBUG
	err =SCOTCH_meshCheck(&meshdat);
	printf("Completed SCOTCH MESH Check, err = %d\n",err);
#endif
	//Permutation tables is what we want, will embed these into SGRID
	SCOTCH_Num *permtab;
	SCOTCH_Num *peritab;
	//should it be nnodes or nvertices of graph? seems to be nnodes
	int nvertices = vnodnbr+velmnbr+1;
	//permtab = (SCOTCH_Num*) malloc(sizeof(SCOTCH_Num)*(nnodes));
	//peritab = (SCOTCH_Num*) malloc(sizeof(SCOTCH_Num)*(nnodes));
	(permtab) = (SCOTCH_Num *) tl_alloc(sizeof(SCOTCH_Num), vnodnbr);
	(peritab) = (SCOTCH_Num *) tl_alloc(sizeof(SCOTCH_Num), vnodnbr);
	//intialize to 0
	sarray_init_int(permtab,vnodnbr);
	sarray_init_int(peritab,vnodnbr);
	//printf("Calling SCOTCH MESH reorder\n");
	//seg faults any thing above 10x10 mesh??
 	err = SCOTCH_meshOrder(&meshdat,&stratdat, permtab, peritab, NULL, NULL, NULL);
 	//printf("SCOTCH MESH reorder completed\n");
 	//printf("Permutation info generated\n");
 	//for(k=0;k<vnodnbr;k++){printf("permtab[%d] = %d, peritab[%d] = %d\n",k,permtab[k],k,peritab[k]);}

 	//use permtab to reorder nodes and then update connectivity
 	//first loop through elements and update id's
 	for(j=0;j<nelems3d;j++){
		nnodes_on_elem = grid->elem3d[j].nnodes;
		for(i=0;i<nnodes_on_elem;i++){
			nd1 = grid->elem3d[j].nodes[i];
			grid->elem3d[j].nodes[i] = permtab[nd1];
		}
	}
	for(j=0;j<nelems2d;j++){
		nnodes_on_elem = grid->elem2d[j].nnodes;
		for(i=0;i<nnodes_on_elem;i++){
			nd1 = grid->elem2d[j].nodes[i];
			grid->elem2d[j].nodes[i] = permtab[nd1];
		}
	}
	for(j=0;j<nelems1d;j++){
		nnodes_on_elem = grid->elem1d[j].nnodes;
		for(i=0;i<nnodes_on_elem;i++){
			nd1 = grid->elem1d[j].nodes[i];
			grid->elem1d[j].nodes[i] = permtab[nd1];
		}
	}
	//store permtab in node.id
	for(i=0;i<vnodnbr;i++){
		grid->node[i].id = permtab[i];
	}
	//store peritab in sgrid
	grid->inv_per_node = peritab;


	//printf("Rank %d out of %d\n",rank,npe);
	//err = SCOTCH_dgraphInit(&grafdat, MPI_COMM_WORLD);

 	//free up memory
 	SCOTCH_meshExit(&meshdat);
 	SCOTCH_stratExit(&stratdat);
 	verttab = (SCOTCH_Num *) tl_free(sizeof(SCOTCH_Num),vnodnbr+velmnbr+1,verttab);
 	edgetab = (SCOTCH_Num *) tl_free(sizeof(SCOTCH_Num), edgenbr,edgetab);
 	permtab = (SCOTCH_Num *) tl_free(sizeof(SCOTCH_Num), vnodnbr, permtab);

	if ( err != 0) {
		printf("ERROR in dgraph INIT!\n");
	}else{
		printf("SCOTCH reorder completed\n");
	}
	return err;
#endif


	

}
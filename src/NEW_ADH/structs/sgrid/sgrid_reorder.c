#include "adh.h"
//can expose to frontend later, but hard coding strategy for now
static char reorder_strat[] = "g{pass=10}";
//static char reorder_strat[] = "c{rat=0.7,cpr=g{pass=10},unc=g{pass=10}}";
static SCOTCH_Mesh sgrid_to_scotch_mesh(SGRID *grid);
static SCOTCH_Graph sgrid_to_scotch_graph(SGRID *grid);
static int reverse_cuthill_mckee(SGRID *grid, int *permtab, int *peritab);
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
int sgrid_reorder(SGRID *grid, int option){
#ifdef _MPI
	return 0;
#else
	int seed=0;
	int err=0;
	int nnodes = grid->nnodes;
	int nelems1d = grid->nelems1d;
	int nelems2d = grid->nelems2d;
	int nelems3d = grid->nelems3d;
	int nd1;
	int i,j,k;
	int nnodes_on_elem;
	//maybe debug only
#ifdef _DEBUG
	assert(sizeof(SCOTCH_Num) == sizeof(int));
	assert (SCOTCH_numSizeof () == sizeof (SCOTCH_Num));
    assert (option == 1 | option == 2 | option == 3);
#endif
    //Permutation tables is what we want, will embed these into SGRID
    SCOTCH_Num *permtab;
    SCOTCH_Num *peritab;
    SCOTCH_Graph grafdat;
    SCOTCH_Strat  stratdat;
    SCOTCH_Mesh meshdat;
    permtab = (SCOTCH_Num *) tl_alloc(sizeof(SCOTCH_Num), nnodes);
    peritab = (SCOTCH_Num *) tl_alloc(sizeof(SCOTCH_Num), nnodes);
	//OPTION 1
    if(option == 1){
        //routine that builds SCOTCHmesh then sends to graph
        err = SCOTCH_graphInit(&grafdat);
        err = SCOTCH_meshInit(&meshdat);
        meshdat = sgrid_to_scotch_mesh(grid);
        err = SCOTCH_meshGraph(&meshdat, &grafdat);
    }
    else if(option == 2){
	    //OPTION 2
        //routine that builds graph directly from mesh
        grafdat = sgrid_to_scotch_graph(grid);

    }
    else if (option ==3){
        //in house cuthill-mckee
        err = reverse_cuthill_mckee(grid, permtab,peritab);
    }
    else{
        //if none of these options, return an error
        return -1;
    }


    //intermediate step if using SCOTCH
    if ( (option == 1) | (option == 2)){
        err = SCOTCH_stratInit(&stratdat);
#ifdef _DEBUG
    assert(err==0);
#endif
        err = SCOTCH_stratGraphOrder(&stratdat, reorder_strat);
#ifdef _DEBUG
    assert(err==0);
#endif
	   //intialize to 0
	   sarray_init_int(permtab,nnodes);
	   sarray_init_int(peritab,nnodes);
	   //a way to do same things more granularly
 	  //SCOTCH_Ordering     ordedat;
 	  //SCOTCH_graphOrderInit (&grafdat, &ordedat, NULL, NULL, NULL, NULL, NULL);
 	  //SCOTCH_graphOrderCompute (&grafdat, &ordedat, &stratdat);
 	  //SCOTCH_graphOrderCheck (&grafdat, &ordedat);
 	  //can write to file
 	  //FILE *fileptr = fopen("reordering.txt", "w");
 	  //SCOTCH_graphOrderSave     (&grafdat, &ordedat, fileptr); 
 	  //a convenience routine that wraps OrderInit and OrderCompute
 	  err = SCOTCH_graphOrder(&grafdat, &stratdat ,permtab, peritab, NULL, NULL, NULL);
    }
#ifdef _DEBUG
    assert(err==0);
    //for(k=0;k<nnodes;k++){
 	//	printf("permtab[%d] = %d, peritab[%d] = %d\n",k,permtab[k],k,peritab[k]);}
#endif

    //FOR ALL OPTIONS	
 	//use permtab to reorder nodes and then update connectivity
 	//first loop through elements and update id's
 	//ACTUALLY WILL BREAK STUFF, WE DONT CHANGE CONNECTIONS BUT NODE ORDERINGS
 	//JUST ADD PERMUTATION IN GRID INSTEAD
 	//OR DO THIS AND SWITCH ACTUAL NODE ORDER, SOMETHING IS CAUSING
 	//PRECISION ISSUES IN THIS BLOCK
 	//what if we  do nothing
 	for(j=0;j<nelems3d;j++){
		nnodes_on_elem = grid->elem3d[j].nnodes;
		for(i=0;i<nnodes_on_elem;i++){
			nd1 = grid->elem3d[j].nodes[i];
			grid->elem3d[j].nodes[i] = permtab[nd1];
		}
	}
	for(j=0;j<nelems2d;j++){
		nnodes_on_elem = grid->elem2d[j].nnodes;
		//printf("Elem [%d] Old nodes: %d, %d, %d\n",j,grid->elem2d[j].nodes[0],grid->elem2d[j].nodes[1],grid->elem2d[j].nodes[2]);
		for(i=0;i<nnodes_on_elem;i++){
			nd1 = grid->elem2d[j].nodes[i];
			grid->elem2d[j].nodes[i] = permtab[nd1];
		}
		//printf("Elem [%d] New nodes: %d, %d, %d\n",j,grid->elem2d[j].nodes[0],grid->elem2d[j].nodes[1],grid->elem2d[j].nodes[2]);
	}
	for(j=0;j<nelems1d;j++){
		nnodes_on_elem = grid->elem1d[j].nnodes;
		for(i=0;i<nnodes_on_elem;i++){
			nd1 = grid->elem1d[j].nodes[i];
			grid->elem1d[j].nodes[i] = permtab[nd1];
		}
	}
	//Reorder node data, must copy data temporarily
	SNODE *node_temp;
	snode_alloc_array(&node_temp,nnodes) ;
	snode_copy_array(node_temp, grid->node, nnodes);
	for(i=0;i<nnodes;i++){
		//overwrite data
		//node permtab[i] is new location of original node i
		//grid->node[permtab[i]] = node_temp[i];
		snode_copy(&(grid->node[permtab[i]]), node_temp[i]);
	}
	//store permtab in node.id? or just separately?
	//not used because ids are stored in elem.nodes
	for(i=0;i<nnodes;i++){
		//grid->node[i].id = i;
		grid->inv_per_node[i] = peritab[i];
	}
	//free node_temp
	snode_free_array(node_temp, nnodes);
 	//free up memory if SCOTCH used
    if (option == 1){
 	  SCOTCH_meshExit(&meshdat);
    }
    if (option == 1 || option ==2){
 	  SCOTCH_stratExit(&stratdat);
 	  SCOTCH_graphExit(&grafdat);
    }
 	permtab = (SCOTCH_Num *) tl_free(sizeof(SCOTCH_Num), nnodes, permtab);
 	peritab = (SCOTCH_Num *) tl_free(sizeof(SCOTCH_Num), nnodes, peritab);
	if ( err != 0) {
		printf("ERROR in dgraph INIT!\n");
	}else{
		printf("SCOTCH reorder completed\n");
	}
	return err;
#endif

}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Renumbers nodes to minimize bandwidth using SCOTCH library, by default
 *             uses Gibbs-Poole algorithm.
 *             Should only be called for serial runs. Upon completion, node ids will
 *             be permuted and an inverse permutation array will be filled out 
 *             for convenient output.
 *             Designed for only serial and initial mesh, will not update be valid if grid is
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
//helper function that creates nodal graph (not the bipartite one)
//just node to node connections
SCOTCH_Graph sgrid_to_scotch_graph(SGRID *grid){
	int i,j,k;
	int err;
	int n_connections;
	SCOTCH_Graph grafdat;
	err = SCOTCH_graphInit(&grafdat);
	int nelems3d = grid->nelems3d;
	int nelems2d = grid->nelems2d;
	int nelems1d = grid->nelems1d;
	int current_node, other_node;
	int count=0;
	//nodes always start at 0
	SCOTCH_Num baseval=0;
	SCOTCH_Num vertnbr=grid->nnodes;
	int nnodes_on_elem;
	int *nodes; //for aliasing
	int Nedges=0;
	int *n_con;
	int *n_con_no_dup;
	int *edgetab;
	int **temp_edgetab;
	int *verttab;
	//quickly compute number of edges in the mesh
	//the nodes will not be as trivial, need to know how many elements each node is connected to
	//loop through elements and add up
	//need nodal array to store this

	//really one in the same as slin_sys_init_sparsity_mono
	//need two nodal arrays
	//these will store the number of nodes each node is connected to
	n_con = (int *) tl_alloc(sizeof(int), vertnbr);
    n_con_no_dup = (int *) tl_alloc(sizeof(int), vertnbr);
    verttab = (int *) tl_alloc(sizeof(int),vertnbr+1);
    sarray_init_int(n_con, vertnbr);
    sarray_init_int(n_con_no_dup, vertnbr);
    //First set of loops is solely to establish how many nodes are connected to each node
    for (i=0;i<nelems3d;i++){
    	//nnodes on the element
        nnodes_on_elem = grid->elem3d[i].nnodes;
        //Could get stuff for node weights from physics mat
       	//for this element, find nodes
        nodes = grid->elem3d[i].nodes;
        //keep sparsity pattern in temp_cols_diag, temp_cols_off_diag
        for (j=0;j<nnodes_on_elem;j++){
            //each node inside an element, add this to the count
            current_node = nodes[j];
            //loop through each node that is not the node itself
            //nodes before current node
            for (k=0;k<j;k++){
            	//other_node = nodes[k];
                n_con[current_node]+=1;
            }
            for(k=j+1;k<nnodes_on_elem;k++){
            	//other_node = nodes[k];
            	n_con[current_node]+=1;
            }
        }
    }
    for (i=0;i<nelems2d;i++){
    	//nnodes on the element
        nnodes_on_elem = grid->elem2d[i].nnodes;
        //Could get stuff for node weights from physics mat
       	//for this element, find nodes
        nodes = grid->elem2d[i].nodes;
        //keep sparsity pattern in temp_cols_diag, temp_cols_off_diag
        for (j=0;j<nnodes_on_elem;j++){
            //each node inside an element, add this to the count
            current_node = nodes[j];
            //loop through each node that is not the node itself
            //nodes before current node
            for (k=0;k<j;k++){
            	//other_node = nodes[k];
                n_con[current_node]+=1;
            }
            for(k=j+1;k<nnodes_on_elem;k++){
            	//other_node = nodes[k];
            	n_con[current_node]+=1;
            }
        }
    }
    for (i=0;i<nelems1d;i++){
    	//nnodes on the element
        nnodes_on_elem = grid->elem1d[i].nnodes;
        //Could get stuff for node weights from physics mat
       	//for this element, find nodes
        nodes = grid->elem1d[i].nodes;
        //keep sparsity pattern in temp_cols_diag, temp_cols_off_diag
        for (j=0;j<nnodes_on_elem;j++){
            //each node inside an element, add this to the count
            current_node = nodes[j];
            //loop through each node that is not the node itself
            //nodes before current node
            for (k=0;k<j;k++){
            	//other_node = nodes[k];
                n_con[current_node]+=1;
            }
            for(k=j+1;k<nnodes_on_elem;k++){
            	//other_node = nodes[k];
            	n_con[current_node]+=1;
            }
        }
    }
    //use nnz_rows to dynamically allocate
    //int temp_cols_diag[nrows][nCon3d];
    temp_edgetab = (int**) tl_alloc(sizeof(int*), vertnbr);
    for(j=0;j<vertnbr;j++){
        temp_edgetab[j] = (int*) tl_alloc(sizeof(int), n_con[j]);
        for(k=0;k<n_con[j];k++){
            temp_edgetab[j][k]=INT_MAX;
        }
    }
    //Seems redundant but must reuse as indexing
    sarray_init_int(n_con, vertnbr);
    //First set of loops is solely to establish how many nodes are connected to each node
    for (i=0;i<nelems3d;i++){
    	//nnodes on the element
        nnodes_on_elem = grid->elem3d[i].nnodes;
        //Could get stuff for node weights from physics mat
       	//for this element, find nodes
        nodes = grid->elem3d[i].nodes;
        //keep sparsity pattern in temp_cols_diag, temp_cols_off_diag
        for (j=0;j<nnodes_on_elem;j++){
            //each node inside an element, add this to the count
            current_node = nodes[j];
            //loop through each node that is not the node itself
            //nodes before current node
            for (k=0;k<j;k++){
            	other_node = nodes[k];
                temp_edgetab[current_node][n_con[current_node]]=other_node;
                n_con[current_node]+=1;
            }
            for(k=j+1;k<nnodes_on_elem;k++){
            	other_node = nodes[k];
            	temp_edgetab[current_node][n_con[current_node]]=other_node;
                n_con[current_node]+=1;
            }
        }
    }
    for (i=0;i<nelems2d;i++){
    	//nnodes on the element
        nnodes_on_elem = grid->elem2d[i].nnodes;
        //Could get stuff for node weights from physics mat
       	//for this element, find nodes
        nodes = grid->elem2d[i].nodes;
        //keep sparsity pattern in temp_cols_diag, temp_cols_off_diag
        for (j=0;j<nnodes_on_elem;j++){
            //each node inside an element, add this to the count
            current_node = nodes[j];
            //loop through each node that is not the node itself
            //nodes before current node
            for (k=0;k<j;k++){
            	other_node = nodes[k];
                temp_edgetab[current_node][n_con[current_node]]=other_node;
                n_con[current_node]+=1;
            }
            for(k=j+1;k<nnodes_on_elem;k++){
            	other_node = nodes[k];
            	temp_edgetab[current_node][n_con[current_node]]=other_node;
                n_con[current_node]+=1;
            }
        }
    }
    for (i=0;i<nelems1d;i++){
    	//nnodes on the element
        nnodes_on_elem = grid->elem1d[i].nnodes;
        //Could get stuff for node weights from physics mat
       	//for this element, find nodes
        nodes = grid->elem1d[i].nodes;
        //keep sparsity pattern in temp_cols_diag, temp_cols_off_diag
        for (j=0;j<nnodes_on_elem;j++){
            //each node inside an element, add this to the count
            current_node = nodes[j];
            //loop through each node that is not the node itself
            //nodes before current node
            for (k=0;k<j;k++){
            	other_node = nodes[k];
                temp_edgetab[current_node][n_con[current_node]]=other_node;
                n_con[current_node]+=1;
            }
            for(k=j+1;k<nnodes_on_elem;k++){
            	other_node = nodes[k];
            	temp_edgetab[current_node][n_con[current_node]]=other_node;
                n_con[current_node]+=1;
            }
        }
    }
    for (i=0;i<vertnbr;i++){
        // sort the column indices (j-entries)
        //use stdlib.h qsort
        qsort(temp_edgetab[i], n_con[i], sizeof(int), compare_ints);
        //this should hopefully remove duplicates?
        n_connections = sarray_unique_int(temp_edgetab[i], n_con[i]);
        //overwrite nnz row with sarray_unique_int?
        n_con_no_dup[i] = n_connections;
        //add nnz in a row to the NNZ
        //maybe overwrire nnz_rows[i] instead if we want this stored, and then sum it
        Nedges+=n_connections;
    }
    //allocate adjacency data
    edgetab = (int *) tl_alloc(sizeof(int), Nedges);
    int *nodal_edges; //alias

    //now use info to fill in indptr and cols
    for(i=0;i<vertnbr;i++){
        //printf("filling in index ptr and column entries, row %d\n",i);
        verttab[i] = count;
        //printf("filling in index ptr and column entries, row %d\n",i);
        nodal_edges = temp_edgetab[i];
        //think about how to do this since each temp_rows may be different size
        for(j=0;j<n_con_no_dup[i];j++){
            edgetab[count] = nodal_edges[j];
            count++;
        }
    }
    //also last entry
    verttab[vertnbr] = count;

    //build graph
    err = SCOTCH_graphBuild(&grafdat,baseval,vertnbr,verttab,verttab+1,NULL,NULL,Nedges,edgetab,NULL);
#ifdef _DEBUG
    assert(err==0);
#endif
    err = SCOTCH_graphCheck(&grafdat);
#ifdef _DEBUG
    assert(err==0);
#endif
    //freeing memory
    n_con_no_dup = (int *) tl_free(sizeof(int), vertnbr, n_con_no_dup);
    for(j=0;j<vertnbr;j++){
        temp_edgetab[j] = (int*) tl_free(sizeof(int), n_con[j],temp_edgetab[j]);
    }
    temp_edgetab = (int**) tl_free(sizeof(int*), vertnbr,temp_edgetab);
    n_con = (int *) tl_free(sizeof(int), vertnbr, n_con);
    //not freeing here, hopefully this gets freed when graph free is called later
    //verttab = (int *) tl_free(sizeof(int), vertnbr, verttab);
    //edgetab = (int *) tl_free(sizeof(int), Nedges,edgetab);
    return grafdat;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Renumbers nodes to minimize bandwidth using SCOTCH library, by default
 *             uses Gibbs-Poole algorithm.
 *             Should only be called for serial runs. Upon completion, node ids will
 *             be permuted and an inverse permutation array will be filled out 
 *             for convenient output.
 *             Designed for only serial and initial mesh, will not update be valid if grid is
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
 //the bipartite graph which I dont think works as intended
 SCOTCH_Mesh sgrid_to_scotch_mesh(SGRID *grid){
 	int err;
 	SCOTCH_Mesh meshdat;
    err = SCOTCH_meshInit(&meshdat);
 	//serial runs
	int i,j,l;
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
	SCOTCH_Num nverttab = vnodnbr+velmnbr+1;
	verttab = (SCOTCH_Num *) tl_alloc(sizeof(SCOTCH_Num),nverttab);
	sarray_init_int(verttab,nverttab);
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
	sarray_init_int(edgetab,edgenbr);
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
	//sarray_printScreen_int(verttab,vnodnbr+velmnbr+1 ,"verttab");
	//sarray_printScreen_int(edgetab,edgenbr ,"edgetab");
	//clear nodal_offsets
	nodal_offsets = (int *) tl_free(sizeof(int), vnodnbr, nodal_offsets);
	// Set seed and reset SCOTCH random number generator to produce
    // deterministic partitions on repeated calls
    //doesnt seem to work
    err = SCOTCH_meshBuild (&meshdat, velmbas, vnodbas, velmnbr, vnodnbr, verttab, NULL,
	NULL, NULL, NULL, edgenbr, edgetab);

    return meshdat;
 }

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Renumbers nodes to minimize bandwidth using SCOTCH library, by default
 *             uses Gibbs-Poole algorithm.
 *             Should only be called for serial runs. Upon completion, node ids will
 *             be permuted and an inverse permutation array will be filled out 
 *             for convenient output.
 *             Designed for only serial and initial mesh, will not update be valid if grid is
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
 int reverse_cuthill_mckee(SGRID *grid, int *permtab, int *peritab){
    //WORKS BUT DOESNT SCALE!!
    //turn grid into graph first
    int i,j,k;
    int err=0; //error code
    int n_connections;
    int nelems3d = grid->nelems3d;
    int nelems2d = grid->nelems2d;
    int nelems1d = grid->nelems1d;
    int current_node, other_node;
    int count=0;
    //nodes always start at 0
    int vertnbr=grid->nnodes; //number of vertices in graph same as nodes on mesh
    int nnodes_on_elem;
    int *nodes; //for aliasing
    int Nedges=0;
    int *n_con;
    int *n_con_no_dup;
    int *edgetab;
    int **temp_edgetab;
    int *verttab;
    //quickly compute number of edges in the mesh
    //the nodes will not be as trivial, need to know how many elements each node is connected to
    //loop through elements and add up
    //need nodal array to store this

    //really one in the same as slin_sys_init_sparsity_mono
    //need two nodal arrays
    //these will store the number of nodes each node is connected to
    n_con = (int *) tl_alloc(sizeof(int), vertnbr);
    n_con_no_dup = (int *) tl_alloc(sizeof(int), vertnbr);
    verttab = (int *) tl_alloc(sizeof(int),vertnbr+1);
    sarray_init_int(n_con, vertnbr);
    sarray_init_int(n_con_no_dup, vertnbr);
    //First set of loops is solely to establish how many nodes are connected to each node
    for (i=0;i<nelems3d;i++){
        //nnodes on the element
        nnodes_on_elem = grid->elem3d[i].nnodes;
        //Could get stuff for node weights from physics mat
        //for this element, find nodes
        nodes = grid->elem3d[i].nodes;
        //keep sparsity pattern in temp_cols_diag, temp_cols_off_diag
        for (j=0;j<nnodes_on_elem;j++){
            //each node inside an element, add this to the count
            current_node = nodes[j];
            //loop through each node that is not the node itself
            //nodes before current node
            for (k=0;k<j;k++){
                //other_node = nodes[k];
                n_con[current_node]+=1;
            }
            for(k=j+1;k<nnodes_on_elem;k++){
                //other_node = nodes[k];
                n_con[current_node]+=1;
            }
        }
    }
    for (i=0;i<nelems2d;i++){
        //nnodes on the element
        nnodes_on_elem = grid->elem2d[i].nnodes;
        //Could get stuff for node weights from physics mat
        //for this element, find nodes
        nodes = grid->elem2d[i].nodes;
        //keep sparsity pattern in temp_cols_diag, temp_cols_off_diag
        for (j=0;j<nnodes_on_elem;j++){
            //each node inside an element, add this to the count
            current_node = nodes[j];
            //loop through each node that is not the node itself
            //nodes before current node
            for (k=0;k<j;k++){
                //other_node = nodes[k];
                n_con[current_node]+=1;
            }
            for(k=j+1;k<nnodes_on_elem;k++){
                //other_node = nodes[k];
                n_con[current_node]+=1;
            }
        }
    }
    for (i=0;i<nelems1d;i++){
        //nnodes on the element
        nnodes_on_elem = grid->elem1d[i].nnodes;
        //Could get stuff for node weights from physics mat
        //for this element, find nodes
        nodes = grid->elem1d[i].nodes;
        //keep sparsity pattern in temp_cols_diag, temp_cols_off_diag
        for (j=0;j<nnodes_on_elem;j++){
            //each node inside an element, add this to the count
            current_node = nodes[j];
            //loop through each node that is not the node itself
            //nodes before current node
            for (k=0;k<j;k++){
                //other_node = nodes[k];
                n_con[current_node]+=1;
            }
            for(k=j+1;k<nnodes_on_elem;k++){
                //other_node = nodes[k];
                n_con[current_node]+=1;
            }
        }
    }
    //use nnz_rows to dynamically allocate
    //int temp_cols_diag[nrows][nCon3d];
    temp_edgetab = (int**) tl_alloc(sizeof(int*), vertnbr);
    for(j=0;j<vertnbr;j++){
        temp_edgetab[j] = (int*) tl_alloc(sizeof(int), n_con[j]);
        for(k=0;k<n_con[j];k++){
            temp_edgetab[j][k]=INT_MAX;
        }
    }
    //Seems redundant but must reuse as indexing
    sarray_init_int(n_con, vertnbr);
    //First set of loops is solely to establish how many nodes are connected to each node
    for (i=0;i<nelems3d;i++){
        //nnodes on the element
        nnodes_on_elem = grid->elem3d[i].nnodes;
        //Could get stuff for node weights from physics mat
        //for this element, find nodes
        nodes = grid->elem3d[i].nodes;
        //keep sparsity pattern in temp_cols_diag, temp_cols_off_diag
        for (j=0;j<nnodes_on_elem;j++){
            //each node inside an element, add this to the count
            current_node = nodes[j];
            //loop through each node that is not the node itself
            //nodes before current node
            for (k=0;k<j;k++){
                other_node = nodes[k];
                temp_edgetab[current_node][n_con[current_node]]=other_node;
                n_con[current_node]+=1;
            }
            for(k=j+1;k<nnodes_on_elem;k++){
                other_node = nodes[k];
                temp_edgetab[current_node][n_con[current_node]]=other_node;
                n_con[current_node]+=1;
            }
        }
    }
    for (i=0;i<nelems2d;i++){
        //nnodes on the element
        nnodes_on_elem = grid->elem2d[i].nnodes;
        //Could get stuff for node weights from physics mat
        //for this element, find nodes
        nodes = grid->elem2d[i].nodes;
        //keep sparsity pattern in temp_cols_diag, temp_cols_off_diag
        for (j=0;j<nnodes_on_elem;j++){
            //each node inside an element, add this to the count
            current_node = nodes[j];
            //loop through each node that is not the node itself
            //nodes before current node
            for (k=0;k<j;k++){
                other_node = nodes[k];
                temp_edgetab[current_node][n_con[current_node]]=other_node;
                n_con[current_node]+=1;
            }
            for(k=j+1;k<nnodes_on_elem;k++){
                other_node = nodes[k];
                temp_edgetab[current_node][n_con[current_node]]=other_node;
                n_con[current_node]+=1;
            }
        }
    }
    for (i=0;i<nelems1d;i++){
        //nnodes on the element
        nnodes_on_elem = grid->elem1d[i].nnodes;
        //Could get stuff for node weights from physics mat
        //for this element, find nodes
        nodes = grid->elem1d[i].nodes;
        //keep sparsity pattern in temp_cols_diag, temp_cols_off_diag
        for (j=0;j<nnodes_on_elem;j++){
            //each node inside an element, add this to the count
            current_node = nodes[j];
            //loop through each node that is not the node itself
            //nodes before current node
            for (k=0;k<j;k++){
                other_node = nodes[k];
                temp_edgetab[current_node][n_con[current_node]]=other_node;
                n_con[current_node]+=1;
            }
            for(k=j+1;k<nnodes_on_elem;k++){
                other_node = nodes[k];
                temp_edgetab[current_node][n_con[current_node]]=other_node;
                n_con[current_node]+=1;
            }
        }
    }
    for (i=0;i<vertnbr;i++){
        // sort the column indices (j-entries)
        //use stdlib.h qsort
        qsort(temp_edgetab[i], n_con[i], sizeof(int), compare_ints);
        //this should hopefully remove duplicates?
        n_connections = sarray_unique_int(temp_edgetab[i], n_con[i]);
        //overwrite nnz row with sarray_unique_int?
        n_con_no_dup[i] = n_connections;
        //add nnz in a row to the NNZ
        //maybe overwrire nnz_rows[i] instead if we want this stored, and then sum it
        Nedges+=n_connections;
    }
    //allocate adjacency data
    edgetab = (int *) tl_alloc(sizeof(int), Nedges);
    int *nodal_edges; //alias

    //now use info to fill in indptr and cols
    for(i=0;i<vertnbr;i++){
        //printf("filling in index ptr and column entries, row %d\n",i);
        verttab[i] = count;
        //printf("filling in index ptr and column entries, row %d\n",i);
        nodal_edges = temp_edgetab[i];
        //think about how to do this since each temp_rows may be different size
        for(j=0;j<n_con_no_dup[i];j++){
            edgetab[count] = nodal_edges[j];
            count++;
        }
    }
    //also last entry
    verttab[vertnbr] = count;

    //get rid of unwanted data
    //freeing memory
    for(j=0;j<vertnbr;j++){
        temp_edgetab[j] = (int*) tl_free(sizeof(int), n_con[j],temp_edgetab[j]);
    }
    temp_edgetab = (int**) tl_free(sizeof(int*), vertnbr,temp_edgetab);
    n_con = (int *) tl_free(sizeof(int), vertnbr, n_con);

    //NOW WE HAVE MESH IN GRAPH FORMAT
    //verifying:
    //sarray_printScreen_int(verttab,vertnbr+1,"verttab");
    //sarray_printScreen_int(n_con_no_dup,vertnbr,"degrees");
    //sarray_printScreen_int(edgetab,Nedges,"edgetab");
    //OVERWRITE FOR DEMONSTRATION
//    sarray_init_value_int(verttab,vertnbr+1,UNSET_INT);
//    sarray_init_value_int(n_con_no_dup,vertnbr,UNSET_INT);
//    sarray_init_value_int(edgetab,Nedges,UNSET_INT);
    //overwrite
//    vertnbr = 5;
//    Nedges = 12;
//    verttab[0] = 0; verttab[1] = 2;
//    verttab[2] = 5;
//    verttab[3] = 7;
//    verttab[4] = 10;
//    verttab[5] = 12;
//    n_con_no_dup[0] = 2;
//    n_con_no_dup[1] = 3;
//    n_con_no_dup[2] = 2;
//    n_con_no_dup[3] = 3;
//    n_con_no_dup[4] = 2;
//    edgetab[0] = 1;
//    edgetab[1] = 2;
//    edgetab[2] = 0;
//    edgetab[3] = 3;
//    edgetab[4] = 4;
//    edgetab[5] = 0;
//    edgetab[6] = 3;
//    edgetab[7] = 1;
//    edgetab[8] = 2;
//    edgetab[9] = 4;
//    edgetab[10] = 1;
//    edgetab[11] = 3;
    //example 2
//    vertnbr = 10;
//    Nedges = 30;
//    n_con_no_dup[0] = 3;
//    n_con_no_dup[1] = 4;
//    n_con_no_dup[2] = 2;
//    n_con_no_dup[3] = 3;
//    n_con_no_dup[4] = 5;
//    n_con_no_dup[5] = 2;
//    n_con_no_dup[6] = 3;
//    n_con_no_dup[7] = 2;
//    n_con_no_dup[8] = 3;
//    n_con_no_dup[9] = 3;
//    verttab[0] = 0; verttab[1] = 3;
//    verttab[2] = 7;
//    verttab[3] = 9;
//    verttab[4] = 12;
//    verttab[5] = 17;
//    verttab[6] = 19;
//    verttab[7] = 22;
//    verttab[8] = 24;
//    verttab[9] = 27;
//    verttab[10] = 30;
//    
//    edgetab[0] = 1;
//    edgetab[1] = 6;
//    edgetab[2] = 8;//

//    edgetab[3] = 0;
//    edgetab[4] = 4;
//    edgetab[5] = 6;
//    edgetab[6] = 9;//

//    edgetab[7] = 4;
//    edgetab[8] = 6;//

//    edgetab[9] = 4;
//    edgetab[10] = 5;
//    edgetab[11] = 8;//

//    edgetab[12] = 1;
//    edgetab[13] = 2;
//    edgetab[14] = 3;
//    edgetab[15] = 5;
//    edgetab[16] = 9;//

//    edgetab[17] = 3;
//    edgetab[18] = 4;//

//    edgetab[19] = 0;
//    edgetab[20] = 1;
//    edgetab[21] = 2;//

//    edgetab[22] = 8;
//    edgetab[23] = 9;//

//    edgetab[24] = 0;
//    edgetab[25] = 3;
//    edgetab[26] = 7;//

//    edgetab[27] = 1;
//    edgetab[28] = 4;
//    edgetab[29] = 7;
    //vertnbr = n vertices
    //verttab = starting indices to edge table for each vertex has
    //n_con_no_dup = has #vertices each vertex is connected to
    //edgetab - the table that has each of the connections in the graph
    //Nedges - length of edgetab
    int nr =0; //keeps track of how many nodes are in permtab, count will keep track of Q
    int val_found = 0;
    //find max degrees in graph
    int n_connections_max = sarray_max_int(n_con_no_dup, vertnbr);
    int temp_nodes[n_connections_max];
    int temp_degrees[n_connections_max];
    int temp_indeces[n_connections_max];
    //NOW RCM begins , returns 0 if succesful
    //Step 0: prepare empty result array permtab and empty queue Q
    //assumes permtab is already provided
    //permtab = tl_alloc(sizeof(int),vertnbr);
    int *Q = tl_alloc(sizeof(int),vertnbr); //shouldnt ever get bigger than this?
    //initalize to 0 and count to 0
    count = 0;
    sarray_init_value_int(Q, vertnbr, UNSET_INT);
    sarray_init_value_int(permtab, vertnbr, UNSET_INT);
    //Step 1, select node with lowest degree using n_con_no_dup
    //For now we will be lazy and pick first one in array
    //but may be wiser to use random number
    i = sarray_argmin_int(n_con_no_dup, vertnbr);
    //Step 2, insert P ino permtab
    permtab[nr] = i;
    nr+=1;
    //Step 3, add all neighboring nodes of i to Q
    n_connections = n_con_no_dup[i];
    int start_ind, end_ind;
    start_ind = verttab[i];
    end_ind = verttab[i+1];
    assert(n_connections == (end_ind-start_ind));
    //add all neighboring nodes of i in increasing degree!
    //need to get degrees of each node first
    for(i=0,j=start_ind;j<end_ind;i++,j++){
        temp_nodes[i] = edgetab[j];
        //use this to get number of connections
        temp_degrees[i] = n_con_no_dup[temp_nodes[i]];
    }
    //sort the temporary degrees on each node
    //temp indeces will hace the order
    sarray_argsort_int(temp_degrees, temp_indeces, n_connections);
    //sarray_reverse_argsort_int(temp_degrees, temp_indeces, n_connections);
    //shuffle the nodes using the index array
    sarray_shuffle_int(temp_nodes, temp_indeces, n_connections);
    //add the new nodes to Q
    for(i=0;i<n_connections;i++){
        Q[count] = temp_nodes[i];
        count++;
    }
    do
    {
        //Step 4.1: Extract first node from Q
        i = Q[0];
        count--;
        //shuffle down elements in Q
        for(j=0;j<count;j++){
            Q[j] = Q[j+1];
        }

        //Step 4.2: Check if i isnt already in permtab
        val_found =  sarray_is_in_int(permtab, nr, i);
        //If it is new add all neighbors of i
        if (!val_found){
            permtab[nr] = i;
            nr+=1;
            //add to Q all neighbors of i that are not in permtab
            //in ascending order of their degree
            n_connections = n_con_no_dup[i];
            start_ind = verttab[i];
            end_ind = verttab[i+1];
            for(i=0,j=start_ind;j<end_ind;i++,j++){
                temp_nodes[i] = edgetab[j];
                //use this to get number of connections
                temp_degrees[i] = n_con_no_dup[temp_nodes[i]];
            }
            //sort the temporary degrees on each node
            //temp indeces will hace the order
            sarray_argsort_int(temp_degrees, temp_indeces, n_connections);
            //sarray_reverse_argsort_int(temp_degrees, temp_indeces, n_connections);
            //shuffle the nodes using the index array
            sarray_shuffle_int(temp_nodes, temp_indeces, n_connections);
            //add the new nodes to Q if they arent already in permtab
            for(i=0;i<n_connections;i++){
                j = temp_nodes[i];
                val_found =  sarray_is_in_int(permtab, nr, j);
                if (!val_found){
                    Q[count] = j;
                    count++;
                }
            }
        }
    //debug assert Q never is too long
//#ifdef _DEBUG
//    assert(count < vertnbr);
//#endif
        //printf("%d Nodes Complete\n",nr);
    //repeat until Q is empty
    }while((nr!=vertnbr));
    //for reverse cuthill mckee reverse order of permtab
    //for regular, dont change
    start_ind = 0;
    end_ind = vertnbr - 1;
    while (start_ind < end_ind) {
        // Swap elements at start and end
        i = permtab[start_ind];
        permtab[start_ind] = permtab[end_ind];
        permtab[end_ind] = i;
        // Move towards the center
        start_ind++;
        end_ind--;
    }
    //use permtab to get peritab (inverse)
    for (i = 0; i < vertnbr; i++) {
        peritab[permtab[i]] = i; // The value arr[i] becomes the index in the inverse, and i is the value at that index
    }
    //freeing memory here
    verttab = (int *) tl_free(sizeof(int), vertnbr+1, verttab);
    edgetab = (int *) tl_free(sizeof(int), Nedges,edgetab);
    n_con_no_dup = (int *) tl_free(sizeof(int), vertnbr,n_con_no_dup);
    Q = tl_free(sizeof(int),vertnbr,Q);
    return 0;
}





   




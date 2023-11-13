/* This is the NSM Wrapper
   GSavant NSM V3.0 */
/* May 18th 2015 */

// #include "assert.h"
// #include "string.h"
// #include "define.h"
// #include "type.h"
// #include "C.Enumerations.h"
// #include "C.Q2.Channel.Cell.h"
#include "nsmv7_0.h"
#include "global_header.h"

#define YES 1       

/* STRUCT DECLARATIONS */
typedef struct {
  unsigned long number_cells;
  NSMChannelCell** cell_list;
  NSMNutrientProperties* nutrientProperties;
  NSM_AQUATIC_ENVIRONMENT* aqEnv;     /* Ptr to an aquatic environment for all channel cells*/
  NSMSedimentEnvironment* sedEnv;     /* Ptr to a sediment environment for all channel cells */
} CHANNEL_SYSTEM;

/* Prototype Declerations */
CHANNEL_SYSTEM*  CHANNEL_SYSTEM_Create();
void CHANNEL_SYSTEM_ApplyKinetics(CHANNEL_SYSTEM*  system, NSMMetEnvironment*  me);
void CHANNEL_SYSTEM_PrintNO3Concs(CHANNEL_SYSTEM* system);
CHANNEL_SYSTEM*  CHANNEL_SYSTEM_Delete(CHANNEL_SYSTEM* system);
NSM_AQUATIC_ENVIRONMENT nsm_ae;

/* Delete System */
CHANNEL_SYSTEM*  CHANNEL_SYSTEM_Delete(CHANNEL_SYSTEM* system) {
	unsigned long m;
	NSMChannelCell*  cell;
    cell = system->cell_list[0];
    cell = NSMChannelCell_Delete(cell);
    assert(cell==0);
	free(system->cell_list);
	system->nutrientProperties = NSMNutrientProperties_Delete(system->nutrientProperties);
	if(system->aqEnv != 0)
	{
		system->aqEnv = NSM_AQUATIC_ENVIRONMENT_Delete(system->aqEnv);	
	}
	if(system->sedEnv != 0)
	{
		system->sedEnv = NSMSedimentEnvironment_Delete(system->sedEnv);
	}
	free(system);
	system = 0;
	return system;
}

void nsmwq (SMODEL *mod)
{
    /* Assume timestep in seconds */
    int c, n, m, e;
	int i, j,k, ne ;
    int eid;        // Element id
    int nid;        // Node id
    int mid;        // Material id
	int nodesonelem;
	double water_head;
	double head;
	double *volume;
    size_t blkSize;
	//aliasing 
	int nmat = mod->grid->nmat;
	int nnode = mod->grid->nnodes;
	double timestep = mod->dt;

    if (mod->flag.SW2_FLOW == TRUE)
    nodesonelem = 3;
	else if (mod->flag.SW3_FLOW == TRUE)
	nodesonelem = 4;

	int nelem;
	if (mod->flag.SW2_FLOW == TRUE)
    nelem = mod->grid->nelems2d;
	else if (mod->flag.SW3_FLOW == TRUE)
	nelem = mod->grid->nelems3d;


	volume = (double *) tl_alloc(sizeof(double), nnode);

    /* NSM */
    double dmdt1;
	static NSMMetEnvironment*  nsm_me;
	static NSMChannelCell*  nsm_cell;
    static NSMNutrientProperties*  nsm_np;
    static NSM_AQUATIC_ENVIRONMENT*  nsm_ae;
    static NSM_AQUATIC_ENVIRONMENT**  nsm_ae_by_mat_id = 0;
    static NSM_AQUATIC_ENVIRONMENT**  nsm_ae_by_node_id = 0;


	for (i = 0; i < nnode; i ++)
		volume[i] = 0.;


	for (i = 0; i < nnode; i ++)
	{
		for (j = 0; j < nelem; j ++)
			for (k = 0; k < nodesonelem; k ++)
			{
				if (nodesonelem == 3){
					if (i == mod->grid->elem2d[j].nodes[k])
				       volume[i] += (1./3.) * mod->grid->elem2d[j].djac;
				}
				if (nodesonelem == 4){
					if (i == mod->grid->elem3d[j].nodes[k])
				       volume[i] += (1./4.) * mod->grid->elem3d[j].djac;
				}
			}
	}



	/* Step 1 INitialize NSM cells and environment */
	 nsm_ae = NSM_AQUATIC_ENVIRONMENT_Create();
     nsm_cell = NSMChannelCell_Create();
	 nsm_me = NSMMetEnvironment_Create();

    /* Step 2: Apply NSM to all nodes */
    /* Update NSM meteorlogical struct */
    nsm_me->airTemp.specified = 1;		/* Alias for "true"*/
	nsm_me->airTemp.value = 20.0;		/* C*/
	nsm_me->windSpeed.specified = 1;	/* Alias for "true"*/
	nsm_me->windSpeed.value = 0.2;  	/* m/s*/
	nsm_me->julianDay = 0;
	nsm_me->timeStepInDays = timestep/86400.0;

	NSMMetEnvironment_Update(nsm_me);

    for (n =0 ; n<nnode; n++)
    {

		/* Transfer Properties to Nodes from materals */
		for (e = 0; e < nelem; e ++)
			for (ne = 0; ne < nodesonelem; ne ++)
			{
				if (nodesonelem == 3)
				if (n == mod->grid->elem2d[e].nodes[ne])
				{
					
					m = mod->grid->elem2d[e].mat;
				
			nsm_ae->alpha0.rc = mod->mat[m].wnsm->alpha0;
            nsm_ae->alpha1.rc = mod->mat[m].wnsm->alpha1;
            nsm_ae->alpha2.rc = mod->mat[m].wnsm->alpha2;
            nsm_ae->alpha3.rc = mod->mat[m].wnsm->alpha3;
            nsm_ae->alpha4.rc = mod->mat[m].wnsm->alpha4;
            nsm_ae->alpha5.rc = mod->mat[m].wnsm->alpha5;
            nsm_ae->alpha6.rc = mod->mat[m].wnsm->alpha6;

            nsm_ae->beta1.rc = mod->mat[m].wnsm->beta1;
            nsm_ae->beta2.rc = mod->mat[m].wnsm->beta2;
            nsm_ae->beta3.rc = mod->mat[m].wnsm->beta3;
            nsm_ae->beta4.rc = mod->mat[m].wnsm->beta4;

            nsm_ae->K1.rc = mod->mat[m].wnsm->k1;
            nsm_ae->K2.rc = mod->mat[m].wnsm->k2;
			
            nsm_ae->K3.rc = mod->mat[m].wnsm->k3;
            nsm_ae->K4.rc = mod->mat[m].wnsm->k4;
            nsm_ae->mu.rc = mod->mat[m].wnsm->mu;

            nsm_ae->rho.rc = mod->mat[m].wnsm->rho;

            nsm_ae->sigma1.rc = mod->mat[m].wnsm->sigma1;
            nsm_ae->sigma2.rc = mod->mat[m].wnsm->sigma2;
            nsm_ae->sigma3.rc = mod->mat[m].wnsm->sigma3;
            nsm_ae->sigma4.rc = mod->mat[m].wnsm->sigma4;
            nsm_ae->sigma5.rc = mod->mat[m].wnsm->sigma5;

            nsm_ae->Kl = mod->mat[m].wnsm->kl;
            nsm_ae->Kn = mod->mat[m].wnsm->kn;
            nsm_ae->Kp = mod->mat[m].wnsm->kp;
            nsm_ae->Pn = mod->mat[m].wnsm->pn;

            nsm_ae->lambda0 = mod->mat[m].wnsm->lambda0;
            nsm_ae->lambda1 = mod->mat[m].wnsm->lambda1;
            nsm_ae->lambda2 = mod->mat[m].wnsm->lambda2;

            nsm_ae->min_do_conc = mod->mat[m].wnsm->min_do_conc;

				}
				if (nodesonelem == 4)
				if (n == mod->grid->elem3d[e].nodes[ne])
				{
					
					m = mod->grid->elem3d[e].mat;
				
			nsm_ae->alpha0.rc = mod->mat[m].wnsm->alpha0;
            nsm_ae->alpha1.rc = mod->mat[m].wnsm->alpha1;
            nsm_ae->alpha2.rc = mod->mat[m].wnsm->alpha2;
            nsm_ae->alpha3.rc = mod->mat[m].wnsm->alpha3;
            nsm_ae->alpha4.rc = mod->mat[m].wnsm->alpha4;
            nsm_ae->alpha5.rc = mod->mat[m].wnsm->alpha5;
            nsm_ae->alpha6.rc = mod->mat[m].wnsm->alpha6;

            nsm_ae->beta1.rc = mod->mat[m].wnsm->beta1;
            nsm_ae->beta2.rc = mod->mat[m].wnsm->beta2;
            nsm_ae->beta3.rc = mod->mat[m].wnsm->beta3;
            nsm_ae->beta4.rc = mod->mat[m].wnsm->beta4;

            nsm_ae->K1.rc = mod->mat[m].wnsm->k1;
            nsm_ae->K2.rc = mod->mat[m].wnsm->k2;
			
            nsm_ae->K3.rc = mod->mat[m].wnsm->k3;
            nsm_ae->K4.rc = mod->mat[m].wnsm->k4;
            nsm_ae->mu.rc = mod->mat[m].wnsm->mu;

            nsm_ae->rho.rc = mod->mat[m].wnsm->rho;

            nsm_ae->sigma1.rc = mod->mat[m].wnsm->sigma1;
            nsm_ae->sigma2.rc = mod->mat[m].wnsm->sigma2;
            nsm_ae->sigma3.rc = mod->mat[m].wnsm->sigma3;
            nsm_ae->sigma4.rc = mod->mat[m].wnsm->sigma4;
            nsm_ae->sigma5.rc = mod->mat[m].wnsm->sigma5;

            nsm_ae->Kl = mod->mat[m].wnsm->kl;
            nsm_ae->Kn = mod->mat[m].wnsm->kn;
            nsm_ae->Kp = mod->mat[m].wnsm->kp;
            nsm_ae->Pn = mod->mat[m].wnsm->pn;

            nsm_ae->lambda0 = mod->mat[m].wnsm->lambda0;
            nsm_ae->lambda1 = mod->mat[m].wnsm->lambda1;
            nsm_ae->lambda2 = mod->mat[m].wnsm->lambda2;

            nsm_ae->min_do_conc = mod->mat[m].wnsm->min_do_conc;

				}
			
			}
			

		if (mod->flag.SW2_FLOW == TRUE)
			head = mod->sw->d2->head[n];
		if (mod->flag.SW3_FLOW == TRUE)
			head = 10.; /* THis is filler for now */

        if (head > 0)    /* Temp. Replace with Partial element integration GS */
        {

            nsm_cell->depth = head;           /* M*/
			if (nodesonelem == 3)
            nsm_cell->vol = volume[n] * head;                              /* M^3*/
			else
			{
			nsm_cell->vol = volume[n];
			nsm_cell->depth = 1.;
			}

				/*nsm_cell->vol = volume[n];*/
            nsm_cell->waterTemp = 25;              /* Deg C*/

            nsm_cell->aqSpecies.list[0]->conc = 0;
			nsm_cell->aqSpecies.list[1]->conc = 0;
			nsm_cell->aqSpecies.list[2]->conc = 0;
			nsm_cell->aqSpecies.list[3]->conc = 0;
			nsm_cell->aqSpecies.list[4]->conc = 0;
			nsm_cell->aqSpecies.list[5]->conc = 0;
			nsm_cell->aqSpecies.list[6]->conc = 0;
			nsm_cell->aqSpecies.list[7]->conc = 0;
			nsm_cell->aqSpecies.list[8]->conc = 0;
			

			nsm_cell->aqSpecies.list[0]->mass = 0;
			nsm_cell->aqSpecies.list[1]->mass = 0;
			nsm_cell->aqSpecies.list[2]->mass = 0;
			nsm_cell->aqSpecies.list[3]->mass = 0;
			nsm_cell->aqSpecies.list[4]->mass = 0;
			nsm_cell->aqSpecies.list[5]->mass = 0;
			nsm_cell->aqSpecies.list[6]->mass = 0;
			nsm_cell->aqSpecies.list[7]->mass = 0;
			nsm_cell->aqSpecies.list[8]->mass = 0;

			for (i = 0; i < mod->ntransport; i ++)
			{
              if (mod->con[i].type == NOO) nsm_cell->aqSpecies.list[0]->conc = MAX(mod->con[i].concentration[n], 0);
			  if (mod->con[i].type == NOOO) nsm_cell->aqSpecies.list[1]->conc = MAX(mod->con[i].concentration[n], 0);
			  if (mod->con[i].type == NHHHH) nsm_cell->aqSpecies.list[2]->conc = MAX(mod->con[i].concentration[n], 0);
			  if (mod->con[i].type == OrN) nsm_cell->aqSpecies.list[3]->conc = MAX(mod->con[i].concentration[n], 0);
			  if (mod->con[i].type == OP) nsm_cell->aqSpecies.list[4]->conc = MAX(mod->con[i].concentration[n], 0);
			  if (mod->con[i].type == DP) nsm_cell->aqSpecies.list[5]->conc = MAX(mod->con[i].concentration[n], 0);
			  if (mod->con[i].type == ALG) nsm_cell->aqSpecies.list[6]->conc = MAX(mod->con[i].concentration[n], 0);
			  if (mod->con[i].type == CBOD) nsm_cell->aqSpecies.list[7]->conc = MAX(mod->con[i].concentration[n], 0);
			  if (mod->con[i].type == DO) nsm_cell->aqSpecies.list[8]->conc = MAX(mod->con[i].concentration[n], 0);
			  if (mod->con[i].type == NOO) nsm_cell->aqSpecies.list[0]->mass = MAX(mod->con[i].concentration[n] * volume[n], 0);
			  if (mod->con[i].type == NOOO) nsm_cell->aqSpecies.list[1]->mass = MAX(mod->con[i].concentration[n] * volume[n], 0);
			  if (mod->con[i].type == NHHHH) nsm_cell->aqSpecies.list[2]->mass = MAX(mod->con[i].concentration[n] * volume[n], 0);
			  if (mod->con[i].type == OrN) nsm_cell->aqSpecies.list[3]->mass = MAX(mod->con[i].concentration[n] * volume[n], 0);
			  if (mod->con[i].type == OP) nsm_cell->aqSpecies.list[4]->mass = MAX(mod->con[i].concentration[n] * volume[n], 0);
			  if (mod->con[i].type == DP) nsm_cell->aqSpecies.list[5]->mass = MAX(mod->con[i].concentration[n] * volume[n], 0);
			  if (mod->con[i].type == ALG) nsm_cell->aqSpecies.list[6]->mass = MAX(mod->con[i].concentration[n] * volume[n], 0);
			  if (mod->con[i].type == CBOD) nsm_cell->aqSpecies.list[7]->mass = MAX(mod->con[i].concentration[n] * volume[n], 0);
			  if (mod->con[i].type == DO) nsm_cell->aqSpecies.list[8]->mass = MAX(mod->con[i].concentration[n] * volume[n], 0);
			}
		


            /* Execute MNSM kinetics on current "cell" =- actually, a node. */
			
            NSMChannelCell_Kinetics(nsm_cell, nsm_me, nsm_ae);
			

	        /* Update all constituent Concs  */
		   {
			   if (mod->bc_mask[n] != YES)
			   for (i =0; i < mod->ntransport ; i++)
			   {
				   if(mod->con[i].type == NOO) mod->con[i].concentration[n] += nsm_cell->aqSpecies.list[0]->dmdt * (timestep / nsm_cell->vol);
				   if(mod->con[i].type == NOOO) mod->con[i].concentration[n] += nsm_cell->aqSpecies.list[1]->dmdt * (timestep / nsm_cell->vol);
				   if(mod->con[i].type == NHHHH) mod->con[i].concentration[n] += nsm_cell->aqSpecies.list[2]->dmdt * (timestep / nsm_cell->vol);
				   if(mod->con[i].type == OrN) mod->con[i].concentration[n] += nsm_cell->aqSpecies.list[3]->dmdt * (timestep / nsm_cell->vol);
				   if(mod->con[i].type == OP) mod->con[i].concentration[n] += nsm_cell->aqSpecies.list[4]->dmdt * (timestep / nsm_cell->vol);
				   if(mod->con[i].type == DP) mod->con[i].concentration[n] += nsm_cell->aqSpecies.list[5]->dmdt * (timestep / nsm_cell->vol);
				   if(mod->con[i].type == ALG) mod->con[i].concentration[n] += nsm_cell->aqSpecies.list[6]->dmdt * (timestep / nsm_cell->vol);
				   if(mod->con[i].type == CBOD) mod->con[i].concentration[n] += nsm_cell->aqSpecies.list[7]->dmdt * (timestep / nsm_cell->vol);
				   /*if (mod->flag.SW3_FLOW == TRUE)
				   if(nsm[n].type == 1)
				       if(mod->con[i].type == DO) mod->con[i].concentration[n] += nsm_cell->aqSpecies.list[8]->dmdt * (timestep / nsm_cell->vol);
				   if (mod->flag.SW2_FLOW == TRUE) */
					   if(mod->con[i].type == DO) mod->con[i].concentration[n] += nsm_cell->aqSpecies.list[8]->dmdt * (timestep / nsm_cell->vol);
			   }


           }
        }
    }
}

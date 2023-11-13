#ifndef H_SSTRUCTURES_
#define H_SSTRUCTURES_

typedef struct
{
  int nodeu,egsu,upnode;                                    /* upstream node and edge string */
  int noded,egsd,dnode;                                    /* downstream node and edge string*/
  double cd,z1,z2,w1,w2,sl1,sl2,hw;
  double fluxu,fluxd;                           /* flux through edge string */
} SWEIR_C; 

typedef struct
{
  int nodeu,egsu,upnode;                                    /* upstream node and edge string */
  int noded,egsd,dnode;                                    /* downstream node and edge string*/
  double A, B, C, D, E, F;
  double GG;   /* Length of Flap Gate*/
  double fluxu,fluxd;                           /* flux through edge string */
  double density, length, width;
  int calc;
} SFLAP_C; 

typedef struct
{
	double a;
	int opening;
	int egsu, egsd;
	double fluxu, fluxd;
} SSLUICE_C;


/*********************************************************/
/* struct methods -------------------------------------- */

void sstructures_init_weirs(SWEIR_C **, int , int);
void sstructures_init_flaps(SFLAP_C **, int , int);
void sstructures_init_sluices(SSLUICE_C **, int , int);
void sstructures_free_weirs(SWEIR_C *, int , int);
void sstructures_free_flaps(SFLAP_C *, int , int);
void sstructures_free_sluices(SSLUICE_C *, int , int);

#endif

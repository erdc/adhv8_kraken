#include "global_header.h"

#ifdef _ADH_BREACH
#include "breach_type.h"
#endif

/********************************************************/
/********************************************************/

void sstructures_init_weirs(SWEIR_C **weirs, int nstring, int nweir) {

  
  int istruct;
  (*weirs) = (SWEIR_C *) tl_alloc(sizeof(SWEIR_C), nweir);
  SWEIR_C *wstr = *(weirs);   // alias
   if (nweir > 0){
  for (istruct = 0; istruct < nweir; istruct++) {
    wstr[istruct].egsu = UNSET_INT;
	wstr[istruct].egsd = UNSET_INT;
	wstr[istruct].w1 = 0.;
	wstr[istruct].hw = 0.;
	wstr[istruct].z1 = 0.;
    }
   }

 }

void sstructures_init_flaps(SFLAP_C **flaps, int nstring, int nflap) {
	
	
	int istruct;
    (*flaps) = (SFLAP_C *) tl_alloc(sizeof(SFLAP_C), nflap);
    SFLAP_C *fstr = *(flaps);   // alias
   if (nflap > 0){
  for (istruct = 0; istruct < nflap; istruct++) {
    fstr[istruct].egsu = UNSET_INT;
	fstr[istruct].egsd = UNSET_INT;
	fstr[istruct].A = 0.;
	fstr[istruct].B = 0.;
	fstr[istruct].C = 0.;
	fstr[istruct].D = 0.;
	fstr[istruct].GG = 0.;
    }
   }
}


void sstructures_init_sluices(SSLUICE_C **sluices, int istring, int nsluice){

    
	int istruct;
    (*sluices) = (SSLUICE_C *) tl_alloc(sizeof(SSLUICE_C), nsluice);
    SSLUICE_C *ssluc = *(sluices);   // alias
   if (nsluice > 0){
  for (istruct = 0; istruct < nsluice; istruct++) {
    ssluc[istruct].egsu = UNSET_INT;
	ssluc[istruct].egsd = UNSET_INT;
	ssluc[istruct].opening = UNSET_INT;
	ssluc[istruct].a = 0.;
    }
   }
}
//void sstr_value_print(STR_VALUE string) {

void sstructures_free_weirs(SWEIR_C *weirs, int nstring, int nweir) {
  if (nweir > 0)
  weirs = (SWEIR_C *) tl_free(sizeof(SWEIR_C), nweir, weirs);
}
void sstructures_free_flaps(SFLAP_C *flaps, int nstring, int nflap){
  if (nflap > 0)
  flaps = (SFLAP_C *) tl_free(sizeof(SFLAP_C), nflap, flaps);
}
void sstructures_free_sluices(SSLUICE_C *sluices, int nstring, int nsluice){
  if (nsluice > 0)
  sluices = (SSLUICE_C *) tl_free(sizeof(SSLUICE_C), nsluice, sluices);
}

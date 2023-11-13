//////////////////////////////////////////////////////////////////////////////////////////////////
// FOLLOWING LINES ADDED BY GAJANAN [ gkc July 2015 ]. These are for ADCIRC wind library usage. //

#ifdef WINDLIB

/* Note: Entire file added by Gajanan */ 
#include "global_header.h"
#ifdef WIND_GNU_NM
void __wind_MOD_timeconv(int *iyr, int *imo, int *iday, int *ihr, int *imin, double *isec,
                         double *timesec, int *MyProc, int *NScreen, int *ScreenUnit);
#endif
#ifdef WIND_INTEL_NM
void wind_mp_timeconv_(int *iyr, int *imo, int *iday, int *ihr, int *imin, double *isec,
                       double *timesec, int *MyProc, int *NScreen, int *ScreenUnit);
#endif
#ifdef WIND_PGI_NM
void wind_timeconv_(int *iyr, int *imo, int *iday, int *ihr, int *imin, double *isec,
                    double *timesec, int *MyProc, int *NScreen, int *ScreenUnit);
#endif
void read_bc_WINDLIB(SMODEL *mod, char *data) {

  // allocate windlib struct
  swindlib_alloc(&(mod->windlib));

  char *subdata = NULL;         /* the data after the second card   is read */
  SIO info = *(mod->io);        // alias

  // open and initialize sedlib input file
  init_adh_in_file(&(mod->io->windfile));
  sprintf(mod->io->windfile.filename, "fort.22");
  open_input_file(&(mod->io->windfile), "windlib file", strlen(mod->io->sup.filename));
  assert(mod->io->windfile.fp);

  parse_card(data, &subdata);

  mod->windlib->ics = 1; // Default, Cartesian! Spherical (ICS=2) NOT ALLOWED yet.
  mod->windlib->nws = read_int_field(info, &subdata);  // 1. NWS type
  if (mod->windlib->nws!=1){
    mod->windlib->statim = read_dbl_field(info, &subdata);  // 2. Wind data starting time
    mod->windlib->wtiminc = read_dbl_field(info, &subdata);  // 3. Wind data time increment
  }
  // cjt :: store in mod->flags
  mod->flag.NWS = mod->windlib->nws;

  switch (mod->windlib->nws) {
  case 1:
#ifdef _MESSG
    printf("NWS = %i is not allowed to run in parallel. Please compile AdH in serial.",mod->windlib->nws);
    exit(0);
#else
// Nothing required here; wind loading is mesh specific, and wind time stepping matches simulation time stepping.
#endif
    return;
  case 2:
#ifdef _MESSG
    printf("NWS = %i is not allowed to run in parallel. Please compile AdH in serial.",mod->windlib->nws);
    exit(0);
#else
// Nothing required here; wind loading is mesh specific, but interpolation of wind in time required.
#endif
    return;
  case 3:
    mod->windlib->nwlon    = read_int_field(info, &subdata);
    mod->windlib->nwlat    = read_int_field(info, &subdata);
    mod->windlib->wlonmin  = read_dbl_field(info, &subdata);
    mod->windlib->wlatmax  = read_dbl_field(info, &subdata);
    mod->windlib->wloninc  = read_dbl_field(info, &subdata);
    mod->windlib->wlatinc  = read_dbl_field(info, &subdata);

    mod->windlib->irefyr   = read_int_field(info, &subdata); 
    mod->windlib->irefmo   = read_int_field(info, &subdata); 
    mod->windlib->irefday  = read_int_field(info, &subdata); 
    mod->windlib->irefhr   = read_int_field(info, &subdata); 
    mod->windlib->irefmin  = read_int_field(info, &subdata); 
    mod->windlib->refsec   = read_dbl_field(info, &subdata); 
#ifdef WIND_GNU_NM
    __wind_MOD_timeconv( &(mod->windlib->irefyr),  &(mod->windlib->irefmo),
                         &(mod->windlib->irefday), &(mod->windlib->irefhr),
                         &(mod->windlib->irefmin), &(mod->windlib->refsec),
                         &(mod->windlib->wreftim), &(mod->windlib->MyProc),
                         &(mod->windlib->NScreen), &(mod->windlib->ScreenUnit) );
#endif
#ifdef WIND_INTEL_NM
    wind_mp_timeconv_( &(mod->windlib->irefyr),  &(mod->windlib->irefmo),
                       &(mod->windlib->irefday), &(mod->windlib->irefhr),
                       &(mod->windlib->irefmin), &(mod->windlib->refsec),
                       &(mod->windlib->wreftim), &(mod->windlib->MyProc),
                       &(mod->windlib->NScreen), &(mod->windlib->ScreenUnit) );
#endif
#ifdef WIND_PGI_NM
    wind_timeconv_( &(mod->windlib->irefyr),  &(mod->windlib->irefmo),
                    &(mod->windlib->irefday), &(mod->windlib->irefhr),
                    &(mod->windlib->irefmin), &(mod->windlib->refsec),
                    &(mod->windlib->wreftim), &(mod->windlib->MyProc),
                    &(mod->windlib->NScreen), &(mod->windlib->ScreenUnit) );
#endif
    break;
  case 4:
#ifdef _MESSG
    printf("NWS = %i is not allowed to run in parallel. Please compile AdH in serial.",mod->windlib->nws);
    exit(0);
#else
// Nothing required here; wind loading is mesh specific, but interpolation of wind in time required.
#endif
    break;
  case 5:
#ifdef _MESSG
    printf("NWS = %i is not allowed to run in parallel. Please compile AdH in serial.",mod->windlib->nws);
    exit(0);
#else
// Nothing required here; wind loading is mesh specific, but interpolation of wind in time required.
#endif
    return;
  case 6:
    mod->windlib->nwlon    = read_int_field(info, &subdata);
    mod->windlib->nwlat    = read_int_field(info, &subdata);
    mod->windlib->wlonmin  = read_dbl_field(info, &subdata);
    mod->windlib->wlatmax  = read_dbl_field(info, &subdata);
    mod->windlib->wloninc  = read_dbl_field(info, &subdata);
    mod->windlib->wlatinc  = read_dbl_field(info, &subdata);
    break;
  case 7:
    mod->windlib->nwlon    = read_int_field(info, &subdata);
    mod->windlib->nwlat    = read_int_field(info, &subdata);
    mod->windlib->wlonmin  = read_dbl_field(info, &subdata);
    mod->windlib->wlatmax  = read_dbl_field(info, &subdata);
    mod->windlib->wloninc  = read_dbl_field(info, &subdata);
    mod->windlib->wlatinc  = read_dbl_field(info, &subdata);
    break;
  case 8:
    printf("AdH reading ERA data\n");
    return;
  default:
    printf("\n NWS = %i not allowed", mod->windlib->nws);
    printf("\n Error reading 'NWS': Only NWS values {1, 2, 3, 4, 5, 6, 7, 8} allowed! Terminating program.\n");
    exit(0);
    break;
  }
  fclose(mod->io->windfile.fp); /* Executed only for NWS 3, 4, 6, 7 */
}

#endif
// ABOVE LINES ADDED BY GAJANAN                                                                 //
//////////////////////////////////////////////////////////////////////////////////////////////////

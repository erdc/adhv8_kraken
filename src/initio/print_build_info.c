#include "global_header.h"
#include "adh_version.h"

/*!
   \brief Print ADH build information to screen
 */
void print_build_info()
{
#ifdef _PRE_ADH
  printf("\n ---------------\n pre-AdH Version\n ---------------\n");
#else
  printf("\n -----------\n AdH Version\n -----------\n");
#endif

  printf(" Version %s\n", ADH_VERSION);

#ifdef _ADH_SEDIMENT
  printf(" Compatible with Sedlib v%s\n", SEDLIB_VERSION);
#endif

#ifdef _PRE_ADH
  printf("\n -------------------------\n pre-AdH Build Information\n -------------------------\n");
#else
  printf("\n ---------------------\n AdH Build Information\n ---------------------\n");
#endif

//  printf(" SVN revision # %s\n", ADHREV);
//  printf(" Build Date/Time: %s / %s\n", ADHREVDATE, ADHREVTIME);

  /* debug */
#ifdef _DEBUG
  printf(" Built with DEBUG enabled\n");
#endif

#ifdef _DEBUG_KITCHEN_SINK
  printf(" Built with DEBUG_KITCHEN_SINK enabled\n");
#endif

#ifdef _DEBUG_TEMPORARY
  printf(" Built with DEBUG_TEMPORARY enabled\n");
#endif

  /* physics */
#ifdef _ADH_GROUNDWATER
  printf(" Built with GW physics enabled\n");
#endif

#ifdef _ADH_HEAT
  printf(" Built with HEAT physics enabled\n");
#endif

#ifdef _ADH_NS
  printf(" Built with NS physics enabled\n");
#endif

#ifdef _ADH_OVERLAND
  printf(" Built with OVERLAND physics enabled\n");
#endif

#ifdef _ADH_REACT
  printf(" Built with REACT physics enabled\n");
#endif

#ifdef _ADH_SW2
  printf(" Built with SW2 physics enabled\n");
#endif

#ifdef _ADH_SW3
  printf(" Built with SW3 physics enabled\n");
#endif

  /* output */
#ifdef _ADH_HDF5
  printf(" Built with XDMF (*.xml with *.h5) output file format\n");
#elif defined _ADH_BINARY
  printf(" Built with XMDF (XMS binary) output file format\n");
#else
  printf(" Built with XMS output file format\n");
#endif

  /* extras modules */
#ifdef _ADH_HEC
  printf(" Built with HEC DSS file support enabled\n");
#endif

#ifdef _SOCKETS
  printf(" Built with SOCKETS enabled\n");
#endif

#ifdef _ADH_STRUCTURES
  printf(" Built with STRUCTURES enabled\n");
#endif

#ifdef _ADH_BREACH
  printf(" Built with BREACH enabled\n");
#endif

#ifdef UNREF_CONSV
  printf(" Built with UNREFINEMENT CONSERVATION enabled\n");
#endif

  /* adh libraries */
#ifdef _ADH_ICM
  printf(" Built with ICM library support enabled\n");
#endif

#ifdef _ADH_NSM
  printf(" Built with NSM library support enabled\n");
#endif

#ifdef _ADH_SEDIMENT
  printf(" Built with Sedlib support enabled\n");
#endif

  /* external libraries */
#ifdef _MPI
  printf(" Built with MPI enabled\n");
#endif

#ifdef _METIS
  printf(" Built with PARMETIS enabled\n");
#endif

#ifdef _UMFPACK
#if defined _UMFPACK_VERSION && defined _UMFPACK_INT_SIZE
  printf(" Built with UMFPACK Version %d (int size: %d bytes) enabled\n",
    _UMFPACK_VERSION, _UMFPACK_INT_SIZE);
#else
  printf(" Built with UMFPACK enabled\n");
#endif
#endif
  
  return;
}

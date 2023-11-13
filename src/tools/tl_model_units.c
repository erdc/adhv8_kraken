#include "global_header.h"

void tl_model_units(SMODEL *mod) {
    if (!mod->flag.UNITS) {
        if (IS_ABOUT(mod->density, 1000., 4.) == YES && IS_ABOUT(mod->gravity, 9.8, 4.) ==   YES) {
            mod->flag.UNITS = MKS;
        }
        else if (IS_ABOUT(mod->density, 1.94, 4.) == YES && IS_ABOUT(mod->gravity, 32.2, 4.) == YES) {
            mod->flag.UNITS = FPS;
            printf(" Calculations are in FPS, temperature is expected to be in Fahrenheit.   \n");
        }
#ifdef _ADH_GROUNDWATER
	else if (mod->flag.GW_FLOW) {
	  printf("Running GW Model. Units could not be determined\n");
	}
#endif
        else {
            tl_error(">> density and gravity don't correspond to MKS or FPS\n");
        }
    }
}

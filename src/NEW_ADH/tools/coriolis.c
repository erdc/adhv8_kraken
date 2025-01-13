#include "adh.h"
#define EARTH_ROTAT 7.2722E-5   /* Earth's rotation in 1/s (coriolis: added 6-02) */
double get_coriolis_angular_speed(double coriolis_factor) {
    return (sin((PI / 180.) * coriolis_factor) * 2. * EARTH_ROTAT);
}

/* this uses a simplified linear routine to calcuate density */
/* this is in Kg-m-s-C units */
/* the salinity is in psu */
/* if either the temperature or the salinity is left out then */
/* 20 degrees C and 0 psu are assumed */
/* EOS: A summary concerning the newly adopted practical salinity scale, 1978 and the international equation of state of seawater, 1980 */
/* Pritchard */
/* Implemented GSavant 02192017*/

#include "global_header.h"

void tl_density_calculator_metric(double density,  double *temperature, double stmp, double *salinity, double ssal, int number_entries, double *rho, int flag) {
  
  int i;
  double salt_temp, tmp_temp;
  double a0, a1, a2, a3, a4, a5;
  double b0, b1, b2, b3, b4, b5;
  double c0, c1, c2;
  double d0;
  
  
  a0 = 999.842594;
  a1 = 0.06793952;
  a2 = -0.009095290;
  a3 = 0.0001001685;
  a4 = -0.000001120083;
  a5 = 0.000000006536332;
  b0 = 0.824493;
  b1 = -0.0040899;
  b2 = 0.000076438;
  b3 = -0.00000082467;
  b4 = 0.0000000053875;
  c0 = -0.00572466;
  c1 = 0.00010227;
  c2 = -0.0000016546;
  d0 = 0.00048314;
  
  if (flag == 0) /* Use linearized EOS */
    {
      if (temperature == NULL) {
	for(i = 0; i < number_entries; i++) {
	  salt_temp = MIN(MAX(salinity[i] * ssal,0.0),50.);
	  rho[i] = (density / 998.2) * (1025.5 + 0.78 * (salt_temp  - 35.));
	}
      } else if(salinity == NULL) {
	for(i = 0; i < number_entries; i++) {
	  tmp_temp = MIN(MAX(temperature[i] * stmp,0.0),50.);
	  rho[i] = (density / 998.2) * (999.7 - 0.15 * (tmp_temp - 10.));
	}
      } else {
	for(i = 0; i < number_entries; i++) {
	  salt_temp = MIN(MAX(salinity[i] * ssal,0.0),50.);
	  tmp_temp = MIN(MAX(temperature[i] * stmp,0.0),50.);
	  rho[i] = (density / 998.2) * (1027. - 0.15 * (tmp_temp  - 10.) + 0.78 * (salt_temp  - 35.));
	}
      }
    }
  if (flag == 1) { /* Use full EOS */
    if (temperature == NULL) {
      for(i = 0; i < number_entries; i++) {
	salt_temp = MIN(MAX(salinity[i] * ssal,0.0),50.);
	rho[i] = a0 + (b0)*salt_temp + (c0)*pow(salt_temp, 3./2.) + (d0)*pow(salt_temp,2.);
	rho[i] = rho[i] * (density / 998.2);
      }
    }
    else if(salinity == NULL) {
      for(i = 0; i < number_entries; i++) {
	tmp_temp = MIN(MAX(temperature[i] * stmp,0.0),50.);
	rho[i] = (density / 998.2) * (999.7 - 0.15 * (tmp_temp - 10.));
      }
    }
    else {
      for(i = 0; i < number_entries; i++) {
	salt_temp = MIN(MAX(salinity[i] * ssal,0.0),50.);
	tmp_temp = MIN(MAX(temperature[i] * stmp,0.0),50.);
	rho[i] = a0 + a1*tmp_temp + a2 * pow(tmp_temp, 2.) + a3 * pow(tmp_temp, 3.) + a4 * pow(tmp_temp, 4.) + a5 * pow(tmp_temp, 5.) +
	  (b0 + b1 * tmp_temp + b2 * pow (tmp_temp, 2.) + b3 * pow(tmp_temp, 3.) + b4 * pow(tmp_temp, 4.)) * salt_temp +
	  (c0 + c1 * tmp_temp + c2 * pow(tmp_temp, 2.)) * pow(salt_temp, 3./2.) + d0 * pow(salt_temp, 2.);
	rho[i] = rho[i] * (density / 998.2);
      }
    }
  }
#ifdef _ADH_GROUNDWATER
  if (flag == 2) { /* GW salinity */
    for (i=0; i < number_entries; i++) {
      rho[i] = density;
    }
    if (salinity != NULL) {
      for (i=0; i < number_entries; i++) {
	rho[i] = density*(1.0 + 0.00078141 * salinity[i]);
      }
    }
  } else if (flag == 3) { /* GW Heat */
    for (i=0; i < number_entries; i++) {
      rho[i] = density;
    }
    if (temperature != NULL) {
      for (i=0; i < number_entries; i++) {
	rho[i] = density*(1.003 - 0.0015027 * temperature[i]);
      }
    }
  }
  
#endif
}

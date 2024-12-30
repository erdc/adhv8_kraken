#include "adh.h"

double l2_error(double *v1, double *v2, int n){
	double err = 0.0;
	double temp;
	int i; 
	for (i = 0; i < n; i++) {
		temp = v1[i] - v2[i];
		err += temp*temp;
	}
	return sqrt(err);
}

double linf_error(double *v1, double *v2, int n){
	double err = 0.0;
	double temp;
	int i; 
	for (i = 0; i < n; i++) {
		temp = fabs(v1[i] - v2[i]);
		if ( temp > err) {
            err = temp;
        }
    }
	return err;


}

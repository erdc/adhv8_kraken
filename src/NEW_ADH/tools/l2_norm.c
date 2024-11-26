#include "adh.h"
double l2_norm(double *v, int size){
	double l2norm = 0.0;
	for(int i=0;i<size;i++){
		l2norm+= v[i]*v[i];
	}
	return sqrt(l2norm);
}

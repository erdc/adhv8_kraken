#include "global_header.h"

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

double read_adh_velocity_frame(FILE *fp, SGRID *grid, SVECT *vel, int unit_conversion, bool normalize) {
    
    assert(fp!=NULL);
    
    int i,idum;
    char *line = NULL, card[10];
    size_t len = 0;
    ssize_t readfile = -1;
    double time = -1.0, vx, vy, vz;
    
    while(1) {
        readfile = getline(&line, &len, fp); // TS 0 time
        sscanf(line,"%s %d %lf",card,&idum,&time);
        if (strcmp(card, "TS") == 0) {
            for (i=0; i<grid->nnodes; i++) {
                readfile = getline(&line, &len, fp);
                //printf("node %d :: line = %s\n",i,line);
                if (readfile == -1) {
                    printf("ERROR :: velocity file end reached\n");
                    exit(-1);
                }
                sscanf(line,"%lf %lf %lf",&vx,&vy,&vz);
                vel[i].x = vx;
                vel[i].y = vy;
                vel[i].z = vz;
                if (normalize) {
                    vel[i].x *= grid->xLinv;
                    vel[i].y *= grid->yLinv;
                }
            }
           return time * unit_conversion;
        }
        if (strcmp(card, "ENDDS") == 0) {
            printf("** read_adh_velocity_frame :: Velocity E.O.F.\n");
            break;
        }
    }
    return -1;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

double read_adh_wq_frame(FILE *fp, SGRID *grid, double *oxy, double *sal, double *sun, int unit_conversion) {
    
    assert(fp!=NULL);
    
    int i,idum;
    char *line = NULL, card[10];
    size_t len = 0;
    ssize_t readfile = -1;
    double time = -1.0, vx, vy, vz;
    
    while(1) {
        readfile = getline(&line, &len, fp); // TS 0 time
        sscanf(line,"%s %d %lf",card,&idum,&time);
        if (strcmp(card, "TS") == 0) {
            for (i=0; i<grid->nnodes; i++) {
                readfile = getline(&line, &len, fp);
                //printf("node %d :: line = %s\n",i,line);
                if (readfile == -1) {
                    printf("ERROR :: velocity file end reached\n");
                    exit(-1);
                }
                sscanf(line,"%lf %lf %lf",&vx,&vy,&vz);
                oxy[i] = vx;
                sal[i] = vy;
                sun[i] = vz;
            }
            return time * unit_conversion;
        }
    }
    return -1;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

double read_adh_dpl_frame(FILE *fp, SGRID *grid, double *dpl, int unit_conversion) {
    
    assert(fp!=NULL);
    
    int i,idum;
    char *line = NULL, card[10];
    size_t len = 0;
    ssize_t readfile = -1;
    double time = -1.0, d;
    
    while(1) {
        readfile = getline(&line, &len, fp); // TS 0 time
        sscanf(line,"%s %d %lf",card,&idum,&time);
        if (strcmp(card, "TS") == 0) {
            for (i=0; i<grid->nnodes; i++) {
                readfile = getline(&line, &len, fp);
                //printf("node %d :: line = %s\n",i,line);
                d=0.;
                if (readfile == -1) {
                    printf("ERROR :: velocity file end reached\n");
                    exit(-1);
                }
                sscanf(line,"%lf",&d);
                dpl[i] = d;
            }
            return time * unit_conversion;
        }
    }
    return -1;
}

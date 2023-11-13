//////////////////////////////////////////////////////////////////////////////////////////////////
// FOLLOWING LINES ADDED BY GAJANAN [ gkc July 2015 ]. These are for ADCIRC wind library usage. //

#include "global_header.h"


/***********************************************************/
/***********************************************************/
/***********************************************************/

void swindlib_alloc(SWINDLIB **windlib) {
    (*windlib) = (SWINDLIB *) tl_alloc(sizeof(SWINDLIB), 1);
}

void swindlib_init(SWINDLIB *wl, int nnodes) {
    int i;

    wl->prevx = (double *) tl_alloc(sizeof(double), nnodes);
    wl->prevy = (double *) tl_alloc(sizeof(double), nnodes);
    wl->prevp = (double *) tl_alloc(sizeof(double), nnodes);
    wl->currx = (double *) tl_alloc(sizeof(double), nnodes);
    wl->curry = (double *) tl_alloc(sizeof(double), nnodes);
    wl->currp = (double *) tl_alloc(sizeof(double), nnodes);
    for(i = 0; i<nnodes; i++){
       wl->prevx[i] = 0.0; wl->prevy[i] = 0.0; wl->prevp[i] = 0.0;
       wl->currx[i] = 0.0; wl->curry[i] = 0.0; wl->currp[i] = 0.0;
    }
}

/***********************************************************/
/***********************************************************/
/***********************************************************/

void swindlib_free(SWINDLIB *windlib, int nnodes) {
    windlib->prevx = (double *) tl_free(sizeof(double), nnodes, windlib->prevx);
    windlib->prevy = (double *) tl_free(sizeof(double), nnodes, windlib->prevy);
    windlib->prevp = (double *) tl_free(sizeof(double), nnodes, windlib->prevp);
    windlib->currx = (double *) tl_free(sizeof(double), nnodes, windlib->currx);
    windlib->curry = (double *) tl_free(sizeof(double), nnodes, windlib->curry);
    windlib->currp = (double *) tl_free(sizeof(double), nnodes, windlib->currp);

    windlib = (SWINDLIB *) tl_free(sizeof(SWINDLIB), 1, windlib);
}

/***********************************************************/
/***********************************************************/
/***********************************************************/

void swindlib_printScreen(SWINDLIB *windlib) {

    printf("\n-----------------\nWind library data\n-----------------\n");
    printf("NWS:          %i\n", windlib->nws);
    printf("ICS:          %i\n", windlib->ics);
    if (windlib->nws !=1){
        printf("STATIM:       %f\n", windlib->statim);
        printf("WTIMINC:      %f\n\n", windlib->wtiminc);
    }

    if (windlib->nws == 3 || windlib->nws == 6 || windlib->nws == 7 ){
      printf("\nAdditional parameters for NWS = %i -----------------\n", windlib->nws);
      printf("NWLON:        %i\n", windlib->nwlon);
      printf("NWLAT:        %i\n", windlib->nwlat);
      printf("WLONMIN:      %f\n", windlib->wlonmin);
      printf("WLATMAX:      %f\n", windlib->wlatmax);
      printf("WLONINC:      %f\n", windlib->wloninc);
      printf("WLATINC:      %f\n\n", windlib->wlatinc);
    }
    if (windlib->nws == 3){
      printf("IREFYR:       %i\n", windlib->irefyr);
      printf("IREFMO:       %i\n", windlib->irefmo);
      printf("IREFDAY:      %i\n", windlib->irefday);
      printf("IREFHR:       %i\n", windlib->irefhr);
      printf("IREFMIN:      %i\n", windlib->irefmin);
      printf("REFSEC:       %f\n", windlib->refsec);
      printf("WREFTIM:      %f\n\n", windlib->wreftim);
    }
}

/***********************************************************/
/***********************************************************/
/***********************************************************/

void swindlib_update(SMODEL *mod, SWINDLIB *wl, double t_curr, int nnodes, double **values) {
    int i,id3d;
    double *x, *y, wtratio;
    double slam = 0; double sfea = 0;
    x = (double *) tl_alloc(sizeof(double), nnodes);
    y = (double *) tl_alloc(sizeof(double), nnodes);
    if (mod->grid->ndim == 2) {
        for (i=0; i<nnodes; i++) {
            x[i] = mod->grid->node[i].x;
            y[i] = mod->grid->node[i].y;
        }
    } else if (mod->grid->ndim == 3) {
        for (i=0; i<nnodes; i++) {
            id3d = mod->grid->nodeID_2d_to_3d_sur[i];
            x[i] = mod->grid->node[id3d].x;
            y[i] = mod->grid->node[id3d].y;
        }
    }
    switch (wl->nws){
        case 1: // Wind STRESSES read in directly at all grid nodes. No spatial / temporal interpolaiton.
            AdH_direct_nwsget(wl->currx, wl->curry, wl->currp, nnodes, mod->io);
            for (i=0;i<nnodes;i++){
                values[0][i] = wl->currx[i];
                values[1][i] = wl->curry[i];
                //			values[2][i] = wl->currp[i];
            }
            x = (double *) tl_free(sizeof(double), nnodes, x);
            y = (double *) tl_free(sizeof(double), nnodes, y);
            return;
            break;
        case 2: // Wind STRESSES read in directly at all grid nodes. Interpolated temporally.
            if (t_curr > wl->wtime_curr){
                wl->wtime_prev = wl->wtime_curr;
                wl->wtime_curr += wl->wtiminc;
                for (i=0;i<nnodes;i++){
                    wl->prevx[i] = wl->currx[i];
                    wl->prevy[i] = wl->curry[i];
                    wl->prevp[i] = wl->currp[i];
                }
                AdH_direct_nwsget(wl->currx, wl->curry, wl->currp, nnodes, mod->io);
            }
            break;
        case 3: // Wind VELOCITIES in US Navy Fleet format interpolated spaciotemporally.
            if (t_curr > wl->wtime_curr){
                wl->wtime_prev = wl->wtime_curr;
                wl->wtime_curr += wl->wtiminc;
                for (i=0;i<nnodes;i++){
                    wl->prevx[i] = wl->currx[i];
                    wl->prevy[i] = wl->curry[i];
                }
                windlib_nws3get( x, y, &slam, &sfea, wl->currx, wl->curry, &(wl->iwtime),
                                   &(wl->iwyr), &(wl->wtimed), &(nnodes), &(wl->nwlon),
                                   &(wl->nwlat), &(wl->wlatmax), &(wl->wlonmin), &(wl->wlatinc),
                                   &(wl->wloninc), &(wl->ics), &(wl->NScreen), &(wl->ScreenUnit) );
            }
            break;
        case 4: // Wind VELOCITIES read in PBL/JAG format directly at all grid nodes. Interpolated temporally.
            if (t_curr > wl->wtime_curr){
                wl->wtime_prev = wl->wtime_curr;
                wl->wtime_curr += wl->wtiminc;
                for (i=0;i<nnodes;i++){
                    wl->prevx[i] = wl->currx[i];
                    wl->prevy[i] = wl->curry[i];
                    wl->prevp[i] = wl->currp[i];
                }
                windlib_nws4get(wl->currx,wl->curry, wl->currp,
                                   &(nnodes), &(mod->density), &(mod->gravity));
            }
            break;
        case 5: // Wind VELOCITIES read in directly at all grid nodes. Interpolated temporally.
            if (t_curr > wl->wtime_curr){
                wl->wtime_prev = wl->wtime_curr;
                wl->wtime_curr += wl->wtiminc;
                for (i=0;i<nnodes;i++){
                    wl->prevx[i] = wl->currx[i];
                    wl->prevy[i] = wl->curry[i];
                    wl->prevp[i] = wl->currp[i];
                }
                AdH_direct_nwsget(wl->currx, wl->curry, wl->currp, nnodes, mod->io);
            }
            break;
        case 6: // Wind VELOCITIES read in from a rectangular grid, and interpolated spaciotemporally.
            if (t_curr > wl->wtime_curr){
                wl->wtime_prev = wl->wtime_curr;
                wl->wtime_curr += wl->wtiminc;
                for (i=0;i<nnodes;i++){
                    wl->prevx[i] = wl->currx[i];
                    wl->prevy[i] = wl->curry[i];
                    wl->prevp[i] = wl->currp[i];
                }
                windlib_nws6get( x, y, &slam, &sfea,
                                   wl->currx, wl->curry, wl->currp, &(nnodes), &(wl->nwlon),
                                   &(wl->nwlat), &(wl->wlatmax), &(wl->wlonmin), &(wl->wlatinc),
                                   &(wl->wloninc), &(wl->ics), &(mod->density), &(mod->gravity) );
            }
            break;
        case 7: // Wind STRESSES read in from a rectangular grid, and interpolated spaciotemporally.
            if (t_curr > wl->wtime_curr){
                wl->wtime_prev = wl->wtime_curr;
                wl->wtime_curr += wl->wtiminc;
                for (i=0;i<nnodes;i++){
                    wl->prevx[i] = wl->currx[i];
                    wl->prevy[i] = wl->curry[i];
                    wl->prevp[i] = wl->currp[i];
                }
                windlib_nws7get( x, y, &slam, &sfea,
                                   wl->currx, wl->curry, wl->currp, &(nnodes), &(wl->nwlon),  
                                   &(wl->nwlat), &(wl->wlatmax), &(wl->wlonmin), &(wl->wlatinc),
                                   &(wl->wloninc), &(wl->ics), &(mod->density), &(mod->gravity) );
            }
            break;
        case 8: // ERA Wind VELOCITIES read in directly at all grid nodes. Interpolated temporally.
            if (t_curr > wl->wtime_curr){
                wl->wtime_prev = wl->wtime_curr;
                wl->wtime_curr += wl->wtiminc;
                for (i=0;i<nnodes;i++){
                    wl->prevx[i] = wl->currx[i];
                    wl->prevy[i] = wl->curry[i];
                    wl->prevp[i] = wl->currp[i];
                }
                AdH_ERA_nwsget(wl->currx, wl->curry, wl->currp, mod->grid, mod->io);
            }
            break;
        default:
            printf("\nNot sure how you got here, but NWS = %i is not allowed!\n",wl->nws);
            printf("\nOnly NWS = {1, 2, 3, 4, 5, 6, 7, 8} are allowed. Terminating program.\n");
            printf("******************************************************************\n");
            exit(0);
            break;
    }
    wtratio = (t_curr - wl->wtime_prev)/wl->wtiminc;
    for (i=0;i<nnodes;i++){
        values[0][i] = wl->prevx[i] + wtratio * (wl->currx[i] - wl->prevx[i]);
        values[1][i] = wl->prevy[i] + wtratio * (wl->curry[i] - wl->prevy[i]);
        //		values[2][i] = wl->prevp[i] + wtratio * (wl->currp[i] - wl->prevp[i]);
    }

    x = (double *) tl_free(sizeof(double), nnodes, x);
    y = (double *) tl_free(sizeof(double), nnodes, y);
}

/***********************************************************/
/***********************************************************/
/***********************************************************/

void swindlib_firstread(SMODEL *mod, SWINDLIB *wl, int nnodes ) {
    printf("\n*************************************************************\n");
    printf("-- reading wind forcing file: %s\n",mod->io->windfile.filename);
    int i,id3d;
    double *x, *y;
    double slam=0; double sfea=0;
    double iwtimep;
    x = (double *) tl_alloc(sizeof(double), nnodes);
    y = (double *) tl_alloc(sizeof(double), nnodes);
    if (mod->grid->ndim == 2) {
        for (i=0; i<nnodes; i++) {
            x[i] = mod->grid->node[i].x;
            y[i] = mod->grid->node[i].y;
        }
    } else if (mod->grid->ndim == 3) {
        for (i=0; i<nnodes; i++) {
            id3d = mod->grid->nodeID_2d_to_3d_sur[i];
            x[i] = mod->grid->node[id3d].x;
            y[i] = mod->grid->node[id3d].y;
        }
    }
    switch (wl->nws) {
        case 1:
            x = (double *) tl_free(sizeof(double), nnodes, x);
            y = (double *) tl_free(sizeof(double), nnodes, y);
            return;
        case 2:
            AdH_direct_nwsget(wl->prevx, wl->prevy, wl->prevp, nnodes, mod->io);
            AdH_direct_nwsget(wl->currx, wl->curry, wl->currp, nnodes, mod->io);
            break;
        case 3:
            do {
                windlib_nws3get( x, y, &slam, &sfea, wl->currx, wl->curry, &(wl->iwtime),
                                   &(wl->iwyr), &(wl->wtimed), &(nnodes), &(wl->nwlon),
                                   &(wl->nwlat), &(wl->wlatmax), &(wl->wlonmin), &(wl->wlatinc),
                                   &(wl->wloninc), &(wl->ics), &(wl->NScreen), &(wl->ScreenUnit));
                if (wl->iwyr != wl->irefyr){
                    iwtimep = wl->iwtime;
                    for (i=0; i<nnodes;i++){
                        wl->prevx[i] = wl->currx[i];
                        wl->prevy[i] = wl->curry[i];
                    }
                    continue;
                }
                if (wl->wtimed <= wl->wreftim){
                    iwtimep = wl->iwtime;
                    for (i=0; i<nnodes;i++){
                        wl->prevx[i] = wl->currx[i];
                        wl->prevy[i] = wl->curry[i];
                    }
                }
            }while (wl->iwyr != wl->irefyr || wl->wtimed <= wl->wreftim);
            iwtimep = wl->iwtime;
            for (i=0; i<nnodes;i++){
                wl->prevx[i] = wl->currx[i];
                wl->prevy[i] = wl->curry[i];
            }
            //            printf("\nFound wind data at time t1 = %f",iwtimep);
            //            printf("\nFound wind data at time t2 = %f\n",wl->iwtime);
            wl->wtime_curr = wl->wtimed - wl->wreftim;
            wl->wtime_prev = wl->wtime_curr - wl->wtiminc;
            x = (double *) tl_free(sizeof(double), nnodes, x);
            y = (double *) tl_free(sizeof(double), nnodes, y);
            printf("\nFound wind data at time t1 = %f",wl->wtime_prev);
            printf("\nFound wind data at time t2 = %f\n",wl->wtime_curr);
            return;
        case 4:
            windlib_nws4get(wl->prevx, wl->prevy, wl->prevp,
                               &(nnodes), &(mod->density), &(mod->gravity));
            windlib_nws4get(wl->currx,wl->curry, wl->currp,
                               &(nnodes), &(mod->density), &(mod->gravity));
            
            break;
        case 5:
            AdH_direct_nwsget(wl->prevx, wl->prevy, wl->prevp, nnodes, mod->io);
            AdH_direct_nwsget(wl->currx, wl->curry, wl->currp, nnodes, mod->io);
            break;
        case 6:
            windlib_nws6get( x, y, &slam, &sfea,
                               wl->prevx,  wl->prevy, wl->prevp, &(nnodes), &(wl->nwlon),
                               &(wl->nwlat), &(wl->wlatmax), &(wl->wlonmin), &(wl->wlatinc),
                               &(wl->wloninc), &(wl->ics), &(mod->density), &(mod->gravity) );
            windlib_nws6get( x, y, &slam, &sfea,
                               wl->currx, wl->curry, wl->currp, &(nnodes), &(wl->nwlon),
                               &(wl->nwlat), &(wl->wlatmax), &(wl->wlonmin), &(wl->wlatinc),
                               &(wl->wloninc), &(wl->ics), &(mod->density), &(mod->gravity) );
            break;
        case 7:
            windlib_nws7get( x, y, &slam, &sfea,
                               wl->prevx,  wl->prevy, wl->prevp, &(nnodes), &(wl->nwlon),
                               &(wl->nwlat), &(wl->wlatmax), &(wl->wlonmin), &(wl->wlatinc),
                               &(wl->wloninc), &(wl->ics), &(mod->density), &(mod->gravity) );
            windlib_nws7get( x, y, &slam, &sfea,
                               wl->currx, wl->curry, wl->currp, &(nnodes), &(wl->nwlon),
                               &(wl->nwlat), &(wl->wlatmax), &(wl->wlonmin), &(wl->wlatinc),
                               &(wl->wloninc), &(wl->ics), &(mod->density), &(mod->gravity) );
            break;
        case 8:
            AdH_ERA_nwsget(wl->prevx, wl->prevy, wl->prevp, mod->grid, mod->io);
            AdH_ERA_nwsget(wl->currx, wl->curry, wl->currp, mod->grid, mod->io);
            break;
        default:
            printf("\nNot sure how you got here, but NWS = %i is not allowed!\n",wl->nws);
            printf("\nOnly NWS = {1, 2, 3, 4, 5, 6, 7, 8} are allowed. Terminating program.\n");
            printf("*************************************************************\n\n");
            exit(0);
            break;
    }

    x = (double *) tl_free(sizeof(double), nnodes, x);
    y = (double *) tl_free(sizeof(double), nnodes, y);

    wl->wtime_prev = wl->statim;
    wl->wtime_curr = wl->wtime_prev + wl->wtiminc;
    printf("Wind forcing starter calls completed successfully for NWS = %i\n",wl->nws);
    printf("*************************************************************\n\n");
}

void AdH_direct_nwsget(double *wvnx, double *wvny, double *press, int nnodes, SIO * io){
    int i, nodenum;
    char line[MAXLINE];
    char *data;
    for (i=0; i<nnodes; i++){
        fgets(line, MAXLINE, io->windfile.fp);
        data = line;
        parse_card(line,&data);
        nodenum          = read_int_field(*io, &data);
        if (nodenum > nnodes) {
            printf("\n\nProblem encountered in wind file fort.22");
            printf("Line from fort.22:  %s",line);
            printf("\nNode number encountered %i exceeds the number of surface nodes %i",nodenum,nnodes);
            printf("Terminating AdH run...");
            exit(0);
        }
        else{
            nodenum--; /* Since arrays start at [0] in C */
        }
        wvnx[nodenum]    = read_dbl_field(*io, &data);
        wvny[nodenum]    = read_dbl_field(*io, &data);
        press[nodenum]   = read_dbl_field(*io, &data);
    }
}

// ABOVE LINES ADDED BY GAJANAN                                                                 //

// cjt :: added for ERA data option || probably won't work with adaption!
void AdH_ERA_nwsget(double *wvnx, double *wvny, double *press, SGRID *grid, SIO * io){
    
    int i, j, found, found_anywhere;
    char line[MAXLINE];
    char *data;
    assert(io->windfile.fp);
    fgets(line, MAXLINE, io->windfile.fp); // TS TIME_UNITS SNAP_TIME
    //printf("line: %s \n",line);
    //printf("grid->orig_macro_nnodes: %d\n",grid->orig_macro_nnodes);
    //printf("grid->nnodes: %d\n",grid->nnodes);
    for (i=0; i<grid->orig_macro_nnodes; i++) {
        fgets(line, MAXLINE, io->windfile.fp);
        data = line;
        parse_card(line,&data);
#ifdef _MESSG
        // there is a faster way to do this I'm sure
        found = 0;
        for (j=0; j<grid->my_nnodes; j++) { // loop over local nodes
            if (grid->node[j].original_id == i) { // is the PE local node the same as the global original
                wvnx[j]  = read_dbl_field(*io, &data);
                wvny[j]  = read_dbl_field(*io, &data);
                press[j] = 0.0; //read_dbl_field(*io, &data);
                found = 1;
                break;
            }
        }
        //printf("global node: %d \t myid: %d \t found: %d\n",i,grid->smpi->myid,found);
        found_anywhere = 0;
        MPI_Reduce(&found, &found_anywhere, 1, MPI_INT, MPI_MAX, 0, grid->smpi->ADH_COMM);
        if (grid->smpi->myid == 0) {
            if (found_anywhere == 0) {
                printf("Global node: %d \n",i);
                tl_error(">> FATAL ERROR: Could not find original wind node %d on any PE subdomains");
            }
        }
#else
        wvnx[i]  = read_dbl_field(*io, &data);
        wvny[i]  = read_dbl_field(*io, &data);
        press[i] = 0.0; //read_dbl_field(*io, &data);
#endif
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////

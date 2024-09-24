/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file read_bc.c This file collects methods to read an AdH SuperModel boundary condition input file       */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
static int DEBUG = OFF;
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Reads a SuperModel boundary condition file
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] *sm (SUPER_MODEL *)  a double pointer to an array of AdH supermodels
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void read_bc(SMODEL_SUPER *mod) {
    char *line = NULL, cdum[20], line_save[100];
    size_t len = 0;
    ssize_t read;
    char *token, *subtoken, *subsubtoken;
    char delim[] = " ";
    int idum;
    double fdum;

    SIO info = *(mod->io);   // alias

    assert(mod->io->bc.fp);

#ifdef _DEBUG
    if (DEBUG) {
        printf("------------------------------------------------------\n");
        printf("reading bc file for debugging first\n");
        printf("------------------------------------------------------\n");
    }
#endif
    while ((read = getline(&line, &len, mod->io->bc.fp)) != -1) {
        get_token(line,&token);
        if (strcmp(token, "NOTERM") == 0) {
            read_bc_NOTERM(mod,&token);
        } else if (strcmp(token, "DEBUG") == 0) {
            read_bc_DEBUG(mod,&token);
        } else if (strcmp(token, "SCREEN_OUTPUT") == 0) {
            read_bc_SCREEN_OUTPUT(mod,&token);
        } else if (strcmp(token, "FILE_OUTPUT") == 0) {
            read_bc_FILE_OUTPUT(mod,&token);
        } else if (strcmp(token, "TEST") == 0) {
            read_bc_TEST(mod,&token);
        }
    }
    rewind(mod->io->bc.fp);

#ifdef _DEBUG
    if (DEBUG) {
        printf("------------------------------------------------------\n");
        printf("reading bc file for series extraction\n");
        printf("------------------------------------------------------\n");
    }
#endif
    int type = UNSET_INT;
    while ((read = getline(&line, &len, mod->io->bc.fp)) != -1) {
        get_token(line,&token);
        if (strcmp(token, "SERIES") == 0) {
            get_next_token(&token);
            if (strcmp(token, "AWRITE") == 0) {
                type = OUTPUT_SERIES;
                read_bc_AWRITE(mod,&token);
            }
            if (strcmp(token, "WIND") == 0) {
                type = WIND_SERIES;
                read_bc_SERIES(mod,&token);
            }
            if (strcmp(token, "WAVE") == 0) {
                type = WAVE_SERIES;
                read_bc_SERIES(mod,&token);
            }
            if (strcmp(token, "DT") == 0) {
                type = DT_SERIES;
                read_bc_SERIES(mod,&token);
            }
            if (strcmp(token, "WRITE") == 0) {
                type = OUTPUT_SERIES;
                read_bc_SERIES(mod,&token);
            }
#ifdef _ADH_GROUNDWATER
//            if (strcmp(token, "GWCONSTITUENT") == 0) {
//                type = CONSTITUITIVE_SERIES;
//                read_bc_SERIES(mod,&token);
//            }
#endif
        }
    }
    rewind(mod->io->bc.fp);

#ifdef _DEBUG
    if (DEBUG) {
        printf("------------------------------------------------------\n");
        printf("reading bc file for string extraction\n");
        printf("------------------------------------------------------\n");
    }
#endif
    while ((read = getline(&line, &len, mod->io->bc.fp)) != -1) {
        get_token(line,&token);
        if (strcmp(token, "NDS") == 0 || strcmp(token, "EGS") == 0 ||
            strcmp(token, "MDS") == 0 || strcmp(token, "MTS") == 0) {
            read_bc_STRINGS(mod,&token);
        }
    }
    rewind(mod->io->bc.fp);

#ifdef _DEBUG
    if (DEBUG) {
        printf("------------------------------------------------------\n");
        printf("reading rest of bc file\n");
        printf("------------------------------------------------------\n");
    }
#endif
    while ((read = getline(&line, &len, mod->io->bc.fp)) != -1) {
        get_token(line,&token);
        if (strcmp(token, "OP") == 0) {
            read_bc_OP(mod,&token);
        }
        if (strcmp(token, "IP") == 0) {
            read_bc_IP(mod,&token);
        }
        if (strcmp(token, "PC") == 0) {
            read_bc_PC(mod,&token);
        }
        if (strcmp(token, "TC") == 0) {
            read_bc_TC(mod,&token);
        }
        if (strcmp(token, "OFF") == 0) {
            read_bc_OFF(mod,&token);
        }
        if (strcmp(token, "MP") == 0) {
            read_bc_MP(mod,&token);
        }
        if (strcmp(token, "DB") == 0) {
            read_bc_DB(mod,&token);
        }
        if (strcmp(token, "NB") == 0) {
            read_bc_NB(mod,&token);
        }
//        if (strcmp(token, "HY") == 0) {
//            read_bc_HY(mod,&token);
//        }
        if (strcmp(token, "CN") == 0) {
            read_bc_CN(mod,&token);
        }
#ifdef _SEDLIB
//        if (strcmp(token, "SEDLIB") == 0) {
//            read_bc_SEDLIB(mod,line);
//        }
#endif
        if (strcmp(token, "OP") == 0) {
            read_bc_OP(mod,&token);
        }
    }
    rewind(mod->io->bc.fp);
}


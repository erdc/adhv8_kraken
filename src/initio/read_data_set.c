/* this routine reads a data set */
#include "global_header.h"

int read_data_set(
  SIO info,  
  SFILE_IN *file,		/* input file with data set */
  SGRID *grid,
  double *dbl_array,		/* the double array where we keep the data set */
  char *name			      /* the entire NAME line of the data set which may
                           include the constituent or layer ID */
)
{
  char line[MAXLINE];		/* the input line */
  char *data = NULL;		/* the data */
  CARD card, ds_type = NO_CARD;
  int inode, iarray, ival;	/* loop counters */
  int ds_nelem = UNSET_INT,	/* the number of data set elements */
      ds_nnode = UNSET_INT,	/* the number of data set nnodes */
      ds_nval = UNSET_INT,  /* the number of values per node */
      ds_parts = 0;     /* the number of data set parts read (total of 6) */
  int res = NO;

  assert(file && file->fp);
  sprintf(name, "");

  while ((res == NO) && (fgets(line, MAXLINE, file->fp) != NULL))
    {
      io_save_line(&info, file->fp, file->filename, line);
      if (strip_comments(line) <= 1) /* ignore empty ('\n') or comment lines */
        {
          continue;
        }
      card = parse_card(line, &data);
      switch (card)
        { 
          case CARD_BEGSCL: /* Regular scalar dataset */
          case CARD_BEGVEC: /* Regular vector dataset */
            ++ds_parts;
            if (ds_type != NO_CARD)
              {
                io_read_error(info, "Found extra data set begin card before the end of the "
                  "data set was reached (most likely missing ENDDS card).", TRUE);
              }
            ds_type = card;
            /* set dataset type specific attribute */
            switch (card)
              {
                case CARD_BEGSCL:
                  ds_nval = 1; /* only 1 value per node */
                  break;
                case CARD_BEGVEC:
                  ds_nval = 3; /* 3 vector componenet values per node */
                  break;
                default:
                  assert(0); /* unhandled case */
                  break;
              }
            break;
          case CARD_NAME:
            ++ds_parts;
            if (strlen(name))
              {
                io_read_error(info, "Found extra NAME card before the end of the data set "
                  "was reached.", TRUE);
              }
            /* copy the entire file line so the constituent or layer ID can be
               parsed later when this dataset is processed */
            strcpy(name, line);
            break;
          case CARD_NC:
            ++ds_parts;
            if (ds_nelem != UNSET_INT)
              {
                io_read_error(info, "Found extra NC card before the end of the data set was "
                  "reached.", TRUE);
              }
            ds_nelem = read_int_field_custom(info, &data, NULL, "number of elements", 1, TRUE);
            if (ds_nelem != grid->macro_nelems3d && ds_nelem != grid->macro_nelems2d && (ds_nelem != grid->nelems2d_sur && grid->smpi->npes ==1))
              {
                char msg[MAXLINE];
                
                sprintf(msg, "The number of geometery elements (2D: %d, 3D: %d, 3D-surface: %d)  "
                  "does not match the number of data set elements (%d).",
                  grid->macro_nelems2d, grid->macro_nelems3d, grid->nelems2d_sur, ds_nelem);
                io_read_error(info, msg, TRUE);
              }
            break;
          case CARD_ND:
            ++ds_parts;
            if (ds_nnode != UNSET_INT)
              {
                io_read_error(info, "Found extra ND card before the end of the data set was "
                  "reached.", TRUE);
              }
            ds_nnode = read_int_field_custom(info, &data, NULL, "number of nnodes", 1, TRUE);
            if ((ds_nnode != grid->macro_nnodes) && (ds_nnode != grid->macro_nnodes_sur))
              {
                char msg[MAXLINE];
                
                sprintf(msg, "The number of geometry nodes (%d) or surface geometry nodes (%d) does not match "
                  "the number of data set nodes (%d)).", grid->macro_nnodes, grid->macro_nnodes_sur,ds_nnode);
                io_read_error(info, msg, TRUE);
              }
            break;
          case CARD_TS:
            ++ds_parts;
            /* before reading data, ensure other information was read first */
            if (ds_type == NO_CARD)
              {
                io_read_error(info, "A data set begin card (e.g. BEGSCL) must be specified "
                  "prior to the data set values.", TRUE);
              }
            if (ds_nnode == UNSET_INT)
              {
                io_read_error(info, "The ND card must be specified prior to the data set "
                  "values.", TRUE);
              }
            if (ds_nelem == UNSET_INT)
              {
                io_read_error(info, "The NC card must be specified prior to the data set "
                  "values.", TRUE);
              }
            if (strlen(name) == 0)
              {
                io_read_error(info, "The NAME card must be specified prior to the data set "
                  "values.", TRUE);
              }
            /* Ignore the TS card information; we are assuming this is an initial condition */
            inode = 0;
            /* Read the nodal values */
            switch (ds_nval)
              {
                case UNSET_INT:
                  assert(0); /* ds_nval was not set above */
                  io_read_error(info, "Internal error occurred reading data set.", TRUE);
                  break;
                /* read one value per node (one value on each line) */
                case 1:
                  //while (inode < grid->nnodes && fgets(line, MAXLINE, file->fp) != NULL)
                  while (inode < ds_nnode && fgets(line, MAXLINE, file->fp) != NULL)
                    {
                      io_save_line(&info, file->fp, file->filename, line);
                      data = &line[0];
                      dbl_array[inode] = read_dbl_field_custom(info, &data, NULL,"data set value", 1, TRUE);
                      ++inode;
                    }
                  break;
                /* read multiple values per node (multiple values on each line) */
                default:
                  iarray = 0; 
                  //while (inode < grid->nnodes && fgets(line, MAXLINE, file->fp) != NULL)
                  while (inode < ds_nnode && fgets(line, MAXLINE, file->fp) != NULL)
                    {
                      
                      //assert(grid->node[inode].string != UNSET_INT);
                      io_save_line(&info, file->fp, file->filename, line);
                      data = &line[0];
                      for (ival = 0; ival < ds_nval; ++ival, ++iarray)
                        {
                          dbl_array[iarray] = read_dbl_field_custom(info, &data, NULL,"data set value", ival + 1, TRUE);
                        }
                      ++inode;
                    }
                  break;
              }
            /* check for premature end of data set */
            //if (inode != grid->nnodes)
            if (inode != ds_nnode)
              {
                io_read_error(info, "Data set ended prematurely; not enough values "
                  "were read.", TRUE);
              }
            break;
          case CARD_ENDDS:
            /* found the end of the data set */
            ++ds_parts;
            res = YES; 
            if (ds_parts != 6)
              {
                /* TS card checks for other required card, so it must be missing */
                io_read_error(info, "The TS card must be specified prior to end of data "
                  "set.", TRUE);
              }
            break;
          default:
            /* do nothing */
            break;
        }
    }

  /* ensure that we read all parts of the data set, if we found a dataset */
  if (res == NO && ds_parts > 0)
    {
      io_read_error(info, "Missing ENDDS card at end of data set.", TRUE);
    }

  return res;
}

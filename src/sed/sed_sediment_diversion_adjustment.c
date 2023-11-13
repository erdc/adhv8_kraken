/* calculates the sediment diversion adjustment */

#include "share_extern.h"

void
sed_sediment_diversion_adjustment ()	
{
  int *nchk, i, k, n, ielem1, iupdate;
  double k_n, tsrh, sd_adj;

  nchk = tl_alloc (sizeof (int), nnode);

  iupdate = 0;
  for (n = 0; n < nstring; n++)
    {
      if (str_values[n].sed_div_flag == TRUE)
	{
          iupdate = 1;
	  for (i = 0; i < my_nnode; i++)
	    nchk[i] = 0;

	  for (ielem1 = 0; ielem1 < nelem1d; ielem1++)
	    {
	      if (elem1d_flags[ielem1] == n)
		{
		  nchk[elem1d[ielem1].nodes[0]] = 1;
		  nchk[elem1d[ielem1].nodes[1]] = 1;
		}
	    }

	  for (i = 0; i < my_nnode; i++)
	    {
	      if (nchk[i] != 0)
		{
                   tsrh = fr_coef[i];
                   if (fr_flag[i] == 1)
                   {
                     k_n = 1.0;
                     if (gravity > 10.0) k_n = 1.486;
                     tsrh = pow ((8.25 * sqrt (gravity) * fr_coef[i] / k_n), 6.);
                   }
                  for (k = 0; k < nsed; k++)
                    {
                      sd_adj = sed_diversion_correction_factor(tsrh, bed_elev[i], ol_head[i], bed_displacement[i], 
                               div_coef[0], div_coef[1], div_coef[2], rouse_coef[k][i]);
		      mfcf[k][i] *= sd_adj;
                    }
		}
	    }
	}

      if (iupdate == 1)
      {

#ifdef _MESSG
                  for (k = 0; k < nsed; k++)
                    comm_update_double (mfcf[k], ONE);
#endif
      }
    }

  nchk = tl_free (sizeof (int), nnode, nchk);
}

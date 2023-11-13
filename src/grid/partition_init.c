/* this routine partitions the initial mesh */

#include "share_extern.h"

void partition_init(
  void
)
{
#ifdef _MESSG
  int i, j;			/* loop counters */
  int ifirst, ilast;		/* bounds on the loops to load the subdomains */
  int nnode_inc;		/* roughly the number of elements in a subdomain */
  int isd;			/* loop counter over the subdomains */

  /* determine how many nodes are on each processor */
  nnode_inc = nnode / npes;
  if(nnode_inc < 1)
    tl_error("Too many subdomains for the number of nodes.");

  /* loops over the subdomains */
  for(isd = 0, ilast = 0; isd < npes; isd++)
    {
      /* sets the first and last nodes of the subdomain */
      ifirst = ilast;
      if(isd == npes - 1)
	ilast = nnode;
      else
	ilast += nnode_inc;

      /* sets the number of nodes on this processor */
      if(isd == myid)
	{
	  my_nnode = ilast - ifirst;
	}

      /* assigns the nodes */
      for(i = ifirst, j = 0; i < ilast; i++, j++)
	{
	  node_pair[i].sd = isd;
	  node_pair[i].rnode = j;
	}
    }
#endif
}

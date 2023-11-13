#include "global_header.h"
#ifdef _MESSG
int test(
  double,
  int,
  double
);
void ROTATE_IT(
  int *,
  int,
  int,
  int
);


void findrank(
  EDGE_RANK *a,
  int n,
  int *rank,
  SMPI *smpi
)
{

  /*
   *  This is a parallel Odd-Even Transpostion sort to compute the rank
   *    of an array of doubles.  No actual data is moved from the original
   *    input array a, but corresponding ranks are returned in rank[].
   *    It is assumed that the elements of a[] are locally sorted.
   *
   *  A simple merge is used to sort data locally after each exchange.
   *    This requires the use of three (3) extra arrays.  Each element in the
   *    new arrays is a structure containing a double (from a[]) and an
   *    integer (initialized to the PE rank which held the value).  Once these
   *    structures have been sorted according to the data value, the global
   *    rank is computed and this rank is communicated back to the PE which
   *    held the data.
   *
   *  Author: Clay P. Breshears <clay@turing.wes.hpc.mil>
   *          CEWES MSRC
   *          Jul 23 1998
   *
   */

#define LT 1
#define GT 2

  int myid, numprocs, max1, max2, amount, p1, p2, itemsize;
  int i, tag, iteration, count;
  int odd, odd_start, odd_inc, odd_indx, odd_compare;
  int even, even_start, even_inc, even_indx, even_compare;
  int *outranks, *outrankr;
  struct sortitem
  {
    double data;
    int where;
  };
  /* cjt */
  struct sortitem *local_odd; 
  struct sortitem *local_even;
  struct sortitem *foreign;
  MPI_Status status, statusa;
  MPI_Request done1, done2;
  

  /*
   *  Find out who I am and how many PEs are involved
   */

  itemsize = sizeof(struct sortitem);
  myid=smpi->myid;
  numprocs=smpi->npes;
  MPI_Comm ADH_COMM;
  ADH_COMM = smpi->ADH_COMM;
  /*
   *  Create local arrays of sort items, one for the ODD phase
   *    and one for the EVEN phases.
   */
  if(n != 0)
    {
      local_odd = (struct sortitem *)tl_alloc(itemsize, n);
      local_even = (struct sortitem *)tl_alloc(itemsize, n);
    }

  /*
   *  Initialize local array, set up neighbor addresses and other
   *    constants to assist in computation.
   */

  for(i = 0; i < n; i++)
    {
      local_odd[i].data = a[i].length;
      local_odd[i].where = myid;
    }

  if(myid % 2)
    {				/* odd # PE */
      odd = myid + 1;
      odd_inc = +1;
      odd_start = 0;
      odd_compare = LT;
      even = myid - 1;
      even_inc = -1;
      even_start = n - 1;
      even_compare = GT;
    }
  else
    {
      odd = myid - 1;
      odd_inc = -1;
      odd_start = n - 1;
      odd_compare = GT;
      even = myid + 1;
      even_inc = +1;
      even_start = 0;
      even_compare = LT;
    }
  if(odd < 0 || odd > numprocs - 1)
    odd = MPI_PROC_NULL;
  if(even > numprocs - 1)
    even = MPI_PROC_NULL;

  /*
   *  Find size of maximum data set residing in neighbor processor and
   *    allocate array to hold data during neighbor items during merge phases
   */

  MPI_Isend(&n, 1, MPI_INT, odd, 4, ADH_COMM, &done1);
  MPI_Isend(&n, 1, MPI_INT, even, 4, ADH_COMM, &done2);
  max1 = max2 = 0;
  MPI_Recv(&max1, 1, MPI_INT, odd, 4, ADH_COMM, &status);
  if(n > max1)
    max1 = n;
  MPI_Recv(&max2, 1, MPI_INT, even, 4, ADH_COMM, &status);
  if(max2 > max1)
    max1 = max2;
  MPI_Wait(&done1, &statusa);
  MPI_Wait(&done2, &statusa);
  max1 *= itemsize;		/* number of bytes outside */
  foreign = (struct sortitem *)tl_alloc(itemsize, max1);

  amount = n * itemsize;	/* number of bytes locally */

  /*
   *  Use exclusive prefix scan to determine number of items in preceding PEs
   *  MPI_Scan counts values of n for processors = 0-myid.  Then, subtract the
   *  n from the value to get the number of items in the preceding PEs
   */

  MPI_Scan(&n, &max2, 1, MPI_INT, MPI_SUM, ADH_COMM);
  max2 -= n;
  /*
   *  Odd-Even Transposition Sort
   */

  for(iteration = 0; iteration < (numprocs + 1) / 2; iteration++)
    {

      /*
       *  Odd Merge Phase
       */

      tag = iteration;
      MPI_Isend(local_odd, amount, MPI_BYTE, odd, tag, ADH_COMM, &done1);
      MPI_Recv(foreign, max1, MPI_BYTE, odd, tag, ADH_COMM, &status);
      MPI_Get_count(&status, MPI_BYTE, &count);
      MPI_Wait(&done1, &statusa);
      count /= itemsize;
      if(odd == MPI_PROC_NULL)
	count = 0;

      /*
       *  Set up local and foreign array pointers
       */

      odd_indx = p1 = p2 = odd_start;
      if(!(myid % 2))
	{
	  p2 = count - 1;
	  odd_indx = n - 1;
	}
      /*
       *   Main merge loop; only n iterations needed to replace local data
       */

      for(i = 0; i < n; i++)
	{

	  /*
	   *   Check to see if foreign values have been exhausted; simple copy
	   *   Else determine which value goes into merge results array
	   *   Update pointers
	   */

	  if(count == 0)
	    {
	      local_even[odd_indx] = local_odd[p1];
	      p1 += odd_inc;
	    }
	  else if(test(local_odd[p1].data, odd_compare, foreign[p2].data))
	    {
	      local_even[odd_indx] = local_odd[p1];
	      p1 += odd_inc;
	    }
	  else
	    {
	      local_even[odd_indx] = foreign[p2];
	      p2 += odd_inc;
	      count--;
	    }
	  odd_indx += odd_inc;

	}

      /*
       *  Even Merge Phase
       */
      tag *= numprocs;
      MPI_Isend(local_even, amount, MPI_BYTE, even, tag, ADH_COMM, &done1);
      MPI_Recv(foreign, max1, MPI_BYTE, even, tag, ADH_COMM, &status);
      MPI_Get_count(&status, MPI_BYTE, &count);
      MPI_Wait(&done1, &statusa);
      count /= itemsize;
      if(even == MPI_PROC_NULL)
	count = 0;

      /*
       *  Set up local and foreign array pointers
       */

      even_indx = p1 = p2 = even_start;
      if(myid % 2)
	{
	  p2 = count - 1;
	  even_indx = n - 1;
	}

      /*
       *   Main merge loop; only n iterations needed to replace local data
       */

      for(i = 0; i < n; i++)
	{

	  /*
	   *   Check to see if foreign values have been exhausted; simple copy
	   *   Else determine which value goes into merge results array
	   *   Update pointers
	   */

	  if(count == 0)
	    {
	      local_odd[even_indx] = local_even[p1];
	      p1 += even_inc;
	    }
	  else if(test(local_even[p1].data, even_compare, foreign[p2].data))
	    {
	      local_odd[even_indx] = local_even[p1];
	      p1 += even_inc;
	    }
	  else
	    {
	      local_odd[even_indx] = foreign[p2];
	      p2 += even_inc;
	      count--;
	    }
	  even_indx += even_inc;

	}
    }				/* end for iteration */
  foreign = (struct sortitem *)tl_free(itemsize, max1, foreign);
  
  /*
   *   Copy all local data to other array; change where field to item rank
   */

  for(i = 0; i < n; i++)
    {
      local_even[i].data = local_odd[i].data;
      local_even[i].where = max2++;
    }

  /* 
   *   Pull out rank of items that began in this PE
   */

  p1 = 0;
  for(i = 0; i < n; i++)
    {
      if(local_odd[i].where == myid)
	{
	  rank[p1++] = local_even[i].where;
	}
    }

  /*
   *  Find Maximum number of items held by any PE.
   */

  MPI_Allreduce(&n, &max2, 1, MPI_INT, MPI_MAX, ADH_COMM);
  max2++;

  /* 
   *   Pull out ranks of other PEs (in reverse order from myid), store
   *     in outrank using [0] element containing the count of elements in 
   *     array, send and receive outrank, place received values in rank
   */

  odd = myid;
  even = myid;
  outranks = (int *)tl_alloc(sizeof(int), max2);
  outrankr = (int *)tl_alloc(sizeof(int), max2);
  for(tag = 1; tag < numprocs; tag++)
    {
      odd = (odd == 0) ? numprocs - 1 : odd - 1;
      even = (even == numprocs - 1) ? 0 : even + 1;
      p2 = 1;
      for(i = 0; i < n; i++)
	{
	  if(local_odd[i].where == odd)
	    {
	      outranks[p2++] = local_even[i].where;
	    }
	}
      outranks[0] = p2 - 1;

      MPI_Isend(outranks, max2, MPI_INT, odd, tag, ADH_COMM, &done1);
      MPI_Recv(outrankr, max2, MPI_INT, even, tag, ADH_COMM, &status);
      MPI_Wait(&done1, &statusa);

      for(i = 1; i <= outrankr[0]; i++)
	{
	  rank[p1++] = outrankr[i];
	}

    }				/* end for tag */

  /* 
   *   Find location of ranks from lowest # PE 
   *   Rotate ranks to place in order
   */

  for(i = 1; i < n; i++)
    if(rank[i - 1] > rank[i])
      break;

  ROTATE_IT(rank, 0, n, i);

  local_odd = (struct sortitem *)tl_free(itemsize, n, local_odd);
  local_even = (struct sortitem *)tl_free(itemsize, n, local_even);
  outranks = (int *)tl_free(sizeof(int), max2, outranks);
  outrankr = (int *)tl_free(sizeof(int), max2, outrankr);

}				/* end findrank */

int test(
  double a,
  int comp,
  double b
)
{
  if(comp == LT)
    return (a < b);
  else
    return (a > b);
}

#define SWAP(x,y) {int t; t = x; x = y; y = t;}

void ROTATE_IT(
  int *a,
  int i,
  int h,
  int l
)
{
  int n, j;

  for(n = i, j = i + h - 1; n < j; n++, j--)
    SWAP(a[n], a[j]);
  for(n = i, j = i + h - l - 1; n < j; n++, j--)
    SWAP(a[n], a[j]);
  for(n = i + h - l, j = i + h - 1; n < j; n++, j--)
    SWAP(a[n], a[j]);
}				/* end ROTATE */
#else
void findrank()
  {
    return;
  }
#endif

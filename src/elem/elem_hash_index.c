/* returns index into element hash table - elements are 
   hashed based on their first node */

#include "global_header.h"

int elem_hash_index(
  int nd1			/* the first node */
)
{
  return nd1 % HASHSIZE;
}

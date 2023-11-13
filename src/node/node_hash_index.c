/* returns index into node hash table */

#include "global_header.h"

#ifdef _MESSG
int node_hash_index(SNODE global,  /* the global node number */
                    int npes
  )
{
  int key;
  int key1;
  int key2;

  key1 = (global.resident_id) %HASHSIZE;
  key2 = npes % HASHSIZE;
  key = (key1 * key2 + global.resident_pe) %HASHSIZE;
  return key;
}
#else
void node_hash_index(void)
{
}
#endif

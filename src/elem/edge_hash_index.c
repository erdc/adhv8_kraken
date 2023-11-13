/* returns index into edge hash table */

#include "global_header.h"

int edge_hash_index(
  int nd1,
  int nd2
)
{
  int key;
  int key1;
  int key2;

  key1 = nd1 % HASHSIZE;
  key2 = INT_MAX % HASHSIZE;
  key = (key1 * key2 + nd2) % HASHSIZE;
  return key;
}

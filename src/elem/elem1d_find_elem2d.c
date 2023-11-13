/* returns pointer to structure containing the 2D element */
/* We have to run the elem_hash_index on every node in the 1D element */
/* since we do not know the first node for the 2D element.  This determines */
/* the hash index, i.e. how the element is stored. */
/* Only two lookups at most are required since the hash table is setup using */
/* the two smallest node numbers in the 2d elements (one of which might not be */
/* in the 1d element that shares the element face, but one of the two is. */
/* We do the look up using the smallest node in the 1d element. */

#include "global_header.h"

ELEM2D_LIST_ITEM *elem1d_find_elem2d(
                                     int nd1,			/* the first node */
                                     int nd2,			/* the second node */
                                     ELEM2D_LIST_ITEM ** elem2d_hashtab	/* the hash table */
)
{
    int elem2d_hashval;		/* the index for the element */
    ELEM2D_LIST_ITEM *ep;   /* the entry in the hash table */
    
    elem2d_hashval = elem_hash_index(nd1); //this goes to the function used to breakup the elements into hashs
    
    for(ep = elem2d_hashtab[elem2d_hashval]; ep != NULL; ep = ep->next) {
        if((nd1 == ep->nd1) && ((nd2 == ep->nd2) || (nd2 == ep->nd3)))
            return ep; /* found */
    }
    return (NULL);
}

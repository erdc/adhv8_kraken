/* ADH Version 2.0.0 6-04 */
/* Takes a 3D element and finds the coordinates of the centroid of the 2D projection of the
 * a tetrahedron which has one vertical edge */

#include "global_header.h"

static int DEBUG = OFF;

/* When determining if we are at the correct location, DIFF_TOL is used
 * to allow some slight fudge in location */
#define DIFF_TOL .0000001

/* SHIFT is a multiplicative factor used in hash functions to convert
 * float values to integer values */
/* #define SHIFT 100000 */
#define SHIFT 1

/* the number of bins in each direction for the hash function */
#define nxBox 1000
#define nyBox 1000

#define QQ 7
#define ALPHA 1
#define BETA 1

/****************************************************************/
/****************************************************************/

void build_columns(SGRID *grid, int first_pass) {
    
    
#ifdef _DEBUG
    if (DEBUG) {
#ifdef _MESSG
      tag(grid->smpi->ADH_COMM);
#else
      tag();
#endif
    }
#endif
    
    /* Set up hash function and initialize */
    if (first_pass) {
        tl_column_init(grid);
#ifdef _DEBUG
    if (DEBUG) {
#ifdef _MESSG
      tag(grid->smpi->ADH_COMM);
#else
      tag();
#endif
    }
#endif
    }
    
    /* Classify all faces according the face file input */
    //classify_2d_elements(&grid); //tag();
    
    /* Build the columns of prismatic elements */
    build_column_hash(&grid);
#ifdef _DEBUG
    if (DEBUG) {
#ifdef _MESSG
      tag(grid->smpi->ADH_COMM);
#else
      tag();
#endif
    }
#endif
    build_column_list(&grid);
#ifdef _DEBUG
    if (DEBUG) {
#ifdef _MESSG
      tag(grid->smpi->ADH_COMM);
#else
      tag();
#endif
    }
#endif
    
    /* identify the 2D and 3D surface elements */
    classify_2d_elements(&grid);
#ifdef _DEBUG
    if (DEBUG) {
#ifdef _MESSG
      tag(grid->smpi->ADH_COMM);
#else
      tag();
#endif
    }
#endif
    //tl_find_elem3d_sur(&grid);
    
    /* build the hash and lists of vertical lines of nodes */
    build_vertical_hash(&grid);
#ifdef _DEBUG
    if (DEBUG) {
#ifdef _MESSG
      tag(grid->smpi->ADH_COMM);
#else
      tag();
#endif
    }
#endif
    build_vertical_list(&grid);
#ifdef _DEBUG
    if (DEBUG) {
#ifdef _MESSG
      tag(grid->smpi->ADH_COMM);
#else
      tag();
#endif
    }
#endif
    
    /* identify the bedom face of each column */
    //find_elems2d_bed(&grid);
    
    /* build the hash and lists of midpts (edges) */
    build_midpt_hash(&grid);
#ifdef _DEBUG
    if (DEBUG) {
#ifdef _MESSG
      tag(grid->smpi->ADH_COMM);
#else
      tag();
#endif
    }
#endif
    build_midpt_list(&grid);
#ifdef _DEBUG
    if (DEBUG) {
#ifdef _MESSG
      tag(grid->smpi->ADH_COMM);
#else
      tag();
#endif
    }
#endif
    
    // only call to first allocate ...
    //find_bedom_nodes(&grid);
    
    /* construct list of 2D sidewall elements */
    build_sidewall_list(&grid);
    /* fill in the elements that lie above/below each midpt (edge) */
    /* this isn't be used as of now, and I'm about to break it anyway */
    
    //  build_midpt_connectivity();
    
#ifdef _DEBUG
    if (DEBUG) {
#ifdef _MESSG
      tag(grid->smpi->ADH_COMM);
#else
      tag();
#endif
    }
#endif
}

/****************************************************************/
/****************************************************************/

#if 0
int tl_column_init() {
    
    /* compute the size of the bounding box */
    xL = x_max - x_min;
    yL = y_max - y_min;
    
    x_factor = ((double) nxBox) / xL;
    y_factor = ((double) nyBox) / yL;
    
    grid->hash_size = (nxBox + 1) * (nyBox + 1);
    
    return 1;
    
}
#endif

int tl_column_init(SGRID *grid) {
#define MAX_SIZE 100000
    int ii;
    double temp;
    int q = 7;
    int max_range = 10;
    
    grid->shift_factor = 1;
    grid->hash_size = 65536;
    
#if 1
    /* Here, I build a scaling factor which will be a power of 2.
     * This seems to work okay, but I'm a bit worried about the case where
     * the cooridinate values are large ( > 10^^5) in magnitude.
     * So, this will do until I work up a hack based on the mantissa
     */
    temp = MAX(fabs(grid->x_min), fabs(grid->x_max));
    if (temp > MAX_SIZE) {
        return 0;
    }
    else {
        while (temp < MAX_SIZE) {
            temp = temp * 2.0;
            grid->shift_factor = grid->shift_factor * 2;
        }
    }
    return 0;
#endif
    
    grid->scale_factor = 1;
    
    /* compute the limiting size for number of decimal digits (q) */
    for (ii = 1; ii < q; ii++) {
        max_range *= 10;
    }
    
    temp = MAX(fabs(grid->x_min), fabs(grid->x_max));
    if (temp > max_range) {
        return 0;
    }
    else {
        while (temp < max_range) {
            temp = 10 * temp;
            grid->scale_factor *= 10;
        }
    }
    return 0;
}

/****************************************************************/
/****************************************************************/

void build_column_hash(SGRID **grid3d) {
    
    int i;
    int ie;
    int ncolumns;
    
    SGRID *grid = (*grid3d); // alias
    
    
    free_column_hash(grid, grid->column_hash);
    free_column_hash2d(grid, grid->column_hash2d);
    
    if (grid->column_hash == NULL) {
        grid->column_hash = (CENT_LIST_ITEM **) tl_alloc(sizeof(CENT_LIST_ITEM *), grid->hash_size);
    }
    if (grid->column_hash2d == NULL) {
        grid->column_hash2d = (CENT_LIST_ITEM **) tl_alloc(sizeof(CENT_LIST_ITEM *), grid->hash_size);
    }
    
    for (i = 0; i < grid->hash_size; i++) {
        grid->column_hash[i] = (CENT_LIST_ITEM *) tl_alloc(sizeof(CENT_LIST_ITEM), ONE);
        grid->column_hash[i]->next = NULL;
        
        grid->column_hash2d[i] = (CENT_LIST_ITEM *) tl_alloc(sizeof(CENT_LIST_ITEM), ONE);
        grid->column_hash2d[i]->next = NULL;
    }
    grid->is_allocated_column_hash = 1;
    grid->is_allocated_column_hash2d = 1;
    
    
    grid->ncolumns = 0;
    for (ie = 0; ie < grid->nelems3d; ie++) {
        Add_column_hash_entry(grid, grid->column_hash, ie);
    }
    ncolumns = grid->ncolumns;
    
    grid->ncolumns = 0;
    for (ie = 0; ie < grid->nelems2d; ie++) {
        if (grid->elem2d[ie].bflag == 0) {
            Add_column_hash2d_entry(grid, grid->column_hash2d, ie);
        }
    }
    
    if (grid->ncolumns != ncolumns) {
        int nsurf=0, nbed=0, nside=0;
        for (i=0; i<grid->nelems2d; i++) {
            if (grid->elem2d[i].bflag == 0) nsurf++;
            if (grid->elem2d[i].bflag == 1) nbed++;
            if (grid->elem2d[i].bflag == 2) nside++;
        }
        printf("2d element counts :: nsurf: %d \t nbed: %d \t nside: %d\n",nsurf,nbed,nside);
        
        printf("NUMBER 2D COLUMNS = %d\n", grid->ncolumns);
        printf("NUMBER 3D COLUMNS = %d\n", ncolumns);
        tl_error("Number of 2D columns and 3D columns does not agree.\n");
    }
    
    /*
     if (grid->nelems2d_sur != grid->ncolumns) {
     printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
     tl_error(">> the number of surface elements should be equal to the number of columns\n");
     }
     if (grid->nelems2d_bed != grid->ncolumns) {
     printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
     tl_error(">> the number of bedom elements should be equal to the number of columns\n");
     }
     */
    
    /* DEBUG Stuff */
    /* hash_histogram(grid, grid->column_hash); */
}

/****************************************************************/
/****************************************************************/

void Add_column_hash_entry(SGRID *grid, CENT_LIST_ITEM **head, int ie)
{
    CENT_LIST_ITEM *ptr;
    SVECT2D center, entry;
    int hash_val;
    
    center = tl_calc_3dcentroid(grid, ie);
    hash_val = column_hash_function(grid, center);
    
    /* first check to see if there is no entry in this bin */
    if (head[hash_val] == NULL) {
        /* Null pointer means that we haven't filled this bin yet. Fire away. */
        push_column_hash_entry(grid, &head[hash_val], center);
        grid->ncolumns++;
        return;
    }
    
    ptr = head[hash_val];
    while (ptr->next != NULL) {
        entry = ptr->vect;
        /* check to see if this center is the same as a pre-existing entry */
        if ((fabs(entry.x - center.x) < DIFF_TOL) && (fabs(entry.y - center.y) < DIFF_TOL)) {
            return;   /* this point has already been hashed */
        }
        ptr = ptr->next;
    }
    
    /* Found a new entry; add it to the hashtable */
    push_column_hash_entry(grid, &head[hash_val], center);
    grid->ncolumns++;
    return;
}

/****************************************************************/
/****************************************************************/

void Add_column_hash2d_entry(SGRID *grid, CENT_LIST_ITEM ** head, int ie)
{
    CENT_LIST_ITEM *ptr;
    SVECT2D center, entry;
    int hash_val;
    
    center = tl_calc_2dcentroid(grid, ie);
    hash_val = column_hash_function(grid, center);
    
    /* first check to see if there is no entry in this bin */
    if (head[hash_val] == NULL) {
        /* Null pointer means that we haven't filled this bin yet. Fire away. */
        push_column_hash_entry(grid, &head[hash_val], center);
        grid->ncolumns++;
        return;
    }
    
    ptr = head[hash_val];
    while (ptr->next != NULL) {
        entry = ptr->vect;
        /* check to see if this center is the same as a pre-existing entry */
        if ((fabs(entry.x - center.x) < DIFF_TOL) && (fabs(entry.y - center.y) < DIFF_TOL)) {
            return;   /* this point has already been hashed */
        }
        ptr = ptr->next;
    }
    /* Found a new entry; add it to the hashtable */
    push_column_hash_entry(grid, &head[hash_val], center);
    grid->ncolumns++;
    return;
}

/****************************************************************/
/****************************************************************/

void hash_histogram(SGRID *grid, CENT_LIST_ITEM ** hash)
{
    FILE *fp;
    int *bins;
    int ii;
    int num_bins = 10;
    int length;
    
    fp = io_fopen("HISTOGRAM", "w", TRUE);
    
    bins = tl_alloc(sizeof(int), num_bins + 1);
    for (ii = 0; ii < num_bins + 1; ii++) {
        bins[ii] = 0;
    }
    
    for (ii = 0; ii < grid->hash_size; ii++) {
        /* get the length of the linked list */
        length = length_column_list(hash[ii]);
        if (length < 0) {
            tl_error("Error in length of list\n");
        }
        if (length > num_bins - 1) {
            bins[num_bins]++;
        }
        else {
            bins[length]++;
        }
    }
    
    fprintf(fp, "Length    No. of bins\n\n");
    for (ii = 0; ii <= num_bins; ii++) {
        fprintf(fp, "%d    %d\n", ii, bins[ii]);
    }
    fclose(fp);
}

/****************************************************************/
/****************************************************************/


int pre_intcompare(const void *p1, const void *p2)
{
    int i = *((int *) p1);
    int j = *((int *) p2);
    
    if (i > j) {
        return (1);
    }
    if (i < j) {
        return (-1);
    }
    return (0);
}

/****************************************************************/
/****************************************************************/

int pre_dblcompare(const void *p1, const void *p2)
{
    double i = *((int *) p1);
    double j = *((int *) p2);
    
    if (i > j) {
        return (1);
    }
    if (i < j) {
        return (-1);
    }
    return (0);
}

/****************************************************************/
/****************************************************************/

SVECT2D tl_calc_3dcentroid(SGRID *grid, int iel)
{
    
    int i, j;
    int flag;
    int count;
    int nnodes = grid->elem3d[iel].nnodes;
    double xval[3], yval[3];
    SVECT2D data[NDONPRISM];
    SVECT2D vect;
    
    for (i = 0; i < nnodes; i++) {
        data[i].x = grid->node[ grid->elem3d[iel].nodes[i] ].x;
        data[i].y = grid->node[ grid->elem3d[iel].nodes[i] ].y;
    }
    
    if (nnodes == NDONTET){
        count = 0;
        xval[count] = data[0].x;
        yval[count] = data[0].y;
        count++;
        
        for (i = 1; i < 4; i++) {
            flag = 1;
            /* check to see that this point doesn't lie beneath any of the previous ones */
            for (j = i - 1; j >= 0; j--) {
                if ((fabs(data[i].x - data[j].x) < DIFF_TOL) && (fabs(data[i].y - data[j].y) < DIFF_TOL)) {
                    flag = 0;
                    break;
                }
            }
            if (flag > 0) {
                if (count > 2) {
                    printf("Error: two points in the 3D element should lie in same vertical line\n");
                }
                xval[count] = data[i].x;
                yval[count] = data[i].y;
                count++;
            }
        }
    }
    else if (nnodes == NDONPRISM){ /* Assumption: Prism numbering follows nodes on triangular_face_1 followed by those on triangular_face_2 */
        for (i = 0; i < NDONTRI; i++) {
            xval[i] = data[i].x;
            yval[i] = data[i].y;
        }
    }
    
    /* Force the order of operations to be the same ... */
    qsort(xval, 3, sizeof(double), pre_dblcompare);
    qsort(yval, 3, sizeof(double), pre_dblcompare);
    
    vect.x = 0.0;
    vect.y = 0.0;
    for (i = 0; i < 3; i++) {
        vect.x += xval[i];
        vect.y += yval[i];
    }
    vect.x = vect.x / 3.0;
    vect.y = vect.y / 3.0;
    
    return (vect);
}

/****************************************************************/
/****************************************************************/

SVECT2D tl_calc_2dcentroid(SGRID *grid, int iel) {
    
    int i;
    int nnodes = grid->elem2d[iel].nnodes;
    double xval[NDONQUAD], yval[NDONQUAD]; /* over-allocating */
    SVECT2D data[NDONQUAD];
    SVECT2D vect;
    
    for (i = 0; i < nnodes; i++) {
        data[i].x = grid->node[ grid->elem2d[iel].nodes[i] ].x;
        data[i].y = grid->node[ grid->elem2d[iel].nodes[i] ].y;
    }
    
    for (i = 0; i < nnodes; i++)
    {
        xval[i] = data[i].x;
        yval[i] = data[i].y;
    }
    
    /* Force the order of operations to be the same ... */
    qsort(xval, nnodes, sizeof(double), pre_dblcompare);
    qsort(yval, nnodes, sizeof(double), pre_dblcompare);
    
    vect.x = 0;
    vect.y = 0;
    
    for (i = 0; i < nnodes; i++) {
        vect.x += xval[i];
        vect.y += yval[i];
    }
    vect.x = vect.x / (double) nnodes;
    vect.y = vect.y / (double) nnodes;
    
    return (vect);
}

/****************************************************************/
/****************************************************************/

int column_hash_function(SGRID *grid, SVECT2D center)
{
    int val;
    float temp1;
    float temp2;
    
    
    /*  temp1 = fabs(grid->scale_factor * center.x);
     temp2 = fabs(grid->scale_factor * center.y);
     val = ((int) (ALPHA * ((int) (temp1)) + BETA * ((int) (temp2)))) % grid->hash_size;
     */
    temp1 = fabs(grid->shift_factor * center.x);
    temp2 = fabs(grid->shift_factor * center.y);
    val = ((int) (ALPHA * ((int) (temp1)) + BETA * ((int) (temp2)))) % grid->hash_size;
    
    return (val % grid->hash_size);
}

/****************************************************************/
/****************************************************************/

void build_column_list(SGRID **grid3d)
{
    int i;
    int ie;
    int ival;
    int foo = 0;
    SGRID *grid = (*grid3d); // alias
    
    free_id_list(grid->column_list, grid->isize_ncolumns);
    if (grid->column_list != NULL) {
        tl_free(sizeof(ID_LIST_ITEM *), grid->isize_ncolumns, (void *) grid->column_list);
    }
    free_id_list(grid->column_list2d, grid->isize_ncolumns);
    if (grid->column_list2d != NULL) {
        tl_free(sizeof(ID_LIST_ITEM *), grid->isize_ncolumns, (void *) grid->column_list2d);
    }
    
    grid->isize_ncolumns = grid->ncolumns;
    
    grid->column_list = tl_alloc(sizeof(ID_LIST_ITEM *), grid->ncolumns);
    grid->column_list2d = tl_alloc(sizeof(ID_LIST_ITEM *), grid->ncolumns);
    for (i = 0; i < grid->ncolumns; i++) {
        grid->column_list[i] = tl_alloc(sizeof(ID_LIST_ITEM), ONE);
        grid->column_list[i]->next = NULL;
    }
    
    for (i = 0; i < grid->ncolumns; i++) {
        grid->column_list2d[i] = tl_alloc(sizeof(ID_LIST_ITEM), ONE);
        grid->column_list2d[i]->next = NULL;
    }
    for (ie = 0; ie < grid->nelems3d; ie++) {
        ival = find_column_from_centroid(grid, ie, grid->column_hash);
        Add_column_entry(&(grid->column_list[ival]), ie);
    }
    
    if (foo == 1) {
        int v1, v2;
        v1 = find_2d_hashvalue(grid, 0);
        v2 = find_3d_hashvalue(grid, 1);
    }
    
    for (ie = 0; ie < grid->nelems2d; ie++) {
        ival = find_column_from_center(grid, ie, grid->column_hash);
        if (ival != -1) {
            Add_column_entry(&(grid->column_list2d[ival]), ie);
        }
    }
}

/****************************************************************/
/****************************************************************/

int find_column_from_centroid(SGRID *grid, int ielem, CENT_LIST_ITEM ** list)
{
    SVECT2D center, entry;
    CENT_LIST_ITEM *head;
    int hash_val;
    
    center = tl_calc_3dcentroid(grid, ielem);
    hash_val = column_hash_function(grid, center);
    
    head = list[hash_val];
    while (head->next != NULL) {
        entry = head->vect;
        if ((fabs(entry.x - center.x) < DIFF_TOL) && (fabs(entry.y - center.y) < DIFF_TOL)) {
            /* found the match */
            return head->index;
        }
        head = head->next;
    }
    tl_error("Did not find a column for this element.\n");
    return -1;
}

int find_column_from_center(SGRID *grid, int ielem, CENT_LIST_ITEM ** column_hash)
{
    SVECT2D center, entry;
    CENT_LIST_ITEM *head;
    int hash_val;
    center = tl_calc_2dcentroid(grid, ielem);
    hash_val = column_hash_function(grid, center);
    
    head = column_hash[hash_val];
    while (head->next != NULL) {
        entry = head->vect;
        if ((fabs(entry.x - center.x) < DIFF_TOL) && (fabs(entry.y - center.y) < DIFF_TOL)) {
            /* found the match */
            return head->index;
        }
        head = head->next;
    }
    return -1;
}

int find_column_from_point(SGRID *grid, SVECT2D center, CENT_LIST_ITEM ** list)
{
    SVECT2D entry;
    CENT_LIST_ITEM *head;
    int hash_val;
    
    hash_val = column_hash_function(grid, center);
    
    head = list[hash_val];
    while (head->next != NULL) {
        entry = head->vect;
        if ((fabs(entry.x - center.x) < DIFF_TOL) && (fabs(entry.y - center.y) < DIFF_TOL)) {
            /* found the match */
            return head->index;
        }
        head = head->next;
    }
    tl_error("Did not find a column for this element.\n");
    return -1;
}

/****************************************************************/
/****************************************************************/

void Add_column_entry(ID_LIST_ITEM ** head, int id)
{
    ID_LIST_ITEM *ptr;
    
    ptr = tl_alloc(sizeof(ID_LIST_ITEM), ONE);
    ptr->id = id;
    ptr->next = *head;
    *head = ptr;
    return;
}

/****************************************************************/
/****************************************************************/

int find_2d_hashvalue(SGRID *grid, int ie)
{
    int ival;
    
    ival = find_column_from_center(grid, ie, grid->column_hash);
    return ival;
}

/****************************************************************/
/****************************************************************/

int find_3d_hashvalue(SGRID *grid, int ie)
{
    
    int ival;
    
    ival = find_column_from_centroid(grid, ie, grid->column_hash);
    return ival;
}

/****************************************************************/
/****************************************************************/

void classify_2d_elements(SGRID **grid3d) {

    
    //tl_check_all_pickets(__FILE__, __LINE__);
    
    int i, ie, inode, ival;
    SGRID *grid = (*grid3d); // alias
    
    //assert(grid->type == COLUMNAR);
    
    grid->nelems2d_sur = 0;
    grid->nelems2d_bed = 0;
    grid->nelems2d_sidewall = 0;
    for (ie = 0; ie < grid->nelems2d; ie++) {
        if      (grid->elem2d[ie].bflag == 0) (grid->nelems2d_sur)++;
        else if (grid->elem2d[ie].bflag == 1) (grid->nelems2d_bed)++;
        else if (grid->elem2d[ie].bflag == 2) (grid->nelems2d_sidewall)++;
        else {
            printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
            tl_error(">> some 2d element does not have valid boundary flag\n");
        }
    }
    
    /* sanity checks */
    if ( (grid->nelems2d - grid->nelems2d_sur - grid->nelems2d_bed - grid->nelems2d_sidewall) != 0) {
        printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
        tl_error(">> some 2d elements do not have correct boundary flags\n");
    }
    
    /* deallocate previous arrays */
    if (grid->isize_sur_elems > 0) {
        tl_free(sizeof(int), grid->isize_sur_elems, (void *) grid->elem3d_sur);
        tl_free(sizeof(int), grid->isize_sur_elems, (void *) grid->elem2d_sur);
    }
    if (grid->isize_bed_elems > 0) {
        tl_free(sizeof(int), grid->isize_bed_elems, (void *) grid->elem3d_bed);
        tl_free(sizeof(int), grid->isize_bed_elems, (void *) grid->elem2d_bed);
    }
    if (grid->isize_sidewall_elems > 0) {
        tl_free(sizeof(int), grid->isize_sidewall_elems, (void *) grid->elem3d_sidewall);
        tl_free(sizeof(int), grid->isize_sidewall_elems, (void *) grid->elem2d_sidewall);
    }
    
    if(grid->nelems2d_sur>0)grid->elem3d_sur = tl_alloc(sizeof(int), grid->nelems2d_sur);
    if(grid->nelems2d_sur>0)grid->elem2d_sur = tl_alloc(sizeof(int), grid->nelems2d_sur);
    if(grid->nelems2d_bed>0)grid->elem3d_bed = tl_alloc(sizeof(int), grid->nelems2d_bed);
    if(grid->nelems2d_bed>0)grid->elem2d_bed = tl_alloc(sizeof(int), grid->nelems2d_bed);
    if(grid->nelems2d_sidewall>0)grid->elem3d_sidewall = tl_alloc(sizeof(int), grid->nelems2d_sidewall);
    if(grid->nelems2d_sidewall>0)grid->elem2d_sidewall = tl_alloc(sizeof(int), grid->nelems2d_sidewall);
    
    grid->isize_sur_elems = grid->nelems2d_sur;
    grid->isize_bed_elems = grid->nelems2d_bed;
    grid->isize_sidewall_elems = grid->nelems2d_sidewall;
    
    for (i = 0; i < grid->nelems2d_sur; i++) {
        grid->elem3d_sur[i] = UNSET_INT;   // 3d surface element id belonging to column i
        grid->elem2d_sur[i] = UNSET_INT;   // 2d surface element id belonging to column i
	/*mwf looks like these should be looped separately from surface? using nelems2d_bed*/
        /* grid->elem3d_bed[i] = UNSET_INT;   // 3d bedom element id belonging to column i */
        /* grid->elem2d_bed[i] = UNSET_INT;   // 2d bedom element id belonging to column i */
    }
    for (i = 0; i < grid->nelems2d_bed; i++) {
         grid->elem3d_bed[i] = UNSET_INT;   // 3d bedom element id belonging to column i
         grid->elem2d_bed[i] = UNSET_INT;   // 2d bedom element id belonging to column i
    }
    for (i = 0; i < grid->nelems2d_sidewall; i++) {
        grid->elem3d_sidewall[i] = UNSET_INT;   // 3d sidewall element belonging to column i
        grid->elem2d_sidewall[i] = UNSET_INT;   // 2d sidewall element belonging to column i
    }
    
    /* Match the free surface elements to columns. Also, mark nodes that are on the free surface. */
    int *sur_nodes = tl_alloc(sizeof(int), grid->nnodes);
    int *bed_nodes = tl_alloc(sizeof(int), grid->nnodes);
    for (inode = 0; inode < grid->nnodes; inode++) {
        sur_nodes[inode] = UNSET_INT;
        bed_nodes[inode] = UNSET_INT;
    }
    
    int elem2d_sur_count = 0, elem2d_bed_count = 0; /* Gajanan gkc May 2020. Only for non-columnar grids */
    for (ie = 0; ie < grid->nelems2d; ie++) {
        if (grid->elem2d[ie].bflag == 0) {
            // surface element (also columns)
            for (inode = 0; inode < NDPRFC; inode++) {
#ifdef _DEBUG
                if (grid->elem2d[ie].nodes[inode] >= grid->nnodes) {
                    printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
                    printf("grid->elem2d[ie].nodes[inode]: %d \t grid->nnodes: %d\n",grid->elem2d[ie].nodes[inode],grid->nnodes);
                    tl_error(">> There are no nodes on this 2d grid.");
                }
#endif
                sur_nodes[grid->elem2d[ie].nodes[inode]] = 1;
            }
            // below only makes sense for colunnar grids
            if (grid->type == COLUMNAR) {
                ival = find_column_from_center(grid, ie, grid->column_hash);
                if (ival != -1) {
                    grid->elem3d_sur[ival] = grid->elem2d[ie].id_3d;
                    grid->elem2d_sur[ival] = grid->elem2d[ie].id;
                }
            }
            else{ /* Gajanan gkc May 2020. Non-columnar grids, just march through serially to number elements */
                grid->elem3d_sur[elem2d_sur_count] = grid->elem2d[ie].id_3d;
                grid->elem2d_sur[elem2d_sur_count] = grid->elem2d[ie].id;
                elem2d_sur_count++;
            }
        }
        else if (grid->elem2d[ie].bflag == 1) {
            // bed elements
            for (inode = 0; inode < NDPRFC; inode++) {
#ifdef _DEBUG
                if (grid->elem2d[ie].nodes[inode] >= grid->nnodes) {
                    printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
                    printf("grid->elem2d[ie].nodes[inode]: %d \t grid->nnodes: %d\n",grid->elem2d[ie].nodes[inode],grid->nnodes);
                    tl_error(">> There are no nodes on this 2d grid.");
                }
#endif
                bed_nodes[grid->elem2d[ie].nodes[inode]] = 1;
            }
            // below only makes sense for colunnar grids
            if (grid->type == COLUMNAR) {
                ival = find_column_from_center(grid, ie, grid->column_hash);
                if (ival != -1) {
                    grid->elem3d_bed[ival] = grid->elem2d[ie].id_3d;
                    grid->elem2d_bed[ival] = grid->elem2d[ie].id;
                }
            }
            else{ /* Gajanan gkc May 2020. Non-columnar grids, just march through serially to number elements */
                grid->elem3d_bed[elem2d_bed_count] = grid->elem2d[ie].id_3d;
                grid->elem2d_bed[elem2d_bed_count] = grid->elem2d[ie].id;
                elem2d_bed_count++;
            }
        }
    }
    
    if (grid->type == COLUMNAR) {
        if (grid->nelems2d_sur != grid->ncolumns) {
            printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
            tl_error(">> the number of surface elements should be equal to the number of columns\n");
        }
        if (grid->nelems2d_bed != grid->ncolumns) {
            printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
            tl_error(">> the number of bedom elements should be equal to the number of columns\n");
        }
    }
    
    /* count the number of surface and bed nodes */
    grid->nnodes_sur = 0;
    grid->my_nnodes_sur = 0.;
    grid->nnodes_bed = 0;
    grid->my_nnodes_bed = 0.;
    for (inode = 0; inode < grid->nnodes; inode++) {
        if (sur_nodes[inode] == 1) {
            grid->nnodes_sur++;
            if(grid->node[inode].resident_pe == grid->smpi->myid) {
                grid->my_nnodes_sur++;
            }
        };
        if (bed_nodes[inode] == 1) {
            grid->nnodes_bed++;
            if(grid->node[inode].resident_pe == grid->smpi->myid) {
                grid->my_nnodes_bed++;
            }
        }
    }
    
    /* allocate maps */
    /*mwf separate out bed and surface logic? */
    //if (grid->type == COLUMNAR) {
        if(grid->max_nnodes_sur == 0){
            grid->nodeID_2d_to_3d_sur = (int *) tl_alloc(sizeof(int), grid->nnodes_sur);
            grid->nodeID_3d_to_2d_sur = (int *) tl_alloc(sizeof(int), grid->nnodes);
        }
        if(grid->max_nnodes_bed == 0){
            grid->nodeID_2d_to_3d_bed = (int *) tl_alloc(sizeof(int), grid->nnodes_bed);
            grid->nodeID_3d_to_2d_bed = (int *) tl_alloc(sizeof(int), grid->nnodes);
        }
        
        /* initialize maps */
        sarray_init_value_int(grid->nodeID_2d_to_3d_sur, grid->nnodes_sur, UNSET_INT);
        sarray_init_value_int(grid->nodeID_2d_to_3d_bed, grid->nnodes_bed, UNSET_INT);
        sarray_init_value_int(grid->nodeID_3d_to_2d_sur, grid->nnodes, UNSET_INT);
        sarray_init_value_int(grid->nodeID_3d_to_2d_bed, grid->nnodes, UNSET_INT);
        
        /* count the number of surface nodes */
        int nsur = 0;
        int nbed = 0;
        for (inode = 0; inode < grid->nnodes; inode++) {
            if (sur_nodes[inode] == 1) {
                grid->nodeID_3d_to_2d_sur[inode] = nsur;
                grid->nodeID_2d_to_3d_sur[nsur] = inode;
                nsur++;
            };
            if (bed_nodes[inode] == 1) {
                grid->nodeID_3d_to_2d_bed[inode] = nbed;
                grid->nodeID_2d_to_3d_bed[nbed] = inode;
                nbed++;
            }
        }
    //}
    
    if (sur_nodes != NULL) sur_nodes = (int*) tl_free(sizeof(int), grid->nnodes, sur_nodes);
    if (bed_nodes != NULL) bed_nodes = (int*) tl_free(sizeof(int), grid->nnodes, bed_nodes);
    
}
//
///* Look for the free surface element in each column of elements */
//void tl_find_elem3d_sur(SGRID **grid3d)
//{
//  int i, j;
//  int ii;
//  int ie;
//  int istr;
//  int test_count;
//  int nd[3];                    /* the surface nodes */
//  int test_nodes[3];
//  int found_elem2d = 0;
//  ID_LIST_ITEM *ptr;
//  int *sur_nodes;
//  int ival;
//
//  SGRID *grid = (*grid3d); // alias
//
//  if (grid->isize_sur_elems > 0) {
//    tl_free(sizeof(int), grid->isize_sur_elems, (void *) grid->elem3d_sur);
//    tl_free(sizeof(int), grid->isize_sur_elems, (void *) grid->elem2d_sur);
//  }
//
//  grid->elem3d_sur = tl_alloc(sizeof(int), grid->ncolumns);
//  grid->elem2d_sur = tl_alloc(sizeof(int), grid->ncolumns);
//  grid->isize_sur_elems = grid->ncolumns;
//
//  for (i = 0; i < grid->ncolumns; i++) {
//    grid->elem3d_sur[i] = UNSET_INT;
//    grid->elem2d_sur[i] = UNSET_INT;
//  }
//
//  /* Match the free surface elements to columns. Also, mark nodes that are on the free surface. */
//  sur_nodes = tl_alloc(sizeof(int), grid->nnodes);
//  for (ii = 0; ii < grid->nnodes; ii++) {
//    sur_nodes[ii] = UNSET_INT;
//  }
//
//  for (ie = 0; ie < grid->nelems2d; ie++) {
//    if (str_values[grid->elem2d[ie].flag].flow.bc_flag == BCT_FRS) {
//      for (ii = 0; ii < NDPRFC; ii++) {
//        sur_nodes[grid->elem2d[ie].nodes[ii]] = 1;
//      }
//      ival = find_column_from_center(grid, ie, grid->column_hash);
//      if (ival != -1) {
//        grid->elem2d_sur[ival] = ie;
//      }
//      else {
//        printf("Ruh roh! Did not find a column for this element!!!!!!\n\n\n");
//      }
//    }
//  }
//
//  /* count the number of surface nodes */
//  grid->nnodes_sur = 0;
//  for (ii = 0; ii < grid->nnodes; ii++) {
//    if (sur_nodes[ii] == 1) {
//      grid->nnodes_sur++;
//    }
//  }
//
//  for (i = 0; i < grid->ncolumns; i++) {
//    ptr = grid->column_list[i];
//    while (ptr->next != NULL) {
//      ie = ptr->id;
//      test_count = 0;
//
//      for (j = 0; j < NDPRELM; j++) {
//        if (sur_nodes[grid->elem3d[ie].nodes[j]] == 1) {
//          test_count++;
//        }
//      }
//
//      if (test_count > 2) {
//        /* We've found an element on the free surface */
//        grid->elem3d_sur[i] = ie;
//        break;
//      }
//      ptr = ptr->next;
//    }
//  }
//
//  /* cjt */
//  if (sur_nodes != NULL) sur_nodes = (int*) tl_free(sizeof(int), grid->nnodes, sur_nodes);
//}


/****************************************************************/
/****************************************************************/

/* Use the hash to build a vertical columns of nodes that are sorted
 * by elevation
 */

/* WARNING: need to check logic if the point to be added will come at the end of the
 * linked list
 */
void Add_sorted_vertical_entry(SGRID *grid, ID_LIST_ITEM ** head, int id)
{
    ID_LIST_ITEM *ptr, *entry;
    ID_LIST_ITEM *next = NULL;
    ID_LIST_ITEM *prev = NULL;
    double z, this_z;
    
    ptr = tl_alloc(sizeof(ID_LIST_ITEM), ONE);
    ptr->id = id;
    
    entry = *head;
    
    /* Check to see if the list is empty; if so, add this entry and return */
    if (entry->next == NULL) {
        ptr->next = *head;
        ptr->prev = ptr;
        *head = ptr;
        return;
    }
    
    z = grid->node[entry->id].z;
    this_z = grid->node[id].z;
    
    /* Entries are sorted by z-depth. The first entry in the list will be the node
     * on the surface (with greatest z value)
     */
    
    /* If the new entry lies above the first in the linked list, then the new entry
     * becomes the head of the linked list
     */
    if (this_z > z) {
        ptr->next = entry;
        ptr->prev = entry->prev;
        *head = ptr;
        entry->prev = ptr;
        return;
    }
    
    /* Follow the linked list until we find location to insert this new entry  */
    while (entry->next != NULL) {
        next = entry->next;
        z = grid->node[entry->id].z;
        if (this_z > z)
        /* Append entries with greater z-depth , this way the linked
         * list starts at the top, and we go down the column as
         * we follow the links
         */
        {
            if (prev == NULL) {
                *head = ptr;
                prev = *head;
            }
            else {
                prev->next = ptr;
            }
            ptr->next = entry;
            ptr->prev = prev;
            entry->prev = ptr;
            return;
        }
        prev = entry;
        entry = entry->next;
    }
    
    /* If we reach this point, then the entry to be added is at the end of the linked list */
    prev->next = ptr;
    ptr->next = entry;
    ptr->prev = prev;
    (*head)->prev = ptr;
    
    return;
}

/****************************************************************/
/****************************************************************/

int length_node_list(ID_LIST_ITEM * head)
{
    int count = 0;
    
    while (head->next != NULL) {
        count++;
        head = head->next;
    }
    return count;
}


/****************************************************************/
/****************************************************************/

void Add_vertical_node_hash_entry(SGRID *grid, CENT_LIST_ITEM ** head, int inode)
{
    CENT_LIST_ITEM *ptr;
    SVECT2D center, entry;
    int hash_val;
    
    center.x = grid->node[inode].x;
    center.y = grid->node[inode].y;
    hash_val = column_hash_function(grid, center);
    
    ptr = head[hash_val];
    while (ptr->next != NULL) {
        entry = ptr->vect;
        /* check to see if this center is the same as a pre-existing entry */
        if ((fabs(entry.x - center.x) < DIFF_TOL) && (fabs(entry.y - center.y) < DIFF_TOL)) {
            return;   /* this point has already been hashed */
        }
        ptr = ptr->next;
    }
    /* Found a new entry; add it to the hashtable */
    push_vertical_hash_entry(grid, &head[hash_val], center);
    grid->num_vert_segments++;
    return;
}

/****************************************************************/
/****************************************************************/


int find_sidewall_column(SGRID *grid, int ielem, CENT_LIST_ITEM ** list)
{
    SVECT2D center, entry, data[3];
    CENT_LIST_ITEM *head;
    int hash_val;
    int i;
    
    for (i = 0; i < 3; i++) {
        data[i].x = grid->node[grid->elem2d[ielem].nodes[i]].x;
        data[i].y = grid->node[grid->elem2d[ielem].nodes[i]].y;
    }
    
    center.x = data[0].x;
    center.y = data[0].y;
    
    if ((fabs(data[1].x - center.x) < DIFF_TOL) && (fabs(data[1].y - center.y) < DIFF_TOL)) {
        center.x += data[1].x;
        center.y += data[1].y;
    }
    else {
        center.x += data[1].x;
        center.y += data[1].y;
    }
    
    center.x = center.x / 2.0;
    center.y = center.y / 2.0;
    
    hash_val = column_hash_function(grid, center);
    
    head = list[hash_val];
    while (head->next != NULL) {
        entry = head->vect;
        if ((fabs(entry.x - center.x) < DIFF_TOL) && (fabs(entry.y - center.y) < DIFF_TOL)) {
            /* found the match */
            return head->index;
        }
        head = head->next;
    }
    return -1;
}

/****************************************************************/
/****************************************************************/

/* Build linked lists of sidewall 2D elements
 * NOTE: The midpt hash must be built first
 */
void build_sidewall_list(SGRID **grid3d)
{
    int i;
    int ie;
    int ival;
    int nd1, nd2, nd3;
    SVECT midpt;
    
    SGRID *grid = (*grid3d); // alias
    
    free_id_list(grid->sidewall_list, grid->isize5);
    /* cjt */
    if (grid->sidewall_list != NULL) {
        tl_free(sizeof(ID_LIST_ITEM *), grid->isize5, (void *) grid->sidewall_list);
    }
    
    
    if(grid->num_midpts>0) grid->sidewall_list = tl_alloc(sizeof(ID_LIST_ITEM *), grid->num_midpts);
    grid->isize5 = grid->num_midpts;
    
    for (i = 0; i < grid->num_midpts; i++) {
        grid->sidewall_list[i] = tl_alloc(sizeof(ID_LIST_ITEM), ONE);
        grid->sidewall_list[i]->next = NULL;
    }
    for (ie = 0; ie < grid->nelems2d; ie++) {
        nd1 = grid->elem2d[ie].nodes[0];
        nd2 = grid->elem2d[ie].nodes[1];
        nd3 = grid->elem2d[ie].nodes[2];
        
        /* Check to see if a pair of nodes is on same vertical line */
        
        if ((fabs(grid->node[nd1].x - grid->node[nd2].x) < DIFF_TOL) && (fabs(grid->node[nd1].y - grid->node[nd2].y) < DIFF_TOL)) {
            midpt.x = (grid->node[nd1].x + grid->node[nd3].x) * 0.5;
            midpt.y = (grid->node[nd1].y + grid->node[nd3].y) * 0.5;
        }
        else if ((fabs(grid->node[nd1].x - grid->node[nd3].x) < DIFF_TOL) && (fabs(grid->node[nd1].y - grid->node[nd3].y) < DIFF_TOL)) {
            midpt.x = (grid->node[nd1].x + grid->node[nd2].x) * 0.5;
            midpt.y = (grid->node[nd1].y + grid->node[nd2].y) * 0.5;
        }
        else if ((fabs(grid->node[nd2].x - grid->node[nd3].x) < DIFF_TOL) && (fabs(grid->node[nd2].y - grid->node[nd3].y) < DIFF_TOL)) {
            midpt.x = (grid->node[nd1].x + grid->node[nd3].x) * 0.5;
            midpt.y = (grid->node[nd1].y + grid->node[nd3].y) * 0.5;
        }
        /* If we get to here, then this is not a sidewall element */
        else {
            continue;
        }
        
        ival = find_midpt_segment(grid, midpt, grid->midpt_hash);
        if (ival >= 0) {
            Add_sidewall_entry(&grid->sidewall_list[ival], ie);
        }
    }
}

/****************************************************************/
/****************************************************************/

void Add_sidewall_entry(ID_LIST_ITEM ** head, int id)
{
    ID_LIST_ITEM *ptr;
    
    
    ptr = tl_alloc(sizeof(ID_LIST_ITEM), ONE);
    ptr->id = id;
    
    ptr->next = *head;
    *head = ptr;
    return;
}

/****************************************************************/
/****************************************************************/

void push_vertical_hash_entry(SGRID *grid, CENT_LIST_ITEM ** head, SVECT2D center)
{
    CENT_LIST_ITEM *ptr;
    
    ptr = tl_alloc(sizeof(CENT_LIST_ITEM), ONE);
    ptr->vect = center;
    ptr->index = grid->num_vert_segments;
    ptr->next = *head;
    *head = ptr;
    
    return;
}

/****************************************************************/
/****************************************************************/

void push_column_hash_entry(SGRID *grid, CENT_LIST_ITEM ** head, SVECT2D center)
{
    CENT_LIST_ITEM *ptr;
    
    ptr = tl_alloc(sizeof(CENT_LIST_ITEM), ONE);
    ptr->vect = center;
    ptr->index = grid->ncolumns;
    ptr->next = *head;
    *head = ptr;
    
    return;
    
}

/****************************************************************/
/****************************************************************/

void free_column_hash(SGRID *grid, CENT_LIST_ITEM ** list)
{
    CENT_LIST_ITEM *ptr, *next;
    int i;
    
    /* check to see if we need to free anything */
    if (grid->is_allocated_column_hash == 0) {
        return;
    }
    
    for (i = 0; i < grid->hash_size; i++) {
        ptr = list[i];
        while (ptr->next != NULL) {
            next = ptr->next;
            tl_free(sizeof(CENT_LIST_ITEM), ONE, (void *) ptr);
            ptr = next;
        }
        tl_free(sizeof(CENT_LIST_ITEM), ONE, (void *) ptr);
    }
}

/****************************************************************/
/****************************************************************/

void free_id_list(ID_LIST_ITEM ** list, int number  /* number of columns (the extent of this array) */
)
{
    ID_LIST_ITEM *ptr, *next;
    int i;
    
    for (i = 0; i < number; i++) {
        ptr = list[i];
        while (ptr->next != NULL) {
            next = ptr->next;
            tl_free(sizeof(ID_LIST_ITEM), ONE, (void *) ptr);
            ptr = next;
        }
        tl_free(sizeof(ID_LIST_ITEM), ONE, (void *) ptr);
    }
}

/****************************************************************/
/****************************************************************/

int find_vertical_segment(SGRID *grid, int inode, CENT_LIST_ITEM ** list)
{
    SVECT2D center, entry;
    CENT_LIST_ITEM *head;
    int hash_val;
    
    center.x = grid->node[inode].x;
    center.y = grid->node[inode].y;
    
    hash_val = column_hash_function(grid, center);
    
    head = list[hash_val];
    while (head->next != NULL) {
        entry = head->vect;
        if ((fabs(entry.x - center.x) < DIFF_TOL) && (fabs(entry.y - center.y) < DIFF_TOL)) {
            /* found the match */
            return head->index;
        }
        head = head->next;
    }
    /* tl_error("Did not find a vertical string for this node.\n"); */
    /* printf("DANGER WILL ROBINSON!\n"); */
    
    /* If we reach here it is not necessarily an error. This can happen if
     * the unrefinement has introduced a node owned by another processor, and
     * thus not currently in the hash.
     *
     * This is not a problem; when this occurs, we are trying to identify nodes
     * that lie beneath the surface node being removed. The return value indicates
     * failure, and the other routine ignores the node (which is correct).
     */
    return -1;
}

/****************************************************************/
/****************************************************************/

/* Look for the 2D element at the bedom of each column of elements */
/* Note: we can't rely on boundary conditions to identify the elements, so we
 * will find the three bedom nodes of the column and then match the 2D element
 */
//void find_elems2d_bed(SGRID **grid3d)
//{
//  int i, j, ii;
//  int ie;
//  int nd[3];
//  int ival;
//  int bot_nd[3];
//  int test_nodes[3];
//  int found_elem2d = 0;
//  ID_LIST_ITEM *ptr;
//  SVECT2D center;
//  double xval[3], yval[3];
//
//  SGRID *grid = (*grid3d); // alias
//
//  if (grid->isize_bed_elems > 0) {
//    tl_free(sizeof(int), grid->isize_bed_elems, (void *) grid->elems2d_bed);
//  }
//
//  grid->elems2d_bed = tl_alloc(sizeof(int), grid->ncolumns);
//  grid->isize_bed_elems = grid->ncolumns;
//
//  for (i = 0; i < grid->ncolumns; i++) {
//    found_elem2d = 0;
//    /* get the 3 surface nodes in the column and find the
//     * bedom nodes in the column */
//    ie = grid->elem2d_sur[i];
//    for (ii = 0; ii < 3; ii++) {
//      nd[ii] = grid->elem2d[ie].nodes[ii];
//      ival = find_vertical_segment(grid, nd[ii], grid->vertical_hash);
//      ptr = grid->vertical_list[ival];
//      bot_nd[ii] = ptr->prev->id;
//    }
//    /* now find the 2D element that has these nodes ... */
//    center.x = 0.0;
//    center.y = 0.0;
//
//    /* Force the order of operations to be the same ... */
//    for (j = 0; j < 3; j++)
//      {
//        xval[j] = grid->node[nd[j]].x;
//        yval[j] = grid->node[nd[j]].y;
//      }
//    qsort(xval, 3, sizeof(double), pre_dblcompare);
//    qsort(yval, 3, sizeof(double), pre_dblcompare);
//
//    center.x = 0.0;
//    center.y = 0.0;
//
//    for (ii = 0; ii < 3; ii++)
//      {
//        center.x += xval[ii];
//        center.y += yval[ii];
//      }
//
///*    for (j = 0; j < 3; j++) {
//      center.x += node[nd[j]].x;
//      center.y += node[nd[j]].y;
//    }
//*/
//
//    center.x = center.x / 3.0;
//    center.y = center.y / 3.0;
//    ival = find_column_from_point(grid, center, grid->column_hash);
//    qsort((void *) bot_nd, 3, sizeof(int), pre_intcompare);
//
//    /* loop over the 2D elements in the column looking for a match */
//    ptr = grid->column_list2d[ival];
//    while (ptr->next != NULL) {
//      ie = ptr->id;
//      test_nodes[0] = grid->elem2d[ie].nodes[0];
//      test_nodes[1] = grid->elem2d[ie].nodes[1];
//      test_nodes[2] = grid->elem2d[ie].nodes[2];
//      qsort((void *) test_nodes, 3, sizeof(int), pre_intcompare);
//      if ((test_nodes[0] == bot_nd[0])
//          && (test_nodes[1] == bot_nd[1])
//          && (test_nodes[2] == bot_nd[2])) {
//        /* we found a match */
//        found_elem2d = 1;
//        grid->elems2d_bed[i] = ie;
//        break;
//      }
//      ptr = ptr->next;
//    }
//    if (found_elem2d == 0) {
//      tl_error("Could not find the bedom surface element.\n");
//    }
//  }
//}
//
/****************************************************************/
/****************************************************************/

//void find_bedom_nodes(SGRID **grid3d)
//{
//  int i, j = 0;
//  ID_LIST_ITEM *ptr;
//
//  SGRID *grid = (*grid3d); // alias
//
//  if (grid->bedom_nodes != NULL) {
//	//if (myid ==3) printf("**** freeing bedom_node from node: %d\n",grid->isize_nnodes3d_prev);
//	tl_free(sizeof(int), grid->isize_nnodes3d_prev, (void *) grid->bedom_nodes);  // cjt
//  }
//
//  //if (myid==3) printf("**** allocating bedom_nodes to nnode: %d\n",nnode);
//  grid->isize_nnodes3d_prev = grid->nnodes;
//  grid->bedom_nodes = (int *) tl_alloc(sizeof(int), grid->isize_nnodes3d_prev);
//
//  for (i = 0; i < grid->nnodes; i++)
//    grid->bedom_nodes[i] = NO;
//
//  for (i = 0; i < grid->nnodes_sur; i++) {
//    ptr = grid->vertical_list[i];
//    grid->bedom_nodes[(ptr->prev)->id] = YES;
//  }
//}

/****************************************************************/
/****************************************************************/

void build_vertical_hash(SGRID **grid3d)
{
    int i;
    
    SGRID *grid = (*grid3d); // alias
    
    free_vertical_hash(grid, grid->vertical_hash);
    grid->num_vert_segments = 0;
    if (grid->vertical_hash == NULL) {
        grid->vertical_hash = (CENT_LIST_ITEM **) tl_alloc(sizeof(CENT_LIST_ITEM *), grid->hash_size);
    }
    for (i = 0; i < grid->hash_size; i++) {
        grid->vertical_hash[i] = (CENT_LIST_ITEM *) tl_alloc(sizeof(CENT_LIST_ITEM), ONE);
        grid->vertical_hash[i]->next = NULL;
    }
    for (i = 0; i < grid->nnodes; i++) {
        Add_vertical_node_hash_entry(grid, grid->vertical_hash, i);
    }
    
    grid->is_allocated_vertical_hash = 1;
    
}

/****************************************************************/
/****************************************************************/

void free_vertical_hash(SGRID *grid, CENT_LIST_ITEM ** list)
{
    CENT_LIST_ITEM *ptr, *next;
    int i;
    
    /* check to see if we need to free anything */
    if (grid->is_allocated_vertical_hash == 0) {
        return;
    }
    for (i = 0; i < grid->hash_size; i++) {
        ptr = list[i];
        while (ptr->next != NULL) {
            next = ptr->next;
            tl_free(sizeof(CENT_LIST_ITEM), ONE, (void *) ptr);
            ptr = next;
        }
        tl_free(sizeof(CENT_LIST_ITEM), ONE, (void *) ptr);
    }
}

/****************************************************************/
/****************************************************************/

void build_vertical_list(SGRID **grid3d)
{
    int i, ival;
    SGRID *grid = (*grid3d); // alias
    
    free_id_list(grid->vertical_list, grid->isize3);
    if (grid->vertical_list != NULL) {
        tl_free(sizeof(ID_LIST_ITEM *), grid->isize3, (void *) grid->vertical_list);
    }
    
    //  grid->nnodes_sur = count_surface_nodes();
    grid->isize3 = grid->nnodes_sur;
    
    grid->vertical_list = tl_alloc(sizeof(ID_LIST_ITEM *), grid->nnodes_sur);
    for (i = 0; i < grid->nnodes_sur; i++) {
        grid->vertical_list[i] = tl_alloc(sizeof(ID_LIST_ITEM), ONE);
        grid->vertical_list[i]->next = NULL;
    }
    for (i = 0; i < grid->nnodes; i++) {
        ival = find_vertical_segment(grid, i, grid->vertical_hash);
        Add_sorted_vertical_entry(grid, &grid->vertical_list[ival], i);
    }
}

/****************************************************************/
/****************************************************************/

int length_column_list(CENT_LIST_ITEM * head)
{
    int count = 0;
    
    while (head->next != NULL) {
        count++;
        head = head->next;
    }
    return count;
}

/****************************************************************/
/****************************************************************/

void free_midpt_hash(SGRID *grid, CENT_LIST_ITEM ** list)
{
    CENT_LIST_ITEM *ptr, *next;
    int i;
    
    /* check to see if we need to free anything */
    if (grid->is_allocated_midpt_hash == 0) {
        return;
    }
    
    for (i = 0; i < grid->hash_size; i++) {
        ptr = list[i];
        while (ptr->next != NULL) {
            next = ptr->next;
            tl_free(sizeof(CENT_LIST_ITEM), ONE, (void *) ptr);
            ptr = next;
        }
        tl_free(sizeof(CENT_LIST_ITEM), ONE, (void *) ptr);
    }
}

/****************************************************************/
/****************************************************************/

void build_midpt_hash(SGRID **grid3d)
{
    int i, j;
    int ie;
    int nd1, nd2;
    SGRID *grid = (*grid3d); // alias
    
    free_midpt_hash(grid, grid->midpt_hash);
    grid->num_midpts = 0;
    
    grid->is_allocated_midpt_hash = 1;
    
    /* num_midpts is the the number of edges (and hence midpoints) on the 2D surface mesh */
    grid->num_midpts = 0;
    if (grid->midpt_hash == NULL) {
        grid->midpt_hash = (CENT_LIST_ITEM **) tl_alloc(sizeof(CENT_LIST_ITEM *), grid->hash_size);
    }
    for (i = 0; i < grid->hash_size; i++) {
        grid->midpt_hash[i] = tl_alloc(sizeof(CENT_LIST_ITEM), ONE);
        grid->midpt_hash[i]->next = NULL;
    }
    
    /* loop over the columns and get the surface face */
    for (i = 0; i < grid->ncolumns; i++) {
        ie = grid->elem2d_sur[i];
        if (grid->elem2d[ie].nnodes == NDONTRI) {
            for (j = 0; j < grid->elem2d[ie].nnodes; j++) {
                nd1 = grid->elem2d[ie].nodes[grid->nd_on_TriEdge[j][0]];
                nd2 = grid->elem2d[ie].nodes[grid->nd_on_TriEdge[j][1]];
                Add_midpt_hash_entry(grid, grid->midpt_hash, nd1, nd2);
#ifdef _DEBUG
                if (nd1 < 0 || nd2 < 0) {
                    printf("icol: %d  2d element id: %d nd1: %d nd2: %d\n",i,ie,nd1,nd2);
                    tl_error("triangle edge nodes IDs jacked up!");
                }
#endif
            }
        } else {
            for (j = 0; j < grid->elem2d[ie].nnodes; j++) {
                nd1 = grid->elem2d[ie].nodes[grid->nd_on_QuadEdge[j][0]];
                nd2 = grid->elem2d[ie].nodes[grid->nd_on_QuadEdge[j][1]];
                Add_midpt_hash_entry(grid, grid->midpt_hash, nd1, nd2);
            }
        }
    }
    
    /* loop over all the surface nodes and add an entry for them */
    for (i = 0; i < grid->nnodes_sur; i++) {
        nd1 = grid->vertical_list[i]->id;
        nd2 = grid->vertical_list[i]->id;
#ifdef _DEBUG
        if (nd1 < 0 || nd2 < 0) {
            printf("surface node: %d nd1: %d nd2: %d\n",i,nd1,nd2);
            tl_error("surface nodes IDs jacked up!");
        }
#endif
        Add_midpt_hash_entry(grid, grid->midpt_hash, nd1, nd2);
    }
    
    //  for(ie=0; ie < nelem2d; ie++)
    //    {
    //      Add_edges = 1;
    //      /* verify that this is a free surface face */
    //      for(i=0; i<3; i++)
    //        {
    //          istr = node_flags[elem2d[ie].nodes[i]];
    //          if(str_values[istr].displacement.bc_flag != BCT_FREE_DIR)
    //            {
    //              printf("NOTE: 2D element #%d is not on free surface\n", ie);
    //              Add_edges = 0;
    //              break;
    //            }
    //        }
    //      if(Add_edges == 1)
    //        {
    //          for(i=0; i<3; i++)
    //            {
    //              nd1 = elem2d[ie].nodes[nd_on_2dedge[i][0]];
    //              nd2 = elem2d[ie].nodes[nd_on_2dedge[i][1]];
    //              Add_midpt_hash_entry(grid->midpt_hash, nd1, nd2);
    //            }
    //        }
    //    }
    
}

/****************************************************************/
/****************************************************************/
// CJT NOTE :: both column_hash and midpt_hash are of type CENT_LIST_ITEM
// CJT NOTE :: they both utilized column_hash_function to store
// CJT NOTE :: here, we are concerned with midpt_hash

void Add_midpt_hash_entry(SGRID *grid, CENT_LIST_ITEM ** head, int nd1, int nd2)
{
    CENT_LIST_ITEM *ptr;
    SVECT2D midpt, entry;
    int hash_val;
    
    
    midpt.x = (grid->node[nd1].x + grid->node[nd2].x) * 0.5;
    midpt.y = (grid->node[nd1].y + grid->node[nd2].y) * 0.5;
    
    hash_val = column_hash_function(grid, midpt);
    /* first check to see if there is no entry in this bin */
    if (head[hash_val] == NULL) {
        /* Null pointer means that we haven't filled this bin yet. Fire away. */
        push_column_hash_entry(grid, &head[hash_val], midpt);
        grid->num_midpts++;
        return;
    }
    
    ptr = head[hash_val];
    while (ptr->next != NULL) {
        entry = ptr->vect;
        /* check to see if this center is the same as a pre-existing entry */
        if ((fabs(entry.x - midpt.x) < DIFF_TOL) && (fabs(entry.y - midpt.y) < DIFF_TOL)) {
            return;   /* this point has already been hashed */
        }
        ptr = ptr->next;
    }
    /* Found a new entry; add it to the hashtable */
    push_midpt_hash_entry(grid, &head[hash_val], midpt);
    grid->num_midpts++;
    return;
}

/****************************************************************/
/****************************************************************/

void push_midpt_hash_entry(SGRID *grid, CENT_LIST_ITEM ** head, SVECT2D center)
{
    CENT_LIST_ITEM *ptr;
    
    ptr = tl_alloc(sizeof(CENT_LIST_ITEM), ONE);
    ptr->vect = center;
    ptr->index = grid->num_midpts;
    ptr->next = *head;
    *head = ptr;
    
    return;
}

/****************************************************************/
/****************************************************************/

void build_midpt_list(SGRID **grid3d)
{
    int i;
    int ie;
    int ival;
    int nd1, nd2;
    SVECT midpt;
    SGRID *grid = (*grid3d); // alias
    
    free_midpt_list(grid->midpt_list, grid->isize4);
    grid->midpt_list = tl_realloc(sizeof(ID_LIST_ITEM *), grid->num_midpts, grid->isize4, grid->midpt_list);
    grid->isize4 = grid->num_midpts;
    
    for (i = 0; i < grid->num_midpts; i++) {
        grid->midpt_list[i] = tl_alloc(sizeof(MIDPT_LIST_ITEM), ONE);
        grid->midpt_list[i]->next = NULL;
    }
    for (ie = 0; ie < grid->nelems3d; ie++) {
        for (i = 0; i < grid->elem3d[ie].nedges; i++) {
            nd1 = grid->elem3d[ie].nodes[grid->elem3d[ie].edges[i][0]];
            nd2 = grid->elem3d[ie].nodes[grid->elem3d[ie].edges[i][1]];
            
#ifdef _DEBUG
            if (nd1 < 0 || nd2 < 0) {
                printf("3d element id: %d nd1: %d nd2: %d\n",ie,nd1,nd2);
                tl_error("triangle edge nodes IDs jacked up!");
            }
#endif
            
            midpt.x = (grid->node[nd1].x + grid->node[nd2].x) * 0.5;
            midpt.y = (grid->node[nd1].y + grid->node[nd2].y) * 0.5;
            midpt.z = (grid->node[nd1].z + grid->node[nd2].z) * 0.5;
            
            ival = find_midpt_segment(grid, midpt, grid->midpt_hash);
            if (ival >= 0) {
                Add_sorted_midpt_entry(grid, &grid->midpt_list[ival], midpt, nd1, nd2);
            }
        }
    }
    
}

/****************************************************************/
/****************************************************************/
// NOTE :: Vertical edges are not hashed
// NOTE :: The column_hash_function is used for hashing columns, midpoints, etc, it's just called "column"
//         Here, however, the values it returns are used to find items for a midpoint list

int find_midpt_segment(SGRID *grid, SVECT midpt, CENT_LIST_ITEM ** list)
{
    SVECT2D center, entry;
    CENT_LIST_ITEM *head;
    int hash_val;
    
    center.x = midpt.x;
    center.y = midpt.y;
    
    hash_val = column_hash_function(grid, center);
    
    head = list[hash_val];
    while (head->next != NULL) {
        entry = head->vect;
        if ((fabs(entry.x - center.x) < DIFF_TOL) && (fabs(entry.y - center.y) < DIFF_TOL)) {
            /* found the match */
            return head->index;  // return the index of this column
        }
        head = head->next;
    }
    /* Didn't find a match; this edge must be a vertical edge (which are not hashed) */
    return -1;
}

/****************************************************************/
/****************************************************************/

void free_midpt_list(MIDPT_LIST_ITEM ** list, int number    /* number of midpt segments (the extent of this array) */
)
{
    MIDPT_LIST_ITEM *ptr, *next;
    int i;
    
    for (i = 0; i < number; i++) {
        ptr = list[i];
        while (ptr->next != NULL) {
            next = ptr->next;
            tl_free(sizeof(MIDPT_LIST_ITEM), ONE, (void *) ptr);
            ptr = next;
        }
        tl_free(sizeof(MIDPT_LIST_ITEM), ONE, (void *) ptr);
    }
}

/****************************************************************/
/****************************************************************/

void Add_sorted_midpt_entry(SGRID *grid, MIDPT_LIST_ITEM ** head, SVECT midpt, int nd1, int nd2)
{
    MIDPT_LIST_ITEM *ptr, *entry;
    MIDPT_LIST_ITEM *next = NULL;
    MIDPT_LIST_ITEM *prev = NULL;
    ID_LIST_ITEM *sptr = NULL;
    int snode1, snode2;
    int ival1, ival2;
    int temp;
    double z, this_z;
    
    ptr = tl_alloc(sizeof(MIDPT_LIST_ITEM), ONE);
    
    ptr->value[0] = UNSET_FLT;
    ptr->value[1] = UNSET_FLT;
    ptr->value[2] = UNSET_FLT;
    ptr->value[3] = UNSET_FLT;
    ptr->value[4] = UNSET_FLT;
    
    ptr->column1 = UNSET_INT;
    ptr->column2 = UNSET_INT;
    ptr->elem_upper[0] = UNSET_INT;
    ptr->elem_upper[1] = UNSET_INT;
    ptr->elem_lower[0] = UNSET_INT;
    ptr->elem_lower[1] = UNSET_INT;
    
    /* Find the surface node that lies above each of the endpoints */
    ival1 = find_vertical_segment(grid, nd1, grid->vertical_hash);
    sptr = grid->vertical_list[ival1];
    snode1 = sptr->id;
    
    ival2 = find_vertical_segment(grid, nd2, grid->vertical_hash);
    sptr = grid->vertical_list[ival2];
    snode2 = sptr->id;
    
    if (ival1 == ival2) {
        ptr->vertical = 1;
        /* put top node first */
        if (grid->node[nd1].z < grid->node[nd2].z) {
            temp = nd2;
            nd2 = nd1;
            nd1 = temp;
        }
        
    }
    else {
        ptr->vertical = 0;
        /* sort nodes by order of the surface nodes */
        if (snode1 > snode2) {
            temp = snode2;
            snode2 = snode1;
            snode1 = temp;
            temp = nd2;
            nd2 = nd1;
            nd1 = temp;
        }
    }
    
    ptr->vect = midpt;
    ptr->node1 = nd1;
    ptr->node2 = nd2;
    ptr->surf_node1 = snode1;
    ptr->surf_node2 = snode2;
    
    entry = *head;
    if (entry->next == NULL) {
        ptr->next = *head;
        ptr->prev = ptr;
        *head = ptr;
        return;
    }
    
    this_z = midpt.z;
    
    while (entry->next != NULL) {
        next = entry->next;
        z = entry->vect.z;
        
        /* Check to see if we've already included this edge */
        if (fabs(this_z - z) < DIFF_TOL) {
            tl_free(sizeof(MIDPT_LIST_ITEM), ONE, (void *) ptr);
            return;
        }
        /* Append entries with greater z-depth , this way the linked
         * list starts at the top, and we go down the column as
         * we follow the links
         */
        if (this_z > z) {
            if (prev == NULL) {
                *head = ptr;
                prev = *head;
            }
            else {
                prev->next = ptr;
            }
            ptr->next = entry;
            ptr->prev = prev;
            return;
        }
        prev = entry;
        entry = entry->next;
    }
    
    /* If we reach this point, then the entry to be added is at the end of the linked list */
    prev->next = ptr;
    ptr->next = entry;
    ptr->prev = prev;
    (*head)->prev = ptr;
    
    return;
}

/****************************************************************/
/****************************************************************/

/* Given the indices of the 2 nodes of an edge, set the value (e.g., pressure) on
 * the corresponding midpoint
 */

/* JLH   NOTE: THIS STILL NEEDS WORK to handle the various pressure values .. */

int set_value_midpt_list(SGRID *grid, MIDPT_LIST_ITEM ** head, int nd1, int nd2, double value, int index /* the index into the value array to fill */
)
{
    SVECT midpt;
    int ival;
    int temp;
    MIDPT_LIST_ITEM *ptr;
    
    if ((index < 0) || (index > 4)) {
        tl_error("set_value_midpt_list error: invalid index");
    }
    
    /* sort note indices, smallest first */
    if (nd1 > nd2) {
        temp = nd2;
        nd2 = nd1;
        nd1 = temp;
    }
    
    midpt.x = (grid->node[nd1].x + grid->node[nd2].x) * 0.5;
    midpt.y = (grid->node[nd1].y + grid->node[nd2].y) * 0.5;
    midpt.z = (grid->node[nd1].z + grid->node[nd2].z) * 0.5;
    
    ival = find_midpt_segment(grid, midpt, grid->midpt_hash);
    
    if (ival < 0) {
        /* Invalid edge; probably a vertical segment */
        return -1;
    }
    
    ptr = head[ival];
    
    while (ptr->next != NULL) {
        if ((ptr->node1 == nd1) && (ptr->node2 == nd2))
        /* Found the match */
        {
            ptr->value[index] = value;
            return 1;
        }
        ptr = ptr->next;
    }
    
    /* Fell through and didn't find a midpoint in the hash */
    return -1;
}

/****************************************************************/
/****************************************************************/

/* Given the indices of the 2 nodes of an edge, get the value (e.g., pressure) on
 * the corresponding midpoint
 * Returns 1 for success, -1 if the edge was not found
 */
int get_value_midpt_list(SGRID *grid, MIDPT_LIST_ITEM ** head, int nd1, int nd2, double *value)
{
    SVECT midpt;
    int ival;
    MIDPT_LIST_ITEM *ptr;
    
    midpt.x = (grid->node[nd1].x + grid->node[nd2].x) * 0.5;
    midpt.y = (grid->node[nd1].y + grid->node[nd2].y) * 0.5;
    midpt.z = (grid->node[nd1].z + grid->node[nd2].z) * 0.5;
    
    ival = find_midpt_segment(grid, midpt, grid->midpt_hash);
    if (ival < 0) {
        /* Couldn't find midpt segment --- most likely because we are on a vertical edge */
        return -1;
    }
    
    ptr = head[ival];
    
    while (ptr->next != NULL) {
        if (   ((ptr->node1 == nd1) && (ptr->node2 == nd2))
            || ((ptr->node1 == nd2) && (ptr->node2 == nd1)) )
        /* Found the match */
        {
            value[0] = ptr->value[0];
            value[1] = ptr->value[1];
            value[2] = ptr->value[2];
            value[3] = ptr->value[3];
            value[4] = ptr->value[4];
            return 1;
        }
        ptr = ptr->next;
    }
    
    /* Fell through and didn't find a midpoint in the hash */
    return -1;
}

/****************************************************************/
/****************************************************************/

/* Associate elements to edges/midpts
 * MUST have already build the midpt list by invoking:
 *     build_midpt_hash();
 *     build_midpt_list();
 * AND must have already built the prism column list by invoking:
 *     build_column_hash();
 *     build_column_list();
 * ALSO need to have built the linked lists of vertical node segments by:
 *     build_vertical_hash();
 *     build_vertical_list();
 *     (this last is used in the function above_below)
 */

void build_midpt_connectivity(SGRID *grid)
{
    int icol, ie, j;
    int ival;
    int nd1, nd2;
    SVECT midpt;
    
    ID_LIST_ITEM *ptr;
    
    /* Loop over the prism columns */
    for (icol = 0; icol < grid->ncolumns; icol++) {
        ptr = grid->column_list[icol];
        while (ptr->next != NULL) {
            ie = ptr->id;
            for (j = 0; j < grid->elem3d[ie].nedges; j++) {
                nd1 = grid->elem3d[ie].nodes[grid->elem3d[ie].edges[j][0]];
                nd2 = grid->elem3d[ie].nodes[grid->elem3d[ie].edges[j][1]];
                
                midpt.x = (grid->node[nd1].x + grid->node[nd2].x) * 0.5;
                midpt.y = (grid->node[nd1].y + grid->node[nd2].y) * 0.5;
                midpt.z = (grid->node[nd1].z + grid->node[nd2].z) * 0.5;
                
                /* make sure that this edge is hashed (if it is a vertical edge it won't be) */
                ival = find_midpt_segment(grid, midpt, grid->midpt_hash);
                if (ival >= 0) {
                    set_midpt_elem(grid, &grid->midpt_list[ival], icol, midpt, nd1, nd2, ie);
                }
            }
        }
        
        ptr = ptr->next;
    }
    
}

/****************************************************************/
/****************************************************************/

void set_midpt_elem(SGRID *grid, MIDPT_LIST_ITEM ** head, int icol, SVECT midpt, int nd1, int nd2, int ie)
{
    int temp;
    MIDPT_LIST_ITEM *ptr;
    int ic;
    int iwhich;
    
    /* sort note indices, smallest first */
    if (nd1 > nd2) {
        temp = nd2;
        nd2 = nd1;
        nd1 = temp;
    }
    
    midpt.x = (grid->node[nd1].x + grid->node[nd2].x) * 0.5;
    midpt.y = (grid->node[nd1].y + grid->node[nd2].y) * 0.5;
    midpt.z = (grid->node[nd1].z + grid->node[nd2].z) * 0.5;
    
    ptr = *head;
    
    while (ptr->next != NULL) {
        if ((ptr->node1 == nd1) && (ptr->node2 == nd2)) {
            /* Found an edge that matches */
            /* Now, make sure that this element has a vertical face on this side of
             * the column. We can do this by checking the nodes of the element
             * and seeing if there is another point above/below one of the segment
             * endpoints.
             */
            iwhich = above_or_below(grid, ie, nd1, nd2);
            if (iwhich == 0) {    /* the element has no vertical face on this side of column */
                return;
            }
            
            /* Now figure out which column we are in */
            /* On the first entry, we set column1 */
            if (ptr->column1 == UNSET_INT) {
                ptr->column1 = icol;
                ic = 0;
            }
            /* check to see if this column is the same as what we've already seen */
            else if (ptr->column1 == icol) {
                ic = 0;
            }
            /* This must be the second column */
            else {
                ptr->column2 = icol;
                ic = 1;
            }
            if (iwhich > 0) { /* this element lies above the midpt */
                /* santity check */
                if (ptr->elem_upper[ic] != UNSET_INT) {
                    tl_error("ERROR. We seem to have two elements above midpoint.\n");
                }
                ptr->elem_upper[ic] = ie;
            }
            else {    /* this element lies below the midpt */
                
                /* santity check */
                if (ptr->elem_lower[ic] != UNSET_INT) {
                    tl_error("ERROR. We seem to have two elements below midpoint.\n");
                }
                ptr->elem_lower[ic] = ie;
            }
            return;
        }
        ptr = ptr->next;
    }
    tl_error("ERROR: fall through trying to attach element to midpt. \n");
}

#if 0
/* Determine which element lies above the other in the column
 * Returns 1 if ie1 is the upper element,
 *        -1 if ie2 is the lower element,
 *         0 if we are confused
 */
int find_upper_element(int ie1, int ie2, int nd1, int nd2)
{
    int i;
    int iseg1, iseg2;
    int is;
    int this_node;
    double z_max = -999999999999999.9;
    
    /* we look at each node of the elements, and determine if they lie in the vertical
     * line of the endpoints of our segment. This gives of the planar face of the
     * vertical column. Then we take the largest z-value to determine
     * which element is on top.
     */
    
    iseg1 = find_vertical_segment(nd1, grid->vertical_hash);
    iseg2 = find_vertical_segment(nd2, grid->vertical_hash);
    
    /* Loop over the first element and find the largest z-value on the face */
    for (i = 0; i < NDPRELM; i++) {
        this_node = elem3d[ie1].nodes[i];
        is = find_vertical_segment(this_node, grid->vertical_hash);
        /* verify that the node is on the proper face */
        if ((is != iseg1) && (is != iseg2)) {
            continue;
        }
        z_max = MAX(z_max, node[this_node].z);
    }
    
    /* Now loop over the second element */
    for (i = 0; i < NDPRELM; i++) {
        this_node = elem3d[ie2].nodes[i];
        is = find_vertical_segment(this_node, grid->vertical_hash);
        /* verify that the node is on the proper face */
        if ((is != iseg1) && (is != iseg2)) {
            continue;
        }
        if (node[this_node].z > z_max) {
            return -1;
        }
    }
    return 1;
}
#endif

/****************************************************************/
/****************************************************************/

/* Given an element and two nodes of on edge, determine if the element has a vertical face
 * that contains the edge, and if so, if the element lies above or below this edge
 *   Returns:
 *      0  if the element doesn't have a vertical face
 *      1  if the element is above the edge
 *     -1  if the element is below the edge
 */
int above_or_below(SGRID *grid, int ie, int nd1, int nd2)
{
    int i, j;
    int hash[2];
    int my_nodes[2];
    int jnode;
    int ihash;
    
    my_nodes[0] = nd1;
    my_nodes[1] = nd2;
    
    hash[0] = find_vertical_segment(grid, my_nodes[0], grid->vertical_hash);
    hash[1] = find_vertical_segment(grid, my_nodes[1], grid->vertical_hash);
    //      data[i].y = node[elem2d[iel].nodes[i]].y;
    
    for (i = 0; i < NDPRELM; i++) {
        jnode = grid->elem3d[ie].nodes[i];
        if ((jnode == nd1) || (jnode == nd2)) {
            continue;
        }
        ihash = find_vertical_segment(grid, jnode, grid->vertical_hash);
        /* loop over the two endpoints */
        for (j = 0; j < 2; j++) {
            if (ihash == hash[j]) {
                /* This node of the element lies above/below the endpoint */
                if (grid->node[jnode].z > grid->node[my_nodes[j]].z) {  /* above */
                    return 1;
                }
                else if (grid->node[jnode].z < grid->node[my_nodes[j]].z) { /* below */
                    return -1;
                }
                else {
                    /* sanity check ... this should not happen! */
                    tl_error("ERROR in above_below\n");
                }
            }
        }
    }
    return 0;
}

/****************************************************************/
/****************************************************************/

/* Given the nodes of an endpoint, return the item in the link list associated with
 * this edge/midpoint
 */
MIDPT_LIST_ITEM *get_midpt_entry(SGRID *grid, MIDPT_LIST_ITEM ** head, int nd1, int nd2)
{
    int temp;
    int ival;
    SVECT midpt;
    MIDPT_LIST_ITEM *ptr;
    
    /* sort note indices, smallest first */
    if (nd1 > nd2) {
        temp = nd2;
        nd2 = nd1;
        nd1 = temp;
    }
    
    midpt.x = (grid->node[nd1].x + grid->node[nd2].x) * 0.5;
    midpt.y = (grid->node[nd1].y + grid->node[nd2].y) * 0.5;
    midpt.z = (grid->node[nd1].z + grid->node[nd2].z) * 0.5;
    
    ival = find_midpt_segment(grid, midpt, grid->midpt_hash);
    
    ptr = head[ival];
    
    while (ptr->next != NULL) {
        if ((ptr->node1 == nd1) && (ptr->node2 == nd2))
        /* Found the match */
        {
            return ptr;
        }
        ptr = ptr->next;
    }
    return NULL;
}

/*!
 * \brief Function to allows us to sort list of ints
 *
 * \param *p1 Pointer to First Entry
 * \param *p2 Pointer to Second Entry
 *  \return 1 If P1>P2, -1 If P1<P2, 0 If P1==P2
 */

/****************************************************************/
/****************************************************************/

void free_column_hash2d(SGRID *grid, CENT_LIST_ITEM ** list)
{
    CENT_LIST_ITEM *ptr, *next;
    int i;
    
    /* check to see if we need to free anything */
    if (grid->is_allocated_column_hash2d == 0) {
        return;
    }
    for (i = 0; i < grid->hash_size; i++) {
        ptr = list[i];
        while (ptr->next != NULL) {
            next = ptr->next;
            tl_free(sizeof(CENT_LIST_ITEM), ONE, (void *) ptr);
            ptr = next;
        }
        tl_free(sizeof(CENT_LIST_ITEM), ONE, (void *) ptr);
    }
}

/****************************************************************/
/****************************************************************/

void global_free_columns(SGRID *grid) {
    
    
    /* free column hashes */
    free_column_hash(grid, grid->column_hash);
    if (grid->column_hash != NULL) {
        tl_free(sizeof(CENT_LIST_ITEM *), grid->hash_size, (void *) grid->column_hash);
    }
    free_column_hash2d(grid, grid->column_hash2d);
    if (grid->column_hash2d != NULL) {
        tl_free(sizeof(CENT_LIST_ITEM *), grid->hash_size, (void *) grid->column_hash2d);
    }
    
    /* free column lists */
    free_id_list(grid->column_list, grid->isize_ncolumns);
    if (grid->column_list != NULL) {
        tl_free(sizeof(ID_LIST_ITEM *), grid->isize_ncolumns, (void *) grid->column_list);
    }
    free_id_list(grid->column_list2d, grid->isize_ncolumns);
    if (grid->column_list2d != NULL) {
        tl_free(sizeof(ID_LIST_ITEM *), grid->isize_ncolumns, (void *) grid->column_list2d);
    }
    /* cjt */
    free_id_list(grid->sidewall_list, grid->num_midpts);
    if (grid->sidewall_list != NULL) {
        tl_free(sizeof(ID_LIST_ITEM *), grid->num_midpts, (void *) grid->sidewall_list);
    }
    
    if (grid->isize_sur_elems > 0) {
        tl_free(sizeof(int), grid->isize_sur_elems, (void *) grid->elem3d_sur);
        tl_free(sizeof(int), grid->isize_sur_elems, (void *) grid->elem2d_sur);
    }
    if (grid->isize_bed_elems > 0) {
        tl_free(sizeof(int), grid->isize_bed_elems, (void *) grid->elem3d_bed);
        tl_free(sizeof(int), grid->isize_bed_elems, (void *) grid->elem2d_bed);
    }
    
    if (grid->isize_sidewall_elems > 0) {
        tl_free(sizeof(int), grid->isize_sidewall_elems, (void *) grid->elem3d_sidewall);
        tl_free(sizeof(int), grid->isize_sidewall_elems, (void *) grid->elem2d_sidewall);
    }
    
    /* free maps */
    if (grid->nodeID_2d_to_3d_sur != NULL) { // this is a 2d grid
        grid->nodeID_2d_to_3d_sur = (int *) tl_free(sizeof(int), grid->old_max_nnodes_sur, grid->nodeID_2d_to_3d_sur);
    }
    
    if (grid->nodeID_2d_to_3d_bed != NULL) { // this is a 2d grid
        grid->nodeID_2d_to_3d_bed = (int *) tl_free(sizeof(int), grid->old_max_nnodes_bed, grid->nodeID_2d_to_3d_bed);
    }
    
    if (grid->nodeID_3d_to_2d_sur != NULL) { // this is a 3d grid
        grid->nodeID_3d_to_2d_sur = (int *) tl_free(sizeof(int), grid->max_nnodes, grid->nodeID_3d_to_2d_sur);
    }
    
    if (grid->nodeID_3d_to_2d_bed != NULL) { // this is a 3d grid
        grid->nodeID_3d_to_2d_bed = (int *) tl_free(sizeof(int), grid->max_nnodes, grid->nodeID_3d_to_2d_bed);
    }
    
    /* free vertical hash and vertical lists */
    free_vertical_hash(grid, grid->vertical_hash);
    if (grid->vertical_hash != NULL) {
        tl_free(sizeof(CENT_LIST_ITEM **), grid->hash_size, (void *) grid->vertical_hash);
    }
    free_id_list(grid->vertical_list, grid->isize3);
    if (grid->vertical_list != NULL) {
        tl_free(sizeof(ID_LIST_ITEM *), grid->isize3, (void *) grid->vertical_list);
    }
    
    /* free midpt hash */
    free_midpt_hash(grid, grid->midpt_hash);
    if (grid->midpt_hash != NULL) {
        tl_free(sizeof(CENT_LIST_ITEM *), grid->hash_size, (void *) grid->midpt_hash);
    }
    
    /* free midpt list */
    free_midpt_list(grid->midpt_list, grid->isize4);
    if (grid->midpt_list != NULL) {
        tl_free(sizeof(ID_LIST_ITEM *), grid->isize4, (void *) grid->midpt_list);
    }
    
#ifdef _MESSG
    // column_adpt_unref_clean();
#endif
    
}


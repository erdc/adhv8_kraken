#include "global_header.h"

int get_node_id(SIO info, char **data, int nnodes) {

    int ind = read_int_field(info,data);
    if (ind > nnodes) {
        io_read_error(info, "Bad node number.", TRUE);
    }
    ind--;
    return ind;
}

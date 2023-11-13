/* Given a list of 4 integers, one of which is identical to one of the others, return a list of the 3 unique integers */
#include "global_header.h"

void tl_find_common_int(int *int_list) {

  if(int_list[0]==int_list[1]) {
	  int_list[1]=int_list[2];
	  int_list[2]=int_list[3];
  }
  else if(int_list[0]==int_list[2]) {
	  int_list[2]=int_list[3];
  }
  else if(int_list[1]==int_list[2]) {
	  int_list[2]=int_list[3];
  }
  else {
	  /* this means that the 4th integer in redundant, so just return the original list */
  }
}

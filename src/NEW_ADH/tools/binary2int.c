#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int fromBinary(const char *s) {
  return (int) strtol(s, NULL, 2);
}

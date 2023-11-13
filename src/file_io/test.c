
#include <stdio.h>

extern void create_fileio(
  int *
);
extern void close_fileio(
  int *
);
extern int write_fileio_dbl(
  int *,
  double *,
  int *
);
extern int read_fileio_dbl(
  int *,
  double *,
  int *
);

int main(
  int argc,
  char **argv
)
{
  int num;
  double mydata[940000];
  double readdata[500];
  int nrecv;
  int datacount = 940000;
  int i;

  num = 20;

  for(i = 0; i < datacount; i++)
    mydata[i] = (double)i;

  create_fileio(&num);

  for(i = 0; i < 10; i++)
    {
      printf("writing data ");
      write_fileio_dbl(&num, mydata, &datacount);
      printf("done\n");

    }

  close_fileio(&num);

}

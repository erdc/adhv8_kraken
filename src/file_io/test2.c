
#include <stdio.h>

extern void create_fileio(
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
  double mydata[100];
  double readdata[950000];
  int nrecv;
  int datacount = 100;
  int i, j;

  num = 20;

  for(i = 0; i < 10; i++)
    {
      read_fileio_dbl(&num, readdata, &nrecv);

      printf("%d: recived %d doubles\n", i, nrecv);

      for(j = 0; j < 5; j++)
	printf("%f\n", readdata[j]);

      printf("\n");
    }
}

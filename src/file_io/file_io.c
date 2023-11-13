
/*
   Purpose:  set of library functions that mimic sockets but using plain
   files.

   Author:   JRB / Sat Jul  1 23:45:42 CDT 2006

   Plan:  Communicating via files will require two files.  One data file
   and the other is a semaphore (thanx Jeff H.) file.  The semaphore
   file contains the status of the data file.

   created file ready for writing             1
   data in file ready for reading             2

   to open  1) create files
   2) change semaphore to 1

   to write 1) check semaphore to see if = 1 (ready for writing)
   2) change semaphre to -2 to indicate data is being written
   3) write data
   4) change semaphore to 2 to indicate data is complete

   to read  1) check semaphore to see if = 2 (reading for reading)
   2) change semaphore to -1 to indicate data is being read
   3) read data
   4) change semaphore to 1 to indicate that data read is complete
 */

/*
   Basic functions
   create fileio  (0)
   close fileio   (0)
   read fileio    (0)
   write fileio   (0)
   status fileio  (0)
 */

#include <stdio.h>
#include <fcntl.h>
/*#include <unistd.h> */
#ifndef WINDOWS
#include <unistd.h>
#else
#include <windows.h>
#endif
#include <sys/types.h>
#include <sys/stat.h>
#include "header_file_io.h"

#define SLEEP_TIME 5
#define TRUE 1
#define FALSE 0

/* forward declarations to avoid including all the headers */
FILE * io_fopen(const char *, const char *, const int);

void create_fileio(
  int *socknum
)
{
  char fname[80];
  char dname[80], semaphore[80];
  FILE *fd_dname, *fd_semaphore;
  struct stat buf;
  int istat;

  sprintf(fname, "../file_sock%d", *socknum);

  sprintf(dname, "%s.fsd", fname);
  sprintf(semaphore, "%s.sem", fname);

  istat = stat(dname, &buf);

  fprintf(stderr, "create_fileio: istat %d\n", istat);

  if (istat == -1)
    {
      fd_dname = io_fopen(dname, "w", TRUE);
      fd_semaphore = io_fopen(semaphore, "w", TRUE);

      fprintf(fd_semaphore, "1");

      fclose(fd_dname);
      fclose(fd_semaphore);
    }
}

void close_fileio(
  int *socknum
)
{
  char fname[80];
  char dname[80], semaphore[80];
  FILE *fd_semaphore;
  int sem_status;

  sprintf(fname, "../file_sock%d", *socknum);

  sprintf(dname, "%s.fsd", fname);
  sprintf(semaphore, "%s.sem", fname);

  /*
   ** check sem status
   **   to write 1) check semaphore to see if = 1 
   **    (finished reading OK to delete)
   */

  sem_status = -99;

  do
    {
      fd_semaphore = io_fopen(semaphore, "r", FALSE);

      if (fd_semaphore != NULL)
        {
          fscanf(fd_semaphore, "%d", &sem_status);
          fclose(fd_semaphore);
        }

      if (sem_status != 1)
        {
          fprintf(stderr, "waiting for semaphore to close\n");
#ifndef WINDOWS
          sleep(SLEEP_TIME);
#else
          Sleep(SLEEP_TIME);
#endif
        }
    }
  while (sem_status != 1);

  unlink(dname);
  unlink(semaphore);
}

int write_fileio_dbl(
  int *socknum,
  double *dataout,
  int *nsend
)
{
  char fname[80];
  char dname[80], semaphore[80];
  FILE *fd_dname, *fd_semaphore;

  int sem_status, i;
  /*int fstatus; */

  create_fileio(socknum);

  sprintf(fname, "../file_sock%d", *socknum);

  sprintf(dname, "%s.fsd", fname);
  sprintf(semaphore, "%s.sem", fname);

  /*
   ** check sem status
   **   to write 1) check semaphore to see if = 1 (ready for writing)
   */

  sem_status = -99;

  do
    {
      fd_semaphore = io_fopen(semaphore, "r", FALSE);

      if (fd_semaphore != NULL)
        {
          fscanf(fd_semaphore, "%d", &sem_status);
          fclose(fd_semaphore);
        }

      if (sem_status != 1)
        {
          fprintf(stderr, "waiting for semaphore to write\n");
#ifndef WINDOWS
          sleep(SLEEP_TIME);
#else
          Sleep(SLEEP_TIME);
#endif
        }
    }
  while (sem_status != 1);

  /*  
   **  Now cleared to write
   **    2) change semaphre to -2 to indicate data is being written
   */

  fd_semaphore = io_fopen(semaphore, "w", TRUE);
  fprintf(fd_semaphore, "-2");
  fclose(fd_semaphore);

  /*
   ** 3) write data
   */
  fd_dname = io_fopen(dname, "w", TRUE);

  /* fstatus = fwrite(dataout, sizeof(double), *nsend, fd_dname); */

  for (i = 0; i < *nsend; i++)
    {
      if (i < 5 || i > *nsend - 5)
        {
          fprintf(stderr, "node %d val %f\n", i, dataout[i]);
        }
      fprintf(fd_dname, "%f\n", dataout[i]);
    }

  fclose(fd_dname);

  /*
   **  4) change semaphore to 2 to indicate data is complete
   */
  fd_semaphore = io_fopen(semaphore, "w", TRUE);
  fprintf(fd_semaphore, "2");
  fclose(fd_semaphore);

  /* return(fstatus); */
  return (1);
}

int read_fileio_dbl(
  int *socknum,
  double *datain,
  int *nrecv
)
{
  char fname[80];
  char dname[80], semaphore[80];
  FILE *fd_dname, *fd_semaphore;

  int sem_status;
  int dblcount;
  size_t readsize;
  double dbl_val;

  fprintf(stderr, "socknum %d\n", *socknum);
  sprintf(fname, "../file_sock%d", *socknum);

  sprintf(dname, "%s.fsd", fname);
  sprintf(semaphore, "%s.sem", fname);

  /*
   ** check sem status
   **   to write 1) check semaphore to see if = 2 (ready for reading)
   */
  sem_status = -99;

  do
    {
      fd_semaphore = io_fopen(semaphore, "r", FALSE);

      if (fd_semaphore != NULL)
        {
          fscanf(fd_semaphore, "%d", &sem_status);
          fclose(fd_semaphore);
        }

      if (sem_status != 2)
        {
          fprintf(stderr, "waiting for semaphore to read\n");
#ifndef WINDOWS
          sleep(SLEEP_TIME);
#else
          Sleep(SLEEP_TIME);
#endif
        }
    }
  while (sem_status != 2);

  fprintf(stderr, "Cleared to read\n");

  /*
   **   cleared to read
   ** 2) change semaphore to -1 to indicate data is being read
   */
  fd_semaphore = io_fopen(semaphore, "w", TRUE);
  fprintf(fd_semaphore, "-1");
  fclose(fd_semaphore);
  fprintf(stderr, "Changed semaphore\n");

  /*
   **  3) read data
   */
  dblcount = 0;
  fd_dname = io_fopen(dname, "r", TRUE);
  fprintf(stderr, "Opened file --- begin read loop\n");

  do
    {
      /* readsize = fread ( datain, sizeof(double), 1, fd_dname ); */
      readsize = fscanf(fd_dname, "%lf", &dbl_val);

      if (dblcount < 5)
        {
          fprintf(stderr, "dblcount: %d %lu %f\n", dblcount, readsize, dbl_val);
        }

      if (readsize == 1)
        {
          datain[dblcount] = dbl_val;
          dblcount++;
        }
    }
  while (readsize == 1);

  fprintf(stderr, "last dblcount: %d %lu %f\n", dblcount, readsize, dbl_val);

  fclose(fd_dname);

  *nrecv = dblcount;

  /*
   ** 4) change semaphore to 1 to indicate that data read is complete
   */
  fd_semaphore = io_fopen(semaphore, "w", TRUE);
  fprintf(fd_semaphore, "1");
  fclose(fd_semaphore);

  return (dblcount);
}

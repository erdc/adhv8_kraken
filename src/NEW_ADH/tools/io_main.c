#include <ctype.h>
#include "adh.h"

/*!
   \file newio_main.c
   \brief Basic utility functions for file I/O.
 */

/*!
   \brief Open a file and gives error if not successful

   Given a filename, this function attempts to open the file.
   If not successful, and told to abort, a error message is provided
   before exiting the program.

   \param[in] filename The path and name of file to open
   \param[in] mode The mode to open the file (e.g. "r" for reading)
   \param[in] abort If abort == TRUE, then execution is terminated, when
                opening the file fails. If abort == FALSE, then you are
                responsible for checking the file pointer before use.
   \return A pointer to the open file
 */
FILE * io_fopen(const char *filename, const char *mode, const int abort)
{
  FILE *fp = NULL;
  char *ptr = NULL;
  
  assert(filename && strlen(filename));
  assert(mode && strlen(mode));

  fp = fopen(filename, mode);
  /* if successful, or failure and NOT aborting */
  if (fp || (!fp && !abort))
    {
      /* return file pointer */
      return fp;
    }

  /* determine if opening file for reading or writing */ 
  ptr = strpbrk(mode, "w");
  fprintf(stderr, "\n Unable to find and open file for %s: %s",
    (ptr ? "writing" : "reading"), filename);
  exit(0);
  return NULL;
}

/*!
   \brief Open a file and prepare dataset

   Given a filename, this function opens the file and then
   prepares dataset by writing the header information.

   \param[in] filename The path and name of file to open
   \param[in] ps_flag The data set flag (PS_FLAG_xxx)
   \return A pointer to the open dataset file
 */
FILE * io_fopen_prep_ds(SMODEL_SUPER *mod, const char *filename, const int ps_flag, const int grid_flag)
{
  FILE *fp = NULL;
  
  assert(filename);
  fp = io_fopen(filename, "w", TRUE);
  print_header(mod, fp, ps_flag);
  return fp;
}

/*!
   \brief Open a file and prepare constituent dataset

   Given a filename, this function opens the file and then
   prepares dataset by writing the header information.

   \param[in] filename The path and name of file to open
   \param[in] ps_flag The data set flag (PS_FLAG_xxx)
   \param[in] index The constituent / variable index or UNSET_INT if not applicable
   \return A pointer to the open dataset file
 */
FILE * io_fopen_prep_ds_index(SMODEL_SUPER *mod, const char *filename, const int ps_flag, const int index, const int grid_flag)
{
  FILE *fp = NULL;
  
  assert(filename);
  fp = io_fopen(filename, "w", TRUE);
  mod->itrns = index;
  print_header(mod, fp, ps_flag);
  return fp;
}

/*!
   \brief Build up a filename from base and component strings
   
   Given a base name and extension, create and return a filename.

   \param[out] res The string buffer to hold the resulting filename
   \param[in] buffer_size The buffer size of the 
   \param[in] base The file base name
   \param[in] part The extra part to concatenate
   \return A pointer to the resulting filename (same as first variable)
 */
char * build_filename(char *res, const int buffer_size, const char *base,
  const char *part)
{
  return build_filename3(res, buffer_size, base, part, UNSET_INT, "", UNSET_INT,
    "", UNSET_INT);
}

/*!
   \brief Build up a filename from base and component strings
   
   Given a base name and various extensions and numbers, create and return
   a filename.

   \param[out] res The string buffer to hold the resulting filename
   \param[in] buffer_size The buffer size of the 
   \param[in] base The file base name
   \param[in] part1 The first extra part to concatenate
   \param[in] part1_num The first number (index) to concatenate or UNSET_INT if not applicable 
   \param[in] part2 The second extra part to concatenate
   \param[in] part2_num The second number (index) to concatenate or UNSET_INT if not applicable 
   \return A pointer to the resulting filename (same as first variable)
 */
char * build_filename2(char *res, const int buffer_size, const char *base,
  const char *part1, const int part1_num, const char *part2, const int part2_num)
{
  return build_filename3(res, buffer_size, base, part1, part1_num, part2, part2_num,
    "", UNSET_INT);
}

/*!
   \brief Build up a filename from base and component strings
   
   Given a base name and various extensions and numbers, create and return
   a filename.

   \param[out] res The string buffer to hold the resulting filename
   \param[in] buffer_size The buffer size of the 
   \param[in] base The file base name
   \param[in] part1 The first extra part to concatenate
   \param[in] part1_num The first number (index) to concatenate or UNSET_INT if not applicable 
   \param[in] part2 The second extra part to concatenate
   \param[in] part2_num The second number (index) to concatenate or UNSET_INT if not applicable 
   \param[in] part3 The third extra part to concatenate
   \param[in] part3_num The third number (index) to concatenate or UNSET_INT if not applicable 
   \return A pointer to the resulting filename (same as first variable)
 */
char * build_filename3(char *res, const int buffer_size, const char *base,
  const char *part1, const int part1_num, const char *part2, const int part2_num,
  const char *part3, const int part3_num)
{
  char num1[MAXLINE] = "", num2[MAXLINE] = "", num3[MAXLINE] = "";
  int size;

  assert(res && base && part1 && part2 && part3);
  size = strlen(base) + strlen(part1) + strlen(part2) + strlen(part3);
  if (part1_num != UNSET_INT)
    {
      sprintf(num1, "%d", part1_num);
      size += strlen(num1);
    }
  if (part2_num != UNSET_INT)
    {
      sprintf(num2, "%d", part2_num);
      size += strlen(num2);
    }
  if (part3_num != UNSET_INT)
    {
      sprintf(num3, "%d", part3_num);
      size += strlen(num3);
    }
  if (size >= buffer_size)
    {
      fprintf(stderr, "Filename exceeds buffer size of string array.");
      exit(0);
    }
  sprintf(res, "%s%s%s%s%s%s%s", base, part1, num1, part2, num2, part3, num3);
  return res;
}

/*!
   \brief Report a warning or error message for invalid input.

   Prints the warning/error message 'msg' and the contents of 'line' (the
   original full input line) to stderr. If abort is true, then program
   execution is terminated.

   \param[in] line The complete original input line being parsed
   \param[in] msg The error message
   \param[in] abort If abort == TRUE, then execution is terminated,
                otherwise continue on and hope for the best.
   \return TRUE on success
 */
int io_read_error(SIO io, const char *msg, const int abort)
{ 
  int myid;
#ifdef _MESSG
  myid = messg_comm_rank(cstorm_comm); // cjt :: only works for 1 grid in CSTORM :: MPI_COMM_WORLD);
#else
  myid = 0;
#endif

  assert(strlen(msg));
  if (myid <= 0)
    {
      if (strlen(io.cur_filename) && strlen(io.cur_line))
        {
          fprintf(stderr, " %s\n The file line (in \"%s\") was:\n   %s", msg, io.cur_filename, io.cur_line);
        }
      else
        {
          /* the filename and line were not saved prior to calling this error */
          fprintf(stderr, " %s", msg);
        }
    }
  if (abort)
    {
      exit(0);
    }
  return TRUE;
}

/*!
   \brief Copy the file name and file line being processed.

   Copy the file name and file line being processed. Information will be
   reported to the user when an error occurs.

   \param[in] fp Pointer to the current file (necessary to check buffer)
   \param[in] filename The name of the current file
   \param[in] line The current line
   \return TRUE on success
 */
int io_save_line(SIO *io, FILE *fp, const char *filename, const char *line)
{
  if (filename)
    {
      strncpy(io->cur_filename, filename, MAXLINE);
    }
  else
    {
      sprintf(io->cur_filename, " "); /* cjt */
    }
  if (line)
    {
      strncpy(io->cur_line, line, MAXLINE);
    }
  else
    {
      sprintf(io->cur_line, " "); /* cjt */
    }

  /* perform a check to ensure that the buffer was big enough to get the
     entire file line */
  if (fp && !feof(fp) && line && strlen(line) && !strchr(line, '\n'))
    {
      /* MAXLINE should be increased if file is valid and we reach this */
      char msg[MAXLINE];
      sprintf(msg, "File line is longer than the internal buffer; Can only read "
        "first %d characters of line.", MAXLINE);
      io_read_error(*io, msg, TRUE);
      return FALSE;
    }
  return TRUE;
}

/*!
   \brief Strip comments from file line

   Strips all characters from a file line after a comment symbol and
   returns resulting string length

   \param[in] fileline The file line to remove comments from
   \return String length of result
 */
int strip_comments(char *fileline)
{
  char *ptr = strpbrk(fileline, COMMENT_CHARACTERS);
  
  if (ptr)
    {
      *ptr = '\0';
    }
  return strlen(fileline);
}

/*!
   \brief Increment string pointer past leading white space

   \param[in] string A pointer to a string array
   \return TRUE on success
 */
int strip_white(char **string)
{
  while (**string != '\0' && isspace(**string))
    {
      (*string)++;
    }
  return TRUE;
}

/*!
   \brief Get text version of index
   
   Text string indicating an index is useful for file read error messages

   \param[in] index Number to convert into text
   \param[in] text String buffer to write result to
   \return TRUE on success, FALSE for failure
 */
int get_index_as_text(const int index, char *text)
{
  assert(text); /* should exist already */
  assert(index > 0); /* indices should be positive */
  switch (index)
    {
      case 1:
        sprintf(text, "1st");
        break;
      case 2:
        sprintf(text, "2nd");
        break;
      case 3:
        sprintf(text, "3rd");
        break;
      default: /* valid for 4 through 20 */
        sprintf(text, "%dth", index);
        break;
    }
  return (index > 0);
}

/*!
   \brief Convert text to uppercase

   Convert a text string to uppercase characters, in place (no new string buffer)

   \param[out] text String buffer to convert
   \return A pointer to the resulting text (same as first variable)
 */
char * convert_to_uppercase(char *text)
{
  assert(text);
  char *c = text, *end = strchr(text, '\0');
  for (; c != end; ++c)
    {
      *c = toupper(*c);
    }
  return text;
}

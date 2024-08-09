#include "adh.h"
#include <ctype.h>
#include <errno.h>
#include <string.h>

/*!
   \brief Read an integer datum value from a string
   
   "Read an integer datum value from a ASCII file line. This version calls the
   read/query version, passing NULL for query."  
   
   \param[in] data The piece of the line that is currently being parsed
   \return the value of the integer datum

 */   
int read_int_field(SIO io, char **data)
{
  return read_query_int_field(io, data, NULL);
}

/*!
   \brief Read an integer datum value from a string
   
   "Read an integer datum value from a ASCII file line."
   
   \param[in] data The piece of the line that is currently being parsed
   \param[out] query The optional flag which denotes whether a valid integer
   was found and read (READ_SUCCESS = successful, READ_FAIL = failed, or
   READ_NONE = no input field). Using this flag silences error messages and
   avoids automatic program termination, allowing the querying of input data
   (is the next input an integer?); pass NULL for regular operation. data is
   only incremented when successful.
   \return the value of the integer datum

 */
int read_query_int_field(SIO io, char **data, int *query)
{
  char *end = NULL;             /* pointer for error trapping */ 
  int token_len;                /* length of token (data input field) */
  int val = 0;                  /* holds the val */

  /* default flag to fail */
  if (query)
    {
      *query = READ_FAIL;
    }

  /* get length of datum input field */
  strip_white(data);
  token_len = strcspn(*data, REJECT COMMENT_CHARACTERS);
  if (token_len < 1)
    {
      if (query)
        {
          *query = READ_NONE;
          return val;
        }
      io_read_error(io, "Did not find expected integer input field.", TRUE);
    }

#ifdef IO_VALIDATION
  /* This provides additional data validation than that built into the strtol or
   * atoi functions. strtol and atoi read "123abc" as "123", and strtol leaves "abc"
   * for later, possibly delaying any error if not purposefully checked
   */
  char *datum = NULL, *sign = NULL;

  /* copy datum into seperate buffer */
  datum = (char *) malloc(sizeof(char) * (token_len + 1));
  strncpy(datum, *data, token_len);
  datum[token_len] = '\0';

  /* check that the input can be interpreted as an integer:
   *   - consists of only decimal digits ("1234567890") and sign ("-+") characters
   *   - must include at least one digit
   *   - sign is optional (only one instance allowed) and must precede all digits
   */
  sign = strpbrk(datum, "-+");  /* strpbrk returns a pointer to 1st instance or NULL */
  if (strspn(datum, "-+1234567890") != token_len || (sign && (sign != datum || token_len == 1 || strpbrk(sign + 1, "-+"))))
    {
      if (query)
        {
          return val;
        }
      char msg[MAXLINE];
      sprintf(msg, "Encountered an invalid integer input: %s", datum);
      io_read_error(io, msg, TRUE);
    }
  free(datum);
#endif

  /* strtol is preferred over atoi because it provides a pointer to after the
   * read token. The following ignores any immediate invalid characters thus
   * preventing any error trapping (bad practice).
   *
   *  val = atoi(*data);
   *  (*data) += token_len + 1;
   */

  /* read integer datum and trap for errors */
  errno = 0;
  val = strtol(*data, &end, 10); /* reads base 10 integers */
  if (errno || *data == end)
    {
      if (query)
        {
          return val;
        }
      char msg[MAXLINE];
      sprintf(msg, "Encountered an invalid integer input: %s", *data);
      io_read_error(io, msg, TRUE);
    }
  /* reset data to point to first character after read token */  
  *data = end;

  /* set flag to success */
  if (query)
    {
      *query = READ_SUCCESS;
    }
  
  return val;
}

/*!
   \brief Read an integer datum value from a string and provide a custom error
   message, if necessary.
   
   "Read an integer datum value from a ASCII file line. This version provides
   custom messages if the value is not found or the input is invalid (wrong
   type)."
   
   \param[in] data The piece of the line that is currently being parsed
   \param[out] query The optional flag to return result
   \param[in] field_name The name of the parameter field
   \param[in] field_idx The parameter field index, or UNSET_INT to ignore
   \param[in] abort Flag to abort if an error is encountered (TRUE/FALSE)

   \return the value of the integer datum

 */   
int read_int_field_custom(SIO io, char **data, int *query, char *field_name,
  int field_idx, int abort)
{
  int val, res, *pquery = NULL;
  
  if (query)
    {
      pquery = query;
    }
  else
    {
      pquery = &res;
    }

  val = read_query_int_field(io, data, pquery);

  if (*pquery == READ_SUCCESS)
    {
      return val;
    }

  char msg[MAXLINE], pos[10], pos_msg[25] = "";

  if (field_idx != UNSET_INT)
    {
      get_index_as_text(field_idx, pos);
      sprintf(pos_msg, " (%s parameter)", pos);
    }

  switch (*pquery)
    {
    case READ_NONE:
      sprintf(msg, "Did not find the %s integer input field%s.",
        field_name, pos_msg);
      break;
    case READ_FAIL:
      sprintf(msg, "Encountered an invalid integer input in the %s field%s.",
        field_name, pos_msg);
      break;
    default:
      assert(0); /* unhandled case */
      break;
    }
  io_read_error(io, msg, abort);

  return val;
}

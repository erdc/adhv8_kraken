#include "global_header.h"
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <errno.h>
#include <string.h>

/*!
   \brief Read a double datum value from a string
   
   "Read a double datum value from a ASCII file line. This version calls the
   read/query version, passing NULL for query."  
   
   \param[in] data The piece of the line that is currently being parsed
   \return the value of the double datum

 */   
double read_dbl_field(SIO io, char **data)
{
  return read_query_dbl_field(io, data, NULL);
}

/*!
   \brief Read a double datum value from a string
   
   "Read a double datum value from a ASCII file line."
   
   \param[in] data The piece of the line that is currently being parsed
   was found and read (READ_SUCCESS = successful, READ_FAIL = failed, or
   READ_NONE = no input field). Using this flag silences error messages and
   avoids automatic program termination, allowing the querying of input data
   (is the next input a double?); pass NULL for regular operation. data is
   only incremented when successful.

 */
double read_query_dbl_field(SIO io, char **data, int *query)
{
  char *end = NULL;             /* pointer for error trapping */ 
  int token_len;                /* length of token (datum input field) */
  double val = 0.0;             /* holds the value */

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
      io_read_error(io,"Did not find the expected double input field.", TRUE);
    }

#ifdef IO_VALIDATION
  /* This provides additional and clearer data validation than that built into the
   * strtod function called below.
   * strtod, for example:
   *   - allows NAN and INFINITY (case-insensitive)
   *   - updates errno only for overflow and underflow (and not for invalid input)
   *   - reads "3.0x10" as "3.0" and leaves "x10" for later, possibly delaying any
   *     error if not purposefully checked
   */
  char *datum = NULL, *sign = NULL, *digit = NULL, *decimal = NULL, *e = NULL;

  /* copy datum into seperate buffer */
  datum = (char *) malloc(sizeof(char) * (token_len + 1));
  strncpy(datum, *data, token_len);
  datum[token_len] = '\0';

  /* check that the input can be interpreted as a double:
   *   - consists of only decimal digits (1234567890), demical point (.),
   *     expontential notation (Ee), and sign (-+) characters.
   *   - must include at least one digit.
   *   - sign is optional with two instances permitted: (1) preceding all characters
   *     and (2) directly after expontential notation (preceding all expontential digits).
   *   - decimal point is optional with only one instance permitted preceding the
   *     exponential notation.
   *   - exponential notation is optional with only one instance permitted preceded and
   *     followed by a digit.
   * this does not distingush an integer from a double (integers are a subset of reals).
   * valid: "[-+]D.D{[Ee][-+]D}" (where D is at least one decimal digit)
   * invalid: ".", "-", "+", "e", "-.e", "De", "eD", ...
   */
  sign = strpbrk(datum, "-+");  /* strpbrk returns a pointer to 1st instance or NULL */
  digit = strpbrk(datum, "1234567890");
  decimal = strpbrk(datum, ".");
  e = strpbrk(datum, "Ee");
  if (strspn(datum, "-+1234567890.Ee") != token_len || !digit ||
      /* check decimal */
      (decimal && strpbrk(decimal + 1, ".")) ||
      /* check exponential notation */
      (e && (e < digit || e < decimal || strpbrk(e + 1, "Ee") || !isdigit(datum[token_len - 1]))))
    {
      if (query)
        {
          return val;
        }
      char msg[MAXLINE];
      sprintf(msg, "Encountered an invalid double input: %s", datum);
      io_read_error(io,msg, TRUE);
    }
  /* check sign */
  while (sign)
    {
      if ((!e && sign != datum) || (e && sign != datum && sign != e + 1))
        {
          if (query)
            {
              return val;
            }
          char msg[MAXLINE];
          sprintf(msg, "Encountered an invalid double input: %s", datum);
          io_read_error(io,msg, TRUE);
        }
      sign = strpbrk(sign + 1, "-+");
    }
  free(datum);
#endif

  /* read double datum and trap for errors */
  errno = 0;
  val = strtod(*data, &end);
  if (errno || *data == end)
    {
      /* try to read long double instead */
      /* strtold is prefered over sscanf and is included in C99, however PGI does not support it */
#ifndef __PGI
      errno = 0;

      long double ldval = strtold(*data, &end);

      if (errno || *data == end)
#else
      long double ldval;
      char format[10];

      /* create the format string which tells sscanf how many characters to read
         and interpret as a long double */
      sprintf(format, "%%%dLf", token_len);
      
      if (sscanf(*data, &format[0], &ldval) != 1)
#endif
        {
          /* still an error */
          if (query)
            {
              return val;
            }
          char msg[MAXLINE];
          sprintf(msg, "Encountered an invalid double input: %s", *data);
          io_read_error(io,msg, TRUE);
        }

      /* need to convert from a long double to a regular double */
      val = ldval;
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
   \brief Read a double datum value from a string and provide a custom error
   message, if necessary.
   
   "Read a double datum value from a ASCII file line. This version provides
   custom messages if the value is not found or the input is invalid (wrong
   type)."
   
   \param[in] data The piece of the line that is currently being parsed
   \param[out] query The optional flag to return result
   \param[in] field_name The name of the parameter field
   \param[in] field_idx The parameter field index, or UNSET_INT to ignore
   \param[in] abort Flag to abort if an error is encountered (TRUE/FALSE)

   \return the value of the double datum

 */   
double read_dbl_field_custom(SIO io, char **data, int *query, char *field_name,
  int field_idx, int abort)
{
  double val;
  int res, *pquery = NULL;

  if (query)
    {
      pquery = query;
    }
  else
    {
      pquery = &res;
    }
  
  val = read_query_dbl_field(io, data, pquery);

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
      sprintf(msg, "Did not find the %s double input field%s.", field_name, pos_msg);
      break;
    case READ_FAIL:
      sprintf(msg, "Encountered an invalid double input in the %s field%s.", field_name, pos_msg);
      break;
    default:
      assert(0); /* unhandled case */
      break;
    }
  io_read_error(io, msg, abort);

  return val;
}

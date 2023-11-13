#include "global_header.h"
#include <ctype.h>
#include <string.h>

/*!
   \brief Read a text field form a string
   
   "Read a text field from a ASCII file line. This version calls the
   read/query version, passing NULL for query."  
   
   \param[in] data The piece of the line that is currently being parsed
   \param[out] field The character string to copy the text field into
   \param[in] buffer_size The allocated size of field character array
   \return The text field (same pointer as input parameter) 

 */
char * read_text_field(SIO io, char **data, char *field, int buffer_size)
{
  return read_query_text_field(io, data, field, buffer_size, NULL);
}

/*!
   \brief Read a text field form a string
   
   "Read a text field from a ASCII file line. A single field is delimited by
   white space unless surrounded by double quotes ("") or single quotes ('')."
   
   \param[in] data The piece of the line that is currently being parsed
   \param[out] field The character string to copy the text field into
   \param[in] buffer_size The allocated size of field character array
   \param[out] query The optional flag which denotes whether a valid text field
   was found and read (READ_SUCCESS = successful, READ_FAIL = failed, or
   READ_NONE = no input field). Using this flag silences error messages and
   avoids automatic program termination, allowing the querying of input data
   (is there another field?); pass NULL for regular operation. data is
   only incremented when successful.
   \return The text field (same pointer as input parameter) 

 */
char * read_query_text_field(SIO io, char **data, char *field, int buffer_size, int *query)
{
  char quote[2] = "\"";
  int idx, nchar = 0, end_quote = FALSE;

  assert(*data);
  assert(field);
  assert(buffer_size > 1);
  if (query)
    {
      *query = READ_SUCCESS;
    }
  strip_white(data);
  if (!strlen(*data))
    {
      field[0] = '\0';
      if (query)
        {
          *query = READ_NONE;
          return field;
        }
      io_read_error(io, "Did not find expected text input field.", TRUE);
    }
  
  /* handle quote ("" or '') delimited text field */
  if (**data == '\"' || **data == '\'')
    {
      /* copy the quote character */
      strncpy(quote, *data, 1);
      (*data)++;
      /* remove any leading white space */ 
      while (isspace(**data))
        {
          (*data)++;
        }
      nchar = strcspn(*data, quote);
      /* check whether an ending quote was found */
      if ( (int)strlen(*data) != nchar)
        {
          end_quote = TRUE;
        }
      else if (query)
        {
          *query = READ_FAIL;
        }
      else
        {
          io_read_error(io, "Text input field missing end quote.", TRUE);
        }
    }
  else
    {
      nchar = strcspn(*data, " \t");
    }

  /* copy text field */
  if (nchar < buffer_size)
    {
      /* do not copy the quotes */
      strncpy(field, *data, nchar);
      field[nchar] = '\0';
      idx = nchar - 1;
    }
  else
    {
      /* copy only the as much text as fits */
      strncpy(field, *data, buffer_size);
      field[buffer_size - 1] = '\0';
      idx = buffer_size - 2;
      if (query)
        {
          *query = READ_FAIL;
        }
      else
        {
          char msg[MAXLINE];

          sprintf(msg, "Test input field exceeds the allowed length of %d "
            "characters.", buffer_size);
          io_read_error(io, msg, TRUE);
        }
    }

  /* remove any ending white space including the '\n' */
  while (isspace(field[idx]) && idx >= 0)
    {
      field[idx] = '\0';
      --idx;
    }

  /* increment the data pass the text field */
  for (idx = (end_quote ? -1 : 0); idx < nchar; ++idx)
    {
      (*data)++;
    }

  return field;
}

/*!
   \brief Read a text field from a string and provide a custom error
    message, if necessary.

   "Read a text field from a ASCII file line. This version provides custom
   messages if the value is not found or the input is invalid (no closing
   quote or the buffer is too small)."

   \param[in] data The piece of the line that is currently being parsed
   \param[out] field The character string to copy the text field into
   \param[in] buffer_size The allocated size of field character array
   \param[out] query The optional flag to return result
   \param[in] field_name The name of the parameter field
   \param[in] field_idx The parameter field index, or UNSET_INT to ignore
   \param[in] abort Flag to abort if an error is encountered (TRUE/FALSE)
   \return The text field (same pointer as input parameter)

*/
char * read_text_field_custom(SIO io, char **data, char *field, int buffer_size,
  int *query, char *field_name, int field_idx, int abort)
{
  int res, *pquery = NULL;

  if (query)
    {
      pquery = query;
    }
  else
    {
      pquery = &res;
    }

  field = read_query_text_field(io, data, field, buffer_size, pquery);

  if (*pquery == READ_SUCCESS)
    {
      return field;
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
        sprintf(msg, "Did not find the %s text input field%s.",
          field_name, pos_msg);
        break;
      case READ_FAIL:
        sprintf(msg, "Encountered an error reading the %s text field%s, "
          "missing end quote or buffer size too small.", field_name,
          pos_msg);
        break;
      default:
        assert(0); /* unhandled case */
        break;
    }
  io_read_error(io, msg, abort);

  return field;
}

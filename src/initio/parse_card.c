#include <global_header.h>
#include <ctype.h>

/* Currently all input card tokens are strings of length <= 4. However, this function
   will work if the the card length is increased in the future */

/* input is the string to be parsed; data will be pointer to the remainder of 
   the line */

/* Cards must begin with an alphabetic character so they can be distinguished from
   numeric input */

CARD parse_card(
  char *string, /* to be parsed */
  char **data   /* pointer to the remainder of the line (after the card) */
)
{
  unsigned int i;
  char *ptr = NULL;
  int card[MAX_CARD_LENGTH];
#if MAX_CARD_LENGTH <= 3
  unsigned int card_value, constant;
#elif MAX_CARD_LENGTH <= 6
  unsigned long card_value, constant;
#endif
  size_t token_length;

  /* strip off leading white space */
  strip_white(&string);

  /* if empty or the first character is not alphabetic then return no card */
  token_length = strcspn(string, REJECT);
  if (token_length < 1 || !isalpha(*string))
    {
      return NO_CARD;
    }

  /* ignore (do not include)  any characters beyond the specified card length */
  /* for example, "OBJTYPE" from a hot start file will be read as "OBJTYP", leaving 
     the "E" off */
  token_length = MIN(token_length, MAX_CARD_LENGTH);

  /* initialize */
  for (i = 0; i < MAX_CARD_LENGTH; ++i)
    {
      card[i] = 0;
    }
  ptr = string;


  /* Calculate integer value equivalence of file card using base 37, where individual
     characters A-Z and a-z are converted to 1-26 and 0-9 are converted to 27-36 and
     combined:
 
     card_value = card[0] + (37 * card[1]) + (37^2 * card[2]) + ...

     For example: "XY1" converts to:

            'X'           'Y'           '1'
     (37^0 * 24) + (37^1 * 25) + (37^2 * 28) = 24 + (37 * 25) + (37^2 * 28) = 39281
      
     This formulation allows for any length card.
   */

  /* Inspect each character of the card and convert into a base 37 value */
  for (i = 0; i < token_length; ++i, ++ptr)
    {
      /* bail out if there is a non-alphanumeric character */
      if (!isalnum(*ptr))
        {
          return NO_CARD;
        }
      /* convert characters to uppercase and save ASCII value */
      card[i] = toupper(*ptr);
      /* convert from ASCII to base 37 */
      if (isalpha(card[i]))
        {
          card[i] = card[i] - 'A' + 1;
        }
      else
        {
          card[i] = card[i] - '0' + 27;
        }
    }

  card_value = 0;
  constant = 1;
  for (i = 0; i < token_length; ++i)
    {
      card_value += constant * card[i];
      constant *= 37;
    }

  /* if pointer is not NULL then set it to point to after the card */
  if (data)
    {
      (*data) = string + token_length;
      strip_white(data);
    }

  return (CARD)card_value;
}

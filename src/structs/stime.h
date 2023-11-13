/*!
   \file type_time.h
   \brief Data Structures Needed by Linear Solver Routine
 */

#ifdef TIME_INIT
#define EXTERN_TIME
#else
#define EXTERN_TIME extern
#endif

/*!
   \struct Time_Info
   \brief General Parameters for Time Stepping
 */
typedef struct
{
  double dt_reduce_factor;	/* Multiplier to Reduce Time Step Size */
  double dt_increase_factor;	/* Multiplier to Increase Time Step Size */
  double dt;			/* Current Proposed Time Step */
  double dt_prev;		/* Previous Accepted Time Step */
  double dt_init;		/* Initial Time Step */
  double dt_min;		/* Minimum Time Step (Thus Far) */
  double dt_max;		/* Maximum Time Step (Thus Far) */
  double t_init;		/* Initial Time */
  double t_prev;		/* Previous Accepted Time */
  double t_final;		/* Final Time */
  double Percent_Done;		/* Percent Done */
  int linearSteps;		/* Previous Linear Steps */
  int nonlinearSteps;		/* Previous Non-Linear Steps */
  int Schema;			/* Schema for Changing Time Steps */
} Time_Info;

/************************************************************/

EXTERN_TIME Time_Info time_info;

/************************************************************/

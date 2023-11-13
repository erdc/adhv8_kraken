C******************************************************************************
C PADCIRC VERSION 45.12 03/17/2006                                            *
C  last changes in this file VERSION 45.08                                    *
C... TCM V50.66.01 --ADDING TIME DEPENDENT BATHYMETRY                         *
C                 ADDED NDDT,BTIMINC,DP1,DP2,BTIME1,BTIME2                    *
C                                                                             *
C This module declares all global variables that are not exclusive to the 3D  *
C routines.  The 3D exclusive variables are declared in global_3DVS           *
C                                                                             *
C******************************************************************************

C
      MODULE GLOBAL
      USE SIZES
      !USE VERSION
      USE KDTREE2_MODULE   !v49.48.01 tcm fast search algorithm
      IMPLICIT NONE
      SAVE
C...
C...SET GLOBAL PHYSICAL CONSTANTS
C...
      REAL(sz) :: g = 9.80665 !default gravitational acceleration

C.....nominal density of water RHOWAT0
      REAL(SZ), PARAMETER ::  RhoWat0=1000.D0

C.....Background Atmospheric Pressure in mb
      REAL(SZ), PARAMETER :: PRBCKGRND = 1013.0d0   !tcm 20100617 v49.16 Added
      REAL(8) :: rhoAir = 1.15d0         ! kg/m^3
      REAL(8) :: windReduction = 0.9d0
      REAL(8) :: one2ten = 0.8928d0

C.....PI etc
      REAL(8), PARAMETER ::  PI=3.141592653589793D0
      REAL(8), PARAMETER ::  DEG2RAD = PI/180.D0  ! degrees to radians
      REAL(8), PARAMETER ::  RAD2DEG = 180.D0/PI  ! radians to degrees
      REAL(8), PARAMETER ::  e=2.718281828459045d0


      ! jgf50.38.03: Additional useful constants
      REAL(8), PARAMETER :: Rearth = 6378206.4d0 ! radius of earth (m)
      REAL(8), PARAMETER :: nm2m = 1852.d0 ! nautical miles to meters
      REAL(8), PARAMETER :: m2nm = 1.d0/nm2m ! meters to nautical miles
      REAL(8), PARAMETER :: kt2ms = nm2m / 3600.0 ! knots to m/s
      REAL(8), PARAMETER :: ms2kt = 1.d0 / kt2ms ! m/s to knots
      REAL(8), PARAMETER :: omega = 2.0d0*pi / 86164.2d0
      REAL(8), PARAMETER :: mb2pa = 100.d0


C... v49.48 tcm -- added for coupling with STWAVE (NRS=4)
      LOGICAL  :: CPL2STWAVE = .FALSE.  !Coupled to STWAVE (NRS=4)
      LOGICAL  :: Flag_ElevError = .FALSE.  

C      real(sz) rampriver,drampriver,rampriver1,rampriver2

      INTEGER :: screenUnit = 6 ! jgf46.19 I/O unit where screen output is sent

      REAL(SZ),ALLOCATABLE, TARGET ::  DP(:)  !TCM V50.66.01 FOR TIME VARYING BATHYMETRY


      REAL(SZ),ALLOCATABLE ::   WVNX1(:),WVNY1(:),PRN1(:)
      REAL(SZ),ALLOCATABLE ::   WVNX2(:),WVNY2(:),PRN2(:)

      INTEGER NP,IHOT

C...
C...DECLARE REAL(8) AND CHAR VARIABLES, EQUIVALENCES
C...
      REAL(8) SLAM0,SFEA0
      REAL(8),ALLOCATABLE :: SLAM(:),SFEA(:),X(:),Y(:)


C...
      INTEGER  NWS

      INTEGER :: NABOUT = 0
      INTEGER :: NSCREEN = 1
      ! COMMUNICATORS FOR GLOBAL FILE WRITING
      INTEGER, ALLOCATABLE ::  COMM_WRITER(:)
      ! COMMUNICATORS FOR HOTSTART FILE WRITING !st3
      INTEGER, ALLOCATABLE ::  COMM_WRITEH(:), COMM_HSLEEP(:)

      INTEGER :: SIG_TERM  = 999  ! signal to terminate adc_writer
!  tcm v49.48.01
C     End of variables shared between messenger and writer modules
C
C     jgf49.44: Log levels, in order from largest amount of log messages
C     written (DEBUG) to fewest log messages written (ERROR). Compared
C     to the value of NABOUT when determining which messages to write
C     to the screen or to log files.
      INTEGER, PARAMETER :: DEBUG = -1  ! write all messages and echo input
      INTEGER, PARAMETER :: ECHO = 0    ! echo input, plus write all non-debug
      INTEGER, PARAMETER :: INFO = 1    ! don't echo input; write all non-debug
      INTEGER, PARAMETER :: WARNING = 2 ! don't echo input; write only warn/err
      INTEGER, PARAMETER :: ERROR = 3   ! don't echo input; only fatal msgs
C
      CHARACTER(len=10), dimension(5) :: logLevelNames
      CHARACTER(len=50), dimension(50) :: messageSources ! subroutine names
      CHARACTER(len=1024) :: scratchMessage ! used for formatted messages
      INTEGER :: sourceNumber ! index into messageSources for current sub


C-------------------end of data declarations----------------------------------C

      CONTAINS

C******************************************************************************
C                                                                             *
C    Transform from lon,lat (lamda,phi) coordinates into CPP coordinates.     *
C    Lon,Lat must be in radians.                                              *
C                                                                             *
C******************************************************************************
      SUBROUTINE CPP(X,Y,RLAMBDA,PHI,RLAMBDA0,PHI0)
      IMPLICIT NONE
      REAL*8 X,Y,RLAMBDA,PHI,RLAMBDA0,PHI0,R
      R=6378206.4d0
      X=R*(RLAMBDA-RLAMBDA0)*COS(PHI0)
      Y=PHI*R
      RETURN
      END SUBROUTINE CPP


C--------------------------------------------------------------------
C     S U B R O U T I N E    O P E N   L O G   F I L E
C--------------------------------------------------------------------
C     jgf50.65: Open the log file; this must be called after make dirname
C     so that we know where to put this log file.
C--------------------------------------------------------------------
      SUBROUTINE openLogFile()
      IMPLICIT NONE
C...
C...  OPEN STATEMENT FOR UNIT 16 OUTPUT FILE (ADCIRC LOG FILE)
C...
#ifdef CMPI
      OPEN(16,FILE=trim(localdir)//'/'//'fort.16', ACTION='WRITE',
     &    STATUS='REPLACE')
#else
      OPEN(16,FILE='fort.16', ACTION='WRITE', STATUS='REPLACE')
#endif
C--------------------------------------------------------------------
      END SUBROUTINE openLogFile
C--------------------------------------------------------------------

C--------------------------------------------------------------------
C     S U B R O U T I N E    S C R E E N
C--------------------------------------------------------------------
C     jgf49.44: General purpose subroutine to write a message to
C     the screen with a certain "logging level", and subject to the
C     user's selection of where to write screen output. The logging
C     level is controlled by NABOUT from the fort.15 file. The actual
C     destination of messages written to the screen is controlled by
C     NSCREEN from the fort.15 file.
C
C     In parallel, only the processor with rank 0 actually writes
C     the message.
C
C     This subroutine assumes that the global variable "caller" has
C     been set to the name of the subroutine calling it. Therefore,
C     the setMessageSource subroutine must be called at the beginning
C     of the subroutine that calls this one, and unsetMessageSource
C     must be called at the end.
C--------------------------------------------------------------------
      SUBROUTINE screenMessage(level, message)
      IMPLICIT NONE
      INTEGER, intent(in) :: level
      CHARACTER(*), intent(in) :: message
      INTEGER j  ! loop counter for stack
C
      IF (myProc.eq.0) THEN
         IF (NSCREEN.ne.0) THEN
            IF (level.ge.NABOUT) THEN
#ifdef FULL_STACK
               write(screenUnit,331,advance="no")
     7             trim(logLevelNames(level+2)),
     &             (trim(messageSources(j)),j=1,sourceNumber)
               write(screenUnit,332) trim(message)
#else
               write(screenUnit,333) trim(logLevelNames(level+2)),
     &             trim(messageSources(sourceNumber)), trim(message)
#endif
#ifdef FLUSH_MESSAGES
               flush(screenUnit)
#endif
            ENDIF
         ENDIF
      ENDIF
331   FORMAT(A,": ",A,50(:,"->",A))
332   FORMAT(": ",A)
333   FORMAT(A,": ",A,": ",A)
C--------------------------------------------------------------------
      END SUBROUTINE screenMessage
C--------------------------------------------------------------------


C--------------------------------------------------------------------
C     S U B R O U T I N E    L O G   M E S S A G E
C--------------------------------------------------------------------
C     jgf49.44: General purpose subroutine to write a message to
C     the fort.16 file. In parallel, processors of all ranks will
C     write the message to their own subdomain fort.16 files.
C
C     This subroutine assumes that the global variable "caller" has
C     been set to the name of the subroutine calling it. Therefore,
C     the setMessageSource subroutine must be called at the beginning
C     of the subroutine that calls this one, and unsetMessageSource
C     must be called at the end.
C--------------------------------------------------------------------
      SUBROUTINE logMessage(level, message)
      IMPLICIT NONE
      INTEGER, intent(in) :: level
      CHARACTER(*), intent(in) :: message
      INTEGER j  ! loop counter for stack
C
      IF (level.ge.NABOUT) THEN
#ifdef FULL_STACK
         write(16,331,advance="no") trim(logLevelNames(level+2)),
     &       (trim(messageSources(j)),j=1,sourceNumber)
         write(16,332) trim(message)
#else
         write(16,333) trim(logLevelNames(level+2)),
     &       trim(messageSources(sourceNumber)), trim(message)
#endif
#ifdef FLUSH_MESSAGES
               flush(16)
#endif
      ENDIF
331   FORMAT(A,": ",A,50(:,"->",A))
332   FORMAT(": ",A)
333   FORMAT(A,": ",A,": ",A)
C--------------------------------------------------------------------
      END SUBROUTINE logMessage
C--------------------------------------------------------------------

C--------------------------------------------------------------------
C     S U B R O U T I N E   A L L    M E S S A G E
C--------------------------------------------------------------------
C     jgf49.44: General purpose subroutine to write a message to
C     both the screen and to the fort.16 log file.
C--------------------------------------------------------------------
      SUBROUTINE allMessage(level, message)
      IMPLICIT NONE
      INTEGER, intent(in) :: level
      CHARACTER(*), intent(in) :: message
C
      call logMessage(level, message)
      call screenMessage(level, message)
C--------------------------------------------------------------------
      END SUBROUTINE allMessage
C--------------------------------------------------------------------

C--------------------------------------------------------------------
C     S U B R O U T I N E   S E T   M E S S A G E   S O U R C E
C--------------------------------------------------------------------
C     jgf49.44: Sets the name of the subroutine that is writing
C     log and/or screen messages. Must use at the start of any subroutine
C     that calls screen, logMessage, or allMessage.
C--------------------------------------------------------------------
      SUBROUTINE setMessageSource(source)
      IMPLICIT NONE
      CHARACTER(*), intent(in) :: source
C
      sourceNumber = sourceNumber + 1
      messageSources(sourceNumber) = source
C--------------------------------------------------------------------
      END SUBROUTINE setMessageSource
C--------------------------------------------------------------------

C--------------------------------------------------------------------
C     S U B R O U T I N E   U N S E T   M E S S A G E   S O U R C E
C--------------------------------------------------------------------
C     jgf49.44: Removes the name of the subroutine that is no longer
C     writing log and/or screen messages. Must use at the end of
C     any subroutine that calls screen, logMessage, or allMessage.
C--------------------------------------------------------------------
      SUBROUTINE unsetMessageSource()
      IMPLICIT NONE
C
      sourceNumber = sourceNumber - 1
C--------------------------------------------------------------------
      END SUBROUTINE unsetMessageSource
C--------------------------------------------------------------------


C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E   O P E N  F I L E  F O R  R E A D
C-----------------------------------------------------------------------
C     jgf50.16 Added general subroutine for opening an existing input
C     file for reading. Includes error checking.
C-----------------------------------------------------------------------
      SUBROUTINE openFileForRead(lun, filename, errorIO)
      IMPLICIT NONE
      INTEGER, intent(in) :: lun   ! fortran logical unit number
      CHARACTER(*), intent(in) :: filename ! full pathname of file
      INTEGER, intent(out) :: errorIO  ! zero if the file opened successfully
      LOGICAL :: fileFound    ! .true. if the file is present
C
      call setMessageSource("openFileForRead")
#if defined(GLOBAL_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C
      errorIO = 0
C
C     Check to see if file exists
      write(scratchMessage,21) lun
  21  format("Searching for file to open on unit ",I5,"...")
      call logMessage(INFO,trim(scratchMessage))
      inquire(FILE=trim(filename),EXIST=fileFound)
      if (fileFound.eqv..false.) then
         write(scratchMessage,23) trim(filename)
  23     format("The file '",A,"' was not found.")
         call allMessage(ERROR,scratchMessage)
         errorIO = 1
#if defined(GLOBAL_TRACE) || defined(ALL_TRACE)
         call allMessage(DEBUG,"Return.")
#endif
         call unsetMessageSource()
         return  ! file not found
      else
         write(scratchMessage,24) trim(filename)
  24     format("The file '",A,"' was found. The file will be opened.")
         call logMessage(INFO,trim(scratchMessage))
      endif
C
C     Open existing file
      OPEN(lun,FILE=trim(filename),STATUS='OLD',
     & ACTION='READ',IOSTAT=errorIO)
      if (errorIO.ne.0) then
          write(scratchMessage,25) trim(filename)
  25      format("Could not open the file '",A,"'.")
          call allMessage(ERROR,trim(scratchMessage))
#if defined(GLOBAL_TRACE) || defined(ALL_TRACE)
          call allMessage(DEBUG,"Return.")
#endif
          call unsetMessageSource()
          return  ! file found but could not be opened
      else
         write(scratchMessage,26) trim(filename)
  26     format("The file '",A,"' was opened successfully.")
         call logMessage(INFO,trim(scratchMessage))
      endif

#if defined(GLOBAL_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      return
C-----------------------------------------------------------------------
      END SUBROUTINE openFileForRead
C-----------------------------------------------------------------------

C     ----------------------------------------------------------------
C      F U N C T I O N   S P H E R I C A L   D I S T A N C E
C     ----------------------------------------------------------------
C     jgf49.1001 Function to get the distance along the surface of
C     a sphere (the earth's surface in this case).
C     ----------------------------------------------------------------
      REAL(SZ) FUNCTION sphericalDistance(dx, dy, y1, y2)
      IMPLICIT NONE
      REAL(SZ), intent(in) :: DX    ! longitude distance in radians
      REAL(SZ), intent(in) :: DY    ! latitude distance in radians
      REAL(SZ), intent(in) :: y1    ! degrees latitude of starting point
      REAL(SZ), intent(in) :: y2    ! degrees latitude of ending point
C
C     compute the distances based on haversine formula for
C     distance along a sphere
      sphericalDistance = Rearth *
     &  ( 2.0d0*ASIN(
     &     sqrt(sin(DY/2.0d0)**2.0d0
     &      + cos(y1*DEG2RAD)*cos(y2*DEG2RAD)*sin(DX/2.0d0)**2.0d0
     &   )))

      RETURN
C     ----------------------------------------------------------------
      END FUNCTION sphericalDistance
C     ----------------------------------------------------------------
      !=================================================================
      END MODULE GLOBAL
C     ------------------------------------------------------------------
C     ------------------------------------------------------------------

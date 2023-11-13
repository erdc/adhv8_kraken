C******************************************************************************
C PADCIRC VERSION 45.12 03/17/2006                                            *
C  last changes in this file VERSION 45.56                                    *
C                                                                             *
C******************************************************************************
C
      MODULE SIZES
      IMPLICIT NONE
C
C...SET NUMBER OF BYTES "SZ" IN REAL(SZ) DECLARATIONS
C...SET "NBYTE" FOR PROCESSING INPUT DATA RECORD LENGTH

#ifdef REAL4
      INTEGER, PARAMETER :: SZ = 4
      !INTEGER, PARAMETER :: NBYTE=4
#endif
#ifdef REAL8
      INTEGER, PARAMETER :: SZ = 8
      !INTEGER, PARAMETER :: NBYTE=8
#endif


      CHARACTER(2048), TARGET, SAVE :: GLOBALDIR, LOCALDIR
      INTEGER,SAVE :: MYPROC
      CHARACTER(2048), SAVE :: GBLINPUTDIR
      INTEGER ::  MNE,MNP,MNEI,MNETA,MNVEL,MNBFR, MNFFR, MNTIF
     & ,MNSTAV, MNSTAE, MNBOU, MNSTAM, MNSTAC, MNOPE
      INTEGER :: MNPROC         ! Number of compute processors
      INTEGER :: MNWPROC        ! Number of writer processors
      INTEGER :: MNWPROH        ! Number of Hwriter processors !st3 100711
      LOGICAL      :: READ_LOCAL_HOT_START_FILES
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      END MODULE SIZES
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

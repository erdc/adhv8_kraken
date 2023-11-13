C******************************************************************************
C PADCIRC VERSION 47.xx 10/13/2006                                             *
C******************************************************************************
C
      MODULE MESSENGER
      USE SIZES
      USE GLOBAL, ONLY: COMM_WRITER, COMM_WRITEH, COMM_HSLEEP, 
     &    SIG_TERM,CPL2STWAVE,Flag_ElevError,DEBUG, 
     &    allMessage, setMessageSource,unsetMessageSource     
     
#ifdef HAVE_MPI_MOD
      USE MPI
#endif
      IMPLICIT NONE
      SAVE


C--------------------------------------------------------------------------
C  This module supplies the MPI Message-Passing Interface for PADCIRC.
C  Uses asynchronous communication with buffer packing as performance
C  enhancement for "cluster" architectures.
C--------------------------------------------------------------------------


C  Message-Passing Array space

      INTEGER ::  MPI_COMM_ADCIRC                     ! Local communicator
      INTEGER ::  IERR
      INTEGER ::  TAG = 100


      !jgf50.82: Create a flag for unrecoverable issue on a subdomain
      LOGICAL :: subdomainFatalError ! true if mpi_abort should be called
C
C
      CONTAINS

C---------------------end of data declarations--------------------------------C


C---------------------------------------------------------------------
C      S U B R O U T I N E   M S G _ F I N I 
C---------------------------------------------------------------------
C  Delete MPI resources and Shutdown MPI library.
C  vjp  8/29/1999
C---------------------------------------------------------------------
      SUBROUTINE MSG_FINI (NO_MPI_FINALIZE)
      IMPLICIT NONE
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
      LOGICAL, OPTIONAL :: NO_MPI_FINALIZE
      INTEGER I
C
      call setMessageSource("msg_fini")
#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.") 
#endif
C
      IF(MNWPROC > 0) THEN
        IF(MYPROC.eq.0) THEN
          DO I=1,MNWPROC
            WRITE(16,*)'PROC ',MYPROC,' IS SENDING SIG_TERM TO WRITER ',
     &           I
            CALL MPI_SEND(SIG_TERM,1,MPI_INTEGER,MNPROC,
     &           TAG,COMM_WRITER(I),IERR)
          ENDDO
        ENDIF
      ENDIF

      IF(MNWPROH > 0) THEN             !st3 100711 for hsfile
        IF(MYPROC.eq.0) THEN
          DO I=1,MNWPROH
            WRITE(16,*)'PROC ',MYPROC,' IS SENDING SIG_TERM TO HSWRITER',
     &           I
            CALL MPI_BARRIER(COMM_HSLEEP(I),IERR)
            CALL MPI_SEND(SIG_TERM,1,MPI_INTEGER,MNPROC,
     &           TAG,COMM_WRITEH(I),IERR)
          ENDDO
        ENDIF
      ENDIF
      
      IF (subdomainFatalError.eqv..true.) THEN
         ! jgf50.82: Return the rank of the offending processor
         ! as the error code        
         CALL MPI_ABORT(MPI_COMM_ADCIRC,MYPROC,IERR)
      ENDIF

      ! tcm v51.32  added a "go nuclear" option for killing
      ! all MPI processes when Elevation Greater than Error 
      ! elevation is exceeded, not just the mpi processes owned
      ! by ADCIRC when coupled with STWAVE via CSTORM-MS
      ! Note: mpi_comm_world is defined in mpich.f
      IF (CPL2STWAVE.eqv..true.) then
         IF (Flag_ElevError.eqv..true.) then
            call mpi_abort(mpi_comm_world,myproc,ierr)
         endif
      endif

      IF (PRESENT(NO_MPI_FINALIZE)) THEN
         IF (.NOT.NO_MPI_FINALIZE) THEN
            CALL MPI_FINALIZE(IERR)
            IF (MYPROC.EQ.0)
     &      PRINT *, "MPI terminated with Status = ",IERR
         ENDIF
      ELSE
         CALL MPI_FINALIZE(IERR)
         IF (MYPROC.EQ.0)
     &   PRINT *, "MPI terminated with Status = ",IERR
      ENDIF

#if defined(MESSENGER_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      RETURN
C---------------------------------------------------------------------      
      END SUBROUTINE MSG_FINI
C---------------------------------------------------------------------

      END MODULE MESSENGER

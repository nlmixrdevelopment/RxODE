      subroutine rprint(msg)
      character (len=*) msg
        call dblepr(msg, -1, 0, 0)
      end subroutine 

      subroutine rprintid(msg, i1, d1)
      character (len=*) msg
      double precision d1
      integer i1
        call dblepr(msg, -1, d1, 1)
        call intpr(" ", -1, i1, 1)
      end subroutine 

      subroutine rprintd1(msg, d1)
      character (len=*) msg
      double precision d1
        call dblepr(msg, -1, d1, 1)
      end subroutine 

      subroutine rprintd2(msg, d1, d2)
      character (len=*) msg
      double precision DBL(2), d1, d2
        DBL(1) = d1
        DBL(2) = d2
        call dblepr(msg, -1, DBL, 2)
      end subroutine 

      subroutine rprinti1(msg, i1)
      character (len=*) msg
      integer i1
        call intpr(msg, -1, i1, 1)
      end subroutine 

      subroutine rprinti2(msg, i1, i2)
      character (len=*) msg
      INTEGER IN(2), i1, i2
        IN(1) = i1
        IN(2) = i2
        call intpr(msg, -1, IN, 2)
      end subroutine 

      subroutine rprinti3(msg, i1, i2, i3)
      character (len=*) msg
      INTEGER IN(3), i1, i2, i3
        IN(1) = i1
        IN(2) = i2
        IN(3) = i3
        call intpr(msg, -1, IN, 3)
      end subroutine 

      subroutine rprint2(msg)
      implicit none
      character (len = *) msg
        call dblepr(msg, 61, 0, 0)
      end subroutine 


*DECK XERRWD
      SUBROUTINE XERRWD (MSG, NMES, NERR, LEVEL, NI, I1, I2, NR, R1, R2)

C***PURPOSE  Write error message with values.
C***original AUTHOR  Hindmarsh, Alan C., (LLNL)
C Rewritten to be used with R by Karline Soetaert
C
C  All arguments are input arguments.
C
C  MSG    = The message (character array).
C  NMES   = The length of MSG (number of characters).
C  NERR   = The error number (not used).
C  LEVEL  = The error level..
C           0 or 1 means recoverable (control returns to caller).
C           2 means fatal (run is aborted--see note below).
C  NI     = Number of integers (0, 1, or 2) to be printed with message.
C  I1,I2  = Integers to be printed, depending on NI.
C  NR     = Number of reals (0, 1, or 2) to be printed with message.
C  R1,R2  = Reals to be printed, depending on NR.
C
C-----------------------------------------------------------------------
C
C  Declare arguments.
C
      DOUBLE PRECISION R1, R2, RVEC(2), Dummy
      INTEGER I, NMES, NERR, LEVEL, NI, I1, I2, NR, Ivec(2)
      CHARACTER(LEN=*) MSG
      INTEGER LUNIT, IXSAV, MESFLG

      dummy = 0.d0
C      call dblepr(MSG, NMES, dummy, 0)
      MSG = MSG(1:NMES) // char(0)
      call rprintf(MSG)


      IF (NI .EQ. 1) THEN
C       call intpr('In above message, I = ', 22, I1, 1)
        MSG = 'In above message, I1 = %d' // char(0)
        call rprintfi1(MSG, I1)
        MSG = ' ' // char(0)
        call rprintf(MSG)
      ENDIF

      IF (NI .EQ. 2) THEN 
       IVEC(1) = I1
       IVEC(2) = I2
C       call intpr('In above message, I = ', 22, IVEC, 2)
        MSG = 'In above message, I1 = %d, I2 = %d' // char(0)
        call rprintfi2(MSG, I1, I2)
        MSG = ' ' // char(0)
        call rprintf(MSG)
      ENDIF

      IF (NR .EQ. 1) THEN
C       call dblepr('In above message, R = ', 22, R1, 1)
        MSG = 'In above message, R1 = %g' // char(0)
        call rprintfd1(MSG, R1)
        MSG = ' ' // char(0)
        call rprintf(MSG)
      ENDIF

      IF (NR .EQ. 2) THEN
       RVEC(1) = R1
       RVEC(2) = R2
C       call dblepr('In above message, R1 = ', 22, RVEC, 2)
        MSG = 'In above message, R1 = %g, R2 = %g' // char(0)
        call rprintfd2(MSG, R1, R2)
        MSG = ' ' // char(0)
        call rprintf(MSG)
      ENDIF

C  Abort the run if LEVEL = 2.
       if (LEVEL .eq. 2) call rexit ("fatal error")
 100  RETURN

      END

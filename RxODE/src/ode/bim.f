C -----------------------------------------------------------------------------------
C     THE CODE BIM NUMERICALLY SOLVES (STIFF) DIFFERENTIAL ODE 
C     PROBLEMS.
C
C     Copyright (C)2002-2007   
C
C     Authors: CECILIA MAGHERINI (cecilia.magherini@ing.unipi.it)
C              LUIGI   BRUGNANO  (brugnano@math.unifi.it) 
C
C
C     This program is free software; you can redistribute it and/or
C     modify it under the terms of the GNU General Public License
C     as published by the Free Software Foundation; either version 2
C     of the License, or (at your option) any later version.
C
C     This program is distributed in the hope that it will be useful,
C     but WITHOUT ANY WARRANTY; without even the implied warranty of
C     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C     GNU General Public License for more details.
C
C     Licensed under The GNU General Public License, Version 2 or later.
C       http://www.gnu.org/licenses/info/GPLv2orLater.html
C
C     You should have received a copy of the GNU General Public License
C     along with this program; if not, write to the Free Software
C     Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
C     USA.
C -----------------------------------------------------------------------------------


      SUBROUTINE BIM(M,FCN,T0,TEND,Y0,H,RTOL,ATOL,
     &               JAC,IJAC,MLJAC,MUJAC,
     &               WORK,LWORK,IWORK,LIWORK,
     &               RPAR,IPAR,IOUT,IDID)

C -----------------------------------------------------------------------------------
C -----------------------------------------------------------------------------------
C
C     PURPOSE:    BIM SOLVES A (STIFF) DIFFERENTIAL ODE PROBLEM,
C     --------
C				Y'    = F(T,Y),     T0<=T<=TEND,
C				Y(T0) = Y0,
C
C                 BY MEANS OF BLENDED IMPLICIT METHODS.
C                 BLENDED IMPLICIT METHODS ARE A CLASS OF BLOCK 
C                 METHODS	PROVIDING A (RELATIVELY) EASY DEFINITION 
C                 OF SUITABLE NONLINEAR SPLITTINGS FOR SOLVING THE
C                 CORRESPONDING DISCRETE PROBLEMS [1,4-6].
C                 THE CODE BIM IMPLEMENTS A VARIABLE STEPSIZE-
C                 VARIABLE ORDER METHOD. ORDERS: 4-6-8-10-12.
C                 IMPLEMENTATION DETAILS ARE IN REFERENCES [1-4].
C
C
C     AUTHORS:    L.BRUGNANO
C     --------    DIPARTIEMNTO DI MATEMATICA "U.DINI"
C                 VIALE MORGAGNI 67/A
C                 50134 FIRENZE
C                 ITALY
C                 E-MAIL: BRUGNANO@MATH.UNIFI.IT
C
C                 C.MAGHERINI, 
C                 DIPARTIMENTO DI MATEMATICA APPLICATA "U.DINI"
C                 VIA BUONARROTI, 1/C 
C                 56127 PISA
C                 ITALY
C                 E-MAIL: CECILIA.MAGHERINI@ING.UNIPI.IT
C
C
C
C     CODE HOME PAGE:   http://www.math.unifi.it/~brugnano/BiM/index.html
C     ---------------
C
C     CODE:       THE CODE IS MADE UP OF THREE FILES:
C     -----        - BIM.F      (I.E. THE PRESENT FILE) WHICH CONTAINS THE MAIN 
C                    INTEGRATION PROCEDURE
C                  - SUBBIM.F   CONTAINING ADDITIONAL AND LINEAR ALGEBRA
C                    PROCEDURES
C                  - PARAM.DAT  CONTAINING VARIOUS PARAMETERS USED IN THE CODE
C
C     CURRENT RELEASE:  2.0, APRIL 4, 2005.
C     ----------------
C
C     REL.HISTORY:  1.0    fall of 2002
C     ------------         - pre-release;
C                   1.1    June 18, 2003
C                          - first version released;
C                   1.1.1  July 20, 2003 
C                          - minor bugs fixed;
C                   1.1.2  October 31, 2003
C                          - criterion to recognize slowly varying solutions modified;
C                   2.0    april 4, 2005
C                          - order 14 method removed, 
C                          - slightly modification of the I/O interface,
C                          - bug in the order reduction recovery fixed,
C                          - some minor changes in the order variation strategy.
C
C     REFERENCES:
C     -----------
C                 [1] L.BRUGNANO, C.MAGHERINI 
C                     The BiM code for the numerical solution of ODEs
C                     Jour. CAM 164-165 (2004) 145-158.
C
C                 [2] L.BRUGNANO, C.MAGHERINI 
C                     Some Linear Algebra issues concerning the implementation 
C                     of Blended Implicit Methods 
C                     Numer. Linear Alg. Appl. 12 (2005) 305-314.
C
C                 [3] L.BRUGNANO, C.MAGHERINI
C                     Economical Error Estimates for Block Implicit Methods for 
C                     ODEs via Deferred Correction.
C                     Appl. Numer. Math. 56 (2006) 608-617.
C
C                 [4] L.BRUGNANO, C.MAGHERINI 
C                     Blended Implementation of Block Implicit Methods for ODEs 
C                     Appl. Numer. Math. 42 (2002) 29-45.
C
C                 [5] L.BRUGNANO, D.TRIGIANTE 
C                     Block Implicit Methods for ODEs
C                     in "Recent Trends in Numerical Analysis", D.Trigiante Ed.
C                     Nova Science Publ. Inc., New York, 2001, pp. 81-105.
C
C                 [6] L.BRUGNANO 
C                     Blended Block BVMs (B$_3$VMs): a Family of economical 
C                     implicit methods for ODEs
C                     Jour. CAM 116 (2000) 41-62.
C                 
C
C -----------------------------------------------------------------------------------
C -----------------------------------------------------------------------------------
C
C     USAGE:
C     ------
C
C      CALL BIM(M,FCN,T0,TEND,Y0,H,RTOL,ATOL,JAC,IJAC,MLJAC,MUJAC,
C     &               WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IOUT,IDID)
C
C
C     NOTE:   IN ORDER TO GAIN THE BEST PERFORMANCE, THE EXECUTABLE HAS TO 
C     -----   BE CREATED WITH THE OPTION ALLOWING  TO CONTINUE THE EXECUTION 
C             AFTER A FLOATING-POINT EXCEPTION  (E.G., BY USING THE OPTION 
C             -FPE; SEE YOUR FORTRAN COMPILER REFERENCE MANUAL).       
C			THE ISNAN LOGICAL FUNCTION IS REQUIRED, TO RECOGNIZE NANs. IF
C			NOT SUPPORTED BY YOUR COMPILER, A STANDARD ONE IS PROVIDED AT
C			THE TOP OF THE SUBBIM.F FILE.
C
C
C -----------------------------------------------------------------------------------
C           INPUT PARAMETERS
C -----------------------------------------------------------------------------------
C
C M		  SIZE OF THE PROBLEM
C
C FCN		  SUBROUTINE WITH THE FUNCTION F(T,Y) TO BE INTEGRATED. IT IS IN THE FORM:
C
C      subroutine fcn(m,t,y,dy,ierr,rpar,ipar)
C      double precision t,y,dy,rpar(*)
C      integer m,ierr,ipar(*)
C      dimension y(m),dy(m)
CC     m      size of the continuous problem
CC     t,y    is the point where f is evaluated
CC     dy     will contain the value of f(t,y)
CC     ierr   is a return code (0 means OK)
CC     rpar   possible external real parameters
CC     ipar   possible external integer parameters
C      ................
C      return
C      end
C
C
C T0-TEND   INTEGRATION INTERVAL
C
C Y0        INITIAL CONDITION
C
C H         INITIAL STEPSIZE
C
C RTOL-ATOL RELATIVE AND ABSOLUTE TOLERANCES
C
C JAC       SUBROUTINE EVALUATING THE JACOBIAN OF F (DUMMY, IF IJAC=0). 
C           IF IJAC .NE. 0, IT IS IN THE FORM:
C
C      subroutine jac(m,t,y,jac,ldjac,ierr,rpar,ipar)
C      double precision t,y,jac,rpar(*)
C      integer neqn,ldjac,ierr,ipar(*)
C      dimension y(m),jac(ldjac,m)
CC     m      size of the continuous problem
CC     t,y	is the point where the Jacobian is evaluated
CC     jac	will contain the value of the Jacobian at (t,y)
CC     ldjac  leading dimension of the array ldjac
CC     ierr	is a return code (0 means OK)
CC     rpar   possible external real parameters
CC     ipar	possible external integer parameters
C      ............
C      return
C      end
C
C
C IJAC      FLAG: 0=NUMERICAL JACOBIAN, ANALYTICAL OTHERWISE
C
C MLJAC-MUJAC	LOWER-UPPER BANDWIDTH OF THE JACOBIAN (MLJAC=M IF FULL JACOBIAN) 				
C
C LWORK     LENGTH OF WORK   ( LWORK >= 14 +KMAX +8*M +4*KMAX*M +M*(LDJAC+LDLU),
C
C            WHERE:
C
C                            LDJAC = LDLU = M,            IN CASE OF A FULL JACOBIAN, 
C            LDJAC = MLJAC+MUJAC+1,  LDLU = LDJAC+MLJAC,  IN CASE OF A SPARSE JACOBIAN;   
C
C                             KMAX = ORDMAX-2,            IF ORDMAX>4,                      
C                                    3,                   IF ORDMAX=4. )                     
C
C WORK(1)   UROUND. MACHINE PRECISION. (DEFAULT = 1.D-16)
C
C WORK(2)   HMAX. MAXIMUM INTEGRATION STEP. (DEFAULT = (TEND-T0)/8)
C
C WORK(3)   FACTOL - SAFETY FACTOR FOR THE STOPPING CRITERION OF THE BLENDLED ITERATION.
C           METHOD OF ORDER 4. (DEFAULT = 1.D-1)
C
C WORK(4)   FACTOL - SAFETY FACTOR FOR THE STOPPING CRITERION OF THE BLENDLED ITERATION.
C           METHOD OF ORDER 6. (DEFAULT = 1.D-1)
C
C WORK(5)   FACTOL - SAFETY FACTOR FOR THE STOPPING CRITERION OF THE BLENDLED ITERATION.
C           METHOD OF ORDER 8. (DEFAULT = 1.D-1)
C
C WORK(6)   FACTOL - SAFETY FACTOR FOR THE STOPPING CRITERION OF THE BLENDLED ITERATION.
C           METHOD OF ORDER 10. (DEFAULT = 1.D-1)
C
C WORK(7)   FACTOL - SAFETY FACTOR FOR THE STOPPING CRITERION OF THE BLENDLED ITERATION.
C           METHOD OF ORDER 12. (DEFAULT = 1.D-1)
C
C WORK(8)   FACTOL - SAFETY FACTOR FOR THE STOPPING CRITERION OF THE BLENDLED ITERATION
C           IN CASE OF SMALL VALUES OF min(abs(y_0)), min(abs(f_0)) AND OF max(abs(f_0)).
C           (DEFAULT = 5D-3)  
C
C WORK(9)   FACTOL - SAFETY FACTOR FOR THE STOPPING CRITERION OF THE BLENDLED ITERATION
C           IN CASE OF SLOWLY VARYING SOLUTIONS (DEFAULT = 5D-2)  
C 
C WORK(10)-WORK(11)  FACL-FACR. THE NEW STEPSIZE MUST SATISFY
C                            FACL<=HNEW/HOLD<= FACR.
C           (DEFAULT: WORK(10)=1.2D-1, WORK(11)=10D0)
C
C WORK(12)  SFTY - SAFETY FACTOR FOR PREDICTING THE NEW STEPSIZE FOR THE CURRENT ORDER 
C           METHOD. (DEFAULT = 1D0/2D1)
C
C WORK(13)  SFTYUP - SAFETY FACTOR FOR PREDICTING THE NEW STEPSIZE FOR THE HIGHER ORDER 
C           METHOD. (DEFAULT = SFTY/2D0)
C
C WORK(14)  SFTYDN - SAFETY FACTOR FOR PREDICTING THE NEW STEPSIZE FOR THE LOWER ORDER
C           METHOD. (DEFAULT = SFTY) 
C
C LIWORK    LENGTH OF IWORK  (LIWORK >= M+37)     
C
C IWORK(1)  MAX NUMBER OF INTEGRATION STEPS (DEFAULT = 100000).
C
C IWORK(2)  ORDMIN, 4<=ORDMIN<=12. (DEFAULT = 4).
C
C IWORK(3)  ORDMAX, ORDMIN<=ORDMAX<=12. (DEFAULT = 12).
C
C IWORK(4)  MAX NUMBER OF BLENDED ITERATIONS PER INTEGRATION STEP, METHOD OF ORDER 4
C           (DEFAULT = 10).
C
C IWORK(5)  MAX NUMBER OF BLENDED ITERATIONS PER INTEGRATION STEP, METHOD OF ORDER 6
C           (DEFAULT = 12).
C
C IWORK(6)  MAX NUMBER OF BLENDED ITERATIONS PER INTEGRATION STEP, METHOD OF ORDER 8
C           (DEFAULT = 14).
C
C IWORK(7)  MAX NUMBER OF BLENDED ITERATIONS PER INTEGRATION STEP, METHOD OF ORDER 10
C           (DEFAULT = 16).
C
C IWORK(8)  MAX NUMBER OF BLENDED ITERATIONS PER INTEGRATION STEP, METHOD OF ORDER 12
C           (DEFAULT = 18).
C
C RPAR,IPAR   REAL AND INTEGER PARAMETERS (OR PARAMETER ARRAYS) WHICH  
C             CAN BE USED FOR COMMUNICATION BETWEEN YOUR CALLING
C             PROGRAM AND THE FCN, JAC SUBROUTINES.    
C
C IOUT      OUTPUT FLAG, SET TO 1, IN CASE OF OUTPUT AT EACH SUCCESFULL STEP.
C           IF NOT SET TO 1, THE SUBROUTINE SOLOUT MAY BE DUMMY. OTHERWISE, 
C           IT MUST BE PROVIDED IN THE FOLLOWING FORM: 
C
C      subroutine solout(m,t,y,f,k,ord,irtrn)   
CC     m			is the size of the problem
CC     (t,y)		is the current point
CC     f			is the derivative of y
CC     k			is the block-size of the method
CC     ord		is the order of the method
CC     irtrn		is a return code (0 means that everything is OK)
C      ................
C      return
C      end
C
C
C -----------------------------------------------------------------------------------
C           OUTPUT PARAMETERS
C -----------------------------------------------------------------------------------
C
C
C T0        VALUE OF T UP TO WHERE THE SOLUTION HAS BEEN COMPUTED
C           (IF THE INTEGRATION HAS BEEN SUCCESFULL,THEN T0=TEND) 
C
C Y0        NUMERICAL SOLUTION AT T0    
C
C IDID      RETURN CODE:
C              0  SUCCESFULL RUN
C             -1  WRONG INPUT PARAMETERS
C             -2  A LARGER NMAX IS NEEDED
C             -3  STEPSIZE TOO SMALL
C             -4  REPEATEDLY SINGULAR MATRIX
C             -5  TOO MANY CONSECUTIVE NEWTON FAILURES
C             -6  ERROR CODE RETURNED BY THE JAC SUBROUTINE 
C                 OR BY THE FCN SUBROUTINE AT THE STARTING
C                 POINT
C
C IWORK( 9) NUMBER OF FUNCTION EVALUATIONS
C
C IWORK(10) NUMBER OF JACOBIAN EVALUATIONS
C
C IWORK(11) NUMBER OF LU DECOMPOSITION
C
C IWORK(12) NUMBER OF LINEAR SYSTEMS SOLVED
C
C IWORK(13)-IWORK(17) NUMBER OF BLENDED ITERATIONS PER METHOD
C
C IWORK(18)-IWORK(22) NUMBER OF STEPS PER METHOD
C
C IWORK(23)-IWORK(27) NUMBER OF ACCEPTED STEPS PER METHOD
C
C IWORK(28)-IWORK(32) NUMBER OF REFUSED STEPS PER METHOD (ERROR TEST)
C
C IWORK(33)-IWORK(37) NUMBER OF REFUSED STEPS PER METHOD (NEWTON'S CONVERGENCE)
C
C -----------------------------------------------------------------------------------
C -----------------------------------------------------------------------------------
C
Cc-----------------------------------------------------------------------
Cc     Sample driver for the code BIM
Cc-----------------------------------------------------------------------
C      program BIMDO
C      implicit none
C      integer MMAX,lwork,liwork
C      parameter(MMAX=3,lwork=24+MMAX*(48+2*MMAX),liwork=MMAX+37)
C      double precision y(MMAX),work(lwork)
C      integer iwork(lwork),ijac,iout
C      external feval, jeval,  solout
C      character problm*8
C      double precision t0,tf,h0,h,rtol,atol,
C     &                 y0(MMAX), rpar(1)
C      integer neqn,mljac,mujac,itol,i,ierr, ipar(1)
C      integer NSTEPS,NACCEPT,NFAILERR,NFAILNEWT,NITER,idid
Cc-----------------------------------------------------------------------
Cc     get the problem dependent parameters
Cc-----------------------------------------------------------------------
C      neqn = 3
C      t0   = 0d0
C      tf   = 1d15
C      ijac = 1
C      mljac = 3
C      mujac = 3
C      y(1) = 1d0
C      y(2) = 0d0
C      y(3) = 0d0
C
Cc-----------------------------------------------------------------------
Cc     read the tolerances and initial stepsize
Cc-----------------------------------------------------------------------
C      write(6,*) 'give the absolute tolerance'
C      read(5,*)   atol
C      write(6,*) 'give the relative tolerance'
C      read(5,*)   rtol
C      write(6,*) 'give the initial  stepsize '
C      read(5,*)   h0
C
C      h = h0
C      do i=1,8
C         iwork(i) = 0
C      end do
C      do i=1,14
C         work(i) = 0d0
C      end do
C      iout = 0
C      idid = 0
C
Cc-----------------------------------------------------------------------
Cc     call of the subroutine BIM
Cc-----------------------------------------------------------------------
C      call BIM(neqn,feval,t0,tf,y,h,rtol,atol,
C     &         jeval,ijac,mljac,mujac,
C     &         work,lwork,iwork,liwork,
C     &         rpar,ipar,iout,idid)
C      if (idid.ne.0) then
C         write(6,*) 'ERROR: returned idid =', idid
C         goto 20
C      endif
Cc-----------------------------------------------------------------------
Cc     print final solution
Cc-----------------------------------------------------------------------
C      write(6,10)
C   10 format(//)
C
C      write(6,11) atol,rtol,h0
C   11 format(/,' we solved the problem with',//,
C     +       '       absolute tolerance = ',d10.4,',',/,
C     +       '       relative tolerance = ',d10.4,',',/,
C     +       '     and initial stepsize = ',d10.4,//)
Cc-----------------------------------------------------------------------
Cc     print error with respect to reference solution
Cc-----------------------------------------------------------------------
C      NSTEPS    = 0
C      NACCEPT   = 0
C      NFAILERR  = 0
C      NFAILNEWT = 0
C      NITER     = 0
C      DO I=1,5
C       NITER     = NITER     + iwork(I+12)
C       NSTEPS    = NSTEPS    + iwork(I+17)
C       NACCEPT   = NACCEPT   + iwork(I+22)
C       NFAILERR  = NFAILERR  + iwork(I+27)
C       NFAILNEWT = NFAILNEWT + iwork(I+32)
C      END DO
C
C      write(6,41) NSTEPS,IWORK(18),IWORK(19),IWORK(20),IWORK(21),
C     &            IWORK(22),
C     &            NACCEPT,NFAILNEWT,NFAILERR,
C     &            IWORK(9),IWORK(10),IWORK(11),IWORK(12),
C     &            NITER
C   41 format(  ///,
C     +         ' # Steps              = ',i8,/
C     +         ' # Composition        = ',i8,i8,i8,i8,i8,/
C     +         ' # Accept             = ',i8,/,
C     +         ' # Failnwt            = ',i8,/,
C     +         ' # Failerr            = ',i8,/,
C     +         ' # F-eval             = ',i8,/,
C     +         ' # Jac-eval           = ',i8,/,
C     +         ' # LU-decomp          = ',i8,/,
C     +         ' # Linear systems     = ',i8,/,
C     +         ' # Newt. iterat.      = ',i8)
C      write(6,*) 'Yf = '
C      do i = 1, neqn
C		write(6,*) y(i)
C	 end do
C20    continue
C      end
C
Cccccccccccccccccccccccccccccccc AUXILIARY SUBROUTINES 
C
C      subroutine feval(neqn,t,y,dy,ierr,rpar,ipar)
C      integer neqn,ierr,ipar
C      double precision t,y(neqn),dy(neqn),rpar
C      dy(1) = -.04d0*y(1) + 1.d4*y(2)*y(3)
C      dy(3) = 3.d7*y(2)*y(2)
C      dy(2) = -dy(1) - dy(3)
C      return
C      end
C
C      subroutine jeval(neqn,t,y,jac,ldim,ierr,rpar,ipar)
C      integer neqn,ldim,ierr,ipar
C      double precision t,y(neqn),jac(ldim,neqn),rpar
C
C      integer i,j
C
C      do 20 j=1,neqn
C         do 10 i=1,neqn
C            jac(i,j)=0d0
C   10    continue
C   20 continue
C
C      jac(1,1) = -.04d0
C      jac(1,2) = 1.d4*y(3)
C      jac(1,3) = 1.d4*y(2)
C      jac(2,1) = .04d0
C      jac(2,3) = -jac(1,3)
C      jac(3,1) = 0.0
C      jac(3,2) = 6.d7*y(2)
C      jac(3,3) = 0.0
C      jac(2,2) = -jac(1,2) - jac(3,2)
C      return
C      end
C
C      subroutine solout(m,t,y,f,k,ord,irtrn)
C      implicit none
C      integer m,k,ord,irtrn
C      double precision t(k),y(m,k),f(m,k)
C      integer i,j
C 20   format (e22.16,a1)
C      write(20,20) t(k)
C      do i=1,m
C         write(20,20) y(i,k),' '
C      end do
C      write(20,30) ' '
C 30   format(A1,/)
C      irtrn=0
C      return
C      end
C
C -----------------------------------------------------------------------------------
C -----------------------------------------------------------------------------------
      IMPLICIT NONE

      EXTERNAL FCN,JAC
      INTEGER M,LWORK,LIWORK,IWORK(LIWORK),
     &        IJAC,MLJAC,MUJAC,LDJAC,
     &        LDLU,IJOB,IPAR(*),IOUT,IDID
      LOGICAL JBAND
      DOUBLE PRECISION T0,TEND,Y0(M),H,RTOL,ATOL,WORK(LWORK),RPAR(*)

      INTEGER NMETH,KMAX
      PARAMETER (NMETH=5,KMAX=10)
      
      INTEGER MAXSTEP,ORDMIN,ORDMAX,ITMAX(NMETH),STEP_ORD(NMETH),
     &        NFEVAL,NJEVAL,NLU,NLINSYS,NITER(NMETH),NSTEP(NMETH),
     &        NACCEPT(NMETH),NFAILERR(NMETH),NFAILNEWT(NMETH)

      DOUBLE PRECISION UROUND, FACNEWTV(NMETH),FACNSMALL,
     &                 FACNRESTR,FACL,FACR,
     &                 SFTY, SFTYUP, SFTYDN, HMAX,
     &                 RHOMUV(NMETH),RHOMLV(NMETH),
     &                 TOLESTRAP(2)

      INTEGER I,INDF0,INDT,INDIPVT,INDEJ0,INDY,INDF,INDTHETA,INDJ0,
     &        INDERR,INDSCAL,INDDELJ0,INDDELJ0OLD,
     &        INDFJ0,indord,INDSCALEXT,IND_DD

      LOGICAL STOPINT


      STEP_ORD(1)= 3
      STEP_ORD(2)= 4
      STEP_ORD(3)= 6
      STEP_ORD(4)= 8
      STEP_ORD(5)= 10


      STOPINT = .FALSE.

C     INITIAL STEP-SIZE AND TOLERANCES
      IF (H.EQ.0D0) THEN 
          H=1.D-6
      ELSEIF(H.LT.0D0) THEN
          WRITE(6,*) 'WRONG INPUT H=',H
          STOPINT=.TRUE.
      END IF

C--------------------------------------------------
C     PARAMETERS INITIALIZATION
C--------------------------------------------------

      IF (IWORK(1).EQ.0) THEN
         MAXSTEP=100000
      ELSE
         MAXSTEP=IWORK(1)
         IF (MAXSTEP.LE.0) THEN
             WRITE(6,*) 'WRONG INPUT IWORK(1)=',IWORK(1)
             STOPINT=.TRUE.
         END IF
      ENDIF

      IF (IWORK(2).EQ.0) THEN
         ORDMIN = 4
      ELSE
         ORDMIN = IWORK(2)
         indord = ORDMIN/2-1
         IF ((indord.LE.0).OR.(indord.GT.NMETH)) THEN
             WRITE(6,*) 'WRONG INPUT IWORK(2)=',IWORK(2)
             STOPINT=.TRUE.
         END IF
         ORDMIN = 2*(indord+1)
      ENDIF

      IF (IWORK(3).EQ.0) THEN
         ORDMAX = 12
         indord = NMETH
      ELSE
         ORDMAX = IWORK(3)
         indord = ORDMAX/2-1
         IF ((indord.LE.0).OR.(indord.GT.NMETH)) THEN
             WRITE(6,*) 'WRONG INPUT IWORK(3)=',IWORK(3)
             STOPINT=.TRUE.
         END IF
         ORDMAX = 2*(indord+1)
      ENDIF

      IF (ORDMIN.GT.ORDMAX) THEN
        WRITE(6,1000) IWORK(2),IWORK(3)
 1000   FORMAT(/,/,'INVALID VALUES FOR IWORK(2)=',i3,' (ORDMIN)',/,
     &             '               AND IWORK(3)=',i3,' (ORDMAX)',/,/)
        STOPINT=.TRUE.
      END IF

      IF (IWORK(4) .EQ. 0) THEN
         ITMAX(1) = 10
      ELSE
         ITMAX(1) = IWORK(4)
         IF (ITMAX(1).LE.0) THEN
             WRITE(6,*) 'WRONG INPUT IWORK(4)=',IWORK(4)
             STOPINT=.TRUE.
         END IF
      END IF         

      IF (IWORK(5) .EQ. 0) THEN
         ITMAX(2) = 12
      ELSE
         ITMAX(2) = IWORK(5)
         IF (ITMAX(2).LE.0) THEN
             WRITE(6,*) 'WRONG INPUT IWORK(5)=',IWORK(5)
             STOPINT=.TRUE.
         END IF
      END IF         

      IF (IWORK(6) .EQ. 0) THEN
         ITMAX(3) = 14
      ELSE
         ITMAX(3) = IWORK(6)
         IF (ITMAX(3).LE.0) THEN
             WRITE(6,*) 'WRONG INPUT IWORK(6)=',IWORK(6)
             STOPINT=.TRUE.
         END IF
      END IF         

      IF (IWORK(7) .EQ. 0) THEN
         ITMAX(4) = 16
      ELSE
         ITMAX(4) = IWORK(7)
         IF (ITMAX(4).LE.0) THEN
             WRITE(6,*) 'WRONG INPUT IWORK(7)=',IWORK(7)
             STOPINT=.TRUE.
         END IF
      END IF         

      IF (IWORK(8) .EQ. 0) THEN
         ITMAX(5) = 18
      ELSE
         ITMAX(5) = IWORK(8)
         IF (ITMAX(5).LE.0) THEN
             WRITE(6,*) 'WRONG INPUT IWORK(8)=',IWORK(8)
             STOPINT=.TRUE.
         END IF
      END IF         


      IF (WORK(1) .EQ. 0D0) THEN
         UROUND = 1.0D-16
      ELSE
         UROUND = WORK(1)
         IF ((UROUND.LE.0D0).OR.(UROUND.GE.1D0)) THEN
             WRITE(6,*) 'WRONG INPUT WORK(1)=',WORK(1)
             STOPINT=.TRUE.
         END IF
      END IF
      
      IF (ATOL.LE.0D0 .OR. RTOL.LE.UROUND) THEN
          WRITE(6,*) 'TOLERANCES ARE TOO SMALL'
          STOPINT = .TRUE.
      END IF

      IF (WORK(2) .LE. 0D0) THEN
         HMAX = (TEND-T0)/8d0
      ELSE
         HMAX = WORK(2)
         IF (HMAX.GT.(TEND-T0)) HMAX = TEND-T0
      END IF

      IF (WORK(3) .EQ. 0D0) THEN
         FACNEWTV(1) = 1D-1
      ELSE
         FACNEWTV(1) = WORK(3)
         IF ((FACNEWTV(1).LE.0D0).OR.(FACNEWTV(1).GE.1D0)) THEN
             WRITE(6,*) 'WRONG INPUT WORK(3)=',WORK(3)
             STOPINT=.TRUE.
         END IF
      END IF

      IF (WORK(4) .EQ. 0D0) THEN
         FACNEWTV(2) = 1D-1
      ELSE
         FACNEWTV(2) = WORK(4)
         IF ((FACNEWTV(2).LE.0D0).OR.(FACNEWTV(2).GE.1D0)) THEN
             WRITE(6,*) 'WRONG INPUT WORK(4)=',WORK(4)
             STOPINT=.TRUE.
         END IF
      END IF

      IF (WORK(5) .EQ. 0D0) THEN
         FACNEWTV(3) = 1D-1
      ELSE
         FACNEWTV(3) = WORK(5)
         IF ((FACNEWTV(3).LE.0D0).OR.(FACNEWTV(3).GE.1D0)) THEN
             WRITE(6,*) 'WRONG INPUT WORK(5)=',WORK(5)
             STOPINT=.TRUE.
         END IF
      END IF

      IF (WORK(6) .EQ. 0D0) THEN
         FACNEWTV(4) = 1D-1
      ELSE
         FACNEWTV(4) = WORK(6)
         IF ((FACNEWTV(4).LE.0D0).OR.(FACNEWTV(4).GE.1D0)) THEN
             WRITE(6,*) 'WRONG INPUT WORK(6)=',WORK(6)
             STOPINT=.TRUE.
         END IF
      END IF

      IF (WORK(7) .EQ. 0D0) THEN
         FACNEWTV(5) = 1D-1
      ELSE
         FACNEWTV(5) = WORK(7)
         IF ((FACNEWTV(5).LE.0D0).OR.(FACNEWTV(5).GE.1D0)) THEN
             WRITE(6,*) 'WRONG INPUT WORK(7)=',WORK(7)
             STOPINT=.TRUE.
         END IF
      END IF

      IF (WORK(8) .EQ. 0D0) THEN
         FACNSMALL = 5d-3
      ELSE
         FACNSMALL = WORK(8)
         IF ((FACNSMALL.LE.0D0).OR.(FACNSMALL.GE.1D0)) THEN
             WRITE(6,*) 'WRONG INPUT WORK(8)=',WORK(8)
             STOPINT=.TRUE.
         END IF
      END IF

      IF (WORK(9) .EQ. 0D0) THEN
         FACNRESTR = 5d-2
      ELSE
         FACNRESTR = WORK(9)
         IF ((FACNRESTR.LE.0D0).OR.(FACNRESTR.GE.1D0)) THEN
             WRITE(6,*) 'WRONG INPUT WORK(9)=',WORK(9)
             STOPINT=.TRUE.
         END IF
      END IF

      IF (WORK(10) .EQ. 0D0) THEN
         FACL = 1.2D-1
      ELSE
         FACL = WORK(10)
         IF (FACL.LT.0D0) THEN
             WRITE(6,*) 'WRONG INPUT WORK(10)=',WORK(10)
             STOPINT=.TRUE.
         END IF
      END IF

      IF (WORK(11) .EQ. 0D0) THEN
         FACR = 1D1
      ELSE
         FACR = WORK(11)
         IF(FACR.LE.0D0) THEN
             WRITE(6,*) 'WRONG INPUT WORK(11)=',WORK(11)
             STOPINT=.TRUE.
         END IF
         IF(FACL.GE.FACR) THEN
            WRITE(6,1010) WORK(10),WORK(11)
 1010       FORMAT(/,/,'INVALID VALUES FOR WORK(10)=',e10.2,' (FACL)',/,
     &             '               AND WORK(11)=',e10.2,' (FACR)',/,/)
            STOPINT=.TRUE.
         END IF
      END IF

      IF (WORK(12) .EQ. 0D0) THEN
         SFTY = 1D0/20D0
      ELSE
         SFTY = WORK(12)
         IF(SFTY.LE.0D0) THEN
             WRITE(6,*) 'WRONG INPUT WORK(12)=',WORK(12)
             STOPINT=.TRUE.
         END IF
      END IF

      IF (WORK(13) .EQ. 0D0) THEN
         SFTYUP = .5d0*SFTY
      ELSE
         SFTYUP = WORK(13)
         IF(SFTYUP.LE.0D0) THEN
             WRITE(6,*) 'WRONG INPUT WORK(13)=',WORK(13)
             STOPINT=.TRUE.
         END IF
      END IF
      
      IF (WORK(14) .EQ. 0D0) THEN
         SFTYDN = SFTY
      ELSE
         SFTYDN = WORK(14)
         IF(SFTYDN.LE.0D0) THEN
             WRITE(6,*) 'WRONG INPUT WORK(14)=',WORK(14)
             STOPINT=.TRUE.
         END IF
      END IF

      IF (STOPINT) THEN
C     INVALID INPUT PARAMETERS
         IDID = -1
         RETURN
      END IF

C---------------------------------------------------------
C     FIXED PARAMETERS
C---------------------------------------------------------

      RHOMUV(1) = 1d-2*DABS(DLOG10(DMIN1(RTOL,ATOL,1D-1)))
      RHOMLV(1) = 5d-1
      DO I=2,NMETH
         RHOMUV(I) = RHOMUV(I-1)**(DBLE(STEP_ORD(I))
     &                             /DBLE(STEP_ORD(I-1)))
         RHOMLV(I) = RHOMLV(I-1)**(DBLE(STEP_ORD(I))
     &                             /DBLE(STEP_ORD(I-1)))
      end do
      TOLESTRAP(1) = DMIN1(1D-2,1D2*RTOL)
      TOLESTRAP(2) = DMIN1(1D-2,1D2*ATOL)
      
C---------------------------------------------------------
C     BANDED JACOBIAN
C---------------------------------------------------------

      JBAND = (MLJAC .LT. M)
      IF (JBAND) THEN
         LDJAC = MLJAC+MUJAC+1
         LDLU  = MLJAC+LDJAC
         IJOB=2
      ELSE
         LDJAC = M
         LDLU  = M
         IJOB  = 1
      END IF

C---------------------------------------------------------
C     COMPUTE THE VECTORS ENTRY-POINT IN IWORK AND WORK
C---------------------------------------------------------

      INDIPVT = 38
      IF ((INDIPVT + M-1) .GT. LIWORK) THEN
        WRITE(6,*) 'INSUFF. STORAGE FOR IWORK, MIN. LIWORK=',INDIPVT+M-1
        IDID = -1
        RETURN
      END IF
      
      INDF0       = 15
      INDT        = INDF0       + M
      INDY        = INDT        + STEP_ORD(indord)
      INDF        = INDY        + M*STEP_ORD(indord)
      INDTHETA    = INDF        + M*STEP_ORD(indord)
      INDERR      = INDTHETA    + M*LDLU
      INDSCAL     = INDERR      + M*STEP_ORD(indord)
      INDSCALEXT  = INDSCAL     + M
      INDJ0       = INDSCALEXT  + M
      INDDELJ0    = INDJ0       + M*LDJAC
      INDDELJ0OLD = INDDELJ0    + M
      IND_DD      = INDDELJ0OLD + M
      INDFJ0      = IND_DD      + M*(STEP_ORD(indord)+1)
      INDEJ0      = INDFJ0      + M
      IF ((INDEJ0 + M-1) .GT. LWORK) THEN
        WRITE(6,*) 'INSUFF. STORAGE FOR WORK, MIN. LWORK=',INDEJ0+M-1
        IDID = -1
        RETURN
      END IF

      CALL  BIM0(M,FCN,JAC,NMETH,STEP_ORD(indord),Y0,WORK(INDF0),
     &           T0,TEND,H,RTOL,ATOL,
     &           MAXSTEP,ORDMIN,ORDMAX,ITMAX,UROUND,HMAX,FACNEWTV,
     &           FACNSMALL,FACNRESTR,FACL,FACR,SFTY,SFTYUP,SFTYDN,
     &           RHOMUV,RHOMLV,
     &           IWORK(9),IWORK(10),IWORK(11),IWORK(12),IWORK(13),
     &           IWORK(18),IWORK(23),IWORK(28),IWORK(33),
     &           IWORK(INDIPVT),STEP_ORD,          
     &           WORK(INDT),WORK(INDY),WORK(INDF),WORK(INDTHETA),
     &           WORK(INDERR),WORK(INDJ0),WORK(INDDELJ0),
     &           WORK(INDDELJ0OLD),WORK(INDFJ0),WORK(INDEJ0),
     &           WORK(INDSCAL),
     &           WORK(IND_DD),TOLESTRAP,WORK(INDSCALEXT),
     &           IJAC,MLJAC,MUJAC,LDJAC,LDLU,JBAND,IJOB,
     &           RPAR,IPAR,IOUT,IDID)

      RETURN
      END

c------------------------------------------------------------------------------------------------------
C    CORE INTEGRATOR
C------------------------------------------------------------------------------------------------------

      SUBROUTINE  BIM0(M,FCN,JAC,NMETH,KMAX,Y0,F0,T0,TEND,H,RTOL,ATOL,
     &                 MAXSTEP,ORDMIN,ORDMAX,ITMAX,UROUND,HMAX,FACNEWTV,
     &                 FACNSMALL,FACNRESTR,FACL,FACR,SFTY,SFTYUP,SFTYDN,
     &                 RHOMUV,RHOMLV,
     &                 NFEVAL,NJEVAL,NLU,NLINSYS,NITER,NSTEP,NACCEPT,
     &                 NFAILERR,NFAILNEWT,
     &                 IPVT,STEP_ORD,T,Y,F,THETA,ERR,J0,
     &                 DELJ0,DELJ00,FJ0,EJ0,SCAL,
     &                 DD,TOLESTRAP,SCALEXTRAP,
     &                 IJAC,MLJAC,MUJAC,LDJAC,LDLU,JBAND,IJOB,
     &                 RPAR,IPAR,IOUT,IDID)
 
      IMPLICIT NONE
C
C INPUT PARAMETERS
C

      EXTERNAL FCN,JAC
      INTEGER M,NMETH,KMAX,MAXSTEP,ORDMIN,ORDMAX,ITMAX(NMETH),
     &        IPVT(M),STEP_ORD(NMETH),
     &        MLJAC,MUJAC,LDJAC,LDLU,IJOB,IJAC,IPAR(*)
      LOGICAL JBAND

      DOUBLE PRECISION TEND,H,RTOL,ATOL,UROUND,HMAX,
     &                 FACNEWTV(NMETH),FACNSMALL,FACNRESTR,FACL,FACR,
     &                 SFTY,SFTYUP,SFTYDN,
     &                 T(KMAX),F0(M),THETA(LDLU,M),J0(LDJAC,M),DELJ0(M),
     &                 DELJ00(M),Y(M,KMAX),F(M,KMAX),ERR(M,KMAX),
     &                 RHOMUV(NMETH),RHOMLV(NMETH),SCAL(M),DD(KMAX+1,M),
     &                 TOLESTRAP(2),FJ0(M),SCALEXTRAP(M),EJ0(M),
     &                 RPAR(*)

C
C I/O PARAMETERS
C

      DOUBLE PRECISION T0,Y0(M)

C
C OUTPUT PARAMETERS
C

      INTEGER NFEVAL,NJEVAL,NLU,NLINSYS,NITER(NMETH),NSTEP(NMETH),
     &        NACCEPT(NMETH),NFAILERR(NMETH),NFAILNEWT(NMETH),
     &        IOUT,IDID

C
C LOCAL VARIABLES
C

      INCLUDE 'param.dat'
 
      INTEGER I,J,IT,MAXIT,ORD,ORD_IND,ORDNEW,K,KOLD,KNEW,KUP,K0,
     &        INDMIN,NFAILCONV,NORD,INFO,IRTRN,IERR,IERR0,
     &        NFAILCONS,NSTEPS,NSING,NERROR
      DOUBLE PRECISION NERR,NERRUP,NERRUP1,NERRDOWN,
     &                 NERR0,NERROLD,NERRSTOP,
     &                 RHO,RHO0,RHOBAD,RHOMU,RHOML,RHOT,RHOTUP,
     &                 RHOOLD,RHOI,RHOIUP,
     &                 RHOMAX,HNEW,HNUP,HNDN,H0,H00,RATH,RATRHO,
     &                 FI,ESP,ESPUP,ESPDN,
     &                 FACNEWT,GAMMA,HGAMMA,
     &                 MINY0,MAXF0,MAXDELTA,
     &                 NU, NUUP,NUDN, NU1, NUUP1, FATERR,HNUP1,FI1,
     &                 SCALJ0_1,NJ0,NJ00,HJ0,RHOJ0,
     &                 DJ0,FATDJ0,FATDJ0I,
     &                 ABSY0,
     &                 DELTAH2,DELTAH1SF,DELTAH1,CFAT1,CFAT2,CFAT3,
     &                 RATHH,HFATT,ALFAFATT,RHOFATT,DISCR,
     &                 tolrhoJ0,RTOLATOL,
     &                 DELT,YSAFE,ITNEW,
     &                 RHOUP,RHONEW,RHODN,
     &                 VUP,SISERR,SISERRDN,SISERRUP,CSIS,CFACT,
     &                 NF0,NF,RTOL0

      LOGICAL JVAI,LAST,EXTRAP,EXTRAP0,CALJAC,CALFACT,
     &        SUCCESS,RESTRICT,
     &        TRUEJAC,STAGNA,
     &        LINEAR,CALJAC0,CALFACT0,
     &        QINF,QINFJ,QINFF,NQINF,SMALLM,
     &        ESTIM,ESTIM1,NODJ0,NODJ00,ISNAN,
     &        ERROR

C-------------------------------------------------------------------------------------------------------------------
C STATISTICS
C
      NFEVAL = 0
      NJEVAL = 0
      NLU    = 0
      NLINSYS= 0 
      NSTEPS = 0
      DO I=1,NMETH
         NITER(I)     = 0
         NSTEP(I)     = 0
         NACCEPT(I)   = 0
         NFAILERR(I)  = 0
         NFAILNEWT(I) = 0
      END DO

C-------------------------------------------------------------------------------------------------------------------
C     OTHER INITIALIZATIONS

      IF (JBAND) THEN
         CSIS  = DBLE(2*M*(MLJAC+MUJAC))
         CFACT = DBLE(2*M*MLJAC*MUJAC)
      ELSE
         CSIS  = DBLE(2*M**2)
         CFACT = DBLE(2*M**3)/3D0
      END IF
      SMALLM   = M.LE.5

C     VECTOR TO BE USED FOR THE ESTIMATE OF JACOBIAN VARIATION
      NERR     = 0.3141592654D0
      NERRUP   = 0D0
      DO I = 1,M
         NERR   = 4D0*NERR*(1D0-NERR)
         EJ0(I) = NERR
         IF (NERR.LT.5D-1) EJ0(I)=EJ0(I)+5D-1
         NERRUP = DMAX1(NERRUP,EJ0(I))
      END DO
      DO I = 1,M
         EJ0(I)=EJ0(I)/NERRUP
      END DO
    
      NFAILCONV  = 0
      NFAILCONS  = 0

      ORD     = ORDMIN
      ORD_IND = INT(ORD/2) - 1
      K       = STEP_ORD(ORD_IND)
      KOLD    = 0

      H       = DMIN1(H,HMAX)
C -------------------------------------------------------
C     INITIAL STEPSIZE TOO SMALL !!!!!
      IF (H .LT. 1D1*UROUND*T0) H = 2D1*UROUND*T0
C -------------------------------------------------------

      RTOL0 = RTOL
      RTOL  = DMAX1(RTOL,5D1*UROUND/SFTY)
      RTOLATOL =  RTOL/ATOL

      RHO      =  0D0
      MAXDELTA =  0D0
      
      LAST     = .FALSE.
      EXTRAP   = .FALSE.
      EXTRAP0  = .FALSE.
      RESTRICT = .FALSE.
      SUCCESS  = .FALSE.
      CALJAC   = .TRUE.      
      CALJAC0  = .TRUE.
      CALFACT0 = .TRUE.
      LINEAR   = .FALSE.
      
      IERR = 0
      CALL FCN(M,T0,Y0,F0,IERR,RPAR,IPAR)
      IF (IERR.NE.0) THEN 
          WRITE(6,1020) IERR
 1020     FORMAT(/,/,'ERROR CODE IERR = ', I3, /
     &               'RETURNED BY THE SUBROUTINE FCN',/,/)
          WRITE(6,900) T0
          IDID = -6
          GOTO 800
      END IF
      NFEVAL = NFEVAL + 1
      SCALJ0_1 = 0D0
      DO I=1,M
         ABSY0 = DABS(Y0(I))
         SCAL(I)=1D0/(1d0+RTOLATOL*ABSY0)
         SCALEXTRAP(I)=1d0/(1d0+ABSY0)
      END DO

      HFATT=1D0

      IF(IOUT.EQ.1) 
     &  CALL SOLOUT(M,T0,Y0,F0,1,ORD,IRTRN)
      
c     MAIN LOOP
100   CONTINUE

      IF (K .EQ. KOLD) GOTO 140
C     THE METHOD'S ORDER HAS BEEN CHANGED
      ESP = 1D0/DBLE(K+1)
      IF (ORD .LT. ORDMAX) 
     &    ESPUP=1D0/DBLE(STEP_ORD(ORD_IND+1)+1D0)
      IF (ORD .GT. ORDMIN) 
     &    ESPDN=1D0/DBLE(STEP_ORD(ORD_IND-1)+1D0)
      MAXIT    = ITMAX(ORD_IND)
      RHOML    = RHOMLV(ORD_IND)
      RHOMU    = RHOMUV(ORD_IND)      
      RHOOLD   = 0D0
      NORD     = 0
      NERROR   = 0
      ERROR    = .FALSE.
C     DEFINE THE PARAMETERS DEPENDING ON THE METHOD'S ORDER
      GOTO(105,110,115,120,125) ORD_IND
105   GAMMA  = gamma4
      RHOBAD = rhobad4
      RHOMAX = rhom4
      RHOT   = rhot4
      RHOTUP = rhot6
      RHOI   = rhoi4
      RHOIUP = rhoi6
      SISERR    = SISERR4
      SISERRUP  = SISERR6
      FATERR = faterr4
      VUP    = gamma4*VUP4
      FATDJ0  = fatdj04
      FATDJ0I = fatdj04i
      tolrhoJ0 = tolrhoJ4
      DELTAH2   = deltah2_4
      DELTAH1SF = deltah1_4sf
      CFAT1  = cfat4_1
      CFAT2  = cfat4_2 
      GOTO 140
110   GAMMA  = gamma6
      RHOBAD = rhobad6
      RHOMAX = rhom6
      RHOT   = rhot6
      RHOTUP = rhot8
      RHOI   = rhoi6
      RHOIUP = rhoi8
      SISERR    = SISERR6
      SISERRUP  = SISERR8
      SISERRDN  = SISERR4
      FATERR = faterr6
      VUP    = gamma6*VUP6
      FATDJ0  = fatdj06
      FATDJ0I = fatdj06i
      tolrhoJ0 = tolrhoJ6
      DELTAH2   = deltah2_6
      DELTAH1SF = deltah1_6sf
      CFAT1  = cfat6_1
      CFAT2  = cfat6_2 
      GOTO 140
115   GAMMA  = gamma8
      RHOBAD = rhobad8
      RHOMAX = rhom8
      RHOT   = rhot8
      RHOTUP = rhot10      
      RHOI   = rhoi8
      RHOIUP = rhoi10
      SISERR    = SISERR8
      SISERRUP  = SISERR10
      SISERRDN  = SISERR6
      FATERR = faterr8
      VUP    = gamma8*VUP8
      FATDJ0  = fatdj08
      FATDJ0i = fatdj08i
      tolrhoJ0 = tolrhoJ8
      DELTAH2   = deltah2_8
      DELTAH1SF = deltah1_8sf
      CFAT1  = cfat8_1
      CFAT2  = cfat8_2 
      GOTO 140
120   GAMMA  = gamma10
      RHOBAD = rhobad10
      RHOMAX = rhom10
      RHOT   = rhot10
      RHOTUP = rhot12
      RHOI   = rhoi10
      RHOIUP = rhoi12
      SISERR    = SISERR10
      SISERRUP  = SISERR12
      SISERRDN  = SISERR8
      FATERR = faterr10
      VUP    = gamma10*VUP10
      FATDJ0  = fatdj010
      FATDJ0i = fatdj010i
      tolrhoJ0 = tolrhoJ10
      DELTAH2   = deltah2_10
      DELTAH1SF = deltah1_10sf
      CFAT1  = cfat10_1
      CFAT2  = cfat10_2 
      GOTO 140
125   GAMMA  = gamma12
      RHOBAD = rhobad12
      RHOMAX = rhom12
      RHOT   = rhot12
      RHOI   = rhoi12
      SISERR    = SISERR12
      SISERRDN  = SISERR10
      FATDJ0  = fatdj012
      FATDJ0i = fatdj012i
      tolrhoJ0 = tolrhoJ12
      DELTAH2   = deltah2_12
      DELTAH1SF = deltah1_12sf
      CFAT1  = cfat12_1
      CFAT2  = cfat12_2 

140   CONTINUE

C-------------------------------------------------------------------------------------------------------------------
C     JACOBIAN EVALUATION
      LINEAR = .FALSE.
      ESTIM =(.NOT.SMALLM).AND. 
     &       (SUCCESS.OR.(.NOT.SUCCESS.AND.CALJAC))
      IF (ESTIM) THEN
        SCALJ0_1 = 0D0
        DO I = 1,M
           SCALJ0_1 = DMAX1(SCALJ0_1,DABS(Y0(I)))
        END DO
        SCALJ0_1  = scalJ0*(1d0+SCALJ0_1)
        DO I=1,M
          DELJ0(I) = Y0(I) + EJ0(I)*SCALJ0_1
        END DO
        IERR=0
        CALL FCN(M,T0,DELJ0,FJ0,IERR,RPAR,IPAR)
        NFEVAL = NFEVAL + 1
        NODJ0 = (IERR.NE.0) 
        NJ00 = NJ0
        IF (.NOT. NODJ0) THEN 
            SCALJ0_1 = 1D0/SCALJ0_1
            NJ0 = 0D0
            DO I=1,M
              DELJ0(I) = (FJ0(I)-F0(I))*SCALJ0_1
              NJ0 = DMAX1(NJ0,DABS(DELJ0(I)))
            END DO
        END IF
      END IF

      IF (SUCCESS) THEN 

        CALJAC=CALJAC0
        IF (CALJAC) GOTO 150

        ESTIM1 = ESTIM.AND..NOT.NODJ00.AND..NOT.NODJ0
        IF (ESTIM1) THEN
           DJ0  = 0D0
           DO I=1,M
             DJ0   = DMAX1(DJ0,
     &                   DABS(DELJ0(I)-DELJ00(I)))
           END DO
           CALJAC=(DJ0.GT.1D2*UROUND*NJ00)
           LINEAR =.NOT.CALJAC
           IF (LINEAR) GOTO 150
        END IF

        CALJAC=(K.NE.KOLD).OR.(.NOT.EXTRAP)
        IF (CALJAC) GOTO 150
      
        CALJAC=(RHO.GE.tolrhoJ0).AND.(IT.GE.itmaxJ0)
        IF (.NOT.CALJAC.OR..NOT.ESTIM1) GOTO 150

        NQINF  = .NOT.(QINF.OR.QINFJ)
        CALJAC=.NOT.
     &      ( ((RHO.LT.5D-2) .OR. (IT .LT. 4)) .AND.
     &      ( (QINF.AND.(DJ0 .LE. fatdJ0i*NJ00)).OR.
     &        ((.NOT.QINF).AND.NQINF.AND.(RHOJ0.LT.1D0).AND.
     &         (HJ0*DJ0 .LE. fatdJ0*RHOJ0/RHOT)) ))
      END IF

150   CONTINUE
      TRUEJAC = CALJAC .OR.(.NOT.SUCCESS)
      IF (CALJAC) THEN
           IF (IJAC.EQ.0) THEN
C NUMERICAL JACOBIAN
               DO I=1,M
                 YSAFE=Y0(I)
                 DELT = DSQRT(UROUND*DMAX1(1.D-5,DABS(YSAFE)))
                 Y0(I) = YSAFE+DELT
                 IERR = 0
                 CALL FCN(M,T0,Y0,FJ0,IERR,RPAR,IPAR)
                 IF (IERR.NE.0) THEN
                     WRITE(6,1030) IERR
 1030                FORMAT(/,/,'ERROR CODE IERR = ', I3, /,
     &               'RETURNED BY THE SUBROUTINE FCN DURING',/,
     &               'THE NUMERICAL EVALUATION OF THE JACOBIAN',/,/)
                    WRITE(6,900) T0
                    GOTO 800
                 END IF
                 IF (JBAND) THEN
                    DO J=MAX(1,I-MUJAC),MIN(M,I+MLJAC)
                       J0(J-I+MUJAC+1,I)=(FJ0(J)-F0(J))/DELT 
                    END DO
                 ELSE
                    DO J=1,M
                       J0(J,I)=(FJ0(J)-F0(J))/DELT
                    END DO
                 END IF
                 Y0(I)=YSAFE
               END DO
           ELSE 
C ANALYTICAL JACOBIAN
              IERR = 0
              CALL JAC(M,T0,Y0,J0,LDJAC,IERR,RPAR,IPAR)
              IF (IERR.NE.0) THEN
                 WRITE(6,1040) IERR
 1040            FORMAT(/,/,'ERROR CODE IERR = ', I3, /,
     &           'RETURNED BY THE SUBROUTINE JAC',/,/)
                 WRITE(6,900) T0
                 IDID = -6
                 GOTO 800
              END IF
           END IF
           NJEVAL = NJEVAL + 1
           DO I=1,M
              DELJ00(I) = DELJ0(I)
           END DO
           NJ00   = NJ0
           NODJ00 = NODJ0
           IF (SMALLM) THEN
            IF (JBAND) THEN
             NJ0 = 0D0
             DO I=1,M
                DELT=0D0
                DO J=MAX(1,I-MLJAC),MIN(M,MUJAC+I)
                   DELT=DELT+DABS(J0(MUJAC+1-J+I,J))
                END DO
                NJ0 = DMAX1(NJ0,DELT)
             END DO
            ELSE
              NJ0 = 0D0
              DO I=1,M
                 DELT = 0D0
                 DO J=1,M
                   DELT = DELT +DABS(J0(I,J))
                 END DO
                 NJ0 = DMAX1(NJ0,DELT)
              END DO
            END IF
            NJ0 = NJ0/DBLE(M)
           END IF
      END IF

155   CONTINUE

C-------------------------------------------------------------------------------------------------------------------
C     COMPUTE AND FACTORIZE THE ITERATION MATRIX THETA = (I -h*gamma*J0)

      RATHH = H/HFATT
      CALFACT=((KOLD.NE.K).OR.(.NOT.SUCCESS).OR.CALJAC .OR.
     &        (M.LT.2*K).OR.
     &        (QINF.AND.(DABS(RATHH-1D0).GT.fatdJ0i))
     &         .OR. ((.NOT.QINF).AND.
     &        ((RATHH.GT.DELTAH2).OR.(RATHH.LT.DELTAH1SF)))
     &        .OR. CALFACT0
     &        .OR.((RHO.GT.5D-2).AND.(IT.GE.4)) )
      IF (CALFACT.OR.QINF) GOTO 160
      CALFACT = QINFF.OR.(RHOFATT.GE.1D0)
      IF (CALFACT.OR.(RATHH.GE.1D0).OR.(IT.EQ.1)) GOTO 160
      
      RHONEW = RHO * H/H0
      ITNEW = DBLE(IT)*DLOG(RHO)/DLOG(RHONEW)
      ALFAFATT = (DBLE(M)/(DBLE(6*K)*ITNEW)) + 1D0
      CFAT3= RHOT/(GAMMA*RHOFATT)*
     &            (DELTAH1SF*RHOFATT)**(1D0/ALFAFATT)
      CFAT3=CFAT2-CFAT3**2D0
      DISCR = CFAT1**2D0 - CFAT3
      IF (DISCR .GE. 0D0) THEN
         DELTAH1 = -(CFAT1 + DSQRT(DISCR))
      ELSE
         DELTAH1 = 1D0
      ENDIF
      CALFACT=RATHH.LT.DELTAH1

160   CONTINUE
      IF (CALFACT) THEN      
         NSING = 0
170      CONTINUE
         HGAMMA = H*GAMMA
         IF (JBAND) THEN
               DO J=1,M
                 DO I=1,LDJAC
                    THETA(I+MLJAC,J)=-HGAMMA*J0(I,J)
                 END DO
                 THETA(LDJAC,J)=1D0 +THETA(LDJAC,J)
               END DO
         ELSE
               DO J=1,M
                 DO I=1,M
                    THETA(I,J) = -HGAMMA*J0(I,J)
                 END DO
                 THETA(J,J) = 1D0 + THETA(J,J)
               END DO
         END IF
         CALL DECLU(M,THETA,LDLU,MLJAC,MUJAC,IPVT,IJOB,INFO)
         NLU = NLU + 1
         IF (INFO.NE.0) THEN
               NSING = NSING + 1
               IF (NSING.GT.5) THEN
                  WRITE(*,*) 'MATRIX IS REPEATEDLY SINGULAR, IER= ',INFO
                  WRITE(6,900) T0
                  IDID = -4
                  GOTO 800
               ELSE
                  H = H*.5D0
                  IF (.1d0*DABS(H) .LE. DABS(T0)*UROUND) THEN
                    WRITE(6,*) 'STEPSIZE TOO SMALL, H=',H
                    WRITE(6,900) T0
                    IDID = -3
                    GOTO 800
                  END IF
                  GOTO 170
               END IF
         END IF
      END IF

      CALJAC0  = .FALSE.
      CALFACT0 = .FALSE.

C -----------------------------------------------------------------------------------------------------------------

      DO I=1,K
         T(I)=T0+I*H
      END DO

C-------------------------------------------------------------------------------------------------------------------
C     STOPPING CRITERION WHEN THE SOLUTION HAS SMALL ENTRIES

      FACNEWT = FACNEWTV(ORD_IND)
      MINY0   = DABS(Y0(1))
      MAXF0   = DABS(F0(1))
      INDMIN  = 1
      DO I=2,M
         IF (DABS(Y0(I)).LT.MINY0) THEN
             MINY0  = DABS(Y0(I))
             INDMIN = I
         ENDIF
         MAXF0 = DMAX1(MAXF0,DABS(F0(I)))
      END DO
      IF (((MINY0 .LT. TOLMINY0).AND.(DABS(F0(INDMIN)).LT.1d-1)))
     &THEN
          FACNEWT = DMIN1(1D-1,FACNEWT)
          IF ((DABS(F0(INDMIN)).LT.TOLMINF0).AND.(MAXF0.LT.TOLMAXF0))
     &        FACNEWT=DMIN1(FACNSMALL,FACNEWT)
      END IF 
      IF (RESTRICT) FACNEWT=DMIN1(FACNEWT,FACNRESTR)

C-------------------------------------------------------------------------------------------------------------------
C     SOLUTION INITIALIZATION      
      IF ( (.NOT.EXTRAP).OR.(NFAILCONV.GT.flmx) ) THEN
C     THE SOLUTION IS INITIALIZED WITH THE CONSTANT PROFILE
         IF (NFAILCONV.GT.flhlt) THEN
             WRITE(*,*) 'TOO MANY CONSECUTIVE NEWTON FAILURES'
             WRITE(6,900) T0
             IDID = -5
             GOTO 800
             RETURN
         END IF
         DO J = 1,K
            DO I=1,M
               Y(I,J)=Y0(I)
            END DO
         END DO
      ELSE
         CALL EXTRAPOLA(M,K0,K,H0,H,Y,DD)
      END IF

C-------------------------------------------------------------------------------------------------------------------
C     NEWTON ITERATION

      IT       = 0
      RHO      = 0D0
      NERR0    = 1D0
      NERRSTOP  = ATOL*DMAX1(FACNEWT,UROUND/RTOL)
200   CONTINUE
      IERR0=0
      DO I=1,K
        IERR = 0
        CALL FCN(M,T(I),Y(1,I),F(1,I),IERR,RPAR,IPAR)
        IERR0 = IERR0+IERR
      END DO
      NFEVAL = NFEVAL + K
      IF (IERR.NE.0) THEN
         NERR = NERRSTOP + 1D0
         JVAI = .FALSE.
         GOTO 305
      END IF

      GOTO (210,220,230,240,250) ORD_IND

210   CALL blendstep4(M,Y0,F0,Y,F,H,THETA,IPVT,ERR,
     &                LDLU,MLJAC,MUJAC,IJOB)      
      GOTO 300

220   CALL blendstep6(M,Y0,F0,Y,F,H,THETA,IPVT,ERR,
     &                LDLU,MLJAC,MUJAC,IJOB)      
      GOTO 300

230   CALL blendstep8(M,Y0,F0,Y,F,H,THETA,IPVT,ERR,
     &                LDLU,MLJAC,MUJAC,IJOB)      
      GOTO 300

240   CALL blendstep10(M,Y0,F0,Y,F,H,THETA,IPVT,ERR,
     &                 LDLU,MLJAC,MUJAC,IJOB)      
      GOTO 300

250   CALL blendstep12(M,Y0,F0,Y,F,H,THETA,IPVT,ERR,
     &                 LDLU,MLJAC,MUJAC,IJOB)      
c      GOTO 300

300   CONTINUE

      IT = IT + 1

      CALL NORM(M,K,SCAL,ERR,NERR,NERRUP)
c CHECK FOR NANs 
      IF (ISNAN(NERR)) THEN
         NERR = NERRSTOP + 1D0
         JVAI = .FALSE.
         GOTO 305
      END IF

C SPECTRAL RADIUS ESTIMATE

      NERROLD = NERR0
      NERR0   = NERR
      RHO0    = RHO
      RHO     = NERR0/NERROLD
      IF (IT.GT.2) RHO = DSQRT(RHO0*RHO)

      JVAI = (NERR .GT. NERRSTOP) .AND. (IT .LE. MAXIT) 
     &       .AND.(IT .LE. 2 .OR. RHO .LE. RHOBAD) 
305   IF (JVAI)  GOTO 200

C END OF NEWTON'S ITERATION
C-------------------------------------------------------------------------------------------------------------------
      NLINSYS  = NLINSYS  + 2*IT*K
      NITER(ORD_IND) = NITER(ORD_IND) + IT
      NSTEP(ORD_IND) = NSTEP(ORD_IND) + 1

      IF (NERR.GT.NERRSTOP) GOTO 308

      IERR0 = 0
      DO I=1,K
        IERR = 0
        CALL FCN(M,T(I),Y(1,I),F(1,I),IERR,RPAR,IPAR)
        IERR0 = IERR0 + IERR
      END DO
      NFEVAL = NFEVAL + K
      IF (IERR0.NE.0)  NERR = NERRSTOP + 1D0

308   CONTINUE

      IF (NERR .GT. NERRSTOP) THEN
C     NEWTON HAS FAILED
          NFAILNEWT(ORD_IND) = NFAILNEWT(ORD_IND) + 1
          NFAILCONV  = NFAILCONV  + 1
          NFAILCONS  = MAX(NFAILCONS + 1,2)
          H       = facnocon * H
          KOLD    = K
          IF (ORD .GT. ORDMIN) THEN
              ORD     = ORD - 2
              ORD_IND = INT(ORD/2) - 1
              K       = STEP_ORD(ORD_IND)
          END IF
          CALJAC   = .NOT.TRUEJAC
          TRUEJAC  = .TRUE.
          SUCCESS  = .FALSE.
          EXTRAP   = .FALSE.
          LAST     = .FALSE.
          RESTRICT = .FALSE.
          GOTO 550
      END IF

C-------------------------------------------------------------------------------------------------------------------
C     LOCAL ERROR ESTIMATION

      GOTO (310,320,330,340,350) ORD_IND
310   CALL localerr4(M,F0,F,H,ERR,SCAL,NERR,NERRUP,NLINSYS,
     &               THETA,IPVT,LDLU,MLJAC,MUJAC,IJOB) 
      GOTO 400

320   CALL localerr6(M,F0,F,H,ERR,SCAL,NERR,NERRUP,NLINSYS,
     &               THETA,IPVT,LDLU,MLJAC,MUJAC,IJOB) 
      GOTO 400

330   CALL localerr8(M,F0,F,H,ERR,SCAL,NERR,NERRUP,NLINSYS,
     &               THETA,IPVT,LDLU,MLJAC,MUJAC,IJOB) 
      GOTO 400

340   CALL localerr10(M,F0,F,H,ERR,SCAL,NERR,NERRUP,NLINSYS,
     &               THETA,IPVT,LDLU,MLJAC,MUJAC,IJOB) 
      GOTO 400

350   CALL localerr12(M,F0,F,H,ERR,SCAL,NERR,NERRUP,NLINSYS,
     &               THETA,IPVT,LDLU,MLJAC,MUJAC,IJOB) 
c      GOTO 400

400   CONTINUE
      IF (ISNAN(NERR)) THEN 
          NERR=NERRSTOP + 1D0
          GOTO 308
      END IF
      NFAILCONV = 0

      IF (NERR .GT. 0d0) THEN
          HNEW = H*(SFTY*ATOL/NERR)**ESP
      ELSE
          HNEW = facr*H
      END IF
      
      IF (NERR .GE. ATOL) THEN
C     FAILURE DUE TO LOCAL ERROR TEST
          HNEW = H*(1D-1*ATOL/NERR)**ESP
          NFAILERR(ORD_IND) = NFAILERR(ORD_IND) + 1
          NFAILCONS = MAX(NFAILCONS + 1,2)
          CALJAC   = .NOT. TRUEJAC
          TRUEJAC  = .TRUE.
          SUCCESS  = .FALSE.
          EXTRAP   =  EXTRAP0
          LAST     = .FALSE.
          RESTRICT = .FALSE.
          KOLD     = K
          IF (ERROR.AND.(ORD.GT.ORDMIN)) THEN
             ORD = ORD - 2
             ORD_IND = INT(ORD/2) - 1
             K = STEP_ORD(ORD_IND)
             H = DMAX1(HNEW/2D0,facl*H)
          ELSE
             H  = DMAX1(HNEW,facl*H)
          END IF
          ERROR = .TRUE.
          GOTO 550 
      END IF

C     STEP ACCEPTED
      NACCEPT(ORD_IND) = NACCEPT(ORD_IND) + 1

      IF (IOUT.EQ.1)
     &   CALL SOLOUT(M,T,Y,F,K,ORD,IRTRN)

      IF (LAST) GOTO 600

      SUCCESS   = .TRUE.
      NORD      =  NORD + 1
      NFAILCONS = MAX(NFAILCONS - 1,0)

      EXTRAP0 = .FALSE.
      MAXF0   = 0D0
      DO I=1,M    
         MAXDELTA = DABS(Y(I,K)-Y0(I))*SCALEXTRAP(I)
         EXTRAP0 = EXTRAP0 .OR.
     &    ((dabs(Y0(I)).LE.1D-1).AND.(MAXDELTA.GE.TOLESTRAP(2))).OR.
     &    ((dabs(Y0(I)).GT.1D-1).AND.(MAXDELTA.GE.TOLESTRAP(1)))
         MAXF0 = DMAX1(MAXF0,DABS(F(I,K)))
      END DO
      EXTRAP0  = EXTRAP0 .OR. (MAXF0.GE.5D-1)
      RESTRICT = .NOT.EXTRAP0
      EXTRAP   = EXTRAP0
      IF (EXTRAP) CALL DIFFDIV(M,K,Y0,Y,DD)

C     ORDER VARIATION
      RATH   = HNEW/H 
      IF (IT .LT. 2) RHO=1D-3
      RATRHO = RHOOLD/RHO
      KNEW   = K
      ORDNEW = ORD
      
      QINF   = (NERRUP.GE.NERR).AND.(NERR.NE.0D0)
      IF (TRUEJAC) THEN
         HJ0   = H
         RHOJ0 = RHO
         IF (IT.LE.1) RHOJ0 = 1D0
         QINFJ = QINF
      END IF 
      IF (CALFACT) THEN
         HFATT   = H
         RHOFATT = RHO
         IF (IT.LE.1) RHOFATT = 1D0
         QINFF   = QINF
      END IF
      
      IF ((ERROR).AND.(ORD.GT.ORDMIN)) THEN
         NERROR = MIN(NERROR+1,2)
         IF (NERROR.GT.1) THEN
            ORDNEW  = ORD - 2
            ORD_IND = ORD_IND - 1
            KNEW    = STEP_ORD(ORD_IND)
            HNEW    = DMIN1(H,HNEW)/2D0
            CALJAC  = .TRUE.
            TRUEJAC = .TRUE.
            GOTO 500
         END IF
      ELSE
         NERROR = 0
      END IF
      ERROR = .FALSE.

      IF (RHO.LE.RHOML) GOTO 445
c     SPECTRAL RADIUS TOO LARGE
      IF ((ORD.EQ.ORDMIN).OR.(IT.LE.3)) GOTO 500

405   CONTINUE      

      GOTO (410,415,420,425) ORD_IND-1
410   CALL errdown6(M,F0,F,H,ERR(1,2),SCAL,NERRDOWN,NLINSYS,
     &              QINF,THETA,IPVT,LDLU,MLJAC,MUJAC,IJOB) 
      GOTO 440

415   CALL errdown8(M,F0,F,H,ERR(1,2),SCAL,NERRDOWN,NLINSYS,
     &              QINF,THETA,IPVT,LDLU,MLJAC,MUJAC,IJOB) 
      GOTO 440

420   CALL errdown10(M,F0,F,H,ERR(1,2),SCAL,NERRDOWN,NLINSYS,
     &               QINF,THETA,IPVT,LDLU,MLJAC,MUJAC,IJOB) 
      GOTO 440

425   CALL errdown12(M,F0,F,H,ERR(1,2),SCAL,NERRDOWN,NLINSYS,
     &               QINF,THETA,IPVT,LDLU,MLJAC,MUJAC,IJOB) 
c      GOTO 440

440   CONTINUE

      IF (NERRDOWN .GT. 0D0) THEN
         HNDN = H*(SFTYDN*ATOL/NERRDOWN)**ESPDN
      ELSE
         HNDN = HNEW
      END IF

      IF (QINF.AND.(HNDN.LT.HNEW)) GOTO 500
      IF (.NOT. QINF) HNDN = DMIN1(HNDN,HNEW)
      ORDNEW  = ORD - 2
      ORD_IND = ORD_IND - 1
      KNEW    = STEP_ORD(ORD_IND)
      HNEW    = HNDN
      GOTO 500

445   CONTINUE

      IF ((NFAILCONS.GT.0).OR.(ORD.GE.ORDMAX)
     &    .OR. (NORD.LT.3) 
     &    .OR.(RATH.GE.FACU1).OR.(RATH.LE.facu2))
     &GOTO 500 

      IF ((NERRUP.LT.1D1*UROUND).AND.(.NOT.QINF).AND.
     &    (IT.LT.3).AND.(RHO.LT.1D1*UROUND))
     &THEN
C       THE PROBLEM IS A QUADRATURE
        CALL ERRUP(M,K,ORD_IND,ERR,H,H0,H00,VUP,NERRUP,SCAL,
     &             THETA,IPVT,LDLU,MLJAC,MUJAC,IJOB)             
        NLINSYS = NLINSYS + 1
        IF (NERRUP .GT. 0D0) THEN
           HNUP = H*(SFTYUP*ATOL/NERRUP)**ESPUP
        ELSE
           HNUP = facr*H
        END IF
        NU   = DBLE(IT)
        NUUP = DBLE(IT)
        KUP  = STEP_ORD(ORD_IND+1)
        FI   = (HNEW/HNUP)*(DBLE(K)/DBLE(KUP))*
     &          (CSIS*(2D0*KUP*NUUP+SISERRUP)+CFACT)/
     &          (CSIS*(2D0*K*NU+SISERR)+CFACT)
        IF (FI.LT.1D0) THEN
         ORDNEW  = ORD + 2
         ORD_IND = ORD_IND + 1
         KNEW    = STEP_ORD(ORD_IND)
         HNEW    = HNUP
        END IF
        GOTO 500
      END IF

      IF (QINF) GOTO 450
      
      IF (RHO.GE.RHOMU) THEN
         STAGNA = (DABS(RATH-1D0).LT.1D-1).AND.
     &            (DABS(RATRHO-1D0).LT.1D-1)
         IF (.NOT.STAGNA.OR.(T(K).EQ.0D0)) GOTO 450
         IERR = 0
         CALL FCN(M,T(K)*(1D0+1D-5),Y(1,K),FJ0,IERR,RPAR,IPAR)
         NFEVAL = NFEVAL + 1
         IF (IERR.NE.0) GOTO 450
         NF0 = 0d0
         NF  = 0d0
         DO I=1,M
            NF0 = DMAX1(NF0,DABS(FJ0(I)-F(I,K)))
            NF  = DMAX1(NF,DABS(F(I,K)))
         END DO
         NF0 = NF0/(1d-5*T(K))    
         IF (NJ0*NF.GE.2D0*NF0) GOTO 450
      END IF
 
      IF (NERRUP .GT. 0D0) THEN
         HNUP = H*(SFTYUP*ATOL/NERRUP)**ESPUP
      ELSE
         HNUP = facr*H
      END IF

      NU   = (DBLE(IT) *DLOG(RHO))/DLOG(RHO*RATH)
      NUUP = (DBLE(IT) *DLOG(RHO))/DLOG(RHO*(RHOTUP/RHOT)*(HNUP/H)) 

      IF ((NUUP.LT.0D0).OR.(NU.LT.0D0)) GOTO 450

      KUP = STEP_ORD(ORD_IND+1)
      FI  = (HNEW/HNUP)*(DBLE(K)/DBLE(KUP))*
     &      (CSIS*(2D0*KUP*NUUP+SISERRUP)+CFACT)/
     &      (CSIS*(2D0*K*NU+SISERR)+CFACT)

      IF (FI.LT.1D0)
     &THEN
         ORDNEW  = ORD + 2
         ORD_IND = ORD_IND + 1
         KNEW    = STEP_ORD(ORD_IND)
         HNEW    = HNUP
         GOTO 500
      END IF
 
450   CONTINUE
C ---------------------------------------------------------------
C ORDER REDUCTION RECOVERY
C ---------------------------------------------------------------
      STAGNA = (RATH .GT. rath1) .AND. (RATH .LT. rath2) .AND.
     &         (RATRHO .GT. ratrho1).AND.(RATRHO.LT.ratrho2)
      IF ( QINF .OR.
     &     (STAGNA.AND.(FATERR*NERRUP.GE.NERR)) )
     &THEN

             IF ((.NOT.CALFACT.OR.(.NOT.TRUEJAC.AND..NOT.LINEAR))
     &          .AND.QINF) THEN
c IN THIS CASE, THE SPECTRAL RADIUS DOESN'T BEHAVES LIKE rho=rhoti/q.
                CALJAC0  = .NOT.LINEAR
                CALFACT0 = .TRUE.
                GOTO 500
             END IF

             CALL ERRUP(M,K,ORD_IND,ERR,H,H0,H00,VUP,NERRUP1,SCAL,
     &                  THETA,IPVT,LDLU,MLJAC,MUJAC,IJOB)             
             NLINSYS = NLINSYS + 1
            
             IF (NERRUP1 .GT. 0D0) THEN
                 HNUP1 = H*(SFTYUP*ATOL/NERRUP1)**ESPUP
             ELSE
                 HNUP1 = facr*H
             END IF
             
             IF (QINF.AND.(HNUP1.GT.H)) THEN
                  RHOUP=RHO*(RHOIUP/RHOI)*(H/HNUP1)
                  IF (RHOUP.GT.RHOMLV(ORD_IND+1)) GOTO 500
             END IF

             NU1    = (DBLE(IT)*DLOG(RHO))/DLOG(RHO/RATH) 
             NUUP1  = (DBLE(IT)*DLOG(RHO))/
     &                DLOG(RHO*(H/HNUP1)*(RHOIUP/RHOI)) 
             IF ((NUUP1.LT.0D0).OR.(NU1.LT.0D0)) GOTO 500

             KUP    = STEP_ORD(ORD_IND+1)
             FI1    = (HNEW/HNUP1)*(DBLE(K)/DBLE(KUP))*
     &                (CSIS*(2D0*KUP*NUUP1+SISERRUP)+CFACT)/
     &                (CSIS*(2D0*K*NU1+SISERR)+CFACT)
             IF   (FI1 .LT. 1D0) 
     &       THEN
                   ORDNEW  = ORD + 2
                   ORD_IND = ORD_IND + 1
                   KNEW    = STEP_ORD(ORD_IND)
                   HNEW    = HNUP1
             END IF
      END IF

500   CONTINUE

      DO I=1,M
         Y0(I)   = Y(I,K)
         ABSY0   = DABS(Y0(I))
         F0(I)   = F(I,K)
         SCAL(I) = 1D0/(1d0+RTOLATOL*ABSY0)
         SCALEXTRAP(I) = 1d0/(1d0+ABSY0)
         IF (ORD.LT.ORDMAX) THEN
            ERR(I,K+2) = ERR(I,K+1)
            ERR(I,K+1) = ERR(I,1)
         END IF
      END DO

      T0   = T(K)

      K0   = K
      KOLD = K
      K    = KNEW
      ORD  = ORDNEW

      H00 = H0
      H0  = H
      IF (NFAILCONS.GT.0) HNEW = DMIN1(H,HNEW)
      HNEW = DMAX1(HNEW,facl*H)
      H    = DMIN1(HNEW,facr*H,HMAX)

      RHOOLD = RHO

550   CONTINUE

      NSTEPS = NSTEPS + 1

      IF (.1d0*DABS(T0-TEND)/DBLE(K).GT.dabs(T0)*UROUND) THEN
         IF (.1d0*DABS(H) .LE. DABS(T0)*UROUND) THEN
            WRITE(6,*) 'STEPSIZE TOO SMALL, H=',H
            WRITE(6,900) T0
            IDID = -3
            GOTO 800
         END IF
         IF (NSTEPS .GE. MAXSTEP) THEN
            WRITE(6,*) 'MORE THAN NMAX = ',MAXSTEP, ' STEPS ARE NEEDED'
            WRITE(6,900) T0
            IDID = -2
            GOTO 800
         END IF
	   IF ((T0 + DBLE(K)*H) .GT. TEND) THEN
            H = (TEND -T0)/DBLE(K)
            LAST = .TRUE.
         END IF

         GOTO 100 
C     END OF MAIN LOOP      
      ELSE
         IDID = 0
      END IF 
            
C     INTEGRATION END
600   CONTINUE
      DO I=1,M
          Y0(I) = Y(I,K)
          F0(I) = F(I,K)
      END DO
      T0 = T(K)

      RTOL = RTOL0
      
900   FORMAT(' EXIT AT T = ',D18.4)
800   RETURN

      END

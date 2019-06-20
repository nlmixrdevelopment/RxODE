*DECK DUMACH
      DOUBLE PRECISION FUNCTION DUMACH ()
C***BEGIN PROLOGUE  DUMACH
C***PURPOSE  Compute the unit roundoff of the machine.
C***CATEGORY  R1
C***TYPE      DOUBLE PRECISION (RUMACH-S, DUMACH-D)
C***KEYWORDS  MACHINE CONSTANTS
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C *Usage:
C        DOUBLE PRECISION  A, DUMACH
C        A = DUMACH()
C
C *Function Return Values:
C     A : the unit roundoff of the machine.
C
C *Description:
C     The unit roundoff is defined as the smallest positive machine
C     number u such that  1.0 + u .ne. 1.0.  This is computed by DUMACH
C     in a machine-independent manner.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DUMSUM
C***REVISION HISTORY  (YYYYMMDD)
C   19930216  DATE WRITTEN
C   19930818  Added SLATEC-format prologue.  (FNF)
C   20030707  Added DUMSUM to force normal storage of COMP.  (ACH)
C***END PROLOGUE  DUMACH
C
      DOUBLE PRECISION U, COMP
C***FIRST EXECUTABLE STATEMENT  DUMACH
      U = 1.0D0
 10   U = U*0.5D0
      CALL DUMSUM(1.0D0, U, COMP)
      IF (COMP .NE. 1.0D0) GO TO 10
      DUMACH = U*2.0D0
      RETURN
C----------------------- End of Function DUMACH ------------------------
      END
      SUBROUTINE DUMSUM(A,B,C)
C     Routine to force normal storing of A + B, for DUMACH.
      DOUBLE PRECISION A, B, C
      C = A + B
      RETURN
      END
*DECK DCFODE
      SUBROUTINE DCFODE (METH, ELCO, TESCO)
C***BEGIN PROLOGUE  DCFODE
C***SUBSIDIARY
C***PURPOSE  Set ODE integrator coefficients.
C***TYPE      DOUBLE PRECISION (SCFODE-S, DCFODE-D)
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C
C  DCFODE is called by the integrator routine to set coefficients
C  needed there.  The coefficients for the current method, as
C  given by the value of METH, are set for all orders and saved.
C  The maximum order assumed here is 12 if METH = 1 and 5 if METH = 2.
C  (A smaller value of the maximum order is also allowed.)
C  DCFODE is called once at the beginning of the problem,
C  and is not called again unless and until METH is changed.
C
C  The ELCO array contains the basic method coefficients.
C  The coefficients el(i), 1 .le. i .le. nq+1, for the method of
C  order nq are stored in ELCO(i,nq).  They are given by a genetrating
C  polynomial, i.e.,
C      l(x) = el(1) + el(2)*x + ... + el(nq+1)*x**nq.
C  For the implicit Adams methods, l(x) is given by
C      dl/dx = (x+1)*(x+2)*...*(x+nq-1)/factorial(nq-1),    l(-1) = 0.
C  For the BDF methods, l(x) is given by
C      l(x) = (x+1)*(x+2)* ... *(x+nq)/K,
C  where         K = factorial(nq)*(1 + 1/2 + ... + 1/nq).
C
C  The TESCO array contains test constants used for the
C  local error test and the selection of step size and/or order.
C  At order nq, TESCO(k,nq) is used for the selection of step
C  size at order nq - 1 if k = 1, at order nq if k = 2, and at order
C  nq + 1 if k = 3.
C
C***SEE ALSO  DLSODE
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791129  DATE WRITTEN
C   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
C   890503  Minor cosmetic changes.  (FNF)
C   930809  Renamed to allow single/double precision versions. (ACH)
C***END PROLOGUE  DCFODE
C**End
      INTEGER METH
      INTEGER I, IB, NQ, NQM1, NQP1
      DOUBLE PRECISION ELCO, TESCO
      DOUBLE PRECISION AGAMQ, FNQ, FNQM1, PC, PINT, RAGQ,
     1   RQFAC, RQ1FAC, TSIGN, XPIN
      DIMENSION ELCO(13,12), TESCO(3,12)
      DIMENSION PC(12)
C
C***FIRST EXECUTABLE STATEMENT  DCFODE
C     GO TO (100, 200), METH
      IF (METH .EQ. 1) GO TO 100
      IF (METH .EQ. 2) GO TO 200
C
 100  ELCO(1,1) = 1.0D0
      ELCO(2,1) = 1.0D0
      TESCO(1,1) = 0.0D0
      TESCO(2,1) = 2.0D0
      TESCO(1,2) = 1.0D0
      TESCO(3,12) = 0.0D0
      PC(1) = 1.0D0
      RQFAC = 1.0D0
      DO 140 NQ = 2,12
C-----------------------------------------------------------------------
C The PC array will contain the coefficients of the polynomial
C     p(x) = (x+1)*(x+2)*...*(x+nq-1).
C Initially, p(x) = 1.
C-----------------------------------------------------------------------
        RQ1FAC = RQFAC
        RQFAC = RQFAC/NQ
        NQM1 = NQ - 1
        FNQM1 = NQM1
        NQP1 = NQ + 1
C Form coefficients of p(x)*(x+nq-1). ----------------------------------
        PC(NQ) = 0.0D0
        DO 110 IB = 1,NQM1
          I = NQP1 - IB
C110      PC(I) = PC(I-1) + FNQM1*PC(I)
          PC(I) = PC(I-1) + FNQM1*PC(I)
 110      CONTINUE
        PC(1) = FNQM1*PC(1)
C Compute integral, -1 to 0, of p(x) and x*p(x). -----------------------
        PINT = PC(1)
        XPIN = PC(1)/2.0D0
        TSIGN = 1.0D0
        DO 120 I = 2,NQ
          TSIGN = -TSIGN
          PINT = PINT + TSIGN*PC(I)/I
C120      XPIN = XPIN + TSIGN*PC(I)/(I+1)
          XPIN = XPIN + TSIGN*PC(I)/(I+1)
 120      CONTINUE
C Store coefficients in ELCO and TESCO. --------------------------------
        ELCO(1,NQ) = PINT*RQ1FAC
        ELCO(2,NQ) = 1.0D0
        DO 130 I = 2,NQ
C130      ELCO(I+1,NQ) = RQ1FAC*PC(I)/I
          ELCO(I+1,NQ) = RQ1FAC*PC(I)/I
 130      CONTINUE
        AGAMQ = RQFAC*XPIN
        RAGQ = 1.0D0/AGAMQ
        TESCO(2,NQ) = RAGQ
        IF (NQ .LT. 12) TESCO(1,NQP1) = RAGQ*RQFAC/NQP1
        TESCO(3,NQM1) = RAGQ
 140    CONTINUE
      RETURN
C
 200  PC(1) = 1.0D0
      RQ1FAC = 1.0D0
      DO 230 NQ = 1,5
C-----------------------------------------------------------------------
C The PC array will contain the coefficients of the polynomial
C     p(x) = (x+1)*(x+2)*...*(x+nq).
C Initially, p(x) = 1.
C-----------------------------------------------------------------------
        FNQ = NQ
        NQP1 = NQ + 1
C Form coefficients of p(x)*(x+nq). ------------------------------------
        PC(NQP1) = 0.0D0
        DO 210 IB = 1,NQ
          I = NQ + 2 - IB
C210      PC(I) = PC(I-1) + FNQ*PC(I)
          PC(I) = PC(I-1) + FNQ*PC(I)
 210      CONTINUE
        PC(1) = FNQ*PC(1)
C Store coefficients in ELCO and TESCO. --------------------------------
        DO 220 I = 1,NQP1
C220      ELCO(I,NQ) = PC(I)/PC(2)
          ELCO(I,NQ) = PC(I)/PC(2)
 220      CONTINUE
        ELCO(2,NQ) = 1.0D0
        TESCO(1,NQ) = RQ1FAC
        TESCO(2,NQ) = NQP1/ELCO(1,NQ)
        TESCO(3,NQ) = (NQ+2)/ELCO(1,NQ)
        RQ1FAC = RQ1FAC/FNQ
 230    CONTINUE
      RETURN
C----------------------- END OF SUBROUTINE DCFODE ----------------------
      END
*DECK DINTDY
      SUBROUTINE DINTDY (T, K, YH, NYH, DKY, IFLAG)
C***BEGIN PROLOGUE  DINTDY
C***SUBSIDIARY
C***PURPOSE  Interpolate solution derivatives.
C***TYPE      DOUBLE PRECISION (SINTDY-S, DINTDY-D)
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C
C  DINTDY computes interpolated values of the K-th derivative of the
C  dependent variable vector y, and stores it in DKY.  This routine
C  is called within the package with K = 0 and T = TOUT, but may
C  also be called by the user for any K up to the current order.
C  (See detailed instructions in the usage documentation.)
C
C  The computed values in DKY are gotten by interpolation using the
C  Nordsieck history array YH.  This array corresponds uniquely to a
C  vector-valued polynomial of degree NQCUR or less, and DKY is set
C  to the K-th derivative of this polynomial at T.
C  The formula for DKY is:
C               q
C   DKY(i)  =  sum  c(j,K) * (T - tn)**(j-K) * h**(-j) * YH(i,j+1)
C              j=K
C  where  c(j,K) = j*(j-1)*...*(j-K+1), q = NQCUR, tn = TCUR, h = HCUR.
C  The quantities  nq = NQCUR, l = nq+1, N = NEQ, tn, and h are
C  communicated by COMMON.  The above sum is done in reverse order.
C  IFLAG is returned negative if either K or T is out of bounds.
C
C***SEE ALSO  DLSODE
C***ROUTINES CALLED  XERRWD
C***COMMON BLOCKS    DLS001
C***REVISION HISTORY  (YYMMDD)
C   791129  DATE WRITTEN
C   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
C   890503  Minor cosmetic changes.  (FNF)
C   930809  Renamed to allow single/double precision versions. (ACH)
C   010418  Reduced size of Common block /DLS001/. (ACH)
C   031105  Restored 'own' variables to Common block /DLS001/, to
C           enable interrupt/restart feature. (ACH)
C   050427  Corrected roundoff decrement in TP. (ACH)
C***END PROLOGUE  DINTDY
C**End
      INTEGER K, NYH, IFLAG
      DOUBLE PRECISION T, YH, DKY
      DIMENSION YH(NYH,*), DKY(*)
      INTEGER IOWND, IOWNS,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     2   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      DOUBLE PRECISION ROWNS,
     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      COMMON /DLS001/ ROWNS(209),
     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     2   IOWND(6), IOWNS(6),
     3   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     4   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     5   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      INTEGER I, IC, J, JB, JB2, JJ, JJ1, JP1
      DOUBLE PRECISION C, R, S, TP
      CHARACTER (LEN=80) MSG
C
C***FIRST EXECUTABLE STATEMENT  DINTDY
      IFLAG = 0
      IF (K .LT. 0 .OR. K .GT. NQ) GO TO 80
      TP = TN - HU -  100.0D0*UROUND*SIGN(ABS(TN) + ABS(HU), HU)
      IF ((T-TP)*(T-TN) .GT. 0.0D0) GO TO 90
C
      S = (T - TN)/H
      IC = 1
      IF (K .EQ. 0) GO TO 15
      JJ1 = L - K
      DO 10 JJ = JJ1,NQ
C10     IC = IC*JJ
        IC = IC*JJ
 10     CONTINUE
 15   C = IC
      DO 20 I = 1,N
C20     DKY(I) = C*YH(I,L)
        DKY(I) = C*YH(I,L)
 20     CONTINUE
      IF (K .EQ. NQ) GO TO 55
      JB2 = NQ - K
      DO 50 JB = 1,JB2
        J = NQ - JB
        JP1 = J + 1
        IC = 1
        IF (K .EQ. 0) GO TO 35
        JJ1 = JP1 - K
        DO 30 JJ = JJ1,J
C30       IC = IC*JJ
          IC = IC*JJ
 30       CONTINUE
 35     C = IC
        DO 40 I = 1,N
C40       DKY(I) = C*YH(I,JP1) + S*DKY(I)
          DKY(I) = C*YH(I,JP1) + S*DKY(I)
 40       CONTINUE
 50     CONTINUE
      IF (K .EQ. 0) RETURN
 55   R = H**(-K)
      DO 60 I = 1,N
C60     DKY(I) = R*DKY(I)
        DKY(I) = R*DKY(I)
 60     CONTINUE
      RETURN
C
 80   MSG = 'DINTDY-  K (=I1) illegal      '
c     CALL XERRWD (MSG, 30, 51, 0, 1, K, 0, 0, 0.0D0, 0.0D0)
      CALL XERRWD (51, 0)
      IFLAG = -1
      RETURN
 90   MSG = 'DINTDY-  T (=R1) illegal      '
c     CALL XERRWD (MSG, 30, 52, 0, 0, 0, 0, 1, T, 0.0D0)
c     MSG='      T not in interval TCUR - HU (= R1) to TCUR (=R2)      '
c     CALL XERRWD (MSG, 60, 52, 0, 0, 0, 0, 2, TP, TN)
      CALL XERRWD (52, 0)
      IFLAG = -2
      RETURN
C----------------------- END OF SUBROUTINE DINTDY ----------------------
      END
*DECK DSOLSY
      SUBROUTINE DSOLSY (WM, IWM, X)
C     SUBROUTINE DSOLSY (WM, IWM, X, TEM)
C***BEGIN PROLOGUE  DSOLSY
C***SUBSIDIARY
C***PURPOSE  ODEPACK linear system solver.
C***TYPE      DOUBLE PRECISION (SSOLSY-S, DSOLSY-D)
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C
C  This routine manages the solution of the linear system arising from
C  a chord iteration.  It is called if MITER .ne. 0.
C  If MITER is 1 or 2, it calls DGESL to accomplish this.
C  If MITER = 3 it updates the coefficient h*EL0 in the diagonal
C  matrix, and then computes the solution.
C  If MITER is 4 or 5, it calls DGBSL.
C  Communication with DSOLSY uses the following variables:
C  WM    = real work space containing the inverse diagonal matrix if
C          MITER = 3 and the LU decomposition of the matrix otherwise.
C          Storage of matrix elements starts at WM(3).
C          WM also contains the following matrix-related data:
C          WM(1) = SQRT(UROUND) (not used here),
C          WM(2) = HL0, the previous value of h*EL0, used if MITER = 3.
C  IWM   = integer work space containing pivot information, starting at
C          IWM(21), if MITER is 1, 2, 4, or 5.  IWM also contains band
C          parameters ML = IWM(1) and MU = IWM(2) if MITER is 4 or 5.
C  X     = the right-hand side vector on input, and the solution vector
C          on output, of length N.
C  TEM   = vector of work space of length N, not used in this version.
C  IERSL = output flag (in COMMON).  IERSL = 0 if no trouble occurred.
C          IERSL = 1 if a singular matrix arose with MITER = 3.
C  This routine also uses the COMMON variables EL0, H, MITER, and N.
C
C***SEE ALSO  DLSODE
C***ROUTINES CALLED  DGBSL, DGESL
C***COMMON BLOCKS    DLS001
C***REVISION HISTORY  (YYMMDD)
C   791129  DATE WRITTEN
C   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
C   890503  Minor cosmetic changes.  (FNF)
C   930809  Renamed to allow single/double precision versions. (ACH)
C   010418  Reduced size of Common block /DLS001/. (ACH)
C   031105  Restored 'own' variables to Common block /DLS001/, to
C           enable interrupt/restart feature. (ACH)
C***END PROLOGUE  DSOLSY
C**End
      INTEGER IWM
      DOUBLE PRECISION WM, X
      DIMENSION WM(*), IWM(*), X(*)
C     DOUBLE PRECISION WM, X, TEM
C     DIMENSION WM(*), IWM(*), X(*), TEM(*)
      INTEGER IOWND, IOWNS,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     2   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      DOUBLE PRECISION ROWNS,
     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      COMMON /DLS001/ ROWNS(209),
     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     2   IOWND(6), IOWNS(6),
     3   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     4   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     5   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      INTEGER I, MEBAND, ML, MU
      DOUBLE PRECISION DI, HL0, PHL0, R
C
C***FIRST EXECUTABLE STATEMENT  DSOLSY
      IERSL = 0
C     GO TO (100, 100, 300, 400, 400), MITER
      IF (MITER .EQ. 1) GO TO 100
      IF (MITER .EQ. 2) GO TO 100
      IF (MITER .EQ. 3) GO TO 300
      IF (MITER .EQ. 4) GO TO 400
      IF (MITER .EQ. 5) GO TO 400
 100  CALL DGESL (WM(3), N, N, IWM(21), X, 0)
      RETURN
C
 300  PHL0 = WM(2)
      HL0 = H*EL0
      WM(2) = HL0
      IF (HL0 .EQ. PHL0) GO TO 330
      R = HL0/PHL0
      DO 320 I = 1,N
        DI = 1.0D0 - R*(1.0D0 - 1.0D0/WM(I+2))
        IF (ABS(DI) .EQ. 0.0D0) GO TO 390
C320    WM(I+2) = 1.0D0/DI
        WM(I+2) = 1.0D0/DI
 320    CONTINUE
 330  DO 340 I = 1,N
C340    X(I) = WM(I+2)*X(I)
        X(I) = WM(I+2)*X(I)
 340    CONTINUE
      RETURN
 390  IERSL = 1
      RETURN
C
 400  ML = IWM(1)
      MU = IWM(2)
      MEBAND = 2*ML + MU + 1
      CALL DGBSL (WM(3), MEBAND, N, ML, MU, IWM(21), X, 0)
      RETURN
C----------------------- END OF SUBROUTINE DSOLSY ----------------------
      END
*DECK DEWSET
      SUBROUTINE DEWSET (N, ITOL, RTOL, ATOL, YCUR, EWT)
C***BEGIN PROLOGUE  DEWSET
C***SUBSIDIARY
C***PURPOSE  Set error weight vector.
C***TYPE      DOUBLE PRECISION (SEWSET-S, DEWSET-D)
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C
C  This subroutine sets the error weight vector EWT according to
C      EWT(i) = RTOL(i)*ABS(YCUR(i)) + ATOL(i),  i = 1,...,N,
C  with the subscript on RTOL and/or ATOL possibly replaced by 1 above,
C  depending on the value of ITOL.
C
C***SEE ALSO  DLSODE
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791129  DATE WRITTEN
C   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
C   890503  Minor cosmetic changes.  (FNF)
C   930809  Renamed to allow single/double precision versions. (ACH)
C***END PROLOGUE  DEWSET
C**End
      INTEGER N, ITOL
      INTEGER I
      DOUBLE PRECISION RTOL, ATOL, YCUR, EWT
      DIMENSION RTOL(*), ATOL(*), YCUR(N), EWT(N)
C
C***FIRST EXECUTABLE STATEMENT  DEWSET
C     GO TO (10, 20, 30, 40), ITOL
      IF (ITOL .EQ. 1) GO TO 10
      IF (ITOL .EQ. 2) GO TO 20
      IF (ITOL .EQ. 3) GO TO 30
      IF (ITOL .EQ. 4) GO TO 40
 10   CONTINUE
      DO 15 I = 1,N
C15     EWT(I) = RTOL(1)*ABS(YCUR(I)) + ATOL(1)
        EWT(I) = RTOL(1)*ABS(YCUR(I)) + ATOL(1)
 15     CONTINUE
      RETURN
 20   CONTINUE
      DO 25 I = 1,N
C25     EWT(I) = RTOL(1)*ABS(YCUR(I)) + ATOL(I)
        EWT(I) = RTOL(1)*ABS(YCUR(I)) + ATOL(I)
 25     CONTINUE
      RETURN
 30   CONTINUE
      DO 35 I = 1,N
C35     EWT(I) = RTOL(I)*ABS(YCUR(I)) + ATOL(1)
        EWT(I) = RTOL(I)*ABS(YCUR(I)) + ATOL(1)
 35     CONTINUE
      RETURN
 40   CONTINUE
      DO 45 I = 1,N
C45     EWT(I) = RTOL(I)*ABS(YCUR(I)) + ATOL(I)
        EWT(I) = RTOL(I)*ABS(YCUR(I)) + ATOL(I)
 45     CONTINUE
      RETURN
C----------------------- END OF SUBROUTINE DEWSET ----------------------
      END
*DECK DSTODA
      SUBROUTINE DSTODA (NEQ, Y, YH, NYH, YH1, EWT, SAVF, ACOR,
     1   WM, IWM, F, JAC, PJAC, SLVS)
      EXTERNAL F, JAC, PJAC, SLVS
      INTEGER NEQ, NYH, IWM
      DOUBLE PRECISION Y, YH, YH1, EWT, SAVF, ACOR, WM
      DIMENSION NEQ(*), Y(*), YH(NYH,*), YH1(*), EWT(*), SAVF(*),
     1   ACOR(*), WM(*), IWM(*)
      INTEGER IOWND, IALTH, IPUP, LMAX, MEO, NQNYH, NSLP,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     2   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      INTEGER IOWND2, ICOUNT, IRFLAG, JTYP, MUSED, MXORDN, MXORDS
      DOUBLE PRECISION CONIT, CRATE, EL, ELCO, HOLD, RMAX, TESCO,
     2   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      DOUBLE PRECISION ROWND2, CM1, CM2, PDEST, PDLAST, RATIO,
     1   PDNORM
      COMMON /DLS001/ CONIT, CRATE, EL(13), ELCO(13,12),
     1   HOLD, RMAX, TESCO(3,12),
     2   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     3   IOWND(6), IALTH, IPUP, LMAX, MEO, NQNYH, NSLP,
     4   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     5   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     6   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      COMMON /DLSA01/ ROWND2, CM1(12), CM2(5), PDEST, PDLAST, RATIO,
     1   PDNORM,
     2   IOWND2(3), ICOUNT, IRFLAG, JTYP, MUSED, MXORDN, MXORDS
      INTEGER I, I1, IREDO, IRET, J, JB, M, NCF, NEWQ
      INTEGER LM1, LM1P1, LM2, LM2P1, NQM1, NQM2
      DOUBLE PRECISION DCON, DDN, DEL, DELP, DSM, DUP, EXDN, EXSM, EXUP,
     1   R, RH, RHDN, RHSM, RHUP, TOLD, DMNORM
      DOUBLE PRECISION ALPHA, DM1,DM2, EXM1,EXM2,
     1   PDH, PNORM, RATE, RH1, RH1IT, RH2, RM, SM1(12)
      SAVE SM1
      DATA SM1/0.5D0, 0.575D0, 0.55D0, 0.45D0, 0.35D0, 0.25D0,
     1   0.20D0, 0.15D0, 0.10D0, 0.075D0, 0.050D0, 0.025D0/
C-----------------------------------------------------------------------
C DSTODA performs one step of the integration of an initial value
C problem for a system of ordinary differential equations.
C Note: DSTODA is independent of the value of the iteration method
C indicator MITER, when this is .ne. 0, and hence is independent
C of the type of chord method used, or the Jacobian structure.
C Communication with DSTODA is done with the following variables:
C
C Y      = an array of length .ge. N used as the Y argument in
C          all calls to F and JAC.
C NEQ    = integer array containing problem size in NEQ(1), and
C          passed as the NEQ argument in all calls to F and JAC.
C YH     = an NYH by LMAX array containing the dependent variables
C          and their approximate scaled derivatives, where
C          LMAX = MAXORD + 1.  YH(i,j+1) contains the approximate
C          j-th derivative of y(i), scaled by H**j/factorial(j)
C          (j = 0,1,...,NQ).  On entry for the first step, the first
C          two columns of YH must be set from the initial values.
C NYH    = a constant integer .ge. N, the first dimension of YH.
C YH1    = a one-dimensional array occupying the same space as YH.
C EWT    = an array of length N containing multiplicative weights
C          for local error measurements.  Local errors in y(i) are
C          compared to 1.0/EWT(i) in various error tests.
C SAVF   = an array of working storage, of length N.
C ACOR   = a work array of length N, used for the accumulated
C          corrections.  On a successful return, ACOR(i) contains
C          the estimated one-step local error in y(i).
C WM,IWM = real and integer work arrays associated with matrix
C          operations in chord iteration (MITER .ne. 0).
C PJAC   = name of routine to evaluate and preprocess Jacobian matrix
C          and P = I - H*EL0*Jac, if a chord method is being used.
C          It also returns an estimate of norm(Jac) in PDNORM.
C SLVS   = name of routine to solve linear system in chord iteration.
C CCMAX  = maximum relative change in H*EL0 before PJAC is called.
C H      = the step size to be attempted on the next step.
C          H is altered by the error control algorithm during the
C          problem.  H can be either positive or negative, but its
C          sign must remain constant throughout the problem.
C HMIN   = the minimum absolute value of the step size H to be used.
C HMXI   = inverse of the maximum absolute value of H to be used.
C          HMXI = 0.0 is allowed and corresponds to an infinite HMAX.
C          HMIN and HMXI may be changed at any time, but will not
C          take effect until the next change of H is considered.
C TN     = the independent variable. TN is updated on each step taken.
C JSTART = an integer used for input only, with the following
C          values and meanings:
C               0  perform the first step.
C           .gt.0  take a new step continuing from the last.
C              -1  take the next step with a new value of H,
C                    N, METH, MITER, and/or matrix parameters.
C              -2  take the next step with a new value of H,
C                    but with other inputs unchanged.
C          On return, JSTART is set to 1 to facilitate continuation.
C KFLAG  = a completion code with the following meanings:
C               0  the step was succesful.
C              -1  the requested error could not be achieved.
C              -2  corrector convergence could not be achieved.
C              -3  fatal error in PJAC or SLVS.
C          A return with KFLAG = -1 or -2 means either
C          ABS(H) = HMIN or 10 consecutive failures occurred.
C          On a return with KFLAG negative, the values of TN and
C          the YH array are as of the beginning of the last
C          step, and H is the last step size attempted.
C MAXORD = the maximum order of integration method to be allowed.
C MAXCOR = the maximum number of corrector iterations allowed.
C MSBP   = maximum number of steps between PJAC calls (MITER .gt. 0).
C MXNCF  = maximum number of convergence failures allowed.
C METH   = current method.
C          METH = 1 means Adams method (nonstiff)
C          METH = 2 means BDF method (stiff)
C          METH may be reset by DSTODA.
C MITER  = corrector iteration method.
C          MITER = 0 means functional iteration.
C          MITER = JT .gt. 0 means a chord iteration corresponding
C          to Jacobian type JT.  (The DLSODA/DLSODAR argument JT is
C          communicated here as JTYP, but is not used in DSTODA
C          except to load MITER following a method switch.)
C          MITER may be reset by DSTODA.
C N      = the number of first-order differential equations.
C-----------------------------------------------------------------------
      KFLAG = 0
      TOLD = TN
      NCF = 0
      IERPJ = 0
      IERSL = 0
      JCUR = 0
      ICF = 0
      DELP = 0.0D0
      IF (JSTART .GT. 0) GO TO 200
      IF (JSTART .EQ. -1) GO TO 100
      IF (JSTART .EQ. -2) GO TO 160
C-----------------------------------------------------------------------
C On the first call, the order is set to 1, and other variables are
C initialized.  RMAX is the maximum ratio by which H can be increased
C in a single step.  It is initially 1.E4 to compensate for the small
C initial H, but then is normally equal to 10.  If a failure
C occurs (in corrector convergence or error test), RMAX is set at 2
C for the next increase.
C DCFODE is called to get the needed coefficients for both methods.
C-----------------------------------------------------------------------
      LMAX = MAXORD + 1
      NQ = 1
      L = 2
      IALTH = 2
      RMAX = 10000.0D0
      RC = 0.0D0
      EL0 = 1.0D0
      CRATE = 0.7D0
      HOLD = H
      NSLP = 0
      IPUP = MITER
      IRET = 3
C Initialize switching parameters.  METH = 1 is assumed initially. -----
      ICOUNT = 20
      IRFLAG = 0
      PDEST = 0.0D0
      PDLAST = 0.0D0
      RATIO = 5.0D0
      CALL DCFODE (2, ELCO, TESCO)
      DO 10 I = 1,5
C10     CM2(I) = TESCO(2,I)*ELCO(I+1,I)
        CM2(I) = TESCO(2,I)*ELCO(I+1,I)
 10     CONTINUE
      CALL DCFODE (1, ELCO, TESCO)
      DO 20 I = 1,12
C20     CM1(I) = TESCO(2,I)*ELCO(I+1,I)
        CM1(I) = TESCO(2,I)*ELCO(I+1,I)
 20     CONTINUE
      GO TO 150
C-----------------------------------------------------------------------
C The following block handles preliminaries needed when JSTART = -1.
C IPUP is set to MITER to force a matrix update.
C If an order increase is about to be considered (IALTH = 1),
C IALTH is reset to 2 to postpone consideration one more step.
C If the caller has changed METH, DCFODE is called to reset
C the coefficients of the method.
C If H is to be changed, YH must be rescaled.
C If H or METH is being changed, IALTH is reset to L = NQ + 1
C to prevent further changes in H for that many steps.
C-----------------------------------------------------------------------
 100  IPUP = MITER
      LMAX = MAXORD + 1
      IF (IALTH .EQ. 1) IALTH = 2
      IF (METH .EQ. MUSED) GO TO 160
      CALL DCFODE (METH, ELCO, TESCO)
      IALTH = L
      IRET = 1
C-----------------------------------------------------------------------
C The el vector and related constants are reset
C whenever the order NQ is changed, or at the start of the problem.
C-----------------------------------------------------------------------
 150  DO 155 I = 1,L
C155    EL(I) = ELCO(I,NQ)
        EL(I) = ELCO(I,NQ)
 155    CONTINUE
      NQNYH = NQ*NYH
      RC = RC*EL(1)/EL0
      EL0 = EL(1)
      CONIT = 0.5D0/(NQ+2)
C     GO TO (160, 170, 200), IRET
      IF (IRET .EQ. 1) GO TO 160
      IF (IRET .EQ. 2) GO TO 170
      IF (IRET .EQ. 3) GO TO 200
C-----------------------------------------------------------------------
C If H is being changed, the H ratio RH is checked against
C RMAX, HMIN, and HMXI, and the YH array rescaled.  IALTH is set to
C L = NQ + 1 to prevent a change of H for that many steps, unless
C forced by a convergence or error test failure.
C-----------------------------------------------------------------------
 160  IF (H .EQ. HOLD) GO TO 200
      RH = H/HOLD
      H = HOLD
      IREDO = 3
      GO TO 175
 170  RH = MAX(RH,HMIN/ABS(H))
 175  RH = MIN(RH,RMAX)
      RH = RH/MAX(1.0D0,ABS(H)*HMXI*RH)
C-----------------------------------------------------------------------
C If METH = 1, also restrict the new step size by the stability region.
C If this reduces H, set IRFLAG to 1 so that if there are roundoff
C problems later, we can assume that is the cause of the trouble.
C-----------------------------------------------------------------------
      IF (METH .EQ. 2) GO TO 178
      IRFLAG = 0
      PDH = MAX(ABS(H)*PDLAST,0.000001D0)
      IF (RH*PDH*1.00001D0 .LT. SM1(NQ)) GO TO 178
      RH = SM1(NQ)/PDH
      IRFLAG = 1
 178  CONTINUE
      R = 1.0D0
      DO 181 J = 2,L
        R = R*RH
        DO 180 I = 1,N
C180      YH(I,J) = YH(I,J)*R
          YH(I,J) = YH(I,J)*R
 180      CONTINUE
 181      CONTINUE
      H = H*RH
      RC = RC*RH
      IALTH = L
      IF (IREDO .EQ. 0) GO TO 690
C-----------------------------------------------------------------------
C This section computes the predicted values by effectively
C multiplying the YH array by the Pascal triangle matrix.
C RC is the ratio of new to old values of the coefficient  H*EL(1).
C When RC differs from 1 by more than CCMAX, IPUP is set to MITER
C to force PJAC to be called, if a Jacobian is involved.
C In any case, PJAC is called at least every MSBP steps.
C-----------------------------------------------------------------------
 200  IF (ABS(RC-1.0D0) .GT. CCMAX) IPUP = MITER
      IF (NST .GE. NSLP+MSBP) IPUP = MITER
      TN = TN + H
      I1 = NQNYH + 1
      DO 215 JB = 1,NQ
        I1 = I1 - NYH
CDIR$ IVDEP
        DO 210 I = I1,NQNYH
C210      YH1(I) = YH1(I) + YH1(I+NYH)
          YH1(I) = YH1(I) + YH1(I+NYH)
 210      CONTINUE
 215    CONTINUE
      PNORM = DMNORM (N, YH1, EWT)
C-----------------------------------------------------------------------
C Up to MAXCOR corrector iterations are taken.  A convergence test is
C made on the RMS-norm of each correction, weighted by the error
C weight vector EWT.  The sum of the corrections is accumulated in the
C vector ACOR(i).  The YH array is not altered in the corrector loop.
C-----------------------------------------------------------------------
 220  M = 0
      RATE = 0.0D0
      DEL = 0.0D0
      DO 230 I = 1,N
C230    Y(I) = YH(I,1)
        Y(I) = YH(I,1)
 230    CONTINUE
      CALL F (NEQ, TN, Y, SAVF)
      NFE = NFE + 1
      IF (IPUP .LE. 0) GO TO 250
C-----------------------------------------------------------------------
C If indicated, the matrix P = I - H*EL(1)*J is reevaluated and
C preprocessed before starting the corrector iteration.  IPUP is set
C to 0 as an indicator that this has been done.
C-----------------------------------------------------------------------
      CALL PJAC (NEQ, Y, YH, NYH, EWT, ACOR, SAVF, WM, IWM, F, JAC)
      IPUP = 0
      RC = 1.0D0
      NSLP = NST
      CRATE = 0.7D0
      IF (IERPJ .NE. 0) GO TO 430
 250  DO 260 I = 1,N
C260    ACOR(I) = 0.0D0
        ACOR(I) = 0.0D0
 260    CONTINUE
 270  IF (MITER .NE. 0) GO TO 350
C-----------------------------------------------------------------------
C In the case of functional iteration, update Y directly from
C the result of the last function evaluation.
C-----------------------------------------------------------------------
      DO 290 I = 1,N
        SAVF(I) = H*SAVF(I) - YH(I,2)
C290    Y(I) = SAVF(I) - ACOR(I)
        Y(I) = SAVF(I) - ACOR(I)
 290    CONTINUE
      DEL = DMNORM (N, Y, EWT)
      DO 300 I = 1,N
        Y(I) = YH(I,1) + EL(1)*SAVF(I)
C300    ACOR(I) = SAVF(I)
        ACOR(I) = SAVF(I)
 300    CONTINUE
      GO TO 400
C-----------------------------------------------------------------------
C In the case of the chord method, compute the corrector error,
C and solve the linear system with that as right-hand side and
C P as coefficient matrix.
C-----------------------------------------------------------------------
 350  DO 360 I = 1,N
C360    Y(I) = H*SAVF(I) - (YH(I,2) + ACOR(I))
        Y(I) = H*SAVF(I) - (YH(I,2) + ACOR(I))
 360    CONTINUE
      CALL SLVS (WM, IWM, Y, SAVF)
      IF (IERSL .LT. 0) GO TO 430
      IF (IERSL .GT. 0) GO TO 410
      DEL = DMNORM (N, Y, EWT)
      DO 380 I = 1,N
        ACOR(I) = ACOR(I) + Y(I)
C380    Y(I) = YH(I,1) + EL(1)*ACOR(I)
        Y(I) = YH(I,1) + EL(1)*ACOR(I)
 380    CONTINUE
C-----------------------------------------------------------------------
C Test for convergence.  If M .gt. 0, an estimate of the convergence
C rate constant is stored in CRATE, and this is used in the test.
C
C We first check for a change of iterates that is the size of
C roundoff error.  If this occurs, the iteration has converged, and a
C new rate estimate is not formed.
C In all other cases, force at least two iterations to estimate a
C local Lipschitz constant estimate for Adams methods.
C On convergence, form PDEST = local maximum Lipschitz constant
C estimate.  PDLAST is the most recent nonzero estimate.
C-----------------------------------------------------------------------
 400  CONTINUE
      IF (DEL .LE. 100.0D0*PNORM*UROUND) GO TO 450
      IF (M .EQ. 0 .AND. METH .EQ. 1) GO TO 405
      IF (M .EQ. 0) GO TO 402
      RM = 1024.0D0
      IF (DEL .LE. 1024.0D0*DELP) RM = DEL/DELP
      RATE = MAX(RATE,RM)
      CRATE = MAX(0.2D0*CRATE,RM)
 402  DCON = DEL*MIN(1.0D0,1.5D0*CRATE)/(TESCO(2,NQ)*CONIT)
      IF (DCON .GT. 1.0D0) GO TO 405
      PDEST = MAX(PDEST,RATE/ABS(H*EL(1)))
      IF (PDEST .NE. 0.0D0) PDLAST = PDEST
      GO TO 450
 405  CONTINUE
      M = M + 1
      IF (M .EQ. MAXCOR) GO TO 410
      IF (M .GE. 2 .AND. DEL .GT. 2.0D0*DELP) GO TO 410
      DELP = DEL
      CALL F (NEQ, TN, Y, SAVF)
      NFE = NFE + 1
      GO TO 270
C-----------------------------------------------------------------------
C The corrector iteration failed to converge.
C If MITER .ne. 0 and the Jacobian is out of date, PJAC is called for
C the next try.  Otherwise the YH array is retracted to its values
C before prediction, and H is reduced, if possible.  If H cannot be
C reduced or MXNCF failures have occurred, exit with KFLAG = -2.
C-----------------------------------------------------------------------
 410  IF (MITER .EQ. 0 .OR. JCUR .EQ. 1) GO TO 430
      ICF = 1
      IPUP = MITER
      GO TO 220
 430  ICF = 2
      NCF = NCF + 1
      RMAX = 2.0D0
      TN = TOLD
      I1 = NQNYH + 1
      DO 445 JB = 1,NQ
        I1 = I1 - NYH
CDIR$ IVDEP
        DO 440 I = I1,NQNYH
C440      YH1(I) = YH1(I) - YH1(I+NYH)
          YH1(I) = YH1(I) - YH1(I+NYH)
 440      CONTINUE
 445    CONTINUE
      IF (IERPJ .LT. 0 .OR. IERSL .LT. 0) GO TO 680
      IF (ABS(H) .LE. HMIN*1.00001D0) GO TO 670
      IF (NCF .EQ. MXNCF) GO TO 670
      RH = 0.25D0
      IPUP = MITER
      IREDO = 1
      GO TO 170
C-----------------------------------------------------------------------
C The corrector has converged.  JCUR is set to 0
C to signal that the Jacobian involved may need updating later.
C The local error test is made and control passes to statement 500
C if it fails.
C-----------------------------------------------------------------------
 450  JCUR = 0
      IF (M .EQ. 0) DSM = DEL/TESCO(2,NQ)
      IF (M .GT. 0) DSM = DMNORM (N, ACOR, EWT)/TESCO(2,NQ)
      IF (DSM .GT. 1.0D0) GO TO 500
C-----------------------------------------------------------------------
C After a successful step, update the YH array.
C Decrease ICOUNT by 1, and if it is -1, consider switching methods.
C If a method switch is made, reset various parameters,
C rescale the YH array, and exit.  If there is no switch,
C consider changing H if IALTH = 1.  Otherwise decrease IALTH by 1.
C If IALTH is then 1 and NQ .lt. MAXORD, then ACOR is saved for
C use in a possible order increase on the next step.
C If a change in H is considered, an increase or decrease in order
C by one is considered also.  A change in H is made only if it is by a
C factor of at least 1.1.  If not, IALTH is set to 3 to prevent
C testing for that many steps.
C-----------------------------------------------------------------------
      KFLAG = 0
      IREDO = 0
      NST = NST + 1
      HU = H
      NQU = NQ
      MUSED = METH
      DO 461 J = 1,L
        DO 460 I = 1,N
C460      YH(I,J) = YH(I,J) + EL(J)*ACOR(I)
          YH(I,J) = YH(I,J) + EL(J)*ACOR(I)
 460      CONTINUE
 461      CONTINUE
      ICOUNT = ICOUNT - 1
      IF (ICOUNT .GE. 0) GO TO 488
      IF (METH .EQ. 2) GO TO 480
C-----------------------------------------------------------------------
C We are currently using an Adams method.  Consider switching to BDF.
C If the current order is greater than 5, assume the problem is
C not stiff, and skip this section.
C If the Lipschitz constant and error estimate are not polluted
C by roundoff, go to 470 and perform the usual test.
C Otherwise, switch to the BDF methods if the last step was
C restricted to insure stability (irflag = 1), and stay with Adams
C method if not.  When switching to BDF with polluted error estimates,
C in the absence of other information, double the step size.
C
C When the estimates are OK, we make the usual test by computing
C the step size we could have (ideally) used on this step,
C with the current (Adams) method, and also that for the BDF.
C If NQ .gt. MXORDS, we consider changing to order MXORDS on switching.
C Compare the two step sizes to decide whether to switch.
C The step size advantage must be at least RATIO = 5 to switch.
C-----------------------------------------------------------------------
      IF (NQ .GT. 5) GO TO 488
      IF (DSM .GT. 100.0D0*PNORM*UROUND .AND. PDEST .NE. 0.0D0)
     1   GO TO 470
      IF (IRFLAG .EQ. 0) GO TO 488
      RH2 = 2.0D0
      NQM2 = MIN(NQ,MXORDS)
      GO TO 478
 470  CONTINUE
      EXSM = 1.0D0/L
      RH1 = 1.0D0/(1.2D0*DSM**EXSM + 0.0000012D0)
      RH1IT = 2.0D0*RH1
      PDH = PDLAST*ABS(H)
      IF (PDH*RH1 .GT. 0.00001D0) RH1IT = SM1(NQ)/PDH
      RH1 = MIN(RH1,RH1IT)
      IF (NQ .LE. MXORDS) GO TO 474
         NQM2 = MXORDS
         LM2 = MXORDS + 1
         EXM2 = 1.0D0/LM2
         LM2P1 = LM2 + 1
         DM2 = DMNORM (N, YH(1,LM2P1), EWT)/CM2(MXORDS)
         RH2 = 1.0D0/(1.2D0*DM2**EXM2 + 0.0000012D0)
         GO TO 476
 474  DM2 = DSM*(CM1(NQ)/CM2(NQ))
      RH2 = 1.0D0/(1.2D0*DM2**EXSM + 0.0000012D0)
      NQM2 = NQ
 476  CONTINUE
      IF (RH2 .LT. RATIO*RH1) GO TO 488
C THE SWITCH TEST PASSED.  RESET RELEVANT QUANTITIES FOR BDF. ----------
 478  RH = RH2
      ICOUNT = 20
      METH = 2
      MITER = JTYP
      PDLAST = 0.0D0
      NQ = NQM2
      L = NQ + 1
      GO TO 170
C-----------------------------------------------------------------------
C We are currently using a BDF method.  Consider switching to Adams.
C Compute the step size we could have (ideally) used on this step,
C with the current (BDF) method, and also that for the Adams.
C If NQ .gt. MXORDN, we consider changing to order MXORDN on switching.
C Compare the two step sizes to decide whether to switch.
C The step size advantage must be at least 5/RATIO = 1 to switch.
C If the step size for Adams would be so small as to cause
C roundoff pollution, we stay with BDF.
C-----------------------------------------------------------------------
 480  CONTINUE
      EXSM = 1.0D0/L
      IF (MXORDN .GE. NQ) GO TO 484
         NQM1 = MXORDN
         LM1 = MXORDN + 1
         EXM1 = 1.0D0/LM1
         LM1P1 = LM1 + 1
         DM1 = DMNORM (N, YH(1,LM1P1), EWT)/CM1(MXORDN)
         RH1 = 1.0D0/(1.2D0*DM1**EXM1 + 0.0000012D0)
         GO TO 486
 484  DM1 = DSM*(CM2(NQ)/CM1(NQ))
      RH1 = 1.0D0/(1.2D0*DM1**EXSM + 0.0000012D0)
      NQM1 = NQ
      EXM1 = EXSM
 486  RH1IT = 2.0D0*RH1
      PDH = PDNORM*ABS(H)
      IF (PDH*RH1 .GT. 0.00001D0) RH1IT = SM1(NQM1)/PDH
      RH1 = MIN(RH1,RH1IT)
      RH2 = 1.0D0/(1.2D0*DSM**EXSM + 0.0000012D0)
      IF (RH1*RATIO .LT. 5.0D0*RH2) GO TO 488
      ALPHA = MAX(0.001D0,RH1)
      DM1 = (ALPHA**EXM1)*DM1
      IF (DM1 .LE. 1000.0D0*UROUND*PNORM) GO TO 488
C The switch test passed.  Reset relevant quantities for Adams. --------
      RH = RH1
      ICOUNT = 20
      METH = 1
      MITER = 0
      PDLAST = 0.0D0
      NQ = NQM1
      L = NQ + 1
      GO TO 170
C
C No method switch is being made.  Do the usual step/order selection. --
 488  CONTINUE
      IALTH = IALTH - 1
      IF (IALTH .EQ. 0) GO TO 520
      IF (IALTH .GT. 1) GO TO 700
      IF (L .EQ. LMAX) GO TO 700
      DO 490 I = 1,N
C490    YH(I,LMAX) = ACOR(I)
        YH(I,LMAX) = ACOR(I)
 490    CONTINUE
      GO TO 700
C-----------------------------------------------------------------------
C The error test failed.  KFLAG keeps track of multiple failures.
C Restore TN and the YH array to their previous values, and prepare
C to try the step again.  Compute the optimum step size for this or
C one lower order.  After 2 or more failures, H is forced to decrease
C by a factor of 0.2 or less.
C-----------------------------------------------------------------------
 500  KFLAG = KFLAG - 1
      TN = TOLD
      I1 = NQNYH + 1
      DO 515 JB = 1,NQ
        I1 = I1 - NYH
CDIR$ IVDEP
        DO 510 I = I1,NQNYH
C510      YH1(I) = YH1(I) - YH1(I+NYH)
          YH1(I) = YH1(I) - YH1(I+NYH)
 510      CONTINUE
 515    CONTINUE
      RMAX = 2.0D0
      IF (ABS(H) .LE. HMIN*1.00001D0) GO TO 660
      IF (KFLAG .LE. -3) GO TO 640
      IREDO = 2
      RHUP = 0.0D0
      GO TO 540
C-----------------------------------------------------------------------
C Regardless of the success or failure of the step, factors
C RHDN, RHSM, and RHUP are computed, by which H could be multiplied
C at order NQ - 1, order NQ, or order NQ + 1, respectively.
C In the case of failure, RHUP = 0.0 to avoid an order increase.
C The largest of these is determined and the new order chosen
C accordingly.  If the order is to be increased, we compute one
C additional scaled derivative.
C-----------------------------------------------------------------------
 520  RHUP = 0.0D0
      IF (L .EQ. LMAX) GO TO 540
      DO 530 I = 1,N
C530    SAVF(I) = ACOR(I) - YH(I,LMAX)
        SAVF(I) = ACOR(I) - YH(I,LMAX)
 530    CONTINUE
      DUP = DMNORM (N, SAVF, EWT)/TESCO(3,NQ)
      EXUP = 1.0D0/(L+1)
      RHUP = 1.0D0/(1.4D0*DUP**EXUP + 0.0000014D0)
 540  EXSM = 1.0D0/L
      RHSM = 1.0D0/(1.2D0*DSM**EXSM + 0.0000012D0)
      RHDN = 0.0D0
      IF (NQ .EQ. 1) GO TO 550
      DDN = DMNORM (N, YH(1,L), EWT)/TESCO(1,NQ)
      EXDN = 1.0D0/NQ
      RHDN = 1.0D0/(1.3D0*DDN**EXDN + 0.0000013D0)
C If METH = 1, limit RH according to the stability region also. --------
 550  IF (METH .EQ. 2) GO TO 560
      PDH = MAX(ABS(H)*PDLAST,0.000001D0)
      IF (L .LT. LMAX) RHUP = MIN(RHUP,SM1(L)/PDH)
      RHSM = MIN(RHSM,SM1(NQ)/PDH)
      IF (NQ .GT. 1) RHDN = MIN(RHDN,SM1(NQ-1)/PDH)
      PDEST = 0.0D0
 560  IF (RHSM .GE. RHUP) GO TO 570
      IF (RHUP .GT. RHDN) GO TO 590
      GO TO 580
 570  IF (RHSM .LT. RHDN) GO TO 580
      NEWQ = NQ
      RH = RHSM
      GO TO 620
 580  NEWQ = NQ - 1
      RH = RHDN
      IF (KFLAG .LT. 0 .AND. RH .GT. 1.0D0) RH = 1.0D0
      GO TO 620
 590  NEWQ = L
      RH = RHUP
      IF (RH .LT. 1.1D0) GO TO 610
      R = EL(L)/L
      DO 600 I = 1,N
C600    YH(I,NEWQ+1) = ACOR(I)*R
        YH(I,NEWQ+1) = ACOR(I)*R
 600    CONTINUE
      GO TO 630
 610  IALTH = 3
      GO TO 700
C If METH = 1 and H is restricted by stability, bypass 10 percent test.
 620  IF (METH .EQ. 2) GO TO 622
      IF (RH*PDH*1.00001D0 .GE. SM1(NEWQ)) GO TO 625
 622  IF (KFLAG .EQ. 0 .AND. RH .LT. 1.1D0) GO TO 610
 625  IF (KFLAG .LE. -2) RH = MIN(RH,0.2D0)
C-----------------------------------------------------------------------
C If there is a change of order, reset NQ, L, and the coefficients.
C In any case H is reset according to RH and the YH array is rescaled.
C Then exit from 690 if the step was OK, or redo the step otherwise.
C-----------------------------------------------------------------------
      IF (NEWQ .EQ. NQ) GO TO 170
 630  NQ = NEWQ
      L = NQ + 1
      IRET = 2
      GO TO 150
C-----------------------------------------------------------------------
C Control reaches this section if 3 or more failures have occured.
C If 10 failures have occurred, exit with KFLAG = -1.
C It is assumed that the derivatives that have accumulated in the
C YH array have errors of the wrong order.  Hence the first
C derivative is recomputed, and the order is set to 1.  Then
C H is reduced by a factor of 10, and the step is retried,
C until it succeeds or H reaches HMIN.
C-----------------------------------------------------------------------
 640  IF (KFLAG .EQ. -10) GO TO 660
      RH = 0.1D0
      RH = MAX(HMIN/ABS(H),RH)
      H = H*RH
      DO 645 I = 1,N
C645    Y(I) = YH(I,1)
        Y(I) = YH(I,1)
 645    CONTINUE
      CALL F (NEQ, TN, Y, SAVF)
      NFE = NFE + 1
      DO 650 I = 1,N
C650    YH(I,2) = H*SAVF(I)
        YH(I,2) = H*SAVF(I)
 650    CONTINUE
      IPUP = MITER
      IALTH = 5
      IF (NQ .EQ. 1) GO TO 200
      NQ = 1
      L = 2
      IRET = 3
      GO TO 150
C-----------------------------------------------------------------------
C All returns are made through this section.  H is saved in HOLD
C to allow the caller to change H on the next step.
C-----------------------------------------------------------------------
 660  KFLAG = -1
      GO TO 720
 670  KFLAG = -2
      GO TO 720
 680  KFLAG = -3
      GO TO 720
 690  RMAX = 10.0D0
 700  R = 1.0D0/TESCO(2,NQU)
      DO 710 I = 1,N
C710    ACOR(I) = ACOR(I)*R
        ACOR(I) = ACOR(I)*R
 710    CONTINUE
 720  HOLD = H
      JSTART = 1
      RETURN
C----------------------- End of Subroutine DSTODA ----------------------
      END
*DECK DPRJA
      SUBROUTINE DPRJA (NEQ, Y, YH, NYH, EWT, FTEM, SAVF, WM, IWM,
     1   F, JAC)
      EXTERNAL F, JAC
      INTEGER NEQ, NYH, IWM
      DOUBLE PRECISION Y, YH, EWT, FTEM, SAVF, WM
      DIMENSION NEQ(*), Y(*), YH(NYH,*), EWT(*), FTEM(*), SAVF(*),
     1   WM(*), IWM(*)
      INTEGER IOWND, IOWNS,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     2   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      INTEGER IOWND2, IOWNS2, JTYP, MUSED, MXORDN, MXORDS
      DOUBLE PRECISION ROWNS,
     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      DOUBLE PRECISION ROWND2, ROWNS2, PDNORM
      COMMON /DLS001/ ROWNS(209),
     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     2   IOWND(6), IOWNS(6),
     3   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     4   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     5   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      COMMON /DLSA01/ ROWND2, ROWNS2(20), PDNORM,
     1   IOWND2(3), IOWNS2(2), JTYP, MUSED, MXORDN, MXORDS
      INTEGER I, I1, I2, IER, II, J, J1, JJ, LENP,
     1   MBA, MBAND, MEB1, MEBAND, ML, ML3, MU, NP1
      DOUBLE PRECISION CON, FAC, HL0, R, R0, SRUR, YI, YJ, YJJ,
     1   DMNORM, DFNORM, DBNORM
C-----------------------------------------------------------------------
C DPRJA is called by DSTODA to compute and process the matrix
C P = I - H*EL(1)*J , where J is an approximation to the Jacobian.
C Here J is computed by the user-supplied routine JAC if
C MITER = 1 or 4 or by finite differencing if MITER = 2 or 5.
C J, scaled by -H*EL(1), is stored in WM.  Then the norm of J (the
C matrix norm consistent with the weighted max-norm on vectors given
C by DMNORM) is computed, and J is overwritten by P.  P is then
C subjected to LU decomposition in preparation for later solution
C of linear systems with P as coefficient matrix.  This is done
C by DGEFA if MITER = 1 or 2, and by DGBFA if MITER = 4 or 5.
C
C In addition to variables described previously, communication
C with DPRJA uses the following:
C Y     = array containing predicted values on entry.
C FTEM  = work array of length N (ACOR in DSTODA).
C SAVF  = array containing f evaluated at predicted y.
C WM    = real work space for matrices.  On output it contains the
C         LU decomposition of P.
C         Storage of matrix elements starts at WM(3).
C         WM also contains the following matrix-related data:
C         WM(1) = SQRT(UROUND), used in numerical Jacobian increments.
C IWM   = integer work space containing pivot information, starting at
C         IWM(21).   IWM also contains the band parameters
C         ML = IWM(1) and MU = IWM(2) if MITER is 4 or 5.
C EL0   = EL(1) (input).
C PDNORM= norm of Jacobian matrix. (Output).
C IERPJ = output error flag,  = 0 if no trouble, .gt. 0 if
C         P matrix found to be singular.
C JCUR  = output flag = 1 to indicate that the Jacobian matrix
C         (or approximation) is now current.
C This routine also uses the Common variables EL0, H, TN, UROUND,
C MITER, N, NFE, and NJE.
C-----------------------------------------------------------------------
      NJE = NJE + 1
      IERPJ = 0
      JCUR = 1
      HL0 = H*EL0
C     GO TO (100, 200, 300, 400, 500), MITER
      IF (MITER .EQ. 1) GO TO 100
      IF (MITER .EQ. 2) GO TO 200
      IF (MITER .EQ. 3) GO TO 300
      IF (MITER .EQ. 4) GO TO 400
      IF (MITER .EQ. 5) GO TO 500
C If MITER = 1, call JAC and multiply by scalar. -----------------------
 100  LENP = N*N
      DO 110 I = 1,LENP
C110    WM(I+2) = 0.0D0
        WM(I+2) = 0.0D0
 110    CONTINUE
      CALL JAC (NEQ, TN, Y, 0, 0, WM(3), N)
      CON = -HL0
      DO 120 I = 1,LENP
C120    WM(I+2) = WM(I+2)*CON
        WM(I+2) = WM(I+2)*CON
 120    CONTINUE
      GO TO 240
C If MITER = 2, make N calls to F to approximate J. --------------------
 200  FAC = DMNORM (N, SAVF, EWT)
      R0 = 1000.0D0*ABS(H)*UROUND*N*FAC
      IF (R0 .EQ. 0.0D0) R0 = 1.0D0
      SRUR = WM(1)
      J1 = 2
      DO 230 J = 1,N
        YJ = Y(J)
        R = MAX(SRUR*ABS(YJ),R0/EWT(J))
        Y(J) = Y(J) + R
        FAC = -HL0/R
        CALL F (NEQ, TN, Y, FTEM)
        DO 220 I = 1,N
C220      WM(I+J1) = (FTEM(I) - SAVF(I))*FAC
          WM(I+J1) = (FTEM(I) - SAVF(I))*FAC
 220      CONTINUE
        Y(J) = YJ
        J1 = J1 + N
 230    CONTINUE
      NFE = NFE + N
 240  CONTINUE
C Compute norm of Jacobian. --------------------------------------------
      PDNORM = DFNORM (N, WM(3), EWT)/ABS(HL0)
C Add identity matrix. -------------------------------------------------
      J = 3
      NP1 = N + 1
      DO 250 I = 1,N
        WM(J) = WM(J) + 1.0D0
C250    J = J + NP1
        J = J + NP1
 250    CONTINUE
C Do LU decomposition on P. --------------------------------------------
      CALL DGEFA (WM(3), N, N, IWM(21), IER)
      IF (IER .NE. 0) IERPJ = 1
      RETURN
C Dummy block only, since MITER is never 3 in this routine. ------------
 300  RETURN
C If MITER = 4, call JAC and multiply by scalar. -----------------------
 400  ML = IWM(1)
      MU = IWM(2)
      ML3 = ML + 3
      MBAND = ML + MU + 1
      MEBAND = MBAND + ML
      LENP = MEBAND*N
      DO 410 I = 1,LENP
C410    WM(I+2) = 0.0D0
        WM(I+2) = 0.0D0
 410    CONTINUE
      CALL JAC (NEQ, TN, Y, ML, MU, WM(ML3), MEBAND)
      CON = -HL0
      DO 420 I = 1,LENP
C420    WM(I+2) = WM(I+2)*CON
        WM(I+2) = WM(I+2)*CON
 420    CONTINUE
      GO TO 570
C If MITER = 5, make MBAND calls to F to approximate J. ----------------
 500  ML = IWM(1)
      MU = IWM(2)
      MBAND = ML + MU + 1
      MBA = MIN(MBAND,N)
      MEBAND = MBAND + ML
      MEB1 = MEBAND - 1
      SRUR = WM(1)
      FAC = DMNORM (N, SAVF, EWT)
      R0 = 1000.0D0*ABS(H)*UROUND*N*FAC
      IF (R0 .EQ. 0.0D0) R0 = 1.0D0
      DO 560 J = 1,MBA
        DO 530 I = J,N,MBAND
          YI = Y(I)
          R = MAX(SRUR*ABS(YI),R0/EWT(I))
C530      Y(I) = Y(I) + R
          Y(I) = Y(I) + R
 530      CONTINUE
        CALL F (NEQ, TN, Y, FTEM)
        DO 550 JJ = J,N,MBAND
          Y(JJ) = YH(JJ,1)
          YJJ = Y(JJ)
          R = MAX(SRUR*ABS(YJJ),R0/EWT(JJ))
          FAC = -HL0/R
          I1 = MAX(JJ-MU,1)
          I2 = MIN(JJ+ML,N)
          II = JJ*MEB1 - ML + 2
          DO 540 I = I1,I2
C540        WM(II+I) = (FTEM(I) - SAVF(I))*FAC
            WM(II+I) = (FTEM(I) - SAVF(I))*FAC
 540        CONTINUE
 550      CONTINUE
 560    CONTINUE
      NFE = NFE + MBA
 570  CONTINUE
C Compute norm of Jacobian. --------------------------------------------
      PDNORM = DBNORM (N, WM(ML+3), MEBAND, ML, MU, EWT)/ABS(HL0)
C Add identity matrix. -------------------------------------------------
      II = MBAND + 2
      DO 580 I = 1,N
        WM(II) = WM(II) + 1.0D0
C580    II = II + MEBAND
        II = II + MEBAND
 580    CONTINUE
C Do LU decomposition of P. --------------------------------------------
      CALL DGBFA (WM(3), MEBAND, N, ML, MU, IWM(21), IER)
      IF (IER .NE. 0) IERPJ = 1
      RETURN
C----------------------- End of Subroutine DPRJA -----------------------
      END
*DECK DMNORM
      DOUBLE PRECISION FUNCTION DMNORM (N, V, W)
C-----------------------------------------------------------------------
C This function routine computes the weighted max-norm
C of the vector of length N contained in the array V, with weights
C contained in the array w of length N:
C   DMNORM = MAX(i=1,...,N) ABS(V(i))*W(i)
C-----------------------------------------------------------------------
      INTEGER N,   I
      DOUBLE PRECISION V, W,   VM
      DIMENSION V(N), W(N)
      VM = 0.0D0
      DO 10 I = 1,N
C10     VM = MAX(VM,ABS(V(I))*W(I))
        VM = MAX(VM,ABS(V(I))*W(I))
 10     CONTINUE
      DMNORM = VM
      RETURN
C----------------------- End of Function DMNORM ------------------------
      END
*DECK DFNORM
      DOUBLE PRECISION FUNCTION DFNORM (N, A, W)
C-----------------------------------------------------------------------
C This function computes the norm of a full N by N matrix,
C stored in the array A, that is consistent with the weighted max-norm
C on vectors, with weights stored in the array W:
C   DFNORM = MAX(i=1,...,N) ( W(i) * Sum(j=1,...,N) ABS(a(i,j))/W(j) )
C-----------------------------------------------------------------------
      INTEGER N,   I, J
      DOUBLE PRECISION A,   W, AN, SUM
      DIMENSION A(N,N), W(N)
      AN = 0.0D0
      DO 20 I = 1,N
        SUM = 0.0D0
        DO 10 J = 1,N
C10       SUM = SUM + ABS(A(I,J))/W(J)
          SUM = SUM + ABS(A(I,J))/W(J)
 10       CONTINUE
        AN = MAX(AN,SUM*W(I))
 20     CONTINUE
      DFNORM = AN
      RETURN
C----------------------- End of Function DFNORM ------------------------
      END
*DECK DBNORM
      DOUBLE PRECISION FUNCTION DBNORM (N, A, NRA, ML, MU, W)
C-----------------------------------------------------------------------
C This function computes the norm of a banded N by N matrix,
C stored in the array A, that is consistent with the weighted max-norm
C on vectors, with weights stored in the array W.
C ML and MU are the lower and upper half-bandwidths of the matrix.
C NRA is the first dimension of the A array, NRA .ge. ML+MU+1.
C In terms of the matrix elements a(i,j), the norm is given by:
C   DBNORM = MAX(i=1,...,N) ( W(i) * Sum(j=1,...,N) ABS(a(i,j))/W(j) )
C-----------------------------------------------------------------------
      INTEGER N, NRA, ML, MU
      INTEGER I, I1, JLO, JHI, J
      DOUBLE PRECISION A, W
      DOUBLE PRECISION AN, SUM
      DIMENSION A(NRA,N), W(N)
      AN = 0.0D0
      DO 20 I = 1,N
        SUM = 0.0D0
        I1 = I + MU + 1
        JLO = MAX(I-ML,1)
        JHI = MIN(I+MU,N)
        DO 10 J = JLO,JHI
C10       SUM = SUM + ABS(A(I1-J,J))/W(J)
          SUM = SUM + ABS(A(I1-J,J))/W(J)
 10       CONTINUE
        AN = MAX(AN,SUM*W(I))
 20     CONTINUE
      DBNORM = AN
      RETURN
C----------------------- End of Function DBNORM ------------------------
      END

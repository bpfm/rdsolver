*     DGEES Example Program Text
*     NAG Copyright 2005.
*     .. Parameters ..
      INTEGER          NIN, NOUT
      PARAMETER        (NIN=5,NOUT=6)
      INTEGER          NB, NMAX
      PARAMETER        (NB=64,NMAX=10)
      INTEGER          LDA, LDVS, LWORK
      PARAMETER        (LDA=NMAX,LDVS=NMAX,LWORK=(2+NB)*NMAX)
*     .. Local Scalars ..
      INTEGER          I, IFAIL, INFO, J, LWKOPT, N, SDIM
*     .. Local Arrays ..
      DOUBLE PRECISION A(LDA,NMAX), VS(LDVS,NMAX), WI(NMAX),
     +                 WORK(LWORK), WR(NMAX)
      LOGICAL          BWORK(NMAX)
*     .. External Functions ..
      LOGICAL          SELECT
      EXTERNAL         SELECT
*     .. External Subroutines ..
      EXTERNAL         DGEES, X04CAF
*     .. Executable Statements ..
      WRITE (NOUT,*) 'DGEES Example Program Results'
      WRITE (NOUT,*)
*     Skip heading in data file
      READ (NIN,*)
      READ (NIN,*) N
      IF (N.LE.NMAX) THEN
*
*        Read the matrix A from data file
*
         READ (NIN,*) ((A(I,J),J=1,N),I=1,N)
*
*        Find the Schur factorization
*
         CALL DGEES('Vectors (Schur)','Sort',SELECT,N,A,LDA,SDIM,WR,WI,
     +              VS,LDVS,WORK,LWORK,BWORK,INFO)
         LWKOPT = WORK(1)
*
         IF (INFO.EQ.0 .OR. INFO.EQ.(N+2)) THEN
*
*           Print solution
*
            WRITE (NOUT,99999)
     +        'Number of eigenvalues for which SELECT is true = ', SDIM
            WRITE (NOUT,*)
            IF (INFO.EQ.(N+2)) THEN
               WRITE (NOUT,99998) '***Note that rounding errors mean ',
     +           'that leading eigenvalues in the Schur form',
     +           'no longer satisfy SELECT = .TRUE.'
               WRITE (NOUT,*)
            END IF
*
*           Print out factors of the Schur factorization
*
            IFAIL = 0
            CALL X04CAF('General',' ',N,N,A,LDA,'Schur matrix T',IFAIL)
*
            WRITE (NOUT,*)
            CALL X04CAF('General',' ',N,N,VS,LDVS,
     +                  'Matrix of Schur vectors Z',IFAIL)
         ELSE
            WRITE (NOUT,99997) 'Failure in DGEES.  INFO = ', INFO
         END IF
*
*        Print workspace information
*
         IF (LWORK.LT.LWKOPT) THEN
            WRITE (NOUT,*)
            WRITE (NOUT,99996) 'Optimum workspace required = ', LWKOPT,
     +        'Workspace provided         = ', LWORK
         END IF
      ELSE
         WRITE (NOUT,*)
         WRITE (NOUT,*) 'NMAX too small'
      END IF
      STOP
*
99999 FORMAT (1X,A,I4)
99998 FORMAT (1X,2A,/1X,A)
99997 FORMAT (1X,A,I4)
99996 FORMAT (1X,A,I5,/1X,A,I5)
      END

      LOGICAL FUNCTION SELECT(AR,AI)
*     .. Scalar Arguments ..
*
*     Logical function SELECT for use with DGEES
*
*     Returns the value .TRUE. if the imaginary part of the eigenvalue
*     (AR + AI*i) is zero, i.e. the eigenvalue is real
*
      DOUBLE PRECISION        AI, AR
*     .. Local Scalars ..
      LOGICAL                 D
*     .. Executable Statements ..
      IF (AI.EQ.0.0D0) THEN
         D = .TRUE.
      ELSE
         D = .FALSE.
      END IF
*
      SELECT = D
*
      RETURN
      END

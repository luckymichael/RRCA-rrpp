C  ------------------------------------------------------------------  C
C  --------------  Compute and Invert Kriging Matrix  ---------------  C
C  ------------------------------------------------------------------  C
      SUBROUTINE KRIGMATINV(X0,Y0,XP,YP,AMAT,AINV,N,IDR)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (MAXN=2048)
      INTEGER*4 INDX(MAXN)
      DIMENSION XP(N),YP(N)
      DIMENSION AMAT(N+2*IDR+1,N+2*IDR+1),AINV(N+2*IDR+1,N+2*IDR+1)

      ND = N+2*IDR+1
      if (ND.gt.MAXN) call ERROR('Kriging matrix too large')
      call KRIGMAT(X0,Y0,XP,YP,AMAT,N,IDR)
C  Invert A by LU Decomposition
      if (LUDCMP(AMAT,ND,INDX,AINV).EQ.0)
     |  call ERROR('Cannot invert kriging matrix')
      do I=1,ND
        do J=1,ND
           AINV(J,I) = 0
        enddo
        AINV(I,I) = 1
        call LUBKSB(AMAT,ND,INDX,AINV(1,I))
      enddo
      END
C  ------------------------------------------------------------------  C
C  --------  Evaluate Kriging given Inverse Kriging Matrix  ---------  C
C  ------------------------------------------------------------------  C
      FUNCTION EVALKRIGE(AINV,X0,Y0,XP,YP,ZP,X,Y,N,IDR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXND=1024)
      DIMENSION XP(N),YP(N),ZP(N)
      DIMENSION AINV(N+2*IDR+1,N+2*IDR+1),S(MAXND)

      ND = N+2*IDR+1
      if (ND.gt.MAXND) call ERROR('N too large in EVALKRIGE')
C  Calculate S
      do K=1,N
        DX = X-XP(K)
        DY = Y-YP(K)
        S(K) = SQRT(DX**2 + DY**2)
      enddo
      S(N+1) = 1
      do K=1,IDR
        S(N+2*K+0) = (X-X0)**K
        S(N+2*K+1) = (Y-Y0)**K
      enddo
C  Evaluate Function
      EVALKRIGE = 0
      do I=1,N
        W = 0
        do J=1,ND
          W = W + AINV(J,I)*S(J)
        enddo
        EVALKRIGE = EVALKRIGE + W*ZP(I)
      enddo
      END
C  ------------------------------------------------------------------  C
C  ---------------------  Compute Kriging Matrix  -------------------  C
C  ------------------------------------------------------------------  C
      SUBROUTINE KRIGMAT(X0,Y0,XP,YP,A,N,IDR)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION  XP(N),YP(N),A(N+2*IDR+1,N+2*IDR+1)

C  Build coefficient matrix A - assumed linear semivariogram
      do I=1,N
        A(I,I)   = 0
        A(N+1,I) = 1
        A(I,N+1) = 1
        do J=1,I-1
          DX = XP(I)-XP(J)
          DY = YP(I)-YP(J)
          H = SQRT(DX**2 + DY**2)
          if (H.eq.0) call ERROR('Repeated point')
          A(I,J) = H
          A(J,I) = A(I,J)
        enddo
      enddo
C  Drift
      do K=1,IDR
        do I=1,N
          XK = (XP(I)-X0)**K
          YK = (YP(I)-Y0)**K
          A(I,N+2*K+0) = XK
          A(I,N+2*K+1) = YK
          A(N+2*K+0,I) = XK
          A(N+2*K+1,I) = YK
        enddo
      enddo
C  Zero block
      do J=N+1,N+2*IDR+1
        do I=N+1,N+2*IDR+1
          A(I,J) = 0.0
        enddo
      enddo
      END
C  ------------------------------------------------------------------  C
C  ------------------------  Krige to Point  ------------------------  C
C  ------------------------------------------------------------------  C
      FUNCTION POINTKRIGE(X,Y,X0,Y0,XK,YK,ZK,N,IDR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXN=2048)
      DIMENSION XK(N),YK(N),ZK(N)
      DIMENSION A(MAXN,MAXN),W(MAXN)

      if (N+2*IDR+1.GT.MAXN) call ERROR('POINTKRIGE too many points')

C  Construct kriging matrix
      call KRIGMAT(X0,Y0,XK,YK,A,N,IDR)
C  Evaluate RHS
      do K=1,N
        DX = X-XK(K)
        DY = Y-YK(K)
        W(K) = SQRT(DX**2 + DY**2)
      enddo
      W(N+1) = 1
      do K=1,IDR
        W(N+2*K+0) = (X-X0)**K
        W(N+2*K+1) = (Y-Y0)**K
      enddo
C  Solve in place
      call LUSOLV(A,W,N+2*IDR+1)
C  Calculate weighted average
      F = 0
      do L=1,N
        F = F + W(L)*ZK(L)
      enddo
      POINTKRIGE = F
      END
C  ------------------------------------------------------------------- C
C  ------  Solve in place by LU decomposition  (Pivoted Crout)  ------ C
C  ------------------------------------------------------------------- C
      SUBROUTINE LUSOLV(A,X,N)
      IMPLICIT  DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXN=2048)
      DIMENSION A(N,N),X(N),W(MAXN)

      if (N.GT.MAXN) call ERROR('LUSOLV buffer overflow')

C  Decompose
      do i=1,N
        dmax = 0
        do j=1,N
          dmax = MAX(dmax,ABS(A(i,j)))
        enddo
        if (dmax.eq.0) RETURN
        W(i) = 1/dmax
      enddo
      do j=1,N
        do i=1,j-1
          sum = A(i,j)
          do k=1,i-1
            sum = sum - A(i,k)*A(k,j)
          enddo
          A(i,j) = sum
        enddo
        imax = j
        dmax = 0
        do i=j,N
          sum = A(i,j)
          do k=1,j-1
            sum = sum - A(i,k)*A(k,j)
          enddo
          A(i,j) = sum
          dum = W(i)*ABS(sum)
          if (dum.gt.dmax) then
            imax  = i
            dmax = dum
          endif
        enddo
        if (j.ne.imax)then
          do k=1,N
            dum       = A(imax,k)
            A(imax,k) = A(j,k)
            A(j,k)    = dum
          enddo
          dum     = X(imax)
          X(imax) = X(j)
          X(j)    = dum
          W(imax) = W(j)
        endif
        if (A(j,j).eq.0) call ERROR('LUSOLV singular')
        dum = 1/A(j,j)
        do i=j+1,N
          A(i,j) = A(i,j)*dum
        enddo
      enddo

C  Backsolve
      do i=1,N
        do j=1,i-1
          X(i) = X(i) - A(i,j)*X(j)
        enddo
      enddo
      do i=N,1,-1
        do j=i+1,N
          X(i) = X(i)-A(i,j)*X(j)
        enddo
        X(i) = X(i)/A(i,i)
      enddo
      END
C  ------------------------------------------------------------------- C
C  -----------  LU decomposition of matrix (Pivoted Crout)  ---------- C
C  ------------------------------------------------------------------- C
      FUNCTION LUDCMP(A,N,INDX,W)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER        (tiny=1.0E-20)
      INTEGER          N,INDX(N)
      DOUBLE PRECISION A(N,N),W(N)

      LUDCMP = 0
      do i=1,N
        dmax = 0
        do j=1,N
          dmax = MAX(dmax,ABS(A(i,j)))
        enddo
        if (dmax.eq.0) RETURN
        W(i) = 1/dmax
      enddo
      do j=1,N
        do i=1,j-1
          sum = A(i,j)
          do k=1,i-1
            sum = sum - A(i,k)*A(k,j)
          enddo
          A(i,j) = sum
        enddo
        imax = j
        dmax = 0
        do i=j,N
          sum = A(i,j)
          do k=1,j-1
            sum = sum - A(i,k)*A(k,j)
          enddo
          A(i,j) = sum
          dum = W(i)*ABS(sum)
          if (dum.gt.dmax) then
            imax  = i
            dmax = dum
          endif
        enddo
        if (j.ne.imax)then
          do k=1,N
            dum       = A(imax,k)
            A(imax,k) = A(j,k)
            A(j,k)    = dum
          enddo
          W(imax) = W(j)
        endif
        INDX(j) = imax
        if (j.ne.N)then
          if (A(j,j).eq.0) A(j,j) = tiny
          dum = 1/A(j,j)
          do i=j+1,N
            A(i,j) = A(i,j)*dum
          enddo
        endif
      enddo
      if (A(N,N).eq.0) A(N,N) = tiny
      LUDCMP = 1
      END
C  ------------------------------------------------------------------- C
C  ------------------  Backsolve LU decomposition  ------------------- C
C  ------------------------------------------------------------------- C
      SUBROUTINE LUBKSB(A,N,INDX,B)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER          N,INDX(N)
      DOUBLE PRECISION A(N,N),B(N)

      ii = 0
      do i=1,N
        k = INDX(i)
        sum = B(k)
        B(k) = B(i)
        if (ii.ne.0) then
          do j=ii,i-1
            sum = sum - A(i,j)*B(j)
          enddo
        else if (sum.ne.0) then
          ii = i
        endif
        B(i) = sum
      enddo
      do i=N,1,-1
        sum = B(i)
        do j=i+1,N
          sum=sum-A(i,j)*B(j)
        enddo
        B(i)=sum/A(i,i)
      enddo
      END

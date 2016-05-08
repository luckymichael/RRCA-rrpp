! This program reads SFI files for a reference and four imact runs
! and then performs the accounting calculations
! Willem A. Schreuder
! May 17, 2003
!
! Modified 10/18/2013 for 5 Run procedure

      
! Modified by Michael Ou 1/29/2016 takes file names for 5 runs
      PROGRAM ACCT
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INCLUDE 'acct.ins'
      INTEGER        CCP
      INTEGER        IREF
      DIMENSION      VOLM(MAXACCT,0:31),SMS(6),SUM(32),VV(32)
      CHARACTER*4    YEAR
      CHARACTER*21   GAGE(MAXACCT,NGAGE)
      CHARACTER*24   NAME(MAXACCT)
      CHARACTER*256  inFILE,accFILE,sfiFILE(0:31),csvFILE
      CHARACTER*256  CHARIN
      CHARACTER*1024 LINE
!  If the first character of NAME is * it is on the mainstem
!  The first character of GAGE must be +/- to decide if it is
!  added or subtracted in the accumulation or ' ' to omit

!  Process command line arguments
      
      K = 1
      NK = 6
100   if (IARGC().lt.K) call ERROR('Usage: rrca-acct32 file.in')
      call GETARG(K,inFILE)
      write(*,'(A,A)') 'Reading file names of 32-Run SFI files and the Output CSV file from ',inFILE
      
      OPEN(8,FILE=inFILE,STATUS='OLD',IOSTAT=IOS)
      if (IOS.ne.0) call ERROR('Cannot open file '//inFILE)
      do
        read(8,"(A)") accFILE
        if (accFILE(1:1) /= '#') exit
      end do
      do L = 0, 31
        read(8,"(A)") sfiFILE(L)
      end do 
      read(8,"(A)") csvFILE
      close(8)
        
!  Read parameter file
      OPEN(8,FILE=accFILE,STATUS='OLD',IOSTAT=IOS)
      if (IOS.ne.0) call ERROR('Cannot open file '//accFILE)
      NSKIP = 1
      NACCT = 0
      READ(8,'(A)',IOSTAT=IOS) LINE
      do while (IOS.eq.0)
        NACCT = NACCT+1
        if (NACCT.gt.MAXACCT) call ERROR('Too many accounting points')
        NAME(NACCT) = CHARIN(LINE,1,IPOS,IOS)
        if (IOS.ne.0) call ERROR('Error reading account name '//LINE)
        do J=1,NGAGE
          GAGE(NACCT,J) = CHARIN(LINE,IPOS,IPOS,IOS)
        enddo
        READ(8,'(A)',IOSTAT=IOS) LINE
      enddo
      CLOSE(8)
!  Read reference case and four impact runs
      do L=0,31
        call READSFI(sfiFILE(L),GAGE,VOLM(:,L),NACCT,NSKIP)   
      end do
!  Initialize sums
      SMS = 0
      SUM = 0
      
!  Write accounting table
      OPEN(8,FILE=csvFILE,STATUS='UNKNOWN')
      WRITE(8,'(A)') 'Location,Nebraska,LRNRD,MRNRD,URNRD,TRINRD,OTHER.NRD,&
          &THETA,LMUTO,LMUT,LMUO,LOMT,LOTU,MUTO,LUO,LMU,LMO,LTO,LMT,LUT,UTO,MTO,& 
          &MUT,MUO,LM,LU,LT,LO,MU,MT,MO,UT,UO,TO,L,M,U,T,O'
!  Entries
      VV=0.d0
      do K=1,NACCT
        !Nebraska
        VV(1)=VOLM(K,0)-VOLM(K,1)
        !Lower Republican NRD Impact
        VV(2)=(VOLM(K,0)-VOLM(K,27))/5+ &
          ((VOLM(K,28)-VOLM(K,17))/4+(VOLM(K,29)-VOLM(K,18))/4+(VOLM(K,30)-VOLM(K,19))/4+(VOLM(K,31)-VOLM(K,20))/4)/5+ &
          ((VOLM(K,25)-VOLM(K,7))/6+(VOLM(K,21)-VOLM(K,8))/6+(VOLM(K,23)-VOLM(K,9))/6+(VOLM(K,26)-VOLM(K,10))/6+(VOLM(K,22)-VOLM(K,11))/6+(VOLM(K,24)-VOLM(K,12))/6)/5+ &
          ((VOLM(K,15)-VOLM(K,2))/4+(VOLM(K,16)-VOLM(K,3))/4+(VOLM(K,14)-VOLM(K,4))/4+(VOLM(K,13)-VOLM(K,5))/4)/5+ &
          (VOLM(K,6)-VOLM(K,1))/5
        !Middle Republican NRD Impact
        VV(3)=(VOLM(K,0)-VOLM(K,28))/5+ &
          ((VOLM(K,27)-VOLM(K,17))/4+(VOLM(K,29)-VOLM(K,21))/4+(VOLM(K,30)-VOLM(K,22))/4+(VOLM(K,31)-VOLM(K,23))/4)/5+ &
          ((VOLM(K,25)-VOLM(K,16))/6+(VOLM(K,18)-VOLM(K,8))/6+(VOLM(K,20)-VOLM(K,9))/6+(VOLM(K,26)-VOLM(K,14))/6+(VOLM(K,19)-VOLM(K,11))/6+(VOLM(K,24)-VOLM(K,15))/6)/5+ &
          ((VOLM(K,12)-VOLM(K,2))/4+(VOLM(K,7)-VOLM(K,3))/4+(VOLM(K,10)-VOLM(K,4))/4+(VOLM(K,13)-VOLM(K,6))/4)/5+ &
          (VOLM(K,5)-VOLM(K,1))/5
        !Upper Republican NRD Impact
        VV(4)=(VOLM(K,0)-VOLM(K,29))/5+ &
          ((VOLM(K,27)-VOLM(K,18))/4+(VOLM(K,28)-VOLM(K,21))/4+(VOLM(K,30)-VOLM(K,24))/4+(VOLM(K,31)-VOLM(K,25))/4)/5+ &
          ((VOLM(K,23)-VOLM(K,16))/6+(VOLM(K,17)-VOLM(K,8))/6+(VOLM(K,19)-VOLM(K,12))/6+(VOLM(K,22)-VOLM(K,15))/6+(VOLM(K,26)-VOLM(K,13))/6+(VOLM(K,20)-VOLM(K,7))/6)/5+ &
          ((VOLM(K,11)-VOLM(K,2))/4+(VOLM(K,9)-VOLM(K,3))/4+(VOLM(K,10)-VOLM(K,5))/4+(VOLM(K,14)-VOLM(K,6))/4)/5+ &
          (VOLM(K,4)-VOLM(K,1))/5
        !Tri Basin NRD Impact
        VV(5)=(VOLM(K,0)-VOLM(K,30))/5+ &
          ((VOLM(K,27)-VOLM(K,19))/4+(VOLM(K,28)-VOLM(K,22))/4+(VOLM(K,29)-VOLM(K,24))/4+(VOLM(K,31)-VOLM(K,26))/4)/5+ &
          ((VOLM(K,20)-VOLM(K,10))/6+(VOLM(K,17)-VOLM(K,11))/6+(VOLM(K,18)-VOLM(K,12))/6+(VOLM(K,21)-VOLM(K,15))/6+(VOLM(K,25)-VOLM(K,13))/6+(VOLM(K,23)-VOLM(K,14))/6)/5+ &
          ((VOLM(K,8)-VOLM(K,2))/4+(VOLM(K,9)-VOLM(K,4))/4+(VOLM(K,7)-VOLM(K,5))/4+(VOLM(K,16)-VOLM(K,6))/4)/5+ &
          (VOLM(K,3)-VOLM(K,1))/5
        !Other NRD Impacts
        VV(6)=(VOLM(K,0)-VOLM(K,31))/5+ &
          ((VOLM(K,27)-VOLM(K,20))/4+(VOLM(K,28)-VOLM(K,23))/4+(VOLM(K,29)-VOLM(K,25))/4+(VOLM(K,30)-VOLM(K,26))/4)/5+ &
          ((VOLM(K,18)-VOLM(K,7))/6+(VOLM(K,17)-VOLM(K,9))/6+(VOLM(K,19)-VOLM(K,10))/6+(VOLM(K,22)-VOLM(K,14))/6+(VOLM(K,24)-VOLM(K,13))/6+(VOLM(K,21)-VOLM(K,16))/6)/5+ &
          ((VOLM(K,8)-VOLM(K,3))/4+(VOLM(K,11)-VOLM(K,4))/4+(VOLM(K,12)-VOLM(K,5))/4+(VOLM(K,15)-VOLM(K,6))/4)/5+ &
          (VOLM(K,2)-VOLM(K,1))/5

        do L=1,NK
    
          if (NAME(K)(1:1).eq.'*') SMS(L) = SMS(L) + VV(L)
          if (NAME(K)(1:1).ne.'0') SUM(L) = SUM(L)+VV(L)
        enddo
        !print '(I4,4F10.0)', K,VV(1:4)
        WRITE(8,'(A,38(A1,I0))') NAME(K)(2:LENGTH(NAME(K))),((',',IREP(VV(L))),L=1,NK),((',',IREP(VOLM(K,L))),L=0,31)
      enddo
!  Mainstem
      WRITE(8,"(A,6(A1,I0),32(','))") 'Mainstem',((',',IREP(SMS(L))),L=1,NK)
!  Total
      WRITE(8,"(A,6(A1,I0),32(','))") 'Total',((',',IREP(SUM(L))),L=1,NK)
      
      CLOSE(8)
      END
!  ------------------------------------------------------------------  !
!  ------------------------  Read HYDMOD file  ----------------------  !
!  ------------------------------------------------------------------  !
!  This subroutine reads the HYDMOD file and calculates volumes for one run
!  The symbolic names are mapped to columns based on directory entries
!  A bunch of steps are skipped before calculations begin
!  Each year is integrated based on actual times for all the stations
!  The gage values are then accumulated over the accounting items
      SUBROUTINE READSFI(FILE,GAGE,VOLM,NACCT,NSKIP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INCLUDE 'acct.ins'
      PARAMETER     (MAXN=256)
      DIMENSION     INDEX(MAXACCT,NGAGE),VOLM(MAXACCT)
      DIMENSION     SUM(MAXN),Q(MAXN)
      CHARACTER*20  DIR(MAXN)
      CHARACTER*21  GAGE(MAXACCT,NGAGE)
      CHARACTER*(*) FILE
      INCLUDE 'openspec.inc'
      

!  Open hydmod file and read directory
      OPEN(8,FILE=FILE,STATUS='OLD',FORM=FORM,ACCESS=ACCESS,IOSTAT=IOS)
      if (IOS.ne.0) call ERROR2('Cannot open file',FILE)
      READ(8,IOSTAT=IOS) N,ITMUNI
      if (IOS.ne.0) call ERROR2('Cannot read count',FILE)
      if (N.gt.MAXN) call ERRORINT('Too many HYDMOD items',N)
      READ(8,IOSTAT=IOS) IDUMMY,(DIR(K),K=1,N)
      if (IOS.ne.0) call ERROR2('Cannot read directory',FILE)
!  Match gages to directory entries
      do I=1,NACCT
        do J=1,NGAGE
          if (GAGE(I,J).eq.' ') then
            INDEX(I,J) = 0
          else
            L = 0
            do K=1,N
              if (GAGE(I,J)(2:21).eq.DIR(K)) L = K
            enddo
            if (L.eq.0) then
              call ERROR2('No match to gage',GAGE(I,J))
            else if (GAGE(I,J)(1:1).eq.'+') then
              INDEX(I,J) = +L
            else if (GAGE(I,J)(1:1).eq.'-') then
              INDEX(I,J) = -L
            else
              call ERROR2('Invalid gage type',GAGE(I,J))
            endif
          endif
        enddo
      enddo
!  Skip initial, steady state and initial year records
      do I=1,NSKIP
        READ(8) T
      enddo
!  Compute annual values
!       Initialize accumulators
      do K=1,N
        SUM(K) = 0
      enddo
!  Read 12 months by 2 time steps
!  Integrate rates into sum
      do L=1,2*12
        T0 = T
        READ(8,IOSTAT=IOS) T,(Q(K),K=1,N)
        if (IOS.ne.0) call ERROR2('Error reading file',FILE)
        FACT = (T-T0)/43560
        do K=1,N
          SUM(K) = SUM(K) + FACT*Q(K)
        enddo
      enddo
!  Map annual totals to gage and sum accounting points
      do I=1,NACCT
        VOLM(I) = 0
        do J=1,NGAGE
          if (INDEX(I,J).gt.0) then
            VOLM(I) = VOLM(I) + SUM(INDEX(I,J))
          else if (INDEX(I,J).lt.0) then
            VOLM(I) = VOLM(I) - SUM(-INDEX(I,J))
          endif
        enddo
      enddo
      CLOSE(8)
      END
!  ------------------------------------------------------------------  !
!  -----------------------  Round for Reporting  --------------------  !
!  ------------------------------------------------------------------  !
      PURE FUNCTION IREP(X)
      DOUBLE PRECISION, INTENT(IN) :: X
      INTEGER :: IREP
      INTEGER, PARAMETER :: REPMIN=10
!  Do not report values less than REPMIN
      if (ABS(X).lt.REPMIN) then
        IREP = 0
!  Round to positive integers
      else if (X.gt.0) then
        IREP = X+0.5
!  Round to negative integers
      else
        IREP = X-0.5
      endif
      END

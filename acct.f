C This program reads SFI files for a reference and four imact runs
C and then performs the accounting calculations
C Willem A. Schreuder
C May 17, 2003
C
C Modified 10/18/2013 for 5 Run procedure
      PROGRAM ACCT
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INCLUDE 'acct.ins'
      INTEGER        CCP
      INTEGER        IREF
      DIMENSION      VOLM(MAXACCT,0:5),SMS(5),SUM(5)
      CHARACTER*4    YEAR
      CHARACTER*21   GAGE(MAXACCT,NGAGE)
      CHARACTER*24   NAME(MAXACCT)
      CHARACTER*32   FILE
      CHARACTER*256  CHARIN
      CHARACTER*1024 LINE
C  If the first character of NAME is * it is on the mainstem
C  The first character of GAGE must be +/- to decide if it is
C  added or subtracted in the accumulation or ' ' to omit

C  Process command line arguments
      FILE = '../data0/acct.12s2'
      IREF = 0
      CCP = 0
      K = 1
      NK = 4
100   if (IARGC().lt.K) call ERROR('Usage: acct year')
      call GETARG(K,YEAR)
C     Alternative control file
      if (YEAR.eq.'-f') then
        call GETARG(K+1,FILE)
        K = K+2
        goto 100
C     5 Run procedure (IREF=d)
      else if (YEAR.eq.'-5') then
        IREF = 4
        K = K+1
        goto 100
C     Kansas Method for CCP
      else if (YEAR.eq.'-c') then
        CCP = 1
        NK = 5
        K = K+1
        goto 100
C     Colorado Method for CCP
      else if (YEAR.eq.'-C') then
        CCP = -1
        call GETARG(K+1,LINE)
        READ(LINE,'(F20.0)') VCPP
        K = K+2
        NK = 5
        goto 100
      endif
C  Read parameter file
      OPEN(8,FILE=FILE,STATUS='OLD',IOSTAT=IOS)
      if (IOS.ne.0) call ERROR('Cannot open file '//FILE)
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
C  Read reference case and four impact runs
      call READSFI(YEAR//'.sfi',GAGE,VOLM(1,0),NACCT,NSKIP)
      call READSFI(YEAR//'a.sfi',GAGE,VOLM(1,1),NACCT,NSKIP)
      call READSFI(YEAR//'b.sfi',GAGE,VOLM(1,2),NACCT,NSKIP)
      call READSFI(YEAR//'c.sfi',GAGE,VOLM(1,3),NACCT,NSKIP)
      call READSFI(YEAR//'d.sfi',GAGE,VOLM(1,4),NACCT,NSKIP)
C  Initialize sums
      do L=1,NK
        SMS(L) = 0
        SUM(L) = 0
      enddo
C  Process CCP
      if (CCP.gt.0) then
        call READSFI(YEAR//'e.sfi',GAGE,VOLM(1,5),NACCT,NSKIP)
      else
        do K=1,NACCT
          VOLM(K,5) = VOLM(K,IREF)
        enddo
        VOLM(6,5) = VOLM(6,IREF)+VCPP
      endif
C  Write accounting table
      OPEN(8,FILE=YEAR//'.htm',STATUS='UNKNOWN')
      WRITE(8,'(A)')  '<html>'
      WRITE(8,'(3A)') '<head><title>Impacts ',YEAR,'</title></head>'
      WRITE(8,'(A)')  '<body bgcolor="#FFFFFF">'
C  Header
      WRITE(8,'(A)') '<table border>'
      WRITE(8,'(A,I1,3A)') '<tr><th align=center colspan=',NK+1,
     |       '><font size="+4">Impacts ',YEAR,
     |       ' (acre-feet)</font></th></tr>'
      if (CCP.ne.0) then
        WRITE(8,'(8A)') '<tr align="center">',
     |      '<th>Location</th>',
     |      '<th>Colorado<br>Pumping</th>',
     |      '<th>Kansas<br>Pumping</th>',
     |      '<th>Nebraska<br>Pumping</th>',
     |      '<th>Nebraska<br>Mound</th>',
     |      '<th>CCP<br>Credit</th>',
     |      '</tr>'
      else
        WRITE(8,'(7A)') '<tr align="center">',
     |      '<th>Location</th>',
     |      '<th>Colorado<br>Pumping</th>',
     |      '<th>Kansas<br>Pumping</th>',
     |      '<th>Nebraska<br>Pumping</th>',
     |      '<th>Nebraska<br>Mound</th>',
     |      '</tr>'
      endif
C  Entries
      do K=1,NACCT
        WRITE(8,'(4A)') '<tr align="right">',
     |    '<td align=left>',NAME(K)(2:LENGTH(NAME(K))),'</td>'
        do L=1,NK
          if (L.eq.4) then
            VOL = VOLM(K,0)-VOLM(K,L)
          else
            VOL = VOLM(K,L)-VOLM(K,IREF)
          endif
          if (NAME(K)(1:1).eq.'*') SMS(L) = SMS(L) + VOL
          SUM(L) = SUM(L)+VOL
          WRITE(8,'(A,I10,A)') '<td>',IREP(VOL),'</td>'
        enddo
        WRITE(8,'(A)')'</tr>'
      enddo
C  Mainstem
      WRITE(8,'(4A)') '<tr align="right">',
     |  '<th align=left>Mainstem</th>'
      do L=1,NK
        WRITE(8,'(A,I10,A)') '<th>',IREP(SMS(L)),'</th>'
      enddo
      WRITE(8,'(A)')'</tr>'
C  Total
      WRITE(8,'(4A)') '<tr align="right">',
     |  '<th align=left>Total</th>'
      do L=1,NK
        WRITE(8,'(A,I10,A)') '<th>',IREP(SUM(L)),'</th>'
      enddo
      WRITE(8,'(A)')'</tr>'
      WRITE(8,'(A)') '</table>'
      if (CCP.gt.0) then
         WRITE(8,'(A)') '<H3>Kansas Method</H3>'
      else if (IREF.eq.0) then
         WRITE(8,'(A)') '<H3>Original Procedure</H3>'
      else
         WRITE(8,'(A)') '<H3>5 Run Procedure</H3>'
      endif
      WRITE(8,'(A)') '</body>'
      WRITE(8,'(A)') '</html>'
      CLOSE(8)
      END
C  ------------------------------------------------------------------  C
C  ------------------------  Read HYDMOD file  ----------------------  C
C  ------------------------------------------------------------------  C
C  This subroutine reads the HYDMOD file and calculates volumes for one run
C  The symbolic names are mapped to columns based on directory entries
C  A bunch of steps are skipped before calculations begin
C  Each year is integrated based on actual times for all the stations
C  The gage values are then accumulated over the accounting items
      SUBROUTINE READSFI(FILE,GAGE,VOLM,NACCT,NSKIP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INCLUDE 'acct.ins'
      PARAMETER     (MAXN=256)
      DIMENSION     INDEX(MAXACCT,NGAGE),VOLM(MAXACCT)
      DIMENSION     SUM(MAXN),Q(MAXN)
      CHARACTER*20  DIR(MAXN)
      CHARACTER*21  GAGE(MAXACCT,NGAGE)
      CHARACTER*(*) FILE

C  Open hydmod file and read directory
      OPEN(8,FILE=FILE,STATUS='OLD',FORM='UNFORMATTED',IOSTAT=IOS)
      if (IOS.ne.0) call ERROR2('Cannot open file',FILE)
      READ(8,IOSTAT=IOS) N,ITMUNI
      if (IOS.ne.0) call ERROR2('Cannot read count',FILE)
      if (N.gt.MAXN) call ERRORINT('Too many HYDMOD items',N)
      READ(8,IOSTAT=IOS) IDUMMY,(DIR(K),K=1,N)
      if (IOS.ne.0) call ERROR2('Cannot read directory',FILE)
C  Match gages to directory entries
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
C  Skip initial, steady state and initial year records
      do I=1,NSKIP
        READ(8) T
      enddo
C  Compute annual values
C       Initialize accumulators
      do K=1,N
        SUM(K) = 0
      enddo
C  Read 12 months by 2 time steps
C  Integrate rates into sum
      do L=1,2*12
        T0 = T
        READ(8,IOSTAT=IOS) T,(Q(K),K=1,N)
        if (IOS.ne.0) call ERROR2('Error reading file',FILE)
        FACT = (T-T0)/43560
        do K=1,N
          SUM(K) = SUM(K) + FACT*Q(K)
        enddo
      enddo
C  Map annual totals to gage and sum accounting points
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
C  ------------------------------------------------------------------  C
C  -----------------------  Round for Reporting  --------------------  C
C  ------------------------------------------------------------------  C
      FUNCTION IREP(X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (REPMIN=10)
C  Do not report values less than REPMIN
      if (ABS(X).lt.REPMIN) then
        IREP = 0
C  Round to positive integers
      else if (X.gt.0) then
        IREP = X+0.5
C  Round to negative integers
      else
        IREP = X-0.5
      endif
      END

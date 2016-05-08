C This program reads SFI files for a reference and four imact runs
C and then performs the accounting calculations
C Willem A. Schreuder
C May 17, 2003
      PROGRAM ACCT
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INCLUDE 'acct.ins'
      DIMENSION      VOLM(MAXACCT,MAXYEAR,5)
      CHARACTER*21   GAGE(MAXACCT,NGAGE)
      CHARACTER*24   NAME(MAXACCT)
      CHARACTER*256  FILE,CHARIN
      CHARACTER*1024 LINE
C  If the first character of NAME is * it is on the mainstem
C  The first character of GAGE must be +/- to dcide if it is
C  added or subtracted in the accumulation or ' ' to omit

      if (IARGC().ne.7)
     |  call ERROR('Usage: acct par ref co ks ne mound ver')
      call GETARG(1,FILE)
      OPEN(8,FILE=FILE,STATUS='OLD',IOSTAT=IOS)
      if (IOS.ne.0) call ERROR('Cannot open file '//FILE)
      READ(8,'(A)',IOSTAT=IOS) LINE
      NSKIP = INTIN(LINE,   1,IPOS,IOS)
      if (NSKIP.lt.0.or.IOS.ne.0) call ERROR('Invalid SKIP '//LINE)
      IYR0  = INTIN(LINE,IPOS,IPOS,IOS)-1
      if (IYR0.le.0.or.IOS.ne.0) call ERROR('Invalid start year '//LINE)
      IYR1  = INTIN(LINE,IPOS,IPOS,IOS)
      if (IYR1.le.0.or.IOS.ne.0) call ERROR('Invalid end year '//LINE)
      NYEAR = IYR1-IYR0
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
C  Read reference case and four impact runs
      do L=1,5
        call GETARG(L+1,FILE)
        call READSFI(FILE,GAGE,VOLM(1,1,L),NACCT,NYEAR,NSKIP)
      enddo
C  Write accountinf table
      call GETARG(7,FILE)
      L = LENGTH(FILE)
      WRITE(6,'(A)')  '<html>'
      WRITE(6,'(3A)') '<head><title>Model Version ',FILE(1:L),
     |                '</title></head>'
      WRITE(6,'(A)')  '<body bgcolor="#FFFFFF">'
      WRITE(6,'(3A)') '<H1>Impacts Version ',FILE(1:L),'</H1>'
      call TABLE('Impact of Colorado Pumping',
     |   VOLM(1,1,1),VOLM(1,1,2),NAME,IYR0,NACCT,NYEAR)
      call TABLE('Impact of Kansas Pumping',
     |   VOLM(1,1,1),VOLM(1,1,3),NAME,IYR0,NACCT,NYEAR)
      call TABLE('Impact of Nebraska Pumping',
     |   VOLM(1,1,1),VOLM(1,1,4),NAME,IYR0,NACCT,NYEAR)
      call TABLE('Impact of Mound',
     |   VOLM(1,1,5),VOLM(1,1,1),NAME,IYR0,NACCT,NYEAR)
      WRITE(6,'(A)') '</body>'
      WRITE(6,'(A)') '</html>'
      END
C  ------------------------------------------------------------------  C
C  ------------------------  Read HYDMOD file  ----------------------  C
C  ------------------------------------------------------------------  C
C  This subroutine reads the HYDMOD file and calculates volumes for one run
C  The symbolic names are mapped to columns based on directory entries
C  A bunch of steps are skipped before calculations begin
C  Each year is integrated based on actual times for all the stations
C  The gage values are then accumulated over the accounting items
      SUBROUTINE READSFI(FILE,GAGE,VOLM,NACCT,NYEAR,NSKIP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INCLUDE 'acct.ins'
      PARAMETER     (MAXN=256)
      DIMENSION     INDEX(MAXACCT,NGAGE),VOLM(MAXACCT,MAXYEAR)
      DIMENSION     SUM(MAXN),Q(MAXN),Q0(MAXN)
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
C  The offset-1M gages are all 'SO'
      do K=1,N
        if (DIR(K)(1:2).eq.'SO') then
          Q0(K) = 1000000
        else
          Q0(K) = 0
        endif
      enddo
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
      do IYR=1,NYEAR
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
            SUM(K) = SUM(K) + FACT*(Q(K)-Q0(K))
          enddo
        enddo
C  Map annual totals to gage and sum accounting points
        do I=1,NACCT
          VOLM(I,IYR) = 0
          do J=1,NGAGE
            if (INDEX(I,J).gt.0) then
              VOLM(I,IYR) = VOLM(I,IYR) + SUM(INDEX(I,J))
            else if (INDEX(I,J).lt.0) then
              VOLM(I,IYR) = VOLM(I,IYR) - SUM(-INDEX(I,J))
            endif
          enddo
        enddo
      enddo
      CLOSE(8)
      END
C  ------------------------------------------------------------------  C
C  ----------------------  Save accounting table  -------------------  C
C  ------------------------------------------------------------------  C
C  This subroutine saves one table to an HTML table
C  The impact is the difference between VOL0 and VOL1
      SUBROUTINE TABLE(TITLE,VOL0,VOL1,NAME,IYR0,NACCT,NYEAR)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INCLUDE 'acct.ins'
      DIMENSION     VOL0(MAXACCT,MAXYEAR),VOL1(MAXACCT,MAXYEAR)
      DIMENSION     AVG(MAXACCT)
      CHARACTER*(*) TITLE,NAME(MAXACCT)

C  Header
      WRITE(6,'(A)') '<table border>'
      WRITE(6,'(A,I2,3A)') '<tr align="center"><th colspan=',NACCT+3,
     |  '>',TITLE,' (acre-feet)</th></tr>'
      WRITE(6,'(A)') '<tr align="center"><th>Year</th>'
      do K=1,NACCT
        AVG(K) = 0
        WRITE(6,'(3A)') '<th>',NAME(K)(2:LENGTH(NAME(K))),'</th>'
      enddo
      WRITE(6,'(A)') '<th>Mainstem Total</th><th>Total</th></tr>'
C  Year by year
      do IYR=1,NYEAR
        WRITE(6,'(A,I4,A)') '<tr align="right"><td>',IYR0+IYR,'</td>'
        SMS = 0
        SUM = 0
        do K=1,NACCT
          VOL = VOL1(K,IYR)-VOL0(K,IYR)
          AVG(K) = AVG(K)+VOL
          if (NAME(K)(1:1).eq.'*') SMS = SMS + VOL
          SUM = SUM + VOL
          WRITE(6,'(A,I10,A)') '<td>',IREP(VOL),'</td>'
        enddo
        WRITE(6,'(A,2(I10,A))')'<td>',IREP(SMS),'</td><td>',
     |                                IREP(SUM),'</td></tr>'
      enddo
C  Average
      WRITE(6,'(A,2(I4,A))') '<tr align="right"><td>Average ',
     |  IYR0+1,'-',IYR0+NYEAR,'</td>'
      SMS = 0
      SUM = 0
      do K=1,NACCT
        VOL = AVG(K)/NYEAR
        if (NAME(K)(1:1).eq.'*') SMS = SMS + VOL
        SUM = SUM + VOL
        WRITE(6,'(A,I10,A)') '<td>',IREP(VOL),'</td>'
      enddo
      WRITE(6,'(A,2(I10,A))') '<td>',IREP(SMS),'</td><td>',
     |                               IREP(SUM),'</td></tr>'
      WRITE(6,'(A)') '</table>'
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

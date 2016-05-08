C  ------------------------------------------------------------------- C
C  -----------------  Convert string to upper case  ------------------ C
C  ------------------------------------------------------------------- C
      SUBROUTINE UPPERCASE(STR)
      CHARACTER*(*) STR

      do I=1,LEN(STR)
        if ('a'.le.STR(I:I).and.STR(I:I).le.'z')
     |    STR(I:I) = CHAR(ICHAR(STR(I:I))-32)
      enddo
      END
C  ------------------------------------------------------------------- C
C  --------------  Test if we are at the end of line  ---------------- C
C  ------------------------------------------------------------------- C
      FUNCTION EOL(LINE,IPOS)
      LOGICAL       EOL
      CHARACTER*(*) LINE

      LN = LEN(LINE)
      do while (IPOS.le.LN .and. LINE(IPOS:IPOS).eq.' ')
        IPOS = IPOS+1
      enddo
      EOL = (IPOS.gt.LN)
      END
C  ------------------------------------------------------------------- C
C  -----------------  Read a string, status in IOS  ------------------ C
C  ------------------------------------------------------------------- C
      FUNCTION CHARIN(LINE,ISTART,IEND,IOS)
      CHARACTER*256 CHARIN
      CHARACTER*(*) LINE

      IOS = -1
      CHARIN = ' '
      LN = LEN(LINE)
      IEND = ISTART
      do while (IEND.le.LN .and. LINE(IEND:IEND).eq.' ')
        IEND = IEND+1
      enddo
      if (IEND.gt.LN) RETURN
      if (LINE(IEND:IEND).eq.'"') then
        IEND = IEND+1
        I0 = IEND
        do while (IEND.le.LN .and. LINE(IEND:IEND).ne.'"')
          IEND = IEND+1
        enddo
        if (IEND.gt.LN) RETURN
        I1 = IEND-1
        IEND = IEND+1
      else
        I0 = IEND
        do while (IEND.le.LN .and. LINE(IEND:IEND).ne.' ')
          IEND = IEND+1
        enddo
        I1 = MIN(IEND,LN)
      endif
      CHARIN = LINE(I0:I1)
      IOS= 0
      END
C  ------------------------------------------------------------------- C
C  ----------------  Read an integer, status in IOS  ----------------- C
C  ------------------------------------------------------------------- C
      FUNCTION INTIN(LINE,ISTART,IEND,IOS)
      CHARACTER*(*) LINE

      IOS = -1
      INTIN = 0
      LN = LEN(LINE)
      IEND = ISTART
      do while (IEND.le.LN .and. LINE(IEND:IEND).eq.' ')
        IEND = IEND+1
      enddo
      if (IEND.gt.LN) RETURN
      I0 = IEND
      do while (IEND.le.LN .and. LINE(IEND:IEND).ne.' ')
        IEND = IEND+1
      enddo
      I1 = MIN(IEND,LN)
      READ(LINE(I0:I1),'(I32)',IOSTAT=IOS) INTIN
      END
C  ------------------------------------------------------------------- C
C  ------------------  Read a real, status in IOS  ------------------- C
C  ------------------------------------------------------------------- C
      FUNCTION REALIN(LINE,ISTART,IEND,IOS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*(*) LINE

      IOS = -1
      REALIN = 0
      LN = LEN(LINE)
      IEND = ISTART
      do while (IEND.le.LN .and. LINE(IEND:IEND).eq.' ')
        IEND = IEND+1
      enddo
      if (IEND.gt.LN) RETURN
      I0 = IEND
      do while (IEND.le.LN .and. LINE(IEND:IEND).ne.' ')
        IEND = IEND+1
      enddo
      I1 = MIN(IEND,LN)
      READ(LINE(I0:I1),'(G32.0)',IOSTAT=IOS) REALIN
      END
C  ------------------------------------------------------------------- C
C  ---------------  Return non-blank length of string  --------------- C
C  ------------------------------------------------------------------- C
      FUNCTION LENGTH0(LINE)
      CHARACTER*(*) LINE

      LN = LEN(LINE)
      do I=LN,1,-1
        if (LINE(I:I).ne.' ') goto 999
      enddo
      I = 0
999   LENGTH0 = I
      END
C  ------------------------------------------------------------------- C
C  ---------------  Return non-blank length of string  --------------- C
C  ------------------------------------------------------------------- C
      FUNCTION LENGTH(LINE)
      CHARACTER*(*) LINE

      LN = LEN(LINE)
      do I=LN,1,-1
        if (LINE(I:I).ne.' ') goto 999
      enddo
      I = 1
999   LENGTH = I
      END
C  ------------------------------------------------------------------- C
C  -----------------  Print error message and stop  ------------------ C
C  ------------------------------------------------------------------- C
      SUBROUTINE ERROR(MSG)
      CHARACTER*(*) MSG

      print *,MSG(1:LENGTH(MSG))
      STOP
      END
C  ------------------------------------------------------------------- C
C  -----------------  Print error message and stop  ------------------ C
C  ------------------------------------------------------------------- C
      SUBROUTINE ERROR2(MSG1,MSG2)
      CHARACTER*(*) MSG1,MSG2

      print *,MSG1(1:LENGTH(MSG1)),' ',MSG2(1:LENGTH(MSG2))
      STOP
      END
C  ------------------------------------------------------------------- C
C  -----------------  Print error message and stop  ------------------ C
C  ------------------------------------------------------------------- C
      SUBROUTINE ERRORINT(MSG,IMSG)
      CHARACTER*(*) MSG

      print *,MSG(1:LENGTH(MSG)),' ',IMSG
      STOP
      END
C  ------------------------------------------------------------------- C
C  -----------------  Print error message and stop  ------------------ C
C  ------------------------------------------------------------------- C
      SUBROUTINE ERRORDBL(MSG,DMSG)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*(*) MSG

      print *,MSG(1:LENGTH(MSG)),' ',DMSG
      STOP
      END
C  ------------------------------------------------------------------- C
C  ----------------  Read flag file with header line  ---------------- C
C  ------------------------------------------------------------------- C
C  This subroutine reads a flag array from the file FILE in format 999I1
C  The flags consists of integers 0 or 1-9.
C  0 means outside the domain, 1-9 indicate up to 9 zones
C  The first line must enumerate the names of the zones that will appear
C  Names are delimited by spaces
C  No more tha "MAX" names may appear (this is usually 9)
C  When MAX.eq.0, IBOUNDS are read in format 999I2 and no line of names
C  If the file name is '*' a dummy zone is created for the entire domain
      SUBROUTINE READFLAGS(FILE,FLAG,NROW,NCOL,LAB,N,MAX)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER       FLAG(NCOL,NROW)
      CHARACTER*(*) FILE,LAB(9)
      CHARACTER*256 LINE,CHARIN,STRING

C  Special case - map everything to 1
      if (FILE.eq.'*') then
        N = 1
        LAB(N) = '*'
        do I=1,NROW
          do J=1,NCOL
            FLAG(J,I) = 1
          enddo
        enddo
C  Read actual file
      else if (MAX.gt.0) then
        OPEN(8,FILE=FILE,STATUS='OLD',IOSTAT=IOS)
        if (IOS.ne.0) call ERROR2('Error opening file',FILE)

C     Read header line that defines zone names
        READ(8,'(A)') LINE
        N = 0
        IPOS = 1
        STRING = CHARIN(LINE,IPOS,IPOS,IOS)
        do while (IOS.eq.0)
          N = N+1
          if (N.gt.MAX) call ERROR2('Too many zones in',FILE)
          LAB(N) = STRING
          STRING = CHARIN(LINE,IPOS,IPOS,IOS)
        enddo

C     Read flag zones
        do I=1,NROW
          READ(8,'(999I1)',IOSTAT=IOS) (FLAG(J,I),J=1,NCOL)
          if (IOS.ne.0) call ERROR2('Error reading file',FILE)
        enddo
        CLOSE(8)

C     Check that all flags are in range
        do I=1,NROW
          do J=1,NCOL
            if (FLAG(J,I).lt.0 .or. FLAG(J,I).gt.N)
     |        call ERROR2('Zone out of range in file',FILE)
          enddo
        enddo
C  Read actual IBOUND file
      else
        OPEN(8,FILE=FILE,STATUS='OLD',IOSTAT=IOS)
        if (IOS.ne.0) call ERROR2('Error opening file',FILE)
        do I=1,NROW
          READ(8,'(999I2)',IOSTAT=IOS) (FLAG(J,I),J=1,NCOL)
          if (IOS.ne.0) call ERROR2('Error reading file',FILE)
        enddo
        CLOSE(8)
      endif
      END
C  ------------------------------------------------------------------- C
C  --------------  Read real file free format with skip  ------------- C
C  ------------------------------------------------------------------- C
      SUBROUTINE READREAL(FILE,X,NROW,NCOL,ISKIP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION     X(NCOL,NROW)
      CHARACTER*(*) FILE

      OPEN(8,FILE=FILE,STATUS='OLD',IOSTAT=IOS)
      if (IOS.ne.0) call ERROR2('Error opening file',FILE)
C  Skip headers
      do I=1,ISKIP
        READ(8,*,IOSTAT=IOS)
        if (IOS.ne.0) call ERROR2('Error reading file header',FILE)
      enddo
C  Read values
      do I=1,NROW
        READ(8,*,IOSTAT=IOS) (X(J,I),J=1,NCOL)
        if (IOS.ne.0) call ERROR2('Error reading file',FILE)
      enddo
      CLOSE(8)
      END
C----------------------------------------------------------------------C
C-------------------  Read annual factors from file  ------------------C
C----------------------------------------------------------------------C
C  This subroutine reads a file in a year-values format.
C  The first line must contain a list of the items to be read
C  These names will be check agains the NAME list to make sure
C  the columns are in the right order.
C  The data format is YEAR Col1 Col2 .... ColN
C  If the years is out of range MINYR..MAXYR it is ignored
C  The entire array is set equal to zero initially so
C  the data may be incomplete and out of order
      SUBROUTINE READYEARS(FILE,X,NAME,N,MINYR,MAXYR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION      X(MINYR:MAXYR,N)
      CHARACTER*32   NAME(N)
      CHARACTER*1024 LINE
      CHARACTER*(*)  FILE

C  Initialize
      do L=1,N
        do IYR=MINYR,MAXYR
          X(IYR,L) = 0
        enddo
      enddo
C  Read
      OPEN(8,FILE=FILE,STATUS='OLD',IOSTAT=IOS)
      if (IOS.ne.0) call ERROR2('Cannot open file',FILE)
C  Read names from first line
      READ(8,'(A)',IOSTAT=IOS) LINE
      if (IOS.ne.0) call ERROR2('Cannot read file',FILE)
      call CHECKHDR(FILE,LINE,NAME,N)
C  Read data
100   READ(8,'(A)',END=199) LINE
      if (LINEFIX(LINE).ne.0) goto 100
      IYR = INTIN(LINE,1,IPOS,IOS)
      if (IOS.ne.0) call ERROR2('Invalid year in file',FILE)
      if (IYR.lt.MINYR .or. IYR.gt.MAXYR) goto 100
      do L=1,N
        X(IYR,L) = REALIN(LINE,IPOS,IPOS,IOS)
        if (IOS.ne.0) call ERROR2('Invalid data in file',FILE)
      enddo
      goto 100
199   CLOSE(8)
      END
C----------------------------------------------------------------------C
C--------------------  Store array to monthly file  -------------------C
C----------------------------------------------------------------------C
C     This subroutine stores array data in a file free format unless the
C     output would consist of all zeroes
C     The array to be saved is scaled by FACTM.
C     The file name is created from the IYR, MO and EXT values in
C     a date-lke format YYYY.MM.EXT in a subdirectory DIR.
C     If IYR<=0 the file name is 'steady'.
C     If MO<= the file name omits the month.
C     A binary file is written if BIN is true.
      SUBROUTINE STOREMONTH(X,FACTM,NROW,NCOL,IYR,MO,DIR,EXT,BIN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(NCOL,NROW)
      LOGICAL       BIN
      CHARACTER*(*) EXT,DIR
      CHARACTER*256 FILE

      if (FACTM.eq.0) RETURN

      SUM = 0
      do I=1,NROW
        do J=1,NCOL
          SUM = SUM + ABS(FACTM*X(J,I))
        enddo
      enddo
      if (SUM.eq.0) RETURN

      L = LENGTH(DIR)
      if (IYR.le.0) then
        WRITE(FILE,'(3A)') DIR(1:L),'/steady.',EXT
      else if (MO.le.0) then
        WRITE(FILE,'(2A,I4,2A)') DIR(1:L),'/',IYR,'.',EXT
      else
        WRITE(FILE,'(2A,F7.2,2A)') DIR(1:L),'/',IYR+MO/1D2,'.',EXT
      endif
      if (BIN) then
        OPEN(8,FILE=FILE,FORM='UNFORMATTED',STATUS='UNKNOWN',IOSTAT=IOS)
        if (IOS.ne.0) call ERROR('Cannot open file '//FILE)
        WRITE(8,IOSTAT=IOS) ((FACTM*X(J,I),J=1,NCOL),I=1,NROW)
        if (IOS.ne.0) call ERROR('Cannot write to file '//FILE)
        CLOSE(8)
      else
        OPEN(8,FILE=FILE,STATUS='UNKNOWN',IOSTAT=IOS)
        if (IOS.ne.0) call ERROR('Cannot open file '//FILE)
        do I=1,NROW
          WRITE(8,*,IOSTAT=IOS) (FACTM*X(J,I),J=1,NCOL)
          if (IOS.ne.0) call ERROR('Cannot write to file '//FILE)
        enddo
        CLOSE(8)
      endif
      END
C  ------------------------------------------------------------------- C
C  ----------------  Remove non-printable characters  ---------------- C
C  ------------------------------------------------------------------- C
C     This function returns
C       1 if a comment line
C      -1 if a blank line
C       0 in all other cases
      FUNCTION LINEFIX(LINE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*(*) LINE

C  Change non-printable characters to spaces - assumes ASCII ordering
      do I=1,LEN(LINE)
        if (LINE(I:I).lt.' ' .or. LINE(I:I).gt.'~') LINE(I:I) = ' '
      enddo
C  Check if blank or comment line
      if (LINE(1:1).eq.'#') then
        LINEFIX = 1
      else if (LENGTH0(LINE).eq.0) then
        LINEFIX = -1
      else
        LINEFIX = 0
      endif
      END
C----------------------------------------------------------------------C
C--------------------  Read column names from line  -------------------C
C----------------------------------------------------------------------C
C  Read column names and compare them against the supplied list
C  Names are converted to upper case before comparison.
      SUBROUTINE CHECKHDR(FILE,LINE,NAME,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*32   STRING,NAME(N)
      CHARACTER*256  CHARIN
      CHARACTER*(*)  FILE,LINE

      STRING = CHARIN(LINE,1,IPOS,IOS)
      if (IOS.ne.0) call ERROR2('Invalid time label',FILE)
      do L=1,N
        STRING = CHARIN(LINE,IPOS,IPOS,IOS)
        if (IOS.ne.0) call ERROR2('Invalid data label',FILE)
        call UPPERCASE(STRING)
        if (STRING.ne.NAME(L))
     |    call ERROR2('Data label out of order '//STRING,FILE)
      enddo
      END
C  ------------------------------------------------------------------- C
C  ------------------  Linear or Quadratic spline  ------------------- C
C  ------------------------------------------------------------------- C
C  Fit a linear or quadratic spline through the data
C  (X,Y) are a set of N knots
C  T is a set of M independent variable values
C  F returns the spline value at each T
C  NSPL is the spline order (1 or 2)
      SUBROUTINE SPLINE(X,Y,N,T,F,M,NSPL)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (MAXN=1024)
      DIMENSION X(N),Y(N),T(M),F(M),C1(MAXN),C2(MAXN)

C
C  Compute spline coefficients
C
      if (N.lt.2) call ERROR('Cannot fit spline to one point')
      if (N.gt.MAXN) call ERROR('Too many points to fit spline')
      if (NSPL.eq.1) then
        do I=1,N-1
          C1(I) = (Y(I+1)-Y(I))/(X(I+1)-X(I))
          C2(I) = 0
        enddo
      else
C       Decide on initial slope
        DX  = X(2)-X(1)
        YP = (Y(2)-Y(1))/DX
        C1(1) = YP
        C2(1) = 0
C       Propagate terminal slope
        do I=2,N-1
          DX  = X(I+1)-X(I)
          C1(I) = YP
          C2(I) = (Y(I+1)-Y(I)-YP*DX)/DX**2
          YP = YP + 2*C2(I)*DX
        enddo
      endif
C     Terminal knot
      C1(N) = 0
      C2(N) = 0

C
C  Evaluate spline
C
      do L=1,M
        if (T(L).lt.X(1) .or. T(L).gt.X(N))
     |    call ERROR('Spline eval data out of range')
        I = 1
        J = N
        do while (J-I.gt.1)
          K = (I+J)/2
          if (X(K).gt.T(L)) then
            J = K
          else
            I = K
          endif
        enddo
        XP = T(L)-X(I)
        F(L) = Y(I) + C1(I)*XP + C2(I)*XP**2
      enddo
      END
C  ------------------------------------------------------------------- C
C  ----  Find entry in a forward sorted list using binary search  ---- C
C  ------------------------------------------------------------------- C
C  Search a set of N ordered (X,Y) values for X0 and return the corresponding Y0
C  The X values must be ordered.  The value of X0 must be in the X range.
C  Use binary search to find the X values adjacent to X0, and use linear
C  interpolation to evaluate Y0.
      FUNCTION BINFIND(X0,X,Y,N)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION X(N),Y(N)

C  Make sure X(1) <= X0 and X0 <= X(N)
      if (N.lt.1 .or. X0.lt.X(1) .or. X0.gt.X(N))
     |  call ERRORDBL('Value out of range in BINFIND',X0)
C  Find X(i) <= X0 and X0 <= X(j) using binary search
      I = 1
      J = N
      do while (J-I.gt.1)
        K = (I+J)/2
        if (X(K).eq.X0) then
          I = K
          J = K
        else if (X(K).gt.X0) then
          J = K
        else
          I = K
        endif
      enddo
C  Found an exact match
      if (I.eq.J) then
        BINFIND = Y(I)
C  Linearly interpolate to X0
      else
        F = (X0-X(I))/(X(J)-X(I))
        BINFIND = (1-F)*Y(I)+F*Y(J)
      ENDIF
      END

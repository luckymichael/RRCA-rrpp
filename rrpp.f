C  Republican River Pre-Processor (RRPP)
C  Republican River Settlement Model Version 12
C  Willem A. Schreuder
C  Last modified: Apr 26, 2012
C
C  Program to generate recharge and well pumping input files for MODFLOW
C  
C  Recharge from precipitation is calculated based on a curve.  Different curves
C  are used for irrigated and non-irrigated land, and for different soil types.
C  The recharge for irrigated and non-irrigated lands are calculated and then added
C  based on the amount of irrigated acreage in the cell.
C  
C  Recharge from applied surface water and groudwater is calculated as an
C  efficiency factor multiplied by the applied amount.  The efficiency factor for
C  groundwater depends on the system, whether flood irrigation or a center pivot
C  sprinkler system.  The efficiency of sprinkler systems also vary over time.
C  Efficiencies by year and by county were estimated and are used in the program
C  mkgw to generate files that are read by this program.  Recharge from surface
C  water and groundwater irrigation applications are read for each state on a cell
C  by cell basis.  For comingled lands, the surface water and groundwater portion
C  of the comingled land are read separately as part of the surface water and
C  groundwater lands.
C
C  All the above calculations are performed on a cell by cell basis.  The area and
C  volume of surface and groundwater irrigation is read for each year from cell by
C  cell files.  Annual precipitation at a set of weather stations is distributed to
C  individual cells by kriging.  A flag array classifies each cell according to one
C  soil type.  The recharge calculations are performed for each land type according
C  to the soil type and then composited to calculate recharge for the cell.  This
C  recharge rate is then adjusted by the terrain factor to yield the final recharge
C  number.  The resulting annual recharge is then distributed to months using a
C  static user supplied distribution.
C
C  Steady state recharge is calculated by averaging the recharge produced for a
C  user supplied range of years.  Since this reflects the pre-development
C  condition, the calculation is done using only the non-irrigated curves.
C
C  For impact runs, i.e. runs where pumping or canals are switched off, that
C  portion of the irrigated lands that relates to the source of water is switched
C  off in the recharge calculation, as is pumping when applicable.
C
C  Modified 9/18/2007 (WAS)
C  Changed addition of ACO so that if neither GW nor SW is applied (NOPUMP and MOUND),
C  the commingled acres are not counted.  Before we tacitly assumed that these options 
C  are mutually exclusive.
C  Correction 4/26/2012 (WAS)
C  Sam Perkins pointed out that the ACO logic was reversed, which
C  effects the mound calculations
      PROGRAM RRPP
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'rrpp.ins'
      DIMENSION RCH(NCOL,NROW),CURV(NCOL,NROW)
      DIMENSION FRCH(MAXRCH,2,0:MAXSOIL),ISOIL(NCOL,NROW)
      DIMENSION IBGT(NCOL,NROW),MOUND(NCOL,NROW),TRN(NCOL,NROW)
      DIMENSION XSTN(MAXSTN),YSTN(MAXSTN),PSTN(MAXSTN,MINYR:MAXYR)
      DIMENSION FACTM(0:12),BGT(0:MAXBGT)
      DIMENSION XFACTOR(NCOL,NROW)
      LOGICAL       GW(3),EOL,APPEND
      INTEGER       NDM(12)
      CHARACTER*4   CYEAR
      CHARACTER*7   ROOT,NAME
      CHARACTER*16  CMD,CBGT(MAXBGT),CSOIL(MAXSOIL),ACCESS
      CHARACTER*256 FILE,RCHFIL(4),WELFIL
      DATA          STATE /'co','ks','ne'/
      DATA          NDM /31,28,31,30,31,30,31,31,30,31,30,31/

C  Generate cell coordinates
      do j=1,NCOL
        X(j) = (j-0.5)*DX+X0
      enddo
      do i=1,NROW
        Y(i) = (NROW-i+0.5)*DY+Y0
      enddo

C  Set defaults
      BIN      = .FALSE.
      NRCH     = 0
      NOUT     = 0
      WELFIL   = ' '
      IYRSS0   = 0
      IYRSS1   = 0
      IYRTR0   = 0
      IYRTR1   = -1
      IYRXF    = 9999
      NSOIL    = 0
      NBGT     = 0
      NSTN     = 0
      NPTS     = 0
      IDR      = 0
      FACTM(0) = 1
      APPEND   = .FALSE.
      do M=1,12
        FACTM(M) = 0
      enddo
      do I=1,NROW
        do J=1,NCOL
          IBOUND(J,I) = 1
          IBGT(J,I)   = 0
          TRN(J,I)    = 1
        enddo
      enddo
C  Switches for impact runs
      do L=1,3
        GW(L) = .TRUE.
      enddo
      do I=1,NROW
        do J=1,NCOL
          MOUND(J,I) = 0
        enddo
      enddo

C======================================================================C
C======================  Read and process inputs  =====================C
C======================================================================C
      if (IARGC().ne.1) call ERROR('Usage: rrpp <parameter file>')
      call GETARG(1,FILE)
      OPEN(1,FILE=FILE,STATUS='OLD',IOSTAT=IOS)
      if (IOS.ne.0) call ERROR('Error opening parameter file')
100   READ(1,'(A)',END=199) LINE
      if (LINEFIX(LINE).ne.0) goto 100
      CMD = CHARIN(LINE,1,IPOS,IOS)
      if (IOS.ne.0) call ERROR('Invalid command '//LINE)
      call UPPERCASE(CMD)

C
C  File directory, extension and type for saving precip data
C    This stores a bunch of files containing intermediate recharge results
C    in the recharge calculation for display purposes.
C
      if (CMD.eq.'FILE') then
        FDIR = CHARIN(LINE,IPOS,IPOS,IOS)
        if (IOS.ne.0) call ERROR('Invalid FILE directory')
        do i=1,MAXEXT
          FEXT(i) = CHARIN(LINE,IPOS,IPOS,IOS)
          if (IOS.ne.0) call ERRORINT('Invalid FILE extention',i)
        enddo
        NOUT = MAXEXT
C
C  Input data file directories
C    Where state pumping and recharge files are read from
C
      else if (CMD.eq.'DIR') then
        do L=1,3
          STATE(L) = CHARIN(LINE,IPOS,IPOS,IOS)
          if (IOS.ne.0) call ERROR('Error reading DIR '//LINE)
        enddo
C
C  APPPEND output file(s)
C
      else if (CMD.eq.'APPEND') then
        APPEND = .TRUE.
C
C  Recharge output file(s)
C    MODFLOW recharge output file name
C    Can split recharge into PPT/GW/SW/Canals in 4 files
C
      else if (CMD.eq.'RECHARGE') then
        NRCH = 0
        do while (.not.EOL(LINE,IPOS))
          NRCH = NRCH+1
          RCHFIL(NRCH) = CHARIN(LINE,IPOS,IPOS,IOS)
          if (IOS.ne.0) call ERROR('Error reading RECHARGE file '//LINE)
        enddo
        if (NRCH.ne.1 .and. NRCH.ne.4)
     |    call ERROR('RECHARGE must specify 1 or 4 files '//LINE)
C
C  Well output file
C    MODFLOW well file
C
      elseif (CMD.eq.'WELL') then
        WELFIL = CHARIN(LINE,IPOS,IPOS,IOS)
        if (IOS.ne.0) call ERROR('Error reading WELL file name '//LINE)
C
C  Binary input files
C    If you want to use binary input files
C
      elseif (CMD.eq.'BIN') then
        BIN = .TRUE.
C
C  Years to use for steady state calculation
C    Specifies start and end year of period used to calculate recharge
C    SS recharge is the average recharge for this period multiplied
C    by the global multiplier specified as the third parameter
C
      else if (CMD.eq.'STEADY') then
        IYRSS0 = INTIN(LINE,IPOS,IPOS,IOS)
        if (IOS.ne.0) call ERROR('Invalid STEADY start year')
        IYRSS1 = INTIN(LINE,IPOS,IPOS,IOS)
        if (IOS.ne.0) call ERROR('Invalid STEADY end year')
        if (IYRSS0.gt.IYRSS1) call ERROR('Invalid STEADY period')
        if (IYRSS0.lt.MINYR) call ERROR('STEADY period start error')
        if (IYRSS1.gt.MAXYR) call ERROR('STEADY period end error')
        SSPPT = REALIN(LINE,IPOS,IPOS,IOS)
        if (IOS.ne.0) call ERROR('Invalid STEADY precip multiplier')
C
C  Years to use for transient calculation
C    Specifies start and end year of transient period
C    output is generated for calendar years in this range.
C
      else if (CMD.eq.'TRANSIENT') then
        IYRTR0 = INTIN(LINE,IPOS,IPOS,IOS)
        if (IOS.ne.0) call ERROR('Invalid TRANSIENT start year')
        IYRTR1 = INTIN(LINE,IPOS,IPOS,IOS)
        if (IOS.ne.0) call ERROR('Invalid TRANSIENT end year')
        if (IYRTR0.gt.IYRTR1) call ERROR('Invalid TRANSIENT period')
        if (IYRTR0.lt.MINYR) call ERROR('TRANSIENT period start error')
        if (IYRTR1.gt.MAXYR) call ERROR('TRANSIENT period end error')
C
C  Read groundwater pumping switch
C    This switch inhibits pumping for this state
C    Used to make impact runs
C
      elseif (CMD.eq.'NOPUMP') then
        L = INSTATE(IOS)
        if (IOS.ne.0) call ERROR('Invalid state '//LINE)
        GW(L) = .FALSE.
C
C  Read mound adjustment
C    This switch inhibits surface water diversions and canal leakage in
C    the mound area.  Takes a flag file designating where the mound is
C    A 0 means outside the mound, 1 means in the mound.
C
      elseif (CMD.eq.'MOUND') then
        FILE = CHARIN(LINE,IPOS,IPOS,IOS)
        if (IOS.ne.0) call ERROR('Error reading MOUND flag name '//LINE)
        call READFLAGS(FILE,MOUND,NROW,NCOL,CMD,N,1)
C
C  Load IBOUNDs from file
C    Designates groundwater model domain
C
      elseif (CMD.eq.'IBOUND') then
        FILE = CHARIN(LINE,IPOS,IPOS,IOS)
        if (IOS.ne.0) call ERROR('Invalid IBOUND file name '//LINE)
        call READFLAGS(FILE,IBOUND,NROW,NCOL,LINE,N,0)
C
C  Load SOIL flags from file
C    Flag array to assign one soil type to each cell
C
      elseif (CMD.eq.'SOIL') then
        FILE = CHARIN(LINE,IPOS,IPOS,IOS)
        if (IOS.ne.0) call ERROR('Invalid SOIL file name '//LINE)
        call READFLAGS(FILE,ISOIL,NROW,NCOL,CSOIL,NSOIL,MAXSOIL)
C
C  Read terrain ZONEs
C    (see subroutine for explanation)
C
      elseif (CMD.eq.'TERRAIN') then
        call READTERRAIN(TRN)
C
C  Read XFACTOR year and multipliers
C    (see subroutine for explanation)
C
      elseif (CMD.eq.'XFACTOR') then
        IYRXF = INTIN(LINE,IPOS,IPOS,IOS)
        if (IOS.ne.0) call ERROR('Invalid XFACTOR year '//LINE)
        FILE = CHARIN(LINE,IPOS,IPOS,IOS)
        if (IOS.ne.0) call ERROR('Invalid XFACTOR file name '//LINE)
        call READREAL(FILE,XFACTOR,NROW,NCOL,0)
C
C  Read budget ZONEs
C    Quick and dirty test - prints budget by these zones to screen
C
      elseif (CMD.eq.'BUDGET') then
        FILE = CHARIN(LINE,IPOS,IPOS,IOS)
        if (IOS.ne.0) call ERROR('Invalid BUDGET file name '//LINE)
        call READFLAGS(FILE,IBGT,NROW,NCOL,CBGT,NBGT,MAXBGT)
        print '(A,99A10)','Time',
     |   (CBGT(L)(1:LENGTH(CBGT(L))),L=1,NBGT),'Total'
C
C  Read precip stations locations, data and drift
C    (see subroutine for explanation)
C
      elseif (CMD.eq.'STATIONS') then
        call READSTATION(XSTN,YSTN,NSTN,PSTN,IDR)
C
C  Read factors for annual -> month recharge conversion
C    Twelve factors (Jan-Dec) used to distribute annual recharge calcs
C    to months (the model uses monthly stress periods)
C
      elseif (CMD.eq.'MONTH') then
        do M=1,12
          FACTM(M) = REALIN(LINE,IPOS,IPOS,IOS)
          if (IOS.ne.0 .or. FACTM(M).LT.0)
     |      call ERROR('Invalid MONTH multiplier '//LINE)
        enddo
C
C  Read recharge function
C    (see subroutine for explanation)
C
      elseif (CMD.eq.'FUNCTION') then
        call READFUNCTION(FRCH,NPTS,NSOIL)
C  Duh!
      else
        call ERROR('Unknown parameter '//CMD)
      endif
      goto 100
199   CLOSE(1)

C======================================================================C
C==================  Calculate Pumping and Recharge  ==================C
C======================================================================C
      if (NSOIL.eq.0) call ERROR('No SOILs specified')
      if (NPTS.eq.0)  call ERROR('No curve specified')

C  Check month multipliers
      SUM = 0
      do M=1,12
        SUM = SUM+FACTM(M)
      enddo
      if (ABS(SUM-1).gt.1e-10)
     |  call ERROR('MONTH multipliers do not sum to 1')

C  Check if terrain kills cells or missing soil
      do I=1,NROW
        do J=1,NCOL
          if (IBOUND(J,I).ne.0) then
             if (TRN(J,I).le.0) print *,'WARNING: ',
     |         'Terrain multiplier error ',I,J,TRN(J,I)
             if (ISOIL(J,I).eq.0) print *,'WARNING: ',
     |         'Missing soil type ',I,J
          endif
        enddo
      enddo

C  Set up precipitation kriging matrix and invert it
      call KRIGMATINV(X0,Y0,XSTN,YSTN,AMAT,AINV,NSTN,IDR)

C
C  Open recharge and pumping output files
C
      if (APPEND) then
        ACCESS = 'APPEND'
      else
        ACCESS = 'SEQUENTIAL'
      endif
      MWEL = 0
      if (WELFIL.eq.' ') call ERROR('No well file specified')
      OPEN(LUOUT,FILE=WELFIL,STATUS='UNKNOWN',ACCESS=ACCESS,IOSTAT=IOS)
      if (IOS.ne.0) call ERROR('Cannot open file '//WELFIL)
      if (.not.APPEND) WRITE(LUOUT,'(2I10,A)') MAXWEL,0,' NOPRINT'
      if (NRCH.eq.0) call ERROR('No recharge file(s) specified')
      do L=1,NRCH
        OPEN(LUOUT+L,FILE=RCHFIL(L),STATUS='UNKNOWN',
     |    ACCESS=ACCESS,IOSTAT=IOS)
        if (IOS.ne.0) call ERROR('Cannot open file '//RCHFIL(L))
        if (.not.APPEND) WRITE(LUOUT+L,'(2I10)') 1,0
      enddo
C
C  Calculate steady state recharge
C
      if (IYRSS0.gt.0) then
        do I=1,NROW
          do J=1,NCOL
            RCH(J,I) = 0
          enddo
        enddo
C       Accumulate recharge for steady state years
C       The SS switch inhibits irrigation for this calculation
        do IYR=IYRSS0,IYRSS1
          call EVALCURVE(CURV,.TRUE.,IYR,GW,MOUND,TRN,ISOIL,NSOIL,
     |      FRCH,NPTS,XSTN,YSTN,NSTN,PSTN,IDR)
          do I=1,NROW
            do J=1,NCOL
              RCH(J,I) = RCH(J,I) + CURV(J,I)
            enddo
          enddo
        enddo
C       Average recharge
        do I=1,NROW
          do J=1,NCOL
            RCH(J,I) = RCH(J,I)/(IYRSS1-IYRSS0+1)
          enddo
        enddo
        call BUDGET('SS  ',RCH,IBGT,BGT,NBGT,NROW,NCOL)
C       Save recharge for annual average period
        NSEC = 365.25*24*60*60
        call PERIOD(RCH,GW,MOUND,'Steady ','steady','steady',
     |              -1,IYRXF,XFACTOR,NSEC,SSPPT)
      endif

C
C  Calculate transient recharge
C
      do IYR=IYRTR0,IYRTR1
C       Evaluate recharge from curve for this year
        call EVALCURVE(CURV,.FALSE.,IYR,GW,MOUND,TRN,ISOIL,NSOIL,
     |    FRCH,NPTS,XSTN,YSTN,NSTN,PSTN,IDR)
        WRITE(CYEAR,'(I4)') IYR
        call BUDGET(CYEAR,CURV,IBGT,BGT,NBGT,NROW,NCOL)
C       Calculate monthly recharge from annual curve and monthly stresses
        do MO=1,12
          if (MOD(IYR,4).eq.0 .and. MO.eq.2) then
            NDAY = 29
          else
            NDAY = NDM(MO)
          endif
          NSEC = NDAY*24*60*60
          WRITE(ROOT,'(F7.2)') IYR+MO/100D0
          WRITE(NAME,'(I2,1H/,I4)') MO,IYR
          call PERIOD(CURV,GW,MOUND,NAME,ROOT,CYEAR,
     |                IYR,IYRXF,XFACTOR,NSEC,FACTM(MO))
        enddo
      enddo
      END
C  ------------------------------------------------------------------- C
C  ------------  Calculate and print state zone budgets  ------------- C
C  ------------------------------------------------------------------- C
      SUBROUTINE BUDGET(DATE,RCH,IBGT,BGT,NBGT,NROW,NCOL)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION RCH(NCOL,NROW),IBGT(NCOL,NROW),BGT(0:NBGT)
      CHARACTER*(*) DATE

      if (NBGT.eq.0) RETURN
C  Calculate recharge by state zone
      do L=0,NBGT
        BGT(L) = 0
      enddo
      do I=1,NROW
        do J=1,NCOL
          BGT(IBGT(J,I)) = BGT(IBGT(J,I)) + RCH(J,I)
        enddo
      enddo
C  Convert inches/cell to acre-feet
      SUM = 0
      do L=1,NBGT
        BGT(L) = BGT(L)/12*640
        SUM = SUM + BGT(L)
      enddo
      print '(A4,99F10.1)',DATE,(BGT(L),L=1,NBGT),SUM
      END
C  ------------------------------------------------------------------- C
C  ---------------------  Read TERRAIN command  ---------------------- C
C  ------------------------------------------------------------------- C
C  Read TERRAIN command.
C  Return a terrain multiplier array TRN
C  An flag array of terrain multipliers are read (up to 9 terrain zones)
C  A list of terrain multipliers are read (X,Y,M1,...,Mn), n=# terrain zones
C  The LIMIT qualifier can be used to limit calculated terrain multipliers.
C
C  For every cell, the terrain multipliers for the zone type corresponding
C  to that cell is kriged and a value for the cell is obtained.  When the
C  zone flag is 0, a terrain multiplier of 1.0 is used.
      SUBROUTINE READTERRAIN(TRN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'rrpp.ins'
      PARAMETER     (MAXTRN=9)
      DIMENSION     XTRN(MAXSTN),YTRN(MAXSTN),FTRN(MAXSTN,MAXTRN)
      DIMENSION     TRN(NCOL,NROW),ITRN(NCOL,NROW)
      CHARACTER*16  CMD,CTRN(MAXTRN)
      CHARACTER*256 FILE

C     Read terrain flags
      FILE = CHARIN(LINE,IPOS,IPOS,IOS)
      if (IOS.ne.0) call ERROR('Invalid TERRAIN flag file '//LINE)
      call READFLAGS(FILE,ITRN,NROW,NCOL,CTRN,NTRN,MAXTRN)
C     Read terrain multipliers
      FILE = CHARIN(LINE,IPOS,IPOS,IOS)
      if(IOS.ne.0)call ERROR('Invalid TERRAIN multiplier file '//LINE)
      OPEN(8,FILE=FILE,STATUS='OLD',IOSTAT=IOS)
      if (IOS.ne.0) call ERROR('Error opening file '//FILE)
C     Get terrain plot file name
      FILE = CHARIN(LINE,IPOS,IPOS,IOS)
      if(IOS.ne.0)call ERROR('Invalid TERRAIN plot file '//LINE)
C     Read optional parameters
      TMIN = 0
      TMAX = 1D100
      CMD = CHARIN(LINE,IPOS,IPOS,IOS)
      do while (IOS.eq.0)
        if (CMD.eq.'LIMIT') then
          TMIN = REALIN(LINE,IPOS,IPOS,IOS)
          if (IOS.ne.0) call ERROR('Invalid TERRAIN LIMIT minimum')
          TMAX = REALIN(LINE,IPOS,IPOS,IOS)
          if (IOS.ne.0) call ERROR('Invalid TERRAIN LIMIT maximum')
          if (TMIN.gt.TMAX)
     |      call ERROR('Invalid TERRAIN LIMIT range '//LINE)
        else
          call ERROR('Invalid TERRAIN parameter '//CMD)
        endif
        CMD = CHARIN(LINE,IPOS,IPOS,IOS)
      enddo
C     Read locations and multipliers
      READ(8,*)
      NMUL = 0
      READ(8,'(A)',IOSTAT=IOS) LINE
      do while (IOS.eq.0 .and. LINE.ne.'END')
        if (LINEFIX(LINE).eq.0) then
          NMUL = NMUL+1
          if (NMUL.gt.MAXSTN) call ERROR('Too many TERRAIN locations')
          XTRN(NMUL) = REALIN(LINE,1,IPOS,IOS)
          if (IOS.ne.0) call ERROR('Invalid terrain X '//LINE)
          YTRN(NMUL) = REALIN(LINE,IPOS,IPOS,IOS)
          if (IOS.ne.0) call ERROR('Invalid terrain Y '//LINE)
          do L=1,NTRN
            FTRN(NMUL,L) = REALIN(LINE,IPOS,IPOS,IOS)
            if (IOS.ne.0)
     |        call ERROR('Invalid terrain multiplier '//LINE)
          enddo
        endif
        READ(8,'(A)',IOSTAT=IOS) LINE
      enddo
      CLOSE(8)
C     Set up kriging matrix and invert it
      call KRIGMATINV(X0,Y0,XTRN,YTRN,AMAT,AINV,NMUL,0)
C     Krige multipliers to cells
      do I=1,NROW
        do J=1,NCOL
          L = ITRN(J,I)
          if (L.eq.0) then
            TRN(J,I) = 1
          else
            TRN(J,I) = EVALKRIGE(AINV,X0,Y0,XTRN,YTRN,FTRN(1,L),
     |                             X(J),Y(I),NMUL,0)
            if (TMIN.gt.TRN(J,I)) TRN(J,I) = TMIN
            if (TMAX.lt.TRN(J,I)) TRN(J,I) = TMAX
          endif
        enddo
      enddo
C     Save terrain multipliers
      OPEN(9,FILE=FILE,STATUS='UNKNOWN',IOSTAT=IOS)
      if (IOS.ne.0) call ERROR('Error opening file '//FILE)
      do I=1,NROW
        WRITE(9,'(999F6.3)',IOSTAT=IOS) (TRN(J,I),J=1,NCOL)
        if (IOS.ne.0) call ERROR('Cannot write file '//FILE)
      enddo
      CLOSE(9)
      END
C  ------------------------------------------------------------------- C
C  ---------------------  Read precip STATIONs  ---------------------- C
C  ------------------------------------------------------------------- C
      SUBROUTINE READSTATION(XSTN,YSTN,NSTN,PSTN,IDR)
C  Read precip station data.
C  Return an array of station locations and one value per year.
C  Also return the kriging drift order (this should be one).
C
C  Station names and locations are read from the location file in
C  free format with NAME, X and Y appearing one per line. Skip the
C  first line as it is considered to be a comment.  Coordinates are
C  model coordinates.
C
C  Precipitation data is read from a free format file.
C  Line one must contain the word YEAR and the station names in the
C  same order as om the locations file
C  Subsequent lines must contain the year and one value per station in
C  inches.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'rrpp.ins'
      CHARACTER*16  CSTN(MAXSTN)
      CHARACTER*256 FILE,STRING
      DIMENSION XSTN(MAXSTN),YSTN(MAXSTN),PSTN(MAXSTN,MINYR:MAXYR)

C     Process input line
      FILE = CHARIN(LINE,IPOS,IPOS,IOS)
      if (IOS.ne.0)
     |  call ERROR('Invalid STATION location file name '//LINE)
      OPEN(8,FILE=FILE,STATUS='OLD',IOSTAT=IOS)
      if (IOS.ne.0) call ERROR('Error opening file '//FILE)
      FILE = CHARIN(LINE,IPOS,IPOS,IOS)
      if (IOS.ne.0)
     |  call ERROR('Invalid STATION precipitation file name '//LINE)
      IDR = INTIN(LINE,IPOS,IPOS,IOS)
      if (IOS.ne.0 .or. IDR.lt.0 .or. IDR.gt.MAXDR)
     |  call ERROR('Invalid STATION drift order '//LINE)
C     Read station names and locations
      READ(8,*)
      NSTN = 0
      READ(8,'(A)',IOSTAT=IOS) LINE
      do while (IOS.eq.0 .and. LINE.ne.'END')
        if (LINEFIX(LINE).eq.0) then
          NSTN = NSTN+1
          if (NSTN.gt.MAXSTN) call ERROR('Too many stations')
          CSTN(NSTN) = CHARIN(LINE,1,IPOS,IOS)
          if (IOS.ne.0) call ERROR('Invalid station name '//LINE)
          XSTN(NSTN) = REALIN(LINE,IPOS,IPOS,IOS)
          if (IOS.ne.0) call ERROR('Invalid station X '//LINE)
          YSTN(NSTN) = REALIN(LINE,IPOS,IPOS,IOS)
          if (IOS.ne.0) call ERROR('Invalid station Y '//LINE)
        endif
        READ(8,'(A)',IOSTAT=IOS) LINE
      enddo
      CLOSE(8)
C     Read precipitation data station names
      OPEN(8,FILE=FILE,STATUS='OLD',IOSTAT=IOS)
      if (IOS.ne.0) call ERROR('Error opening file '//FILE)
      READ(8,'(A)') LINE
      STRING = CHARIN(LINE,1,IPOS,IOS)
      do L=1,NSTN
        STRING = CHARIN(LINE,IPOS,IPOS,IOS)
        if (STRING.ne.CSTN(L))
     |    call ERROR('Precip data out of order '//STRING)
      enddo
C     Initialize precip array to zero
      do I=1,NSTN
        do J=MINYR,MAXYR
          PSTN(I,J) = 0
        enddo
      enddo
C     Read precipitation data - discard data out of range
      READ(8,'(A)',IOSTAT=IOS) LINE
      do while (IOS.eq.0 .and. LINE.ne.'END')
        if (LINEFIX(LINE).eq.0) then
          IYR = INTIN(LINE,1,IPOS,IOS)
          if (IOS.ne.0) call ERROR('Invalid precip year '//LINE)
          if (IYR.gt.MINYR .and. IYR.lt.MAXYR) then
            do L=1,NSTN
              PSTN(L,IYR) = REALIN(LINE,IPOS,IPOS,IOS)
              if (IOS.ne.0) call ERROR('Invalid precip data '//LINE)
            enddo
          endif
        endif
        READ(8,'(A)',IOSTAT=IOS) LINE
      enddo
      CLOSE(8)
      END
C  ------------------------------------------------------------------- C
C  ---------------------  Read precip FUNCTION  ---------------------- C
C  ------------------------------------------------------------------- C
C  Read the FUNCTION command and calculate the curve values
C  Returns a curve in FRCH for every soil with values every 0.1 inches
C  A set of two recharges (non-irrigated and irrigated) is read for every
C  soil type in that order.
C  The curve is set up by interpolating between the user supplied knots
C  with a linear or quadratic spline (NSPL=1 or NSPL=2)
      SUBROUTINE READFUNCTION(FRCH,NPTS,NSOIL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'rrpp.ins'
      DIMENSION     FRCH(MAXRCH,2,0:MAXSOIL)
      DIMENSION     PPT(MAXFUN),RCH(MAXFUN,2,MAXSOIL)
      CHARACTER*256 FILE

      if (NSOIL.eq.0) call ERROR('No SOILs specified before FUNCTION')
      FILE = CHARIN(LINE,IPOS,IPOS,IOS)
      if (IOS.ne.0)
     |  call ERROR('Invalid FUNCTION output file name '//LINE)
      NSPL = INTIN(LINE,IPOS,IPOS,IOS)
      if (IOS.ne.0 .or. NSPL.lt.1 .or. NSPL.gt.2)
     |  call ERROR('Invalid FUNCTION spline type '//LINE)
      NFUN = 0
100   READ(1,'(A)',IOSTAT=IOS) LINE
      if (IOS.ne.0) goto 199
      if (LINEFIX(LINE).ne.0) goto 100
      VAL = REALIN(LINE,1,IPOS,IOS)
      if (IOS.ne.0) goto 198
      NFUN = NFUN+1
      if (NFUN.gt.MAXFUN) call ERROR('Too many FUNCTION data')
      PPT(NFUN) = VAL
      if (NFUN.eq.1) then
        if (VAL.ne.0)
     |    call ERROR('FUNCTION data must start at zero')
      else if (PPT(NFUN-1).ge.VAL) then
        call ERROR('FUNCTION data not monotone increasing')
      endif
      do L=1,NSOIL
        do K=1,2
          RCH(NFUN,K,L) = REALIN(LINE,IPOS,IPOS,IOS)
          if (IOS.ne.0) call ERROR('Invalid FUNCTION data '//LINE)
        enddo
      enddo
      goto 100
198   BACKSPACE(1)
199   CONTINUE
C     Evaluate spines at 1/10 inch increments
      NPTS = 10*PPT(NFUN)+1
      if (NPTS.gt.MAXRCH) call ERROR('FUNCTION range too large')
      do L=1,NPTS
        FRCH(L,1,0) = (L-1)/10D0
        FRCH(L,2,0) = (L-1)/10D0
      enddo
C     Fill in gaps with linear or quadratic splines
      do L=1,NSOIL
        do K=1,2
          call SPLINE(PPT,RCH(1,K,L),NFUN,FRCH,FRCH(1,K,L),NPTS,NSPL)
        enddo
      enddo
C     Save recharge curves to file
      OPEN(9,FILE=FILE,STATUS='UNKNOWN',IOSTAT=IOS)
      if (IOS.ne.0) call ERROR('Error opening file '//FILE)
      WRITE(9,'(A6,999(1X,A,I3))')
     |  'Precip',(('Soil',L,K=1,2),L=1,NSOIL)
      do I=1,NPTS
        WRITE(9,'(F6.1,999F8.3)')
     |    FRCH(I,1,0),((FRCH(I,K,L),K=1,2),L=1,NSOIL)
      enddo
      CLOSE(9)
      END
C----------------------------------------------------------------------C
C---------------------  Read state code from line  --------------------C
C----------------------------------------------------------------------C
      FUNCTION INSTATE(IOS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'rrpp.ins'
      CHARACTER*8 STRING

      STRING = CHARIN(LINE,IPOS,IPOS,IOS)
      call UPPERCASE(STRING)
      if (IOS.ne.0) then
        INSTATE = 0
      else if (STRING.eq.'CO ') then
        INSTATE = CO
      else if (STRING.eq.'KS ') then
        INSTATE = KS
      else if (STRING.eq.'NE ') then
        INSTATE = NE
      else
        INSTATE = 0
        IOS = -1
      endif
      END
C----------------------------------------------------------------------C
C-------------  Calculate recharge from curve for year  ---------------C
C----------------------------------------------------------------------C
      SUBROUTINE EVALCURVE(CURV,SS,IYR,GW,MOUND,TRN,ISOIL,NSOIL,
     |  FRCH,NPTS,XSTN,YSTN,NSTN,PSTN,IDR)
C  This subroutine calculates annual precipitation recharge (inches/year) in CURV
C  Recharge is calculated based on supplied water.
C  To calculate irrigated area, the area of irrigation for each state on 
C  a cell by cell basis is read from the area files stored by mkgw, mksw and mknedat.
C  When SS is true, this is the predevelopment steady state so ignore irrigation
C  Honor GW and MOUND to selectively disable irrigation for impact runs
C  Whatever area is left over in a cell is non-irrigated.
C  If the irrigated area is greater than the cell area, assume the whole cell is irrigated
C  Precipitation for each cell is obtained by kriging met stations data
C  The curve yields a recharge rate in inches/year.
C  Recharge rate is weighted by area fraction of irrigated and non-irrigated in the cell.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'rrpp.ins'
      PARAMETER (CELL=640)
      DIMENSION CURV(NCOL,NROW),RPP(NCOL,NROW),WORK(NCOL,NROW)
      DIMENSION AREA(NCOL,NROW),FIRR(NCOL,NROW)
      DIMENSION FRCH(MAXRCH,2,0:MAXSOIL),MOUND(NCOL,NROW)
      DIMENSION TRN(NCOL,NROW),ISOIL(NCOL,NROW)
      DIMENSION XSTN(MAXSTN),YSTN(MAXSTN),PSTN(MAXSTN,MINYR:MAXYR)
      LOGICAL   SS,GW(3)
      CHARACTER*4 CYEAR

C  Set up string as year for use in file names
      WRITE(CYEAR,'(I4)') IYR

C  Initialize well pumping and recharge arrays
      do I=1,NROW
        do J=1,NCOL
          AREA(J,I) = 0
          CURV(J,I) = 0
          RPP(J,I)  = 0
        enddo
      enddo

C  For steady state there is by definition no irrigated area
C  For impact runs, adjust area for GW and MOUND
      if (.not.SS) then
C  Read groundwater irrigated area
        do L=1,3
          if (GW(L)) then
            call READFILE(WORK,STATE(L),CYEAR,'agw',.FALSE.)
            call ADD(AREA,1D0,WORK)
          endif
        enddo
C  Read surface water irrigated area
        do L=1,3
          call READFILE(WORK,STATE(L),CYEAR,'asw',.FALSE.)
          call CONDADD(AREA,1D0,WORK,MOUND)
        enddo
C  Read comingled irrigated area
        do L=1,3
          call READFILE(WORK,STATE(L),CYEAR,'aco',.FALSE.)
C    With no pumping, add only outside the mound
C    This only makes a difference if no pumping and no mound is simulated
          if (.not.GW(L)) then
            call CONDADD(AREA,1D0,WORK,MOUND)
C    With pumping, add everywhere regardless of mound
          else
            call ADD(AREA,1D0,WORK)
          endif
        enddo
      endif

C  Calculate recharge from curves
      do I=1,NROW
        do J=1,NCOL
          if (IBOUND(J,I).gt.0) then
            L = ISOIL(J,I)
C           Compute fraction of cell that is irrigated
            FIRR(J,I) = MIN(1,AREA(J,I)/CELL)
C           Calculate precipitation in this cell (inches)
            RPP(J,I) = EVALKRIGE(AINV,X0,Y0,XSTN,YSTN,PSTN(1,IYR),
     |                           X(J),Y(I),NSTN,IDR)
            if (RPP(J,I).lt.0) call ERRORDBL('Negative PPT',RPP(J,I))
C           Compute composite precipitation recharge from curve for this cell
            C = (1-FIRR(J,I))*BINFIND(RPP(J,I),FRCH,FRCH(1,1,L),NPTS) +
     |              FIRR(J,I)*BINFIND(RPP(J,I),FRCH,FRCH(1,2,L),NPTS)
C           Adjust for terrain
            CURV(J,I) = TRN(J,I)*C
          endif
        enddo
      enddo
C  Save the results for plotting if requested and not steady state
      if (.not.SS .and. NOUT.gt.0) then
        call STOREMONTH(RPP,1D0,NROW,NCOL,IYR,-1,FDIR,FEXT(1),BIN)
        call STOREMONTH(FIRR,640D0,NROW,NCOL,IYR,-1,FDIR,FEXT(2),BIN)
        call STOREMONTH(CURV,1D0,NROW,NCOL,IYR,-1,FDIR,FEXT(3),BIN)
      endif
      END
C----------------------------------------------------------------------C
C--------  Calculate recharge and pumping for stress period  ----------C
C----------------------------------------------------------------------C
      SUBROUTINE PERIOD(CURV,GW,MOUND,NAME,ROOT,CYEAR,
     |                  IYR,IYRXF,XFACTOR,NSEC,PPTADJ)
C  This subroutine stores the pumping and recharge to the MODFLOW files.
C  This subroutine is called once for every stress period (months).
C  Pumping is read from the ".mi" cell by cell files and ".pmp" files
C  prepared by mkgw and mknedat.  Skip pumping by state for impact runs.
C  After adding all the pumping, write it to the MODFLOW file after 
C  converting pumping from af to cfs.
C  Recharge is calculated from the annual recharge value previously calculated 
C  in EVALCURVE distributed to the month using the MONTH factors.
C  For the Eff+Curve method, the monthly return flows is also read from the
C  surface and groundwater cell by cell return flow files  calculated by
C  mkgw, mksw and mknedat and added to the precip recharge.
C  Finally recharge from canal leakage is added to recharge for both methods,  
C  qualified as needed by the mound.
C  Recharge is then saved in U2DREL format with the elements of the array in
C  inches and a constant multiplier to make the units ft/sec. 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE       'rrpp.ins'
      DIMENSION     PUMP(NCOL,NROW),RECH(NCOL,NROW,4)
      DIMENSION     CURV(NCOL,NROW),WORK(NCOL,NROW),MOUND(NCOL,NROW)
      DIMENSION     XFACTOR(NCOL,NROW)
      LOGICAL       GW(3)
      CHARACTER*(*) NAME,ROOT,CYEAR

C  Fraction of the year
      FRAC = NSEC/(365.25D0*24*60*60)

C  Initialize well pumping and recharge arrays
      do I=1,NROW
        do J=1,NCOL
          PUMP(J,I)   = 0
          RECH(J,I,1) = 0
          RECH(J,I,2) = 0
          RECH(J,I,3) = 0
          RECH(J,I,4) = 0
        enddo
      enddo

C  Recharge from curve (convert inches to AF)
      call ADD(RECH(1,1,1),CIN2AF*PPTADJ,CURV)

C  Read arrays by state
      do L=1,3
C       Well pumping
C         Adjust for impact runs
        if (GW(L)) then
C         M&I (annuals)
          call READFILE(WORK,STATE(L),CYEAR,'mi',.FALSE.)
          call ADD(PUMP,FRAC,WORK)
C         Irrigation Wells
          call READFILE(WORK,STATE(L),ROOT,'pmp',.FALSE.)
          call ADD(PUMP,1D0,WORK)
        endif
C       Groundwater Irrigation Returns
C         None without pumping
        if (GW(L)) then
          call READFILE(WORK,STATE(L),ROOT,'rcg',.FALSE.)
          call ADD(RECH(1,1,2),1D0,WORK)
        endif
C       Surface Water Irrigation Returns
C         Adjust for mound
        call READFILE(WORK,STATE(L),ROOT,'rcs',.FALSE.)
        call CONDADD(RECH(1,1,3),1D0,WORK,MOUND)
C       Canal Leakage
C         Adjust for mound
        call READFILE(WORK,STATE(L),ROOT,'rcc',.FALSE.)
        call CONDADD(RECH(1,1,4),1D0,WORK,MOUND)
      enddo

C  X Factor adjustment to irrigation returns
C    Only when year > initial X factor year
      if (IYR.gt.IYRXF) then
        do I=1,NROW
          do J=1,NCOL
            RECH(J,I,2) = XFACTOR(J,I)*RECH(J,I,2)
            RECH(J,I,3) = XFACTOR(J,I)*RECH(J,I,3)
          enddo
        enddo
      endif

C  Converts inches (over the period) to feet/sec
      CONST  = 1/12D0/NSEC
C  Acre-feet over 640 acres (over the period) to
C  Inches per unit area (over the period)
      AF2IN  = 12.0/640
C  Acre-feet (over the period) to cfs
      AF2CFS = 43560D0/NSEC
C  Write recharge array(s) in inches
      if (NRCH.eq.1) then
        WRITE(LUOUT+1,'(2I10)') 1,0
        WRITE(LUOUT+1,'(A,1P,E15.8,0P,A,I3,1X,A)')
     |    'INTERNAL ',CONST,' (20F9.5) ',-1,NAME
        do I=1,NROW
          WRITE(LUOUT+1,'(20F9.5)') (AF2IN*(RECH(J,I,1)+RECH(J,I,2)+
     |      RECH(J,I,3)+RECH(J,I,4)),J=1,NCOL)
        enddo
      else
        do L=1,NRCH
          WRITE(LUOUT+L,'(2I10)') 1,0
          WRITE(LUOUT+L,'(A,1P,E15.8,0P,A,I3,1X,A)')
     |      'INTERNAL ',CONST,' (20F9.5) ',-1,NAME
          do I=1,NROW
            WRITE(LUOUT+L,'(20F9.5)') (AF2IN*RECH(J,I,L),J=1,NCOL)
          enddo
        enddo
      endif
C  Count pumping wells
      NWEL = 0
      do I=1,NROW
        do J=1,NCOL
          if (PUMP(J,I).ne.0) NWEL = NWEL+1
        enddo
      enddo
      if (NWEL.gt.MAXWEL) then
        print *,'Increase MAXWEL to ',NWEL
        STOP
      endif
      MWEL = MAX(MWEL,NWEL)
C  Write pumping wells
      WRITE(LUOUT,'(I5,I3,1X,A)') NWEL,-1,NAME
      do I=1,NROW
        do J=1,NCOL
          if (PUMP(J,I).ne.0)
     |      WRITE(LUOUT,'(I1,2I4,F11.6)') 1,I,J,-AF2CFS*PUMP(J,I)
        enddo
      enddo
      END
C----------------------------------------------------------------------C
C--------------------  Read data file into array  ---------------------C
C----------------------------------------------------------------------C
      SUBROUTINE READFILE(WORK,DIR,ROOT,EXT,REQUIRED)
C  Read a file of cell by cell values into a work array.
C  If the file is required generate an error if it does not exist.
C  If the file does not exist set array to zero.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE       'rrpp.ins'
      LOGICAL       REQUIRED,EXIST
      DIMENSION     WORK(NCOL,NROW)
      CHARACTER*256 FILE
      CHARACTER*(*) DIR,ROOT,EXT

      WRITE(FILE,'(A,1H/,A,1H.,A)') DIR(1:LENGTH(DIR)),ROOT,EXT
      INQUIRE(FILE=FILE,EXIST=EXIST)
      if (.not.EXIST) then
        if (REQUIRED) call ERROR('File does not exist '//FILE)
        do I=1,NROW
          do J=1,NCOL
            WORK(J,I) = 0
          enddo
        enddo
      else if (BIN) then
        OPEN(8,FILE=FILE,FORM='UNFORMATTED',STATUS='OLD',IOSTAT=IOS)
        if (IOS.ne.0) call ERROR('Cannot open file '//FILE)
        READ(8,IOSTAT=IOS) WORK
        if (IOS.ne.0) call ERROR('Cannot read file '//FILE)
        CLOSE(8)
      else
        call READREAL(FILE,WORK,NROW,NCOL,0)
      endif
      END
C----------------------------------------------------------------------C
C--------------------------  Array A = A+gB  --------------------------C
C----------------------------------------------------------------------C
      SUBROUTINE ADD(A,g,B)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE   'rrpp.ins'
      DIMENSION A(NCOL,NROW),B(NCOL,NROW)

      do I=1,NROW
        do J=1,NCOL
          A(J,I) = A(J,I) + g*B(J,I)
        enddo
      enddo
      END
C----------------------------------------------------------------------C
C--------------------------  Array A = A+gB  --------------------------C
C----------------------------------------------------------------------C
      SUBROUTINE CONDADD(A,g,B,M)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE   'rrpp.ins'
      DIMENSION A(NCOL,NROW),B(NCOL,NROW),M(NCOL,NROW)

      do I=1,NROW
        do J=1,NCOL
          if (M(J,I).eq.0) A(J,I) = A(J,I) + g*B(J,I)
        enddo
      enddo
      END

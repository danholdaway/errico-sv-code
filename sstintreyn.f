      SUBROUTINE JULIAN( MDATE, nyrp, nwkp, WT )
 
      INTEGER LENMON(12), midwk1(81:97), JPREV(12)
 
      DATA LENMON / 31, 28, 31, 30, 31, 30
     A            , 31, 31, 30, 31, 30, 31 /
c Julian day of center date of first week of data for each year
c first line;  COADS  data 1981-1989, Sunday    centered
c second line; NESDIS data 1990-1997, Wednesday centered
c                   D0  D1 D2 D3 D4 D5 D6 D7 D8 D9
      DATA midwk1 /    305, 3, 2, 1, 6, 5, 4, 3, 1
     A             , 3,  2, 1, 6, 5, 4, 3, 1 /
 
C ******           INITIALIZE nwkp, nyrp
      nwkp = 1
      nyrp = -1
 
      IDATE = MDATE / 100
      IYR = IDATE / 10000
      IMO = ( IDATE - IYR*10000 ) / 100
      IDAY = MOD( IDATE, 100 )

c ****** current (5/98) CAS Reynolds weekly (DS277.0)
c        sst dataset is thru dec. 97
      if (iyr.gt.97) then
        print *, '*** SST data only thru Dec 1997.'
        stop 1997
      endif
 
      ILEAP = MOD( IYR, 4 )
      IF (IYR.EQ.0) ILEAP=1
      LENMON(2) = 28
      IF(ILEAP.EQ.0) LENMON(2) = 29
 
      JPREV(1) = 0
      DO 10 J=2,12
         JPREV(J)  = JPREV(J-1) + LENMON(J-1)
   10 CONTINUE
      JULDAY = IDAY + JPREV(IMO)
 
      PRINT *, 'MDATE, IYR, IMO, IDAY, JULDAY = '
     A       ,  MDATE, IYR, IMO, IDAY, JULDAY
 

      iwk=(JULDAY-midwk1(IYR)+7)/7
      ILEAP = 0 
      IF (iwk.EQ.0) THEN
         nyrp=IYR-1
         IF (MOD( nyrp, 4 ).EQ.0 .AND. nyrp.NE.0) ILEAP=1
         nwkp=(365+ILEAP+JULDAY -midwk1(nyrp)+7)/7
      ELSE
         nyrp=IYR
         nwkp=iwk
      ENDIF
      ndwkp=midwk1(nyrp)+(nwkp-1)*7

c Julian day we want minus the center of week Julian day preceeding it
      FNUMER = FLOAT( JULDAY- ndwkp )
c correct it if previous week center is in previous year.
      IF(FNUMER.LT.0.) FNUMER = FNUMER + 365.+ILEAP
      FDENOM = 7.0 -1.0
      IF(nyrp.EQ.89 .AND. nwkp.EQ.53  ) FDENOM=4. -1.0
      WT = FNUMER / FDENOM

      PRINT *, ' nwkp,nyrp,ndwkp, JULDAY, WT = '
     &          ,nwkp,nyrp,ndwkp,JULDAY,WT
 
c Years start one weekday later each year, except after leap years,
c when they advance 2 days.
c
C     calender
c yr st end   file dates
c------------------------
c NESDIS centered on Wednesday
c 98  Th-Th   
c 97  W - W
cl96  M -Tu
c 95  Su-Su
c 94  Sa-Sa
c 93  F - F
cl92  W -Th
c 91  Tu-Tu
c 90  M - M
c COADS  centered on Sunday
c 89  Su-Su
cl88  F -Sa   87/12/31-88/12/28
c 87  Th-Th
c 86  W - W
c 85  Tu-Tu
cl84  Su- M
c 83  Sa-Sa
c
c Aha!  web page 
c     ftp://ncardata.ucar.edu/datasets/ds277.0/oi/wkly/oiweek.info
c     says 
c     1981-89 data are from COADS  data with weeks centered on Sunday. 
c     1990-97 data are from NESDIS data with weeks centered on Wednsday
c ?
c  It also says a bias correction to the OI method introduced some 
c  temporal noise, which they *strongly* recommend filtering; does our 
c  weighting from bracketing data times accomplish this?
c ?
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE MKSST( IUNIT, TSCCM, SST1, SST2, TOPOGM, LANDUS
     A                , NLON2, NLAT2, NYRP, NWKP, WT )
 
      INTEGER LANDUS(NLON2,NLAT2)
      REAL TSCCM(NLON2,NLAT2), TOPOGM(NLON2,NLAT2)
      REAL SST1(NLON2,NLAT2), SST2(NLON2,NLAT2)
 
C ******           INITIALIZE SST1, SST2 (NEEDED FOR 82 JAN CASE)
      DO LON=1,NLON2
      DO LAT=1,NLAT2
         SST1(LON,LAT) = 0.
         SST2(LON,LAT) = 0.
      ENDDO
      ENDDO
 
      IF(NYRP.EQ.-1) THEN
         WT = 1.
         GOTO 15
      ENDIF
 
C ******           READ IN MM4 weekly SST DATASET
   10 READ(IUNIT,END=998) NWK, NYEAR, ((SST1(I,J),J=1,NLAT2),I=1,NLON2)
      PRINT*,'Read Reynolds data for week, year ',NWK,NYEAR
      IF( (NYEAR.NE.NYRP) .OR. (NWK.NE.NWKP) ) GOTO10
      PRINT *, 'Keeping Reynolds SST DATA:', NWK, NYEAR
 
C ******           READ IN MM4 weekly SST DATASET
   15 READ(IUNIT) NWK, NYEAR, ((SST2(I,J),J=1,NLAT2),I=1,NLON2)
      PRINT *, 'READING MM4 SST DATA:', NWK, NYEAR
 
      PRINT *, 'SST1, SST2, WT = '
     A       ,  SST1(1,1), SST2(1,1), WT
      DO 20 I=1,NLON2-1
      DO 20 J=1,NLAT2-1
      IF( (TOPOGM(I,J).LT.0.) .AND. (LANDUS(I,J).EQ.7) ) THEN
         TSCCM(I,J) = (1.-WT) * SST1(I,J) + WT * SST2(I,J)
      ENDIF
   20 CONTINUE
 
      REWIND(IUNIT)
 
      RETURN
  998 STOP 12
      END

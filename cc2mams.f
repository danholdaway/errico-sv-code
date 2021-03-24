      PROGRAM cc2mams

c kdr change increment on xtimef if using 6. hour data !!!
c analyses before 1992 may have only 14 levels; check tech note or Web
c NLEV1 is the number of pressure levels in the analysis

C
C A      SET PARAMETERS
      LOGICAL LTEST, LSNOWC
C
C  X X X X X X X                                      X X X X X X X X X
C  X X X X X X X    USER DEFINED PARAMETERS FOLLOW    X X X X X X X X X
C
C  X X X X  SET 1 :   TRUNCATION PARAMETERS FOR CCM1 DATASET   X X X X
C A1 PARAMETER (NLEV1=14,MMAX=42,NMAX=42,KMAX=42,NLAT1=64,NLON1=128)
c kdr  data > 12/31/91 has 15 levels
C
      PARAMETER (NLEV1=15,MMAX=42,NMAX=42,KMAX=42,NLAT1=64,NLON1=128)
c kdr      PARAMETER (NLEV1= 7,MMAX=42,NMAX=42,KMAX=42,NLAT1=64,NLON1=128)
C
C  NLEV1  = NUMBER OF FULL-SIGMA (DATA) LEVELS IN CCM DATASET.
C  MMAX   = MAXIMUM ZONAL WAVENUMBER.
C  NMAX   = MAXIMUM ORDER OF LEGENDRE POLYNOMIAL FOR WAVENUMBER 0.
C  KMAX   = MAXIMUM ORDER OF LEGENDRE POLYNOMIAL FOR ANY WAVENUMBER.
C  NLAT1  = NUMBER OF LATITUDES ON CCM GRID.
C  NLON1  = NUMBER OF LONGITUDES ON CCM GRID.
C
C  X X X X X   SET 2 :   TRUNCATION PARAMETERS FOR MM4 DATASET   X X X X
C
c                                NJ       NI
c kdr      PARAMETER (NLEV2=14,NLON2=41,NLAT2=36)
c Can low 80km    PARAMETER (NLEV2=10,NLON2=54,NLAT2=43)
c Hammill US 30km PARAMETER (NLEV2=13,NLON2=160,NLAT2=96)
c Gulf 50 and 80 PARAMETER (NLEV2=13,NLON2=86,NLAT2=67)
c Gulf 50x.02  (20 full levels) PARAMETER (NLEV2=19,NLON2=86,NLAT2=67)
c Gulf 50x.01  (20 full levels)
        PARAMETER (NLEV2=20,NLON2=86,NLAT2=67)
c Expl2 PARAMETER (NLEV2=10,NLON2=122,NLAT2=92)
c Rolf 40x60 60km PARAMETER (NLEV2=11,NLON2=60,NLAT2=40)
c   MISTAKE; he used 13 full sigma levels (12 half); the .9 was dropped 
c            from the output I read to set up this case
c Jrain W2(10lvls) W4(20lvls)(from Pred D4T1; 50x80x10)
C      PARAMETER (NLEV2=20,NLON2=80,NLAT2=50)
c Jrain 80km S2(20lvls) 10 level was done in Resol1  
c       same case as J89/June89
C      PARAMETER (NLEV2=20,NLON2=54,NLAT2=43)
C Jrain 80km S3(20lvls)  CONVLOW
C      PARAMETER (NLEV2=20,NLON2=70,NLAT2=50)
c J89EAST      same case as J89/June89
C      PARAMETER (NLEV2=16,NLON2=65,NLAT2=65)
c Pred D3T1; 60x75x20 PARAMETER (NLEV2=20,NLON2=75,NLAT2=60)
C      PARAMETER (NLEV2=13,NLON2=54,NLAT2=42)

C    ---  MUST ALSO CHANGE DIMENSIONS IN INTV1, ETC.  ---

C  NLEV2  = NUMBER OF HALF-SIGMA (DATA) LEVELS IN MM4 DATASET.
C  NLON2  = NUMBER OF LONGITUDES ON OUTPUT GRID.
C  NLAT2  = NUMBER OF LATITUDES ON OUTPUT GRID.
C
C  X X X  SET INTERMEDIATE REGULAR LAT/LON GRID DESIGNATED GRID 4  X X X
C
      PARAMETER (IFCT=1)
      PARAMETER (NLON4=NLON2*IFCT,NLAT4=NLAT2*IFCT)
C
C  NLON4 = NUMBER OF E-W POINTS ON THE INTERMEDIATE LAT/LON GRID.
C  NLAT4 = NUMBER OF E-W POINTS ON THE INTERMEDIATE LAT/LON GRID.
C  IFCT  = RATIO OF LAT-LON GRID DIMENSIONS TO THOSE OF MM4 GRID.
C
C  X X X  SET NUMBER OF MASS STORE VOLUMES OF CCM DATA TO BE READ  X X X
C
c kdr      PARAMETER (MAXVOLS=22)
      PARAMETER (MAXVOLS=1)
C
C  X X X  SET NUMBERS OF FIELDS TO BE SPECTRALLY PROJECTED  X X X
C
      PARAMETER (IP3DSC=3, IP2DSC=3, IP3DVE=2, IP2DVE=0)
C
C  IP3DSC = NUMBER OF 3-DIM. SCALAR FIELDS TO BE INTERP. BY PROJECTION.
C           E.G., T, Q
C  IP2DSC = NUMBER OF 2-DIM. SCALAR FIELDS TO BE INTERP. BY PROJECTION.
C           E.G., PHIS, PS, TS
C  IP3DVE = NUMBER OF 3-DIM. VECTOR FIELDS TO BE INTERP. BY PROJECTION.
C           E.G., U, V
C  IP2DVE = NUMBER OF 2-DIM. VECTOR FIELDS TO BE INTERP. BY PROJECTION.
C
C
C  X X X  SET NUMBER OF SINGLE LEVEL FIELDS FROM CCM USED IN SNOW  X X X
C
      PARAMETER (LSNOWC=.F., JFIELD=4) ! LSNOWC=.T. NOT WORKING YET
C
C  X X X X X X X    END OF USER DEFINED PARAMETERS    X X X X X X X X X
C  X X X X X X X                                      X X X X X X X X X
C
C A2     PARAMETERS FOR POINTERS TO DATA FIELDS.
c kdr for ccr1time 
      PARAMETER (NSIGS=NLEV1*2+1)
c
      PARAMETER (NVECTA=IP2DSC+IP3DSC*NLEV1, NVECTB=IP2DVE+IP3DVE*NLEV1
     A,          NVECTC=IP2DSC+IP3DSC*NLEV2, NVECTD=IP2DVE+IP3DVE*NLEV2
     B,          NVECTE=NVECTC+NVECTD, NVECTS=NVECTA+NVECTB
     C,          IFIELD=IP2DSC+IP3DSC+IP2DVE+IP3DVE  )
C
C A3     PARAMETERS FOR CCM GRID TO SPECTRA TRANSFORMS.
      PARAMETER (NLONP2=NLON1+2, MMAXP1=MMAX+1, NMAXP1=NMAX+1
     A,    NSPECS=MMAXP1*NMAXP1-(MMAX+NMAX-KMAX)*(MMAX+NMAX-KMAX+1)/2
     B,    NCOEFS=NSPECS*NVECTS*2, NPMN=MMAXP1*(NMAX+3)
     C,    NBDAT1=NLONP2*(30*NLEV1+50)+7  )
C
C A4     PARAMETERS FOR SPECTRA TO MM4 GRID TRANSFORMS.
      PARAMETER (NBDAT2=NLON2*NLAT2*NVECTS, NWORK2=2*NVECTS*MMAXP1)
C
C A5     DIMENSION VARIABLES FOR CCM GRID TO SPECTRA TRANSFORM.
      COMMON /COMFFT/ IFAX(13),INC,JUMP,NLONT,TRIGS(3*NLON1/2+1)
      DIMENSION COEFS(2,NSPECS,NVECTS), PMN(NPMN), HMN(NPMN)
c kdr NVECTS is the number of levels in all the 3D and 2D fields
     &,         FIELDS(NLONP2,NLAT1,NVECTS)
     A,         B1(NBDAT1), B1R(NLONP2,NVECTS)
     B,         MTRUNC(NMAXP1), NTRUNC(MMAXP1), WORK1(NVECTS*NLONP2)
     C,         GAUSWT(NLAT1), GLAT(NLAT1), SIGMA1(NLEV1)
c kdr for ccr1time
     &,         ASIGMA(NSIGS), BSIGMA(NSIGS),GLATCCM(NLAT1)
     D,         NPOINT(IFIELD,2), RCOEF(2,NMAXP1,MMAXP1,NVECTS)
     E,         OPOINT(IFIELD,3),SIGMAR(NLEV1)
      CHARACTER*8  CHRD(NBDAT1), CFIELD(IFIELD)
     A,         CFLDSC(IP2DSC+IP3DSC), CFLDVE(IP2DVE+IP3DVE)
      EQUIVALENCE (CFIELD,CFLDSC), (CFIELD(IP2DSC+IP3DSC+1),CFLDVE)
C
C A6     DIMENSION VARIABLES FOR LAT-LON GRID
      DIMENSION B2(NLON4,NLAT4,NVECTS),B2R(NLON4,NVECTS),EPSI(NPMN)
     A,         G4LAT(NLAT4),G4LON(NLON4),WORK2(NWORK2)
C
C A7     DIMENSION VARIABLES FOR MM4 HORIZONTAL GRID (CCM VERTICAL)
      DIMENSION B3(NLON2,NLAT2,NVECTS),B2PD(NLON2,NLAT2)
     A,         XLONS(NLON2,NLAT2,2), XLATS(NLON2,NLAT2,2)
C
C A8     DIMENSION VARIABLES FOR MM4 INPUT FILE
c kdr                                   hydrost. geop. hght
      DIMENSION B4(NLON2,NLAT2,NVECTE), B3H(NLON2,NLAT2,NLEV2)
c kdr                               1&2 are mapscale factors
c kdr           dot pnt pressure    see EQUIV for others
     A,         B3PD(NLON2,NLAT2), XMAPS(NLON2,NLAT2,5)
     B,         SIGMA2(NLEV2), SIGMAF(NLEV2+1), DSIGMA(NLEV2)
     C,         CORIOL(NLON2,NLAT2), XLANDU(NLON2,NLAT2)
     D,         SNOWCV(NLON2,NLAT2), LANDUS(NLON2,NLAT2)
     E,         TOPOGM(NLON2,NLAT2)
      EQUIVALENCE (CORIOL,XMAPS(1,1,3)), (XLANDU,XMAPS(1,1,4))
     A,           (SNOWCV,XMAPS(1,1,5))
      CHARACTER CGTYPE*8
C
C A9     DIMENSION VARIABLES FOR SNOW COVER ROUTINE
      CHARACTER*8 CFLDB5(JFIELD)
      DIMENSION KPOINT(JFIELD,2), B5(NBDAT1), B5R(NLON1,JFIELD)
C
C A10    NAMELISTS, DATA DEFINITIONS, ETC.
      EQUIVALENCE (B1,B3),(B1R,B2R,B5) ,(COEFS,RCOEF,B4)
      DIMENSION NTIMES(MAXVOLS), NREADS(300), IERROR(2)
 
C DIMENSION SURFACE TEMPERATURE ON CCM SURFACE; NOT GIVEN BY ECMWF DATA
CCC   DIMENSION TSCCM(NLON2,NLAT2)
      REAL TSCCM(NLON2,NLAT2), SST1(NLON2,NLAT2), SST2(NLON2,NLAT2)
      COMMON /CDATEI/ IDATE
 
C ******           ARRAYS NEEDED FOR NEW CALCUATION OF P*
      REAL PA(NLON2,NLAT2), ZA(NLON2,NLAT2), TLAYER(NLON2,NLAT2)
 
C DIMENSION THESE 'B' ARRAYS FOR TESTING/PLOTTING SSTS
      REAL  BLON(NLAT2,NLON2), BLAT(NLAT2,NLON2), BCOR(NLAT2,NLON2)
      REAL  BSST(NLAT2,NLON2)
 
C DIMENSION DUMMY ARRAYS NEEDED IN INTV1, INTV2, INTV3, INTGTB
      REAL DUM1(NLON2,NLAT2,NLEV1), DUM2(NLON2,NLAT2,NLEV1,2)
      REAL C1(NLEV1), C2(NLEV2)
 
      COMMON /CPOINT/ NPHIS1,NPS1,NTSP1,NTP1,NQP1,NZA1,NUV1
     A               ,NPHIS2,NPS2,NTSP2,NTP2,NQP2,NZA2,NUV2
      NAMELIST /LFILES/ NVOLUS,NTIMES,NREADS
      NAMELIST /LGRID2/ PTOPCB,SIGMAF,LGTYPE,CLAT,CLON,DELX
      DATA IUNITA/21/, IUNITB/11/, IUNITS/12/, IUNITC/20/, IUNITE/33/
      DATA NREADS/300*0/,NTIMES/MAXVOLS*1/
c kdr
      DATA LGTYPE/8HLAMBERT /
c gulf50 &     ,PTOPCB/10./,CLON/-97./,CLAT/31./,DELX/5.0E4/
c gulf80     &     ,PTOPCB/10./,CLON/-97./,CLAT/31./,DELX/5.0E4/
c Expl2 &     ,PTOPCB/10./,CLON/-75./,CLAT/37./,DELX/6.0E4/
c wrong hammill?     &     ,PTOPCB/10./,CLON/-90./,CLAT/40./,DELX/8.0E4/
c Rolf 4060 60km &     ,PTOPCB/2.5/,CLON/-115./,CLAT/40./,DELX/6.0E4/
c Pred D3T1 &     ,PTOPCB/1.0/,CLON/-115./,CLAT/41./,DELX/4.0E4/
      DATA PMN/NPMN*1./,HMN/NPMN*1./
C
      DATA LTEST/.FALSE./
      DATA CFLDSC/'PHIS    ','PS      ','T       ','Q       ',
     A            'ORO     ','ZA      '/
      DATA CFLDVE/'U       ','V       ' /
C
C  B1 IS FOR DATA RECORDS FROM CCM GAUSSIAN GRID
C  B2 IS FOR LAT-LON GRID WITH CCM VERTICAL STRUCTURE
C  B3 IS FOR MM4 HORIZONTAL GRID, BUT CCM VERTICAL GRID
C  B4 IS FOR MM4 3-DIMENSIONAL GRID
C  B5 IS FOR 2-D FIELDS OF CCM DATA FOR SPECIAL INTERPOLATION ROUTINE.
C
C  CFLDSC IS LIST OF SCALAR FIELDS TO BE INTERPOLATED BY SPECTRAL PROJ.
C  CFLDVE IS LIST OF VECTOR FIELDS TO BE INTERPOLATED BY SPECTRAL PROJ.
C  OPOINT IS POINTER INFO FROM CCM HEADER (FOR B1 ARRAY).
C  NPOINT IS POINTER INFO FOR REFORMATTED DATA (B1R).
C  CFLDB5 IS LIST OF 2-D FIELDS TO BE INTERPOLATED IN SPECIAL ROUTINE.
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X


c NLM parameters and arrays
      
      parameter (ni=NLAT2, nj=NLON2, nk=NLEV2) 
      parameter (imoist=1,  nij=ni*nj, nijk=ni*nj*nk, nkp1=nk+1
     &          ,ni1=ni-1, nj1=nj-1
     &          ,nfieldh=4, nfieldq=3+imoist
     &          ,nfieldr2=2, nfield3w=nfieldq 
     &          ,nfield2w=nfieldr2 
     &          ,nfieldsw=nfield3w+nfield2w+nfieldh
c                           4       2          4
     &          ,nfieldsx=30, n1charsw=128, ndcharsw=6
     &          ,nnftimes=100, nrvars=32+2*nk, nlvars=30
     &          ,nnvars=50+nk+nnftimes+2*nfieldsw)
      parameter (nwork1=nij*nkp1, nfield3m=12, nworkr=2+nij*nfield3m)
      parameter (nsplitmx=3, npkappa=1400, nqsatvp=844)
c
      logical lrain,lground ,lprint,lcoefout
      logical lvarsw(nlvars),  lmoist,lcloud,lsnow,ldrycon,
     &               lbdry,lbdsame,lsplit,llinear
      dimension sigf(nkp1),tg(ni,nj),tdeepg(ni,nj)
      real rvarsw(nrvars), tpmeanb(nkp1),work(ni,nj,nk)
     &,    fieldh(nij,nfieldh), fields3(nijk,nfield3w)
     &,    fields2(nij,nfield2w), t(nijk)
      integer nvarsw(nnvars),  igtype(nfieldsw), iflds(3),msteps(nk)
      character*(8) cfieldsw(nfieldsx)
      character*(n1charsw) charsw(ndcharsw)
      equivalence  (tg,fields2(1,2))
      equivalence (tdeepg,fieldh(1,4))
      data cfieldsw/'      pu','      pv','      pt','      pq'
     &             ,'      ps','      tg','   topog','  sfctyp'
     &             ,'   snowc','  tdeepg',20*'        '/
     &    ,charsw /6*'  '/
     &    ,igtype /2*0,8*1/
     &    ,xtimef  /-1./

C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C                            S T A R T
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
C B      COMPUTE VARIOUS CONSTANTS USED AT ALL NTIME
C
C B1     READ NAMELISTS
      READ (13,LFILES)
      WRITE(6,LFILES)
      READ (13,LGRID2)
      WRITE(6,LGRID2)
C
C B2     SET CONSTANTS FOR CCM GRID TO SPECTRA TRANSFORM.
C        (FFTFAX IS NCAR CRAY LIBRARY ROUTINE).
      PIR180=ATAN(1.)/45.
      INC=1
      JUMP=NLONP2
      NLONT=NLON1
      CALL FFTFAX (NLON1,IFAX,TRIGS)
      CALL MTRA (MTRUNC,NTRUNC,MMAX,NMAX,KMAX)
      CALL GAUAW (GLAT,GAUSWT,WORK1,NLAT1)
C
C B3     SET CONSTANTS FOR SPECTRA TO LAT-LON GRID TRANSFORM.
      WRITE (CGTYPE,'(A8)') LGTYPE
      NMAXPX=NMAX+2
      IF (MOD(NMAXPX,2).NE.0) NMAXPX=NMAX+3
      CALL EPSIL2 (EPSI,NMAXPX,MMAXP1)
      CALL GRIDML (XLONS,XLATS,XMAPS,CORIOL,TOPOGM,XLANDU,LANDUS,SNOWCV
     A,            CLON,CLAT,DELX,GRDFAC,NLON2,NLAT2,IUNITB,CGTYPE)

      CALL GRID4  (XLONS(1,1,2),XLATS(1,1,2),NLON2,NLAT2
     A,            G4LON,G4LAT,NLON4,NLAT4)
      CALL SETSIG (SIGMAF,SIGMA2,NLEV2)
C
C B4     DETERMINE POINTERS FOR FIELDS REQUIRED BY SNOW COVER
      IF (LSNOWC) CALL READPT (OPOINT,KPOINT,CFLDB5,JFIELD,NLON1)
C
C
C C  BEGIN LOOP OVER CCM HISTORY-FILE VOLUMES
C
C
      NTCNTA=0
      IUNITA=IUNITA-1
      DO 99 NVOL=1,NVOLUS
      IUNITA=IUNITA+1
      PRINT 66,NVOL,IUNITA
C
C D      BEGIN LOOP OVER NTIMES
C
      DO 98 NTIME=1,NTIMES(NVOL)
      NTCNTA=NTCNTA+1
      IF (NREADS(NTCNTA).NE.1) THEN
         READ (IUNITA) LENHD,MTYPE
         PRINT*,'MTYPE = ',MTYPE
         MTYPE = MOD(MTYPE,10)
         PRINT*,'MTYPE = ',MTYPE
         IF (MTYPE.EQ.1) THEN
            NSKIPR=NLAT1         ! NUMBER OF RECORDS AFTER FIRST RECORD
         ELSEIF (MTYPE.EQ.2 .OR. MTYPE.EQ.3) THEN
            NSKIPR=NLAT1+2       ! NUMBER OF RECORDS AFTER FIRST RECORD
         ELSE
            PRINT 63,MTYPE
            CALL ABORT ('NUMBER OF RECORDS TO SKIP UNKNOWN')
         ENDIF
         CALL SKIPR ( IUNITA,NSKIPR,IERROR )
         PRINT 64,NTIME,IUNITA
      ELSE         ! NREADS(NTCNTA) = 1
         CALL ZEROA (COEFS,NCOEFS)
         CALL CCR1TIME (FIELDS,ASIGMA,BSIGMA,GLATCCM,OPOINT,
     $        NDATE,NTTIME,NSTEP,CFIELD,NLONP2,NVECTS,IFIELD,NSIGS,
     $        NLON1,NLAT1,NLEV1,MMAX,NMAX,KMAX,IUNITA)
c ps in Pa here
 678   FORMAT(A,1p,10(/10(E11.4,1X)))
      WRITE(*,*) 'ASIGMA = ',ASIGMA

C CHANGE ORDER OF VERTICAL INDEXES FOR PRESSURE LEVELS
c SIGMAR goes from top to bottom, just like the p-levels and data
c on the CCM tapes, but opposite of the p-levels and data in
c the WMO archive.
          DO N=1,NLEV1


C This omits the top (sigma(?) = 0.0).  It picks out the number of half
C levels that there are, but from amoung the full levels.
C Is this correct?
C  It picks out the right numbers( the pressure levels, in bars, of the data), 
C but maybe coincidentally , and uses them as sigmas
C This may not be so bad if ptop(ccm)=0. and ps is taken to be 100cb
C as it seems to be in INTGTB



             SIGMAR(N) = ASIGMA(2*N+1)
             SIGMA1(NLEV1-N+1) = SIGMAR(N)
          END DO
          WRITE(*,*) SIGMA1
c
 
C        SET POINTERS
C
c kdr NPOINT filled based on info in OPOINT
         CALL NEWPNT (NPOINT,OPOINT,CFIELD,IFIELD,NLONP2,0)
c CPOINT common block filled from info in NPOINT
         CALL SETPNT (NPOINT,CFIELD,IFIELD,NLEV1,NLEV2)
C
C E      BEGIN LOOP OVER LATITUDES FOR GAUSSIAN GRID TO SPECTR TRANFORM
C
         DO 11 KLAT=1,NLAT1
            JSIGN=1
            IF (KLAT.LE.NLAT1/2) JSIGN=-1

c these calculate Legendre polynomials and derivs of them.
            ABGLAT=ABS(GLAT(KLAT))
            CALL ALPNM2(PMN,NMAXPX,MMAXP1,ABGLAT,EPSI)
            CALL ALPDR2(HMN,PMN,NMAXPX,MMAXP1,EPSI)
c            IF (LSNOWC) CALL SFILE1 (B1,B5,B5R,KPOINT,OPOINT,NBSIZE
c     A         ,NLON1*JFIELD,JFIELD,NLON1,NLONW,KLAT,NLAT1,IUNITE)

            B1R(1:NLONP2,1:NVECTS) = FIELDS(1:NLONP2,KLAT,1:NVECTS)

c ps (B1R(NPS1) changed to log(ps/10^5)
            CALL GRD2SC (B1R,B1R(1,NPHIS1),B1R(1,NPS1),B1R(1,NUV1)
     A,            PMN,HMN,COEFS,WORK1,MTRUNC,NVECTS,JSIGN,KLAT
     B,            NLONP2,NLEV1,MMAXP1,NMAXPX,NMAXP1,NSPECS,NUV1
     C,            GAUSWT(KLAT),GLAT(KLAT),NTRUNC,NVECTA)
   11    CONTINUE
C
C        END OF TRANSFORM FROM DATA TO SPECTRAL COEFFICIENTS.
C
C
C
C F      BEGINNING OF SPECTRAL TO LAT-LON GRID PROJECTION
C
         DO 50 J=1,NLAT4
            XLAT4=SIN(G4LAT(J)*PIR180)
            CALL ALPNM2  (PMN,NMAXPX,MMAXP1,XLAT4,EPSI)
c ps B2R is equived to B1R, so ps is log (ps/10^5) here
c ps (B2R (NPS1)) is reset to 100*exp (ps), which is ps(Pa)/1000. = cb
            CALL SC2GRD4 (B2R,WORK2,COEFS,PMN,1,MTRUNC,NTRUNC,1
     A,             G4LON,NLON4,B2R(1,NPS1),NSPECS,NVECTA
     B,             MMAXP1,NMAXP1,NMAXPX,1)
            CALL ALPDR2  (HMN,PMN,NMAXPX,MMAXP1,EPSI)
c ps possible no changes to ps here; I haven't traced it all though
            CALL SC2GRD4 (B2R(1,NVECTA+1),WORK2,COEFS(1,1,NUV1),PMN,HMN
     A,             MTRUNC,NTRUNC,XLAT4,G4LON,NLON4,1,NSPECS,NVECTB
     B,             MMAXP1,NMAXP1,NMAXPX,2)
            DO 40 I=1,NLON4
               DO 30 L=1,NVECTA
                  B2(I,J,L)=B2R(I,L)
   30          CONTINUE
               DO 20 L=1,NVECTB
                  B2(I,J,L+NVECTA)=B2R(I,L+NVECTA)
   20          CONTINUE
   40       CONTINUE
   50    CONTINUE
C
C HORIZONTAL INTERPOLATION OF BOTH THE SCALAR AND VECTOR FIELDS
C
c scalar
         CALL BILINX (B3,B2,XLONS,XLATS,G4LON,G4LAT
     A,            NLON2,NLAT2,NLON4,NLAT4,NVECTA,NLON2-1,NLAT2-1,1)
         PRINT*,'After BILINX '
         WRITE(*,678) 'B2(10,10,NPS1) ps = ',B2(10,10,NPS1)
         WRITE(*,678) 'B3(10,10,NPS1) ps = ',B3(10,10,NPS1)
         WRITE(*,678) 'B2(10,10,1:15) T = '
     &               ,(B2(10,10,NTP1+K-1),K=1,NLEV1)
         WRITE(*,678) 'B3(10,10,1:15) T = '
     &               ,(B3(10,10,NTP1+K-1),K=1,NLEV1)
         CALL BOUNDF (B3,0.,NLON2,NLAT2,NVECTA)
c vector
         CALL BILINX (B3(1,1,NVECTA+1),B2(1,1,NVECTA+1),XLONS(1,1,2)
     A,            XLATS(1,1,2),G4LON,G4LAT,NLON2,NLAT2,NLON4,NLAT4
     B,            NVECTB,NLON2,NLAT2,2)
C
C ROTATE U-V FIELDS AFTER HORIZONTAL INTERPOLATION
C
         CALL UVROT4 (B3(1,1,NVECTA+1),XLONS(1,1,2),CLON,GRDFAC
     A            ,NLON2,NLAT2,NVECTB/2)
C
C
C
C X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
C                  V E R T I C A L   I N T E R P O L A T I O N
C
C X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
C
C F1     COPY MM4 TOPOGRAPHY
         DO 60 I=1,NLON2
         DO 60 J=1,NLAT2
            B4(I,J,NPHIS2)=TOPOGM(I,J)
   60    CONTINUE

C ******           NEW CALCULATION OF P* ON MM4 TOPOGRAPHY.
         CALL INTGTB( PA, ZA, TLAYER, B4(1,1,NPHIS2)
     A           , B3(1,1,NTP1), B3(1,1,NZA1)
     A           , SIGMAR, NLON2, NLAT2, NLEV1, DUM1, DUM2 )
 
c kdr B4(1,1,NPS2) calculated here  PA*exp(B4(PHI)-ZA)/TLAYER - PTOPCB
c PA is in cbar
         CALL INTPSN ( B4(1,1,NPS2), B4(1,1,NPHIS2)
     A            , PA, ZA, TLAYER, PTOPCB, NLON2, NLAT2 )
         WRITE(*,678)'B4(NPHIS2) phi mm4 = ',B4(10,10,NPHIS2)
         WRITE(*,678)'B4(NPS2) ps mm4 = ',B4(10,10,NPS2)

         CALL P1P2 (B3PD,B4(1,1,NPS2),NLON2,NLAT2)
C
C F0    DETERMINE SURFACE TEMPS ON MM4 TOPOGRAPHY.
C       INTERPOLATION FROM PRESSURE LEVELS AS IN INTV2
c kdr
        PRINT*,'DON"T FORGET;  B4(NJ,NI,)'
              WRITE(*,679)'Before INTV3 B4(i,i',NPS2,') = '
     &                    ,(B4(i,i,NPS2),i=1,NLAT2-1)
 679   FORMAT(A,I3,A,10(/20(F5.1,1X)))
c B4(,,NPS2) is in cbar here
c                     returned       input

c !? Why is B4(NTSP2)(surface temp?)  created from a level of B3(NTP1) 
c    (Temp)  NTSP2 seems to be labelled 'ORO', and there's no surf T
c            in the lists.  But TOPOGM is stored in B4(NPHIS2)
c            and topog (fieldh()) is set from TOPOGM directly
c    Then B4(NTSP2) is passed to MKSST and reassigned new values!

         CALL INTV3 ( B4(1,1,NTSP2), B3(1,1,NTP1),B4(1,1,NPS2),SIGMAR
     A,             NLON2,NLAT2,NLEV1,DUM1)
         PRINT*,'NTSP2 = ',NTSP2
         WRITE(*,678) 'B4(10,10,1:10)after INTV3  ' ,B4(10,10,NTSP2)
 
C F1    CALCULATE SSTS FOR DATE FROM OBSERVED SSTS (1982-1988)
         PRINT *, 'INPUT DAY FOR SST DATA ACQUISITION:', IDATE
         CALL JULIAN( IDATE, NYRP, NMOP, WT )
         CALL MKSST( IUNITS, B4(1,1,NTSP2), SST1, SST2, TOPOGM, LANDUS
     A          , NLON2, NLAT2, NYRP, NMOP, WT )
         WRITE(*,678) 'B4(10,10,1:10)after MKSST  ',B4(10,10,NTSP2)
 
C F3     INTERPOLATE U, V, T, AND Q.
         CALL INTV1 (B4(1,1,NUV2),B3(1,1,NUV1),B3PD,B2PD
     B,                 SIGMA2,SIGMAR,PTOPCB,NLON2,NLAT2,NLEV2
     B,                 NLEV1,NLON2,NLAT2,2,1,DUM2,C1,C2)
C
c B3(NPS1) is set = 100. here, and is used in conjuction with PTOPCB
c and the MM4 pressures B4(NPS2)
         PRINT*,'Before INTV2'
         WRITE(*,678) 'B3(10,10,1:15) T = '
     &               ,(B3(10,10,NTP1+K-1),K=1,NLEV1)
         CALL INTV2 (B4(1,1,NTP2),B3(1,1,NTP1),B4(1,1,NPS2),B3(1,1,NPS1)
     C,                 SIGMA2,SIGMAR,PTOPCB,NLON2,NLAT2,NLEV2,NLEV1,1
     C,                 DUM2,C1,C2)
         PRINT*,'After INTV2'
         WRITE(*,678) 'B4(10,10,NTP2) T = '
     &               ,(B4(10,10,NTP2+K-1),K=1,NLEV2)
         WRITE(*,678) 'B4(10,10,1:10) T = '
     &               ,(B4(10,10,3+K-1),K=1,NLEV2)
 
 
         CALL HUMID1 (B3(1,1,NTP1),B3(1,1,NQP1),B3(1,1,NPS1)
     D,            0.,SIGMA1,NLON2,NLAT2,NLEV1)
c B3(NPS1) is set = 100. here, and is used in conjuction with PTOPCB
c and the MM4 pressures B4(NPS2)
         CALL INTV1 (B4(1,1,NQP2),B3(1,1,NQP1),B4(1,1,NPS2),B3(1,1,NPS1)
     E,                 SIGMA2,SIGMAR,PTOPCB,NLON2,NLAT2,NLEV2
     E,                 NLEV1,NLON2-1,NLAT2-1,1,0,DUM2,C1,C2)
         CALL HUMID2 (B4(1,1,NTP2),B4(1,1,NQP2),B4(1,1,NPS2)
     F,            PTOPCB,SIGMA2,NLON2,NLAT2,NLEV2)
C
C F4     DETERMINE geopot. hght H = B3H
         CALL HYDROST(B3H,B4(1,1,NTP2),B4(1,1,NPHIS2),B4(1,1,NPS2)
     &,        PTOPCB
     A,        SIGMAF,SIGMA2,DSIGMA,NLON2,NLAT2,NLON2-1,NLAT2-1,NLEV2)

c save uncoupled T
         PRINT*,'B4 before MULTP'
         M = 0
         DO 169 K=1,NLEV2
         DO 169 I=1,NLON2
         DO 169 J=1,NLAT2
            M = M+1
            t(M) = B4(I,J,NTP2+K-1)
            IF (I.EQ.10 .AND. J.EQ.10) PRINT*,t(M)
 169     END DO

C F5     MULTIPLY FIELDS BY P*
         CALL MULTP (B4(1,1,NUV2),B3PD,        NLON2,NLAT2
     H,              NLEV2*2,NLON2  ,NLAT2  )
         CALL MULTP (B4(1,1,NTP2),B4(1,1,NPS2),NLON2,NLAT2
     I,              NLEV2  ,NLON2-1,NLAT2-1)
         WRITE(*,678) 'B4(10,10,1:10) p*T = '
     &               ,(B4(10,10,NTP2+K-1),K=1,NLEV2)
 
         CALL MULTP (B4(1,1,NQP2),B4(1,1,NPS2),NLON2,NLAT2
     J,              NLEV2  ,NLON2-1,NLAT2-1)
         CALL MULTP (B3H         ,B4(1,1,NPS2),NLON2,NLAT2
     K,              NLEV2  ,NLON2-1,NLAT2-1)
C
C G      WRITE AN INITIAL FILE FOR THE NLM

         xhour0 = MOD (IDATE,100)
         iday0 = MOD((IDATE/100),100)
         imonth0 =MOD((IDATE/10000),100)
         iyear0 = 1900 + (IDATE/1000000)
         PRINT*,'IDATE, xhour0, iday0, imonth0, iyear0 = '
     &          ,IDATE, xhour0, iday0, imonth0, iyear0

c kdr these loops fill the data arrays, flipping I and J in the process
c as well as converting I,J, and K into a single index M
         M = 0
         DO 179 K=1,NLEV2
         DO 179 J=1,NLAT2
         DO 179 I=1,NLON2
            M = (K-1)*NLAT2*NLON2 + (I-1)*NLAT2 + J
            fields3(M,1) = B4(I,J,NUV2+K-1)
            fields3(M,2) = B4(I,J,NUV2+K-1+nk)
            fields3(M,3) = B4(I,J,NTP2+K-1)
            fields3(M,4) = B4(I,J,NQP2+K-1)
 179     END DO
         M = 0
         DO 189 J=1,NLAT2
         DO 189 I=1,NLON2
            M = (I-1)*NLAT2 + J
            fields2(M,1) = B4(I,J,NPS2)
            fields2(M,2) = B4(I,J,NTSP2)
 189     END DO

         WRITE(*,678) 'pu,pv,pt,pq ',(fields3(100,IF),IF=1,4)
         WRITE(*,678) 'ps,tg  ',(fields2(100,IF),IF=1,2)

c fill extraneous edge values of cross fields with 0s
c NOT TESTED YET
         DO nfld=3,nfield3w
            CALL Cbdrycg (fields3(1,nfld),ni,nj,nk)
         END DO
         DO nfld=1,nfield2w
            CALL Cbdrycg (fields2(1,nfld),ni,nj, 1)
         END DO

         IF ( NVOL.EQ.1 .AND. xtimef.LT.0.) THEN
c  compute mean values
c kdr t is decoupled
c              call Bdecup (t,pt,ps,ni,nj,nk,1)
c             PRINT*,'t(:,:,10) = '
c             WRITE(*,97) (t(m),m=ni*nj*9+1,ni*nj*10)
c 97          FORMAT (9(/10E10.3,X),/E10.3)
c 97          FORMAT (7(/10E10.3,X))
             call Sfavg (tpmeanb,t,ni,nj,nk,1,1)
             call Sfavg (tpmeanb(nkp1),fields2(1,1),ni,nj,1,1,1)

            DO 100 J=1,NLAT2
            DO 100 I=1,NLON2
c               M = (NLON2-1)*J + I
               M = (I-1)*NLAT2 + J
c tg is filled through equivalence to fields2( ,2) above
               tdeepg(J,I) = tg(J,I)
               fieldh(M,1) = TOPOGM(I,J)
               fieldh(M,2) = XLANDU(I,J)
               fieldh(M,3) = SNOWCV(I,J)
 100        END DO

            ldrycon=.f.
            lmoist=.t.
            lground=.t.
            lrain=.F.
            lcloud=.F.
            lsnow=.F.
            lbdry=.F.
            lbdsame=.F.
            lsplit=.F.
            llinear=.F.
            lcoefout=.F.
            lprint=.T.
      
            iflds(1) = nfield3w
            iflds(2) = nfield2w
            iflds(3) = nfieldh
            ipackw = 1
            nspong = 0
            nsplit=0
            msteps(1:nk)=0
      
            sigf(1:nkp1) = SIGMAF(1:nkp1)

            toutreg = 12.*60.
            tbsinreg = 0.
            time0 = 0.
            ktauout = 0
            ktaubs=0
            dt = 0.
            dtbdr = 720.
            dx = DELX
            radfreq = 0.
            iformt = 1

            xtimef = 0.
            PRINT*,'xtimef = ',xtimef
c  write out header and initial data

            call Nsethead (rvarsw,lvarsw,nvarsw,toutreg,ptopcb,dx,
     &               clat,clon,xhour0,time0,dt,dtbdr,radfreq,
     &               tbsinreg,sigf,tpmeanb,lmoist,lground,lrain,
     &               lcloud,lsnow,lbdry,lbdsame,lsplit,llinear,lcoefout,
     &               linitial,ldrycon, igtype,iyear0,imonth0,iday0,
     &               nsplit,msteps,nnmis,iflds,npkappa,nqsatvp,nfieldsw,
     &               ni,nj,nk,nkp1,nspong,n1charsw,ndcharsw,nnvars,  
     &               nlvars,nrvars,nnftimes,1,iformt,ipackw,tpratio)

            ifp1=nvarsw(11)
            nvarsw(ifp1)=nnftimes
c  total number of times put out new.7; in Nsetoutp  1+ktimax/ktimout
            NCNT = 0
            nvarsw(ifp1+1) = 0
            DO 199 NV=1,NVOLUS
            DO 199 NT=1,NTIMES(NV)
               NCNT=NCNT+1
               IF (NREADS(NCNT).EQ.1) nvarsw(ifp1+1)= nvarsw(ifp1+1) +1
 199        END DO
            PRINT*,'nvarsw(ifp1+1) = ',nvarsw(ifp1+1)
            nvarsw(ifp1+2)= nvarsw(ifp1+1)  
            nvarsw(ifp1+3)= nvarsw(ifp1+1)

            PRINT*,'nvarsw(2),IUNITC,lprint = ',nvarsw(2),IUNITC,lprint
            call Nopen (nvarsw(2),IUNITC,lprint)
            PRINT*,'Must have called Nopen'
C
c      SUBROUTINE Nwriteh (nvars,lvars,rvars,cfields,chars,work,nij,
c     &                    ir,iunit,fieldh,cfldsrw,ifldsrw,lprint,
c     &                    nnvars,nrvars,nlvars,nworkr)
c
            call Nwriteh (nvarsw,lvarsw,rvarsw,cfieldsw,charsw,work,nij,
     &                    irw,IUNITC,fieldh,cfieldsw,iflds,lprint,
     &                    nnvars,nrvars,nlvars,nworkr)
         END IF      !NVOL=1 and xtimef<0.

         call Nwritet (fields3,fields2,work,xtimef,nvarsw,cfieldsw,
     &                 cfieldsw,iflds,iflds,
     &                 nij,nk,ktauout,irw,IUNITC,lprint,nnvars,nworkr)
         PRINT 61,NTIME,IUNITA,IUNITC
C X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
C                  E N D   O F   L O O P   O V E R   N T I M E
C
C X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      END IF      ! NREADS(NTCNTA) block begins at beginning of loop 98

c kdr increment file time 
      IF (xtimef.GE.0.) THEN
         xtimef = xtimef + 12.
         PRINT*,'xtimef = ',xtimef
      END IF

   98 CONTINUE
C ******      RELEASE NOT FOUND IN LIBRARIES: COMMENTED OUT, GTB 9/90
CCC   CALL RELEASE (IERROR,2HDN,FTNAME(IUNITA))
   99 CONTINUE
C
C H      FILE CREATION
      PRINT 65,IUNITC
      CALL Nclose (IUNITC,nvarsw(ifp1+1)+1,.T.)
C
   61 FORMAT('0 FIELDS FOR NTIME =',I3,' ON CCM UNIT ',I3,' INTERPOLATED
     A AND WRITTEN TO MM4 UNIT ',I3)
   63 FORMAT('0 MTYPE=',I2,'  CCM HEADER FORMAT UNKNOWN')
   64 FORMAT('0 FIELDS FOR NTIME =',I3,' ON CCM UNIT ',I3,' SKIPPED')
   65 FORMAT('0 JOB COMPLETED.  OUTPUT FILE ON UNIT',I3)
   66 FORMAT('0 CCM VOLUME NUMBER',I3,' ON UNIT',I3,' BEING PROCESSED')
C
      STOP 999
      END
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE ALPDR2 (DALP,ALP,LALP,LM,EPSI)
C
C     *****   OCT 1975  -  JOHN D. HENDERSON  ****
C     *****  SEPT 1984  -  ROGER DALEY        ****
C    * CALCULATES N-S DERIVATIVES OF EACH ASSOCIATED LEGENDRE POLYNOMIAL
C     * DALP(LALP,LM)   WILL CONTAIN N-S DERIVATIVE OF ALP.
C     * ALP(LALP,LM)   CONTAINS LEGENDRE POLYNOMIALS FOR ONE LATITUDE.
C     * EPSI(LALP,LM)   CONTAINS PREVIOUSLY CALCULATED CONSTANTS.
C
C     * WARNING - LALP MUST BE EVEN.
C     *         - LAST ELEMEMT OF EACH ROW IS SET TO ZERO.
C
C     LEVEL 2,DALP,ALP,EPSI
      DIMENSION DALP(LALP,1),ALP(LALP,1),EPSI(LALP,1)
C-----------------------------------------------------------------------
      LALPM=LALP-1
C
      DO 30 M=1,LM
C
      DO 20 N=1,LALPM
      FNS=FLOAT(M+N-2)
      ALPILM=0.
      IF(N.GT.1) ALPILM=ALP(N-1,M)
      DALP(N,M)=(FNS+1.)*EPSI(N,M)*ALPILM - FNS*EPSI(N+1,M)*ALP(N+1,M)
   20 CONTINUE
C
   30 DALP(LALP,M)=0.
C
      RETURN
      END
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE ALPNM2 (ALP,LALP,LM,SINLAT,EPSI)
      DIMENSION ALP(LALP,LM),EPSI(LALP,LM)
C
C  COMPUTES LEGENDRE POLYNOMIALS FOR GIVEN SINE OF LATITUDE.
C  FROM ROGER DALEY 16/9/83
C  POLYNOMIALS ARE NORMALIZED SUCH THAT ?
C
      LALPH=LALP/2
      LALPE=LALPH*2
C
      COS2=1.-SINLAT*SINLAT
      PROD=1.
      A=1.
      B=0.
C
C  LOOP 20 OVER ZONAL WAVE NUMBER 0 TO LM-1
C
      DO 30 M=1,LM
      FM=FLOAT(M-1)
      IF(M.EQ.1) GO TO 12
      A=A+2.
      B=B+2.
      PROD=PROD*COS2*A/B
C
C  COMPUTE THE FIRST TWO DIAGONALS
C
   12 ALP(1,M)=SQRT(.5*PROD)
      ALP(2,M)=SQRT(2.*FM+3.)*SINLAT*ALP(1,M)
C
C  NEXT, COMPUTE ELEMENTS 3 TO LR IN PAIRS
C
      DO 20 N=3,LALPE,2
      ALP(N,M)=(SINLAT*ALP(N-1,M)-EPSI(N-1,M)*ALP(N-2,M))/EPSI(N,M)
   20 ALP(N+1,M)=(SINLAT*ALP(N,M)-EPSI(N,M)*ALP(N-1,M))/EPSI(N+1,M)
C
      IF(LALPE.NE.LALP) THEN
        N=LALP
        ALP(N,M)=(SINLAT*ALP(N-1,M)-EPSI(N-1,M)*ALP(N-2,M))/EPSI(N,M)
      ENDIF
C
   30 CONTINUE
      RETURN
      END
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE BILINX (B3,B2,XLONS,XLATS,G4LON,G4LAT,NLON2,NLAT2
     A,                 NLON4,NLAT4,NVECTS,NIN,NJN,INDX)
C
C  PERFORMING BI-LINEAR INTERPOLATION USING 4 GRID POINTS FROM A BIGGER
C  RECTANGULAR GRID TO A GRID DESCRIBED BY XLONS AND XLATS OF GRID2.
C  A POINT ON GRID2 IS TRAPPED WITHIN FOUR GRID POINTS ON GRID4.THE
C  GRID POINTS ARE ALWAYS TO THE NORTH AND EAST OF THE TRAPPED POINT.
C  THE ALGORITHM COMPUTES THE FRACTIONAL DISTANCES IN BOTH X AND Y
C  DIRECTION OF THE TRAPPED GRID POINT AND USES THE INFORMATION
C  AS WEIGHTING FACTORS IN THE INTERPOLATION.
C  THERE IS ONE LESS ROW AND COLUMN WHEN THE SCALAR FIELDS ARE
C  INTERPOLATED BECAUSE XLATS AND XLONS ARE NOT DEFINED FOR
C  THE CROSS POINTS IN THE MM4 MODEL.
C
C  B2(NLON4,NLAT4,NVECTS) IS THE INPUT FIELD ON REGULAR LAT/LON GRID.
C  B3(NLON2,NLAT2,NVECTS) IS THE OUTPUT FIELD ON LAMBERT CONFORMAL GRID.
C  G4LON.....LONGITUDE VALUES IN DEGREES OF THE INTERMEDIATE GRID4.
C  G4LAT.....LATITUDE VALUES IN DEGREES OF THE INTERMEDIATE GRID4.
C  P.........EAST-WEST WEIGHTING FACTOR.
C  Q.........NORTH-SOUTH WEIGHTING FACTOR.
C  IP........GRID POINT LOCATION IN EAST-WEST OF TRAPPED GRID POINT.
C  IQ........GRID POINT LOCATION IN NORTH-SOUTH OF TRAPPED GRID POINT.
C
      DIMENSION B3(NLON2,NLAT2,NVECTS),B2(NLON4,NLAT4,NVECTS)
      DIMENSION XLONS(NLON2,NLAT2),XLATS(NLON2,NLAT2)
     A,         G4LON(NLON4),G4LAT(NLAT4)
      DO 120 J=1,NJN
         DO 110 I=1,NIN
            YIND=(((XLATS(I,J)-G4LAT(1))/(G4LAT(NLAT4)-G4LAT(1)))
     A           *FLOAT(NLAT4-1))+1.
            JQ=INT(YIND)
            JQP1=MIN0(JQ+1,NLAT4)
            Q=YIND-JQ
            XIND=(((XLONS(I,J)-G4LON(1))/(G4LON(NLON4)-G4LON(1)))
     A           *FLOAT(NLON4-1))+1.
            IP=INT(XIND)
            IPP1=MIN0(IP+1,NLON4)
            P=XIND-IP
            DO 100 L=1,NVECTS
               B3(I,J,L)=(1.-Q)*((1.-P)*B2(IP,JQ,L)+P*B2(IPP1,JQ,L))+
     A                    Q*(P*B2(IPP1,JQP1,L)+(1.-P)*B2(IP,JQP1,L))
  100       CONTINUE
  110    CONTINUE
  120 CONTINUE
      RETURN
      END
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE BOUNDF (F,C,NI,NJ,NK)
      DIMENSION F(NI,NJ,NK)
C
C  THIS ROUTINE SETS THE BOUNDARY VALUES OF A FIELD EQUAL TO A CONSTANT.
C
      DO 2 K=1,NK
      DO 1 I=1,NI
    1 F(I,NJ,K)=C
      DO 2 J=1,NJ-1
    2 F(NI,J,K)=C
C
      RETURN
      END
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE BSSLZR (BES,N)
      DIMENSION BES(N)
      DIMENSION BZ(50)
C
C  THIS ROUTINE PROVIDES FIRST GUESSES FOR THE ROUTINE GAUAW.
C
      DATA PI/3.14159265358979/
      DATA BZ / 2.4048255577, 5.5200781103, 8.6537279129, 11.7915344391,
     * 14.9309177086, 18.0710639679, 21.2116366299, 24.3524715308, 27.49
     *34791320, 30.6346064684, 33.7758202136, 36.9170983537, 40.05842576
     *46, 43.1997917132, 46.3411883717, 49.4826098974, 52.6240518411, 55
     *.7655107550, 58.9069839261, 62.0484691902, 65.1899648002, 68.33146
     *93299, 71.4729816036, 74.6145006437, 77.7560256304, 80.8975558711,
     * 84.0390907769, 87.1806298436, 90.3221726372, 93.4637187819, 96.60
     *52679510, 99.7468198587, 102.8883742542, 106.0299309165, 109.17148
     *96498, 112.3130502805, 115.4546126537, 118.5961766309, 121.7377420
     *880, 124.8793089132, 128.0208770059, 131.1624462752, 134.304016638
     *3, 137.4455880203, 140.5871603528, 143.7287335737, 146.8703076258,
     * 150.0118824570, 153.1534580192, 156.2950342685/
      NN=N
      IF (.NOT.(N.LE.50))GOTO 23014
      GO TO 12
23014 CONTINUE
      BES(50)=BZ(50)
      DO 23016 J=51,N
      BES(J)=BES(J-1)+PI
23016 CONTINUE
23017 CONTINUE
      NN=49
12    CONTINUE
      DO 23018 J=1,NN
      BES(J)=BZ(J)
23018 CONTINUE
23019 CONTINUE
      RETURN
      END
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE DIVIDG (PHIS,NN)
      DIMENSION PHIS(NN)
      DATA G/9.8/
C
C  THIS ROUTINE REMOVES A FACTOR OF G FROM THE SURFACE GEOPOTENTIAL.
C
      DO 1 N=1,NN
    1 PHIS(N)=PHIS(N)/G
C
      RETURN
      END
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE EPSIL2 (EPSI,LALP,LM)
      DIMENSION EPSI(LALP,LM)
C
C  COMPUTE CONSTANTS USED FOR RECURSIVE CALCULATION OF LEGENDRE
C  POLYNOMIALS.  FROM ROGER DALEY 9/16/83
C  CALCULATES EPSILON(N,M)=SQRT((N**2-M**2)/(4*N**2-1))
C  FOR N=0,LALP-1 AND M=0 TO LM-1
C
      DO 20 M=1,LM
      MS=M-1
      N1=1
      IF(MS.EQ.0) N1=2
      DO 20 N=N1,LALP
      NS=MS+N-1
      FNUM=FLOAT(NS**2-MS**2)
      FDEN=FLOAT(4*NS**2-1)
   20 EPSI(N,M)=SQRT(FNUM/FDEN)
      EPSI(1,1)=0.
C
      RETURN
      END
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE FINDPT (NPOINT,CFIELD,IFIELD,CNAME,NPNT1)
      CHARACTER*8 CFIELD(IFIELD), CNAME
      INTEGER NPOINT(IFIELD,2)
C
C  THIS ROUTINE SETS POINTER NPNT1 FOR FIELD CNAME TO VALUE DETERMINED
C  FROM LIST OF POINTERS AND FIELD NAMES.
C
      I=0
    1 I=I+1
      IF (CFIELD(I).EQ.CNAME) THEN
         NPNT1=NPOINT(I,1)
         GO TO 2
      ENDIF
      IF (I.LT.IFIELD) GO TO 1
      PRINT 3,CNAME
      CALL ABORT ('FIELD NAME NOT FOUND IN FINDPT')
    2 RETURN
    3 FORMAT('0 FIELD ',A8,' NOT FOUND IN LIST')
      END
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      FUNCTION FTNAME (IUNIT)
      ENCODE (4,1,FNAME) IUNIT
      FTNAME=FNAME
    1 FORMAT(*FT*,I2)
      RETURN
      END
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE GAUAW (A,W,WORKE,K)
      DIMENSION A(K),W(K)
      DIMENSION WORKE(1)
C
C  COMPUTE GAUSSIAN LATITUDES AND WEIGHTS.  OUTPUT IS FOR ALTERNATING
C  N., S. LATITUDES IF ROUTINE NSNSNS IS CALLED.
C  A IS ACTUALLY ARRAY OF SINES OF LATITUDES.
C
      EPS=1.E-14
      C=(1.-(2./3.14159265358979)**2)*0.25
      FK=K
      KK=K/2
      CALL BSSLZR(A,KK)
      DO 23000 IS=1,KK
         XZ=COS(A(IS)/SQRT((FK+0.5)**2+C))
         ITER=0
10       PKM2=1.
         PKM1=XZ
         ITER=ITER+1
         IF (.NOT.(ITER.GT.10))GOTO 23002
         GO TO 70
23002    CONTINUE
         DO 23004 N=2,K
            FN=N
            PK=((2.*FN-1.)*XZ*PKM1-(FN-1.)*PKM2)/FN
            PKM2=PKM1
            PKM1=PK
23004    CONTINUE
23005    CONTINUE
         PKM1=PKM2
         PKMRK=(FK*(PKM1-XZ*PK))/(1.-XZ**2)
         SP=PK/PKMRK
         XZ=XZ-SP
         AVSP=ABS(SP)
         IF (.NOT.(AVSP.GT.EPS))GOTO 23006
         GO TO 10
23006    CONTINUE
         A(IS)=XZ
         W(IS)=(2.*(1.-XZ**2))/(FK*PKM1)**2
23000 CONTINUE
23001 CONTINUE
      IF (.NOT.(K.EQ.KK*2))GOTO 23008
      GO TO 50
23008 CONTINUE
      A(KK+1)=0.
      PK=2./FK**2
      DO 23010 N=2,K,2
      FN=N
      PK=PK*FN**2/(FN-1.)**2
23010 CONTINUE
23011 CONTINUE
      W(KK+1)=PK
50    CONTINUE
      DO 23012 N=1,KK
      L=K+1-N
      A(L)=-A(N)
      W(L)=W(N)
23012 CONTINUE
23013 CONTINUE
      CALL NSNSNS(A,WORKE,K)
      CALL NSNSNS(W,WORKE,K)
      RETURN
70    PRINT 100
100   FORMAT(//,5X,14Herror in gauaw )
      STOP
      END
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE GRD2SC (B,PHIS,PS,UV,PMN,HMN
     A,                  COEFS,WORK,MTRUNC,NVECTS,JSIGN,KLAT,NLONP2,NLEV
     B,                  MMAXP1,NMAXPX,NMAXP1,NSPECS,NUV,GAUSWT
     C,                  GLAT,NTRUNC,NVECT1)
      COMMON /COMFFT/ IFAX(13),INC,JUMP,NLON,TRIGS(1)
      DIMENSION COEFS(2,NSPECS,NVECTS)
C
C  THIS SUBROUTINE CALLS ALL THE ROUTINES THAT CARRY OUT THE SUCCESSIVE
C  OPERATIONS NEEDED TO CALCULATE THE CONTIBUTION TO THE SPECTRAL COEF-
C  FICIENTS BY ONE LATITUDE OF DATA.
C
      NLON=NLONP2-2
      NLOND2=NLONP2/2
      NBP=NLONP2*NVECTS
C
      CALL LOGPS1 (PS,NLON)
      CALL DIVIDG (PHIS,NLONP2)
      CALL FFT991 (B,WORK,TRIGS,IFAX,INC,JUMP,NLON,NVECTS,-1)
      CALL SPHER1 (B,COEFS,PMN,MTRUNC,NMAXP1,JSIGN,GAUSWT
     A,  MMAXP1,NMAXPX,NVECT1,NLOND2,NSPECS,1)
      CALL SPHER2 (UV,COEFS(1,1,NUV),PMN,HMN,MTRUNC,JSIGN,GLAT,GAUSWT
     B,  NMAXP1,NMAXPX,MMAXP1,NLEV,NLOND2,NSPECS,2)
C
      RETURN
      END
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE GRIDML (XLONS,XLATS,XMAPS,CORIOL,ZS,XLANDU
     A,                LANDUS,SNOWCV
     B,                CLON,CLAT,DELX,GRDFAC,NLONS,NLATS,IUNIT,CGTYPE)
      CHARACTER CGTYPE*8, CFNAME*8
      DIMENSION XLONS(NLONS,NLATS,2), XLATS(NLONS,NLATS,2)
     A,         XMAPS(NLONS,NLATS,2), CORIOL(NLONS,NLATS)
     B,         ZS(NLONS,NLATS), XLANDU(NLONS,NLATS)
     C,         SNOWCV(NLONS,NLATS), LANDUS(NLONS,NLATS), IER(2)
C
C  THIS SUBROUTINE CALLS ROUTINES TO PRODUCE THE MAP FACTORS
C  IT ALSO READS A FILE OF TOPOGRAPHY AND LANDUSE APPROPRIATE FOR GRID
C  FOR EXPLANATION OF VARIABLES SEE SUBROUTINE MAPRON.
C

      IF (CGTYPE.EQ.'LAMBERT ') THEN  ! LAMBERT CONFORMAL MAP SELECTED
         CALL MAPRON (XLATS,XLONS,XMAPS,CORIOL
     A               ,NLONS,NLATS,CLON,CLAT,DELX,0)
         CALL MAPRON (XLATS(1,1,2),XLONS(1,1,2),XMAPS(1,1,2),CORIOL
     B               ,NLONS,NLATS,CLON,CLAT,DELX,1)
         GRDFAC=.716
      ELSE                          ! UNKNOWN MAP SELECTED
         PRINT 10,CGTYPE
         CALL ABORT ('INVALID OPTION REQUESTED: SEE OUTPUT')
      ENDIF
C
C  READ APPROPRIATE FILE OF TERRAIN AND LANDUSE FOR THIS GRID
C
      READ (IUNIT) CFNAME,LENREC
      READ (IUNIT) ((ZS(J,I),I=1,NLATS),J=1,NLONS)
c kdr
c        DO 541 IR=1,NLATS-1
c        DO 541 JR=1,NLONS-1
c           IF (ZS(JR,IR).GE.100000.) THEN
c              PRINT*,'In GRIDML ZS(',JR,IR,') = ',ZS(JR,IR)
c           END IF
c 541    END DO
c end kdr
      READ (IUNIT) ((LANDUS(J,I),I=1,NLATS),J=1,NLONS)

C
      DO 1 J=1,NLONS
      DO 1 I=1,NLATS
      XLANDU(J,I) = FLOAT(LANDUS(J,I))
      SNOWCV(J,I) = 0.0    ! THIS INDICATES NO SNOW
    1 CONTINUE
      PRINT 11, IUNIT, CFNAME, LENREC
C
   10 FORMAT('0 OPTION REQUESTED =',A8,' UNKNOWN IN SUBROUTINE GRIDML')
   11 FORMAT('0MM4 TERRAIN FILE READ FROM UNIT #',I2,': NAME, LENGTH=',
     A        A8,I5,'.')
      RETURN
      END
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE GRID4 (XLONS,XLATS,NLON2,NLAT2
     A,                G4LON,G4LAT,NLON4,NLAT4)
C
C SUBROUTINE TO CREATE AN INTERMEDIATE LAT/LON GRID THAT WILL CONTAIN
C THE MM4 GRID OF ANY DIMENSION AND AT ANY LOCATION ON THE GLOBE WITH
C THE EXCEPTION OF THE POLES.
C
      DIMENSION XLONS(NLON2,NLAT2),XLATS(NLON2,NLAT2)
     A,         G4LON(NLON4),G4LAT(NLAT4)
      SMLON4=180.
      BGLON4=-180.
      SMLAT4=90.
      BGLAT4=-90.
      DO 10 I=1,NLAT2
      DO 10 J=1,NLON2
      BGLON4=AMAX1(XLONS(J,I),BGLON4)
      SMLON4=AMIN1(XLONS(J,I),SMLON4)
      BGLAT4=AMAX1(XLATS(J,I),BGLAT4)
      SMLAT4=AMIN1(XLATS(J,I),SMLAT4)
   10 CONTINUE
C
C AFTER THE MAXIMUM AND MINIMUM OF BOTH LATITUDE AND LONGITUDE ARE
C DETERMINED AN EXTRA DEGREE IS ADDED TO ALL SIDES OF THE GRID FOR
C INTERPOLATION (GRID4) TO ENSURE THAT THE ORIGINAL GRID (GRID2) LIES
C WITHIN THE NEW GENERATED REGULAR GRID4.
C
      IF(BGLAT4 .GE. 89.9999 .OR. SMLAT4 .LE. -89.9999) GO TO 99
      BGLAT4=BGLAT4+1.
      IF(BGLAT4 .GE. 90.)BGLAT4=89.99
      SMLAT4=SMLAT4-1.
      IF(SMLAT4 .LE. -90.)SMLAT4=-89.99
      BGLON4=BGLON4+1.
      SMLON4=SMLON4-1.
      DX=(BGLON4-SMLON4)/FLOAT(NLON4-1)
      DY=(BGLAT4-SMLAT4)/FLOAT(NLAT4-1)
      DO 20 J=1,NLON4
      G4LON(J)=SMLON4+(FLOAT(J-1)*DX)
   20 CONTINUE
      DO 30 I=1,NLAT4
      G4LAT(I)=SMLAT4+(FLOAT(I-1)*DY)
   30 CONTINUE
      WRITE(6,111)
  111 FORMAT(/42X,' DOMAIN OF THE INTERPOLATED GRID ',)
      WRITE(6,112)
  112 FORMAT(/30X,'CORNERS OF GRID FOR INTERPOLATION (LAT,LON) FOLLOWS:'
     A)
      WRITE(6,113) G4LAT(1),G4LON(1),G4LAT(NLAT4),G4LON(1)
     A,            G4LAT(1),G4LON(NLON4),G4LAT(NLAT4),G4LON(NLON4)
  113 FORMAT(5X/,4(3X,'( ',F7.2,' , ',F7.2,' )'))
      RETURN
   99 PRINT 100
  100 FORMAT(10X,' STOP IN SUBROUTINE GRID4. EITHER THE MAXIMUM
     A LATITUDE IS GREATER THAN 90 DEGREES OR THE MINIMUM IS
     B LESS THAN -90 DEGREES ' )
      STOP 111
      END
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE HUMID1 (T,Q,PS,PT,SIGMA,NI,NJ,NK)
      PARAMETER (TR=1./273.16)
      DIMENSION T(NI,NJ,NK),Q(NI,NJ,NK),PS(NI,NJ),SIGMA(NK)
C
C  THIS ROUTINE REPLACES SPECIFIC HUMIDITY BY RELATIVE HUMIDITY
C
      DO 1 I=1,NI-1
      DO 1 J=1,NJ-1
      DO 1 K=1,NK
      P=(PT+SIGMA(K)*PS(I,J))*10.               ! PRESSURE AT LEVEL K
      HL=597.3-.566*(T(I,J,K)-273.16)           ! LATENT HEAT OF EVAP.
      SATVP=6.11*EXP(9.045*HL*(TR-1./T(I,J,K))) ! SATURATION VAP PRESS.
      QS=.622*SATVP/(P-SATVP)                   ! SAT. MIXING RATIO
      Q(I,J,K)=Q(I,J,K)/QS
    1 CONTINUE
      RETURN
C
      ENTRY HUMID2 (T,Q,PS,PT,SIGMA,NI,NJ,NK)
C
C  THIS ROUTINE REPLACES RELATIVE HUMIDITY BY SPECIFIC HUMIDITY
C
      DO 2 I=1,NI-1
      DO 2 J=1,NJ-1
      DO 2 K=1,NK
      P=(PT+SIGMA(K)*PS(I,J))*10.
C??   P=PT*10.+SIGMA(K)*1000.          ! CHANGED BY GB, APR 90
      HL=597.3-.566*(T(I,J,K)-273.16)
      SATVP=6.11*EXP(9.045*HL*(TR-1./T(I,J,K)))
      IF(P.LT.300.) P = 300.           ! GB MOD: KEEP Q SMALL FOR P<300
      QS=.622*SATVP/(P-SATVP)
      IF (Q(I,J,K).LT.0.1) Q(I,J,K)=0.1
      Q(I,J,K)=Q(I,J,K)*QS
    2 CONTINUE
C
      RETURN
      END
C
C X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE HYDROST (H,T,PHIS,PS,PT,SIGMAF,SIGMAH,DSIGMA
     A,                   NI,NJ,NI1,NJ1,NK)
C
C ROUTINE TO COMPUTE HEIGHT USING THE HYDROSTATIC RELATION.
C THE METHOD UTILIZED HERE IS CONSISTENT WITH THE WAY THE
C HEIGHT IS COMPUTED IN THE MM4 MODEL.
C
      PARAMETER (RGAS=287.04, GRAV=9.8)
      DIMENSION H(NI,NJ,NK),T(NI,NJ,NK),PHIS(NI,NJ),PS(NI,NJ)
      DIMENSION SIGMAF(NK+1),SIGMAH(NK),DSIGMA(NK)
      RG=RGAS/GRAV
      DO 1 K=1,NK
    1 DSIGMA(K)=SIGMAF(K+1)-SIGMAF(K)
C
C SET BOUNDARY VALUES TO ZERO AT ALL LEVELS SINCE THE HEIGHT IS
C DEFINED AT CROSS POINTS AND AT HALF LEVELS.
C
      CALL BOUNDF(H,0.,NI,NJ,NK)
      DO 2 I=1,NI1
      DO 2 J=1,NJ1
      PF=PT/PS(I,J)
      H(I,J,NK)=PHIS(I,J)+RG*T(I,J,NK)*ALOG((1.+PF)/(SIGMAH(NK)+PF))
      DO 2 K2=1,NK-1
      K=NK-K2
      K1=K+1
      TBAR=(T(I,J,K)*DSIGMA(K)+T(I,J,K1)*DSIGMA(K1))/(DSIGMA(K)
     A     +DSIGMA(K1))
      H(I,J,K)=H(I,J,K1)+RG*TBAR*ALOG((SIGMAH(K1)+PF)/
     A         (SIGMAH(K)+PF))
    2 CONTINUE
      RETURN
      END
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE INTPS  (PSMM4,PSCCM,T,ZMM4,ZCCM,PT,NI,NJ)
      PARAMETER (RGAS=287.04, RLAPSE=-6.5E-03, GRAV=9.8)
      DIMENSION PSMM4(NI,NJ), PSCCM(NI,NJ), T(NI,NJ)
     A,          ZMM4(NI,NJ), ZCCM(NI,NJ)
C
C  INTERPOLATE SURFACE PRESSURE. THIS VERTICAL INTERPOLATION IS
C  NECESSARY BECAUSE THE TERRAIN HEIGHTS IN THE CCM AND MM4 DIFFER.
C  THIS INTERPOLATION USES THE HYDROSTATIC RELATION WITH A VERTICAL-
C  MEAN T DETERMINED FROM THE LOWEST LEVEL CCM T AND AN ASSUMED LAPSE
C  RATE RLAPSE.  PSMM4 = SURFACE PRESSURE - PTOP
C
      RL2=RLAPSE/2.
      GDRM=-GRAV/RGAS
      DO 1 I=1,NI-1
      DO 1 J=1,NJ-1
      TB=T(I,J)+RL2*(ZMM4(I,J)-ZCCM(I,J))
      PSMM4(I,J)=PSCCM(I,J)*EXP(GDRM*(ZMM4(I,J)-ZCCM(I,J))/TB)-PT
    1 CONTINUE
      PRINT *, 'PSMM4(45,5) = ', PSMM4(45,5)
 
      RETURN
      END
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE INTV1 (FMM4,FCCM,PSMM4,PSCCM
     A,                 SMM4,SCCM,PT,NI,NJ,KMM4,KCCM,NI1,NJ1
     B,                 IFIELD,KREV,FCCMR,C1,C2)
      PARAMETER (RGAS=287.04, RLAPSE=-6.5E-03, GRAV=9.8 )
      PARAMETER (RGAS2=RGAS/2., B1=GRAV/RLAPSE )

      DIMENSION PSMM4(NI,NJ),PSCCM(NI,NJ),SMM4(KMM4),SCCM(KCCM)
     A,          FMM4(NI,NJ,KMM4,IFIELD), FCCM(NI,NJ,KCCM,IFIELD)
 
C  ARRAYS USED IN REVERSING VERTICAL INDEXES OF CCM DATA.
C  INPUT ARRAY(S) NOT AFFECTED WITH REVERSAL OF INDEXES
      DIMENSION FCCMR(NI,NJ,KCCM,2)
 
C DIMENSION DUMMY ARRAYS TO PRINT FOR CHECK
      DIMENSION C1(KCCM),C2(KMM4)


C  INTV1 IS FOR VERTICAL INTERPOLATION OF U, V, AND RELATIVE HUMIDITY.
C        THE INTERPOLATION IS LINEAR IN P.  WHERE EXTRAPOLATION IS
C        NECESSARY, FIELDS ARE CONSIDERED TO HAVE 0 VERTICAL DERIVATIVE.
C  INTV2 IS FOR VERTICAL INTERPOLATION OF T.  THE INTERPOLATION IS
C        LINEAR IN LOG P.  WHERE EXTRAPOLATION UPWARD IS NECESSARY,
C        THE T FIELD IS CONSIDERED TO HAVE 0 VERTICAL DERIVATIVE.
C        WHERE EXTRAPOLATION DOWNWARD IS NECESSARY, THE T FIELD IS
C        CONSIDERED TO HAVE A LAPSE RATE OF RLAPSE (K/M), AND THE
C        THICKNESS IS DETERMINED HYDROSTATICALLY FROM THE MEAN OF THE
C        TWO EXTREME TEMPERATUES IN THE LAYER.
 
 
      DO 1 I=1,NI1
      DO 1 J=1,NJ1
      PSCCM(I,J)=100.
      DO 2 K=1,KCCM
      KR=KCCM-K+1
      DO 2 IFLD=1,IFIELD
      FCCMR(I,J,KR,IFLD)=FCCM(I,J,K,IFLD)
    2 CONTINUE
    1 CONTINUE
C
      DO 20 I=1,NI1
      DO 20 J=1,NJ1
      DP1=PSMM4(I,J)/PSCCM(I,J)
      PT1=PT/PSCCM(I,J)
      DO 20 N=1,KMM4
      SC=SMM4(N)*DP1+PT1
      K1=0
      DO 10 K=1,KCCM
      IF (SC.GT.SCCM(K)) K1=K
   10 CONTINUE
 
c SCCM looks like pressure levels of data, in bars, even
c though it's read from something called ASIGMA
 
C  CONDITION FOR SC .LT. SCCM(1) FOLLOWS
      IF (K1.EQ.0) THEN
         DO 12 IFLD=1,IFIELD
   12    FMM4(I,J,N,IFLD)=FCCMR(I,J,1,IFLD)
         GO TO 18
      ENDIF
C
C  CONDITION FOR SCCM(1) .LT. SC .LT. SCCM(KCCM) FOLLOWS
      IF (K1.NE.KCCM) THEN
         K1P=K1+1
         RC=(SC-SCCM(K1))/(SCCM(K1)-SCCM(K1P))
         RC1=RC+1.
         DO 14 IFLD=1,IFIELD
   14    FMM4(I,J,N,IFLD)=RC1*FCCMR(I,J,K1,IFLD)-RC*FCCMR(I,J,K1P,IFLD)
         GO TO 18
      ENDIF
C
C  CONDITION FOR SC .GT. SCCM(KCCM) FOLLOWS
         DO 16 IFLD=1,IFIELD
   16    FMM4(I,J,N,IFLD)=FCCMR(I,J,KCCM,IFLD)
C
   18 CONTINUE
   20 CONTINUE
 
      RETURN
 
 
      ENTRY INTV2 (FMM4,FCCM,PSMM4,PSCCM
     A,            SMM4,SCCM,PT,NI,NJ,KMM4,KCCM,IFIELD,FCCMR,C1,C2)
 
 
      DO 101 I=1,NI
      DO 101 J=1,NJ
         PSCCM(I,J) = 100.
         DO 102 K=1,KCCM
            KR=KCCM-K+1
         DO 102 IFLD=1,IFIELD
            FCCMR(I,J,KR,IFLD)=FCCM(I,J,K,IFLD)
            IF(I.EQ.26 .AND. J.EQ.29) THEN
               C1(K)=FCCM(I,J,K,IFLD)
            END IF
  102    CONTINUE
  101 CONTINUE
C      PRINT 39,(C1(K),K=1,KCCM)
C
      DO 40 I=1,NI-1
      DO 40 J=1,NJ-1
         DP1=PSMM4(I,J)/PSCCM(I,J)
         PT1=PT/PSCCM(I,J)
      DO 40 N=1,KMM4
         SC=SMM4(N)*DP1+PT1
         K1=0
         DO 30 K=1,KCCM
            IF (SC.GT.SCCM(K)) K1=K
   30    CONTINUE
 
C  CONDITION FOR SC .LT. SCCM(1) FOLLOWS
         IF (K1.EQ.0) THEN
            FMM4(I,J,N,1)=FCCMR(I,J,1,1)
            GO TO 38
         ENDIF
C
C  CONDITION FOR SCCM(1) .LT. SC .LT. SCCM(KCCM) FOLLOWS
         IF (K1.NE.KCCM) THEN
               K1P=K1+1
            RC=ALOG(SC/SCCM(K1))/ALOG(SCCM(K1)/SCCM(K1P))
            RC1=RC+1.
            FMM4(I,J,N,1)=RC1*FCCMR(I,J,K1,1)-RC*FCCMR(I,J,K1P,1)
 
            IF( (I.EQ.5) .AND. (J.EQ. 5) .AND. (N.EQ.KMM4) ) THEN
               PRINT *, ' RC1, FCCMR(I,J,K1,1), RC, FCCMR(I,J,K1P,1) = '
     A          ,   RC1, FCCMR(I,J,K1,1), RC, FCCMR(I,J,K1P,1)
               PRINT *, 'FMM4(I,J,N,1) = ', FMM4(I,J,N,1)
            ENDIF
 
            GO TO 38
         ENDIF
C
C  CONDITION FOR SC .GT. SCCM(KCCM) FOLLOWS
         A1=RGAS2*ALOG(SC/SCCM(KCCM))
         FMM4(I,J,N,1)=FCCMR(I,J,KCCM,1)*(B1-A1)/(B1+A1)
C
   38    CONTINUE
   40 CONTINUE
      DO 201 K=1,KMM4
       IF(K.LE.KCCM) THEN
       C1(K)=FCCMR(5,5,K,1)
       END IF
  201 C2(K)=FMM4(5,5,K,1)
      PRINT 39,(C1(K),K=1,KCCM)
      PRINT 39,(C2(K),K=1,KMM4)
   39 FORMAT('0 T',10E10.3,'1CCM 2MM4')
 
      RETURN
      END
C
C XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C
      SUBROUTINE INTV3 (FSCCM,FCCM,PSRCCM,SCCM,NI,NJ,KCCM,FCCMR)
c psrccm is the mm4 hor grid, but ccm vert grid
c sccm is all mm4 (and should be in mbar)
 
      PARAMETER (RGAS=287.04, RLAPSE=-6.5E-03, GRAV=9.8 )
      PARAMETER (RGAS2=RGAS/2., B1=GRAV/RLAPSE )

      REAL FSCCM(NI,NJ), FCCM(NI,NJ,KCCM), PSRCCM(NI,NJ), SCCM(KCCM)
 
C** ARRAYS USED IN REVERSING VERTICAL INDEXES OF CCM DATA.
C** INPUT ARRAY(S) NOT AFFECTED WITH REVERSAL OF INDEXES
      REAL FCCMR(NI,NJ,KCCM)
 
C** INTV3 IS FOR VERTICAL INTERPOLATION OF TSCCM.  THE INTERPOLATION IS
C        LINEAR IN LOG P.  WHERE EXTRAPOLATION UPWARD IS NECESSARY,
C        THE T FIELD IS CONSIDERED TO HAVE 0 VERTICAL DERIVATIVE.
C        WHERE EXTRAPOLATION DOWNWARD IS NECESSARY, THE T FIELD IS
C        CONSIDERED TO HAVE A LAPSE RATE OF RLAPSE (K/M), AND THE
C        THICKNESS IS DETERMINED HYDROSTATICALLY FROM THE MEAN OF THE
C        TWO EXTREME TEMPERATUES IN THE LAYER.
C
      DO 1 I=1,NI
      DO 1 J=1,NJ
      DO 1 K=1,KCCM
      KR=KCCM-K+1
      FCCMR(I,J,KR)=FCCM(I,J,K)
    1 CONTINUE
C
C
      DO 40 I=1,NI-1
      DO 40 J=1,NJ-1
      SC=PSRCCM(I,J)/100.
      DO 30 K=1,KCCM-1
         IF (SC.LE.SCCM(K+1) .AND. SC.GE.SCCM(K)) K1=K
   30 CONTINUE
 
C  CONDITION FOR SC .GT. SCCM(KCCM) FOLLOWS
         IF(SC.GT.SCCM(KCCM)) THEN
         A1=RGAS2*ALOG(SC/SCCM(KCCM))
         FSCCM(I,J)=FCCMR(I,J,KCCM)*(B1-A1)/(B1+A1)
         GO TO 38
         END IF
C
C  CONDITION FOR SC .LT. SCCM(KCCM) FOLLOWS
         K1P=K1+1
         RC=ALOG(SC/SCCM(K1))/ALOG(SCCM(K1)/SCCM(K1P))
         RC1=RC+1.
         FSCCM(I,J)=RC1*FCCMR(I,J,K1)-RC*FCCMR(I,J,K1P)
C
   38 CONTINUE
   40 CONTINUE
 
      RETURN
      END
 
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE LOGPS1 (PS,NLON)
      PARAMETER (PBAR=1.E5)
      DIMENSION PS(NLON)
C
C  THIS ROUTINE REPLACES THE SURFACE PRESSURE AT GRID POINTS AT ONE
C  LATITUDE WITH THE NATURAL LOGARITHM OF THE SURFACE PRESSURE.  THE
C  LOG OF A MEAN PRESSURE FIELD PBAR IS REMOVED.
C
      DO 1 N=1,NLON
    1 PS(N)=ALOG(PS(N)/PBAR)
      RETURN
      END
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE MAPRON (XLAT,XLON,XMAP,CORIOL,IX,JX,CLON,CLAT
     A                  ,DELX,ITYPE)
      DIMENSION XLAT(IX,JX), XLON(IX,JX) ,XMAP(IX,JX), CORIOL(IX,JX)
C
C  COMPUTE LATS, LONS, MAP-SCALE FACTORS, AND CORIOLIS PARAMETER FOR
C  LAMBERT CONFORMAL MAP CENTERED AT CLON,CLAT. TRUE AT 30.N AND 60.N.
C
C  IX IS NUMBER OF E-W POINTS.  JX IS NUMBER OF N-S POINTS.
C  CLON, CLAT IS LAT, LON OF CENTER OF GRID (DEGREES EAST, NORTH).
C  DELX IS GRID SPACING IN METERS.
C  ITYPE=1 FOR DOT GRID, 0 FOR CROSS GRID.
C  CORIOLIS FACTOR DEFINED ON DOT GRID ONLY.
C
      XOMEGA=7.2722E-5              ! ANG. ROT OF EARTH IN S**-1
      GRDFAC=.716                   ! FACTOR FOR LAMBERT CONFORMAL
      EARTHR=6.37E6                 ! METERS
      PSI1=30.                      ! LAT AT WHICH PROJ. IS TRUE
C                                   !   TRUE AT COLAT OF PSI1 TOO
      PIFAC=ATAN(1.)/45.            ! CONVERT DEGREES TO RADIANS
      PSI0=90.-CLAT                 ! COLAT OF CENTER
C
      PSI1=PIFAC*PSI1               ! CONVERT TO RADIANS
      PSI0=PIFAC*PSI0               ! CONVERT TO RADIANS
      RGRDF=1./GRDFAC               ! RECIPROCAL OF GRID FACTOR
      ANS=EARTHR*SIN(PSI1)*RGRDF    ! A*SIN(PS1)/N
      C1=-CLON-90.*RGRDF
      C2=ANS*(TAN(PSI0*0.5)/TAN(PSI1*0.5))**GRDFAC
      PF=TAN(PSI1*0.5)/(ANS**RGRDF)
      XMF=SIN(PSI1)/(TAN(0.5*PSI1)**GRDFAC)
      PI90=2.*ATAN(1.)              ! PI OVER TWO (90 DEGREES)
C
      XC=0.5*FLOAT(IX+ITYPE)
      YC=0.5*FLOAT(JX+ITYPE)
      DO 10 I=1,IX
      X=(I-XC)*DELX
      DO 10 J=1,JX
      Y=(J-YC)*DELX
      IF (X.NE.0.) THEN
         XLP=ATAN2(Y-C2,X)
         R=ABS(X/COS(XLP))
         ELSE
         XLP=-PI90
         R=ABS(Y-C2)
      ENDIF
      PSI=ABS(2.*ATAN(PF*(R**RGRDF)))        ! CO LATITUDE (RADIANS)
      PFI=PI90-PSI
      XLON(I,J)=XLP*RGRDF/PIFAC-C1           ! LONGITUDE DEGREES EAST
      XLAT(I,J)=PFI/PIFAC                    ! LATITUDE DEGREES NORTH
      XMAP(I,J)=(XMF/SIN(PSI))*(TAN(0.5*PSI)**GRDFAC)
   10 CONTINUE
C
      IF (ITYPE.EQ.0) RETURN
      XOMEG2=2.*XOMEGA
      DO 20 I=1,IX
      DO 20 J=1,JX
      CORIOL(I,J)=XOMEG2*SIN(XLAT(I,J)*PIFAC)
   20 CONTINUE
C
      RETURN
      END
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE MTRA (MTRUNC,NTRUNC,MMAX,NMAX,KMAX)
      DIMENSION MTRUNC(NMAX+1),NTRUNC(MMAX+1)
C
C  THIS ROUTINE CALCULATES THE TWO VECTOR POINTERS MTRUNC AND NTRUNC
C     MTRUNC IS THE LENGTH OF A DIAGONAL IN N,M SPACE FOR GIVEN N-M+1
C     NTRUNC IS THE LENGTH OF A COLUMN   IN N,M SPACE FOR GIVEN M+1
C  VALUES CONSISTENT WITH A GENERAL PENTAGONAL SPECTRAL TRUNCATION.
C
      KSUBM=KMAX-MMAX+1
      KSUBN=KMAX-NMAX+1
C
      DO 1 N=1,NMAX+1
      MTRUNC(N)=MMAX+1
      IF(N.GT.KSUBM) MTRUNC(N)=MTRUNC(N-1)-1
    1 CONTINUE
C
      DO 2 M=1,MMAX+1
      NTRUNC(M)=NMAX+1
      IF(M.GT.KSUBN) NTRUNC(M)=NTRUNC(M-1)-1
    2 CONTINUE
C
      RETURN
      END
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE MULTP (A,P,NI,NJ,NK,NI1,NJ1)
      DIMENSION A(NI,NJ,NK),P(NI,NJ)
C
C  THIS ROUTINE REPLACES A BY A * P , WHERE  A  IS  U, V, T, OR Q.
C
      DO 1 K=1,NK
      DO 1 J=1,NJ1
      DO 1 I=1,NI1
    1 A(I,J,K)=A(I,J,K)*P(I,J)
C
      RETURN
      END
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE NEWPNT (NPOINT,OPOINT,CFIELD,IFIELD,NLONP2,IFUNCT)
      CHARACTER*8 CFIELD(IFIELD)
      INTEGER NPOINT(IFIELD,2), OPOINT(IFIELD,3)
C
C  THIS ROUTINE DETERMINES THE POINTERS FOR THE RE-ORDERD FIELDS.
C
      NA=0
      DO 1 I=1,IFIELD
      NLEV  = OPOINT(I,2)
      NPOINT(I,1)=NA+1
      IF (IFUNCT.EQ.0) THEN   ! POINTERS FOR SINGLE OR MULTILEVEL FIELDS
         NSKIP=NLONP2*NLEV
         NPOINT(I,2)=NLEV
      ELSE                    ! POINTERS FOR SINGLE LEVELS FIELDS ONLY
         NSKIP=NLONP2
         NPOINT(I,2)=-NPOINT(I,2)
      ENDIF
      NA=NA+NSKIP
    1 CONTINUE
C
C  PRINT TABLE OF NEW AND OLD VALUES
C
      PRINT 3
      DO 2 I=1,IFIELD
      PRINT 4,I,CFIELD(I),(OPOINT(I,J),J=1,3),(NPOINT(I,J),J=1,2)
    2 CONTINUE
C
    3 FORMAT('0',50X,'TABLE OF POINTERS',/,1X,' I',2X
     A,'   FIELD','  OLD PT','  LEVS','  DENS',3X,'NEW PT','  LEVS')
    4 FORMAT(1X,I2,2X,A8,I8,2I6,3X,2I6)
      RETURN
      END
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE NSNSNS (G,W,NN)
      DIMENSION G(NN),W(NN)
C
C  RE-ORDER SO THAT SOUTH TO NORTH
C
      DO 1 N=1,NN,1
    1 W(N)=G(NN+1-N)
      DO 2 N=1,NN
    2 G(N)=W(N)
C
C     N1=0
C     DO 1 N=1,NN,2
C     N1=1+N1
C     W(N)=G(N1)
C   1 W(N+1)=G(NN-N1+1)
C
C     DO 2 N=1,NN
C   2 G(N)=W(N)
C
      RETURN
      END
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE NZONALY (D,BM,G4LON,NLON4,MMAXP1,NVECTS)
C
C THIS SUBROUTINE PERFORMS AN INVERSE FOURIER TRANSFORM.
C IT IS CALLED FOR EVERY LATITUDE OVER ALL LONGITUDES AND FOR ALL
C FIELDS.
C
      DIMENSION D(NLON4,NVECTS),BM(2,MMAXP1,NVECTS),G4LON(NLON4)
      PIR180=ATAN(1.)/45.
      DO 1 M=1,MMAXP1
      M1=M-1
      XF=2.
      IF(M .EQ. 1) XF=1.
      DO 2 I=1,NLON4
      CS=COS(M1*G4LON(I)*PIR180)*XF
      SI=SIN(M1*G4LON(I)*PIR180)*XF
      DO 2 L=1,NVECTS
      D(I,L)=D(I,L)+CS*BM(1,M,L)-SI*BM(2,M,L)
    2 CONTINUE
    1 CONTINUE
      RETURN
      END
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE P1P2 (PD,PX,NI,NJ)
      DIMENSION PD(NI,NJ),PX(NI,NJ)
C
C  THIS ROUTINE DETERMINES P(.) FROM P(X) BY A 4-POINT INTERPOLATION.
C  ON THE X-GRID, A P(X) POINT OUTSIDE THE GRID DOMAIN IS ASSUMED TO
C  SATISFY P(0,J)=P(1,J); P(NI,J)=P(NI-1,J); AND SIMILARLY FOR THE I'S.
C
      NI1=NI-1
      NJ1=NJ-1
C
      DO 1 J=2,NJ1
      DO 1 I=2,NI1
    1 PD(I,J)=0.25*(PX(I,J)+PX(I-1,J)+PX(I,J-1)+PX(I-1,J-1))
C
      DO 2 I=2,NI1
      PD(I,1)=0.5*(PX(I,1)+PX(I-1,1))
    2 PD(I,NJ)=0.5*(PX(I,NJ1)+PX(I-1,NJ1))
C
      DO 3 J=2,NJ1
      PD(1,J)=0.5*(PX(1,J)+PX(1,J-1))
    3 PD(NI,J)=0.5*(PX(NI1,J)+PX(NI1,J-1))
C
      PD(1,1)=PX(1,1)
      PD(1,NJ)=PX(1,NJ1)
      PD(NI,1)=PX(NI1,1)
      PD(NI,NJ)=PX(NI1,NJ1)
C
      RETURN
      END
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE PMN2X (PMN,NSPECTS)
      DIMENSION PMN(NSPECTS)
      RT2=SQRT(2.)
      DO 1 N=1,NSPECTS
    1 PMN(N)=RT2*PMN(N)
      RETURN
      END
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE PRINTB (B,NN,NX,AN)
      CHARACTER*8 AN
      DIMENSION B(1)
      PRINT 1,AN,(B(N),N=1,NN,NX)
    1 FORMAT('0 ',A8,1P6E18.7,5(/,10X,6E18.7))
      RETURN
      END
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE RDATES (NDATE,NTIME,NSTEP)
      COMMON /CDATEI/ IDATE
C
C  UNPACK DATE AND TIME AND PRINT IN APPROPRIATE FORMAT
C
      IDATE=NDATE*100+NTIME/3600
      IYEAR=NDATE/10000
      IMONTH=(NDATE-IYEAR*10000)/100
      IDAY=NDATE-IYEAR*10000-IMONTH*100
      IHOUR=NTIME/3600
      IHOUR1=IHOUR/10
      IHOUR2=IHOUR-IHOUR1*10
      IMINUT=(NTIME-IHOUR*3600)/60
      IMIN1=IMINUT/10
      IMIN2=IMINUT-IMIN1*10
C
      PRINT 10,IYEAR,IMONTH,IDAY
      PRINT 11,IHOUR1,IHOUR2,IMIN1,IMIN2
      PRINT 12,NSTEP
   10 FORMAT('0 DATE  OF ANALYZED DATA  YR/MO/DA = ',I2,'/',I2,'/',I2)
   11 FORMAT('0 TIME  OF ANALYZED DATA           = ',4I1,'Z',//)
   12 FORMAT('0 NSTEP OF ANALYZED DATA           = ',I8)
C
      RETURN
      END
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE READGD (XLONS,XLATS,WORK,ZS,XMAP,DMAP,CORIOL
     A,                  XLANDU,SNOW,NLONS,NLATS,IUNIT)
      DIMENSION XLONS(NLONS,NLATS,2), XLATS(NLONS,NLATS,2)
     A,         WORK(NLATS,4), ZS(NLONS,NLATS), XMAP(NLONS,NLATS)
     B,         DMAP(NLONS,NLATS), CORIOL(NLONS,NLATS)
     C,         XLANDU(NLONS,NLATS), SNOW(NLONS,NLATS), IER(2)
C
C  THIS ROUTINE READS THE FILE FILGRD ON IUNIT AND FILLS THE ARRAYS
C  XLONS AND XLATS.
C  THE THREE INDICES OF XLONS AND XLATS ARE RESPECTIVELY : J,I, AND
C  AN INDEX WHICH DENOTES THE CROSS (=1) OR DOT (=2) GRIDS.
C
      PRINT 79,IUNIT
   79 FORMAT('0MM4 GRID INPUT ON UNIT',I3)
C
      READ (IUNIT)
      DO 4 J=1,NLONS
      CALL SKIPR (IUNIT,6,IER)
      READ (IUNIT) (CORIOL(J,I),I=1,NLATS)
      READ (IUNIT) (XMAP(J,I),I=1,NLATS)
      READ (IUNIT) (DMAP(J,I),I=1,NLATS)
      READ (IUNIT) (WORK(I,1),I=1,NLATS)  ! READ LATS ON X-GRID
      READ (IUNIT) (WORK(I,2),I=1,NLATS)  ! READ LONS ON X-GRID
      READ (IUNIT) (WORK(I,3),I=1,NLATS)  ! READ LATS ON .-GRID
      READ (IUNIT) (WORK(I,4),I=1,NLATS)  ! READ LONS ON .-GRID
      CALL SKIPR (IUNIT,1,IER)
      READ (IUNIT) (ZS(J,I),I=1,NLATS)    ! READ SURFACE HEIGHT
      READ (IUNIT) (XLANDU(J,I),I=1,NLATS)
      READ (IUNIT) (SNOW(J,I),I=1,NLATS)
C
      IF (J.EQ.NLONS) GO TO 2
      DO 1 I=1,NLATS-1
      XLONS(J,I,1)=WORK(I,2)
    1 XLATS(J,I,1)=WORK(I,1)
    2 CONTINUE
      DO 3 I=1,NLATS
      XLONS(J,I,2)=WORK(I,4)
    3 XLATS(J,I,2)=WORK(I,3)
C
    4 CONTINUE
C
      RETURN
      END
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE READHD (IH,XH,CH,SIGMA,CFIELD
     A,                  NPOINT,NLONW,NBSIZE,NDIM,NDIMB,IFIELD
     B,                  NLON,NLAT,NLEV,MMAX,NMAX,KMAX,IUNIT)
      PARAMETER (LSIZE= 6, LNLON= 9, LLONW=10, LNLAT=11
     A,          LNLEV=12, LMMAX=13, LNMAX=14, LKMAX=15
     B,          LNFLD=16, LSTEP=18, LDATE=26, LTIME=27
     C,          LLHDC=31, LLHDR=32, LSIGS=33, LGAUS=34
     D,          LGAUW=35, LFLDS=36, LCFLD=37, LVNAM=64 )
      DIMENSION IH(NDIM), XH(NDIM), SIGMA(NLEV), NPOINT(IFIELD,3)
      CHARACTER*8 CH(NDIM), CFIELD(IFIELD)
C
C  THIS ROUTINE READ A CCM1 HEADER.
C  OUTPUT ARE VALUES FOR NPOINT, SIGMA, NLONW, NBSIZE.
C  VALUES FOR NLON,NLEV,NLAT,MMAX,NMAX,KMAX,NDIMB ARE CHECKED.
C
C
C  DETERMINE HEADER TYPE
C
      READ (IUNIT) LENHD,MTYP10
      MTYPE=MOD(MTYP10,10)
C      IF (MTYPE.NE.2) THEN
         PRINT 100,MTYPE,MTYP10
C         CALL ABORT ('HEADER NOT IN POSSIBLE FORMAT: SEE OUTPUT')
C      ENDIF
C
C  CHECK HEADER RECORD LENGTHS AND READ ALL 3 HEADER RECORDS
C
      IF (LENHD.GT.NDIM) THEN
         PRINT 101,LENHD,NDIM
         CALL ABORT('INSUFFICIENT DIMENSION FOR HEADER RECORDS')
      ENDIF
      BACKSPACE IUNIT
      READ (IUNIT) (IH(J),J=1,LENHD)
      IF ((IH(LLHDC).GT.NDIM).OR.(IH(LLHDR).GT.NDIM)) THEN
         PRINT 102,IH(LLHDC),IH(LLHDR),NDIM
         CALL ABORT('INSUFFICIENT DIMENSION FOR HEADER RECORDS')
      ENDIF
      READ (IUNIT) (CH(J),J=1,IH(LLHDC))
      READ (IUNIT) (XH(J),J=1,IH(LLHDR))
C
C  CHECK THAT TRUNCATION PARAMETERS AGREE.
C
      IERROR=0
C      IF(MMAX.NE.IH(LMMAX))     IERROR=IERROR+1
C      IF(NMAX.NE.IH(LNMAX))     IERROR=IERROR+1
C      IF(KMAX.NE.IH(LKMAX))     IERROR=IERROR+1
      IF(NLEV.NE.IH(LNLEV))     IERROR=IERROR+1
      IF(NLON.NE.IH(LNLON))     IERROR=IERROR+1
      IF(NLAT.NE.IH(LNLAT))     IERROR=IERROR+1
      IF(IERROR.NE.0) THEN
         PRINT 103,MMAX,NMAX,KMAX,NLEV,NLON,NLAT
         PRINT 104,IH(LMMAX),IH(LNMAX),IH(LKMAX),IH(LNLEV)
     A,            IH(LNLON),IH(LNLAT)
         CALL ABORT('TRUNCATION PARAMETERS DISAGREE: SEE OUTPUT')
      ENDIF
C
C  DETERMINE AND PRINT FILE DATE
C
      CALL RDATES (IH(LDATE),IH(LTIME),IH(LSTEP))
C
C  DETERMINE AND PRINT FILE NAME
C
      PRINT 105,(CH(I),I=LVNAM,LVNAM+9)
C
C  COPY SIGMA AT DATA LEVELS
C
      N1=IH(LSIGS)+1
      DO 2 N=1,NLEV
      SIGMA(N)=XH(N1)
      N1=N1+2
    2 CONTINUE
      PRINT 106,(SIGMA(N),N=1,NLEV)
C
C  DETERMINE FIELD POINTERS
C
      NLONW= IH(LLONW)
      MPF=   IH(LFLDS)
      NPF=   IH(LNFLD)
      MPCFLD=IH(LCFLD)
      IDF=3
      CALL RPOINT(MTYPE,IH(MPF),CH(MPCFLD),IDF,NPF,IFIELD,NPOINT,
     A                 CFIELD,NLON,NLEV)
C
C  CHECK DIMENSION SIZE OF B1 ARRAY
C
      IMAX=0
      DO 3 N=1,IFIELD
      IF (IMAX.LT.NPOINT(N,1)) THEN
         IMAX=NPOINT(N,1)
         M=N
      ENDIF
    3 CONTINUE
      KSIZE=NLONW
      NDENS=NPOINT(M,3)
      IF (NDENS.NE.1) KSIZE=2+(NLON+NDENS-1)/NDENS
      NBSIZE=IMAX-1+KSIZE*NPOINT(M,2)
      IF(NBSIZE .GT. NDIMB) THEN
         PRINT 107,NBSIZE,NDIMB
         CALL ABORT (' B1 DIMENSION TOO SMALL: SEE OUTPUT ')
      ENDIF
C
      RETURN
C
  100 FORMAT(' X X X X X X MTYPE, MTYP10 =',2I5,' X X X X X X ')
  101 FORMAT(' X X X X X X LENHD, NDIM =',2I6,' X X X X X X ')
  102 FORMAT(' X X X X X X LENHC, LENHR, NDIM =',3I6,' X X X X X X ')
  103 FORMAT('0   PROGRAM VALUES OF MMAX,NMAX,KMAX,NLEV,NLON,NLAT ARE  '
     A,6I6)
  104 FORMAT('0  CCM FILE VALUES OF MMAX,NMAX,KMAX,NLEV,NLON,NLAT ARE  '
     A,6I6)
  105 FORMAT('0 ANALYSIS DATA FILE NAME ON HEADER = ',10A8)
  106 FORMAT('0(FULL) SIGMA LEVELS ON HEADR = ',15F6.3)
  107 FORMAT(' X X X X X X SIZE OF LATITUDE RECORDS=',I7,' IS TOO LONG',
     A' FOR B1 AS DIMENSIONED, NDIMB =',I7,' X X X X X X ')
      END
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE READHD1 (IH,XH,CH,SIGMA,CFIELD
     A,                  NPOINT,NLONW,NBSIZE,NDIM,NDIMB,IFIELD
     B,                  NLON,NLAT,NLEV,MMAX,NMAX,KMAX,IUNIT)
      PARAMETER (LSIZE= 6, LNLON=10, LLONW=11, LNLAT=12
     A,          LNLEV=13, LMMAX=14, LNMAX=15, LKMAX=16
     B,          LNFLD=17, LSTEP=18, LDATE=27, LTIME=28
     C,          LLHDC=31, LLHDR=32, LSIGS=64, LGAUS=65
     D,          LGAUW=66, LFLDS=67, LCFLD=37, LVNAM=44 )
      DIMENSION IH(NDIM*2), XH(NDIM*2), SIGMA(NLEV), NPOINT(IFIELD,3)
      CHARACTER*8 CH(NDIM), CFIELD(IFIELD)
C
C  THIS ROUTINE READ A CCM1 HEADER.
C  OUTPUT ARE VALUES FOR NPOINT, SIGMA, NLONW, NBSIZE.
C  VALUES FOR NLON,NLEV,NLAT,MMAX,NMAX,KMAX,NDIMB ARE CHECKED.
C
C
C  DETERMINE HEADER TYPE
C
      READ (IUNIT) LENHD,MTYP10
      MTYPE=MOD(MTYP10,10)
      IF (MTYPE.NE.1) THEN
         PRINT 100,MTYPE,MTYP10
         CALL ABORT ('HEADER NOT IN POSSIBLE FORMAT: SEE OUTPUT')
      ENDIF
C
C  CHECK HEADER RECORD LENGTHS AND READ ALL 3 HEADER RECORDS
C
      IF (LENHD.GT.NDIM*2) THEN
         PRINT 101,LENHD,NDIM
         CALL ABORT('INSUFFICIENT DIMENSION FOR HEADER RECORDS')
      ENDIF
      BACKSPACE IUNIT
      READ (IUNIT) (IH(J),J=1,LENHD)
C
C  CHECK THAT TRUNCATION PARAMETERS AGREE.
C
      IERROR=0
C      IF(MMAX.NE.IH(LMMAX))     IERROR=IERROR+1
C      IF(NMAX.NE.IH(LNMAX))     IERROR=IERROR+1
C      IF(KMAX.NE.IH(LKMAX))     IERROR=IERROR+1
      IF(NLEV.NE.IH(LNLEV))     IERROR=IERROR+1
      IF(NLON.NE.IH(LNLON))     IERROR=IERROR+1
      IF(NLAT.NE.IH(LNLAT))     IERROR=IERROR+1
      IF(IERROR.NE.0) THEN
         PRINT 103,MMAX,NMAX,KMAX,NLEV,NLON,NLAT
         PRINT 104,IH(LMMAX),IH(LNMAX),IH(LKMAX),IH(LNLEV)
     A,            IH(LNLON),IH(LNLAT)
         CALL ABORT('TRUNCATION PARAMETERS DISAGREE: SEE OUTPUT')
      ENDIF
C
C  DETERMINE AND PRINT FILE DATE
C
      CALL RDATES (IH(LDATE),IH(LTIME),IH(LSTEP))
C
C  DETERMINE AND PRINT FILE NAME
C
      PRINT 105,(CH(I),I=LVNAM,LVNAM)
C
C  COPY SIGMA AT DATA LEVELS
C
C?      N1=IH(LSIGS)+1
      N1=IH(LSIGS)
      DO 2 N=1,NLEV
         SIGMA(N)=XH(N1)
C?         N1=N1+2
         N1=N1+1
    2 CONTINUE
      PRINT 106,(SIGMA(N),N=1,NLEV)
C
C  DETERMINE FIELD POINTERS
C
      NLONW= IH(LLONW)
      MPF=   IH(LFLDS)
      NPF=   IH(LNFLD)
      MPCFLD=IH(LCFLD)
      IDF=5
      CALL RPOINT(MTYPE,IH(MPF),CH,IDF,NPF,IFIELD,NPOINT,
     A                 CFIELD,NLON,NLEV)
C
C  CHECK DIMENSION SIZE OF B1 ARRAY
C
      IMAX=0
      DO 3 N=1,IFIELD
      IF (IMAX.LT.NPOINT(N,1)) THEN
         IMAX=NPOINT(N,1)
         M=N
      ENDIF
    3 CONTINUE
      KSIZE=NLONW
      NDENS=NPOINT(M,3)
      IF (NDENS.NE.1) KSIZE=2+(NLON+NDENS-1)/NDENS
      NBSIZE=IMAX-1+KSIZE*NPOINT(M,2)
      IF(NBSIZE .GT. NDIMB) THEN
         PRINT 107,NBSIZE,NDIMB
         CALL ABORT (' B1 DIMENSION TOO SMALL: SEE OUTPUT ')
      ENDIF
C
      RETURN
C
  100 FORMAT(' X X X X X X MTYPE, MTYP10 =',2I5,' X X X X X X ')
  101 FORMAT(' X X X X X X LENHD, NDIM =',2I6,' X X X X X X ')
  102 FORMAT(' X X X X X X LENHC, LENHR, NDIM =',3I6,' X X X X X X ')
  103 FORMAT('0   PROGRAM VALUES OF MMAX,NMAX,KMAX,NLEV,NLON,NLAT ARE  '
     A,6I6)
  104 FORMAT('0  CCM FILE VALUES OF MMAX,NMAX,KMAX,NLEV,NLON,NLAT ARE  '
     A,6I6)
  105 FORMAT('0 ANALYSIS DATA FILE NAME ON HEADER = ',10A8)
  106 FORMAT('0(FULL) SIGMA LEVELS ON HEADR = ',10F12.3)
  107 FORMAT(' X X X X X X SIZE OF LATITUDE RECORDS=',I7,' IS TOO LONG',
     A' FOR B1 AS DIMENSIONED, NDIMB =',I7,' X X X X X X ')
      END
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE READPT (OPOINT,NPOINT,CFIELD,IFIELD,NLONP2)
      CHARACTER*8 CFIELD(IFIELD)
      INTEGER NPOINT(IFIELD,2), OPOINT(IFIELD,3)
C
C  THIS ROUTINE READS A LIST OF FIELD NAMES AND LEVELS FROM THE INPUT
C  FILE AND ESTABLISHES POINTERS.
C
      IERROR=0
      DO 1 I=1,IFIELD
      READ (5,2) CFIELD(I), NPOINT(I,2)
      IF (OPOINT(I,2).GT.NPOINT(I,2)) THEN
         IERROR=IERROR+1
         PRINT 3,IERROR,NPOINT(I,2),CFIELD(I),OPOINT(I,2)
      ENDIF
    1 CONTINUE
      CALL NEWPNT (NPOINT,OPOINT,CFIELD,IFIELD,NLONP2,1)
C
      IF (IERROR.NE.0) CALL ABORT ('INPUT FIELD LEVELS WRONG')
      RETURN
    2 FORMAT(A8,I5)
    3 FORMAT('0 ERROR NO.',I3,5X,'REQUESTED LEVEL =',I3,' FOR FIELD'
     A,A8,' NOT LESS THAN NUMBER OF FIELD LEVELS =',I3)
      END
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE REFORM (BOLD,BNEW,OPOINT,NPOINT
     A,                  NBOLD,NBNEW,NLONP2,NLONW,IFIELD)
      INTEGER NPOINT(IFIELD,2), OPOINT(IFIELD,3)
      DIMENSION BOLD(NBOLD), BNEW(NBNEW)
C
C  THIS COPIES DESIRED DATA INTO A NEW BUFFER WITH A NEW ORDER.
C
C  LOOP OVER THE NUMBER OF DESIRED FIELDS
C
      DO 13 I=1,IFIELD
      NDENS = OPOINT(I,3)
      NLEV  = NPOINT(I,2)
      N1    = OPOINT(I,1)-1
      NA    = NPOINT(I,1)-1
C
      NSKIP=NLONW
      IF (NDENS.NE.1) NSKIP=2+(NLONW+NDENS-1)/NDENS
      IF (NLEV.LT.0) THEN    ! PROCESS ONLY 1 LEVEL OF MULTILEVEL FIELD
         N1=N1+NSKIP*(-(1+NLEV))
         NLEV=1
      ENDIF
C
      DO 12 L=1,NLEV
C
      IF (NDENS.EQ.1) THEN         ! DATA IN UNPACKED FORMAT
         DO 10 N=1,NLONW
         BNEW(NA+N) = BOLD(N1+N)
   10    CONTINUE
      ELSE                         ! DATA IN PACKED FORMAT
         CALL UNPKA (BNEW(NA+1), NLONW, BOLD(N1+1), NDENS )

      ENDIF
C
      IF( NLONW .LT. NLONP2) THEN
         DO 11 N=NLONW+1,NLONP2
         BNEW(NA+N) = 0.
   11    CONTINUE
      ENDIF
      N1 = N1+NSKIP
      NA = NA+NLONP2
   12 CONTINUE             ! END LOOP OVER LEVELS FOR 1 FIELD
C
   13 CONTINUE             ! END LOOP OVER FIELDS
C
      RETURN
      END
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE RPOINT (MTYPE,IH,MCFLD,IDF,II,IFIELD,NPOINT,CFIELD,
     A                   NLON,NLEV)
      INTEGER IH(IDF,II), MTYPE, NPOINT(IFIELD,3)
      CHARACTER*8 MCFLD(2,II),  CFIELD(IFIELD)
C
C  THIS ROUTINE SEARCHES THE ARRAY MCFLD FOR THE NAMES STORED IN THE
C  ARRAY CFIELD.  THE ARRAY IH THEN CONTAINS THE CORRESPONDING POINTER
C  INFORMATION.  MCFLD AND IH ARE LISTS OBTAINED FROM THE CCM HEADER.
C  CFIELD IS THE LIST OF NAMES OF FIELDS DESIRED BY THE USER.  ON       T
C  OUTPUT, THE FIRST VALUE OF NPOINT CONTAINS THE POINTER FOR THAT FIELD
C  WITHIN A CCM DATA BUFFER.  THE SECOND VALUE CONTAINS THE NUMBER OF
C  BUFFER.  THE SECOND VALUE WILL CONTAIN THE NUMBER OF LEVELS (1 OR
C  LEVELS (1 OR NLEV) ON WHICH THAT FIELD IS DEFINED. THE THIRD INDEX
C  INDICATES THE DATA PACKING DENSITY.
C
C  KADD IS A POINTER TO THE FIRST FIELD VALUE FOR PACKED DATA OR THE
C  ADDRESS FOR UNPACKED DATA. KSIZ IS THE NUMBER OF PACKED WORDS
C  CORRESPONDING TO NLONW UNPACKED WORDS. 'UNPKA' IS A SYSTEM ROUTINE.
C
C  SET INDICES FOR FORMATS OF PARTICULAR CCM HEADER
      IF (MTYPE .EQ. 1) THEN        ! CCM0B FORMAT
         IPLEV = 2
         IPADR = 3
         IPPAK = 4
         KADD1 = 7      ! LOCATION OF FIRST WORD ON CCM DATA RECORD
         DO 1 I=1,II    ! CHANGE HOLLARITH TO CHARACTER
         WRITE (MCFLD(1,I),'(A8)') IH(1,I)
    1    CONTINUE
      ELSEIF (MTYPE .EQ. 2) THEN    ! CCM1 FORMAT
         IPLEV = 1
         IPADR = 2
         IPPAK = 3
         KADD1 = 1  ! LOCATION OF FIRST WORD ON CCM DATA RECORD
      ELSEIF (MTYPE .GT. 2) THEN    ! UNKNOWN HEADER TYPE
         PRINT 10,MTYPE
         CALL ABORT ('UNKNOWN HEADER TYPE IN SUB. RPOINT')
      ENDIF
C
C  BEGIN LOOP OVER REQUESTED FIELD NAMES
      DO 777 I=1,II
  777  PRINT 778,I,(IH(K,I),K=1,5)
  778     FORMAT(I4,2X,A8,2X,3I10,2X,A8)
C
      IERROR=0
      DO 4 J=1,IFIELD
      NPOINT(J,1)=0
      KADD=KADD1
C
C  BEGIN LOOP OVER FIELDS NAMED ON HEADER
C
      I=0
    2 I=I+1
      NDENS=IH(IPPAK,I)
      KSIZ=2+(NLON+NDENS-1)/NDENS
C
C  DETERMINE NUMBER OF LEVELS OF FIELD NAMED ON HEADER
C
      IF ( MOD(IH(IPLEV,I),10 ) .EQ. 0 ) THEN
         LEVS=1
      ELSE
         LEVS=NLEV
      ENDIF
C
C  COMPARE NAMES OF FIELDS
C
      IF ( MCFLD(1,I) .EQ. CFIELD(J) ) THEN
C
         IF ( NDENS .EQ. 1 )THEN
            NPOINT(J,1) = IH(IPADR,I)
         ELSE
            NPOINT(J,1) = KADD
         ENDIF
C
         NPOINT(J,2) = LEVS
         NPOINT(J,3) = NDENS
      ENDIF
C
      KADD = KADD + LEVS * KSIZ
      IF (NPOINT(J,1).NE.0) GOTO 3      ! POINTER FOUND SUCCESSFULLY
      IF (I.NE.II) GOTO 2               ! END LOOP OVER FIELDS ON HEADER
C
      IF (NPOINT(J,1).EQ.0) THEN
         PRINT 11,J,CFIELD(J)
         IERROR=IERROR+1
      ENDIF
    3 CONTINUE
    4 CONTINUE            ! END LOOP OVER REQUESTED FIELD NAMES
C
      IF (IERROR.NE.0) CALL ABORT(' CCM POINTERS NOT FOUND. SEE OUTPUT')
C
C
   10 FORMAT('0 X X X X   HEADER TYPE = ',I2,' NOT AN OPTION IN RPOINT')
   11 FORMAT('0 X X X X  FIELD ',I2,'   NAMED ',A8,'  NOT FOUND  X X X')
      RETURN
      END
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE RPOINT1 (MTYPE,IH,MCFLD,IDF,II,IFIELD,NPOINT,CFIELD,
     A                   NLON,NLEV)
      INTEGER IH(IDF,II), MTYPE, NPOINT(IFIELD,3)
      CHARACTER*8 MCFLD(2,II),  CFIELD(IFIELD)
C
C  THIS ROUTINE SEARCHES THE ARRAY MCFLD FOR THE NAMES STORED IN THE
C  ARRAY CFIELD.  THE ARRAY IH THEN CONTAINS THE CORRESPONDING POINTER
C  INFORMATION.  MCFLD AND IH ARE LISTS OBTAINED FROM THE CCM HEADER.
C  CFIELD IS THE LIST OF NAMES OF FIELDS DESIRED BY THE USER.  ON       T
C  OUTPUT, THE FIRST VALUE OF NPOINT CONTAINS THE POINTER FOR THAT FIELD
C  WITHIN A CCM DATA BUFFER.  THE SECOND VALUE CONTAINS THE NUMBER OF
C  BUFFER.  THE SECOND VALUE WILL CONTAIN THE NUMBER OF LEVELS (1 OR
C  LEVELS (1 OR NLEV) ON WHICH THAT FIELD IS DEFINED. THE THIRD INDEX
C  INDICATES THE DATA PACKING DENSITY.
C
C  KADD IS A POINTER TO THE FIRST FIELD VALUE FOR PACKED DATA OR THE
C  ADDRESS FOR UNPACKED DATA. KSIZ IS THE NUMBER OF PACKED WORDS
C  CORRESPONDING TO NLONW UNPACKED WORDS. 'UNPKA' IS A SYSTEM ROUTINE.
C
C  SET INDICES FOR FORMATS OF PARTICULAR CCM HEADER
      IF (MTYPE .EQ. 1) THEN        ! CCM0B FORMAT
         IPLEV = 2
         IPADR = 3
         IPPAK = 4
         KADD1 = 7      ! LOCATION OF FIRST WORD ON CCM DATA RECORD
         DO 1 I=1,II    ! CHANGE HOLLARITH TO CHARACTER
         WRITE (MCFLD(1,I),'(A8)') IH(1,I)
    1    CONTINUE
      ELSEIF (MTYPE .EQ. 2) THEN    ! CCM1 FORMAT
         IPLEV = 1
         IPADR = 2
         IPPAK = 3
         KADD1 = 1  ! LOCATION OF FIRST WORD ON CCM DATA RECORD
      ELSEIF (MTYPE .GT. 2) THEN    ! UNKNOWN HEADER TYPE
         PRINT 10,MTYPE
         CALL ABORT ('UNKNOWN HEADER TYPE IN SUB. RPOINT')
      ENDIF
C
C  BEGIN LOOP OVER REQUESTED FIELD NAMES
      DO 777 I=1,II
  777  PRINT 778,I,MCFLD(1,I),MCFLD(2,I),IH(1,I),IH(2,I),IH(3,I)
  778     FORMAT(I4,2X,A8,2X,A8,3I10)
      STOP
C
      IERROR=0
      DO 4 J=1,IFIELD
      NPOINT(J,1)=0
      KADD=KADD1
C
C  BEGIN LOOP OVER FIELDS NAMED ON HEADER
C
      I=0
    2 I=I+1
      NDENS=IH(IPPAK,I)
      KSIZ=2+(NLON+NDENS-1)/NDENS
C
C  DETERMINE NUMBER OF LEVELS OF FIELD NAMED ON HEADER
C
      IF ( MOD(IH(IPLEV,I),10 ) .EQ. 0 ) THEN
         LEVS=1
      ELSE
         LEVS=NLEV
      ENDIF
C
C  COMPARE NAMES OF FIELDS
C
      IF ( MCFLD(1,I) .EQ. CFIELD(J) ) THEN
C
         IF ( NDENS .EQ. 1 )THEN
            NPOINT(J,1) = IH(IPADR,I)
         ELSE
            NPOINT(J,1) = KADD
         ENDIF
C
         NPOINT(J,2) = LEVS
         NPOINT(J,3) = NDENS
      ENDIF
C
      KADD = KADD + LEVS * KSIZ
      IF (NPOINT(J,1).NE.0) GOTO 3      ! POINTER FOUND SUCCESSFULLY
      IF (I.NE.II) GOTO 2               ! END LOOP OVER FIELDS ON HEADER
C
      IF (NPOINT(J,1).EQ.0) THEN
         PRINT 11,J,CFIELD(J)
         IERROR=IERROR+1
      ENDIF
    3 CONTINUE
    4 CONTINUE            ! END LOOP OVER REQUESTED FIELD NAMES
C
      IF (IERROR.NE.0) CALL ABORT(' CCM POINTERS NOT FOUND. SEE OUTPUT')
C
C
   10 FORMAT('0 X X X X   HEADER TYPE = ',I2,' NOT AN OPTION IN RPOINT')
   11 FORMAT('0 X X X X  FIELD ',I2,'   NAMED ',A8,'  NOT FOUND  X X X')
      RETURN
      END
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE SC2GRD4 (D,BM,COEFS,PMN,HMN,MTRUNC,NTRUNC
     A,                   GLAT,G4LON,NLON4,PS,NSPECS,NVECTS
     B,                   MMAXP1,NMAXP1,NMAXPX,IGTYPE)
      PARAMETER (PBAR=1.E2)
      DIMENSION PS(NLON4)
C
C THIS SUBROUTINE CALLS ALL THE SUBROUTINES THAT ARE REQUIRED TO CARRY
C OUT THE NECCESSARY OPERATIONS TO PROJECT THE SPECTRAL FUNCTIONS ON
C TO ONE LATITUDE OF DATA.
C
      CALL ZEROA(D,NVECTS*NLON4)
      CALL ZEROA(BM,MMAXP1*NVECTS*2)
      IF(IGTYPE .EQ. 2)THEN
                      NVECTZ=NVECTS/2
                      CALL SPHER4(BM,COEFS,PMN,HMN,MTRUNC,GLAT
     A,                           NMAXP1,NMAXPX,MMAXP1,NVECTZ,MMAXP1
     B,                           NSPECS,2)
                      CALL NZONALY(D,BM,G4LON,NLON4,MMAXP1,NVECTS)
C                     CALL UVROT4(D,G4LON,NLON4,NVECTZ)
      ELSE
                      CALL SPHER3(BM,COEFS,PMN,MTRUNC,NMAXP1,MMAXP1
     A,                           NMAXPX,NVECTS,MMAXP1,NSPECS,1)
                      CALL NZONALY(D,BM,G4LON,NLON4,MMAXP1,NVECTS)
                      DO 10 I=1,NLON4
   10                 PS(I)=EXP(PS(I))*PBAR
      ENDIF
      RETURN
      END
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE SETPNT (NPOINT,CFIELD,IFIELD,NLEV1,NLEV2)
C
C  THIS ROUTINE CALLS SETPT1 TO FILL A COMMON BLOCK /CPOINT/ OF
C  POINTERS FOR INDICATED FIELDS.  THESE POINTERS ARE FOR THE THIRD
C  INDEX OF ARRAYS SUCH AS B2(NLON,NLAT,NVECT).
C
      CHARACTER*8 CFIELD(IFIELD)
      INTEGER NPOINT(IFIELD,2)
      COMMON /CPOINT/ NPHIS1,NPS1,NTSP1,NTP1,NQP1,NZA1,NUV1
     A               ,NPHIS2,NPS2,NTSP2,NTP2,NQP2,NZA2,NUV2
C
      CALL SETPT1 (NPOINT,CFIELD,IFIELD,NLEV1,'PHIS    ',NPHIS1)
      CALL SETPT1 (NPOINT,CFIELD,IFIELD,NLEV1,'PS      ',NPS1)
      CALL SETPT1 (NPOINT,CFIELD,IFIELD,NLEV1,'T       ',NTP1)
      CALL SETPT1 (NPOINT,CFIELD,IFIELD,NLEV1,'Q       ',NQP1)
      CALL SETPT1 (NPOINT,CFIELD,IFIELD,NLEV1,'ORO     ',NTSP1)
      CALL SETPT1 (NPOINT,CFIELD,IFIELD,NLEV1,'ZA      ',NZA1)
      CALL SETPT1 (NPOINT,CFIELD,IFIELD,NLEV1,'U       ',NUV1)
 
      CALL SETPT1 (NPOINT,CFIELD,IFIELD,NLEV2,'PHIS    ',NPHIS2)
      CALL SETPT1 (NPOINT,CFIELD,IFIELD,NLEV2,'PS      ',NPS2)
      CALL SETPT1 (NPOINT,CFIELD,IFIELD,NLEV2,'T       ',NTP2)
      CALL SETPT1 (NPOINT,CFIELD,IFIELD,NLEV2,'Q       ',NQP2)
      CALL SETPT1 (NPOINT,CFIELD,IFIELD,NLEV2,'ORO     ',NTSP2)
      CALL SETPT1 (NPOINT,CFIELD,IFIELD,NLEV2,'ZA      ',NZA2)
      CALL SETPT1 (NPOINT,CFIELD,IFIELD,NLEV2,'U       ',NUV2)
C
      PRINT*,'CCM PHIS,PS,T,Q,ORO,ZA,U = '
     &      ,NPHIS1,NPS1,NTP1,NQP1,NTSP1,NZA1,NUV1
      PRINT*,'MM4 PHIS,PS,T,Q,ORO,ZA,U = '
     &      ,NPHIS2,NPS2,NTP2,NQP2,NTSP2,NZA2,NUV2
      RETURN
      END
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE SETPT1 (NPOINT,CFIELD,IFIELD,NLEV,CNAME,NPNT1)
      CHARACTER*8 CFIELD(IFIELD), CNAME
      INTEGER NPOINT(IFIELD,2)
C
C  THIS ROUTINE SETS POINTER NPNT1 FOR FIELD CNAME TO VALUE DETERMINED
C  FROM LIST OF FIELD NAMES AND WHETHER THE FIELD IS SINGLE OR MULTI
C  LEVEL (IN THE LATTER CASE, THE NUMBER OF LEVELS IS NLEV).
C
      NPNT1=1
      I=0
    1 I=I+1
      IF (CFIELD(I).EQ.CNAME) RETURN
      IF (NPOINT(I,2).EQ.1) THEN
         NPNT1=NPNT1+1
      ELSE
         NPNT1=NPNT1+NLEV
      ENDIF
      IF (I.LT.IFIELD) GO TO 1
      PRINT 3,CNAME
      CALL ABORT ('FIELD NAME NOT FOUND IN SETPTV')
    2 RETURN
    3 FORMAT('0 FIELD ',A8,' NOT FOUND IN LIST')
      END
C
C X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE SETSIG (SIGMAF,SIGMAH,NK)
C
C THIS ROUTINE DEFINES THE HALF LEVELS WHERE THE VARIABLES
C ARE DEFINED IN THE MM4 MODEL
C
      DIMENSION SIGMAF(NK+1),SIGMAH(NK)
      DO 1 K=1,NK
    1 SIGMAH(K)=0.5*(SIGMAF(K)+SIGMAF(K+1))
      RETURN
      END
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE SFILE1 (B,BCOPY,BR,NPOINT,OPOINT,NBOLD,NBNEW
     A,                  IFIELD,NLONP2,NLONW,KLAT,NLATS,IUNIT)
      INTEGER NPOINT(IFIELD,2), OPOINT(IFIELD,3)
      DIMENSION B(NBOLD), BR(NBNEW), BCOPY(NBOLD)
C
C  THIS ROUTINE SAVES SOME FIELDS ON A FILE
C
      DO 1 N=1,NBOLD
      BCOPY(N)=B(N)
    1 CONTINUE
C
      KXLAT=(KLAT+1)/2                            ! NORTHERN HEMISPHERE
      IF (KXLAT*2-1.NE.KLAT) KXLAT=NLATS-KXLAT+1  ! SOUTHERN HEMISPHERE
      CALL REFORM (BCOPY,BR,OPOINT,NPOINT,NBOLD,NBNEW,NLONP2
     A            ,NLONW,IFIELD)
      WRITE (IUNIT) KXLAT,BR
C
      RETURN
      END
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE SFILE2 (B,NLONP2,IFIELD,NLATS,IUNIT)
      DIMENSION B(NLONP2,NLATS,IFIELD)
C
C  THIS ROUTINE READS A FILE OF SAVED DATA
C
      DO 1 J=1,NLATS
      READ (IUNIT) KLAT,((B(N,KLAT,I),N=1,NLONP2),I=1,IFIELD)
    1 CONTINUE
C
      RETURN
      END
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE SPHER1 (B,COEFS,PMN,MTRUNC,NMAXP1,JSIGN,GAUSWT
     A,           MMAXP1,NMAXPX,NVECTS,NB1,NSPECS,IUV)
C
      PARAMETER (EARTHR=6.37E6)
      COMPLEX B(NB1,NVECTS),COEFS(NSPECS,NVECTS,IUV),UV(NB1,NVECTS,2),CM
      DIMENSION MTRUNC(NMAXP1), PMN(NMAXPX,MMAXP1), HMN(NMAXPX,MMAXP1)
C
C
C  M-1       IS THE ZONAL WAVE NUMBER
C  NM+M-2    IS THE DEGREE OF THE LEGENDRE POLYNOMIAL
C  NM        IS THE INDEX OF THE DIAGONAL ON AN ORDER,DEGREE DIAGRAM
C
      XSIGN=FLOAT(JSIGN)
      J=0
C
      DO 1 NM=1,NMAXP1
      IF(JSIGN.NE.1) XSIGN=-XSIGN
      GE=GAUSWT*XSIGN
      DO 1 M=1,MTRUNC(NM)
      GEPMN=GE*PMN(NM,M)
      J=J+1
      DO 1 L=1,NVECTS
    1 COEFS(J,L,1)=COEFS(J,L,1)+GEPMN*B(M,L)
C
      RETURN
C
C -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
C
      ENTRY SPHER2 (UV,COEFS,PMN,HMN,MTRUNC,JSIGN,GLAT,GAUSWT
     A,    NMAXP1,NMAXPX,MMAXP1,NVECTS,NB1,NSPECS,IUV)
C
C  USE WILLIAMSON'S EQS. 4.2, 3.121, 3.124.
C
      XSIGN=FLOAT(JSIGN)
      J=0
      CSR=GAUSWT/(EARTHR*SQRT(1.-GLAT*GLAT))
      DO 2 NM=1,NMAXP1
      IF (JSIGN.NE.1) XSIGN=-XSIGN
      DO 2 M=1,MTRUNC(NM)
      CM=CMPLX(0.,FLOAT(M-1))
      J=J+1
      FPMN=PMN(NM,M)*XSIGN
      FHMN=HMN(NM,M)*XSIGN*FLOAT(JSIGN)
      DO 2 L=1,NVECTS
      COEFS(J,L,1)=COEFS(J,L,1)
     A     +CSR*(CM*FPMN*UV(M,L,2)+FHMN*UV(M,L,1))
      COEFS(J,L,2)=COEFS(J,L,2)
     B     +CSR*(CM*FPMN*UV(M,L,1)-FHMN*UV(M,L,2))
    2 CONTINUE
C
      RETURN
C
C -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
C
      ENTRY SPHER3 (B,COEFS,PMN,MTRUNC,NMAXP1
     A,           MMAXP1,NMAXPX,NVECTS,NB1,NSPECS,IUV)
C
      J=0
      DO 3 NM=1,NMAXP1
      DO 3 M=1,MTRUNC(NM)
      GE=PMN(NM,M)
      J=J+1
      DO 3 L=1,NVECTS
    3 B(M,L)=B(M,L)+GE*COEFS(J,L,1)
C
      RETURN
C
C -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
C
      ENTRY SPHER4 (UV,COEFS,PMN,HMN,MTRUNC,GLAT
     A,           NMAXP1,NMAXPX,MMAXP1,NVECTS,NB1,NSPECS,IUV)
C
C  USE WILLIAMSON'S EQS. 3.180, 3.181.  HMN IS NEGATIVE OF HMN IN SPHER2
C
      CSR=EARTHR/SQRT(1.-GLAT*GLAT)
      J=0
      DO 4 NM=1,NMAXP1
      DO 4 M=1,MTRUNC(NM)
      J=J+1
      CM=CMPLX(0.,FLOAT(M-1))
      IMN=(NM+M-2)*(NM+M-1)
      IF(IMN.EQ.0) THEN
         XFAC=0.
         ELSE
         XFAC=CSR/FLOAT(IMN)
      ENDIF
      DO 4 L=1,NVECTS
      UV(M,L,1)=UV(M,L,1)
     A    -XFAC*(CM*PMN(NM,M)*COEFS(J,L,2)-HMN(NM,M)*COEFS(J,L,1))
      UV(M,L,2)=UV(M,L,2)
     B    -XFAC*(CM*PMN(NM,M)*COEFS(J,L,1)+HMN(NM,M)*COEFS(J,L,2))
    4 CONTINUE
C
      RETURN
      END
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE UVROT4 (UV,XLONS,CLON,GRIDFC,NLON2,NLAT2,LL)
      DIMENSION UV(NLON2,NLAT2,LL,2),XLONS(NLON2,NLAT2)
C
C CHANGE U AND V FROM TRUE (N,E) TO MAP VALUES (X,Y)
C
      PIR180=ATAN(1.)/45.
      XLONC=CLON*PIR180
      DO 10 J=1,NLAT2
      DO 10 I=1,NLON2
         X=(XLONC-(XLONS(I,J)*PIR180))*GRIDFC
         XS=SIN(X)
         XC=COS(X)
         DO 11 L=1,LL
            D=UV(I,J,L,2)*XS+UV(I,J,L,1)*XC
            UV(I,J,L,2)=UV(I,J,L,2)*XC-UV(I,J,L,1)*XS
            UV(I,J,L,1)=D
   11    CONTINUE
   10 CONTINUE
      RETURN
      END
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE ZEROA (X,NN)
      DIMENSION X(NN)
C
C  THIS ROUTINE ZEROES A VECTOR
C
      DO 1 N=1,NN
    1 X(N)=0.
C
      RETURN
      END
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE INTGTB( PA, ZA, TLAYER, ZMM4
     A                 , TP, ZP, SCCM, NI, NJ, KCCM, TPR, ZPR )
c kdr     NI and NJ have reversed meaning here.
 
      PARAMETER (RGAS=287.04, RLAPSE=-6.5E-03, GRAV=9.8 )
      PARAMETER (RGAS2=RGAS/2., B1=GRAV/RLAPSE )
 
      REAL PA(NI,NJ), ZA(NI,NJ), TLAYER(NI,NJ), ZMM4(NI,NJ)
      REAL TP(NI,NJ,KCCM), ZP(NI,NJ,KCCM), SCCM(KCCM)
 
C  ARRAYS USED IN REVERSING VERTICAL INDEXES OF CCM DATA.
C  INPUT ARRAY(S) NOT AFFECTED WITH REVERSAL OF INDEXES
      REAL TPR(NI,NJ,KCCM), ZPR(NI,NJ,KCCM)
 
C  INTGTB CALCULATES ALL VARIABLES NEEDED TO COMPUTE P* ON THE MM4
C        TOPOGRAPHY.  THE MEAN TEMPERATURE IN THE LAYER BETWEEN
C        THE TOPOGRAPHY AND THE PRESSURE LEVEL ABOVE IS CALULATED
C        BY LINEARLY INTERPOLATING WITH HEIGHT THE TEMPS ON
C        PRESSURE LEVELS.
C        INPUT:    TP        TEMPS ON ECMWF PRESSURE LEVELS
C                  ZP        HEIGHTS OF ECMWF PRESSURE LEVELS
C                  ZMM4      MM4 TOPOGRAPHY
C                  SCCM      ECMWF PRESSURE LEVELS (DIVIDED BY 1000.)
C        OUTPUT:   TLAYER    MEAN LAYER TEMP ABOVE MM4 SURFACE
C                  PA        PRESSURE AT TOP OF LAYER
C                  ZA        HEIGHT AT PRESSURE PA
C
      DO 1 I=1,NI
      DO 1 J=1,NJ
      DO 1 K=1,KCCM
         KR=KCCM-K+1
         TPR(I,J,KR)=TP(I,J,K)
         ZPR(I,J,KR)=ZP(I,J,K)
    1 CONTINUE
C the data seem to be coming from CCR1TIME with the top first
c this puts the bottom first in XXR arrays
 
      DO 40 I=1,NI-1
      DO 40 J=1,NJ-1
 
         KT = 0
c work our way up  THIS CAN NEVER BE SATISFIED!
         DO 30 K=1,KCCM-1
            IF (ZMM4(I,J).LE.ZPR(I,J,K) .AND.
     &          ZMM4(I,J).GT.ZPR(I,J,K+1)    )KT=K
   30    CONTINUE
         KB = KT + 1
 
         IFLAG=0
         IF(KT.NE.0) THEN
            TLAYER(I,J) = ( TPR(I,J,KT) * (ZMM4(I,J)  -ZPR(I,J,KB))
     A                    + TPR(I,J,KB) * (ZPR(I,J,KT)-ZMM4(I,J))  )
     B                    / (ZPR(I,J,KT)-ZPR(I,J,KB))
            TLAYER(I,J) = ( TPR(I,J,KT)+TLAYER(I,J) ) / 2.
            ZA(I,J) = ZPR(I,J,KT)
c kdr ECMWF pressure level in cb
            PA(I,J) = 100. * SCCM(KT)
            IFLAG=IFLAG+1
         ELSE
            TLAYER(I,J) = TPR(I,J,KCCM)
            ZA(I,J) = ZPR(I,J,KCCM)
            PA(I,J) = 100.
         ENDIF
 
   40 CONTINUE
c kdr This section seems to be never used
       PRINT*,'SIGMAR was used in INGTB ',IFLAG,' times'
 
      PRINT *, 'ZMM4, ZPR(KCCM)   =', ZMM4(20,20), ZPR(20,20,KCCM)
      PRINT *, '      TPR(20)   =',            TPR(20,20,KCCM)
      PRINT *, 'TLAYER, ZA, PA =', TLAYER(20,20), ZA(20,20), PA(20,20)
 
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE INTPSN ( PSMM4, ZMM4, PA, ZA, TLAYER,PT,NI,NJ)
c kdr NI and NJ are reversed in meaning here
 
      PARAMETER (RGAS=287.04, RLAPSE=-6.5E-03, GRAV=9.8)
      REAL PSMM4(NI,NJ), ZMM4(NI,NJ)
      REAL PA(NI,NJ), ZA(NI,NJ), TLAYER(NI,NJ)
 
C  EXTRAPOLATE SURFACE PRESSURE FROM CLOSEST PRESSURE LEVEL ABOVE.
C        USE TLAYER CALCULATED IN INTGTB.
C        PSMM4 = SURFACE PRESSURE - PTOP
C
      RL2=RLAPSE/2.
      GDRM=-GRAV/RGAS
      DO 1 I=1,NI-1
      DO 1 J=1,NJ-1
CCC   TB=T(I,J)+RL2*(ZMM4(I,J)-ZCCM(I,J))
      TB = TLAYER(I,J)
      PSMM4(I,J) = PA(I,J) * EXP(GDRM*(ZMM4(I,J)-ZA(I,J))/TB)-PT
c kdr
c      IF (PSMM4(I,J).LE.0.) THEN
      IF (I.EQ.20 .AND. J.EQ.20) THEN
         PRINT *, 'In INTPSN I,J,ZMM4, ZA, PA, PT ='
     &            ,I,J, ZMM4(I,J), ZA(I,J), PA(I,J), PT
         PRINT *, 'TLAYER, PSMM4 = ', TLAYER(I,J), PSMM4(I,J)
      END IF
c end kdr
    1 CONTINUE
 
 
      RETURN
      END
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE Sfavg (favg,f,ni,nj,nk,n1,ng)
      real favg(nk), f(ni,nj,nk)
c
c  compute horizontal averages of fields
c
      ni1=ni-ng
      nj1=nj-ng
      num=(ni1-n1+1)*(nj1-n1+1)
      do 2 k=1,nk
c      IF (k.EQ.10) THEN
c         PRINT*,'t(:,:,10) in Sfavg = '
c         WRITE(*,9) ((f(i,j,10),i=1,ni),j=1,nj)
c      END IF
c 9    FORMAT (//7(/10E10.3,X))
         fg=0.
         do 1 j=n1,nj1
         do 1 i=n1,ni1
            IF (j.EQ.10 .AND. i.eq.10) PRINT*,'f(',i,j,k,') = ',f(i,j,k)
            fg=fg+f(i,j,k)
    1    continue
         favg(k)=fg/num
    2 continue
c
      return
      end
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE CCR1TIME (FIELDS,ASIGMA,BSIGMA,glat,NPOINT,
     $     NDATE,NTIME,NSTEP,CFIELD,NLONPX,NFNL,IFIELD,NSIGS,
     $     NLON,NLAT,NLEV,MMAX,NMAX,KMAX,IUNIT)
      PARAMETER (NDIM=3000, NDIMB=50000, NDIMX=NDIM+1, NDIMY=2*NDIM+1)
      DIMENSION BCCM(NDIMB), IH(NDIM), XH(NDIM)
      CHARACTER*8 CH(NDIM)
C
C  READ 1 TIME OF DATA FROM CCM HISTORY TAPE (INCLUDING HEADER)
C  INSERT SELECTED FIELDS INTO ARRAY FIELDS.
C
c input arguments
      character*(*) cfield(ifield)  ! field name to be read
      integer nlonpx,  ! number of longitudes with wrap-around points
     $        nfnl,    ! number of slices of whole fields  (NVECTS)
     $        ifield,  ! number of fields
     $        nsigs,   ! number of hybrid coefficients
     $        nlon,    ! number of longitudes
     $        nlat,    ! number of latitudes
     $        nlev,    ! number of vertical levels
     $        mmax,    ! m truncation parameter
     $        nmax,    ! n truncation parameter
     $        kmax,    ! k truncation parameter
     $        iunit    ! unit number of a history tape
c output arguments
      real fields(nlonpx,nlat,nfnl), ! fields values
     $     asigma(nsigs),            ! hybrid coefficient a
     $     bsigma(nsigs),            ! hybrid coefficient b
     $     glat(nlat)                ! Gaussia latitudes
      integer npoint(ifield,3)       ! field information
      integer ndate,                 ! date (yymmdd)
     $        ntime,                 ! seconds
     $        nstep                  ! iteration number
c
      CALL CCRHEAD (IH,XH,CH,ASIGMA,BSIGMA,glat,NDATE,NTIME,NSTEP
     A,             CFIELD,NPOINT,NLONW,NBSIZE,NDIM,NDIMB,IFIELD
     B,             NSIGS,NLON,NLAT,NLEV,MMAX,NMAX,KMAX,IUNIT)
      PRINT*,'NLONW from CCRHEAD = ',NLONW
C
      NDIM1=NLONPX*NLAT
      DO 2 L1=1,NLAT
         READ (IUNIT) (BCCM(I),I=1,NBSIZE)
         NPF=1
         L = INT( BCCM(1) )  ! LATITUDE INDEX RUNS FROM SOUTH TO NORTH POLE
         DO 1 I=1,IFIELD
            NPCCM=NPOINT(I,1)
            NLEVS=NPOINT(I,2)
            NDENS=NPOINT(I,3)
            CALL CCRBUFF (FIELDS(1,L,NPF),BCCM(NPCCM)
     A             ,NDENS,NDIM1,NLONW,NLONPX,NLEVS)
            NPF=NPF+NLEVS
    1    CONTINUE
    2 CONTINUE
      PRINT*,'PSCCM(5,5) = ',FIELDS(5,5,2)
C
      RETURN
      END
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE CCRHEAD (IH,XH,CH,ASIGMA,BSIGMA,glat,NDATE,NTIME,NSTEP
     A,                   CFIELD,NPOINT,NLONW,NBSIZE,NDIM,NDIMB,IFIELD
     B,                   NSIGS,NLON,NLAT,NLEV,MMAX,NMAX,KMAX,IUNIT)
      PARAMETER (LSIZE= 6, LNLON= 9, LLONW=10, LNLAT=11
     A,          LNLEV=12, LMMAX=13, LNMAX=14, LKMAX=15
     B,          LNFLD=16, LSTEP=18, LDATE=26, LTIME=27
     C,          LLHDC=31, LLHDR=32, LSIGS=33, LGAUS=34
     D,          LGAUW=35, LFLDS=36, LCFLD=37, LVNAM=38 )
      DIMENSION IH(NDIM), XH(NDIM), NPOINT(IFIELD,3), ASIGMA(NSIGS)
     A,         BSIGMA(NSIGS), glat(nlat)
      CHARACTER*8 CH(NDIM), CFIELD(IFIELD)
C
C  THIS ROUTINE READS A CCM2 HEADER.  OUTPUT ARE VALUES FOR: NPOINT,
C  ASIGMA, BSIGMA, NLONW, NBSIZE, NDATE, NTIME, AND NSTEP.
C  VALUES FOR NLON, NLEV, NLAT, MMAX, NMAX, KMAX, AND NDIMB ARE CHECKED.
C
C
C  READ HEADER RECORDS
C
      READ (IUNIT) LENHD,(IH(J),J=2,LENHD)
      MTYPE=MOD(IH(2),10)
      IF (MTYPE.NE.3) THEN
         PRINT 100,IH(2)
         STOP 'CCRHEAD0'
      ENDIF
      IF ((IH(LLHDC).GT.NDIM).OR.(IH(LLHDR).GT.NDIM)) THEN
         PRINT 101,IH(LLHDC),IH(LLHDR),NDIM
         STOP 'CCRHEAD1'
      ENDIF
      READ (IUNIT) (CH(J),J=1,IH(LLHDC))
      READ (IUNIT) (XH(J),J=1,IH(LLHDR))
C
C  CHECK THAT TRUNCATION PARAMETERS AGREE.
C
      IERROR=0
      IF (MMAX.NE.IH(LMMAX)) IERROR=IERROR+1
      IF (NMAX.NE.IH(LNMAX)) IERROR=IERROR+1
      IF (KMAX.NE.IH(LKMAX)) IERROR=IERROR+1
      IF (NLEV.NE.IH(LNLEV)) IERROR=IERROR+1
      IF (NLON.NE.IH(LNLON)) IERROR=IERROR+1
      IF (NLAT.NE.IH(LNLAT)) IERROR=IERROR+1
      IF(IERROR.NE.0) THEN
         PRINT 102,MMAX,NMAX,KMAX,NLEV,NLON,NLAT
     A      ,IH(LMMAX),IH(LNMAX),IH(LKMAX),IH(LNLEV),IH(LNLON),IH(LNLAT)
         STOP 'CCRHEAD2'
      ENDIF
C
C  DETERMINE AND PRINT FILE DATE
C
      NDATE=IH(LDATE)
      NTIME=IH(LTIME)
      NSTEP=IH(LSTEP)
      CALL CCRDATE (NDATE,NTIME,NSTEP)
C
C  DETERMINE AND PRINT FILE NAME
C
      PRINT 103,(CH(I),I=LVNAM,LVNAM+9)
C
C  A AND B ARRAYS FOR HYBRID COORDINATE
C
      N1=IH(LSIGS)-1 + NSIGS
      N2 = N1 + NSIGS
      DO 1 N=1,NSIGS
         ASIGMA(N)=XH(N1+N)
         BSIGMA(N)=XH(N2+N)
    1 CONTINUE
c
      n1=ih(lgaus)-1
      do lat=1,nlat
         glat(lat)=xh(n1+lat)
      end do
c
C
C  DETERMINE FIELD POINTERS
C
      NLONW= IH(LLONW)
      MPF=   IH(LFLDS)
      NPF=   IH(LNFLD)
      MPCFLD=IH(LCFLD)
      IDF=3
      CALL CCRPOINT (MTYPE,IH(MPF),CH(MPCFLD),IDF,NPF,IFIELD,NPOINT,
     A               CFIELD,NLON,NLEV)
C
C  CHECK DIMENSION SIZE OF B1 ARRAY
C
      IMAX=0
      DO 3 N=1,IFIELD
         IF (IMAX.LT.NPOINT(N,1)) THEN
            IMAX=NPOINT(N,1)
            M=N
         ENDIF
    3 CONTINUE
      KSIZE=NLONW
      NDENS=NPOINT(M,3)
      IF (NDENS.NE.1) KSIZE=2+(NLON+NDENS-1)/NDENS
      NBSIZE=IMAX-1+KSIZE*NPOINT(M,2)
      IF(NBSIZE .GT. NDIMB) THEN
         PRINT 104,NBSIZE,NDIMB
         STOP 'CCRHEAD4'
      ENDIF
C
      RETURN
C
  100 FORMAT('0 TYPE OF HEADER ON CCM FILE CANNOT BE READ BY CCRHEAD',/
     A,5X,'HEADER TYPE ON FILE IS ',I6)
  101 FORMAT('0 LENGTH OF HEADER RECORDS ON CCM FILE TOO LONG FOR XH,CH'
     A,/,5X,'LENGTHS OF RECORDS ARE:',2I7,5X,'LENGTH OF ARRAY IS:',I7)
  102 FORMAT('0 PARAMETER VALUES ON CCM FILE AND THIS PROGRAM DIFFER',/
     A,5X,'PROGRAM  VALUES OF MMAX,NMAX,KMAX,NLEV,NLON,NLAT ARE ',6I6,/
     B,5X,'CCM FILE VALUES OF MMAX,NMAX,KMAX,NLEV,NLON,NLAT ARE ',6I6)
  103 FORMAT('  CCM FILE NAME ON HEADER = ',10A8)
  104 FORMAT('0 SIZE OF DATA ARRAY TOO SHORT FOR CCM FILE RECORDS',/
     A,5X,'REQUIRED SIZE IS:',I7,';   ARRAY SIZE IS:',I7)
      END
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE CCRDATE (NDATE,NTIME,NSTEP)
      COMMON /CDATEI/ IDATE
C
C  UNPACK DATE AND TIME AND PRINT IN APPROPRIATE FORMAT
C
      IDATE=NDATE*100+NTIME/3600
      IYEAR=NDATE/10000
      IMONTH=(NDATE-IYEAR*10000)/100
      IDAY=NDATE-IYEAR*10000-IMONTH*100
      IHOUR=NTIME/3600
      IHOUR1=IHOUR/10
      IHOUR2=IHOUR-IHOUR1*10
      IMINUT=(NTIME-IHOUR*3600)/60
      IMIN1=IMINUT/10
      IMIN2=IMINUT-IMIN1*10
C
      PRINT 10,IYEAR,IMONTH,IDAY,IHOUR1,IHOUR2,IMIN1,IMIN2,NSTEP
   10 FORMAT('0 DATE OF ANALYZED DATA YR/MO/DA = ',I2,'/',I2,'/',I2
     A,5X,'TIME =',4I1,'Z',5X,'NSTEP =',I8)
C
      RETURN
      END
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE CCRPOINT (MTYPE,IH,MCFLD,IDF,II,IFIELD,NPOINT,CFIELD,
     A                     NLON,NLEV)
      INTEGER IH(IDF,II), MTYPE, NPOINT(IFIELD,3)
      CHARACTER*(*) MCFLD(2,II),  CFIELD(IFIELD)
C
C  THIS ROUTINE SEARCHES THE ARRAY MCFLD FOR THE NAMES STORED IN THE
C  ARRAY CFIELD.  THE ARRAY IH THEN CONTAINS THE CORRESPONDING POINTER
C  INFORMATION.  MCFLD AND IH ARE LISTS OBTAINED FROM THE CCM HEADER.
C  CFIELD IS THE LIST OF NAMES OF FIELDS DESIRED BY THE USER.  ON
C  OUTPUT, THE FIRST VALUE OF NPOINT CONTAINS THE POINTER FOR THAT FIELD
C  WITHIN A CCM DATA BUFFER.  THE SECOND VALUE CONTAINS THE NUMBER OF
C  BUFFER.  THE SECOND VALUE WILL CONTAIN THE NUMBER OF LEVELS (1 OR
C  LEVELS (1 OR NLEV) ON WHICH THAT FIELD IS DEFINED. THE THIRD INDEX
C  INDICATES THE DATA PACKING DENSITY.
C
C  SET INDICES FOR FORMATS OF PARTICULAR CCM HEADER
      IF (MTYPE .EQ. 3) THEN    ! CCM2 FORMAT
         IPLEV = 1
         IPADR = 2
         IPPAK = 3
      ELSE                      ! UNKNOWN HEADER TYPE
         PRINT 10,MTYPE
         STOP 'CCRPOINT'
      ENDIF
C
C  BEGIN LOOP OVER REQUESTED FIELD NAMES
C
      IERROR=0
      DO 4 J=1,IFIELD
      NPOINT(J,1)=0
C
C  BEGIN LOOP OVER FIELDS NAMED ON HEADER
C
      I=0
    2 I=I+1
      NDENS = IH(IPPAK,I)
C
C  COMPARE NAMES OF FIELDS
C
      IF ( MCFLD(1,I) .EQ. CFIELD(J) ) THEN
C
         NPOINT(J,1) = IH(IPADR,I)
C
         IF ( MOD(IH(IPLEV,I),10 ) .EQ. 0 ) THEN
            NPOINT(J,2) = 1
         ELSE
            NPOINT(J,2) = NLEV
         ENDIF
C
         NPOINT(J,3) = NDENS
      ENDIF
C
      IF (NPOINT(J,1).NE.0) GOTO 3      ! POINTER FOUND SUCCESSFULLY
      IF (I.NE.II) GOTO 2               ! END LOOP OVER FIELDS ON HEADER
C
      IF (NPOINT(J,1).EQ.0) THEN
         PRINT 11,J,CFIELD(J)
         IERROR=IERROR+1
      ENDIF
    3 CONTINUE
    4 CONTINUE            ! END LOOP OVER REQUESTED FIELD NAMES
C
      IF (IERROR.NE.0) STOP 'CCRPOINT'
      RETURN
C
   10 FORMAT('0 X X X X   HEADER TYPE = ',I2,' NOT AN OPTION IN ',
     A'ROUTINE CCRPOINT   X X X X X')
   11 FORMAT('0 X X X X  FIELD ',I2,'   NAMED ',A8,'  NOT FOUND IN',
     A'ROUTINE CCRPOINT   X X X X X')
      END
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE CCRBUFF (FIELD,BCCM,NDENS,NDIMF1,NLONW,NLONPX,NLEV)
      DIMENSION FIELD(NDIMF1,1), BCCM(1)
C
C  THIS COPIES DESIRED DATA FROM CCM LATITUDE DATA BUFFER INTO FIELD.
C  IF NLEV.LT.0 THEN COPY ONLY ONE LEVEL (=-NLEV) OF MULTILEVEL FIELD,
C  OTHERWISE NLEV IS NUMBER OF LEVELS OF FIELD.
C  NDENS IS PACKING DENSITY ON CCM FILE (NDENS=1 FOR NO PACKING)
C
      IF (NDENS.NE.1) THEN
         NSKIP=2+(NLONW+NDENS-1)/NDENS
      ELSE
         NSKIP=NLONW
      ENDIF
C
      IF (NLEV.LT.0) THEN
         N1=NSKIP*(-(1+NLEV))
         NLEVS=1
      ELSE
         N1=0
         NLEVS=NLEV
      ENDIF
C
C  COPY (NDENS=1) OR UNPACK AND COPY (NDENS.NE.1) DATA BUFFER TO FIELD
C
c kdr
      IF (NLEVS.EQ.1) THEN
         DO 2 L=1,NLEVS
            IF (NDENS.EQ.1) THEN         ! DATA IN UNPACKED FORMAT
               DO 1 N=1,NLONW
                  FIELD(N,L) = BCCM(N1+N)
    1          CONTINUE
            ELSE                         ! DATA IN PACKED FORMAT
               CALL UNPKA (FIELD(1,L),NLONW,BCCM(N1+1),NDENS)
            ENDIF
            N1 = N1+NSKIP
    2    CONTINUE
      ELSE
         DO 5 L=NLEVS,1,-1
            IF (NDENS.EQ.1) THEN         ! DATA IN UNPACKED FORMAT
               DO 4 N=1,NLONW
                  FIELD(N,L) = BCCM(N1+N)
    4          CONTINUE
            ELSE                         ! DATA IN PACKED FORMAT
               CALL UNPKA (FIELD(1,L),NLONW,BCCM(N1+1),NDENS)
            ENDIF
            N1 = N1+NSKIP
    5    CONTINUE
      END IF
C
C  FILL WRAP AROUND POINTS
C
      IF (NLONPX.GT.NLONW) THEN
         DO 3 N=NLONW+1,NLONPX
         DO 3 L=1,NLEVS
         FIELD(N,L)=FIELD(N-NLONW,L)
    3    CONTINUE
      ENDIF
C
      RETURN
      END
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE CCRGETP (NPOINT,CFIELD,IFIELD,NFDIM2,NLEV,NPNT1,CNAME)
      CHARACTER*(*) CFIELD(IFIELD), CNAME
      INTEGER NPOINT(IFIELD,2)
C
C  THIS ROUTINE SETS POINTER NPNT1 FOR FIELD CNAME TO VALUE DETERMINED
C  FROM LIST OF FIELD NAMES AND WHETHER THE FIELD IS SINGLE OR MULTI
C  LEVEL (IN THE LATTER CASE, THE NUMBER OF LEVELS IS NLEV).
C
      NPNT1=1
      DO 1 I=1,IFIELD
      NLEV=NPOINT(I,2)
      IF (CFIELD(I).EQ.CNAME) RETURN
      NPNT1=NPNT1+NLEV*NFDIM2
    1 CONTINUE
C
      PRINT 2,CNAME
      STOP 'CCRGETNP'
    2 FORMAT('0 FIELD ',A8,' NOT FOUND IN LIST BY CCRGETP')
      END
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      subroutine ccrtest(f,nlonpx,nlon,nlat,nlev,cfield)
      implicit none
      integer nlonpx,nlon,nlat,nlev
      real f(nlonpx,nlat,nlev)
      character cfield*8
c
c  Compute max, min, and mean for complete field not including
c  wraparound points.
c
      integer nsum, i, j, k
      real fmax, fmin, fsum
c
      nsum = nlon*nlat*nlev
      fmax = f(1,1,1)
      fmin = fmax
      fsum = 0.
c
      do k = 1, nlev
         do j = 1, nlat
            do i = 1, nlon
               if( f(i,j,k) .lt. fmin ) fmin = f(i,j,k)
               if( f(i,j,k) .gt. fmax ) fmax = f(i,j,k)
               fsum = fsum + f(i,j,k)
            end do
         end do
      end do
c
      fsum = fsum/float( nsum )
c
      print 10, cfield
      print 12, fmax
      print 13, fmin
      print 14, fsum
c
   10 format('0  ',a8,'  field statistics')
   12 format(' max          ', 1pe22.15 )
   13 format(' min          ', 1pe22.15 )
   14 format(' mean abs val ', 1pe22.15 )
c
      return
      end
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      subroutine ccmodht( q, nlonpx, nlat, nlev, iunitin, iunitout )
c
c Read a ccm2 history tape, and write a new tape containing the modifications
c to the q field.
c Files should be opened in the calling routine.
c
      implicit none
c
      integer nlonpx, nlat, nlev, iunitin, iunitout
      real q(nlonpx,nlat,nlev)
c
      integer ifield
      parameter( ifield=1 )
      character*8 cfield(ifield)
      data cfield / 'Q       ' /
c
c Buffers for reading history tape.
c
      integer ndimi, ndimc, ndimr, ndimd
      parameter( ndimi=500, ndimc=500, ndimr=500, ndimd=60000 )
      integer hdi(ndimi)
      real hdr(ndimr), datrec(ndimd)
      character*8 hdc(ndimc)
c
c Variables contained in the integer header, and pointers to those
c values.
c
      integer lenhdi, mtype, maxsize, nlon, nlonw, nfldh,
     *        lenhdc, lenhdr, mpflds, mpcfld
      integer imftyp, imaxsize, inlon, inlonw,
     *        infldh, ilenhdc, ilenhdr, impflds, impcfld
      data imftyp/2/, imaxsize/6/, inlon/9/, inlonw/10/,
     *     infldh/16/, ilenhdc/31/, ilenhdr/32/, impflds/36/,
     *     impcfld/37/
c
c Pointers into ccm2 data records for fields to be replaced.
c
      integer npoint(ifield,3)
c
c Misc. local variables.
c
      integer i, j, jlat
      real lat
c
c     Read the history tape header to obtain the pointers to the data
c     fields.
c
      rewind( iunitin )
      read( iunitin ) lenhdi, (hdi(i),i=2,lenhdi)
c
      mtype = mod( hdi(imftyp), 10 )
      if( mtype .ne. 3 )then
         write(*,*)'ccmodht: wrong history tape header type'
         stop
      endif
c
      lenhdc = hdi(ilenhdc)
      lenhdr = hdi(ilenhdr)
c
      if( ( lenhdc .gt. ndimc )  .or.  ( lenhdr .gt. ndimr) )then
         write(*,*)'ccmodht: arrays for reading header files too small'
         stop
      endif
      read( iunitin ) (hdc(i),i=1,lenhdc)
      read( iunitin ) (hdr(i),i=1,lenhdr)
c
c     Assign variables out of the integer header.
c
      maxsize = hdi(imaxsize)
      nlonw = hdi(inlonw)
      nfldh = hdi(infldh)
      mpflds = hdi(impflds)
      mpcfld = hdi(impcfld)
c
c     Set npoint.
c
      call ccrpoint( mtype, hdi(mpflds), hdc(mpcfld), 3, nfldh, ifield,
     *               npoint, cfield, nlon, nlev )
c
c     Check that none of the fields is packed.
c
      do i = 1, ifield
         if( npoint(i,3) .ne. 1 ) then
            write(*,*)'ccmodht: packed fields not yet implimented'
            stop
         end if
      end do
c
c     Update character header record with info about dataset being made.
c
cxxxxxxxxxxxxxxx not yet implimented xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c
c     Write header records for new dataset.
c
      write( iunitout ) lenhdi, (hdi(i),i=2,lenhdi)
      write( iunitout ) (hdc(i),i=1,lenhdc)
      write( iunitout ) (hdr(i),i=1,lenhdr)
c
c     Modify the data records.
c
      if( maxsize .gt. ndimd )then
         write(*,*)'ccmodht: array for data records too small'
         stop
      end if
c
      do j = 1, nlat
c
         read( iunitin ) lat, (datrec(i),i=2,maxsize)
         jlat = int( lat )
c
         call ccmhtrec( datrec, maxsize, q, nlonw, nlonpx, nlat,
     *                  nlev, jlat, npoint(1,1) )
c
         write( iunitout ) lat, (datrec(i),i=2,maxsize)
c
      end do    ! loop over latitude records in history tape
c
      return
      end
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      subroutine ccmodht1(fs,cfield,nlonpx,nlat,nlev,iunitin,iunitout)
c
c Read a ccm2 history tape, and write a new tape containing the modifications
c to one of the contained fields.
c Files should be opened in the calling routine.
c
      implicit none
c
      integer nlonpx, nlat, nlev, iunitin, iunitout
      real fs(nlonpx,nlat,nlev)
      character*(*) cfield
c
c Buffers for reading history tape.
c
      integer ndimi, ndimc, ndimr, ndimd
      parameter( ndimi=500, ndimc=500, ndimr=500, ndimd=60000 )
      integer hdi(ndimi)
      real hdr(ndimr), datrec(ndimd)
      character*8 hdc(ndimc)
c
c Variables contained in the integer header, and pointers to those
c values.
c
      integer lenhdi, mtype, maxsize, nlon, nlonw, nfldh,
     *        lenhdc, lenhdr, mpflds, mpcfld
      integer imftyp, imaxsize, inlon, inlonw,
     *        infldh, ilenhdc, ilenhdr, impflds, impcfld
      data imftyp/2/, imaxsize/6/, inlon/9/, inlonw/10/,
     *     infldh/16/, ilenhdc/31/, ilenhdr/32/, impflds/36/,
     *     impcfld/37/
c
c Pointers into ccm2 data records for fields to be replaced.
c
      integer ifield
      parameter (ifield=1)
      integer npoint(ifield,3)
c
c Misc. local variables.
c
      integer i, j, jlat
      real lat
c
c     Read the history tape header to obtain the pointers to the data
c     fields.
c
      rewind( iunitin )
      read( iunitin ) lenhdi, (hdi(i),i=2,lenhdi)
c
      mtype = mod( hdi(imftyp), 10 )
      if( mtype .ne. 3 )then
         write(*,*)'ccmodht1: wrong history tape header type'
         stop
      endif
c
      lenhdc = hdi(ilenhdc)
      lenhdr = hdi(ilenhdr)
c
      if( ( lenhdc .gt. ndimc )  .or.  ( lenhdr .gt. ndimr) )then
         write(*,*)'ccmodht1: arrays for reading header files too small'
         stop
      endif
      read( iunitin ) (hdc(i),i=1,lenhdc)
      read( iunitin ) (hdr(i),i=1,lenhdr)
c
c     Assign variables out of the integer header.
c
      maxsize = hdi(imaxsize)
      nlonw = hdi(inlonw)
      nfldh = hdi(infldh)
      mpflds = hdi(impflds)
      mpcfld = hdi(impcfld)
c
c     Set npoint.
c
      call ccrpoint( mtype, hdi(mpflds), hdc(mpcfld), 3, nfldh, ifield,
     *               npoint, cfield, nlon, nlev )
c
c     Check that none of the fields is packed.
c
      do i = 1, ifield
         if( npoint(i,3) .ne. 1 ) then
            write(*,*)'ccmodht1: packed fields not yet implimented'
            stop
         end if
      end do
c
c     Update character header record with info about dataset being made.
c
cxxxxxxxxxxxxxxx not yet implimented xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c
c     Write header records for new dataset.
c
      write( iunitout ) lenhdi, (hdi(i),i=2,lenhdi)
      write( iunitout ) (hdc(i),i=1,lenhdc)
      write( iunitout ) (hdr(i),i=1,lenhdr)
c
c     Modify the data records.
c
      if( maxsize .gt. ndimd )then
         write(*,*)'ccmodht1: array for data records too small'
         stop
      end if
c
      do j = 1, nlat
c
         read( iunitin ) lat, (datrec(i),i=2,maxsize)
         jlat = int( lat )
c
         call ccmhtrec( datrec, maxsize, fs, nlonw, nlonpx, nlat,
     *                  nlev, jlat, npoint(ifield,1) )
c
         write( iunitout ) lat, (datrec(i),i=2,maxsize)
c
      end do    ! loop over latitude records in history tape
c
      return
      end
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      subroutine ccmhtrec( datrec, maxsize, fld, nlonw, nlonpx, nlat,
     *                      nlev, jlat, npoint )
c
c Write data from multi-dimensional array into history tape data record.
c
      implicit none
c
      integer maxsize, nlonw, nlonpx, nlat, nlev, jlat
      integer npoint
      real datrec(maxsize), fld(nlonpx,nlat,nlev)
c
      integer i, ii, k
c
      ii = npoint - 1
      do k = 1, nlev
         do i = 1, nlonw
            datrec(ii+i) = fld(i,jlat,k)
         end do
         ii = ii + nlonw
      end do
c
      return
      end

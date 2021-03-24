	PROGRAM LNCINI
c 
c search for 'user set' to find user set parameters
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C	THIS IS INTENDED TO BE THE MAIN PROGRAM TO 
C	COMPUTE OPTIMAL MODES FOR THE MM4 MODEL 
C	USING THE NAG - LANCZOS ROUTINES. 
C
C	MARTIN EHRENDORFER, 1 AUGUST 1993.
C
C	THIS DRIVER ORIGINATES FROM LANCBL10.F, FOR THE 
C	BATES/LI MODEL. IT CALLS THE 
C	ROUTINES CLANC AND VRETR (WHICH I WROTE, AND WHICH 
C	BY THEMSELVES CALL THE LANCZOS ROUTINES). 
C
C	ALSO CLANC REQUIRES A ROUTINE NAMED ATX93, WHICH 
C	IS TAKEN FROM TLMADJ1 PROVIDED BY RON ERRICO. 
C
C	SEVERAL CHANGES MADE ON 2 SEPTEMBER 1993 FOR EASIER READING. 
C
C	MAIN CHANGES ON 30 SEPTEMBER 1993, TO SAVE AND 
C	UTILIZED INFORMATION FROM PREVIOUS RUNS. 
C
C	RATHER FINAL VERSION OF INTERFACE ROUTINE ATX93 ON 5 OCT 93. 
C	FURTHER IMPROVEMENTS 25 OCTOBER 1993. 
C
C	PUT THE INITIALIZATION INTO CATX TOGETHER WITH RON, 
C	18 NOVEMBER 1993 (STARTING FROM LANCMM5). 
C
C	IN ADDITION TO THE INITIALIZATION, THIS ALSO HAS THE OPTION 
C	TO DO THE POTENTIAL ENSTROPHY NORM. 12 DECEMBER 1993. 
C
C       For restart, acquire file of old vectors that were outout at
C       run time (name given by APTH below) and change ISAVE to non 0
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C	STARTING FROM LNCINI3.F I AM TRYING TO MAKE ONE UNIFIED 
C	DRIVER FOR ALL DIFFERENT NORMS AS WELL AS TO PRODUCE THE 
C	DIMENSIONED VECTORS. 14 DECEMBER 1993. M.E. 
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C	
C	This is the driver to use the Lanczos routines in connection 
C	with the TLM of MM4. In this version there are five different
C	options for the inner product (see parameter list below), 
C	and also, in this version, there is the option to run 
C	initialized or not. The first version of this was 
C	completed 15 December 1993.
C
C	What follows is a short description of the rule how to use 
C	this driver. Basically, all switches to control the Lanczos
C 	routines and the choice of the inner product need to be set 
C	in this driver. But, switches to control the TLM must be set
C	in CATX (like initialization yes/no, or optimization time). 
C	Never set anything in ATX93. However, Tomis routines
C	need some adjustment in CATX. 
C
C	The choice of the inner product is controlled through the 
C	ISET parameters. It is important that only one of these is
C	set to one for one run, the rest must be zero. LANMAX should 
C	be set to the number of iterations to be performed. 
C
C	This driver will perform a Lanczos iteration, and if that 
C	is successful, also produce a file with the dimensioned 
C	eigenvectors (named vectordim.dat). Also, a test of selected
C	eigenvectors can be performed. Further, this driver has the
C	full restart facility (mainly controlled by ISAVE). 
C
C	An important variable to set before a run if APTH, which is
C	the path in the mass store, where results from the computations
C	are written to during the run. 
C
C	This is the rewritten Lanczos driver for the new version of
C	tlmadj (with xnorms-routines). 1 July 1994, M. Ehrendorfer. 
c
c	adapted for moist R-norm on 12/1/1994, M. Ehrendorfer. 
c
c	option 7 included, 12/20/94, m.e.
c
c	starting from lnfull.f.37 -> new hack scheme, 6/13/95, m.e.
c       1/22/96: adapt this driver to the new library (with RAS
c       scheme and vertical differencing scheme). we use energy
c       norm with q ignored (iset1=1; can also work dry), and
c       energy norm with qegt (iset8=1). m.e. remove iset5=1.
C
c	1/17/97 ... adapt this driver for MAMS-2. m.e. with the moist norms
c	1/30/97 ... driver adapted for new R-norm (rme + me)
c       2/3/97  ... TE-norm adapted for # grid points 
C       8/27/98 ... result.dat.X and rm Ooutput96 changes
C                   to avoid problems with restart due to overwriting
C                   record 1 of first file with record LANMAX+1.
C
C
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
c user set
c grid dimensions
        PARAMETER ( NI = 43 , NJ = 54 , NK = 10 ) ! set in TLM/ADJ TOO
C norm choice
	PARAMETER ( ISET1=0,ISET2=0,ISET3=0,ISET4=1,ISET5=0,ISET6=0,
     1              ISET7=0,ISET8=0,ISET9=0 )
C norm truncation parameters
  	PARAMETER ( ITRC = NI/2 , JTRC = NJ/2 , KTRC = 7 )    
c nontruncated (excludes edge values)
C  	PARAMETER ( ITRC = NI-3 , JTRC = NJ-3 , KTRC = NK )    
c moisture truncation
        PARAMETER ( ITRCQ = NI/2 , JTRCQ = NJ/2 , KTRCQ = 4)
c end user set
	PARAMETER ( NI2 = NI-2, NJ2 = NJ-2 )
	PARAMETER ( NI3 = NI-3, NJ3 = NJ-3 ) 
	PARAMETER ( NFULL = 2*NK*NI2*NJ2 + (NK+1)*NI3*NJ3 +
     1                      NK*NI3*NJ3 )
	PARAMETER ( ITYPE = 1*ISET1+2*ISET2+3*ISET3+4*ISET4+5*ISET5 
     1                     +6*ISET6+7*ISET7+8*ISET8+9*ISET9 )
C
c******************************************************************
C	ISET1 = 1 TO USE THE TOTAL ENERGY NORM ROUTINES
c          2/8/97: TE-norm with q-ignored (with # grid points)
c
C	ISET2 = 1 TO USE OPTION 1B OF ENERGY NORM (PROJECT AT T1)
c          2/8/97: ** i n a c t i v e ** at this point
c
C	ISET3 = 1 TO USE OPTION 2A OF ENERGY NORM (PROJECT AT T0, T1)
c          2/8/97: ** i n a c t i v e ** at this point
c
C	ISET4 = 1 TO USE ROTATIONAL MODES NORM (RON) NEW FEB 97
c          2/8/97: R-norm compatible with TE-norm
c
C	ISET5 = 1 TO USE THE SET OF TOMIS NORM ROUTINES
c          2/8/97: ** i n a c t i v e ** removed from this driver
c
C	ISET6 = 1 TO USE MOIST R-NORM (M.E., 12/1/94) 
c          2/8/97: uses new R-norm, qwgt adaption
c
C	ISET7 = 1 TO USE Q ONLY IN THE NORM (M.E., 12/20/94)
c          2/8/97: moisture growth
c
c       ISET8 = 1 TO USE TOTAL ENERGY PLUS (QWGT * Q)^2
C          2/8/97: qwgt^2 = L^2 / (c_p T_r)
c
c       ISET9 = R NORM PLUS SPECTRAL Q NORM
C          3/5/97: qwgt^2 = L^2 / (c_p T_r)
c******************************************************************
C
	PARAMETER ( N = ( NFULL - NK*NI3*NJ3 ) * ISET1 + 
     2                  ( NFULL )              * ISET2 + 
     3                  ( 2*NK*NI2*NJ2 )       * ISET3 + 
     4                  ( ITRC * JTRC * KTRC ) * ISET4 + 
     5                  ( NFULL )              * ISET5 +
     6         ( ITRC*JTRC*KTRC + NK*NI3*NJ3 ) * ISET6 +
     7                  ( NK*NI3*NJ3 )         * ISET7 + 
     8                  ( NFULL )              * ISET8 +
     9  ( ITRC*JTRC*KTRC + ITRCQ*JTRCQ*KTRCQ ) * ISET9    )
C
c user set
c LANMAX should be "three times the number of vectors you want to
c have accurately converged (e.g., set to 30 if you need 10 vectors)."
c This ratio can be smaller for large numbers of iterations, 
c and should be larger for smaller numbers of iterations (<20)
c     Also; the first batch of iters will consist of LANMAX+1 iters, 
c due to the (extra) random initial condition iteration at the beginning.
c To fill tapes evenly set LANMAX = (#TAPES*INUM) -1
c So with INUM=100, tape 1; Random...1...99=100  LANMAX= 99
c                   tape 2  100...199      =100  LANMAX=199
c
 	PARAMETER ( LANMAX = 199 )  ! number of iterations 
C
  	PARAMETER ( JEVEC =  20 )  ! max number of vects to save  
c end user set
C
	PARAMETER ( NW = 6*N+1+4*LANMAX+LANMAX*LANMAX , NP2=N+2 )
C
	PARAMETER ( NFP2 = NFULL + 2 )
C
C	SET THE UNITS FOR THE FILES TO BE USED
C	LUNIT1 FOR PROTOCOL FILE IN THIS DRIVER 
C	LUNIT2 FOR DIRECT ACCESS FILE TO CONTAIN EIGENVECTORS
C	LUNIT3 FOR INTERNAL SEQUENTIAL ACCESS FILE IN LANCZOS ROUTINES
C	NEVER USE UNIT = 94, SINCE THIS IS RESERVED FOR L-ROUTINES
C	NEVER USE UNIT = 95, IT IS USED FOR STORING/READING INFORMATION
C
	PARAMETER ( LUNIT2 = 92, LUNIT3 = 93 ) 
C
C	KAPPA DETERMINES WHETHER E-VECTOR IS WRITTEN OR NOT
C
	REAL KAPPA
	PARAMETER ( KAPPA = 1.0E-02 )
C
	REAL EV (LANMAX), BNDEV (LANMAX), W (NW), EVEC (N,JEVEC)
	INTEGER INDEV (LANMAX) , INDVEC (LANMAX) 
	REAL X1 (N) , X2 (N) , X3 (N) , X4 (N) , XFULL ( NFULL )
C
	CHARACTER*128 TNAME
	CHARACTER*12  FNAME
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C	SET MDIM  NUMBER OF VECTORS TO BE DIMENSIONED
C	SET MTEST NUMBER OF VECTORS TO BE TESTED 
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
 	PARAMETER ( MDIM =  60 , MTEST = 0 ) 
C	PARAMETER ( MDIM = 200 , MTEST = 0 ) 
C
C	LUNIT6 FOR DIRECT ACCESS FILE TO CONTAIN DIMENSIONED VECTORS
C
	PARAMETER ( LUNIT6 = 96 ) 
C
	CHARACTER*12 APTH
	COMMON /IFIRST/ 
     1  IFIRST,LUNIT1,ISCALE,NCALLS,INUM,ISAVE,IOUT,ITYPET,APTH
C
	COMMON /ISECND/ XFULL 
C
	LUNIT1 = 91
	OPEN (UNIT=LUNIT1,FILE='Ooutput91')
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C	PREPARING CALL TO ATX93 TO ACQUIRE AND OPEN FILES 
C
	IFIRST = 0
	CALL ATX93 ( X1 , N )
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C	PREPARING CALL TO xnorms0 to set norm-parameters 
C
        CALL xnorms0 (ni,nj,nk,itrc,jtrc,ktrc,itrcq,jtrcq,ktrcq,itype)
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C
C	PREPARE TO CALL LANCZOS ROUTINES
C
	IFIRST = 1    ! must be 1 to call TLM/ADJ
	ISCALE = 1    ! if 0 just does dimensioning of vectors
	NCALLS = 0    ! =0 at beginning
c user set
c
c INUM can be as large as #words per tape/(nijk*4+nij*2)
c MAMS2B1 has #words/tape as 700 Mwords
c smaller # for INUM may be more convenient; smaller files
c to acquire, when only looking at a few vectors
	INUM   = 200  ! max number of vecs per file (see LANMAX notes)
c
c ISAVE = 0 when no restart is used
c Keep in mind that LANMAX+1 records are available from a
c previous run of  "LANMAX"  iters
	ISAVE  = 10    ! 0, or vector to restart from 
c
c IOUT controls the mid-run output of the restart files, 
c which is useful for preserving data in the case of a system failure
c during a run.
	IOUT   = 20   ! dispose file every iout iteration
c 
c Don't forget to change the filename corresponding to APTH in 
c the script when moving output files to different mass store
c locations via svs.1 instead
 	APTH   = 'TEMP/J89DRY/' ! path name for run-time output
c end user set
C
	ITYPET = ITYPE 
C
	WRITE (LUNIT1,'(A,I10)') 'IFIRST TO OPEN         = ', IFIRST
	WRITE (LUNIT1,'(A,I10)') 'LUNIT1 FOR OUTPUT      = ', LUNIT1
	WRITE (LUNIT1,'(A,I10)') 'ISCALE FOR DIMENSION   = ', ISCALE
	WRITE (LUNIT1,'(A,I10)') 'NCALLS # CALLS TO CATX = ', NCALLS
	WRITE (LUNIT1,'(A,I10)') 'INUM # VECTORS P. FILE = ', INUM
	WRITE (LUNIT1,'(A,I10)') 'ISAVE # VECTORS SAVED  = ', ISAVE
	WRITE (LUNIT1,'(A,I10)') 'IOUT TO CALL NDISPOSE  = ', IOUT
	WRITE (LUNIT1,'(A)')     'APTH PATHNAME FOR FILE = ', APTH
C
	WRITE (LUNIT1,'(A,I10)') 'NI      = ', NI
	WRITE (LUNIT1,'(A,I10)') 'NJ      = ', NJ
	WRITE (LUNIT1,'(A,I10)') 'NK      = ', NK
	WRITE (LUNIT1,'(A,I10)') 'ITRC    = ', ITRC
	WRITE (LUNIT1,'(A,I10)') 'JTRC    = ', JTRC
	WRITE (LUNIT1,'(A,I10)') 'KTRC    = ', KTRC
        WRITE (LUNIT1,'(A,I10)') 'ITRCQ   = ', ITRCQ
        WRITE (LUNIT1,'(A,I10)') 'JTRCQ   = ', JTRCQ
        WRITE (LUNIT1,'(A,I10)') 'KTRCQ   = ', KTRCQ
	WRITE (LUNIT1,'(A,I10)') 'ISET1   = ', ISET1
	WRITE (LUNIT1,'(A,I10)') 'ISET2   = ', ISET2
	WRITE (LUNIT1,'(A,I10)') 'ISET3   = ', ISET3
	WRITE (LUNIT1,'(A,I10)') 'ISET4   = ', ISET4
	WRITE (LUNIT1,'(A,I10)') 'ISET5   = ', ISET5
	WRITE (LUNIT1,'(A,I10)') 'ISET6   = ', ISET6
	WRITE (LUNIT1,'(A,I10)') 'ISET7   = ', ISET7
        WRITE (LUNIT1,'(A,I10)') 'ISET8   = ', ISET8
        WRITE (LUNIT1,'(A,I10)') 'ISET9   = ', ISET9
	WRITE (LUNIT1,'(A,I10)') 'ITYPE   = ', ITYPE
	WRITE (LUNIT1,'(A,I10)') 'N       = ', N
	WRITE (LUNIT1,'(A,I10)') 'NFULL   = ', NFULL
C
	WRITE (LUNIT1,'(/A/)') 'BEFORE GOING INTO CLANC ... ' 
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
	CALL CLANC (N,LANMAX,LUNIT3,KAPPA,W,NW,
     1              LSTEPS,NEIG,EV,BNDEV,NVEC,INDEV,INDVEC,IFAIL)
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
	WRITE (LUNIT1,100)
100	FORMAT (/2X,'AFTER RETURN FROM CLANC ... '//)
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C	DECIDE WHETHER TO DISPOSE Ooutput96 OR NOT AND ACT 
C	ACCORDINGLY. THIS FILE SHOULD BE DISPOSED HERE IF NEW 
C	ITERATIONS HAVE BEEN COMPUTED DURING CLANC AND TESTING
C	WILL BE PERFORMED BELOW. IN THIS CASE THE DISPOSE HERE IS
C	NECESSARY, BECAUSE EVCHK WILL MESS UP THE FILE Ooutput96 
C	DO NOT DISPOSE THIS FILE HERE, IF HISTORY INFORMATION HAS 
C	BEEN READ ONLY. THIS IF-BLOCK REPLACES THE NECESSITY TO 
C	DISPOSE Ooutput96 AT THE END OF THE SCRIPT. 
C	SO, NO LONGER DISPOSE Ooutput96 AT THE END OF THE SCRIPT. 
C	IF IDISP=1 IS USED, SET FNAME CORRECTLY.
C
	IDISP = 1
	IF ( IDISP.EQ.1 ) THEN
	ITAPE = 96
C
c if LANMAX > number of iters that fit on a tape (#?), 
c     if the last pre-existing  tape is not "full" (containing INUM records)
c        then additional iterations will be combined with the existing
c        iterations, and the old version over-written on the mass store (only)
c     if it's still not full, then the old version, still on the
c        CRAY, must be removed before the next run, so that the current
c        version on the mass store is acquired.
c 
	DO I = 1 , 25
	   NC1 = INUM * (I-1) + 1 
	   NC2 = INUM * I 
	   IF ( NCALLS.GE.NC1 .AND. NCALLS.LE.NC2 ) I1 = I
     	ENDDO
	FNAME = 'result.dat.'//CHAR(96+I1) 
C
	TNAME = APTH//FNAME
	CALL NDISPOSE (TNAME,'O',ITAPE,.FALSE.,365,N,' ',.FALSE.)
c	CALL NDISPOSE (TNAME,'O',ITAPE,.FALSE.,365,N,'no',.FALSE.)
	WRITE (LUNIT1,'(//A,A1,A,A24,A,I7//)')
     1  ' WRITE Ooutput96 TO MASS STORE AS ',TNAME
     &  ,' AT NCALLS..',NCALLS	
c kdr
        ierr = ISHELL ('rm Ooutput96')
        WRITE(LUNIT1,'(//A,2I4//)') 
     &       'removed Ooutput96 after NCALLS,LANMAX = ',NCALLS,LANMAX
	ENDIF
C
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
	WRITE (LUNIT1,'(A,I10)')    'DIMENSION OF PROBLEM N  = ', N
	WRITE (LUNIT1,'(A,I10)')    'ITERS-1 PERFORMD LANMAX = ', LANMAX
	WRITE (LUNIT1,'(A,I10)')    'UNIT F. LANCZOS LUNIT3  = ', LUNIT3
	WRITE (LUNIT1,'(A,E22.15)') 'ACC. WRITE VECTOR KAPPA = ', KAPPA
	WRITE (LUNIT1,'(A,I10)')    'STEPS ACT. TAKEN LSTEPS = ', LSTEPS
	WRITE (LUNIT1,'(A,I10)')    'E-VALUES ACCEPTED NEIG  = ', NEIG
	WRITE (LUNIT1,'(A,I10)')    'E-VECTORS ACCEPTED NVEC = ', NVEC
	WRITE (LUNIT1,'(A,I10)')    'ERROR INDICT. (0) IFAIL = ', IFAIL
C
	WRITE (LUNIT1,201) ( EV     (I) , I = 1 , LANMAX )
	WRITE (LUNIT1,202) ( BNDEV  (I) , I = 1 , LANMAX )
201     FORMAT (/2X,'E-VALUS',5E15.7/40(9X,5E15.7/))
202     FORMAT (/2X,'E-BNDS ',5E15.7/40(9X,5E15.7/))
C
	WRITE (LUNIT1,204) ( INDEV   (I) , I = 1 , LANMAX )
	WRITE (LUNIT1,205) ( INDVEC  (I) , I = 1 , LANMAX )
204     FORMAT (2X,'INDEV  ',5I15/40(9X,5I15/))
205     FORMAT (2X,'INDVEC ',5I15/40(9X,5I15/))
C
	CSUMHEV = 0.0             ! compute sum of eigenvalues
	DO 20 I = 1 , LANMAX
	HEV    = EV (I)
	HBNDEV = BNDEV (I)
	IF ( ABS(HEV)    .LT. 1.0E-90 ) HEV    = 0.0
	IF ( ABS(HBNDEV) .LT. 1.0E-90 ) HBNDEV = 0.0
	CSUMHEV = CSUMHEV + HEV
 	WRITE (LUNIT1,'(2X,I7,4D15.7)') I , HEV , HBNDEV , CSUMHEV
20	CONTINUE
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C	IF THE CALL TO CLANC WAS SUCCESSFUL AND EIGENVECTORS
C	HAVE BEEN COMPUTED THEN GO THROUGH THE NEXT IF BLOCK
C	FOR RETRIEVAL AND STORAGE OF EIGENVECTORS. 
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
	IF ( IFAIL.EQ.0 .AND. LUNIT3.GT.0 ) THEN
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
	DO 30  I  = 1 , LANMAX
	EV    (I) = 0.0
30	BNDEV (I) = 0.0
C
c user set
C
C	RETRIEVE SELECTED EIGENVECTORS WRITTEN DURING CLANC
C       L2-L1+1 MUST BE <= JEVEC;  L2 should be <= # converged vects
c       (see README note about convergence)
C
 	L1    = 1
 	L2    = MIN ( NVEC , JEVEC )
C       L1    = 110
C       L2    = 129
c end user set
	IEVEC = N
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
	CALL VRETR ( L1,L2,LUNIT3,NVEC,N,
     1               EV,BNDEV,EVEC,IEVEC,JEVEC,IFAIL )
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
	WRITE (LUNIT1,105)
105	FORMAT (/2X,'AFTER RETURN FROM VRETR ... '//)
C
	WRITE (LUNIT1,'(A,I10)')    'FIRST E-VEC RETR.  L1   = ', L1
	WRITE (LUNIT1,'(A,I10)')    'LAST  E-VEC RETR.  L2   = ', L2
	WRITE (LUNIT1,'(A,I10)')    'UNIT F. LANCZOS LUNIT3  = ', LUNIT3
	WRITE (LUNIT1,'(A,I10)')    'E-VECTORS ACCEPTED NVEC = ', NVEC
	WRITE (LUNIT1,'(A,I10)')    'DIMENSION OF PROBLEM N  = ', N
	WRITE (LUNIT1,'(A,I10)')    '1ST DIM. OF EVEC IEVEC  = ', IEVEC
	WRITE (LUNIT1,'(A,I10)')    '2ND DIM. OF EVEC JEVEC  = ', JEVEC
	WRITE (LUNIT1,'(A,I10)')    'ERROR INDICT. (0) IFAIL = ', IFAIL
C
	WRITE (LUNIT1,201) ( EV     (I) , I = 1 , LANMAX )
	WRITE (LUNIT1,202) ( BNDEV  (I) , I = 1 , LANMAX )
C
C	WRITE EIGENVALUES, BOUNDS, VECTORS TO LUNIT2
C
	OPEN (UNIT=LUNIT2,FILE='vector.dat',FORM='UNFORMATTED',
     1        ACCESS='DIRECT',RECL= NP2 * 8 )
 	DO 10 J = 1 , L2-L1+1
10	WRITE (LUNIT2,REC=J)  EV(J), BNDEV(J),( EVEC (I,J),I=1,N )
C
	ENDIF
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C	SO FAR GOES THE USUAL LANCZOS RUN WITH COMPUTATION AND 
C	STORAGE OF EIGENVECTORS. THE FOLLOWING BLOCK IS FOR 
C	DIMENSIONING AND TESTING OF THE RESULTS. 
C	
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C	DO THE DIMENSIONING ONLY IF CONDITIONS ARE MET 
C
	IF ( IFAIL.EQ.0 .AND. MDIM.GE.1 ) THEN 
C
C	OPEN THE DIMENSIONED PLOT FILE 
C
	OPEN (UNIT=LUNIT6,FILE='vectordim.dat',FORM='UNFORMATTED',
     1        ACCESS='DIRECT',RECL= NFP2 * 8 )
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C	READ EIGENVECTORS FROM LUNIT2, PUT DIMENSIONS BACK ON THEM 
C	AND WRITE THEM OUT TO LUNIT6. THIS IS DONE BY A CALL TO 
C	ATX93 WITH ISCALE SET TO ZERO. IN THIS CASE ATX93 IS LEFT
C	RIGHT AFTER DIMENSIONS HAVE BEEN PUT ON THE EIGENVECTOR. 
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
	WRITE (LUNIT1,50)
50	FORMAT (/2X,'BEFORE PUTTING DIMENSIONS ON EIGENVECTORS ... '//)
C
C	READ EIGENVECTORS, PUT DIMENSIONS BACK ON THEM AND WRITE THEM
C
	MDIM1 = MIN ( NVEC , MDIM )
 	MDIM1 = L2 - L1 + 1
C
	DO 111 J = 1 , MDIM1
	   READ  (LUNIT2,REC=J) X2(1),X2(2),(X1(I),I=1,N)
C
	   IFIRST = 1
	   ISCALE = 0  ! do scaling only 
	   NCALLS = 0
	   INUM   = 200
	   ISAVE  = 0
	   IOUT   = 800
C
	   CALL ATX93 ( X1 , N ) 
 	   WRITE (LUNIT6,REC=J) X2(1),X2(2),(XFULL(I),I=1,NFULL)
C
	   WRITE (LUNIT1,60) J 
60	   FORMAT (/2X,'DIMENSIONING OF EIGENVECTOR DONE ... ',I10//)
C
111	CONTINUE
C
	ELSE
C
	WRITE (LUNIT1,65) MDIM 
65	FORMAT (/2X,'DIMENSIONING OF EIGENVECTORS NOT DONE ... ',I10//)
C
	ENDIF	
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C	NOW COMES THE TESTING OF THE DIMENSIONLESS EIGENVECTORS. 
C	TESTING IS ONLY PERFORMED IF MTEST .GE. 1, OTHERWISE THE 
C	ROUTINE EVCHK IS LEFT IMMEDIATELY. 
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C
	IF ( IFAIL.EQ.0 .AND. MTEST.GE.1 ) THEN 
C
 	   WRITE (LUNIT1,305) MTEST 
 305 	   FORMAT (//2X,'RESULTS FROM TESTING THE EIGEN-PAIRS ...'/
     1             2X,'WITH MTEST = ',I5/)
C
C	SET ISCALE TO 1 TO HAVE ATX93 DO ANYTHING
C
	   IFIRST = 1
	   ISCALE = 1
	   NCALLS = 0
	   INUM   = 200
	   ISAVE  = 0
	   IOUT   = 100
C
	   ACC = 1.0E-10
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	   CALL EVCHK ( N,MTEST,LUNIT1,LUNIT2,ACC,
     1                 X1,X2,X3,X4,W,NW,IFAIL )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
	ELSE
C
	WRITE ( LUNIT1,310 ) 
310	FORMAT (//2X,'NO TESTING OF EIGENVECTORS PERFORMED ...'//)
C
	ENDIF
C	
 	WRITE (LUNIT1,'(A,I10)')    'DIMENSION OF PROBLEM N  = ', N
 	WRITE (LUNIT1,'(A,I10)')    '# VECTORS TESTED  MTEST = ', MTEST
 	WRITE (LUNIT1,'(A,I10)')    'LUNIT1  (OUTPUT)        = ', LUNIT1
 	WRITE (LUNIT1,'(A,I10)')    'LUNIT2  (EIGENVECTORS)  = ', LUNIT2
 	WRITE (LUNIT1,'(A,E22.15)') 'ACCURACY CRIT.   (ACC)  = ', ACC
 	WRITE (LUNIT1,'(A,I10)')    'DIMENSION OF WORK (NW)  = ', NW
 	WRITE (LUNIT1,'(A,I10)')    'ERROR INDICT. (0) IFAIL = ', IFAIL
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C	END OF TESTING. REMEMBER TO DISPOSE FILES.
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
	STOP
	END
C
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C	END OF MAIN PROGRAM TO DRIVE LANCZOS ROUTINES AND 
C	TO PRODUCE PLOT FILE. 14 DECEMBER 1993. 
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C
C
        SUBROUTINE ATX93 ( X , N )
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C	THIS ROUTINE IS DESIGNED TO SERVE AS AN INTERFACE 
C	BETWEEN THE LANCZOS ROUTINES AND THE ACTUAL 
C	COMPUTATIONS FOR THE EVR-TLM-ADJOINT CODE. 
C
C	IN PARTICULAR, THIS ROUTINE WILL - IF TOLD APPROPRIATELY - 
C	ACCESS INFORMATION FROM PREVIOUS COMPUTATIONS. THIS INFORMATION
C	IS SUPPLIED BY THE COMMON BLOCK IFIRST.
C
C	THE BASIC IDEA IS THAT X IS TRANSFORMED EITHER BY 
C	READING THE TRANSFORM FROM RESULTS OF A PREVIOUS RUN (THIS
C	IS CHEAP), OR TO COMPUTE THE SOLUTION BY A CALL TO CATX
C	(THIS IS EXPENSIVE, BUT THE RESULTS WILL BE SAVED). 
C
C	M. EHRENDORFER, 30 SEPTEMBER 1993 AND 25 OCTOBER 1993. 
C
C	MODIFICATIONS, 1 JULY 1994
c
c	adapted for moist r-norm, 12/1/94, M.E. 
C	q-only norm, 12/20/94, m.e.
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
c user set
        PARAMETER ( NI = 43 , NJ = 54 , NK = 10 )
c end user set
	PARAMETER ( NI2 = NI-2, NJ2 = NJ-2 )
	PARAMETER ( NI3 = NI-3, NJ3 = NJ-3 ) 
C
	PARAMETER ( NFULL = 2*NK*NI2*NJ2 + (NK+1)*NI3*NJ3 
     1              + NK*NI3*NJ3 )
C
	REAL X ( N ) , XFULL ( NFULL ) 
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C	THE FOLLOWING COMMON BLOCK IS USED TO CONNECT THIS ROUTINE 
C	WITH THE LANCZOS ROUTINES. ITS ARGUMENTS ARE: 
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C	IFIRST  ...  TO CONTROL THE ACQUIRE OF THE BASIC STATE.
C		     IF SET TO ZERO, ACQUIRE + OPEN BASIC STATE, RETURN. 
C		     IF UNEQUAL ZERO, AQUIRE + OPEN SKIPPED, EXECUTION. 
C
C	LUNIT1  ...  AN OUTPUT UNIT TO WHICH THE CHECK OF THE INNER
C		     PRODUCTS AND OTHER INFORMATION IS WRITTEN. 
C
C	ISCALE  ...  TO CONTROL THE USE OF THIS ROUTINE FOR SIMPLE
C		     DIMENSIONING OF THE FIELDS. 
C		     IF SET TO ZERO THEN A DIMENSIONLESS VECTOR XSTAR
C		     IS DIMENSIONED, AND A RETURN OCCURS. 
C   		     IF UNEQUAL ZERO THE DIMENSIONING BLOCK IS SKIPPED. 
C
C	NCALLS  ...  COUNTER FOR THE NUMBER OF CALLS
C
C		     NCALL = 1 ATX93 IS CALLED FOR THE FIRST TIME
C		     NCALL = 2 ATX93 IS CALLED FOR THE SECOND TIME
C		     NCALL = 3 ATX93 IS CALLED FOR THE THIRD TIME 
C	     	     ETC. 
C		     THIS COUNTING DOES NOT INCLUDE THE CALL 
C		     TO OPEN THE FILES ( IFIRST = 0 ) OR THE 
C		     DIMENSIONING OF THE STATE VECTOR ( ISCALE = 0 ).
C
C	INUM    ...  THE TOTAL NUMBER OF VECTORS TO BE WRITTEN TO 
C		     A SINGLE FILE. THIS NUMBER MIGHT USUALLY BE 200.
C
C	ISAVE   ...  THE NUMBER OF VECTORS RELEVANT FOR THIS 
C		     RUN THAT HAVE BEEN SAVED AND SHOULD BE ACCESSED. 
C		     SET ISAVE TO ZERO, IF NO INFORMATION IS 
C		     AVAILABLE. 
C
C	IOUT    ...  IF NCALLS IS DIVISIBLE BY IOUT, AND NEW 
C		     INFORMATION HAS BEEN COMPUTED THE FILE IS
C		     WRITTEN TO MASS-STORE. IOUT MUST BE LE INUM, 
C		     TO MAKE SURE THAT A FULL FILE IS ALWAYS 
C		     DISPOSED DURING THE RUN. FOR INUM = 200, 
C		     SET IOUT TO 20, FOR EXAMPLE. 
C
C	ITYPE   ...  INDICATOR FOR NORM TO USE, THIS IS THE SUM OF ISETS
C
C	ISET1 = 1 TO USE THE TOTAL ENERGY NORM ROUTINES
C	ISET2 = 1 TO USE OPTION 1B OF ENERGY NORM (PROJECT AT T1)
C	ISET3 = 1 TO USE OPTION 2A OF ENERGY NORM (PROJECT AT T0, T1)
C	ISET4 = 1 TO USE POTENTIAL ENSTROPHY NORM (RON)
C	ISET5 = 1 TO USE THE SET OF TOMIS NORM ROUTINES
C	ISET6 = 1 TO USE MOIST R-NORM (MARTIN AND RON) 
C	ISET7 = 1 TO USE Q ONLY NORM (MARTIN AND RON)
C       ISET8 = 1 TO USE TE-NORM PLUS Q
C
C	APTH    ...  THE PATH NAME TO BE PUT BEFORE THE FILES
C		     WHEN THEY ARE WRITTEN OUT TO MASS-STORE. 
C		     
C	XFULL   ...  A VECTOR NEEDED TO TRANSFER THE DIMENSIONED
C		     VECTOR OUT TO THE CALLING ROUTINE.
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
	CHARACTER*12  FNAME , APTH 
	CHARACTER* 9  FNAM1
	CHARACTER*128 TNAME
C
	COMMON /IFIRST/ 
     1  IFIRST,LUNIT1,ISCALE,NCALLS,INUM,ISAVE,IOUT,ITYPE,APTH
	COMMON /ISECND/ XFULL 
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C	DO THE OPENING OF THE FILES 
C
	IF ( IFIRST.EQ.0 ) THEN 
	CALL CATX (X,N,XFULL,NFULL,IFIRST,LUNIT1,ISCALE,ITYPE)
	WRITE (LUNIT1,'(/A/)') 'OPENING OF FILES SUCCESSFUL ...'
	RETURN
	ENDIF
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C	DO THE DIMENSIONING OF THE STATE VECTOR IF REQUESTED
C
	IF ( ISCALE.EQ.0 ) THEN 
	CALL CATX (X,N,XFULL,NFULL,IFIRST,LUNIT1,ISCALE,ITYPE)
	WRITE (LUNIT1,'(/A/)') 'DIMENSIONING OF VECTOR SUCCESSFUL ...'
	RETURN
	ENDIF
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
	NCALLS = NCALLS + 1
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c kdr
        IF (NCALLS.GT.1 .AND. MOD(NCALLS,INUM).EQ.1) THEN
           ierr = ISHELL ('rm Ooutput96')
           WRITE(LUNIT1,'(//A,2I4//)') 
     &          'removed Ooutput96 after NCALLS,INUM = ',NCALLS,INUM
        ENDIF
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
        IF ( NCALLS .LE. ISAVE ) THEN
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C	IN THIS CASE INFORMATION IS AVAILABLE FOR THIS RUN AND
C	FOR THIS PARTICULAR CALL (IDENTIFIED BY NCALLS) TO ATX93. 
C
C	DETERMINE THE NAME OF THE HISTORY FILE FROM THE VALUE 
C	OF NCALLS. RESTRICTION: AT MOST 25 HISTORY FILES. 
C	THE IDEA IS THAT 
C	THE FIRST  INUM VECTORS ARE ON THE FILE result.dat.a, 
C	THE SECOND INUM VECTORS ARE ON THE FILE result.dat.b,
C	ETC. FOR INUM = 200, THIS WILL WORK UP TO NCALLS = 5000. 
C
C	NOTE THE FILES result.dat.a, result.dat.b, MUST BE AQUIRED
C	BEFORE EXECUTION, AND MUST BE IN TMPDIR, IF ISAVE.GT.ZERO. 
C	
	DO 7010 I = 1 , 25
	NC1 = INUM * (I-1) + 1
	NC2 = INUM * I
	IF ( NCALLS.GE.NC1 .AND. NCALLS.LE.NC2 ) I1 = I
7010 	CONTINUE
C
C	READ THE INFORMATION FROM THE CORRECT HISTORY FILE
C
	FNAME = 'result.dat.'//CHAR(96+I1) 
	OPEN(95,FILE=FNAME,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=N*8)
	IREC = NCALLS - (I1-1) * INUM 	
	READ (95,REC=IREC) ( X (I) , I = 1 , N )
	CLOSE (95) 
	WRITE (LUNIT1,'(/A,I5,A,A12,A,I7/)')  
     1	'NCALLS = ',NCALLS,' READ FROM FILE ',FNAME,' RECORD = ',IREC
C
C	WRITE THE INFORMATION JUST READ TO THE FILE Ooutput96
C	THIS IS DONE BECAUSE IF THE CURRENT HISTORY FILE IS NOT FULL, 
C	Ooutput96 IS DISPOSED LATER BY NDISPOSE UNDER THE NAME 
C	OF THE CURRENT HISTORY FILE. 
C
	FNAM1 = 'Ooutput96'
	OPEN(96,FILE=FNAM1,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=N*8)
	WRITE (96,REC=IREC) ( X (I) , I = 1 , N )
	CLOSE (96) 
	WRITE (LUNIT1,'(A,I5,A,A12,A,I7/)')  
     1	'NCALLS = ',NCALLS,' WRITE TO FILE ',FNAM1,' RECORD = ',IREC
C
C	LEAVE THE ROUTINE WITH THE CORRECT ANSWER
C
	RETURN
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C	OTHERWISE THE RESULT MUST BE COMPUTED BY A CALL TO CATX.
C	IT WILL BE SAVED FOR LATER USE.
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
	ELSE 
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C	MAKE THE EXPENSIVE CALL 
C
C kdr        PRINT*,'      CATX called for NCALLS = ',NCALLS
	CALL CATX (X,N,XFULL,NFULL,IFIRST,LUNIT1,ISCALE,ITYPE)
C
	WRITE (LUNIT1,'(/A,I5,A)')  
     1	'NCALLS = ',NCALLS,' SUCCESSFUL COMPUTATION IN CATX '
C
C	DETERMINE THE NAME OF THE HISTORY FILE TO WRITE TO
C
	DO 7020 I = 1 , 25
	NC1 = INUM * (I-1) + 1
	NC2 = INUM * I
	IF ( NCALLS.GE.NC1 .AND. NCALLS.LE.NC2 ) I1 = I
7020 	CONTINUE
C
	IREC = NCALLS - (I1-1) * INUM 	
	FNAME = 'result.dat.'//CHAR(96+I1) 
C
C	WRITE THE RESULT OF THE COMPUTATION TO Ooutput96
C
	FNAM1 = 'Ooutput96'
	OPEN(96,FILE=FNAM1,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=N*8)
	WRITE (96,REC=IREC) ( X (I) , I = 1 , N )
	CLOSE (96) 
	WRITE (LUNIT1,'(A,I5,A,A12,A,I7/)')  
     1	'NCALLS = ',NCALLS,' WRITE TO FILE ',FNAM1,' RECORD = ',IREC
c PROBLEM;
c     writes to second tape are written to the same UNIT as before,
c     overwriting the first tape iters.  All the remaining first
c     tape iters remain in this file, and are disposed with it!
c FIX; remove Ooutput96 after tape is filled and disposed
C
C	DISPOSE Ooutput96 AS THE CURRENT HISTORY FILE TO MASS STORE
C	ALSO DISPOSE THE PROTOCOL FILE Ooutput91 TO MASS STORE
C
	IF ( MOD (NCALLS,IOUT) .EQ. 0 )  THEN 
C
	ITAPE = 96
	TNAME = APTH//FNAME
c 
c	CALL NDISPOSE (TNAME,'O',ITAPE,.FALSE.,365,N,'no',.FALSE.)
	CALL NDISPOSE (TNAME,'O',ITAPE,.FALSE.,365,N,' ',.FALSE.)
	WRITE (LUNIT1,'(//I7,A,A24//)') 
     1     NCALLS-ISAVE,' RECORDS WRITTEN TO MASS STORE  ',APTH//FNAME
C
c I don't know why NDISPOSE has been partially commented out
c I'll comment out the rest so that it's not misleading.
        ier = ISHELL('cp lanc.out ${HOME}/j89.out')
C	ITAPE = LUNIT1
CC	TNAME = APTH//'lanc.dat'
C	TNAME = APTH//'Ooutput91'
C	CALL NDISPOSE (TNAME,'O',ITAPE,.FALSE.,365,N,' ',.FALSE.)
C	WRITE (LUNIT1,'(//A,A24,A,I7//)') 
C     1  ' WRITE TO MASS STORE  ',TNAME,' AT NCALLS ... ',NCALLS
C	
	ENDIF
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
	ENDIF
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
	RETURN
	END
C
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C   inserted the routine catx from tlmadj.8, 6 July 1994
C
C
c
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C
	SUBROUTINE CATX (XSTAR,NXSTAR,XFULL,NFULL,IFIRST,LUNIT1,ISCALE,
     1                   ITYPE)
C
C
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This is a a combination of the Tangent Linear Model (TLM) and 
c     Adjoint Model (ADJ) component of Version 2 of
c     the Mesoscale Adjoint Modeling System (MAMS2).
c
C     THIS ROUTINE IS DESIGNED TO PERFORM THE MATRIX VECTOR
C     MULTIPLICATION NEEDED IN COURSE OF THE LANCZOS ITERATION.
C
C     IT IS ALWAYS CALLED FROM THE HIGHER LEVEL ROUTINE
C     ATX93.
C
C     M. Ehrendorfer. September 1993, December 1993, July 1994. 
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This is the driver for the Tangent Linear Model (TLM) component
c     of version 1D of the Mesoscale Adjoint Modeling System (MAMS1).
c
c     This code (vesion 1D) was developed by Ronald Errico
c     at NCAR (January 1997).  Prospective users should inform
c     Ronald Errico by e-mail or in writing of any intentions to use
c     this software and agree to reference any appropriate
c     publications of the developers and acknowledge use of this
c     software in any papers.  Furthermore, pre-prints and reprints of
c     any papers produced using this software should be sent to R.
c     Errico at NCAR.
c
c     For a description of MAMS version 1A see NCAR Technical Note
c     #TN/410+IA, 1994, by R. M. Errico, K. Raeder, and T. Vukicevic.
c     Version 2 is very differeant than 1A.  For a summary of the 
c     changes, contact R. Errico.
c
c---------------------------------------------------------------------
c
      logical lmoist,lground,lrain,lcloud,lsnow,lbdry,lbdryx,lsplit
     &,       lbdryp,lbdsame,llinear,lcoefout,linitial,linitrd,llast
     &,       ldrycon
c 
c---------------------------------------------------------------------
c                  USER DEFINED PARAMETERS :
c---------------------------------------------------------------------
c
c
c  grid dimensions
c user set
       parameter ( ni=43, nj=54, nk=10, nspong=5)
       
c
c  timing information
       parameter (dt=60.*4., timax=24.*60., dtbdr=60.*12.)
c
c  output file information
      parameter (toutreg=24.*60, ipackw=1)
c
c  basic state input file information (units=mins)
      parameter (ftimed0X=0., tupdated=12.)
      parameter (ftimem0X=0., tupdatem=12.)
c
c  moisture (actually Jacobian; don't forget dry conv adj) selection
      parameter (lmoist=.t., nqpoints=1000)
c
c  select dry convective adjustment
      parameter (ldrycon=.t.)
c
c  ground temperature physics selection
      parameter (lground=.t.,  lcloud=.t., lsnow=.f.)
c
c  consider boundary perturbations?
      parameter (lbdryp=.f., lbdsame=.f., ftimeb0=0.)
c
c  selection of initial conditions in routine = Ainit
      parameter (ipert=3,  pert=1.0,  tskips1=0., tskips2=0. )
c
c  specification of implicit nonlinear normal mode initialization
c    linitrd to be set to true if basic state of initialization is to
c    be read from a separate file and varies with iteration number.
c       The rule to follow is this: If you want to do an initialization 
c       and (1) an initial file was created along with the basic state 
c       file, you must set linitrd=.t.; (2) an initial file was not 
c       created along with the basic state file, you must set 
c       linitrd=.f. If you set linitial=.f., it doesn't matter what 
c       you set for linitrd.
c Essentially, if an `I' file was created during the run, the basic
c state file will signal the code to read from it when the
c initialization is performed, so it must be acquired along with the
c basic state file.  Run in its normal fashion (outside your code),
c these decisions and the `I' acquire will be done automatically.
c
       parameter (linitial=.t., linitrd=.t.)
c end user set
c
c---------------------------------------------------------------------
c                END OF USER DEFINED PARAMETERS
c---------------------------------------------------------------------
c
c                     OTHER PARAMETERS
c
      parameter (nsplitmx=3, npkappa=1400, nqsatvp=844)
c
c                     COMPUTED PARAMETERS
c
      parameter (dth=dt/3600., dtm=dt/60., dt2=2.*dt )
      parameter (nkp=nk+1, nkp1=nkp, nkm1=nk-1, nkpx=nk+6
     &,          nij=ni*nj, nijk=nij*nk, nij2=(ni+nj)*2, nkp2=nk+2
     &,          nij2n2=nij2*nspong*2, nij2n2k=nij2n2*nk)
c
c  parameters for split explicit mode determination
      parameter (nijs=nij*nsplitmx, nix=nk+4, njx=5, nijx=nix*njx
     &,          nijkx=nijx*nkp2, nkmnsplt=nk-nsplitmx, llinear=.t.)
c
c  parameters for precipitation and dry convection schemes
      parameter (nq1dim=nqpoints, nq2dim=nq1dim*nk
     &,          nq3dim=nq2dim*nk, nqdimw=(nq2dim+2)*2+(nq3dim+2)*4
     &,          nq4dim=nq2dim*2+nq3dim*4, nq5dim=nq2dim*3
     &,          nq6dim=nq2dim+nq3dim)
c
c  parameters set for TLM
      parameter (lcoefout=.f.)  
c
c
c 
c               VARIABLE AND PARAMETER LIST:
c
c parameters subject to user change :
c
c     ni ............ south-north number of grid points
c     nj ............ west-east number of grid pts.
c     nk ............ vertical (top to bottom) no. of grid pts.
c     nspong..........width of sponge domain for lateral bounds
c     dt ............ time step (sec)
c     timax ......... maximum integration time (minutes)
c     dtbdr ..........time interval of boundary value definition
c                     (only necessary if boundaries are perturbed)
c     toutreg ....... if .gt. 0. then periodic output
c     ipackw ........ packing density for output
c     ftimed0 ....... time (mins) to skip on b.s. data file
c     tupdated ...... periodic update interval for b.s. data file
c     ftimem0 ....... time (mins) to skip on moist coefficient files
c     tupdatem ...... periodic update interval for moist coeff files
c
c     lbdryp .........true if lateral boundaries are perturbed
c     lbdsame.........true if initial and lateral boundary pert files same
c     ftimeb0 ....... time (mins) to start on lateral bound pert file
c     lmoist .........true if moisture considered
c     lground ........true if ground temperature considered
c     lcloud .........true if cloud effects on radiative impact on
c                     ground temperature is to be considered 
c     ldrycon.........perform dry convective adjustment
c     nqpoints........must be greater than or equal to the maximum
c                     number of moist-convecting i,j points (for
c                     storing Jacobians).  Also must be greater than
c                     the number of dry convection points.  Also, the
c                     number of nonconvective precipitating points
c                     (counting  each k-level separately) must not be
c                     greater than nk*this.  If any of these limits
c                     are exceeded, execution will stop with an
c                     explanation.  These limits do not apply if the
c                     physical process is turned off (set to .false.).
c     ipert ..........selects method for definition of initial field:
c           (1) gravity wave perturbation; (2) random (all fields);
c           (3) pt at center, bottom; (4) read and differ 2 input files;
c           (5) read from input file; (6) ps at center of grid
c     ipert ..........selects method for definition of initial field
c     pert ...........scaling factor for initial condition
c     tskips1.........selects time on file to be used if initial
c                     conditions to be read from a file on iuniti(1)
c     tskips2.........selects time on file to be used if initial
c                     conditions to be determined as differences
c                     between values on two files iuniti(1)-iuniti(2)
c     linitial....... .t. if adjopint of implicit nonlinear normal 
c                     mode initialization to be done 
c
c  The terms half and full levels respectively refer to data and 
c     interface levels on the vertical grid.
c  The terms dot and cross grids refer to the horizontal grids on
c     which pu, pv and pt, pq, ps, tg are defined, repectively
c
c  Principal dynamic fields (gradients)
c
c  pu......... eastward wind coupled to ps at current time
c  pv......... northward wind coupled to ps at current time
c  pt......... temperature coupled to ps at current time
c  ps......... surface pressure minus ptop at current t (units=cb)
c  pq......... specific humidity coupled to ps at current time
c  tg......... surface temperature
c  prefix 'p' absent denotes decoupled field
c  prefix 'bdr' refers to lateral boundary relaxation fields
c  suffix 'm' indicates values at previous time step
c  suffix 'n' indicates values at time step being forecast
c  suffix 'ten' indicates time tendency of field
c  suffix 'x' denotes dot grid values interpolated to cross grid
c  suffix 'd' denotes cross grid values interpolated to dot grid
c  suffix 'f' on pu or pv denotes values divided by map scale factor
c  a final suffix 'b' indicates a basic state (i.e. reference) field
c
c  Derived fields  
c 
c  omega  ..omega (dp/dt) field
c  sdot  ...sigma-dot (d(sig)/dt) field
c  xflux ...factor for surface and vertical diffusion fluxes
c  dzab ....thickness for vertical diffusion flux
c  div .....horizontal mass flux divergence
c  pkf .....factor for converting t to theta on full levels
c  pkh .....factor for converting t to theta on half levels
c  phi .....hydrostatic geopotential
c  chi .....vertical intergral of div
c  cpmr ....recipricol of Cp for moist air
c  gsw......net short-wave radiative flux at surface
c  glw......downward long-wave radiative flux at surface
c  cldfr ...cloud fractional cover (in layers)
c  ptv......virtual temperature coupled to ps
c  qs.......saturation mixing ratio
c  prw......precipitable water
c  raincv...convective precipitation (units cm)
c  rainncv...non convective precipitation (units cm)
c  qcoefs1, etc are arrays for storing 
c     the Jacobian matrix for the precipitation processes
c
c  Arrays for split explicit scheme and initialization
c
c  edepth...... the vector of values of grav * equivalent depths
c  einverse.... the inverse of ematrixi
c  ematrixi.... matrix whose columns are the vertical eigenvectors
c  gmatrixt.... the hydrostatic matrix for t perturbations
c  gmatrixp.... the hdrostatic vector for ps perturbations
c  pmatrix..... the vector describing response of ps to div perts.
c  tmatrix..... the matrix describing response of pt to div perts.
c  tpmeanb..... the means of t on sigma surfaces and the mean of ps
c  msteps...... the numbers of fine steps for each equivalent depth 
c 
c  Arrays describing the physical nature of the grid
c
c  cf .................... Coriolis coefficient
c  topog ................. basic state topography
c  fmapx.................. cross-point mapscale factor
c  fmapd.................. dot-point mapscale factor
c  fmapx2................. square of fmapx
c  fmapd2................. square of fmapd
c  sigf .................. full-sigma level values 
c  sigh .................. half-sigma level values
c  dsig .................. veritical grid spacing
c  xkc ................... horizontal diffuson coefficient
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c                     DIMENSION STATEMENTS
c
ctlmadj
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c	define the vector that is transformed by this routine
c
        REAL XSTAR ( NXSTAR ) , XFULL ( NFULL ) 
c
c	fields for tomi's routines ( itype = 5 )
c
cx	character*2 normtype
cx	real ynorm (ni,nj,nk,4)
cx	real w1 (nk), w2 (nk), dsigf (nk), works (ni,nj,nk), dsigh (nk)
c
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
c dimension grid arrays (see routine = Nsetup for definitions)
c
      real fmapx(ni,nj), fmapd(ni,nj), fmapx2(ni,nj)
     &,    fmapd2(ni,nj), cf(ni,nj), xland(ni,nj), topog(ni,nj)
     &,    sfctyp(ni,nj), snowc(ni,nj), tdeepg(ni,nj)
     &,    sigh(nk), sigf(nkp), dsig(nk)
     &,    xlat(ni,nj), xlon(ni,nj), dlat(ni,nj), dlon(ni,nj)
     &,    galbedo(ni,nj), gthrmin(ni,nj), gemissv(ni,nj)
     &,    groughl(ni,nj), gmavail(ni,nj), gheatcpr(ni,nj)
     &,    cdragx(ni,nj), cdragd(ni,nj), pkappa(npkappa)
     &,    bfact(nspong)
      integer iland(nj), klevcld(4)
c
c
c diffusion and radiation coefficients
c
      dimension xfluxb(ni,nj,nkpx), xflux(ni,nj,nkpx), xkc(ni,nj)
     &,         glw(ni,nj), gsw(ni,nj)
c
c dimension basic state variables
c
      real pub(ni,nj,nk), pufb(ni,nj,nk), pufxb(ni,nj,nk)
     &,    pvb(ni,nj,nk), pvfb(ni,nj,nk), pvfxb(ni,nj,nk)
     &,    ub(ni,nj,nk), uxb(ni,nj,nk), vb(ni,nj,nk), vxb(ni,nj,nk)
     &,    ptb(ni,nj,nk), tb(ni,nj,nk), divb(ni,nj,nk)
     &,    psb(ni,nj), pdb(ni,nj)
     &,    sdotb(ni,nj,nkp1), sdotdb(ni,nj,nkp1), omegab(ni,nj,nk)
     &,    phib(ni,nj,nk), cpmrb(ni,nj,nk), tdb(ni,nj,nk)
     &,    pkhb(ni,nj,nkp1), tgb(ni,nj), tvb(ni,nj,nk)
     &,    alogpb(ni,nj,nkp1), ptvb(ni,nj,nk), thetab(ni,nj,nkp1)
      real pumb(ni,nj,nk), pvmb(ni,nj,nk), ptmb(ni,nj,nk), psmb(ni,nj)
     &,    umb(ni,nj,nk), vmb(ni,nj,nk), tmb(ni,nj,nk), pdmb(ni,nj)
     &,    tgmb(ni,nj)
c
c dimension perturbation variables 
c
      real pu(ni,nj,nk), puf(ni,nj,nk), pufx(ni,nj,nk), puten(ni,nj,nk)
     &,    pv(ni,nj,nk), pvf(ni,nj,nk), pvfx(ni,nj,nk), pvten(ni,nj,nk)
     &,    u(ni,nj,nk), ux(ni,nj,nk), v(ni,nj,nk), vx(ni,nj,nk)
     &,    pt(ni,nj,nk), ptten(ni,nj,nk), t(ni,nj,nk)
     &,    ps(ni,nj), psten(ni,nj), pd(ni,nj), div(ni,nj,nk)
     &,    sdot(ni,nj,nkp1), sdotd(ni,nj,nkp1), omega(ni,nj,nk)
     &,    phi(ni,nj,nk), tg(ni,nj), tgten(ni,nj), cpmr(ni,nj,nk)
     &,    ptv(ni,nj,nk), chi(ni,nj), tv(ni,nj,nk), td(ni,nj,nk)
     &,    alogp(ni,nj,nkp1), theta(ni,nj,nkp1), xpprp(ni,nj,nkp1)
      real pum(ni,nj,nk), pvm(ni,nj,nk), ptm(ni,nj,nk), psm(ni,nj)
     &,    um(ni,nj,nk), vm(ni,nj,nk), tm(ni,nj,nk), pdm(ni,nj)
     &,    pun(ni,nj,nk), pvn(ni,nj,nk), ptn(ni,nj,nk), psn(ni,nj)
     &,    tgm(ni,nj), tgn(ni,nj), field3c(nijk,5), field2c(nij,2)
c
c  The following are equivalenced for proper use of the file
c  differencing option in routine = Finit 
      equivalence (field3c,     pum), (field3c(1,2),pvm)
     &,           (field3c(1,3),ptm), (field2c,     psm, field3c(1,5))
     &,           (field2c(1,2),tgm)  
c 
c dimension moisture fields 
c
      real prwb(ni,nj), prw(ni,nj), cldfr(ni,nj,3), cldfrb(ni,nj,3)
     &,    qsatvp(nqsatvp)
      real pq(ni,nj,nk), pqm(ni,nj,nk), pqn(ni,nj,nk), pqten(ni,nj,nk)
     &,    q(ni,nj,nk), qm(ni,nj,nk), pqb(ni,nj,nk), pqmb(ni,nj,nk)
     &,    qb(ni,nj,nk), qmb(ni,nj,nk), qs(ni,nj,nk), qsb(ni,nj,nk)
     &,    raincv(ni,nj), rainncv(ni,nj), raincv2(ni,nj)
     &,    rainncv2(ni,nj), qcoefs1(nq4dim), qcoefs3(nq2dim,3)
     &,    workq(nqdimw), qcoefsw(nq4dim)
      integer ipnoncv(ni,nj,nk), ipconvt(ni,nj)
      data raincv/nij*0./, rainncv/nij*0./
      data raincv2/nij*0./, rainncv2/nij*0./
      equivalence (field3c(1,4),pqm)
c
c  The following are for the moist coefficient file headers
      parameter (nnvarsq=30, nrvarsq=10)
      real rvarsq1(nrvarsq), rvarsq2(nrvarsq)
      integer nvarsq1(nnvarsq), nvarsq2(nnvarsq)
c
c  The following are for dry convection coefficients and files
      real dcoefs1(nq6dim), rvarsa1(nrvarsq)
      integer ipdconv(ni,nj), nvarsa1(nnvarsq)
c
c  The following are equivalenced because neither the TLM or ADJ codes
c  properly treat multiple times of the basic state in the time 
c  stepping procedure
c
      equivalence (umb,ub), (vmb,vb), (tmb,tb), (pdmb,pdb)
     &,           (qmb,qb)
c
c dimension arrays for split-explicit and initialization scheme
      real ematrix(nk,nk), einverse(nk,nk), gmatrixt(nk,nk)
     &,    gmatrixp(nk), tmatrix(nk,nk), pmatrix(nk), edepth(nk)
     &,    tpmeanb(nkp1), ematrixi(nk,nk), gmatrixr(nk,nk)
     &,    pweight(nkp2,3)
     &,    ematrixT(nk,nk), einversT(nk,nk), gmatrxtT(nk,nk)
     &,    tmatrixT(nk,nk), ematrxiT(nk,nk), gmatrxrT(nk,nk)
      integer msteps(nk), nnmis(4)
c
c dimension arrays for R-norm
      real Rematrix(nk,nk), Rinverse(nk,nk), Rmatrixt(nk,nk)
     &,    Rmatrixp(nk), Rtmatrix(nk,nk), Rpmatrix(nk), Redepth(nk)
     &,    Rtpmeanb(nkp1), Rmatrixi(nk,nk), Rmatrixr(nk,nk)
     &,    Rpweight(nkp2,3), Rmsteps(nk)
     &,    RinversT(nk,nk), RmatrxtT(nk,nk)
     &,    RmatrxiT(nk,nk), RmatrxrT(nk,nk)
c
c  dimension arrays for spectral q norm
      real qmatrix(nk,nk), qinverse(nk,nk)
      real qmatrixT(nk,nk), qinversT(nk,nk)
c
c
c dimension boundary fields
c
      real bdru(nij2,nk,nspong,2), bdrv(nij2,nk,nspong,2)
     &,    bdrt(nij2,nk,nspong,2), bdrq(nij2,nk,nspong,2)
     &,    bdrp(nij2,nspong,2), field3b(nijk,4), field2b(nij,2)
c
c dimension work arrays 
c (nfield3m is maximum size of nfield3 on any input file)
      parameter (nwork1=nij*nkp1, nfield3m=12
     &,          nworkr=2+(nij+2)*nfield3m )
      dimension work1(ni,nj,nkp1), work2(ni,nj,nkp1), work3(ni,nj,nkp1)
     &,         work4(ni,nj,nk), work5(ni,nj,nk), work6(ni,nj,nk)
     &,         worka(nijs,2), workb(nijs,6), workc(nijkx,12)
     &,         work7b(ni,nk,7), work8b(ni,nj,13), workr(nworkr)
     &,         work7(ni,nk,7), work8(ni,nj,13), work9(ni,nj,nkp2,2)
     &,         work10(ni,nkp1,14)
c
c  The following equivalence is optional to reduce memory required
c  If any of these fields are to be output, they should not be
c  equivalenced here.
c
      equivalence (work1,workr,workq,sdotd,ptv,work7b,work8b,work9
     &,            work10)
     &,           (work2,workc      ,sdotdb,tv,omega)
     &,           (work3,workb      ,pufb,pufx,ux,td,theta)
     &,           (work4,worka      ,pvfb,pvfx,vx,tdb,thetab)
     &,           (work5            ,pufxb,uxb,puf,cpmr,xflux)
     &,           (work6,work7,work8,pvfxb,vxb,pvf,xfluxb)
     &,           (pun,puten),(pvn,pvten),(ptn,ptten),(pqn,pqten)
c
      equivalence (field3b,u,um,qcoefsw), (field3b(1,2),v,vm)
     &,           (field3b(1,3),t,tm), (field3b(1,4),q,qm)
     &,           (field2b,pd,pdm)
c 
c dimension arrays for reading and writing file information
c
      parameter (nfieldh=4, nfldsrw3=8, nfldsrw2=6,
     &           nfldsr3=4, nfldsr3d=3, nfldsr2=2, nfldsr2d=2,
     &           nfldsw3=8, nfldsw3d=6, nfldsw2=6, nfldsw2d=4, 
     &           nfldsp3=4, nfldsp3d=3, nfldsp2=2, nfldsp2d=2 )
c
c  nfieldh  = number of 2-d fields on header
c  nfldsrw3 = number of 3-d fields to be read and/or written
c  nfldsrw2 = number of 2-d fields to be read and/or written
c  nfldsr3  = number of 3-d fields to be read from input b.s. file
c  nfldsr2  = number of 2-d fields to be read from input b.s. file
c  nfldsr3d = number of 3-d dry fields to be read from input bs file
c  nfldsr2d = number of 2-d dry fields to be read from input bs file
c  nfldsw3  = number of 3-d fields to be written to output file
c  nfldsw2  = number of 2-d fields to be written to output file
c  nfldsw3d  = number of 3-d dry fields to be written to output file
c  nfldsw2d  = number of 2-d dry fields to be written to output file
c  nfldsp3  = number of 3-d fields to be read from pert input file
c  nfldsp2  = number of 2-d fields to be read from pert input file
c  nfldsp3d = number of 3-d dry fields to be read from pert input file
c  nfldsp2d = number of 2-d dry fields to be read from pert input file
c  Also change field lists (cflds) below
c
c  nvars  = array of integer variables on header
c  rvars  = array of real variables on header
c  lvars  = array of logical variables on header
c
c  Suffixes indicate variable uses as follows:
c  `r`   for selection of fields to be read as basic state fields;
c        NOTE THAT THE ABOVE IS DIFFERENT THAN IN THE NLM 
c  `w`   for selection of fields to be written;
c  `i`   for information read from input BASIC STATE file;
c  `b`   for information read from input file used as boundary cond.;
c  `rw`  for information regarding the equivalenced arrays fields3,2;
c  `2`   for information regarding 2-d fields;
c  `3`   for information regarding 3-d fields;
c  `t`   for total (a sum of 2-d and 3-d fields).
c  `p`   for selection of initial perturbation fields from a file
c
      parameter (nfldswt=nfldsw3+nfldsw2+nfieldh,
     &           nfldsrt=nfldsr3+nfldsr2+nfieldh,
     &           nfldspt=nfldsp3+nfldsp2,
     &           nfldsrwt=nfldsrw3+nfldsrw2+nfieldh,
     &           ndimfld3=nijk*nfldsrw3, ndimfld2=nij*nfldsrw2)
      parameter (nfldmax=30, n1charsw=128, ndcharsw=6)
      parameter (nnftimeo=100, nnftimes=4+nnftimeo, nrvars=32+2*nk, 
     &           nlvars=30, nnvars=50+nk+nnftimes+2*nfldswt)
      data ibsheadw/1/, ibsformw/1/
c
      logical lvarsw(nlvars), lvarsi(nlvars), lvarsb(nlvars)
     &,       lvarsp(nlvars)
      real rvarsw(nrvars), rvarsi(nrvars), rvarsb(nrvars)
     &,    fieldh(nij,nfieldh), fields3(nijk,nfldsrw3)
     &,    fields2(nij,nfldsrw2), rvarsp(nrvars), tout(nnftimeo)
      integer nvarsw(nnvars), nvarsi(nnvars), nvarsb(nnvars)
     &,       ifldsr(5), ifldsw(5), ifldsrw(3), igtypew(nfldswt)
     &,       nvarsp(nnvars), ifldsp(5), nvarsiw(nnvars)
      character*(8) cfldsw(nfldswt), cfldsr(nfldsrt), cfldsrw(nfldsrwt)
     &,             cfldsi(nfldmax), cfldsb(nfldmax), cfldsp(nfldspt)
      character*(n1charsw) charsw(ndcharsw), charsi(ndcharsw)
      character*65 cusercom
      data ifldsr/nfldsr3,nfldsr2,nfieldh,nfldsr3d,nfldsr2d/
      data ifldsw/nfldsw3,nfldsw2,nfieldh,nfldsw3d,nfldsw2d/
      data ifldsp/nfldsp3,nfldsp2,0,nfldsp3d,nfldsp2d/
      data ifldsrw/nfldsrw3,nfldsrw2,nfieldh/
      equivalence (nvarsb,nvarsi), (lvarsb,lvarsi),  (rvarsb,rvarsi)
     &,           (cfldsi,cfldsb)
      equivalence (pu,fields3(1,1)), (pub,fields3(1,2),pumb)
     &,           (pv,fields3(1,3)), (pvb,fields3(1,4),pvmb)
     &,           (pt,fields3(1,5)), (ptb,fields3(1,6),ptmb)
     &,           (pq,fields3(1,7)), (pqb,fields3(1,8),pqmb)
      equivalence (ps,fields2(1,1)), (psb,fields2(1,2),psmb)
     &,           (tg,fields2(1,3)), (tgb,fields2(1,4),tgmb)
     &,           (raincv,fields2(1,5)), (rainncv,fields2(1,6))
      equivalence (topog,fieldh(1,1)), (sfctyp,fieldh(1,2))
     &,           (snowc,fieldh(1,3)), (tdeepg,fieldh(1,4))
      data cfldsrw/'      pu','     pub','      pv','     pvb',
     &             '      pt','     ptb','      pq','     pqb',
     &             '      ps','     psb','      tg','     tgb',
     &             '  raincv',' rainncv',
     &             '   topog','  sfctyp','   snowc','  tdeepg'/
      data  cfldsp/'      pu','      pv','      pt','      pq',
     &             '      ps','      tg'/
      data  cfldsr/'     pub','     pvb','     ptb','     pqb',
     &             '     psb','     tgb',
     &             '   topog','  sfctyp','   snowc','  tdeepg'/
      data  cfldsw/'      pu','     pub','      pv','     pvb',
     &             '      pt','     ptb','      pq','     pqb',
     &             '      ps','     psb','      tg','     tgb',
     &             '  raincv',' rainncv',
     &             '   topog','  sfctyp','   snowc','  tdeepg'/
      data igtypew/4*0,4*1,10*1/
      data tout/nnftimeo*0./
c
c  cfldsr  = list of names of fields to read from pert. input file
c  cfldsd  = list of names of fields to read from basic state data file
c  cfldsw  = list of names of fields to write to output file,
c  cfldsrw = list of names of which are the union of read/write lists
c  igtypew = a list of indicators for dot (0) or cross (1) grid
c            corresponding to the list in cfldsw
c  tout    = list of time (minutes) for irregular output (toutreg
c            less than 0.). Time 0. will not be output unless
c            specified. Values must be in ascending order.
c
c  The order of the read and write lists are: (1) 3-d dry fields;
c  (2) 3-d moist fields; (3) 2-d dry fields ; (4) 2-d moistfields;
c  (5) 2-d header fields.
c  The moist fields will be automatically removed from the lists
c  if lmoist=.f. is selected by the user.
c  The order of the list cfldsrw must be the order of the equivalencing
c  of: (1) fields3; (2) fields2; (3) fieldh.  For cfldsrw there is no
c  required order for dry vs. moist fields.
c
c
c                        DATA STATEMENTS
c
c *DATA 1
c  lcheck is for checking that two runs produce the same result
c  lprint=.f. removes most printing of run-time information
c  llast should be preset to .f. (see call to Annmi)
      logical lcheck,lprint
c user set
      data lcheck/.f./, lprint/.f./, llast/.f./
c end user set
c
c---------------------------------------------------------------------
c *DATA 2
c unit numbers and file names
      logical ldaccess,ldaccesi,lpshell,lcopy
      integer iuniti(2)
      real tskips(2)
      data ldaccess/.t./, ldaccesi/.t./, lpshell/.t./, lcopy/.t./
      data iuniti/10,11/, iunitb/15/, iunitw/20/, iunitm1/40/
     &,    iunitm2/60/, itapnumd/1/, keepdays/365/
     &,    iunitiw/80/, iunitir/81/,  iunita1/90/
     &,    iunitd/50/, iformhw/1/, iformtw/1/, iunitt/70/
     &,    tskips/tskips1,tskips2/
      parameter (nbsized=500+nfield3m*nij/2
     &,          nbsizec=500+nij+nqdimw/2
     &,          nbsizen=500+nijk+nq2dim*3/2
     &,          nbsizea=500+nij+(nq2dim+nq3dim)/2 )
c
c  ldaccess must = .t. (only direct access fiels may be written)
c  ldaccesi = .t. indicates all files are input as direct access
c  lpshell = .t. indicates UNIX Pshell commands are invoked using
c            Fortran calls to assign, acquire, and dispose files
c  lcopy = .t. if copy done in Pshell Dispose script
c  keepdays = retention period (days) indicated in Pshell Dispose
c  iformhw is format type for headers of data files written (=1 only)
c  iformtw is format type for t-slots of data files written (=1 only)
c
c  iuniti is unit for initial data file
c  iunitb is unit for lateral boundary adjoint file 
c  iunitw is unit for output of selected model fields
c  iunitiw is unit for output of innmi solution
c  iunitir is unit of iterates of basic state initialized fields
c  iunitm1 is unit for input of convectv precip coefs for TLM/ADJ
c  iunitm2 is unit for input of non conv precip coefs for TLM/ADJ
c  iunita1 is unit for input of dry conv. coefs for TLM/ADJ
c  itapnumd is a tape number counter for input data file
c  iunitd is unit for input of basic state data file
c  iunitt is unit for temporary save file
c  nbsized specifies i/o buffer size (words) for a data file
c  nbsizec specifies i/o buffer size (words) for a moist convec. file
c  nbsizen specifies i/o buffer size (words) for a non-conv prec. file
c
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      data charsw(1)/'ADJ  execution'/
      data charsw(2)/'DUMMY/NLMC3/'/
      data charsw(3)/'DUMMY/NLMP4/'/
      data charsw(4)/'no separate boundary file'/
c user set
      data charsw(5)/'/RAEDER/MAMS/J89DRY/070388/MST/H0-24/BS/'/
c end user set
      data charsw(6)/'DUMMY/ADJ3/'/
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
c
c  charsw(1) is type of program being executed
c  charsw(2) is name of initial data file
c  charsw(3) is path name for perturbed initial condition in TLM/ADJ
c  charsw(4) is name of lateral boundary file (unless lbdsame=.t.)
c  charsw(5) is path name for basic state files in TLM/ADJ
c  charsw(6) is path name for output files (must end in '/')
c  cusercom  is comment included in mass store field (no blanks)
c
c 
ccccccccccccccccccc  END OF DIMENSIONS  cccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                     BEGIN EXECUTION
c
c
      if (lmoist) then
         lrain=.t.
      else
         lrain=.f.
      endif
      cusercom(1:65)=charsw(1)(1:65)
c
ctlmadj
      FTIMED0=FTIMED0X
      FTIMEM0=FTIMEM0X
      deps=0.5*dtm                    ! value for rounding times
      ktimax=int((timax+deps)/dtm)    ! maximim number of times steps
c
c  Print some run parameters
      if (IFIRST.eq.0) then
         print *,' '
         print *,' P1: User-specified parameter list:'
         print *,'   ni=',ni,'  nj=',nj,'  nk=',nk,'  nspong=',nspong
     &          ,'  dt=',dt,'  timax=',timax
         print *,'   dtbdr=',dtbdr,'  toutreg=',toutreg,'  ipackw='
     &          ,ipackw,'  ftimed0=',ftimed0,'  ftimem0=',ftimem0
         print *,'   tupdated=',tupdated,'  tupdatem=',tupdatem
     &          ,'  lmoist=',lmoist,'  nqpoints=',nqpoints
         print *,'   lground=',lground,'  lcloud=',lcloud,'  lsnow='
     &          ,lsnow,'   lbdryp=',lbdryp,'  lbdsame=',lbdsame
         print *,'   ipert=',ipert,'  pert=',pert,'  tskips1=',tskips1
     &          ,'  tskips2=',tskips2
         print *,'   linitial=',linitial,'  ldrycon=',ldrycon
         print *,' '
      endif
c
c  Adjust and check lists of fields to read or write
      call Ncfldsm (cfldsr,ifldsr,cfldsrw,ifldsrw,lmoist,lprint,
     &              igtypew,1,'BS')
      call Ncfldsm (cfldsp,ifldsp,cfldsrw,ifldsrw,lmoist,lprint,
     &              igtypew,1,'IP')
      call Ncfldsm (cfldsw,ifldsw,cfldsrw,ifldsrw,lmoist,lprint,
     &              igtypew,nfldswt,'WR')
c
c
c                      INITIALIZE MODEL
c
c  Acquire and Open basic state file, Read and Get required
c  information from header.
c
ctlmadj BLOCKIF
c
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c	do acquire of dry files if ifirst=0, modification necessary
c	for moist case 
c
      IF (IFIRST.eq.0) then
c
c        acquire for dry runs
c
         if (lpshell) call Nacquire (charsw(5),'D',1,lprint,nbsized,
     &                            iunitd)
         if (ldaccesi) call Nopen (0,iunitd,lprint)
         if (linitial.and.linitrd) then  ! read b.s.innmi file
            if (lpshell) call Nacquire (charsw(5),'I',1,lprint,nbsized,
     &                            iunitir)
            if (ldaccesi) call Nopen (0,iunitir,lprint)
         endif
c
c        acquire for moist runs (basic state on 1 file only) 
c
         if (lmoist) then 
               if (lpshell) call Nacquire 
     &                      (charsw(5),'C',1,lprint,nbsizec,iunitm1)
               if (ldaccesi) call Nopen (0,iunitm1,lprint)
               if (lpshell) call Nacquire 
     &                      (charsw(5),'N',1,lprint,nbsizen,iunitm2)
               if (ldaccesi) call Nopen (0,iunitm2,lprint)
         endif
c
c        acquire for dry convective  (basic state on 1 file only) 
c
         if (ldrycon) then 
               if (lpshell) call Nacquire 
     &                      (charsw(5),'A',1,lprint,nbsizea,iunita1)
               if (ldaccesi) call Nopen (0,iunita1,lprint)
         endif
c
      RETURN
      ENDIF
c
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
      call Nreadh (nvarsi,lvarsi,rvarsi,cfldsi,charsi,fieldh,workr,
     &             cfldsrw,cfldsr,ifldsrw,ifldsr,ni,nj,nk,
     &             nij,ir1,iunitd,ldaccesi,.t.,lprint,
     &             nnvars,nlvars,nrvars,nworkr)
      ir1i=ir1  ! save for use in Annmi routine    
      call Ngethead (sigf,tpmeanb,ptop,dx,clat,clon,xhour0,time0,
     &               iyear0,imonth0,iday0,nsplit,radfreq,msteps,
     &               lsplit,lbdry,nnmis,tpratio,cfldsi,
     &               charsi,nk,nkp1,rvarsi,lvarsi,nvarsi,lprint,
     &               nnvars,nrvars,nlvars,' input')
      nbsdataf=int((tupdated+0.5*dtm)/dtm)
      if (abs(timax-ftimed0).gt.dtm*0.5) nnmis(4)=0
c
c  Set various constants
      call Nsetup (Cp,R,xkap,grav,ep1,Rv,ecpm,pqmin,hlatent,
     &             rhdmax,asselin0,asselin1,balpha,bbeta,bfact,
     &             hdeepfac,solcon,stbolt,klevrad,klevcld,degrad,
     &             dx,ptop,nspong,pkappa,qsatvp,npkappa,nqsatvp,
     &             galbedo,gthrmin,gemissv,groughl,gmavail,
     &             gheatcpr,cdragx,cdragd,xland,iland,sfctyp,
     &             snowc,iyear0,imonth0,iday0,xhour0,time0,0.,
     &             xlat,xlon,dlat,dlon,fmapx,fmapd,fmapx2,fmapd2,
     &             cf,xkc,sigf,sigh,dsig,clat,clon,ematrix,
     &             einverse,gmatrixt,gmatrixp,tmatrix,pmatrix,
     &             edepth,tpmeanb,msteps,ptten,psten,puten,
     &             pvten,sdot,work3,work2,work5,dt,nix,njx,nijx,
     &             nijkx,nk,ni,nj,nkp1,nsplit,ematrixi,gmatrixr,
     &             pweight,lprint,linitial,lsplit)
c
c  these are new matrices for r-norm (isothermal T) 
         call Nrsetc (Rtpmeanb,nk,270.)
         Rtpmeanb(nkp1)=100.0 
         call Sinitial (Rematrix,Rinverse,Rmatrixt,Rmatrixp,Rtmatrix,
     &               Rpmatrix,Redepth,Rtpmeanb,Rmsteps,
     &               ptten,psten,puten,pvten,sdot,work3,work2,
     &               Rmatrixi,Rmatrixr,Rpweight,work5,sigf,sigh,dsig,
     &               R,Cp,ptop,dx,dt,nix,njx,nijkx,nk,nkp1,
     &               1,.false.,.true.,.false.)  !  no write 
         call Atransp (RinversT,Rinverse,nk,nk)
         call Atransp (RmatrxtT,Rmatrixt,nk,nk)
         call Atransp (RmatrxiT,Rmatrixi,nk,nk)
         call Atransp (RmatrxrT,Rmatrixr,nk,nk)
c
c  compute matrices for spectral q norm
         call QNmatrix (qmatrix,qinverse,qtotvar,tpmeanb,sigh,dsig,
     &               ptop,qsatvp,work1,work1(1,1,2),work2(1,1,2),work2,
     &               nqsatvp,nk,lprint)
         call Atransp (qinversT,qinverse,nk,nk)
         call Atransp (qmatrixT,qmatrix,nk,nk)
c
c  initialize constants and arrays for lateral boundary adjoint
      if (lbdryp) then
         ktaubd=0
         nbdsteps=int((dtbdr*60.+0.5*dt)/dt)
         call Nrsetc (bdru,nij2n2k,0.)
         call Nrsetc (bdrv,nij2n2k,0.)
         call Nrsetc (bdrt,nij2n2k,0.)
         if (lmoist) call Nrsetc (bdrq,nij2n2k,0.)
         call Nrsetc (bdrp,nij2n2,0.)
         idum=1
         if (lpshell) call Nassign (lprint,iunitb,idum,'B')
         open (iunitb,form='unformatted',iostat=iostatb)
         if (iostatb.eq.0) then       
            if (lprint) then
               print *,' '
               print *,' Open sequential file on unit ',iunitb,
     &                 ' for output of boundary adjoint'
            endif
         else
            print *,' '
            print *,'CAN NOT OPEN BOUNDARY ADJ FILE ON UNIT ',iunitb
            print *,' ERROR=',iostatb,' RETURNED FROM OPEN STATEMENT'
         endif
      endif
c
c  transpose matrices for split-explicit and innmi calculations
         call Atransp (einversT,einverse,nk,nk)
         call Atransp (ematrixT,ematrix,nk,nk)
         call Atransp (tmatrixT,tmatrix,nk,nk)
         call Atransp (gmatrxtT,gmatrixt,nk,nk)
         call Atransp (ematrxiT,ematrixi,nk,nk)
         call Atransp (gmatrxrT,gmatrixr,nk,nk)
c 
c  read basic state data and required coefficients for initial time
c
      ktaubs=0
      nfiled=1
      call Bsdatad (pub,pvb,ptb,psb,ub,vb,tb,pdb,pufb,pvfb,
     &              pkhb,sdotb,omegab,phib,chi,divb,uxb,vxb,
     &              pqb,qb,tvb,pkappa,qsatvp,
     &              alogpb,work9,ptvb,cpmrb,sigf,sigh,dsig,
     &              fmapx,fmapd,fmapx2,topog,workr,fields3,fields2,
     &              rvarsi,nrvars,grav,R,Cp,ptop,dx,ep1,ecpm,
     &              tupdated,ftimed0,nvarsi,ifldsrw,ifldsr,0,
     &              -1,ktaubs,npkappa,nqsatvp,ni,nj,nk,nkp1,
     &              iunitd,0.,nfiled,nbsized,charsw(5),
     &              cfldsrw,cfldsr,cfldsi,lmoist,ldaccesi,lprint,
     &              lpshell,qsb,nnvars,nworkr,ndimfld3,ndimfld2,
     &              pqmin,'ADJ')
      if (lmoist) then
         call Bsdataq0 (nvarsq1,rvarsq1,'C',iunitm1,ktaum1,charsw(5),
     &                  nnvarsq,nrvarsq,nq2dim,nq3dim,2,4,
     &                  qcoefs1,ipnoncv,nbsizec,lprint,.f.,.f.,
     &                  nnvarsq,nrvarsq,nij,nq4dim,nfilem1)
         call Bsdataq0 (nvarsq2,rvarsq2,'N',iunitm2,ktaum2,charsw(5),
     &                  nnvarsq,nrvarsq,nq2dim,nq2dim,2,1,
     &                  qcoefs3,ipnoncv,nbsizen,lprint,.f.,.f.,
     &                  nnvarsq,nrvarsq,nijk,nq5dim,nfilem2)
         call Bsdataq1 (qcoefs1,ipconvt,rvarsq1,nvarsq1,iunitm1,'C',
     &                  ktaum1,nfilem1,tupdatem,ftimem0,workq,
     &                  qcoefsw,work6,nqmaxxx1,
     &                  charsw(5),0.,lprint,lpshell,ldaccesi,'TLM',
     &                  nnvarsq,nrvarsq,nij,nq4dim,nqdimw)
         call Bsdataq1 (qcoefs3,ipnoncv,rvarsq2,nvarsq2,iunitm2,'N',
     &                  ktaum2,nfilem2,tupdatem,ftimem0,workq,
     &                  qcoefsw,work6,nqmaxxx2,
     &                  charsw(5),0.,lprint,lpshell,ldaccesi,'TLM',
     &                  nnvarsq,nrvarsq,nijk,nq5dim,nqdimw)
      endif 
      if (ldrycon) then
         call Bsdataq0 (nvarsa1,rvarsa1,'A',iunita1,ktaua1,charsw(5),
     &                  nnvarsq,nrvarsq,nq2dim,nq3dim,1,1,
     &                  dcoefs1,ipdconv,nbsizea,lprint,.f.,.f.,
     &                  nnvarsq,nrvarsq,nij,nq6dim,nfilea1)
         call Bsdataq1 (dcoefs1,ipdconv,rvarsa1,nvarsa1,iunita1,'A',
     &                  ktaua1,nfilea1,tupdatem,ftimem0,workq,
     &                  qcoefsw,work6,ndmaxxx1,
     &                  charsw(5),0.,lprint,lpshell,ldaccesi,'TLM',
     &                  nnvarsq,nrvarsq,nij,nq6dim,nqdimw)
      endif ! test on ldrycon
         nbsdataq=int((tupdatem+0.5*dtm)/dtm)
c
ctlmadj
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c                  INITIALIZE TLM
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
	IF ( ISCALE .EQ. 0 ) THEN 
c
        call Xnorms4 (u,v,t,ps,q,tg,xstar,xfull,lmoist,
     &              sigf,sigh,dsig,grav,R,Cp,clat,clon,dx,
     &              Redepth,Rmatrixi,Rmatrixp,Rmatrixr,Rpweight,
     &              qmatrix,
     &              Rtpmeanb,work1,work2,work3,work4,work5,work6,
     &              cf,fmapx,fmapd,ni,nj,nk,nxstar,nfull,itype) 
c	
	RETURN
	ENDIF
	
c
c	provide here a call to do the dimensioning option (iscale=0)
c
c 	set u,v,t,q,ps,tg from xstar
c
      call Xnorms1 (u,v,t,ps,q,tg,xstar,lmoist,
     &              sigf,sigh,dsig,grav,R,Cp,clat,clon,dx,
     &              Redepth,Rmatrixi,Rmatrixp,Rmatrixr,Rpweight,
     &              qmatrix,
     &              Rtpmeanb,work1,work2,work3,work4,work5,work6,
     &              cf,fmapx,fmapd,ni,nj,nk,nxstar,itype) 
c
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      call PDECUPL (pu,pv,pt,pq,ps,u,v,t,q,
     &              pub,pvb,ptb,pqb,psb,ub,vb,tb,qb,
     &              pd,pdb,lmoist,ni,nj,nk,2)
c
c  balance initial condition using implicit nnmi
      if (linitial) then
         if (ftimed0.ne.0.) nnmis(4)=0
         call Fnnmi (pu,pv,pt,pq,ps,u,v,t,q,pd,puten,pvten,ptten,chi,
     &               psten,ux,vx,td,tv,ptv,alogp,work9,sdot,sdotd,
     &               puf,pvf,pufx,pvfx,cpmr,omega,phi,fmapx,fmapd,
     &               fmapx2,fmapd2,topog,cf,sigf,sigh,dsig,grav,R,Cp,
     &               ep1,ptop,div,pum,phi,pvm,ptm,work1,work2,
     &               work3,chi,einverse,ematrixi,gmatrixp,gmatrixt,
     &               gmatrixr,edepth,pweight,work1,tpmeanb(nkp1),clat,
     &               ecpm,dx,tpratio,nnmis,ni,nj,nk,nkp1,lmoist,lprint,
     &               nvarsi,cfldsi,charsw,workr,iunitir,cfldsrw,ifldsrw,
     &               lpshell,ldaccess,nnvars,nworkr,fields3,fields2,
     &               ndimfld3,ndimfld2,puten,pvten,pum,pvm,ptm,
     &               psm,pub,pvb,ptb,pqb,psb,pdb,sdotb,
     &               sdotdb,omegab,pufb,pufxb,pvfb,pvfxb,ub,vb,uxb,
     &               vxb,tb,cpmrb,qb,tvb,phib,ptvb,tdb,alogpb,divb,
     &               nbsized,pqmin,ir1i,work8,cfldsr,ifldsr,.true.)
c
c
c  compute some derived basic state fields and coefficients
         if (nnmis(4).gt.0) then  
            call Btendc (pub,pvb,ptb,psb,ub,vb,tb,pdb,pufb,pvfb,
     &                   pkhb,sdotb,omegab,phib,chi,divb,uxb,vxb,tvb,
     &                   ptvb,pqb,qb,cpmrb,alogpb,work9,qsatvp,
     &                   sigf,sigh,dsig,fmapx,fmapd,fmapx2,
     &                   topog,pkappa,grav,R,Cp,ep1,ecpm,ptop,dx,
     &                   qsb,ni,nj,nk,nkp1,npkappa,nqsatvp,lmoist)
         endif  ! end block on setting basic state derived fields
      endif     ! end block on initialization
c
c  Set fields at t=-dt equal to those at time 0
      call Ncopyv (pum,pu,nijk)
      call Ncopyv (pvm,pv,nijk)
      call Ncopyv (ptm,pt,nijk)
      call Ncopyv (psm,ps, nij)
      if (lmoist)  call Ncopyv (pqm,pq,nijk)
      call Ncopyv (tgm,tg, nij)
c
ctlmadj
c  REMOVE WRITE OF INITIAL DATA
c
c
c
c                   BEGIN TLM TIME INTEGRATION
c
c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
c
      ttime0=time0+ftimed0/60.
ctlmadj
      IF (KTIMAX.GT.1) THEN      ! DO LOOP OVER TIME
      do 1 nstep=1,ktimax
c
c  Set times (hours from start)
      xtime=dth*float(nstep)    ! time predicted at end of step
      xtime1=dth*float(nstep-1) ! current time (for pu, pv etc.)
c
      if (lcheck) then
         print *,'nstep=',nstep,'  xtime1=',xtime1
         i=ni/2
         j=nj/2
         print 444,nstep,pu(i,j,10),pq(i,j,10),pt(i,j,10),ps(i,j)
         print 444,nstep,pub(i,j,10),pqb(i,j,10),ptb(i,j,10),psb(i,j)
  444    format(i3,1p4e21.13)
      endif
c
      if (nstep.eq.1) then
         dti=dt
      else
         dti=dt2
      endif
c
c  reset tendencies to zero
      call Nrsetc (puten,nijk,0.)
      call Nrsetc (pvten,nijk,0.)
      call Nrsetc (ptten,nijk,0.)
      call Nrsetc (psten, nij,0.)
      if (lmoist) call Nrsetc (pqten,nijk,0.)
      if (lground) call Nrsetc (tgten,nij,0.)
c
c  Copy edge values from boundary arrays and compute sponge forcing
c  or set edges to zero 
      if (lbdryp) then
         call NLbforce (puten,pvten,ptten,pqten,psten,
     &                  bdru,bdrv,bdrt,bdrq,bdrp,
     &                  bfact,balpha,bbeta,nstep,nbdsteps,
     &                  lmoist,ni,nj,nk,nspong,nij2)
      else
         call Cbdryx ( pu, pv, pt, ps,ni,nj,nk)
         call Cbdryx (pum,pvm,ptm,psm,ni,nj,nk)
         if (lmoist) then 
            call Cbdry0 ( pq,ni,nj,nk,1)
            call Cbdry0 (pqm,ni,nj,nk,1)
         endif
      endif
c
c  Compute diffusion, boundary damping, and PBL forcing
      call Ftendm (puten,pvten,ptten,pqten,psten,tgten,pum,pvm,
     &             ptm,pqm,um,vm,tm,qm,psm,tgm,pdm,umb,vmb,tmb,
     &             qmb,psmb,tgmb,pkhb,pdmb,xflux,xfluxb,xkc,
     &             gheatcpr,work1,work2,sigh,sigf,dsig,bfact,
     &             qsatvp,pqmb,theta,thetab,work1,cdragx,cdragd,
     &             pumb,pvmb,ptmb,xpprp,gmavail,xkap,Cp,
     &             dx,dti,ptop,R,grav,hlatent,balpha,bbeta,lground,
     &             lmoist,nqsatvp,ni,nj,nk,nkp1,nkpx,nspong,lbdry)
c
c  Compute some contributions to the ground temperature tendency
      if (lground) 
     &   call Ftgten (tgten,tgm,pqm,ptm,glw,gsw,cldfr,psm,tm,qm,prw,
     &                qs,work1,tgmb,pqmb,cldfrb,psmb,tmb,qmb,prwb,qsb,
     &                work2,gemissv,gheatcpr,iland,work3,stbolt,
     &                hdeepfac,dsig,sigh,ptop,grav,lmoist,lcloud,
     &                qsatvp,nqsatvp,ni,nj,nk,nstep,radfreq,klevcld,
     &                klevrad,xlat,xlon,galbedo,degrad,solcon,lprint,
     &                xtime1,iyear0,imonth0,iday0,xhour0,ttime0,dtm)
c
c  Compute adiabatic tendencies
      call Ftendc (psten,puten,pvten,ptten,sdot,pd,pu,pv,pt,ps,u,
     &             v,t,puf,pufx,pvf,pvfx,ux,vx,div,sdotd,omega,
     &             alogp,pub,pvb,ptb,psb,alogpb,sdotb,divb,
     &             omegab,pufb,pufxb,pvfb,pvfxb,ub,vb,uxb,vxb,tb,
     &             sdotdb,pdb,chi,fmapx,fmapd,fmapx2,fmapd2,cf,
     &             sigf,sigh,dsig,ptv,ptvb,cpmr,cpmrb,q,pq,
     &             qb,lmoist,ep1,ecpm,R,Cp,ptop,dx,grav,ni,nj,nk,
     &             pqb,pqten,nkp1,td,tdb,tv,tvb,phi,phib,work8,
     &             work9,.true.)
c
c  Compute new values determined from explicit scheme
      call Nstepf1 (pun,pum,puten,dti,ni,nj,nk,0)
      call Nstepf1 (pvn,pvm,pvten,dti,ni,nj,nk,0)
      call Nstepf1 (ptn,ptm,ptten,dti,ni,nj,nk,1)
      call Nstepf1 (psn,psm,psten,dti,ni,nj, 1,1)
      if (lmoist) call Nstepf1 (pqn,pqm,pqten,dti,ni,nj,nk,1)
      if (lground) call Nsteptg (tgn,tgm,tgten,snowc,tgb,dti,ni,nj,
     &                           lsnow,iland,xland)
c
c  Update edges on boundary
      if (lbdryp) then
         if ((nstep.gt.1).and.(mod(nstep,nbdsteps).eq.1))
     &      call NLbread (bdru,bdrv,bdrt,bdrq,bdrp,field3b,field2b,
     &                    workr,xtime,xtimef,lmoist,cfldsb,rvarsb,
     &                    nvarsb,nnvars,nrvars,dtbdr,ktaubd,nij2,ni,nj,
     &                    nk,nijk,nworkr,nspong,ftimeb0,ldaccesi,
     &                    lprint,iunitb)
         call NLbsete (pun,pvn,ptn,pqn,psn,bdru,bdrv,bdrt,bdrq,bdrp,
     &                 nstep,nbdsteps,lmoist,ni,nj,nk,nspong,nij2,-1)
      else
         call Cbdryx (pun,pvn,ptn,psn,ni,nj,nk)
         if (lmoist) call Cbdry0 (pqn,ni,nj,nk,1)
      endif
c
c  Compute split corrections
      if (lsplit.and.(nstep.ne.1))
     &   call Fseqnce (pun,pvn,ptn,psn,pu,pv,pt,ps,pum,pvm,ptm,psm,
     &                 einverse,ematrix,tmatrix,pmatrix,gmatrixt,
     &                 gmatrixp,fmapx2,fmapd,edepth,msteps,
     &                 work1,worka,worka(1,2),
     &                 workb,workb(1,4),worka,worka(1,2),
     &                 tb,psb,dx,dt,ni,nj,nk,nsplit,nij)
c
c  Moist convection
      if (lmoist) then
         nqcoefs2=nqmaxxx1*2*nk+1
         call Fconvp (pqn,ptn,psn,raincv,work1,grav,dsig,qcoefs1,
     &                qcoefs1(nqcoefs2),ipconvt,nqmaxxx1,ni,nj,nk,
     &                dt,dti,raincv2)
      endif
c
c  Adjust for dry convection
      if (ldrycon) then
         ndcoefs2=ndmaxxx1*nk+1
         call Fdryconv (ptn,psn,work1,dcoefs1,dcoefs1(ndcoefs2),
     &                 ipdconv,ndmaxxx1,ni,nj,nk,nspong)
      endif ! test on ldrycon
c
c  Adjust for nonconvective precipitation
      if (lmoist) call Fnconvp (pqn,ptn,psn,rainncv,pqb,dsig,grav,
     &              pqmin,ipnoncv,qcoefs3,hlatent,Cp,nqmaxxx2,
     &              ni,nj,nk,dt,dti,rainncv2)
c
c  Apply Asselin time filter
      call Nwavg (pum,pum,pu,pun,asselin0,asselin1,ni,nj,nk,0,1)
      call Nwavg (pvm,pvm,pv,pvn,asselin0,asselin1,ni,nj,nk,0,1)
      call Nwavg (ptm,ptm,pt,ptn,asselin0,asselin1,ni,nj,nk,1,1)
      call Nwavg (psm,psm,ps,psn,asselin0,asselin1,ni,nj, 1,1,1)
      if (lmoist)
     &   call Nwavg (pqm,pqm,pq,pqn,asselin0,asselin1,ni,nj,nk,1,1)
      if (lground)
     &   call Nwavg (tgm,tgm,tg,tgn,asselin0,asselin1,ni,nj, 1,1,1)
c
c  Shift time values
      if (lbdryp) then
         call NLbcope (pum,pu,ni,nj,nk,0)
         call NLbcope (pvm,pv,ni,nj,nk,0)
         call NLbcope (ptm,pt,ni,nj,nk,1)
         call NLbcope (psm,ps,ni,nj, 1,1)
         if (lmoist) call NLbcope (pqm,pq,ni,nj,nk,1)
      endif
c
      call Ncopyv (pu,pun,nijk)
      call Ncopyv (pv,pvn,nijk)
      call Ncopyv (pt,ptn,nijk)
      call Ncopyv (ps,psn, nij)
      if (lmoist)  call Ncopyv (pq,pqn,nijk)
      if (lground) call Ncopyv (tg,tgn, nij)
c
c read basic state data if update time
      if (mod(nstep,nbsdataf).eq.0)
     &   call Bsdatad (pub,pvb,ptb,psb,ub,vb,tb,pdb,pufb,pvfb,
     &                 pkhb,sdotb,omegab,phib,chi,divb,uxb,vxb,
     &                 pqb,qb,tvb,pkappa,qsatvp,
     &                 alogpb,work9,ptvb,cpmrb,sigf,sigh,dsig,
     &                 fmapx,fmapd,fmapx2,topog,workr,fields3,fields2,
     &                 rvarsi,nrvars,grav,R,Cp,ptop,dx,ep1,ecpm,
     &                 tupdated,ftimed0,nvarsi,ifldsrw,ifldsr,nstep,
     &                 ktimax,ktaubs,npkappa,nqsatvp,ni,nj,nk,nkp1,
     &                 iunitd,xtime,nfiled,nbsized,charsw(5),
     &                 cfldsrw,cfldsr,cfldsi,lmoist,ldaccesi,lprint,
     &                 lpshell,qsb,nnvars,nworkr,ndimfld3,ndimfld2,
     &                 pqmin,'TLM')
      if (lmoist.and.(mod(nstep,nbsdataq).eq.0).and.
     &               (nstep.ne.ktimax)) then
         call Bsdataq1 (qcoefs1,ipconvt,rvarsq1,nvarsq1,iunitm1,'C',
     &                  ktaum1,nfilem1,tupdatem,ftimem0,workq,
     &                  qcoefsw,work6,nqmaxxx1,
     &                  charsw(5),xtime,lprint,lpshell,ldaccesi,'TLM',
     &                  nnvarsq,nrvarsq,nij,nq4dim,nqdimw)
         call Bsdataq1 (qcoefs3,ipnoncv,rvarsq2,nvarsq2,iunitm2,'N',
     &                  ktaum2,nfilem2,tupdatem,ftimem0,workq,
     &                  qcoefsw,work6,nqmaxxx2,
     &                  charsw(5),xtime,lprint,lpshell,ldaccesi,'TLM',
     &                  nnvarsq,nrvarsq,nijk,nq5dim,nqdimw)
      endif 
      if (ldrycon.and.(mod(nstep,nbsdataq).eq.0).and.
     &               (nstep.ne.ktimax))
     &   call Bsdataq1 (dcoefs1,ipdconv,rvarsa1,nvarsa1,iunita1,'A',
     &                  ktaua1,nfilea1,tupdatem,ftimem0,workq,
     &                  qcoefsw,work6,ndmaxxx1,
     &                  charsw(5),xtime,lprint,lpshell,ldaccesi,'TLM',
     &                  nnvarsq,nrvarsq,nij,nq6dim,nqdimw)
ctlmadj
c  REMOVE WRITE OF DATA
c
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
    1 continue ! end time loop
ctlmadj
      ENDIF ! TEST ON KTIMAX.GT.1
c
      FTIMED0=FTIMED0X+timax
      FTIMEM0=FTIMEM0X+timax
      ktaubs=0
      nfiled=1
      call Bsdatad (pub,pvb,ptb,psb,ub,vb,tb,pdb,pufb,pvfb,
     &              pkhb,sdotb,omegab,phib,chi,divb,uxb,vxb,
     &              pqb,qb,tvb,pkappa,qsatvp,
     &              alogpb,work9,ptvb,cpmrb,sigf,sigh,dsig,
     &              fmapx,fmapd,fmapx2,topog,workr,fields3,fields2,
     &              rvarsi,nrvars,grav,R,Cp,ptop,dx,ep1,ecpm,
     &              tupdated,ftimed0,nvarsi,ifldsrw,ifldsr,0,
     &              -1,ktaubs,npkappa,nqsatvp,ni,nj,nk,nkp1,
     &              iunitd,0.,nfiled,nbsized,charsw(5),
     &              cfldsrw,cfldsr,cfldsi,lmoist,ldaccesi,lprint,
     &              lpshell,qsb,nnvars,nworkr,ndimfld3,ndimfld2,
     &              pqmin,'ADJ')
c
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c             END OF TLM INTEGRATION
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
C       THIS BLOCK APPLIES THE OPERATOR A TRANS A TO THE RESULT
C       OF THE TLM INTEGRATION. FOR ITEN=1, THIS IS DONE IN MSTAR2.
C       ALSO WITHIN MSTAR2 PART OF THE AQ-CHECK IS DONE. FOR IPEN=1,
C       THE ROUTINES PVfs2pv AND PVfs2pvA ARE USED. (12/12/93).
C
      call PDECUPL (pu,pv,pt,pq,ps,u,v,t,q,
     &              pub,pvb,ptb,pqb,psb,ub,vb,tb,qb,
     &              pd,pdb,lmoist,ni,nj,nk,1)
c
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
        call Nrsetc (raincv,nij,0.)
        call Nrsetc (raincv2,nij,0.)
        call Nrsetc (rainncv,nij,0.)
        call Nrsetc (rainncv2,nij,0.)

	if ( itype.le.4 .or. itype.eq.6 .or. itype.eq.7 
     &                  .or. itype.eq.8 .or. itype.eq.9 ) then 
c
      call Xnorms2a (u,v,t,ps,q,tg,xfull,lmoist,xnorm1,
     &               sigf,sigh,dsig,grav,R,Cp,clat,clon,dx,
     &               Redepth,Rinverse,Rmatrixp,Rmatrixt,Rpweight,
     &               qinverse,
     &               Rtpmeanb,work1,work2,work3,work4,work5,work6,
     &               cf,fmapx,fmapd,ni,nj,nk,nfull,nxstar,itype)
c
c  Initialize adjoint
c
      call Xnorms2b (u,v,t,ps,q,tg,xfull,lmoist,
     &               sigf,sigh,dsig,grav,R,Cp,clat,clon,dx,
     &               Redepth,RinversT,Rmatrixp,RmatrxtT,Rpweight,
     &               qinversT,
     &               Rtpmeanb,work1,work2,work3,work4,work5,work6,
     &               cf,fmapx,fmapd,ni,nj,nk,nfull,nxstar,itype)
c
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c	
	endif
c
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      pu=0.
      pv=0.
      pt=0.
      if (lmoist) pq=0.
      call PDECUPL (pu,pv,pt,pq,ps,u,v,t,q,
     &              pub,pvb,ptb,pqb,psb,ub,vb,tb,qb,
     &              pd,pdb,lmoist,ni,nj,nk,3)
c
      pum=0.
      pvm=0.
      ptm=0.
      psm=0.
      tgm=0.
      if (lmoist) pqm=0.
c
c
c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
c
c
      if (lmoist) then
         ktaum1=1
         ktaum2=1
      endif
      if (ldrycon) then
         ktaua1=1
      endif
c
      ttime0=time0+ftimed0/60.
      iawrite=2
      dti = dt2
      lbdryx = .f.
c
        if ( ktimax.gt.0 ) then
        do 2 nstep = 1 , ktimax
c
        if (nstep.eq.ktimax) then
           iawrite=1
           dti = dt
        endif
c
c  set factor for rainfall adjoint
      if (lmoist) then
         if (nstep.eq.1) then
            if (ktimax.eq.1) then
               precipfc=1.
            else
               precipfc=.75
            endif
         else
            precipfc=.5
         endif
      endif
c
c Set times
      xtime=-dth*float(nstep)    ! time in hours
      xtime1=-dth*float(nstep-1) ! time in hours for previous time step
c
      if (lcheck) then
         print *,'nstep=',nstep,'  xtime1=',xtime1
         i=ni/2
         j=nj/2
         print 445,nstep,pu(i,j,10),pq(i,j,10),pt(i,j,10),ps(i,j)
         print 445,nstep,pub(i,j,10),pqb(i,j,10),ptb(i,j,10),psb(i,j)
  445    format(i3,1p4e21.13)
      endif
c
ctlmadj  remove write
c
c read basic state data if update time
      if ((iawrite.ne.1).and.mod(nstep-1,nbsdataf).eq.0) then
         if (linitial.and.(nnmis(4).gt.0).and.
     &        (ktimax-nstep.le.nbsdataf)) then  ! read b.s. innmi file
            iunitx=iunitir
            ftimed0y=0.   ! ctlmadj change from ftimed0y here and below 
            tupdatex=1.
            ktaubsx=nnmis(4)
         else                              ! read from b.s. data file
            iunitx=iunitd
            ftimed0y=ftimed0
            tupdatex=tupdated
            ktaubsx=ktaubs
         endif
         call Bsdatad (pub,pvb,ptb,psb,ub,vb,tb,pdb,pufb,pvfb,
     &              pkhb,sdotb,omegab,phib,chi,divb,uxb,vxb,
     &              pqb,qb,tvb,pkappa,qsatvp,
     &              alogpb,work9,ptvb,cpmrb,sigf,sigh,dsig,
     &              fmapx,fmapd,fmapx2,topog,workr,fields3,fields2,
     &              rvarsi,nrvars,grav,R,Cp,ptop,dx,ep1,ecpm,
     &              tupdatex,ftimed0y,nvarsi,ifldsrw,ifldsr,nstep,
     &              ktimax,ktaubsx,npkappa,nqsatvp,ni,nj,nk,nkp1,
     &              iunitx,xtime1,nfiled,nbsized,charsw(5),
     &              cfldsrw,cfldsr,cfldsi,lmoist,ldaccesi,lprint,
     &              lpshell,qsb,nnvars,nworkr,ndimfld3,ndimfld2,
     &              pqmin,'ADJ')
         if (iunitx.eq.iunitd) ktaubs=ktaubsx
      endif
      if ((iawrite.ne.1).and.(mod(nstep-1,nbsdataq).eq.0)
     &               .and.(nstep.ne.ktimax)) then
         if (lmoist) then
            call Bsdataq1 (qcoefs1,ipconvt,rvarsq1,nvarsq1,iunitm1,'C',
     &                  ktaum1,nfilem1,tupdatem,ftimem0,workq,
     &                  qcoefsw,work6,nqmaxxx1,
     &                  charsw(5),xtime1,lprint,lpshell,ldaccesi,'ADJ',
     &                  nnvarsq,nrvarsq,nij,nq4dim,nqdimw)
            call Bsdataq1 (qcoefs3,ipnoncv,rvarsq2,nvarsq2,iunitm2,'N',
     &                  ktaum2,nfilem2,tupdatem,ftimem0,workq,
     &                  qcoefsw,work6,nqmaxxx2,
     &                  charsw(5),xtime1,lprint,lpshell,ldaccesi,'ADJ',
     &                  nnvarsq,nrvarsq,nijk,nq5dim,nqdimw)
         endif
         if (ldrycon)
     &      call Bsdataq1 (dcoefs1,ipdconv,rvarsa1,nvarsa1,iunita1,'A',
     &                  ktaua1,nfilea1,tupdatem,ftimem0,workq,
     &                  qcoefsw,work6,ndmaxxx1,
     &                  charsw(5),xtime1,lprint,lpshell,ldaccesi,'ADJ',
     &                  nnvarsq,nrvarsq,nij,nq6dim,nqdimw)
      endif ! test on iawrite and nstep
c
c  Copy values
      call Ncopyv (pun,pu,nijk)
      call Ncopyv (pvn,pv,nijk)
      call Ncopyv (ptn,pt,nijk)
      call Ncopyv (psn,ps, nij)
      if (lmoist)  call Ncopyv (pqn,pq,nijk)
      if (lground) call Ncopyv (tgn,tg, nij)
c
c  Adjoint of Asselin filter
      call Awavg (pum,pum,pu,pun,asselin0,asselin1,ni,nj,nk,0,1)
      call Awavg (pvm,pvm,pv,pvn,asselin0,asselin1,ni,nj,nk,0,1)
      call Awavg (ptm,ptm,pt,ptn,asselin0,asselin1,ni,nj,nk,1,1)
      call Awavg (psm,psm,ps,psn,asselin0,asselin1,ni,nj, 1,1,1)
      if (lmoist)
     &   call Awavg (pqm,pqm,pq,pqn,asselin0,asselin1,ni,nj,nk,1,1)
      if (lground)
     &   call Awavg (tgm,tgm,tg,tgn,asselin0,asselin1,ni,nj, 1,1,1)
c
c  Adjoint of nonconvective precipitation
      if (lmoist) call Anconvp (pqn,ptn,psn,pqb,pqmin,ipnoncv,qcoefs3,
     &                          dsig,grav,hlatent,Cp,nqmaxxx2,ni,nj,nk,
     &                          rainncv,precipfc)
c
c  Adjoint of dry convection
      if (ldrycon) then
         ndcoefs2=ndmaxxx1*nk+1
         call Adryconv (ptn,psn,work1,dcoefs1,dcoefs1(ndcoefs2),
     &                  ipdconv,ndmaxxx1,ni,nj,nk,nspong)
      endif ! test on ldrycon
c
c  Adjoint for moist convection
      if (lmoist) then
         nqcoefs2=nqmaxxx1*2*nk+1
         call Aconvp (pqn,ptn,psn,work1,qcoefs1,qcoefs1(nqcoefs2),
     &                ipconvt,raincv,precipfc,dsig,nqmaxxx1,grav,
     &                ni,nj,nk)
      endif
c
c  Adjoint of time-split corrections
      if (lsplit.and.(iawrite.ne.1)) 
     &   call Aseqnce (pun,pvn,ptn,psn,pu,pv,pt,ps,pum,pvm,ptm,psm,
     &                 einversT,ematrixT,tmatrixT,pmatrix,gmatrxtT,
     &                 gmatrixp,fmapx2,fmapd,edepth,msteps,
     &                 work1,worka,worka(1,2),
     &                 workb,workb(1,4),worka,worka(1,2),
     &                 tb,psb,dx,dt,ni,nj,nk,nsplit,nij)
c
c  increment boundary adjoint and write out adjoint boundary file
      if (lbdryx) then 
         if (lsplit.and.(iawrite.ne.1)) 
     &      call ALbsete (pun,pvn,ptn,pqn,psn,bdru,bdrv,bdrt,bdrq,bdrp,
     &                    nstep,nbdsteps,lmoist,ni,nj,nk,nspong,nij2,
     &                    ktimax,1,-1)
         if (mod(ktimax-nstep,nbdsteps).eq.0) 
     &      call ALbwrite (bdru,bdrv,bdrt,bdrq,bdrp,xtime,dtbdr,lmoist,
     &                     nstep,nbdsteps,ktimax,nij2,nk,nspong,
     &                     iunitb,ktaubd,0)
      endif
c
c  Compute new values determined from explicit scheme
      call Astepf1 (pun,pum,puten,dti,ni,nj,nk,0)
      call Astepf1 (pvn,pvm,pvten,dti,ni,nj,nk,0)
      call Astepf1 (ptn,ptm,ptten,dti,ni,nj,nk,1)
      call Astepf1 (psn,psm,psten,dti,ni,nj, 1,1)
      if (lmoist) call Astepf1 (pqn,pqm,pqten,dti,ni,nj,nk,1)
      if (lground) call Asteptg (tgn,tgm,tgten,snowc,tgb,dti,ni,nj,
     &                           lsnow,iland,xland)
c
c  Adjoint of adiabatic tendencies
      call Atendc (ps,pu,pv,pt,sdot,pd,puten,pvten,ptten,psten,u,
     &             v,t,puf,pufx,pvf,pvfx,ux,vx,div,sdotd,omega,
     &             alogp,pub,pvb,ptb,psb,alogpb,divb,sdotb,omegab,
     &             pufb,pufxb,pvfb,pvfxb,ub,vb,uxb,vxb,tb,pdb,
     &             sdotdb,chi,fmapx,fmapd,fmapx2,fmapd2,cf,sigf,
     &             sigh,dsig,ptv,ptvb,cpmr,cpmrb,q,pqten,qb,
     &             work1,work2,lmoist,lbdryp,ep1,ecpm,R,Cp,ptop,
     &             dx,ni,nj,nk,pqb,pq,nkp1,td,tdb,tv,tvb,phi,
     &             phib,grav,work9,.true.)
c
c  Adjoint of some contributions to the ground temperature tendency
      if (lground) 
     &   call Atgten (tgten,tgm,pqm,ptm,glw,gsw,cldfr,psm,tm,qm,prw,
     &                qs,work1,tgmb,pqmb,cldfrb,psmb,tmb,qmb,prwb,qsb,
     &                work2,gemissv,gheatcpr,iland,work3,stbolt,
     &                hdeepfac,dsig,sigh,ptop,grav,lmoist,lcloud,
     &                qsatvp,nqsatvp,ni,nj,nk,nstep,radfreq,klevcld,
     &                klevrad,xlat,xlon,galbedo,degrad,solcon,lprint,
     &                xtime,iyear0,imonth0,iday0,xhour0,ttime0,dtm,
     &                iawrite)
c
c  Adjoint of diffusion, boundary damping, and PBL forcing
      call Atendm (pum,pvm,ptm,pqm,psm,tgm,puten,pvten,
     &             ptten,pqten,um,vm,tm,qm,psten,tgten,pdm,umb,
     &             vmb,tmb,
     &             qmb,psmb,tgmb,pkhb,pdmb,xflux,xfluxb,xkc,
     &             gheatcpr,work1,work2,sigh,sigf,dsig,bfact,
     &             qsatvp,pqmb,theta,thetab,work1,cdragx,cdragd,
     &             pumb,pvmb,ptmb,xpprp,gmavail,xkap,Cp,
     &             dx,dti,ptop,R,grav,hlatent,balpha,bbeta,lground,
     &             lmoist,nqsatvp,ni,nj,nk,nkp1,nkpx,nspong,lbdry,
     &             lbdryp)
c
c  Adjoint of sponge forcing
      if (lbdryx) then
         if (lbdry) 
     &      call ALbforce (puten,pvten,ptten,pqten,psten,
     &                     bdru,bdrv,bdrt,bdrq,bdrp,
     &                     bfact,balpha,bbeta,nstep,nbdsteps,ktimax,
     &                     lmoist,ni,nj,nk,nspong,nij2)
         call ALbsete (pu,pv,pt,pq,ps,bdru,bdrv,bdrt,bdrq,bdrp,
     &                 nstep,nbdsteps,lmoist,ni,nj,nk,nspong,nij2,
     &                 ktimax,1,0)
         call ALbsete (pum,pvm,ptm,pqm,psm,bdru,bdrv,bdrt,bdrq,bdrp,
     &                 nstep,nbdsteps,lmoist,ni,nj,nk,nspong,nij2,
     &                 ktimax,0,0)
      endif
c
c  Reset boundaries to zero
      call Cbdryx (pu,pv,pt,ps,ni,nj,nk)
      call Cbdryx (pum,pvm,ptm,psm,ni,nj,nk)
      if (lmoist) then
         call Cbdry0 (pq,ni,nj,nk,1)
         call Cbdry0 (pqm,ni,nj,nk,1)
      endif
c  
      if (iawrite.eq.1) then
         call Naddv (pu,pu,pum,nijk)
         call Naddv (pv,pv,pvm,nijk)
         call Naddv (pt,pt,ptm,nijk)
         call Naddv (ps,ps,psm, nij)
         if (lmoist) call Naddv (pq,pq,pqm,nijk) 
         if (lground) call Naddv (tg,tg,tgm,nij)
         if (linitial) then 
            if  (nstep.eq.ktimax) llast=.t.
            call Annmi (pu,pv,pt,pq,ps,u,v,t,q,pd,puten,pvten,ptten,
     &                  psten,chi,ux,vx,td,tv,ptv,sdot,sdotd,
     &                  puf,pvf,pufx,pvfx,cpmr,omega,phi,fmapx,fmapd,
     &                  fmapx2,fmapd2,topog,cf,sigf,sigh,dsig,
     &                  alogp,alogpb,divb,grav,R,Cp,ep1,ptop,div,
     &                  pum,phi,pvm,work1,work2,work3,work9,chi,
     &                  einversT,ematrxiT,gmatrixp,gmatrxtT,gmatrxrT,
     &                  edepth,pweight,work1,tpmeanb(nkp1),clat,ecpm,
     &                  dx,tpratio,nnmis,ni,nj,nk,nkp1,lmoist,lprint,
     &                  nvarsi,cfldsi,workr,iunitir,cfldsrw,ifldsrw,
     &                  ldaccess,nnvars,nworkr,fields3,fields2,
     &                  ndimfld3,ndimfld2,puten,pvten,pum,pvm,ptm,
     &                  psm,pub,pvb,ptb,pqb,psb,pdb,sdotb,
     &                  sdotdb,omegab,pufb,pufxb,pvfb,pvfxb,ub,vb,uxb,
     &                  vxb,tb,cpmrb,qb,tvb,phib,ptvb,tdb,
     &                  pqmin,rvarsi,nrvars,cfldsr,ifldsr,lbdryp,llast)
             if (llast.and.lbdryp)    
     &       call ALbsete (pu,pv,pt,pq,ps,bdru,bdrv,bdrt,bdrq,bdrp,
     &                     nstep,nbdsteps,lmoist,ni,nj,nk,nspong,nij2,
     &                     ktimax,1,0)
             call Cbdryx (pu,pv,pt,ps,ni,nj,nk)
             if (lmoist) call Cbdry0 (pq,ni,nj,nk,1)
ctlmadj remove write
         endif
      endif 
c
c
c
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
    2 continue
      endif ! check on ktimax
    
c
c             END OF ONE TIME STEP INTEGRATION
c 
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c	END OF ADJOINT
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C	in this last step of catx, the result is made nondimensional, 
c	and the result of the AQ check is printed out. 
c 
      call PDECUPL (pu,pv,pt,pq,ps,u,v,t,q,
     &              pub,pvb,ptb,pqb,psb,ub,vb,tb,qb,
     &              pd,pdb,lmoist,ni,nj,nk,4)
c
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      call Nrsetc (raincv,nij,0.)
      call Nrsetc (raincv2,nij,0.)
      call Nrsetc (rainncv,nij,0.)
      call Nrsetc (rainncv2,nij,0.)
      call Xnorms3 (u,v,t,ps,q,tg,xstar,xfull,xnorm0,lmoist,
     &              sigf,sigh,dsig,grav,R,Cp,clat,clon,dx,
     &              Redepth,RmatrxiT,Rmatrixp,RmatrxrT,Rpweight,
     &              qmatrixT,
     &              Rtpmeanb,work1,work2,work3,work4,work5,work6,
     &              cf,fmapx,fmapd,ni,nj,nk,nxstar,nfull,itype) 
c
c
        xnormd = xnorm0 - xnorm1
C
	WRITE (LUNIT1,'(A,E40.20)') '==== CHECK OF INNER PRODUCTS ===='
	WRITE (LUNIT1,'(A,E40.20)') ' XNORM1 (Lx)^2     = ', xnorm1
	WRITE (LUNIT1,'(A,E40.20)') ' XNORM0 x^T L^T Lx = ', xnorm0
	WRITE (LUNIT1,'(A,E40.20)') ' XNORMD (0-1)      = ', xnormd
	WRITE (LUNIT1,'(A,E40.20)') '================================='
c
c user set
c replace return with stop here to really do just 1 iteration
	RETURN
	END
c
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c

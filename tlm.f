#QSUB -eo -lt 300 -lT 330 -lM 6.0Mw -q reg


#
# This deck is TLM set to run from SV file input.
# It includes extra diagnostics regarding partitioning of the
#  E norm into components.
# It also includes options on filtering R or G modes or Gamma modes.
#

# run TLM on case W1; look at Var norm
#   set lmoist
#   set charsw(5 and 6)
#   select SV vdim file
#   set printed output name

cd $TMPDIR
csh -x << 'JOBEND' >& log
setenv NCPUS 1

rm lnf.* 
ja
 
cat << 'EOFD' > lnf.f


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
      PROGRAM MM4TLM
c 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This is the driver for the Tangent Linear Model (TLM) component 
c     of version 2A of the Mesoscale Adjoint Modeling System (MAMS2).
c
c     This code (vesion 2A) was developed by Ronald Errico
c     and Kevin Raeder at NCAR (January 1997).  Prospective users
c     should inform Ronald Errico by e-mail or in writing of any
c     intentions to use this software and agree to reference any appro-
c     priate publications of the developers and acknowledge use of this
c     software in any papers.  Furthermore, pre-prints and reprints of
c     any papers produced using this software should be sent to R.
c     Errico at NCAR.
c
c     For a description of MAMS version 1A see NCAR Technical Note
c     #TN/410+IA, 1994, by R. M. Errico, K. Raeder, and T. Vukicevic.
c     Version 2A is very differeant than 1A.  For a summary of the 
c     changes, contact R. Errico.
c
c---------------------------------------------------------------------
c
      logical lmoist,lground,lrain,lcloud,lsnow,lbdry,lbdryx,lsplit
     &,       lbdryp,lbdsame,llinear,lcoefout,linitial,linitwr,ldrycon
c 
c---------------------------------------------------------------------
c                  USER DEFINED PARAMETERS :
c---------------------------------------------------------------------
c
cqqqqqqqqqqq
      logical lrfilt, lgfilt, lmode0
c      parameter (lrfilt=.false., lgfilt=.false., lmode0=.false.) ! full
c       parameter (lrfilt=.true., lgfilt=.false., lmode0=.false.) ! G only
       parameter (lrfilt=.false., lgfilt=.true., lmode0=.true.) ! R only 
c  lgfilt is for removing g-waves from pert
c  lrfilt is for removing r-waves from pert
c  lmode0 is for removing gamma mode from pert (should change tpmeanb)
c
c
c  grid dimensions
      parameter (ni=46, nj=61, nk=10, nspong=5)
c
c  timing information (units: seconds, minutes, minutes)
      parameter (dt=60.*5., timax=24.*60., dtbdr=60.*12.)
c
c  output file information (minutes)
      parameter (toutreg=12.*60., ipackw=2)
c
c  basic state input file information (units=mins)
      parameter (ftimed0=0., tupdated=15.)
      parameter (ftimem0=0., tupdatem=15.)
c
c  moisture selection
      parameter (lmoist=.f., nqpoints=960)
c
c  perform dry convective adjustment
      parameter (ldrycon=.t.)
c
c  ground temperature physics selection
      parameter (lground=.t.,  lcloud=.f., lsnow=.f.)
c
c  consider boundary perturbations?
      parameter (lbdryp=.f., lbdsame=.f., ftimeb0=0.)
c
c  selection of initial conditions in routine = Finit
      parameter (ipert=1,  pert=1.,  tskips1=0., tskips2=0. )
c
c  specification of implicit nonlinear normal mode initialization
      parameter (linitial=.t., linitwr=.f.) 
c
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
     &,          nij=ni*nj, nijk=nij*nk, nij2=(ni+nj)*2, nkp2=nk+2)
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
cc
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
c                     ground temperature are to be considered
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
c     pert ...........scaling factor for initial condition
c     tskips1.........selects time on file to be used if initial
c                     conditions to be read from a file on iuniti(1)
c     tskips2.........selects time on file to be used if initial
c                     conditions to be determined as differences
c                     between values on two files iuniti(1)-iuniti(2)
c     linitial....... .t. if implicit nonlinear normal mode 
c                     initialization to be done at start of forecast.
c     linitwr........ .t. if innmi result to be written to file (only if
c                     linitial=.t.)
c
c
c  The terms half and full levels respectively refer to data and 
c     interface levels on the vertical grid.
c  The terms dot and cross grids refer to the horizontal grids on
c     which pu, pv and pt, pq, ps, tg are defined, repectively
c
c  Principal dynamic fields (perturbations)
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
c  div .....horizontal mass flux divergence
c  theta ...theta field on half levels and surface
c  xpprp....factor for perturbing pkh or 1/pkh
c  pkh .....factor for converting t to theta on half levels and surface
c  phi .....hydrostatic geopotential
c  alogp ...natural log of p at each data level and at the surface
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c                     DIMENSION STATEMENTS
c
c dimension grid arrays (see subroutine = Nsetup for definitions)
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
      integer msteps(nk), nnmis(4)
c
c dimension arrays for R-norm
      PARAMETER ( ITRC = ni-3 , JTRC = Nj-3, KTRC = nk )
      PARAMETER ( ITRCQ = 1 , JTRCQ = 1 , KTRCQ = 1)
      parameter (nxstar=itrc*jtrc*ktrc, nfull=nxstar, normtype=4)
      real Rematrix(nk,nk), Rinverse(nk,nk), Rmatrixt(nk,nk)
     &,    Rmatrixp(nk), Rtmatrix(nk,nk), Rpmatrix(nk), Redepth(nk)
     &,    Rtpmeanb(nkp1), Rmatrixi(nk,nk), Rmatrixr(nk,nk)
     &,    Rpweight(nkp2,3), Rmsteps(nk)
     &,    XSTAR (nxstar), XFULL (nfull)
      real qmatrix(nk,nk), qinverse(nk,nk)
      common /cctsmatx/ ctsmatx(11), bnbmatx(10) !qqqqq
c
c dimension boundary fields
c
      real bdru(nij2,nk,nspong,2), bdrv(nij2,nk,nspong,2)
     &,    bdrt(nij2,nk,nspong,2), bdrq(nij2,nk,nspong,2)
     &,    bdrp(nij2,nspong,2)
     &,    field3b(nijk,4), field2b(nij,2)
c
c dimension work arrays 
c (nfield3m is maximum size of nfield3 on any input file)
      parameter (nwork1=nij*nkp1, nfield3m=12
     &,          nworkr=2+(nij+2)*nfield3m )
      real work1(ni,nj,nkp1), work2(ni,nj,nkp1), work3(ni,nj,nkp1)
     &    ,work4(ni,nj,nk), work5(ni,nj,nk), work6(ni,nj,nk)
     &    ,work7(ni,nk,7), work8(ni,nj,13), workr(nworkr)
     &    ,worka(nijs,2), workb(nijs,6), workc(nijkx,12)
     &    ,work9(ni,nj,nkp2,2), work10(ni,nkp1,14)
c
c
c  The following equivalence is optional to reduce memory required.
c  If any of these fields are to be output, they should not be
c  equivalenced here.
c
      equivalence (work1,workr,work7,workq,sdotd,ptv,work9,work10)
     &,           (work2,workc      ,sdotdb,tv)  ! qqq remove omega
     &,           (work3,workb      ,pufx,pufb,ux,td,theta)
     &,           (work4,worka      ,pvfx,pvfb,vx,tdb,thetab)
     &,           (work5            ,pufxb,uxb,puf,cpmr,xflux)
     &,           (work6,work8      ,pvfxb,vxb,pvf,xfluxb)
     &,           (pun,puten),(pvn,pvten),(ptn,ptten),(pqn,pqten)
c
c  If the above optional equivalencing is not used, then use
c
c      equivalence (work6,work8), (work7,work1,workr,workq),  
c     &            (work2,workc), (work3,workb), (work4,worka)
c
      equivalence (field3b,u,um,qcoefsw)
     &,           (field3b(1,2),v,vm), (field3b(1,3),t,tm) 
     &,           (field3b(1,4),q,qm), (field2b,pd,pdm)
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
      logical lcheck,lprint
      data lcheck/.t./, lprint/.t./
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
     &,    iunitiw/80/, iunitir/81/, iunita1/90/
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
c  iunitb is unit for lateral boundary file (unless lbdsame=.t.)
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
c  nbsizea specifies i/o buffer size (words) for a dry conv. file
c
      data charsw(1)/'Lib.2A3'/
      data charsw(2)/'no initial file'/
      data charsw(3)/'no perturbation file'/
      data charsw(4)/'no separate boundary file'/
      data charsw(5)
     &     /'/RAEDER/ADJOINT/MSTACC/EXPL/H48-72/L2B1/NLMBS5Mx3/'/
      data charsw(6)/'SV/SPECTS/TLM/VNORM/W1/DRY/SV1/RONLYPERT/'/
c      data charsw(6)/'SV/SPECTS/TLM/VNORM/W1/DRY/SV1/FULLPERT/'/
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
c  BLOCK FOR TIME SERIES OUTPUT
c  set arrays for time series output
      logical ltseries
      parameter (nfieldts=9, npointts=25, numnorms=4, maxstepe=4000)
      parameter (ntcharsw=n1charsw+8)
      character*(ntcharsw) charswts
      character*8 cnorms(numnorms)
      character*8 cfieldts(nfieldts), cfunitts(nfieldts)
      integer ilocts(npointts), jlocts(npointts), klocts(npointts)
      real savenorm(maxstepe,numnorms)
      data cfieldts /'   ps   ','  omega ','    T   ','    T   ',
     &               '   T    ','   q    ','  Tg    ',' raincv ',
     &               ' rainncv'/
      data cfunitts /'   cb   ','  cb/s  ','    K   ','    K   ',
     &               '    K   ','   g/g  ','    K   ','   cm   ',
     &               '   cm   '/
      data cnorms   /'  EDRY  ','   EQ   ','  ESUM  ','  RNORM '/
      data ltseries /.true./
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
c  Print some run parameters
      if (lprint) then
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
         print *,'   linitial=',linitial,'  linitwr=',linitwr
     &          ,'  ldrycon=',ldrycon
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
c  Check some parameters
      if (toutreg.gt.0.) call Nchktime (toutreg,dtm,'toutreg')
      call Nchktime (tupdated,dtm,'tupdated')
      if (lmoist) call Nchktime (tupdatem,dtm,'tupdatem')
c
c
c
c                      INITIALIZE MODEL
c
c  Acquire and Open basic state file, Read and Get required
c  information from header.
c
      if (lpshell) call Nacquire (charsw(5),'D',1,lprint,nbsized,
     &                            iunitd)
      if (ldaccesi) call Nopen (0,iunitd,lprint)
      call Nreadh (nvarsi,lvarsi,rvarsi,cfldsi,charsi,fieldh,workr,
     &             cfldsrw,cfldsr,ifldsrw,ifldsr,ni,nj,nk,
     &             nij,ir1,iunitd,ldaccesi,.t.,lprint,
     &             nnvars,nlvars,nrvars,nworkr)
      ir1i=ir1  ! save for use in Fnnmi routine    
      call Ngethead (sigf,tpmeanb,ptop,dx,clat,clon,xhour0,time0,
     &               iyear0,imonth0,iday0,nsplit,radfreq,msteps,
     &               lsplit,lbdry,nnmis,tpratio,cfldsi,
     &               charsi,nk,nkp1,rvarsi,lvarsi,nvarsi,lprint,
     &               nnvars,nrvars,nlvars,' input')
      nbsdataf=int((tupdated+0.5*dtm)/dtm)
c
cqqqqqqqqqqqq 
c      ptop=.01
c      tpmeanb=270.
c      tpmeanb(nk+1)=100.
cqqqqqqqqqq

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
c read basic state data and required coefficients for initial time
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
     &              pqmin,'TLM')
      if (lmoist) then
         call Bsdataq0 (nvarsq1,rvarsq1,'C',iunitm1,ktaum1,charsw(5),
     &                  nnvarsq,nrvarsq,nq2dim,nq3dim,2,4,
     &                  qcoefs1,ipnoncv,nbsizec,lprint,ldaccesi,lpshell,
     &                  nnvarsq,nrvarsq,nij,nq4dim,nfilem1)
         call Bsdataq0 (nvarsq2,rvarsq2,'N',iunitm2,ktaum2,charsw(5),
     &                  nnvarsq,nrvarsq,nq2dim,nq2dim,2,1,
     &                  qcoefs3,ipnoncv,nbsizen,lprint,ldaccesi,lpshell,
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
      endif ! test on lmoist
      if (ldrycon) then
         call Bsdataq0 (nvarsa1,rvarsa1,'A',iunita1,ktaua1,charsw(5),
     &                  nnvarsq,nrvarsq,nq2dim,nq3dim,1,1,
     &                  dcoefs1,ipdconv,nbsizea,lprint,ldaccess,lpshell,
     &                  nnvarsq,nrvarsq,nij,nq6dim,nfilea1)       
         call Bsdataq1 (dcoefs1,ipdconv,rvarsa1,nvarsa1,iunita1,'A',
     &                  ktaua1,nfilea1,tupdatem,ftimem0,workq,
     &                  qcoefsw,work6,ndmaxxx1,
     &                  charsw(5),0.,lprint,lpshell,ldaccesi,'TLM',
     &                  nnvarsq,nrvarsq,nij,nq6dim,nqdimw)
      endif ! test on ldrycon
      nbsdataq=int((tupdatem+0.5*dtm)/dtm)
c
c read or set initial data 
c      call Finit (pu,pum,pv,pvm,pt,ptm,pq,pqm,ps,psm,tg,tgm,
c     &            pub,pvb,ptb,pqb,psb,pdb,tgb,workr,work2,
c     &            fields3,fields2,field3c,field2c,ifldsrw,ifldsp,
c     &            cfldsrw,cfldsp,charsw(2),lmoist,lground,lpshell,
c     &            ldaccesi,lprint,nvarsp,charsi,ipert,pert,tskips,
c     &            ni,nj,nk,nbsized,iuniti,nnvars,nworkr,ndimfld3,
c     &            ndimfld2,raincv,rainncv,4,2)
c       pertx=1. !QQQQQQQQQq
ccQQQQQQ
      call FinitGG (pu,pv,pt,ps,pub,pvb,ptb,psb,u,v,t,pd,q,qb,
     &                  ub,vb,tb,pdb,tg,pq,field3c,work1,pert,
     &                  work6,ipert,iuniti,ni,nj,nk,lmoist,pertx)
      tg=0.
      tgm=0.
c      write (13) pu
c      write (13) pv
c      write (13) pt
c      write (13) pq
c      write (13) ps


cqqqqqqqqqqqq
c
c  these are new matrices for r-norm (isothermal T)
      CALL xnorms0 (ni,nj,nk,itrc,jtrc,ktrc,itrcq,jtrcq,ktrcq,normtype)
      call Nrsetc (Rtpmeanb,nk,270.)
      Rtpmeanb(nkp1)=100.0
      call Sinitial (Rematrix,Rinverse,Rmatrixt,Rmatrixp,Rtmatrix,
     &               Rpmatrix,Redepth,Rtpmeanb,Rmsteps,
     &               ptten,psten,puten,pvten,sdot,work3,work2,
     &               Rmatrixi,Rmatrixr,Rpweight,work5,sigf,sigh,dsig,
     &               R,Cp,ptop,dx,dt,nix,njx,nijkx,nk,nkp1,
     &               1,.false.,.true.,.false.)  !  no write
c
c  compute matrices for spectral q norm
      if (lmoist) then 
        call QNmatrix (qmatrix,qinverse,qtotvar,tpmeanb,sigh,dsig,
     &            ptop,qsatvp,work1,work1(1,1,2),work2(1,1,2),work2,
     &            nqsatvp,nk,lprint)
      endif
c
cqqqqqqqqqqq
c  filter g or r modes 
      if (lrfilt.or.lgfilt) then  
        if (lrfilt) then ! save initial fields, to subtract later 
          pum=pu
          pvm=pv
          ptm=pt
          psm=ps
        endif
        call PDECUPL (pu,pv,pt,pq,ps,u,v,t,q,
     &                    pub,pvb,ptb,pqb,psb,ub,vb,tb,qb,
     &                    pd,pdb,lmoist,ni,nj,nk,1)
      call Varnorm (u,v,T,ps,ni,nj,nk,'before filter')
        call PVfs2pv (xfull,u,v,t,ps,work4,work5,work3,
     &                 work1,work2,Redepth,Rinverse,Rmatrixp,
     &                 Rmatrixt,Rtpmeanb(nk+1),clat,dx,ni,nj,nk,
     &                 itrc,jtrc,ktrc)

c       print *,' '
c       print *,'TEST XFULL 1'
c       print ('(1x,1p10e12.2)'),(xfull(i),i=1,nfull,20)
     
        call PVpv2fs (u,v,t,ps,xfull,work1,work2,work3,work4,
     &                 work5,work6,Redepth,Rmatrixi,Rmatrixp,
     &                 Rmatrixr,Rpweight,Rtpmeanb(nk+1),clat,dx,
     &                 ni,nj,nk,itrc,jtrc,ktrc)
        call PDECUPL (pu,pv,pt,pq,ps,u,v,t,q,
     &                    pub,pvb,ptb,pqb,psb,ub,vb,tb,qb,
     &                    pd,pdb,lmoist,ni,nj,nk,2)

c       print *,' '
c       print *,'TEST XFULL 2'
c       call PVfs2pv (xfull,u,v,t,ps,work4,work5,work3,
c     &                 work1,work2,Redepth,Rinverse,Rmatrixp,
c     &                 Rmatrixt,Rtpmeanb(nk+1),clat,dx,ni,nj,nk,
c     &                 itrc,jtrc,ktrc)
c       print ('(1x,1p10e12.2)'),(xfull(i),i=1,nfull,20)
        if (lrfilt) then ! remove r waves from initial fields
          pu=pum-pu
          pv=pvm-pv
          pt=ptm-pt
          ps=psm-ps
c          write (13) pu
c          write (13) pv
c          write (13) pt
c          write (13) ps
c          stop 'TEST'
        endif  ! test on lrfilt
      endif  ! test on filter


c TIME SERIES BLOCK
      if (ltseries) then
cqqqqq
      if (lmode0) then 
      do i=1,ni-1
      do j=1,nj-1
c        ps(i,j)=0.
      do k=1,nk
c        pt(i,j,k)=psb(i,j)*gmatrixr(k,4)
      enddo
      enddo
      enddo
c
      rxxx=1./(ctsmatx(nkp1)+1.)
      call Fdecup (t,pt,ps,tb,psb,ni,nj,nk,1)
      do i=2,ni-2
      do j=2,nj-2
        rxxx1=-ps(i,j)
        do k=1,nk
          rxxx1=rxxx1+ctsmatx(k)*t(i,j,k)
        enddo
        deltaps=rxxx1*rxxx
          ps(i,j)=ps(i,j)+deltaps
        do k=1,nk
          deltat=-bnbmatx(k)*deltaps
          pt(i,j,k)=pt(i,j,k)+deltat*psb(i,j)+deltaps*tb(i,j,k)
          pu(i,j,k)=pu(i,j,k)+deltaps*ub(i,j,k)
          pv(i,j,k)=pv(i,j,k)+deltaps*vb(i,j,k)
        enddo
      enddo
      enddo

      call Fdecup (t,pt,ps,tb,psb,ni,nj,nk,1)
      do i=2,ni-2,20
      do j=2,nj-2,20
        do k=1,nk
          work1(i,j,k)=gmatrixp(k)*ps(i,j)
          do n=1,nk
            work1(i,j,k)=work1(i,j,k)+gmatrixt(k,n)*t(i,j,n)   
          enddo
        enddo
        print *,'ZZ ',i,j,(work1(i,j,k),k=1,nk)
      enddo
      enddo

cqqqqq
      print *,' gmatrixr'
      gggsum=0.
      do k=1,nk
       ggg=gmatrixr(k,1)*gmatrixr(k,1)
       gggsum=gggsum+ggg
       print *,k,gmatrixr(k,1),gggsum
      enddo 
      gggsum=0.1*0.5*gggsum*Cp/270.
      edepr=0.5/edepth(1)
      print *,'gggsum=',gggsum,edepth(1),edepr,redepth(1)

c
      print *,' '
      print *,'edepth, 0.5/edepth'
      do k=1,nk
        edepr=0.5/edepth(k)
        print *,k,edepth(k),edepr
      enddo
      print *,' '
c
      do i=1,nk
      do j=1,nk
      qmatrix(i,j)=gmatrixt(i,j)/edepth(i)
      enddo
      enddo  
      do i=1,nk
      do j=1,nk
      qinverse(i,j)=0.
      do k=1,nk
      qinverse(i,j)=qinverse(i,j)+gmatrixt(k,i)*qmatrix(k,j)
      enddo
      enddo  
      enddo
      do i=1,nk
      do j=1,nk
        qinverse(i,j)=qinverse(i,j)*10.*270/(Cp*tpmeanb(nk+1)**2)
      enddo
      enddo
         call Sprntm (qinverse,nk,nk,'X energy')
c
      do i=1,nk
      do j=1,nk
      qinverse(i,j)=0.
      do k=1,nk
      qinverse(i,j)=qinverse(i,j)+gmatrixr(i,k)*gmatrixt(k,j)
      enddo
      enddo  
      enddo
         call Sprntm (qinverse,nk,nk,'ONE     ')
      endif ! test on lmode0

c

         call TENORMZ (pu,pv,pt,ps,u,v,t,pd,ni,nj,nk,e,
     &                    pub,pvb,ptb,psb,ub,vb,tb,pdb,Cp,R,dsig,
     &                    pqb,qb,pq,q,eq,esum,lmoist)
         call Xnorms2a (u,v,t,ps,q,tg,xfull,lmoist,rnormval,
     &               sigf,sigh,dsig,grav,R,Cp,clat,clon,dx,
     &               Redepth,Rinverse,Rmatrixp,Rmatrixt,Rpweight,
     &               qinverse,
     &               Rtpmeanb,work1,work2,work3,work4,work5,work6,
     &               cf,fmapx,fmapd,ni,nj,nk,nfull,nxstar,normtype)
         print *,'EDRY,EQ,ESUM,Rnorm= ',e,eq,esum,rnormva
      call Varnorm (u,v,T,ps,ni,nj,nk,'before init')
         call Emodesub (xfull,u,v,t,ps,work3
     &                    ,work1,work2,Redepth,Rinverse,Rmatrixp
     &                    ,Rmatrixt,Rtpmeanb(nk+1),grav,ni,nj,nk
     &                    ,itrc,jtrc,ktrc)

         nesaves=0
         savenorm(1,1)=ipert
         savenorm(2,1)=pertx
         savenorm(3,1)=e
         savenorm(1,2)=ipert
         savenorm(2,2)=pertx
         savenorm(3,2)=eq
         savenorm(1,3)=ipert
         savenorm(2,3)=pertx
         savenorm(3,3)=esum
         savenorm(1,4)=ipert
         savenorm(2,4)=pertx
         savenorm(3,4)=rnormval
      endif
c  balance initial condition using implicit nnmi
c
cqqqqqqqqq
C COMPUTE AND PRINT PVORT, DIVG, AGEO STATS FOR COMPUTING TIME SCALE
c BASED ON MODE TENDENCY TO MEASURE BALANCE
c
c  compute vertical mode amplitudes of divergence and vorticity
      call PDECUPL (pu,pv,pt,pq,ps,u,v,t,q,
     &                    pub,pvb,ptb,pqb,psb,ub,vb,tb,qb,
     &                    pd,pdb,lmoist,ni,nj,nk,1)


      call Sproject (work2,u,einverse,ni,nj,nk,nk,nk,0,2,.true.)
      call Sproject (work3,v,einverse,ni,nj,nk,nk,nk,0,2,.true.)
      call Bpuvm (work2,work3,work2,work3,fmapd,ni,nj,nk)
      call Bdiv (work4,work2,work3,fmapx2,dx,ni,nj,nk)
      call Bvort (work5,work2,work3,fmapx2,dx,ni,nj,nk)
c
c  compute vertical mode amplitudes of pseudo-geopotential
      call Sproject (work2,ps,gmatrixp,ni,nj,nk,1,nk,1,1,.true.)
      call Sproject (work2,t,gmatrixt,ni,nj,nk,nk,nk,1,1,.false.)
      call Nmultv (work2,nijk,1./tpmeanb(nk+1))
c
c  compute vertical modes of avort tendency
      pi4=atan(1.)
      f0=16.*pi4*sin(clat*pi4/45.)/(60.*60.*24)
      print *,' '
      print *,' f0,clat=',f0,clat
      print *,' STATS BEFORE FNNMI '
      print *,' '
      call Iagvort  (work1,work5,work2,work3,fmapx2,dx,f0,ni,nj,nk)
c
c  print diagnostics of tendencies of pvort, divg, avort in mode space
      if (lprint) then
         call Ipotvort (work3,work5,work2,edepth,f0,ni,nj,nk)
         call Idiagns (work3,work4,work1,edepth,work2,ni,nj,nk,2)
      endif
cqqqqqqqqqq end
       call Fmoddiag (pu,pv,pt,pq,ps,u,v,t,q,pd,puten,pvten,ptten,chi,
     &               psten,ux,vx,td,tv,ptv,alogp,work9,sdot,sdotd,
     &               puf,pvf,pufx,pvfx,cpmr,omega,phi,fmapx,fmapd,
     &               fmapx2,fmapd2,topog,cf,sigf,sigh,dsig,grav,R,Cp,
     &               ep1,ptop,div,pum,phi,pvm,ptm,work1,work2,
     &               work3,chi,rinverse,ematrixi,rmatrixp,rmatrixt,
     &             gmatrixr,redepth,Rpweight,work1,rtpmeanb(nkp1),clat,
     &               ecpm,dx,tpratio,nnmis,ni,nj,nk,nkp1,lmoist,lprint,
     &               nvarsi,cfldsi,charsw,workr,iunitir,cfldsrw,ifldsrw,
     &               lpshell,ldaccess,nnvars,nworkr,fields3,fields2,
     &               ndimfld3,ndimfld2,puten,pvten,pum,pvm,ptm,
     &               psm,pub,pvb,ptb,pqb,psb,pdb,sdotb,
     &               sdotdb,omegab,pufb,pufxb,pvfb,pvfxb,ub,vb,uxb,
     &               vxb,tb,cpmrb,qb,tvb,phib,ptvb,tdb,alogpb,divb,
     &               nbsized,pqmin,ir1i,work8,cfldsr,ifldsr,.true.,
     &                  Rtpmeanb,xfull,nfull,itrc,jtrc,ktrc,normtype)


c
c  Read perturbed lateral boundary perturbations if requested
c
      if (lbdryp) then
         ktaubd=0
         nbdsteps=int((dtbdr*60.+0.5*dt)/dt)
         if (lbdsame) then   ! Initial and boundary files identical
            iunitb=iuniti(1)
            charsw(4)=charsw(2)
         else  ! Additional file needs to be acquired and read
            if (lpshell) call  Nacquire (charsw(4),'D',1,lprint,
     &                                   nbsized,iunitb)
            if (ldaccesi) call Nopen (0,iunitb,lprint)
            call Nreadh (nvarsb,lvarsb,rvarsb,cfldsb,charsi,work2,workr,
     &                   cfldsrw,cfldsr,ifldsrw,ifldsr,ni,nj,nk,
     &                   nij,ir1,iunitb,ldaccesi,.f.,lprint,
     &                   nnvars,nlvars,nrvars,nworkr)
         endif
         call NLbread (bdru,bdrv,bdrt,bdrq,bdrp,field3b,field2b,
     &                 workr,-1.,xtimef,lmoist,cfldsb,rvarsb,
     &                 nvarsb,nnvars,nrvars,dtbdr,ktaubd,nij2,ni,nj,
     &                 nk,nijk,nworkr,nspong,ftimeb0,ldaccesi,
     &                 lprint,iunitb)
         call NLbread (bdru,bdrv,bdrt,bdrq,bdrp,field3b,field2b,
     &                 workr,0.0,xtimef,lmoist,cfldsb,rvarsb,
     &                 nvarsb,nnvars,nrvars,dtbdr,ktaubd,nij2,ni,nj,
     &                 nk,nijk,nworkr,nspong,ftimeb0,ldaccesi,
     &                 lprint,iunitb)
         call NLbsete (pu,pv,pt,pq,ps,bdru,bdrv,bdrt,bdrq,bdrp,
     &                 0,nbdsteps,lmoist,ni,nj,nk,nspong,nij2,0)
      endif
c
c write output data header and initial fields
      call Nsethead (rvarsw,lvarsw,nvarsw,toutreg,ptop,dx,clat,clon,
     &               xhour0,time0,dt,dtbdr,radfreq,tupdated,sigf,
     &               tpmeanb,lmoist,lground,lrain,lcloud,lsnow,
     &               lbdryp,lbdsame,lsplit,llinear,lcoefout,linitial,
     &               ldrycon,igtypew,iyear0,imonth0,iday0,nsplit,
     &               msteps,nnmis,ifldsw,npkappa,nqsatvp,nfldswt,
     &               ni,nj,nk,nkp1,nspong,n1charsw,ndcharsw,nnvars,
     &               nlvars,nrvars,nnftimes,iformhw,iformtw,ipackw,
     &               tpratio)
      call Nsetoutp (ktimout,ktimax,nvarsw,dt,toutreg,timax,tout,
     &               nnvars,ntoutpt,nnftimes,nnftimeo)
      if (lpshell) call Nassign (lprint,iunitw,itapnumd,'D')
      call Nopen (nvarsw(2),iunitw,lprint)
      call Nwriteh (nvarsw,lvarsw,rvarsw,cfldsw,charsw,workr,nij,
     &              irw1,iunitw,fieldh,cfldsrw,ifldsrw,lprint,
     &              nnvars,nrvars,nlvars,nworkr)
      irw=irw1
      ktauout=0
      if ( (toutreg.gt.0).or.
     &    ((toutreg.lt.0).and.(0.eq.nvarsw(ntoutpt+ktauout)))) 
     &   call Nwritet (fields3,fields2,workr,0.,nvarsw,cfldsrw,
     &                 cfldsw,ifldsrw,ifldsw,nij,nk,ktauout,irw,
     &                 iunitw,lprint,nnvars,nworkr)
c 
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
cqqqqqqqqqqq
          write (13) pu
          write (13) pv
          write (13) pt
          write (13) ps


cqqqqqqqqq
C COMPUTE AND PRINT PVORT, DIVG, AGEO STATS FOR COMPUTING TIME SCALE
c BASED ON MODE TENDENCY TO MEASURE BALANCE
c
c  compute vertical mode amplitudes of divergence and vorticity
      call PDECUPL (pu,pv,pt,pq,ps,u,v,t,q,
     &                    pub,pvb,ptb,pqb,psb,ub,vb,tb,qb,
     &                    pd,pdb,lmoist,ni,nj,nk,1)

      call Varnorm (u,v,T,ps,ni,nj,nk,'after init')
      call Sproject (work2,u,einverse,ni,nj,nk,nk,nk,0,2,.true.)
      call Sproject (work3,v,einverse,ni,nj,nk,nk,nk,0,2,.true.)
      call Bpuvm (work2,work3,work2,work3,fmapd,ni,nj,nk)
      call Bdiv (work4,work2,work3,fmapx2,dx,ni,nj,nk)
      call Bvort (work5,work2,work3,fmapx2,dx,ni,nj,nk)
c
c  compute vertical mode amplitudes of pseudo-geopotential
      call Sproject (work2,ps,gmatrixp,ni,nj,nk,1,nk,1,1,.true.)
      call Sproject (work2,t,gmatrixt,ni,nj,nk,nk,nk,1,1,.false.)
      call Nmultv (work2,nijk,1./tpmeanb(nk+1))
c
c  compute vertical modes of avort tendency
      pi4=atan(1.)
      f0=16.*pi4*sin(clat*pi4/45.)/(60.*60.*24)
      print *,' '
      print *,' f0,clat=',f0,clat
      print *,' STATS AFTER FNNMI '
      print *,' '
      call Iagvort  (work1,work5,work2,work3,fmapx2,dx,f0,ni,nj,nk)
c
c  print diagnostics of tendencies of pvort, divg, avort in mode space
      if (lprint) then
         call Ipotvort (work3,work5,work2,edepth,f0,ni,nj,nk)
         call Idiagns (work3,work4,work1,edepth,work2,ni,nj,nk,2)
      endif
cqqqqqqqqqq end
cqqqqqqqqqq end
       call Fmoddiag (pu,pv,pt,pq,ps,u,v,t,q,pd,puten,pvten,ptten,chi,
     &               psten,ux,vx,td,tv,ptv,alogp,work9,sdot,sdotd,
     &               puf,pvf,pufx,pvfx,cpmr,omega,phi,fmapx,fmapd,
     &               fmapx2,fmapd2,topog,cf,sigf,sigh,dsig,grav,R,Cp,
     &               ep1,ptop,div,pum,phi,pvm,ptm,work1,work2,
     &               work3,chi,rinverse,ematrixi,rmatrixp,rmatrixt,
     &             gmatrixr,redepth,Rpweight,work1,rtpmeanb(nkp1),clat,
     &               ecpm,dx,tpratio,nnmis,ni,nj,nk,nkp1,lmoist,lprint,
     &               nvarsi,cfldsi,charsw,workr,iunitir,cfldsrw,ifldsrw,
     &               lpshell,ldaccess,nnvars,nworkr,fields3,fields2,
     &               ndimfld3,ndimfld2,puten,pvten,pum,pvm,ptm,
     &               psm,pub,pvb,ptb,pqb,psb,pdb,sdotb,
     &               sdotdb,omegab,pufb,pufxb,pvfb,pvfxb,ub,vb,uxb,
     &               vxb,tb,cpmrb,qb,tvb,phib,ptvb,tdb,alogpb,divb,
     &               nbsized,pqmin,ir1i,work8,cfldsr,ifldsr,.true.,
     &                  Rtpmeanb,xfull,nfull,itrc,jtrc,ktrc,normtype)


c
c TIME SERIES BLOCK
      if (ltseries) then
         call TENORMZ (pu,pv,pt,ps,u,v,t,pd,ni,nj,nk,e,
     &                    pub,pvb,ptb,psb,ub,vb,tb,pdb,Cp,R,dsig,
     &                    pqb,qb,pq,q,eq,esum,lmoist)
         call Xnorms2a (u,v,t,ps,q,tg,xfull,lmoist,rnormval,
     &               sigf,sigh,dsig,grav,R,Cp,clat,clon,dx,
     &               Redepth,Rinverse,Rmatrixp,Rmatrixt,Rpweight,
     &               qinverse,
     &               Rtpmeanb,work1,work2,work3,work4,work5,work6,
     &               cf,fmapx,fmapd,ni,nj,nk,nfull,nxstar,normtype)
         print *,'EDRY,EQ,ESUM,Rnorm= ',e,eq,esum,rnormval
         call Emodesub (xfull,u,v,t,ps,work3
     &                    ,work1,work2,Redepth,Rinverse,Rmatrixp
     &                    ,Rmatrixt,Rtpmeanb(nk+1),grav,ni,nj,nk
     &                    ,itrc,jtrc,ktrc)
         savenorm(4,1)=e
         savenorm(4,2)=eq
         savenorm(4,3)=esum
         savenorm(4,4)=rnormval
      endif
c
         if (linitwr) then ! write out last iterate of perturbation
c
c  write out header for initialized (perturbation) file
            call Ncopyi (nvarsiw,nvarsw,nnvars)
            call Ncopyv (work2,rvarsw,nrvars)
            if (lpshell) call Nassign (lprint,iunitiw,1,'I')
            call Nopen (nvarsiw(2),iunitiw,lprint)
            iftp=nvarsiw(11)          ! pointer to nftimes array
            nvarsiw(iftp+1)=1         ! number of tslots on tape
            nvarsiw(iftp+3)=1         ! number of tslots on tape
            work2(1,1,1)=60.          ! regular output time = 1 hr
            call Nwriteh (nvarsiw,lvarsw,work2,cfldsw,charsw,workr,nij,
     &                    irwi,iunitiw,fieldh,cfldsrw,ifldsrw,lprint,
     &                    nnvars,nrvars,nlvars,nworkr)
            xtime=0.
            ktauouti=0
c
c  write out fields for initialized (perturbation) file
            call Nwritet (fields3,fields2,workr,xtime,nvarsiw,cfldsrw,
     &                    cfldsw,ifldsrw,ifldsw,nij,nk,ktauouti,irwi,
     &                    iunitiw,lprint,nnvars,nworkr)
c
c  dispose initialized file
            call Nclose (iunitiw,ktauouti,lprint)
            itapeni=1
            if (lpshell) call Ndispose (charsw(6),'I',itapeni,lcopy,
     &                         keepdays,nvarsiw(2),cusercom,lprint)
         endif  ! end block on initialization output
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
c 
c                   BEGIN TIME INTEGRATION
c
c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
      ttime0=time0+ftimed0/60.
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
         print 444,nstep,pu(i,j,nk),pq(i,j,nk),pt(i,j,nk),ps(i,j)
         print 444,nstep,pub(i,j,nk),pqb(i,j,nk),ptb(i,j,nk),psb(i,j)
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
c TIME SERIES BLOCK
      if (ltseries) then
        if (nstep.eq.1) then
          iutsrs=68
          open (unit=iutsrs,file='Ooutput1',form='UNFORMATTED')
          n=0
          do j=1,5
            j1=5+10*j
            do i=1,5
              i1=0+8*i
              n=n+1
              ilocts(n)=i1
              jlocts(n)=j1
              klocts(n)=0
            enddo ! loop over i
          enddo ! loop over j
          write (iutsrs) dt,npointts,nfieldts,ktimax,nk
          write (iutsrs) cfieldts
          write (iutsrs) cfunitts
          write (iutsrs) ilocts,jlocts,klocts
          write (iutsrs) (dlon(ilocts(n),jlocts(n)),n=1,npointts)
     &,                  (dlat(ilocts(n),jlocts(n)),n=1,npointts)
     &,                  sigh,sigf
        endif
        k=10
        write (iutsrs) k,(ps(ilocts(n),jlocts(n)),n=1,npointts)
        k=5
        write (iutsrs) k,(omega(ilocts(n),jlocts(n),k),n=1,npointts)
        k=1
        write (iutsrs) k,(t(ilocts(n),jlocts(n),k),n=1,npointts)
        k=5
        write (iutsrs) k,(t(ilocts(n),jlocts(n),k),n=1,npointts)
        k=8
        write (iutsrs) k,(t(ilocts(n),jlocts(n),k),n=1,npointts)
        k=8
        write (iutsrs) k,(q(ilocts(n),jlocts(n),k),n=1,npointts)
        k=10
        write (iutsrs) k,(tg(ilocts(n),jlocts(n)),n=1,npointts)
      endif
c  End time series block
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
c
c  Write data if time selected
      if (((toutreg.gt.0).and.(mod(nstep,ktimout).eq.0)).or.
     &    ((toutreg.lt.0).and.(nstep.eq.nvarsw(ntoutpt+ktauout)))) 
     &   call Nwritet (fields3,fields2,workr,xtime,nvarsw,cfldsrw,
     &                 cfldsw,ifldsrw,ifldsw,nij,nk,ktauout,irw,
     &                 iunitw,lprint,nnvars,nworkr)
c
c  correct for double counting of precip due to centered time scheme
      if (lmoist) then
         call Ndiffv (rainncv,rainncv,rainncv2,nij)
         call Ndiffv (raincv, raincv, raincv2, nij)
      endif
c
c  BLOCK FOR TIME SERIES
      if (ltseries) then
         k=0
         write (iutsrs) k,(raincv(ilocts(n),jlocts(n)),n=1,npointts)
         write (iutsrs) k,(rainncv(ilocts(n),jlocts(n)),n=1,npointts)
         call TENORMZ (pu,pv,pt,ps,u,v,t,pd,ni,nj,nk,e,
     &                 pub,pvb,ptb,psb,ub,vb,tb,pdb,Cp,R,dsig,
     &                 pqb,qb,pq,q,eq,esum,lmoist)
         call Xnorms2a (u,v,t,ps,q,tg,xfull,lmoist,rnormval,
     &               sigf,sigh,dsig,grav,R,Cp,clat,clon,dx,
     &               Redepth,Rinverse,Rmatrixp,Rmatrixt,Rpweight,
     &               qinverse,
     &               Rtpmeanb,work1,work2,work3,work4,work5,work6,
     &               cf,fmapx,fmapd,ni,nj,nk,nfull,nxstar,normtype)
         print *,'EDRY,EQ,ESUM,Rnorm= ',e,eq,esum,rnormval
         if (nstep+4.le.maxstepe) then
            nesaves=nstep+4
            savenorm(nstep+4,1)=e
            savenorm(nstep+4,2)=eq
            savenorm(nstep+4,3)=esum
            savenorm(nstep+4,4)=rnormval
         endif
      endif
c
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
    1 continue
c
c             END OF ONE TIME STEP INTEGRATION
c 
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c 
      print *, 'end of integration forecast time,step = ',xtime,nstep
c
cqqqqqqqqq
c      write (13) pu
c      write (13) pv
c      write (13) pt
c      write (13) ps
         call Emodesub (xfull,u,v,t,ps,work3
     &                    ,work1,work2,Redepth,Rinverse,Rmatrixp
     &                    ,Rmatrixt,tpmeanb(nk+1),grav,ni,nj,nk
     &                    ,itrc,jtrc,ktrc)
cqqqqqqqqq
C COMPUTE AND PRINT PVORT, DIVG, AGEO STATS FOR COMPUTING TIME SCALE
c BASED ON MODE TENDENCY TO MEASURE BALANCE
c
c  compute vertical mode amplitudes of divergence and vorticity
      call PDECUPL (pu,pv,pt,pq,ps,u,v,t,q,
     &                    pub,pvb,ptb,pqb,psb,ub,vb,tb,qb,
     &                    pd,pdb,lmoist,ni,nj,nk,1)
      call Sproject (work2,u,einverse,ni,nj,nk,nk,nk,0,2,.true.)
      call Sproject (work3,v,einverse,ni,nj,nk,nk,nk,0,2,.true.)
      call Bpuvm (work2,work3,work2,work3,fmapd,ni,nj,nk)
      call Bdiv (work4,work2,work3,fmapx2,dx,ni,nj,nk)
      call Bvort (work5,work2,work3,fmapx2,dx,ni,nj,nk)
c
c  compute vertical mode amplitudes of pseudo-geopotential
      call Sproject (work2,ps,gmatrixp,ni,nj,nk,1,nk,1,1,.true.)
      call Sproject (work2,t,gmatrixt,ni,nj,nk,nk,nk,1,1,.false.)
      call Nmultv (work2,nijk,1./tpmeanb(nk+1))
c
c  compute vertical modes of avort tendency
      pi4=atan(1.)
      f0=16.*pi4*sin(clat*pi4/45.)/(60.*60.*24)
      print *,' '
      print *,' f0,clat=',f0,clat
      print *,' STATS BEFORE FNNMI '
      print *,' '
      call Iagvort  (work1,work5,work2,work3,fmapx2,dx,f0,ni,nj,nk)
c 
c
c  print diagnostics of tendencies of pvort, divg, avort in mode space
      if (lprint) then
         call Ipotvort (work3,work5,work2,edepth,f0,ni,nj,nk)
         call Idiagns (work3,work4,work1,edepth,work2,ni,nj,nk,2)
      endif
cqqqqqqqqqq end
c
cqqqqqqqqqq end
       call Fmoddiag (pu,pv,pt,pq,ps,u,v,t,q,pd,puten,pvten,ptten,chi,
     &               psten,ux,vx,td,tv,ptv,alogp,work9,sdot,sdotd,
     &               puf,pvf,pufx,pvfx,cpmr,omega,phi,fmapx,fmapd,
     &               fmapx2,fmapd2,topog,cf,sigf,sigh,dsig,grav,R,Cp,
     &               ep1,ptop,div,pum,phi,pvm,ptm,work1,work2,
     &               work3,chi,rinverse,ematrixi,rmatrixp,rmatrixt,
     &             gmatrixr,redepth,Rpweight,work1,rtpmeanb(nkp1),clat,
     &               ecpm,dx,tpratio,nnmis,ni,nj,nk,nkp1,lmoist,lprint,
     &               nvarsi,cfldsi,charsw,workr,iunitir,cfldsrw,ifldsrw,
     &               lpshell,ldaccess,nnvars,nworkr,fields3,fields2,
     &               ndimfld3,ndimfld2,puten,pvten,pum,pvm,ptm,
     &               psm,pub,pvb,ptb,pqb,psb,pdb,sdotb,
     &               sdotdb,omegab,pufb,pufxb,pvfb,pvfxb,ub,vb,uxb,
     &               vxb,tb,cpmrb,qb,tvb,phib,ptvb,tdb,alogpb,divb,
     &               nbsized,pqmin,ir1i,work8,cfldsr,ifldsr,.true.,
     &                  Rtpmeanb,xfull,nfull,itrc,jtrc,ktrc,normtype)
         nnmis(4)=0
         call Fnnmi (pu,pv,pt,pq,ps,u,v,t,q,pd,puten,pvten,ptten,chi,
     &               psten,ux,vx,td,tv,ptv,alogp,work9,sdot,sdotd,
     &               puf,pvf,pufx,pvfx,cpmr,omega,phi,fmapx,fmapd,
     &               fmapx2,fmapd2,topog,cf,sigf,sigh,dsig,grav,R,Cp,
     &               ep1,ptop,div,pum,phi,pvm,ptm,work1,work2,
     &               work3,chi,einverse,ematrixi,gmatrixp,gmatrixt,
     &               gmatrixr,edepth,pweight,work1,tpmeanb(nkp1),clat,
     &  ecpm,dx,tpratio,nnmis,ni,nj,nk,nkp1,lmoist,lprint,
     & nvarsi,cfldsi,charsw,workr,iunitir,cfldsrw,ifldsrw,
     &               lpshell,ldaccess,nnvars,nworkr,fields3,fields2,
     &               ndimfld3,ndimfld2,puten,pvten,pum,pvm,ptm,
     &               psm,pub,pvb,ptb,pqb,psb,pdb,sdotb,
     &               sdotdb,omegab,pufb,pufxb,pvfb,pvfxb,ub,vb,uxb,
     &               vxb,tb,cpmrb,qb,tvb,phib,ptvb,tdb,alogpb,divb,
     &               nbsized,pqmin,ir1i,work8,cfldsr,ifldsr,.true.)







c
      if (lcheck) then
         nic=ni/2
         njc=nj/2
         print *,' ps(',nic,',',njc,') = ',ps(nic,njc)
      endif
c
c  close and dispose output file
      print *,' '
      call Nclose (iunitw,ktauout,lprint)
      if (lpshell) call Ndispose (charsw(6),'D',itapnumd,lcopy,
     &                            keepdays,nvarsw(2),cusercom,lprint)
c
c
c  BLOCK FOR TIME SERIES
      if (ltseries) then
         print *,' '
         print *,' Time Series Output'
         if (nesaves.gt.0) then
            write (iutsrs) numnorms
            write (iutsrs) (cnorms(nnorm),nnorm=1,numnorms)
            do nnorm=1,numnorms
               write (iutsrs) nesaves,(savenorm(n,nnorm),n=1,nesaves)
            enddo
         endif 
         call Nclose (iutsrs,ktimax,lprint)
         nch = LNGTH(charsw(6))
         charswts=charsw(6)(1:nch)//'TSERIES'
         itaptser=1
         inumtser=nfieldts*npointts
         if (lpshell) call Ndispose (charswts,'O',itaptser,lcopy,
     &                           keepdays,inumtser,cusercom,lprint)
      endif
c
      stop 'MAIN_OK'
      end
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE FinitGG (pu,pv,pt,ps,pub,pvb,ptb,psb,u,v,t,pd,q,qb,
     &                  ub,vb,tb,pdb,tg,pq,field3c,work,pert,
     &                  work6,ipert,iuniti,ni,nj,nk,lmoist,pertx)
c
c  Read requested optimal mode and scale
c
      logical lmoist
      real pu(ni,nj,nk),   pv(ni,nj,nk),  pt(ni,nj,nk),  ps(ni,nj),
     &     pub(ni,nj,nk), pvb(ni,nj,nk), ptb(ni,nj,nk), psb(ni,nj),
     &     u(ni,nj,nk), v(ni,nj,nk), t(ni,nj,nk), pd(ni,nj),
     &     q (ni,nj,nk), qb (ni,nj,nk) ,
     &     ub(ni,nj,nk), vb(ni,nj,nk), tb(ni,nj,nk), pdb(ni,nj),
     &     tg(ni,nj), pq(ni,nj,nk), field3c(*), work(nk), work6(ni,nj)
c
      print *,'FINIT ni,nj,nk=',ni,nj,nk
      ni2=ni-2
      nj2=nj-2
      ni3=ni-3
      nj3=nj-3
      nfull=2*nk*ni2*nj2+(nk+1)*ni3*nj3
      nfull=nfull+nk*ni3*nj3
      nfullp2=nfull+2
      nij=ni*nj
      nijk=nij*nk
c
      call Nrsetc (pu,nijk,0.)
      call Nrsetc (pv,nijk,0.)
      call Nrsetc (pt,nijk,0.)
      call Nrsetc (pq,nijk,0.)
      call Nrsetc (ps,nij,0.)
c
      call Nopen (nfullp2,iuniti,.t.)
c
c  read in scaled values of u,v,t,ps
      read (iuniti,rec=ipert) E1,B1,(field3c(k),k=1,nfull)
c
c  copy field3c into u,v,t,ps 
      work=1.0  
      call MSTAR1 (u,v,t,ps,q,field3c,work,
     &                   1.,1.,1.,1.,1.,ni,nj,nk,nfull,
     &             0,0,0,0,0)
      work6=ps  ! value of ps prior to scaling
c
c couple perturbation
      call Bx2d (pd,ps,ni,nj,1)  
      call Fcouple (pu,u,pd,ub,pdb,ni,nj,nk,0,2) 
      call Fcouple (pv,v,pd,vb,pdb,ni,nj,nk,0,2) 
      call Fcouple (pt,t,ps,tb,psb,ni,nj,nk,1,2) 
      if (lmoist) call Fcouple (pq,q,ps,qb,psb,ni,nj,nk,1,2) 
c
c
	i = 20
	j = 20
 	k = 5
	print *, 'uvtpq',u(i,j,k),t(i,j,k),ps(i,j),q(i,j,k)	
	print *, 'uvtpqb',ub(i,j,k),tb(i,j,k),psb(i,j),qb(i,j,k)
c
      tmax=0.
      qmax=0.
      vmax=0.
      do 4 k=1,nk
      do 2 j=2,nj-1
      do 2 i=2,ni-1
      if (abs(u(i,j,k)).gt.vmax) vmax=abs(u(i,j,k))
      if (abs(v(i,j,k)).gt.vmax) vmax=abs(v(i,j,k))
    2 continue
      do 3 j=2,nj-2
      do 3 i=2,ni-2
      if (abs(t(i,j,k)).gt.tmax) tmax=abs(t(i,j,k))
      if (lmoist.and.(abs(q(i,j,k)).gt.qmax)) qmax=abs(q(i,j,k))
    3 continue
    4 continue
      pmax=0.
      do 5 j=2,nj-2
      do 5 i=2,ni-2
      if (abs(ps(i,j)).gt.pmax) pmax=abs(ps(i,j))
    5 continue
      rpmax=pmax/0.1
      rtmax=tmax/2.
      rqmax=qmax/0.001
      rvmax=vmax/4.
      rmax1=amax1(rpmax,rtmax)
      rmax2=amax1(rmax1,rqmax)
      rmax=amax1(rmax2,rvmax)
c
      call cbdryx (pu,pv,pt,ps,ni,nj,nk)
      if (lmoist) call Cbdry0 (pq,ni,nj,nk,1)
c
c  rescale inititial values by pert
c
      pertx=pert/rmax
      pu=pertx*pu
      pv=pertx*pv
      pt=pertx*pt
      ps=pertx*ps
      if (lmoist) pq=pertx*pq
c
c
	i = 20
	j = 20
 	k = 5
	print *, 'uvtpq',pu(i,j,k),pt(i,j,k),ps(i,j),pq(i,j,k)	
c
c     
      print *,' P**: mode number ',ipert,' read: eigenvalue=',E1,B1
      print *,' P19: initial values scaled by pertx=',pertx,pert
      print *,' '
c
      return
      end
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE MSTAR1 (u,v,t,ps,q,xstar,dsig,
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     &                   R,Cp,Tnorm,Pnorm,Anorm,ni,nj,nk,nxstar,
     &             ISET1,ISET2,ISET3,ISET4,ISET5)
C
      real u(ni,nj,nk), v(ni,nj,nk), t(ni,nj,nk), ps(ni,nj)
     &,    xstar(nxstar), dsig(nk) , q (ni,nj,nk) 
c
c  transform from vector format to model format
c  Determine initial (uncoupled) TLM fields from xstar
c  (include scaling factors Matrix A)
c
      ni1=ni-1
      nj1=nj-1
      ni2=ni-2
      nj2=nj-2
      nxu=0
      nxv=nxu+ni2*nj2*nk
      nxt=nxv+ni2*nj2*nk
      nxp=nxt+(ni2-1)*(nj2-1)*nk
      nxq=nxp+(ni2-1)*(nj2-1)
      call Cbdryx (u,v,t,ps,ni,nj,nk)
      call Cbdry0 (q,ni,nj,nk,1)
      afac=sqrt(2.)*Anorm
c
c	fill q
c
	nx = nxq
        do 4 k = 1 , nk
        do 4 j = 2 , nj2
        do 4 i = 2 , ni2
        nx = nx+1
	q (i,j,k) = xstar (nx)
4       continue
c
c  Fill and unscale ps 
      nx=nxp
      fac=afac*Pnorm/sqrt(R*Tnorm)
      do 1 j=2,nj2
      do 1 i=2,ni2
      nx=nx+1
      if ( iset3.eq.1 ) then
      		ps (i,j) = 0.0
	else
      		ps(i,j)=xstar(nx)
      endif
    1 continue
c  
c  Fill and unscale pt
      nx=nxt
      fac1=afac*sqrt(Tnorm/Cp)
      do 2 k=1,nk
      fac=fac1/sqrt(dsig(k)) 
      do 2 j=2,nj2
      do 2 i=2,ni2
      nx=nx+1
      if ( iset3.eq.1 ) then 
      		t(i,j,k) = 0.0
	else
      		t(i,j,k)=xstar(nx)
      endif
    2 continue
c  
c  Fill and unscale pu and pv
      do 3 k=1,nk
      fac=afac/sqrt(dsig(k)) 
      do 3 j=2,nj1
      do 3 i=2,ni1
      nxu=nxu+1
      nxv=nxv+1
      u(i,j,k)=xstar(nxu)
      v(i,j,k)=xstar(nxv)
    3 continue
c
      return
      end
c
      SUBROUTINE TENORMZ (pu,pv,pt,ps,u,v,t,pd,ni,nj,nk,e,
     &                    pub,pvb,ptb,psb,ub,vb,tb,pdb,Cp,R,dsig,
     &                    pqb,qb,pq,q,eq,esum,lmoist)
      logical lmoist
      real u(ni,nj,nk), v(ni,nj,nk), t(ni,nj,nk), ps(ni,nj), dsig(nk)
      real q(ni,nj,nk)
      call Bx2d (pdb,psb,ni,nj,1)
      call Bdecup (ub,pub,pdb,ni,nj,nk,0)
      call Bdecup (vb,pvb,pdb,ni,nj,nk,0)
      call Bdecup (tb,ptb,psb,ni,nj,nk,1)
      if (lmoist) call Bdecup (qb,pqb,psb,ni,nj,nk,1)
      call Bx2d (pd,ps,ni,nj,1)
      call Fdecup (u,pu,pd,ub,pdb,ni,nj,nk,0)
      call Fdecup (v,pv,pd,vb,pdb,ni,nj,nk,0)
      call Fdecup (t,pt,ps,tb,psb,ni,nj,nk,1)
      if (lmoist) call Fdecup (q,pq,ps,qb,psb,ni,nj,nk,1)
c
c  compute dry norm contribution 
c
      e=0.
      do 1 k=1,nk
      do 1 j=2,nj-1
      do 1 i=2,ni-1 
      e=e+(u(i,j,k)*u(i,j,k)+v(i,j,k)*v(i,j,k))*dsig(k)
    1 continue
      e=e/((ni-2)*(nj-2))  

cqqqq
      eprint=e/2.
      print *,'KE EPRINT=',EPRINT
      eppp=0.
      etttt=0. !qqqqqqqq

      tbar=270.
      pbar=100.
      pfac=R*tbar/((pbar**2)*(ni-3)*(nj-3))
      tfac=Cp/(tbar*(ni-3)*(nj-3))
      do 2 j=2,nj-2
      do 2 i=2,ni-2 
      e=e+ps(i,j)*ps(i,j)*pfac
      eppp=eppp+ps(i,j)*ps(i,j)*pfac  !qqqqqqqqq


      do 2 k=1,nk
      etttt=etttt+t(i,j,k)*t(i,j,k)*tfac*dsig(k) !qqqqqqq
      e=e+t(i,j,k)*t(i,j,k)*tfac*dsig(k)
    2 continue  
      e=e/2.  


      eprint=eppp/2.
      print *,'PP EPRINT=',EPRINT

      eprint=etttt/2. + eprint
      print *,'PT EPRINT=',EPRINT

c
c  compute moist norm contribution, 1/29/96
c
      if (lmoist) then
      eq = 0.0
      qwgt = (2.5104e+06 / sqrt ( tbar * 1005.7 )) 
      do 10 k = 1 , nk
      fac = sqrt ( 0.5 * dsig (k) ) * qwgt 
      fac = fac * fac/((ni-3)*(nj-3)) 
      do 10 j = 2 , nj - 2
      do 10 i = 2 , ni - 2
      eq = eq + fac * q (i,j,k) * q (i,j,k) 
10    continue
      else
      eq=0.
      endif
c
      esum=eq+e
cqqq      print *,' E-norm: e,eq,esum = ',e,eq,esum
c
c
      return
      end
c
c
c
      SUBROUTINE Emodesub (xfull,u,v,t,ps,geop
     &                    ,work1,work2,edepth,einverse,gmatrixp
     &                    ,gmatrixt,pmeanb,grav,ni,nj,nk
     &                    ,itrunc,jtrunc,ktrunc)
c
c  compute energy as function of type and vertical mode
c
      real xfull(itrunc,jtrunc,ktrunc)
     &,    geop(ni,nj,nk)
     &,    work1(ni,nj,nk), work2(ni,nj,nk)
     &,    u(ni,nj,nk), v(ni,nj,nk), t(ni,nj,nk), ps(ni,nj)
     &,    edepth(nk), e(nk+1,5)
     &,    einverse(nk,nk), gmatrixp(nk), gmatrixt(nk,nk)
      character*13 ctype(5)
      data ctype(1)/'        KE   '/
      data ctype(2)/'        APE  '/
      data ctype(3)/'        ROT  '/
      data ctype(4)/'        GRAV '/
      data ctype(5)/'        TOT  '/
c             
      call Sproject (work1,u,einverse,ni,nj,nk,nk,nk,0,2,.t.)
      call Sproject (work2,v,einverse,ni,nj,nk,nk,nk,0,2,.t.)
      call Sproject (geop,ps,gmatrixp,ni,nj,nk,1,nk,1,2,.t.)
      call Sproject (geop,t,gmatrixt,ni,nj ,nk,nk,nk,1,2,.f.) 
      rpmeanb=1./pmeanb
      call Nmultv (geop,ni*nj*nk,rpmeanb)
c
      ndgrid=(ni-2)*(nj-2)
      nxgrid=(ni-3)*(nj-3)
c
c  compute KE      
      fac=0.5/ndgrid  
c
      do k=1,nk
        e(k,1)=0.
        do j=2,nj-1  
          do i=2,ni-1
            e(k,1)=e(k,1)+fac*(
     &         work1(i,j,k)*work1(i,j,k)+work2(i,j,k)*work2(i,j,k)) 
          enddo
        enddo
      enddo
c
c  compute APE
      do k=1,nk
        e(k,2)=0.
        fac=0.5/(nxgrid*edepth(k))
        do j=2,nj-2
          do i=2,ni-2
            e(k,2)=e(k,2)+fac*geop(i,j,k)*geop(i,j,k)
          enddo
        enddo
      enddo
c
c  compute total E 
      do k=1,nk
        e(k,3)=0.
        e(k,5)=e(k,1)+e(k,2)
      enddo
c
c  compute ROT E
      do k=1,ktrunc
        do j=1,jtrunc
          do i=1,itrunc
            e(k,3)=e(k,3)+xfull(i,j,k)*xfull(i,j,k)
          enddo
        enddo
      enddo
c
c  compute GRAV E as residual; determine totals 
      do k=1,nk
        e(k,4)=e(k,5)-e(k,3)
      enddo
      do i=1,5
        e(nk+1,i)=0.
        do k=1,nk
          e(nk+1,i)=e(nk+1,i)+e(k,i)
        enddo
      enddo
c
c  print table
c
      print *,' '
      print ('(a,a,5a13)'), 'MODE',' EQ.DEPTH',(ctype(i),i=1,5)
      do k=1,nk
        edx=edepth(k)/grav
        print ('(i4,f9.1,1p5e13.5)'), k,edx,(e(k,i),i=1,5)
      enddo
      print ('(a,9x,1p5e13.5)'),' TOT',(e(nk+1,i),i=1,5) 
      print *,' SUM=',e(nk+1,5)
      print *,' APE=',e(nk+1,2)
      print *,'  KE=',e(nk+1,1)
c
      return
      end

C 
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE Sinitial (ematrix,einverse,gmatrixt,gmatrixp,tmatrix,
     &                     pmatrix,edepth,tpmeanb,msteps,ptten,
     &                     psten,phten,hydros,psfacp,hmatrix,xmatrix,
     &                     ematrixi,gmatrixr,pweight,w,sigf,sigh,dsig,
     &                     R,Cp,ptop,dx,dt,nix,njx,nw,nk,nkp1,
     &                     nsplit,lprint,linitial,lsplit)
c
      logical lprint,linitial,lsplit
      real ematrix(nk,nk), einverse(nk,nk), gmatrixt(nk,nk)
     &,    gmatrixp(nk), tmatrix(nk,nk), pmatrix(nk), edepth(nk)
     &,    tpmeanb(nkp1), ptten(nix,njx,nk)
     &,    psten(nix,njx), phten(nix,njx,nk), hydros(nw), psfacp(nw)
     &,    hmatrix(nk,nk), xmatrix(nk,nk), ematrixi(nk,nk) 
     &,    gmatrixr(nk,nk), pweight(nk+2,3), w(nw,13), sigf(nkp1)
     &,    sigh(nk), dsig(nk)
      integer msteps(nsplit)
c
      common /cctsmatx/ ctsmatx(11), bnbmatx(10) !qqqqq
c
c  determine normal mode projection operators for split-explicit 
c  and implicit normal mode initialization schemes
c
      nk2=nk*nk
      pmeanb=tpmeanb(nkp1)
c
c  compute arrays required to fill matrices
      call Smatrix0 (phten,ptten,psten,hydros,tpmeanb,psfacp,
     &               w(1,1),w(1,2),w(1,3),w(1,4),
     &               w(1,5),w(1,6),w(1,7),w(1,8),w(1,9),
     &               w(1,10),w(1,11),w(1,12),w(1,13),sigf,sigh,dsig,
     &               R,Cp,ptop,dx,grav,nix,njx,nk)
c
c  fill matrices
      call Smfill (hmatrix,phten,nix,njx,nk,nk)
      call Smfill (tmatrix,ptten,nix,njx,nk,nk)
      call Smfill (gmatrixt,hydros,nix,njx,nk,nk)
      call Smfill (pmatrix,psten,nix,njx,nk,1)
      call Smfill (gmatrixp,psfacp,nix,njx,1,nk)
c
c  compute eigenmodes (rg calls eispack routines on the cray)
      numerr=0
      call Ncopyv (xmatrix,hmatrix,nk2)
      call Nmultv (xmatrix,nk2,-1.)
      call Rg (nk,nk,xmatrix,edepth,w,1,ematrixi,w(1,2),w(1,3),ier)
      call Scheki (ier,numerr,'ematrixi','Rg      ')
      call Scheke (edepth,w,nk,numerr,'hmatrix ')
      call Sorder (ematrixi,edepth,w(1,1),w(1,2),nk)
      call Snorml (ematrixi,dsig,nk)
c
c  compute inverse of ematrixi
      call Invmtx (ematrixi,nk,einverse,nk,nk,det,w,ier)
      call Scheki (ier,numerr,'einverse','Invmtx  ')
c
cqqqqqqqqq begin
c
         call Smfill (w(1,2),w(1,11),nix,njx,nk,nk)      !  S / pmeanb
         call Nmultv (w(1,2),nk2,pmeanb)                 !  S
         call Invmtx (w(1,2),nk,xmatrix,nk,nk,det,w,ier) ! S**-1
         call Smultm (ctsmatx,dsig,xmatrix,1,nk,nk)    !-c^T/psmeanb * S**-1   
         pmeanbng=-pmeanb
         call Nmultv (ctsmatx,nk,pmeanbng)             ! c^T * S**-1
         call Ncopyv (xmatrix,gmatrixt,nk2)
         call Invmtx (xmatrix,nk,w(1,2),nk,nk,det,w,ier) ! B**-1
         call Ncopyv (xmatrix,gmatrixp,nk)
         call Sgmatrxp (xmatrix,tpmeanb,pmeanb,sigh,ptop,R,nk)
         call Smultm (bnbmatx,w(1,2),xmatrix,nk,nk,1)    ! B**-1 * b    
         call Smultm (ctsmatx(nk+1),ctsmatx,bnbmatx,1,nk,1) ! c^T S- B- b 
         call Sprntv (ctsmatx,nk+1,'ctsmatx ') 
         call Sprntv (bnbmatx,nk,'bnbmatx ') 
         call Sprntv (xmatrix,nk,'xmatrix ') 
         call Smultm (w(1,2),gmatrixt,bnbmatx,nk,nk,1)    
         call Sprntv (w(1,2),nk,'xmatrix2') 
c
cqqqqqqqqq end
c
c  print some information
c
cqqq      if (lprint.or.(numerr.ne.0)) then
         print *,' '
         print *,' Vertical modes computed'
         if (numerr.ne.0) print *,' WARNING: ',numerr,' ERRORS detected'
         print *,' ' 
         call Sprntv (tpmeanb,nk,'tmeanb  ')
         call Sprntv (pmeanb,1,'pmeanb  ')
         call Sprntv (edepth,nk,'edepth  ')
         call Sprntm (ematrixi,nk,nk,'ematrixi')
         call Sprntm (einverse,nk,nk,'einverse')
         call Sprntm (hmatrix,nk,nk,'hmatrix ')
         call Sprntm (gmatrixt,nk,nk,'gmatrixt')
         call Sprntv (gmatrixp,nk,'gmatrixp') 
         call Sprntv (pmatrix,nk,'pmatrix ') 
c
c  compute and print orthogonality of vertical modes E^T DSIG E    
         call Sdiagmat (w(1,1),dsig,nk)           
         call Atransp (w(1,2),ematrixi,nk,nk)     
         call Smultm (w(1,3),w(1,1),ematrixi,nk,nk,nk)  
         call Smultm (w(1,4),w(1,2),w(1,3),nk,nk,nk)
         call Sprntm (w(1,4),nk,nk,'orthogon')
c
c  compute and print energy matrix for T: -B^T DSIG S^-1
c  S/pmeanb is the matrix multiplying div to get dT/dt, 
c  but div is mass divergence, not wind divergence 
         call Smfill (w(1,2),w(1,11),nix,njx,nk,nk)   !  S / pmeanb 
         pmeanbng=-pmeanb
         call Nmultv (w(1,2),nk2,pmeanbng)            ! -S
         call Invmtx (w(1,2),nk,w(1,3),nk,nk,det,w(1,4),ier) 
         call Atransp (w(1,4),gmatrixt,nk,nk)         ! B^T
         call Smultm (w(1,2),w(1,1),w(1,3),nk,nk,nk)  ! DSIG * S^-1
         call Smultm (w(1,3),w(1,4),w(1,2),nk,nk,nk)  
         call Sprntm (w(1,3),nk,nk,'t energy')       
c
cqqq      endif 
c
c  modify matrices for either scheme
c  after call Sgmatrxp, gmatrixp takes ps to pseudogeop  
c  after call Nmultv, gmatrixp takes ps to modes of psmean*pseudogeop  
      call Sgmatrxp (gmatrixp,tpmeanb,pmeanb,sigh,ptop,R,nk)
      pweight(nk+2,1)=gmatrixp(nk)
      call Ncopyv (xmatrix,gmatrixp,nk)
      call Smultm (gmatrixp,einverse,xmatrix,nk,nk,1)
      call Nmultv (gmatrixp,nk,pmeanb)
c 
c  after call Nmultv, gmatrixt takes t to modes of psmean*pseudogeop  
      pweight(nk+1,1)=gmatrixt(nk,nk)
      call Ncopyv (xmatrix,gmatrixt,nk2)
      call Smultm (gmatrixt,einverse,xmatrix,nk,nk,nk)
      call Nmultv (gmatrixt,nk2,pmeanb)
c
c  modify matrices for split explicit scheme
      if (lsplit) then
         call Ncopyv (xmatrix,tmatrix,nk2)
         call Smultm (tmatrix,xmatrix,ematrixi,nk,nk,nk)
         call Ncopyv (xmatrix,pmatrix,nk)
         call Smultm (pmatrix,xmatrix,ematrixi,1,nk,nk)
         do 1 n=1,nsplit
         rxm=2.*dt/(2.*msteps(n) + 1.)
         pmatrix(n)=rxm*pmatrix(n)
         do 1 k=1,nk
         tmatrix(k,n)=rxm*tmatrix(k,n)
         ematrix(k,n)=rxm*ematrixi(k,n)
    1    continue
      endif
c
c  determine additional matrices for nnmi initialization scheme
      if (linitial) then
c
c  determine weight vector used to determine changes in ps
c  using Machenhauer's scheme.  hmatrix over-written 
c  pweight takes modes of change of h into change of ps
c  pweight= -(pmeanb*dsig)^T H**-1 E
         call Invmtx (hmatrix,nk,xmatrix,nk,nk,det,w,ier)
         call Scheki (ier,numerr,'hmatrix ','Invmtx  ')
         call Smultm (hmatrix,xmatrix,ematrixi,nk,nk,nk)
         do 2 k=1,nk
         pweight(k,1)=0.
         do 2 n=1,nk
         pweight(k,1)=pweight(k,1)-pmeanb*dsig(n)*hmatrix(n,k)
    2    continue
c
c  determine gmatrixr to take modes of phi to t in physical space
         call Ncopyv (xmatrix,gmatrixt,nk2)
         call Invmtx (xmatrix,nk,gmatrixr,nk,nk,det,w,ier)
         call Scheki (ier,numerr,'gmatrixr','Invmtx  ')
         call Nmultv (gmatrixr,nk2,pmeanb)
c
c  determine factors for minimizing vertical gradients of t
         pwsum=0.
         do 3 k=1,nk-1
         pweight(k,3)=1./(sigh(k+1)-sigh(k))
         pwsum=pwsum+pweight(k,3)
    3    continue  
         pweight(nk+1,2)=0.
         call Smultm (pweight(1,2),gmatrixr,gmatrixp,nk,nk,1)
         do 4 k=1,nk-1
         pweight(k,3)=pweight(k,3)/pwsum
         dc=(pweight(k+1,2)-pweight(k,2))/pmeanb
         pweight(k,2)=dc*pweight(k,3)
         pweight(nk+1,2)=pweight(nk+1,2)+dc*pweight(k,2)
    4    continue
c
      endif
c
      return
      end
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE Idiagn1 (f,d,ni,nj,nk,ng,n1)
      real f(ni,nj,nk), d(nk)
c
c  This routine computes the maximum absolute value of one field f.
c
      ni1=ni-ng-n1
      nj1=nj-ng-n1
c
      do 2 k=1,nk
      fs=0.
      nn=0
      do 1 j=n1,nj1
      do 1 i=n1,ni1
      fs=fs+f(i,j,k)*f(i,j,k)
      nn=nn+1
    1 continue
      d(k)=sqrt(fs/nn)
    2 continue
c
      return
      end
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE Fmoddiag 
     &                 (pu,pv,pt,pq,ps,u,v,t,q,pd,puten,pvten,ptten,chi,
     &                  psten,ux,vx,td,tv,ptv,alogp,work9,sdot,sdotd,
     &                  puf,pvf,pufx,pvfx,cpmr,omega,phi,fmapx,fmapd,
     &                  fmapx2,fmapd2,topog,cf,sigf,sigh,dsig,grav,R,Cp,
     &                  ep1,ptop,divg,vort,geop,avort,pvort,work1,work2,
     &                  work3,delp,einverse,ematrixi,gmatrixp,gmatrixt,
     &                  gmatrixr,edepth,pweight,dum1,pmeanb,clat,ecpm,
     &                  dx,tpratio,nnmis,ni,nj,nk,nkp1,lmoist,lprint,
     &                  nvars,cfldsi,charsw,workr,iunit,cfldsrw,ifldsrw,
     &                  lpshell,ldaccess,nnvars,nworkr,fields3,fields2,
     &                  ndimfld3,ndimfld2,dum,dvm,putenb,pvtenb,pttenb,
     &                  pstenb,pub,pvb,ptb,pqb,psb,pdb,sdotb,
     &                  sdotdb,omegab,pufb,pufxb,pvfb,pvfxb,ub,vb,uxb,
     &                  vxb,tb,cpmrb,qb,tvb,phib,ptvb,tdb,alogpb,divb,
     &                  nbsized,pqmin,irb1,work8,cfldsr,iflds,lfirst,
     &                  tpmeanb,xfull,nfull,itrc,jtrc,ktrc,normtype)
c
      real tpmeanb(*), xfull(nfull)
c
c  perform an implicit nonlinear normal mode initialization
c
      logical lmoist,lprint,lpshell,ldaccess,lfirst
      real pu(ni,nj,nk), pv(ni,nj,nk), pt(ni,nj,nk), pq(ni,nj,nk)
     &,    u(ni,nj,nk), v(ni,nj,nk), t(ni,nj,nk), q(ni,nj,nk)
     &,    puten(ni,nj,nk), pvten(ni,nj,nk), ptten(ni,nj,nk)
     &,    ux(ni,nj,nk), vx(ni,nj,nk), td(ni,nj,nk), tv(ni,nj,nk)
     &,    ps(ni,nj), pd(ni,nj), psten(ni,nj), ptv(ni,nj,nk)
     &,    sdot(ni,nj,nkp1), alogp(ni,nj,nkp1)
     &,    sdotd(ni,nj,nkp1), puf(ni,nj,nk), pvf(ni,nj,nk)
     &,    pufx(ni,nj,nk), pvfx(ni,nj,nk), cpmr(ni,nj,nk)
     &,    omega(ni,nj,nk), phi(ni,nj,nk), chi(ni,nj), dum1(ni,nj,nk)
      real pub(ni,nj,nk), pvb(ni,nj,nk), ptb(ni,nj,nk), pqb(ni,nj,nk)
     &,    ub(ni,nj,nk), vb(ni,nj,nk), tb(ni,nj,nk), qb(ni,nj,nk)
     &,    putenb(ni,nj,nk), pvtenb(ni,nj,nk), pttenb(ni,nj,nk)
     &,    uxb(ni,nj,nk), vxb(ni,nj,nk), tvb(ni,nj,nk)
     &,    psb(ni,nj), pdb(ni,nj), pstenb(ni,nj), ptvb(ni,nj,nk)
     &,    sdotb(ni,nj,nkp1)
     &,    pufb(ni,nj,nk), pvfb(ni,nj,nk), pufxb(ni,nj,nk)
     &,    pvfxb(ni,nj,nk), cpmrb(ni,nj,nk), omegab(ni,nj,nk)
     &,    phib(ni,nj,nk),  alogpb(ni,nj,nkp1), divb(ni,nj,nk)
     &,    sdotdb(ni,nj,nkp1), tdb(ni,nj,nk), work8(ni,nj,13)
      real fmapx(ni,nj), fmapd(ni,nj), fmapx2(ni,nj), fmapd2(ni,nj)
     &,    topog(ni,nj), cf(ni,nj), sigf(nkp1), sigh(nk), dsig(nk)
      real workr(nworkr), fields3(ndimfld3), fields2(ndimfld2)
      real divg(ni,nj,nk), vort(ni,nj,nk), geop(ni,nj,nk)
     &,    avort(ni,nj,nk), pvort(ni,nj,nk), delp(ni,nj), edepth(nk)
     &,    work1(ni,nj,nk), work2(ni,nj,nk+1) ,work3(ni,nj,nk)
     &,    einverse(nk,nk), ematrixi(nk,nk), gmatrixp(nk)
     &,    gmatrixt(nk,nk), gmatrixr(nk,nk), pweight(nk+2,3)
     &,    dum(ni,nj,nk), dvm(ni,nj,nk), work9(ni,nj,nk+2,2)
      integer nnmis(4), nvars(nnvars), ifldsrw(3), iflds(5)
      character*(*) cfldsi(*), cfldsrw(*), cfldsr(*)
      character*(*) charsw(*)
c
c  nnmis(1): number of modes to be initialized
c  nnmis(2): number of nonlinear iterations
c  nnmis(3): number of iterations of helmholtz solver 
c  nnmis(4): number of time slots on nnmi-type file
c
c  The following are derived fields in vertical mode space (they may be
c  the fields themselves, their tendencies, or adjustments):
c    divg  = horizontal velocity divergence
c    vort  = vertical component of vorticity
c    geop  = perturbation pseudo-geopotential
c    avort = ageostrophic vort = vort - (1/f0) del**2 geop
c    pvort = linearized potential vorticity = vort - (f0/gH) geop
c
      nij=ni*nj
      nijk=nij*nk
      nk1=nk+1
      irb=irb1+1
c
      if (lprint) print *,' Diagnostic on mode tendencies follows'
c
c  set some constants
      nmodes=nnmis(1)
      pi4=atan(1.)
      f0=16.*pi4*sin(clat*pi4/45.)/(60.*60.*24)
      rpmeanb=1./pmeanb
      if (lprint) then
         nkmd=nk     ! compute results for all vertical modes 
      else
         nkmd=nmodes ! only compute results for vertical modes needed
      endif
      nijkmd=nij*nkmd
c
c  compute adiabatic tendencies of basic state (used for decoupling)
      call Nrsetc (putenb,nijk,0.)
      call Nrsetc (pvtenb,nijk,0.)
      call Nrsetc (pttenb,nijk,0.)
      call Nrsetc (pstenb,nij,0.)
      call Ntendc (pub,pvb,ptb,psb,ub,vb,tb,pdb,alogpb,work9,sdotb,
     &             omegab,phib,chi,divb,uxb,vxb,tvb,ptvb,pqb,
     &             qb,cpmrb,pufb,pvfb,sigf,sigh,dsig,fmapx,fmapd,
     &             fmapx2,fmapd2,topog,grav,R,Cp,putenb,pvtenb,
     &             pttenb,pufxb,pvfxb,sdotdb,tdb,cf,dum1,ep1,
     &             ecpm,ptop,dx,ni,nj,nk,nkp1,pstenb,lmoist,.false.)
c
c  compute adiabatic tendencies of perturbations
      call Nrsetc (puten,nijk,0.)
      call Nrsetc (pvten,nijk,0.)
      call Nrsetc (ptten,nijk,0.)
      call Nrsetc (psten,nij,0.)
      call Ftendc (psten,puten,pvten,ptten,sdot,pd,pu,pv,pt,ps,u,
     &             v,t,puf,pufx,pvf,pvfx,ux,vx,divg,sdotd,omega,
     &             alogp,pub,pvb,ptb,psb,alogpb,sdotb,divb,
     &             omegab,pufb,pufxb,pvfb,pvfxb,ub,vb,uxb,vxb,tb,
     &             sdotdb,pdb,chi,fmapx,fmapd,fmapx2,fmapd2,cf,
     &             sigf,sigh,dsig,ptv,ptvb,cpmr,cpmrb,q,pq,
     &             qb,lmoist,ep1,ecpm,R,Cp,ptop,dx,grav,ni,nj,nk,
     &             pqb,dum1,nkp1,td,tdb,tv,tvb,phi,phib,work8,
     &             work9,.false.)
c
c  Decouple basic state tendencies.   work1 = pstenb on dot grid.
c  Overwite coupled tendencies by corresponding uncoupled ones
      call Bx2d (work1,pstenb,ni,nj,1)
      call Fdecup (putenb,putenb,work1,ub,pdb,ni,nj,nk,0)
      call Fdecup (pvtenb,pvtenb,work1,vb,pdb,ni,nj,nk,0)
      call Fdecup (pttenb,pttenb,pstenb,tb,psb,ni,nj,nk,1)
c
c  Decouple perturbation tendencies.  work2 = psten on dot grid.
c  Overwite coupled tendencies by corresponding uncoupled ones
      call Bx2d (work2,psten,ni,nj,1)
      call Fdecupt (puten,puten,work2,u,pd,putenb,work1,ub,pdb,
     &                    ni,nj,nk,0)
      call Fdecupt (pvten,pvten,work2,v,pd,pvtenb,work1,vb,pdb,
     &                    ni,nj,nk,0)
      call Fdecupt (ptten,ptten,psten,t,ps,pttenb,pstenb,tb,psb,
     &                    ni,nj,nk,1)
c
c  compute perturbed tendency of virtual temperature (assume dq/dt=0)
      call Btmoist (work1,pttenb,qb,ep1,ni,nj,nk,lmoist)
      call Ftmoist (ptten,ptten,q,work1,pttenb,ep1,ni,nj,nk,lmoist,2)
c
c compute and print diagnostic table:
c   ub, vb, tb are work arrays
c 
         call Xnorms2a (u,v,t,ps,q,tg,xfull,lmoist,rnormval,
     &               sigf,sigh,dsig,grav,R,Cp,clat,clon,dx,
     &               edepth,einverse,gmatrixp,gmatrixt,pweight,
     &               qinverse,
     &               tpmeanb,work1,work2,work3,ub,vb,tb,
     &               cf,fmapx,fmapd,ni,nj,nk,nfull,nfull,normtype)
         print *,' '
         print *,' Energy diagnostics'
         call Emodesub (xfull,u,v,t,ps,work3
     &                    ,work1,work2,edepth,einverse,gmatrixp
     &                    ,gmatrixt,pmeanb,grav,ni,nj,nk
     &                    ,itrc,jtrc,ktrc)
         call Xnorms2a (puten,pvten,ptten,psten,q,tg,xfull,
     &               lmoist,rnormval,
     &               sigf,sigh,dsig,grav,R,Cp,clat,clon,dx,
     &               edepth,einverse,gmatrixp,gmatrixt,pweight,
     &               qinverse,
     &               tpmeanb,work1,work2,work3,ub,vb,tb,
     &               cf,fmapx,fmapd,ni,nj,nk,nfull,nfull,normtype)
         print *,' '
         print *,' Tendency diagnostics'
         call Emodesub (xfull,puten,pvten,ptten,psten,work3
     &                    ,work1,work2,edepth,einverse,gmatrixp
     &                    ,gmatrixt,pmeanb,grav,ni,nj,nk
     &                    ,itrc,jtrc,ktrc)
c
      return
      end


      Subroutine Varnorm (u,v,T,ps,ni,nj,nk,cname)
      real u(ni,nj,nk), v(ni,nj,nk), t(ni,nj,nk), ps(ni,nj)
      character*(*) cname
c
c variances crudley interpolated from analysis error variances estimated by
c R. Daley  (April 1999).
      parameter (nkvars=10)
      real tvars(nkvars), uvars(nkvars), pvars
      data tvars/5.,3.,2.,1.7,1.7,2.,1.6,1.2,5.,9./
      data uvars/100.,30.,25.,20.,16.,12.,8.,4.,4.,8./
      data pvars/0.2/
c
      xnorm=0.
      ni1=ni-1
      nj1=nj-1
      ni2=ni-2
      nj2=nj-2
      anorm2=2.*ni2*nj2
c
      do j=2,nj2
      do i=2,ni2
        xnorm=xnorm+ps(i,j)*ps(i,j)/pvars     
      do k=1,nk
        xnorm=xnorm+t(i,j,k)*t(i,j,k)/tvars(k)        
      enddo 
      enddo 
      enddo 
c      
      do j=2,nj1
      do i=2,ni1
      do k=1,nk
        xnorm=xnorm+u(i,j,k)*u(i,j,k)/uvars(k)        
        xnorm=xnorm+v(i,j,k)*v(i,j,k)/uvars(k)
      enddo 
      enddo 
      enddo 
c      
      xnorm=xnorm/anorm2
c
      print *,' '
      print *,' VARNORM = ',xnorm,'  ',cname
      print *,' '
c
      return
      end


'EOFD'
rm /usr/tmp/rmerrico/RAEDER/MAMS/EXPL/H24-48/LANC/DRYTLM/ISVar/I0-10/V1-2/vdim.dat.1
acquire /RAEDER/MAMS/EXPL/H24-48/LANC/DRYTLM/ISVar/I0-10/V1-2/vdim.dat.1 0 xvector

cp ~/Decks/lib2A3.o libadj.a
cp ~/Decks/xnorms.o libxnorms.a

assign -a xvector fort.10
assign -a xinit fort.13

f90raed lnf.f llu "-evi -r2" adj,xnorms,ncarm,newncaro exec
ls -l

#dispose xinit SV/BAL/EXPLOS/TLM/DRYE/C1/NNMITEST/UVTPS copy 365 tseries

rcp lnf.out $TOSUN/vnormtlmR.P

ja -cst
exit

'JOBEND'
rsh faith.cgd.ucar.edu -l ron beep
rcp log $TOSUN/vnormtlmR.L

exit











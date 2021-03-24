c  This file changes routine PVpv2gp and PVpv2gpA so that rotational modes
c  are weighted according to equivalent depth.  This changes the norm only at
c  the initial time.  
c
c
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
         SUBROUTINE Xnorms0 (ni,nj,nk,itrc,jtrc,ktrc,itrcq,jtrcq,ktrcq
     &                      ,itype)
         real qwgt (10) 
         common /cnorm1/ Tnorm,Pnorm,Anorm
         common /cnorm2/ itrunc,jtrunc,ktrunc
         common /cnorm3/ qwgt
         common /cnorm4/ itruncq,jtruncq,ktruncq,qnorm
c
c  Set some weights or variables required by norms, 7/1/94
c
         Tnorm=270.
         Pnorm=100.
         Anorm=1.
c
         cp   = 1005.7 
         clat = 2.5104e+06
c
c  account for number of grid points (in horizontal), 2/3/97 m.e.
c
         xn1 = float ( (ni-3) * (nj-3) )   !  T, p_s, q
         xn2 = float ( (ni-2) * (nj-2) )   !  u, v
         Tnorm = Tnorm * ( xn1 / xn2 )
         Pnorm = Pnorm * ( xn1 / xn2 )
         Anorm = Anorm * sqrt ( xn2 )
c
         itrunc = itrc
         jtrunc = jtrc
         ktrunc = ktrc
c
c  factor of 8 for def of "energy" and integral(sin**2) in x and y
         qnorm=clat/sqrt(8.*Cp*tnorm)
c
         itruncq = itrcq
         jtruncq = jtrcq
         ktruncq = ktrcq
c
         epsilon = 1.0   !  sqrt (10.0) or 1.0/sqrt(10.0)
         do 10 k = 1 , nk
         qwgt (k) = clat / sqrt ( Tnorm * cp )
         qwgt (k) = qwgt (k) * epsilon  !  change rel.  q-wgt
 10      qwgt (k) = qwgt (k) / Anorm    !  # grd pts in q-term 2/4/97
c
         if ( itype.gt.9 ) then
         print *,'NORM OPTION GREATER THAN 9   itype= ',itype
         stop 'Xnorms0'
         endif
c
         return
         end      
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE Xnorms1 (u,v,t,ps,q,tg,xstar,lmoist,
     &                    sigf,sigh,dsig,grav,R,Cp,clat,clon,dx,
     &                    edepth,ematrixi,gmatrixp,gmatrixr,pweight,
     &                    qmatrix,
     &                    tpmeanb,work1,work2,work3,work4,work5,work6,
     &                    cf,fmapx,fmapd,ni,nj,nk,nxstar,itype)
      logical lmoist
      real u(ni,nj,nk), v(ni,nj,nk), t(ni,nj,nk), ps(ni,nj)
     &,    q(ni,nj,nk), tg(ni,nj), xstar(nxstar)
     &,    sigh(nk), sigf(nk+1), dsig(nk), edepth(nk), pweight(nk+2,3)
     &,    ematrixi(nk,nk), gmatrixp(nk), gmatrixr(nk,nk)
     &,    work1(ni,nj,nk+1), work2(ni,nj,nk+1), work3(ni,nj,nk+1) 
     &,    work4(ni,nj,nk), work5(ni,nj,nk), work6(ni,nj,nk)
     &,    tpmeanb(nk+1), cf(ni,nj), fmapx(ni,nj), fmapd(ni,nj)
      real qwgt(10), qmatrix(nk,nk)
      common /cnorm1/ Tnorm,Pnorm,Anorm
      common /cnorm2/ itrunc,jtrunc,ktrunc
      common /cnorm3/ qwgt
      common /cnorm4/ itruncq,jtruncq,ktruncq,qnorm
c
c	In the notation of E/E-94, this routine does
c	Q (A P Q)^(-1) x. 7/1/94
c
c  Determine initial uncoupled TLM fields from xstar
c
	if ( itype.le.3 .or. itype.eq.5 ) then 
	   q=0.0
	   tg=0.0
           call Mstar1 (u,v,t,ps,xstar,dsig,
     &     R,Cp,Tnorm,Pnorm,Anorm,ni,nj,nk,nxstar,itype)
c
        elseif (itype.eq.4) then 
           q = 0.0
           tg = 0.0
           call PVpv2fs (u,v,t,ps,xstar,work1,work2,work3,work4,
     &                 work5,work6,edepth,ematrixi,gmatrixp,
     &                 gmatrixr,pweight,tpmeanb(nk+1),clat,dx,
     &                 ni,nj,nk,itrunc,jtrunc,ktrunc)
c
        elseif (itype.eq.6) then 
           q = 0.0
           tg = 0.0
           call PVpv2fs (u,v,t,ps,xstar,work1,work2,work3,work4,
     &                 work5,work6,edepth,ematrixi,gmatrixp,
     &                 gmatrixr,pweight,tpmeanb(nk+1),clat,dx,
     &                 ni,nj,nk,itrunc,jtrunc,ktrunc)
           call Mstar1q (q,ni,nj,nk,dsig,
     &                   xstar,nxstar,itrunc,jtrunc,ktrunc,qwgt,itype)
c 
        elseif (itype.eq.7) then
           q=0.
           call Mstar1q (q,ni,nj,nk,dsig,
     &                   xstar,nxstar,itrunc,jtrunc,ktrunc,qwgt,itype)
           u  = 0.0
           v  = 0.0
           t  = 0.0
           ps = 0.0
           tg = 0.0
c
        elseif (itype.eq.8) then
           call Mstar1 (u,v,t,ps,xstar,dsig,
     &     R,Cp,Tnorm,Pnorm,Anorm,ni,nj,nk,nxstar,itype)
           call Mstar1q (q,ni,nj,nk,dsig,
     &                   xstar,nxstar,itrunc,jtrunc,ktrunc,qwgt,itype)       
           tg = 0.0              
c
        elseif (itype.eq.9) then
           call PVpv2fs (u,v,t,ps,xstar,work1,work2,work3,work4,
     &                 work5,work6,edepth,ematrixi,gmatrixp,
     &                 gmatrixr,pweight,tpmeanb(nk+1),clat,dx,
     &                 ni,nj,nk,itrunc,jtrunc,ktrunc)
           nx1=itrunc*jtrunc*ktrunc+1  
           call QNs2f (xstar(nx1),q,qmatrix,work1,work2,qnorm,
     &                  ni,nj,nk,itruncq,jtruncq,ktruncq)
           tg = 0.0
c
        else
         print *,'NORM OPTION GREATER THAN 9   itype= ',itype
         stop 'Xnorms1'
c
      endif
c
      return
      end      
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE Xnorms2a (u,v,t,ps,q,tg,xfull,lmoist,xnorm,
     &                     sigf,sigh,dsig,grav,R,Cp,clat,clon,dx,
     &                     edepth,einverse,gmatrixp,gmatrixt,pweight,
     &                     qinverse,
     &                     tpmeanb,work1,work2,work3,work4,work5,work6,
     &                     cf,fmapx,fmapd,ni,nj,nk,nfull,nxstar,itype)
      logical lmoist
      real u(ni,nj,nk), v(ni,nj,nk), t(ni,nj,nk), ps(ni,nj)
     &,    q(ni,nj,nk), tg(ni,nj), xfull(nfull)
     &,    sigh(nk), sigf(nk+1), dsig(nk), edepth(nk), pweight(nk+2,3)
     &,    einverse(nk,nk), gmatrixp(nk), gmatrixt(nk,nk)
     &,    work1(ni,nj,nk+1), work2(ni,nj,nk+1), work3(ni,nj,nk+1) 
     &,    work4(ni,nj,nk), work5(ni,nj,nk), work6(ni,nj,nk)
     &,    tpmeanb(nk+1), cf(ni,nj), fmapx(ni,nj), fmapd(ni,nj)
      real qwgt(10), qinverse(nk,nk)
      common /cnorm1/ Tnorm,Pnorm,Anorm
      common /cnorm2/ itrunc,jtrunc,ktrunc
      common /cnorm3/ qwgt
      common /cnorm4/ itruncq,jtruncq,ktruncq,qnorm
c
c	in the notation of E/E-94, this routine does the 
c	operator A P. 7/1/94
c
c
c  Multiply uncoupled TLM output fields by scaling Matrix A
c
      if (itype.le.3) then
         call Mstar2 (u,v,t,ps,dsig,xnorm,
     &        R,Cp,Tnorm,Pnorm,Anorm,ni,nj,nk,itype)
c
      elseif (itype.eq.4) then
         call PVfs2pv (xfull,u,v,t,ps,work4,work5,work3,
     &                 work1,work2,edepth,einverse,gmatrixp,
     &                 gmatrixt,tpmeanb(nk+1),clat,dx,ni,nj,nk,
     &                 itrunc,jtrunc,ktrunc) 
         xnorm = 0.0 
         do 10 icount = 1 , nxstar
10	 xnorm = xnorm + xfull (icount) * xfull (icount)
c
      elseif (itype.eq.6) then
         call PVfs2pv (xfull,u,v,t,ps,work4,work5,work3,
     &                 work1,work2,edepth,einverse,gmatrixp,
     &                 gmatrixt,tpmeanb(nk+1),clat,dx,ni,nj,nk,
     &                 itrunc,jtrunc,ktrunc) 
         call Mstar2q (q,ni,nj,nk,dsig,
     1                 xfull,nfull,itrunc,jtrunc,ktrunc,qwgt,itype)
c
         xnorm = 0.0 
         do 20 icount = 1 , nxstar
20	 xnorm = xnorm + xfull (icount) * xfull (icount)
c
      elseif (itype.eq.7) then 
         call Mstar2q (q,ni,nj,nk,dsig,
     1                 xfull,nfull,itrunc,jtrunc,ktrunc,qwgt,itype)
         xnorm = 0.0
         do 30 icount = 1 , nxstar
 30      xnorm = xnorm + xfull (icount) * xfull (icount) 
c 
      elseif (itype.eq.8) then
         call Mstar2 (u,v,t,ps,dsig,xnorm,
     &        R,Cp,Tnorm,Pnorm,Anorm,ni,nj,nk,itype)
         call Mstar2q (q,ni,nj,nk,dsig,
     1                 xfull,nfull,itrunc,jtrunc,ktrunc,qwgt,itype)
         nq = 2*nk*(ni-2)*(nj-2) + (nk+1)*(ni-3)*(nj-3)
         do 40 icount = nq+1 , nxstar
 40      xnorm = xnorm + xfull (icount) * xfull (icount)        
c
      elseif (itype.eq.9) then
         call PVfs2pv (xfull,u,v,t,ps,work4,work5,work3,
     &                 work1,work2,edepth,einverse,gmatrixp,
     &                 gmatrixt,tpmeanb(nk+1),clat,dx,ni,nj,nk,
     &                 itrunc,jtrunc,ktrunc) 
         nx1=itrunc*jtrunc*ktrunc+1  
         call QNf2s (xfull(nx1),q,qinverse,work1,work2,qnorm,
     &                  ni,nj,nk,itruncq,jtruncq,ktruncq)
c
         xnorm = 0.0 
         do 50 icount = 1 , nxstar
50	 xnorm = xnorm + xfull (icount) * xfull (icount)
c
      else
         print *,'NORM OPTION GREATER THAN 9   itype= ',itype
         stop 'Xnorms2a'
c
      endif
c
      return
      end      
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE Xnorms2b (u,v,t,ps,q,tg,xfull,lmoist,
     &                     sigf,sigh,dsig,grav,R,Cp,clat,clon,dx,
     &                     edepth,einversT,gmatrixp,gmatrxtT,pweight,
     &                     qinversT,
     &                     tpmeanb,work1,work2,work3,work4,work5,work6,
     &                     cf,fmapx,fmapd,ni,nj,nk,nfull,nxstar,itype)
      logical lmoist
      real u(ni,nj,nk), v(ni,nj,nk), t(ni,nj,nk), ps(ni,nj)
     &,    q(ni,nj,nk), tg(ni,nj), xfull(nfull)
     &,    sigh(nk), sigf(nk+1), dsig(nk), edepth(nk), pweight(nk+2,3)
     &,    einversT(nk,nk), gmatrixp(nk), gmatrxtT(nk,nk)
     &,    work1(ni,nj,nk+1), work2(ni,nj,nk+1), work3(ni,nj,nk+1) 
     &,    work4(ni,nj,nk), work5(ni,nj,nk), work6(ni,nj,nk)
     &,    tpmeanb(nk+1), cf(ni,nj), fmapx(ni,nj), fmapd(ni,nj)
      real qwgt(10), qinversT(nk,nk)
      common /cnorm1/ Tnorm,Pnorm,Anorm
      common /cnorm2/ itrunc,jtrunc,ktrunc
      common /cnorm3/ qwgt 
      common /cnorm4/ itruncq,jtruncq,ktruncq,qnorm
c
c	in the notation of E/E-94, this routine does
c	the operator P^T A^T. It is therefore the adjoint
c	of xnorms2a.  7/1/94
c
c
c  Determine initial uncoupled adjoint fields from uncoupled) TLM 
c  output fields scaled by scaling matrix A
c
      if ( itype.le.3 ) then
         q=0.
         tg=0. 
c
      else if (itype.eq.4) then
         call PVfs2pvA (xfull,u,v,t,ps,work4,work5,work3,
     &                  work1,work2,edepth,einversT,gmatrixp,
     &                  gmatrxtT,tpmeanb(nk+1),clat,dx,ni,nj,nk,
     &                  itrunc,jtrunc,ktrunc) 
         q=0.
         tg=0.
c
      elseif (itype.eq.6) then
         call PVfs2pvA (xfull,u,v,t,ps,work4,work5,work3,
     &                  work1,work2,edepth,einversT,gmatrixp,
     &                  gmatrxtT,tpmeanb(nk+1),clat,dx,ni,nj,nk,
     &                  itrunc,jtrunc,ktrunc) 
         q=0.
         tg=0.
         call Mstar2qA (q,ni,nj,nk,dsig,
     1                  xfull,nfull,itrunc,jtrunc,ktrunc,qwgt,itype)
c
      elseif (itype.eq.7) then
         u  = 0.0
         v  = 0.0
         t  = 0.0
         ps = 0.0
         q  = 0.0
         tg = 0.0
         call Mstar2qA (q,ni,nj,nk,dsig,
     1                  xfull,nfull,itrunc,jtrunc,ktrunc,qwgt,itype)
c
      elseif (itype.eq.8) then
         q=0.0
         tg=0.0
         call Mstar2qA (q,ni,nj,nk,dsig,
     1                  xfull,nfull,itrunc,jtrunc,ktrunc,qwgt,itype)
c
      elseif (itype.eq.9) then
         call PVfs2pvA (xfull,u,v,t,ps,work4,work5,work3,
     &                  work1,work2,edepth,einversT,gmatrixp,
     &                  gmatrxtT,tpmeanb(nk+1),clat,dx,ni,nj,nk,
     &                  itrunc,jtrunc,ktrunc) 
         q=0.
         tg=0.
         nx1=itrunc*jtrunc*ktrunc+1  
         call QNf2sA (xfull(nx1),q,qinversT,work1,work2,qnorm,
     &                  ni,nj,nk,itruncq,jtruncq,ktruncq)
c
      else
         print *,'NORM OPTION GREATER THAN 9   itype= ',itype
         stop 'Xnorms2b'
      endif
c
      return
      end      
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE Xnorms3 (u,v,t,ps,q,tg,xstar,xfull,xnorm,lmoist,
     &                    sigf,sigh,dsig,grav,R,Cp,clat,clon,dx,
     &                    edepth,ematrxiT,gmatrixp,gmatrxrT,pweight,
     &                    qmatrixT,
     &                    tpmeanb,work1,work2,work3,work4,work5,work6,
     &                    cf,fmapx,fmapd,ni,nj,nk,nxstar,nfull,itype)
      logical lmoist
      real u(ni,nj,nk), v(ni,nj,nk), t(ni,nj,nk), ps(ni,nj)
     &,    q(ni,nj,nk), tg(ni,nj), xstar(nxstar) , xfull (nfull)
     &,    sigh(nk), sigf(nk+1), dsig(nk), edepth(nk), pweight(nk+2,3)
     &,    ematrxiT(nk,nk), gmatrixp(nk), gmatrxrT(nk,nk)
     &,    work1(ni,nj,nk+1), work2(ni,nj,nk+1), work3(ni,nj,nk+1) 
     &,    work4(ni,nj,nk), work5(ni,nj,nk), work6(ni,nj,nk)
     &,    tpmeanb(nk+1), cf(ni,nj), fmapx(ni,nj), fmapd(ni,nj)
      real qwgt(10), qmatrixT(nk,nk)
      common /cnorm1/ Tnorm,Pnorm,Anorm
      common /cnorm2/ itrunc,jtrunc,ktrunc
      common /cnorm3/ qwgt
      common /cnorm4/ itruncq,jtruncq,ktruncq,qnorm
c
c
c	in the notation of E/E-94, this routine does
c	[ ( A P Q ) ^ (-1) ] ^T Q^T. It is therefore the 
c	adjoint of xnorms1.  7/1/94
c
c  File xstar from uncoupled fields
c
      if ( itype.le.3 .or. itype.eq.5 ) then
         call MSTAR3 (u,v,t,ps,xstar,dsig,xnorm,
     &        R,Cp,Tnorm,Pnorm,Anorm,ni,nj,nk,nxstar,itype)
c
      elseif (itype.eq.4) then
         call PVpv2fsA (u,v,t,ps,xfull,work4,work5,work2,work3,work6,
     &                  work1,work2,edepth,ematrxiT,gmatrixp,
     &                  gmatrxrT,pweight,tpmeanb(nk+1),clat,dx,
     &                  ni,nj,nk,itrunc,jtrunc,ktrunc)
c
	xnorm = 0.0
	do 10 icount = 1 , nxstar
	xnorm = xnorm + xstar ( icount ) * xfull ( icount )
10      xstar ( icount ) = xfull ( icount )
c
      elseif (itype.eq.6) then
         call PVpv2fsA (u,v,t,ps,xfull,work4,work5,work2,work3,work6,
     &                  work1,work2,edepth,ematrxiT,gmatrixp,
     &                  gmatrxrT,pweight,tpmeanb(nk+1),clat,dx,
     &                  ni,nj,nk,itrunc,jtrunc,ktrunc)
c
         call mstar3q (q,ni,nj,nk,dsig,
     1                 xfull,nfull,itrunc,jtrunc,ktrunc,qwgt,itype)
c
	xnorm = 0.0
	do 20 icount = 1 , nxstar
	xnorm = xnorm + xstar ( icount ) * xfull ( icount )
20      xstar ( icount ) = xfull ( icount )
c
      elseif ( itype.eq.7 ) then
         call mstar3q (q,ni,nj,nk,dsig,
     1                 xfull,nfull,itrunc,jtrunc,ktrunc,qwgt,itype)
        xnorm = 0.0
        do 30 icount = 1 , nxstar
        xnorm = xnorm + xstar ( icount ) * xfull ( icount )
 30     xstar ( icount ) = xfull ( icount ) 
c
      elseif ( itype.eq.8 ) then
         call MSTAR3 (u,v,t,ps,xstar,dsig,xnorm,
     &        R,Cp,Tnorm,Pnorm,Anorm,ni,nj,nk,nxstar,itype)
         call mstar3q (q,ni,nj,nk,dsig,
     1                 xfull,nfull,itrunc,jtrunc,ktrunc,qwgt,itype)
         nq = 2*nk*(ni-2)*(nj-2) + (nk+1)*(ni-3)*(nj-3)
         do 40 icount = nq + 1 , nxstar
         xnorm = xnorm +  xstar ( icount ) * xfull ( icount )
 40      xstar ( icount ) = xfull ( icount ) 
c         
      elseif (itype.eq.9) then
         call PVpv2fsA (u,v,t,ps,xfull,work4,work5,work2,work3,work6,
     &                  work1,work2,edepth,ematrxiT,gmatrixp,
     &                  gmatrxrT,pweight,tpmeanb(nk+1),clat,dx,
     &                  ni,nj,nk,itrunc,jtrunc,ktrunc)
c
         nx1=itrunc*jtrunc*ktrunc+1  
         call QNs2fA (xfull(nx1),q,qmatrixT,work1,work2,qnorm,
     &                  ni,nj,nk,itruncq,jtruncq,ktruncq)
c
	xnorm = 0.0
	do 50 icount = 1 , nxstar
	xnorm = xnorm + xstar ( icount ) * xfull ( icount )
50      xstar ( icount ) = xfull ( icount )
c
c
      else
         print *,'NORM OPTION GREATER THAN 9   itype= ',itype
         stop 'Xnorms3'
      endif
c
      return
      end   
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE Xnorms4 (u,v,t,ps,q,tg,xstar,xfull,lmoist,
     &                    sigf,sigh,dsig,grav,R,Cp,clat,clon,dx,
     &                    edepth,ematrixi,gmatrixp,gmatrixr,pweight,
     &                    qmatrix,
     &                    tpmeanb,work1,work2,work3,work4,work5,work6,
     &                    cf,fmapx,fmapd,ni,nj,nk,nxstar,nfull,itype)
      logical lmoist
      real u(ni,nj,nk), v(ni,nj,nk), t(ni,nj,nk), ps(ni,nj)
     &,    q(ni,nj,nk), tg(ni,nj), xstar (nxstar) , xfull (nfull)
     &,    sigh(nk), sigf(nk+1), dsig(nk), edepth(nk), pweight(nk+2,3)
     &,    ematrixi(nk,nk), gmatrixp(nk), gmatrixr(nk,nk)
     &,    work1(ni,nj,nk+1), work2(ni,nj,nk+1), work3(ni,nj,nk+1) 
     &,    work4(ni,nj,nk), work5(ni,nj,nk), work6(ni,nj,nk)
     &,    tpmeanb(nk+1), cf(ni,nj), fmapx(ni,nj), fmapd(ni,nj)
      real qwgt(10), qmatrix(nk,nk) 
      common /cnorm1/ Tnorm,Pnorm,Anorm
      common /cnorm2/ itrunc,jtrunc,ktrunc
      common /cnorm3/ qwgt
      common /cnorm4/ itruncq,jtruncq,ktruncq,qnorm
c
c	In the notation of E/E-94, this routine does
c	Q (A P Q)^(-1) x and stores on a vector. 7/1/94
c
c  Determine initial uncoupled TLM fields from xstar
c
	if ( itype.le.3 .or. itype.eq.5 ) then 
           xfull = 0.0
           call Mstar4 (u,v,t,ps,xstar,dsig,
     &     R,Cp,Tnorm,Pnorm,Anorm,ni,nj,nk,nxstar,
     &     xfull,nfull,itype)
c
        elseif (itype.eq.4) then 
           xfull = 0.0   ! 6/14/95 to set q-field to zero, m.e.
           call PVpv2fsd (u,v,t,ps,xstar,work1,work2,work3,work4,
     &                 work5,work6,edepth,ematrixi,gmatrixp,
     &                 gmatrixr,pweight,tpmeanb(nk+1),clat,dx,
     &                 ni,nj,nk,itrunc,jtrunc,ktrunc,xfull,nfull)
c
        elseif (itype.eq.6) then 
           call PVpv2fsd (u,v,t,ps,xstar,work1,work2,work3,work4,
     &                 work5,work6,edepth,ematrixi,gmatrixp,
     &                 gmatrixr,pweight,tpmeanb(nk+1),clat,dx,
     &                 ni,nj,nk,itrunc,jtrunc,ktrunc,xfull,nfull)
           call mstar4q (q,ni,nj,nk,dsig,xstar,nxstar,
     1                   xfull,nfull,itrunc,jtrunc,ktrunc,qwgt,itype)
c
        elseif (itype.eq.7) then
           xfull = 0.0
           call mstar4q (q,ni,nj,nk,dsig,xstar,nxstar,
     1                   xfull,nfull,itrunc,jtrunc,ktrunc,qwgt,itype)
c
        elseif (itype.eq.8) then
           xfull = 0.0
           call Mstar4 (u,v,t,ps,xstar,dsig,
     &     R,Cp,Tnorm,Pnorm,Anorm,ni,nj,nk,nxstar,
     &     xfull,nfull,itype)           
           call mstar4q (q,ni,nj,nk,dsig,xstar,nxstar,
     1                   xfull,nfull,itrunc,jtrunc,ktrunc,qwgt,itype)
c
        elseif (itype.eq.9) then 
           call PVpv2fsd (u,v,t,ps,xstar,work1,work2,work3,work4,
     &                 work5,work6,edepth,ematrixi,gmatrixp,
     &                 gmatrixr,pweight,tpmeanb(nk+1),clat,dx,
     &                 ni,nj,nk,itrunc,jtrunc,ktrunc,xfull,nfull)
           nx1=itrunc*jtrunc*ktrunc+1
           call QNs2f (xstar(nx1),q,qmatrix,work1,work2,qnorm,
     &                  ni,nj,nk,itruncq,jtruncq,ktruncq)
           call QNdimv (xfull,q,ni,nj,nk,nfull)
c
        else
         print *,'NORM OPTION GREATER THAN 9   itype= ',itype
         stop 'Xnorms4'
c
      endif
c
      return
      end      
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
c
c   what follows are the 'energy.f' routines in the modified 
c   form to work with itype, 8 july 94, m.e.
c
C
c
c
c
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE MSTAR1 (u,v,t,ps,xstar,dsig,
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     &                   R,Cp,Tnorm,Pnorm,Anorm,ni,nj,nk,nxstar,itype)
C
      real u(ni,nj,nk), v(ni,nj,nk), t(ni,nj,nk), ps(ni,nj)
     &,    xstar(nxstar), dsig(nk)
c
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
      call Cbdryx (u,v,t,ps,ni,nj,nk)
      afac=sqrt(2.)*Anorm
c
c  Fill and unscale ps 
      nx=nxp
      fac=afac*Pnorm/sqrt(R*Tnorm)
      do 1 j=2,nj2
      do 1 i=2,ni2
      nx=nx+1
c      if ( iset3.eq.1 ) then
       if ( itype.eq.3 ) then
      		ps (i,j) = 0.0
	else
      		ps(i,j)=fac*xstar(nx)
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
c      if ( iset3.eq.1 ) then 
       if ( itype.eq.3 ) then
      		t(i,j,k) = 0.0
	else
      		t(i,j,k)=fac*xstar(nx)
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
      u(i,j,k)=fac*xstar(nxu)
      v(i,j,k)=fac*xstar(nxv)
    3 continue
c
      return
      end
C
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE MSTAR2 (u,v,t,ps,dsig,xnorm,
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     &                   R,Cp,Tnorm,Pnorm,Anorm,ni,nj,nk,itype)
      real u(ni,nj,nk), v(ni,nj,nk), t(ni,nj,nk), ps(ni,nj)
     &,    dsig(nk)
c
c  Multiply TLM result by scale factor matrix A**2
c  Also compute norm of result after multiplication by A alone
c  
      ni1=ni-1
      nj1=nj-1
      ni2=ni-2
      nj2=nj-2
      afac=1./(Anorm*sqrt(2.))
      xnorm=0.
c
	zero = 1.0
c	if ( iset2.eq.1 ) zero = 0.0
c	if ( iset3.eq.1 ) zero = 0.0
        if ( itype.eq.2 ) zero = 0.0
        if ( itype.eq.3 ) zero = 0.0
c
c  Scale ps 
      fac=afac*sqrt(R*Tnorm)/Pnorm
      do 1 j=2,nj2
      do 1 i=2,ni2
      ps (i,j) = zero * ps (i,j)
      f=fac*ps(i,j)
      xnorm=xnorm+f*f
      ps(i,j)=fac*f
    1 continue
c  
c  Multiply t 
      fac1=afac*sqrt(Cp/Tnorm)
      do 2 k=1,nk
      fac=fac1*sqrt(dsig(k))
      do 2 j=2,nj2
      do 2 i=2,ni2
      t (i,j,k) = zero * t (i,j,k) 
      f=fac*t(i,j,k)
      xnorm=xnorm+f*f
      t(i,j,k)=fac*f
    2 continue
c  
c  Multiply u and v
      do 3 k=1,nk
      fac=afac*sqrt(dsig(k)) 
      do 3 j=2,nj1
      do 3 i=2,ni1
      fu=fac*u(i,j,k)
      fv=fac*v(i,j,k)
      xnorm=xnorm+fu*fu+fv*fv
      u(i,j,k)=fac*fu
      v(i,j,k)=fac*fv
    3 continue
c
      return
      end
C
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE MSTAR3 (u,v,t,ps,xstar,dsig,xnorm,
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     &                   R,Cp,Tnorm,Pnorm,Anorm,ni,nj,nk,nxstar,itype)
c
      real u(ni,nj,nk), v(ni,nj,nk), t(ni,nj,nk), ps(ni,nj)
     &,    xstar(nxstar), dsig(nk)
c
c  Scale (uncoupled) results and overwrite xstar
c  Also compute norm of old dot new xstar values
c
      ni1=ni-1
      nj1=nj-1
      ni2=ni-2
      nj2=nj-2
      nxu=0
      nxv=nxu+ni2*nj2*nk
      nxt=nxv+ni2*nj2*nk
      nxp=nxt+(ni2-1)*(nj2-1)*nk
      afac=sqrt(2.)*Anorm
      xnorm=0.
c
c  Fill xstar for ps
      nx=nxp
      fac=afac*Pnorm/sqrt(R*Tnorm)
      do 1 j=2,nj2
      do 1 i=2,ni2
      nx=nx+1
      f=fac*ps(i,j)
c      if ( iset3.eq.1 ) then 
       if ( itype.eq.3 ) then
      		continue
	else
      		xnorm=xnorm+f*xstar(nx)
      		xstar(nx)=f
      endif
    1 continue
c
c  Fill xstar for t
      nx=nxt
      fac1=afac*sqrt(Tnorm/Cp)
      do 2 k=1,nk
      fac=fac1/sqrt(dsig(k)) 
      do 2 j=2,nj2
      do 2 i=2,ni2
      nx=nx+1
      f=t(i,j,k)*fac
c      if ( iset3.eq.1 ) then 
       if ( itype.eq.3 ) then
      		continue
	else
      		xnorm=xnorm+f*xstar(nx)
      		xstar(nx)=f
      endif
    2 continue
c  
c  Fill xstar for u and v
      do 3 k=1,nk
      fac=afac/sqrt(dsig(k)) 
      do 3 j=2,nj1
      do 3 i=2,ni1
      nxu=nxu+1
      nxv=nxv+1
      fu=u(i,j,k)*fac
      fv=v(i,j,k)*fac
      xnorm=xnorm+fu*xstar(nxu)+fv*xstar(nxv)
      xstar(nxu)=fu
      xstar(nxv)=fv
    3 continue
c
      return
      end
C
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C 	SUBROUTINE MSTAR4 IS A COPY OF MSTAR1 WITH THE EXCEPTION 
C	THAT IT PUTS THE DIMENSIONED VECTOR BACK ON XSTAR 
C
C	M. EHRENDORFER, 15 SEPTEMBER 1993. 
C
C	MODIFIED ON 15 DECEMBER 1993.
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE MSTAR4 (u,v,t,ps,xstar,dsig,
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     &                   R,Cp,Tnorm,Pnorm,Anorm,ni,nj,nk,nxstar,
     &                   xfull,nfull,itype)
c
      real u(ni,nj,nk), v(ni,nj,nk), t(ni,nj,nk), ps(ni,nj)
     &,    xstar(nxstar), dsig(nk)
C
	REAL XFULL ( NFULL )
c
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
      call Cbdryx (u,v,t,ps,ni,nj,nk)
      afac=sqrt(2.)*Anorm
c  
c  Fill and unscale ps 
      nx=nxp
      fac=afac*Pnorm/sqrt(R*Tnorm)
      do 1 j=2,nj2
      do 1 i=2,ni2
      nx=nx+1
CCC      ps(i,j)=fac*xstar(nx)
c      if ( iset3.eq.1 ) then 
       if ( itype.eq.3 ) then
      		xfull (nx) = 0.0
      	else
      		XFULL(nx)=fac*xstar(nx)
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
CCC      t(i,j,k)=fac*xstar(nx)
c      if ( iset3.eq.1 ) then
       if ( itype.eq.3 ) then
      		xfull (nx) = 0.0
	else
      		XFULL(nx)=fac*xstar(nx)
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
CCC      u(i,j,k)=fac*xstar(nxu)
CCC      v(i,j,k)=fac*xstar(nxv)
      XFULL(nxu)=fac*xstar(nxu)
      XFULL(nxv)=fac*xstar(nxv)
    3 continue
c
      return
      end
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C
c
c
c
c
c
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE PDECUPL (pu,pv,pt,pq,ps,u,v,t,q,
     &                    pub,pvb,ptb,pqb,psb,ub,vb,tb,qb,
     &                    pd,pdb,lmoist,ni,nj,nk,itype)
      logical lmoist
      real  pu(ni,nj,nk),  pv(ni,nj,nk),  pt(ni,nj,nk),  pq(ni,nj,nk)
     &,       u(ni,nj,nk),   v(ni,nj,nk),   t(ni,nj,nk),   q(ni,nj,nk)
     &,     pub(ni,nj,nk), pvb(ni,nj,nk), ptb(ni,nj,nk), pqb(ni,nj,nk)
     &,      ub(ni,nj,nk),  vb(ni,nj,nk),  tb(ni,nj,nk),  qb(ni,nj,nk)
     &,      ps(ni,nj), psb(ni,nj), pd(ni,nj), pdb(ni,nj) 
c
c  Decouple basic state
      call Bx2d (pdb,psb,ni,nj,1)
      call Bdecup (ub,pub,pdb,ni,nj,nk,0)
      call Bdecup (vb,pvb,pdb,ni,nj,nk,0)
      call Bdecup (tb,ptb,psb,ni,nj,nk,1)
      if (lmoist) call Bdecup (qb,pqb,psb,ni,nj,nk,1)
c
c  TLM of Decouple
      if (itype.eq.1) then
         call Bx2d (pd,ps,ni,nj,1)
         call Fdecup (u,pu,pd,ub,pdb,ni,nj,nk,0)
         call Fdecup (v,pv,pd,vb,pdb,ni,nj,nk,0)
         call Fdecup (t,pt,ps,tb,psb,ni,nj,nk,1)
         if (lmoist) call Fdecup (q,pq,ps,qb,psb,ni,nj,nk,1)
c
c  TLM of coupling
      elseif (itype.eq.2) then
         call Bx2d (pd,ps,ni,nj,1)
         call Fcouple (pu,u,pd,ub,pdb,ni,nj,nk,0,1) 
         call Fcouple (pv,v,pd,vb,pdb,ni,nj,nk,0,1)
         call Fcouple (pt,t,ps,tb,psb,ni,nj,nk,1,1)
         if (lmoist) call Fcouple (pq,q,ps,qb,psb,ni,nj,nk,1,1)
c
c  ADJ of Decouple
      elseif (itype.eq.3) then
         pd=0.
         call Adecup (u,pu,pd,ub,pdb,ni,nj,nk,0)
         call Adecup (v,pv,pd,vb,pdb,ni,nj,nk,0)
         call Adecup (t,pt,ps,tb,psb,ni,nj,nk,1)
         if (lmoist) call Adecup (q,pq,ps,qb,psb,ni,nj,nk,1)
         call Ax2d (pd,ps,ni,nj,1,.f.)
c
c  ADJ of coupling
      elseif (itype.eq.4) then
         pd=0.
         call Acouple (pu,u,pd,ub,pdb,ni,nj,nk,0,1) 
         call Acouple (pv,v,pd,vb,pdb,ni,nj,nk,0,1)
         call Acouple (pt,t,ps,tb,psb,ni,nj,nk,1,1)
         if (lmoist) call Acouple (pq,q,ps,qb,psb,ni,nj,nk,1,1)
         call Ax2d (pd,ps,ni,nj,1,.f.)
      endif
c
      return
      end
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE Ndispose (cMSdrct,cflttr,itapnum,lcopy,keepdays
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     &,                    ireclw,cusrcom,lprint)

c cflttr = 'D'ata,'C'onvective, 'N'onconv, 
c          'O'ther, or 'R'adiation file signifier


      character cMSdrct*(*),cflttr*1,ccomm*80,cflnm*128,cusrcom*65
     &,         ccopy*8,ckeep*3,creclw*6,clocnam*9,ctapnum*2
      logical lprint, lcopy
 
      write(ctapnum,'(I2)') itapnum
      if (itapnum.lt.10) ctapnum = ctapnum(2:2)//' '

      clocnam = cflttr//'output'//ctapnum

      nch = LNGTH(cMSdrct)
      if (cflttr.NE.'O') then
         cflnm = cMSdrct(1:nch)//'NLM'//cflttr//ctapnum
      else if (cflttr.EQ.'O') then
         cflnm = cMSdrct(1:nch)
      else
         print*,' '
         print*,'ERROR: invalid file type=',cflttr,
     &          ' specified in call to Ndispose'
      end if
         
      if (lcopy) then
         ccopy = ' copy   '
      else
         ccopy = ' nocopy '
      end if

      write(ckeep,'(I3)') keepdays
      write(creclw,'(I6)') ireclw

      nbl=0
      do i=1,6
         if (creclw(i:i).eq.' ') then
            nbl=nbl+1
         else
            go to 11
         end if
      end do
 11   ccomm = ' recsiz_'//creclw(nbl+1:6)//cusrcom
      
      ierrdis = ISHELL('dispose '//clocnam//' '//cflnm
     &                //ccopy//ckeep//ccomm)

      if (ierrdis.eq.0 .and. lprint) then
         print *,'File disposed: local name =',clocnam
         print *,'File name for dispose = ',cflnm
      else if (ierrdis.ne.0) then
         print*,'ERROR: FAILURE OF DISPOSE OF ',clocnam,
     &                                        ' ierrdis = ',ierrdis
         stop 'Ndispose'
      end if

      itapnum = itapnum + 1

      return
      end
C 
C 
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      FUNCTION LNGTH(S)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      CHARACTER*(*) S
      LNGTH=LEN(S)
    1 IF (S(LNGTH:LNGTH).NE.' ') RETURN
      IF (LNGTH.EQ.1) RETURN
      LNGTH=LNGTH-1
      GOTO 1
      END
C 
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
c
c
c     moisture routines for moist r-norm, 12/1/94, m.e.
c
c     modified to account for q-only norm, 12/20/94, m.e. 
c
c     modified to add option 8 (moist TE-norm), 1/23/96, m.e.
c
c
c
      subroutine mstar1q (q,ni,nj,nk,dsig,
     1                    xstar,nxstar,itrunc,jtrunc,ktrunc,qwgt,itype)
c
c     moist r-norm, m.e. 12/1/94
c
      real q (ni,nj,nk) , xstar (nxstar) , dsig (nk) 
      real qwgt (nk) 
c
      ni2 = ni - 2
      nj2 = nj - 2
      if (itype.eq.6) nq = itrunc * jtrunc * ktrunc
      if (itype.eq.7) nq = 0
      if (itype.eq.8) nq = 2*nk*(ni-2)*(nj-2) + (nk+1)*(ni-3)*(nj-3)
c
      call Cbdry0 (q,ni,nj,nk,1)
c
      do 10 k = 1 , nk
      fac = sqrt ( 0.5 * dsig (k) ) * qwgt (k) 
      fac = 1.0 / fac
      do 10 j = 2 , nj2
      do 10 i = 2 , ni2
      nq = nq + 1
      q (i,j,k) = fac * xstar ( nq )
 10   continue
c
      return
      end
c
c
      subroutine mstar2q (q,ni,nj,nk,dsig,
     1                    xfull,nfull,itrunc,jtrunc,ktrunc,qwgt,itype)
c
c     moist r-norm, m.e. 12/1/94
c
      real q (ni,nj,nk), xfull (nfull), dsig (nk) 
      real qwgt (nk) 
c
      ni2 = ni - 2
      nj2 = nj - 2
      if ( itype.eq.6 ) nq  = itrunc * jtrunc * ktrunc
      if ( itype.eq.7 ) nq = 0
      if ( itype.eq.8 ) nq = 2*nk*(ni-2)*(nj-2) + (nk+1)*(ni-3)*(nj-3)
c
      do 10 k = 1 , nk
      fac =  sqrt ( 0.5 * dsig (k) ) * qwgt (k) 
      do 10 j = 2 , nj2
      do 10 i = 2 , ni2
      nq = nq + 1
      xfull (nq) = fac * q (i,j,k) 
 10   continue
c
      return
      end
c
c
c
      subroutine mstar2qA (q,ni,nj,nk,dsig,
     1                     xfull,nfull,itrunc,jtrunc,ktrunc,qwgt,itype)
c
c     moist r-norm, m.e. 12/1/94
c
      real q (ni,nj,nk), xfull (nfull), dsig (nk) 
      real qwgt (nk) 
c
      ni2 = ni - 2
      nj2 = nj - 2
      if ( itype.eq.6 ) nq  = itrunc * jtrunc * ktrunc
      if ( itype.eq.7 ) nq = 0
      if ( itype.eq.8 ) nq = 2*nk*(ni-2)*(nj-2) + (nk+1)*(ni-3)*(nj-3)
c
      call Cbdry0 (q,ni,nj,nk,1)
c
      do 10 k = 1 , nk
      fac =  sqrt ( 0.5 * dsig (k) ) * qwgt (k) 
      do 10 j = 2 , nj2
      do 10 i = 2 , ni2
      nq = nq + 1
      q (i,j,k) = fac * xfull (nq) 
 10   continue
c
      return
      end
c
c
c
      subroutine mstar3q (q,ni,nj,nk,dsig,
     1                    xfull,nfull,itrunc,jtrunc,ktrunc,qwgt,itype)
c
c     moist r-norm, m.e. 12/1/94
c
      real q (ni,nj,nk) , xfull (nfull) , dsig (nk) 
      real qwgt (nk) 
c
      ni2 = ni - 2
      nj2 = nj - 2
      if ( itype.eq.6 ) nq = itrunc * jtrunc * ktrunc
      if ( itype.eq.7 ) nq = 0
      if ( itype.eq.8 ) nq = 2*nk*(ni-2)*(nj-2) + (nk+1)*(ni-3)*(nj-3)
c
      do 10 k = 1 , nk
      fac = sqrt ( 0.5 * dsig (k) ) * qwgt (k) 
      fac = 1.0 / fac
      do 10 j = 2 , nj2
      do 10 i = 2 , ni2
      nq = nq + 1
      xfull (nq) = fac * q (i,j,k) 
 10   continue
c
      return
      end
c
c
c
      subroutine mstar4q (q,ni,nj,nk,dsig,xstar,nxstar,
     1                    xfull,nfull,itrunc,jtrunc,ktrunc,qwgt,itype)
c
c     moist r-norm, m.e. 12/1/94
c
      real q (ni,nj,nk) , xfull (nfull) , dsig (nk) 
      real qwgt (nk) , xstar (nxstar) 
c
      ni2 = ni - 2
      nj2 = nj - 2
c
      nxu = 0
      nxv = nxu + ni2 * nj2 * nk
      nxt = nxv + ni2 * nj2 * nk
      nxp = nxt + (ni2-1) * (nj2-1) * nk
      nq  = nxp + (ni2-1) * (nj2-1)
      if ( itype.eq.6 ) nq1 = itrunc*jtrunc*ktrunc
      if ( itype.eq.7 ) nq1 = 0
      if ( itype.eq.8 ) nq1 = nq
c
      do 10 k = 1 , nk
      fac = sqrt ( 0.5 * dsig (k) ) * qwgt (k) 
      fac = 1.0 / fac
      do 10 j = 2 , nj2
      do 10 i = 2 , ni2
      nq  = nq  + 1
      nq1 = nq1 + 1
      xfull (nq) = fac * xstar ( nq1 ) 
  10  continue
c
      return
      end
C
C     version 3 Feb 97 ron errico R-norm
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE PVpv2fs (u,v,t,ps,pvort,geop,psi,um,vm,
     &                    work1,work2,edepth,ematrixi,gmatrixp,
     &                    gmatrixr,pweight,pmeanb,clat,dx,
     &                    ni,nj,nk,itrunc,jtrunc,ktrunc) 
C
C  Compute u, v, T, ps from pvort in spectral space
C  Pvort is field of rotational modes.
C
      parameter (tpratio=0.)
      real pvort(itrunc,jtrunc,ktrunc)
     &,    geop(ni,nj,nk), psi(ni,nj,nk), um(ni,nj,nk)
     &,    work1(ni,nj,nk), work2(ni,nj,nk), vm(ni,nj,nk)
     &,    u(ni,nj,nk), v(ni,nj,nk), t(ni,nj,nk), ps(ni,nj)
     &,    edepth(nk), ematrixi(nk,nk), pweight(nk+2,3)
     &,    gmatrixp(nk), gmatrixr(nk,nk) 
c
c
      nij=ni*nj
      nijk=nij*nk
      pi4=atan(1.)
      f0=16.*pi4*sin(clat*pi4/45.)/(60.*60.*24)
      rf0=1./f0
      rpmeanb=1./pmeanb
      ni1=ni-1
      nj1=nj-1
c
c  compute geop in vertical mode space from pvort in spectral space 
      call Nrsetc (geop,nijk,0.)
      call PVpv2gp (geop,pvort,edepth,work1,f0,dx,
     &              ni,nj,itrunc,jtrunc,ktrunc)
c
c  compute psi in vertical mode space
      call Ncopyv (psi,geop,nijk)
      call Nmultv (psi,nijk,rf0)
c
c  compute change in ps (work1=geop change on sigma)
      call Nrsetc (ps,nij,0.)
      call Fchangep (ps,work2,ps,geop,pweight,ematrixi,
     &               gmatrixr,work1,tpratio,ni,nj,nk,ktrunc,
     &               work2,work2,.f.,.true.)
c
c  compute new t 
      call Nrsetc (work2,nijk,0.)
      call Nrsetc (t,nijk,0.)
      call Fchanget (t,work1,t,work2,geop,ps,gmatrixp,gmatrixr,
     &               work2,work2,0.,rpmeanb,ni,nj,nk,ktrunc,.f.)
c
c  compute u and v in vertical mode space (store in um, vm)
      call Nrsetc (work1,nijk,0.)  ! set divg=0.
      call Nrsetc (um,nijk,0.)  
      call Nrsetc (vm,nijk,0.)  
      call Ichanguv (um,vm,psi,work1,dx,ni,nj,ktrunc)
c
c  compute u and v in real space
      call Sproject (u,um,ematrixi,ni,nj,nk,ktrunc,nk,0,2,.t.)
      call Sproject (v,vm,ematrixi,ni,nj,nk,ktrunc,nk,0,2,.t.)
c
      return
      end
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE PVpv2fsA (u,v,t,ps,pvort,geop,psi,um,vm,delp,
     &                     work1,work2,edepth,ematrixi,gmatrixp,
     &                     gmatrixr,pweight,pmeanb,clat,dx,ni,nj,nk,
     &                     itrunc,jtrunc,ktrunc) 
C
C  Compute adjoint of routine = PVpv2fs
C
      parameter (tpratio=0.)
      real pvort(itrunc,jtrunc,ktrunc)
     &,    geop(ni,nj,nk), psi(ni,nj,nk), um(ni,nj,nk)
     &,    work1(ni,nj,nk), work2(ni,nj,nk), vm(ni,nj,nk)
     &,    u(ni,nj,nk), v(ni,nj,nk), t(ni,nj,nk), ps(ni,nj)
     &,    edepth(nk), ematrixi(nk,nk), pweight(nk+2,3)
     &,    gmatrixp(nk), gmatrixr(nk,nk), delp(ni,nj) 
c
c
      nij=ni*nj
      nijk=nij*nk
      pi4=atan(1.)
      f0=16.*pi4*sin(clat*pi4/45.)/(60.*60.*24)
      rf0=1./f0
      nijkmd=nij*ktrunc
      rpmeanb=1./pmeanb
      ni1=ni-1
      nj1=nj-1
c
c  compute u^ and v^ in mode space
      call Sproject (um,u,ematrixi,ni,nj,ktrunc,nk,nk,0,2,.t.)
      call Sproject (vm,v,ematrixi,ni,nj,ktrunc,nk,nk,0,2,.t.)
c
c  compute psi^ in vertical mode space 
      call Achanguv (um,vm,psi,work1,dx,ni,nj,ktrunc)
c
c  compute geop(change)^ and ps^
      call Nrsetc (work2,nijk,0.)
      call Achanget (t,work1,t,work2,geop,delp,gmatrixp,gmatrixr,
     &               work2,work2,0.,rpmeanb,ni,nj,nk,ktrunc,.f.)
c
c  increment geop(change)^
      call Achangep (ps,delp,ps,geop,pweight,ematrixi,
     &               gmatrixr,work1,tpratio,ni,nj,nk,ktrunc,
     &               work1,work1)
c
c  add contribution to geop^ by psi^ in vertical mode space
      call Nmultv (psi,nijkmd,rf0)
      call Naddv (geop,geop,psi,nijkmd)
c
c  compute pvort^ (in spectral space)
      call PVpv2gpA (geop,pvort,edepth,work1,f0,dx,ni,nj,
     &               itrunc,jtrunc,ktrunc)
c
c
      return
      end
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE PVfs2pv (pvort,u,v,t,ps,geop,vort,fmap1,
     &                    work1,work2,edepth,einverse,gmatrixp,
     &                    gmatrixt,pmeanb,clat,dx,ni,nj,nk,
     &                    itrunc,jtrunc,ktrunc) 
C
C  Compute pvort in spectral space from u,v,t,ps
C
      real pvort(itrunc,jtrunc,ktrunc)
     &,    geop(ni,nj,nk), vort(ni,nj,nk)
     &,    work1(ni,nj,nk), work2(ni,nj,nk), fmap1(ni,nj)
     &,    u(ni,nj,nk), v(ni,nj,nk), t(ni,nj,nk), ps(ni,nj)
     &,    edepth(nk)
     &,    einverse(nk,nk), gmatrixp(nk), gmatrixt(nk,nk) 
c
c
      nij=ni*nj
      nijk=nij*nk
      pi4=atan(1.)
      f0=16.*pi4*sin(clat*pi4/45.)/(60.*60.*24)
      rf0=1./f0
      rpmeanb=1./pmeanb
      nkmd=ktrunc
      nijkmd=nij*nkmd
      call Nrsetc (fmap1,nij,1.)
c
c  compute vertical mode amplitude vorticity
      call Sproject (work1,u,einverse,ni,nj,nkmd,nk,nk,0,2,.t.)
      call Sproject (work2,v,einverse,ni,nj,nkmd,nk,nk,0,2,.t.)
      call Bvort (vort,work1,work2,fmap1,dx,ni,nj,nkmd)
c
c  compute vertical mode amplitude of pseudo-geopotential
      call Sproject (geop,ps,gmatrixp,ni,nj,nkmd,1,nk,1,1,.t.)
      call Sproject (geop,t,gmatrixt,ni,nj,nkmd,nk,nk,1,1,.f.)
      call Nmultv (geop,nijkmd,rpmeanb) 
c
c  compute horizontal spectra of pvort (work1=pvort in v.m. space)
      call Ipotvort (work1,vort,geop,edepth,f0,ni,nj,nkmd)
      call Cbdry0 (work1,ni,nj,nkmd,1)
      call PVp2spc (pvort,work1,work2,edepth,f0,dx,ni,nj,
     &              itrunc,jtrunc,ktrunc)
c
      return
      end
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE PVfs2pvA (pvort,u,v,t,ps,geop,vort,fmap1,
     &                     work1,work2,edepth,einverse,gmatrixp,
     &                     gmatrixt,pmeanb,clat,dx,ni,nj,nk,
     &                     itrunc,jtrunc,ktrunc) 
C
C  Comput adjoint of routine = PVfs2pv
C
      real pvort(itrunc,jtrunc,ktrunc)
     &,    geop(ni,nj,nk), vort(ni,nj,nk)
     &,    work1(ni,nj,nk), work2(ni,nj,nk), fmap1(ni,nj)
     &,    u(ni,nj,nk), v(ni,nj,nk), t(ni,nj,nk), ps(ni,nj)
     &,    edepth(nk)
     &,    einverse(nk,nk), gmatrixp(nk), gmatrixt(nk,nk) 
c
c
      nij=ni*nj
      nijk=nij*nk
      pi4=atan(1.)
      f0=16.*pi4*sin(clat*pi4/45.)/(60.*60.*24)
      rf0=1./f0
      rpmeanb=1./pmeanb
      nkmd=ktrunc
      nijkmd=nij*nkmd
      call Nrsetc (fmap1,nij,1.)
c
c  compute vort^ and geop^ in v.m. space from pvort^ in spectral space
c  work1 is pvort^ in vertical mode space 
      call PVp2spcA (pvort,work1,work2,edepth,f0,dx,ni,nj,
     &               itrunc,jtrunc,ktrunc)
      call Apotvort (work1,vort,geop,edepth,f0,ni,nj,nkmd)
c
c  compute t^ and ps^ in physical space from geop^ in v.m. space
      call Nmultv (geop,nijkmd,rpmeanb) 
      call Sproject (t,geop,gmatrixt,ni,nj,nk,nkmd,nk,1,2,.t.)
      call Sproject (ps,geop,gmatrixp,ni,nj,1,nkmd,1,1,2,.t.)
c
c  compute u^ and v^ in physical space from vort^ in v.m. space
      call Nrsetc (work1,nijkmd,0.)
      call Nrsetc (work2,nijkmd,0.)
      call Aavort (vort,work1,work2,fmap1,dx,ni,nj,nkmd)
      call Sproject (u,work1,einverse,ni,nj,nk,nkmd,nk,0,2,.t.)
      call Sproject (v,work2,einverse,ni,nj,nk,nkmd,nk,0,2,.t.)
c
      call Cbdryx (u,v,t,ps,ni,nj,nk)
c
      return
      end
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE PVpv2gp (geop,pvorts,edepth,work,f0,dx,
     &                    ni,nj,itrunc,jtrunc,ktrunc)
      real pvorts(itrunc,jtrunc,ktrunc), geop(ni,nj,ktrunc),
     &     edepth(ktrunc), work(*)
c
c  compute spectral geop from spectral pvort by inverting operator
c  (1/f0)*(del**2 - xlambda) in spectral space
c  then project into vertical mode space on the grid
c
c  sfact is sqrt of 2x2x2; one factor for definition of energy;
c    and a factors of 2 for integrals of sin**2 in each direction
c
      rdx2=2./(dx*dx)
      pi=4.*atan(1.)
      ni1=ni-1
      nj1=nj-1
      nki=ni1-2
      nkj=nj1-2
      xki=pi/(ni1-1)
      xkj=pi/(nj1-1)
      f02=f0*f0
      sfac1=-f0*sqrt(8.)  !qqqqqqqqqqq sfac
c
      do 2 kl=1,ktrunc
      xlambda=f02/edepth(kl)
cqqqqqqqqq
      xxkk=k1
      xxkk=xxkk**.75
      sfac=sfac1/xxkk
cqqqqqqqqqqqq
      do 1 kj=1,jtrunc
      cosj=cos(kj*xkj)
      do 1 ki=1,itrunc
      cosi=cos(ki*xki)
      fac=sqrt(xlambda-(cosi*cosj-1.)*rdx2)
      geop(ki,kj,kl)=sfac*pvorts(ki,kj,kl)/fac 
    1 continue
c
      call Iftproj (geop(1,1,kl),geop(1,1,kl),work,ni,ni1,nj1,.t.)
    2 continue
c
      return
      end
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE PVpv2gpA (geop,pvorts,edepth,work,f0,dx,ni,nj,
     &                     itrunc,jtrunc,ktrunc)
      real pvorts(itrunc,jtrunc,ktrunc), geop(ni,nj,ktrunc),
     &     edepth(ktrunc), work(*)
c
c  Compute adjoint of routine = PVpv2gp
c  fgrid removes factor in Ifttran
c
      rdx2=2./(dx*dx)
      pi=4.*atan(1.)
      ni1=ni-1
      nj1=nj-1
      xki=pi/(ni1-1)
      xkj=pi/(nj1-1)
      f02=f0*f0
      fgrid=(ni-2)*(nj-2)*0.25
      f0fgrid1=-f0*fgrid*sqrt(8.) !qqqqqqqqqqqqqqqq
c
      do 1 kl=1,ktrunc
      xlambda=f02/edepth(kl)
cqqqqqqqqq
      xxkk=k1
      xxkk=xxkk**.75
      f0fgrid=f0fgrid1/xxkk
cqqqqqqqqqqqq
      call Ifttran (geop(1,1,kl),geop(1,1,kl),work,ni,ni1,nj1)
      do 1 kj=1,jtrunc
      cosj=cos(kj*xkj)
      do 1 ki=1,itrunc
      cosi=cos(ki*xki)
      fac=sqrt(xlambda-(cosi*cosj-1.)*rdx2)
      pvorts(ki,kj,kl)=f0fgrid*geop(ki,kj,kl)/fac 
    1 continue
c
      return
      end
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE PVp2spc (pvorts,pvort,work,edepth,f0,dx,
     &                    ni,nj,itrunc,jtrunc,ktrunc)
      real pvorts(itrunc,jtrunc,ktrunc), pvort(ni,nj,ktrunc),
     &     work(ni,nj,2), edepth(ktrunc)
c
c  compute spectral pv from pv in vertical mode space on grid
c
c  sfact is sqrt of 2x2x2; one factor for definition of energy;
c    and a factors of 2 for integrals of sin**2 in each direction
c
      pi=4.*atan(1.)
      ni1=ni-1
      nj1=nj-1
      xki=pi/(ni1-1)
      xkj=pi/(nj1-1)
      rdx2=2./(dx*dx)
      f02=f0*f0
      sfac=sqrt(8.) 
c
      do 1 kl=1,ktrunc
      call Ifttran (work,pvort(1,1,kl),work(1,1,2),ni,ni1,nj1)
      vfac=f02/edepth(kl)
      do 1 kj=1,jtrunc
      cosj=cos(kj*xkj)
      do 1 ki=1,itrunc
      cosi=cos(ki*xki)
      fac=sfac*sqrt(vfac+(1.-cosi*cosj)*rdx2)
      pvorts(ki,kj,kl)=work(ki,kj,1)/fac
    1 continue
c
      return
      end
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE PVp2spcA (pvorts,pvort,work,edepth,f0,dx,
     &                     ni,nj,itrunc,jtrunc,ktrunc)
      real pvorts(itrunc,jtrunc,ktrunc), pvort(ni,nj,ktrunc),
     &     work(ni,nj,2), edepth(ktrunc)
c
c  Compute adjoint of routine = PVp2spc
c  gridfac accounts for the absence of this factor in Iftproj
c
      pi=4.*atan(1.)
      ni1=ni-1
      nj1=nj-1
      nij=ni*nj
      gridfac=4./((ni-2)*(nj-2))
      xki=pi/(ni1-1)
      xkj=pi/(nj1-1)
      rdx2=2./(dx*dx)
      f02=f0*f0
      sfac=sqrt(8.)/gridfac 
c
      call Nrsetc (work,nij,0.)
      do 2 kl=1,ktrunc
      vfac=f02/edepth(kl)
      do 1 kj=1,jtrunc
      cosj=cos(kj*xkj)
      do 1 ki=1,itrunc
      cosi=cos(ki*xki)
      fac=sfac*sqrt(vfac+(1.-cosi*cosj)*rdx2)
      work(ki,kj,1)=pvorts(ki,kj,kl)/fac
    1 continue
      call Iftproj (pvort(1,1,kl),work,work(1,1,2),ni,ni1,nj1,.t.)
    2 continue
c
      return
      end

C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE Apotvort (pvort,vort,geop,edepth,f0,ni,nj,nk)
      real pvort(ni,nj,nk), vort(ni,nj,nk), geop(ni,nj,nk), edepth(nk)
c
c  This routine computes the adjoint of routine = Ipotvort
c
      ni2=ni-2
      nj2=nj-2
c
      do 1 k=1,nk
      fgh=-f0/edepth(k)
      do 1 j=2,nj2
      do 1 i=2,ni2
      vort(i,j,k)=pvort(i,j,k)
      geop(i,j,k)=fgh*pvort(i,j,k)
    1 continue
c
c  zero edges
      call cbdry0 (vort,ni,nj,nk,1)
      call cbdry0 (geop,ni,nj,nk,1)
c
      return
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c	The routine PVpv2fsd is a basically a copy of PVpv2fs
c	that is used in catx to dimension the fields. The only 
c	difference is that it returns the physical fields on 
c	the vector xfull (nfull). 15 december 1993.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE PVpv2fsd (u,v,t,ps,pvort,geop,psi,um,vm,
     &                    work1,work2,edepth,ematrixi,gmatrixp,
     &                    gmatrixr,pweight,pmeanb,clat,dx,
     &                    ni,nj,nk,itrunc,jtrunc,ktrunc,XFULL,NFULL) 
C
C  Compute u, v, T, ps from pvort in spectral space
C  Pvort is field of rotational modes.
C
      parameter (tpratio=0.)
      real pvort(itrunc,jtrunc,ktrunc)
     &,    geop(ni,nj,nk), psi(ni,nj,nk), um(ni,nj,nk)
     &,    work1(ni,nj,nk), work2(ni,nj,nk), vm(ni,nj,nk)
     &,    u(ni,nj,nk), v(ni,nj,nk), t(ni,nj,nk), ps(ni,nj)
     &,    edepth(nk), ematrixi(nk,nk), pweight(nk+2,3)
     &,    gmatrixp(nk), gmatrixr(nk,nk) 
c
	real xfull (nfull) 
c
      nij=ni*nj
      nijk=nij*nk
      pi4=atan(1.)
      f0=16.*pi4*sin(clat*pi4/45.)/(60.*60.*24)
      rf0=1./f0
      rpmeanb=1./pmeanb
      ni1=ni-1
      nj1=nj-1
c
c  compute geop in vertical mode space from pvort in spectral space 
      call Nrsetc (geop,nijk,0.)
      call PVpv2gp (geop,pvort,edepth,work1,f0,dx,
     &              ni,nj,itrunc,jtrunc,ktrunc)
c
c  compute psi in vertical mode space
      call Ncopyv (psi,geop,nijk)
      call Nmultv (psi,nijk,rf0)
c
c  compute change in ps (work1=geop change on sigma)
      call Nrsetc (ps,nij,0.)
      call Fchangep (ps,work2,ps,geop,pweight,ematrixi,
     &               gmatrixr,work1,tpratio,ni,nj,nk,ktrunc,
     &               work2,work2,.f.,.true.)
c
c  compute new t 
      call Nrsetc (work2,nijk,0.)
      call Nrsetc (t,nijk,0.)
      call Fchanget (t,work1,t,work2,geop,ps,gmatrixp,gmatrixr,
     &               work2,work2,0.,rpmeanb,ni,nj,nk,ktrunc,.f.)
c
c  compute u and v in vertical mode space (store in um, vm)
      call Nrsetc (work1,nijk,0.)  ! set divg=0.
      call Nrsetc (um,nijk,0.)  
      call Nrsetc (vm,nijk,0.)  
      call Ichanguv (um,vm,psi,work1,dx,ni,nj,ktrunc)
c
c  compute u and v in real space
      call Sproject (u,um,ematrixi,ni,nj,nk,ktrunc,nk,0,2,.t.)
      call Sproject (v,vm,ematrixi,ni,nj,nk,ktrunc,nk,0,2,.t.)
c
c	now the physical fields are available, put them on the 
c	vector xfull, compare also mstar1 and mstar4
c
      ni1=ni-1
      nj1=nj-1
      ni2=ni-2
      nj2=nj-2
      nxu=0
      nxv=nxu+ni2*nj2*nk
      nxt=nxv+ni2*nj2*nk
      nxp=nxt+(ni2-1)*(nj2-1)*nk
c  
c  ps 
c
      nx=nxp
      do 1 j=2,nj2
      do 1 i=2,ni2
      nx=nx+1
      xfull(nx)=ps(i,j)
    1 continue
c  
c  pt
c
      nx=nxt
      do 2 k=1,nk
      do 2 j=2,nj2
      do 2 i=2,ni2
      nx=nx+1
      xfull(nx)=t(i,j,k)
    2 continue
c  
c  pu and pv
c
      do 3 k=1,nk
      do 3 j=2,nj1
      do 3 i=2,ni1
      nxu=nxu+1
      nxv=nxv+1
      xfull(nxu)=u(i,j,k)
      xfull(nxv)=v(i,j,k)
    3 continue
c
      return
      end
C
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE QNmatrix (qmatrix,qinverse,qtotvar,tpmeanb,sigh,dsig,
     &            ptop,qsatvp,qs,qvar,eofiw,evarri,nqsatvp,nk,lprint)
      parameter (stdev=.2, rhfac=0.1)
      common /cnorm4/ itruncq,jtruncq,ktruncq,qnorm
      logical lprint
      real tpmeanb(nk+1)
      real qs(nk)
      real sigh(nk), dsig(nk), qsatvp(nqsatvp)
      real qvar(nk,nk)
      real evarri(nk,2),qmatrix(nk,nk),eofiw(nk,nk,2),qinverse(nk,nk)
c
      nkp1=nk+1
      nk2=nk*nk
c
c  compute saturation mixing ratio
      do k=1,nk
        call Nqsat (qs(k),tpmeanb(k),tpmeanb(nkp1),qsatvp,
     &              sigh(k),ptop,nqsatvp,1,1,1,1,1,0,0)
      enddo
c
c  compute q covariance matrix based on relative humidity
      sfac=0.5/(stdev**2)
      rhfac2=rhfac**2
      do k1=1,nk
      do k2=k1,nk
      sx=-sfac*(sigh(k1)-sigh(k2))*(sigh(k1)-sigh(k2))
      qvar(k1,k2)=rhfac2*qs(k1)*qs(k2)*exp(sx)
      enddo
      enddo
c
      do k1=2,nk
      do k2=1,k1-1
      qvar(k1,k2)=qvar(k2,k1)
      enddo
      enddo
c
c  pre and post multiply by diag (dsig**0.5)
      do k2=1,nk
      do k1=1,nk
      ds=sqrt(dsig(k1)*dsig(k2))
      qvar(k1,k2)=ds*qvar(k1,k2)
      enddo
      enddo
c
c  compute vertical EOFs
      call Rg (nk,nk,qvar,evarri,evarri(1,2),1,qmatrix,
     &         eofiw,eofiw(1,1,2),ier)
      call Sorder (qmatrix,evarri,eofiw,eofiw(1,1,2),nk)
      call Nrsetc (evarri(1,2),nk,1.)
      call Snorml (qmatrix,evarri(1,2),nk)
c
c  compute total variance (but no longer used)
      qtotvar=0.
      do k=1,nk
      qtotvar=qtotvar+evarri(k,1)
      enddo
c
c  rescale EOFs with dsig weight
      do k1=1,nk
      fac=1./sqrt(dsig(k1))
      do k2=1,nk
      qmatrix(k1,k2)=fac*qmatrix(k1,k2)
      enddo
      enddo
c
c  compute inverse 
      do k2=1,nk
      rfac=dsig(k2)
      do k1=1,nk
      qinverse(k1,k2)=rfac*qmatrix(k2,k1)
      enddo
      enddo
c
      if (lprint) then
      print *,' '
      print *,'QMATRIX (VERTICAL EOFS) '
      do k1=1,nk
      print 1,k1,evarri(k1,1),(qmatrix(k1,k2),k2=1,nk) 
      enddo
c
      print *,' '
      print *,'QINVERSE '
      do k1=1,nk
      print 1,k1,evarri(k1,1),(qinverse(k1,k2),k2=1,nk) 
      enddo
c
    1 format(i3,1p1e10.1,3x,10e10.1)      
c
      call Nrsetc (eofiw,nk2,0.)
      do k3=1,nk
      do k2=1,nk
      do k1=1,nk
      eofiw(k1,k3,1)=eofiw(k1,k3,1)+qmatrix(k1,k2)*qinverse(k2,k3)
      enddo
      enddo
      enddo
c
      print *,' '
      print *,'QMATRIX*QINVERSE '
      do k1=1,nk
      print 1,k1,evarri(k1,1),(eofiw(k1,k2,1),k2=1,nk) 
      enddo
      endif
c
      return
      end

C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE QNf2s (qspecs,q,qinverse,qvmode,work,qnorm,
     &                  ni,nj,nk,itrunc,jtrunc,ktrunc)
c
c  Compute q in vertical EOF, horizontal sine series, scaled from grid q
c
      real qspecs(itrunc,jtrunc,ktrunc)
     &,    q(ni,nj,nk), qvmode(ni,nj,ktrunc), work(ni,nj,2)
     &,    qinverse(nk,nk)
c
      qfac=qnorm
c
c  compute vertical EOF amplitudes on grid
      call Sproject (qvmode,q,qinverse,ni,nj,ktrunc,nk,nk,1,1,.t.)
c
c  compute horizontal spectra of vertical EOFs
      call QNqv2qs (qspecs,qvmode,work,qfac,
     &              ni,nj,itrunc,jtrunc,ktrunc)
c
      return
      end
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE QNqv2qs (qspecs,qvmode,work,qfac,
     &                    ni,nj,itrunc,jtrunc,ktrunc)
      real qspecs(itrunc,jtrunc,ktrunc), qvmode(ni,nj,ktrunc),
     &     work(ni,nj,2)
c
c  compute scaled spectral q from q in vertical EOF space on grid
c
      ni1=ni-1
      nj1=nj-1
c
      do kl=1,ktrunc
        call Ifttran (work,qvmode(1,1,kl),work(1,1,2),ni,ni1,nj1)
        do kj=1,jtrunc
          do ki=1,itrunc
            qspecs(ki,kj,kl)=qfac*work(ki,kj,1)
          enddo
        enddo
      enddo
c
      return
      end
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE QNs2f (qspecs,q,qmatrix,qvmode,work,qnorm,
     &                  ni,nj,nk,itrunc,jtrunc,ktrunc)
c
c  Compute q on grid from scaled q in vert. EOF, horiz. sine series
c
      real qspecs(itrunc,jtrunc,ktrunc)
     &,    q(ni,nj,nk), qvmode(ni,nj,ktrunc), work(ni,nj,2)
     &,    qmatrix(nk,nk)
c
      nij=ni*nj
      nijk=nij*nk
      nkmd=ktrunc
      nijkmd=nij*nkmd
      qfac=qnorm
c
c  compute vertical EOFs from their horizontal spectra
      call Nrsetc (qvmode,nijkmd,0.)
      call QNqs2qv (qspecs,qvmode,work,qfac,
     &              ni,nj,itrunc,jtrunc,ktrunc)
c
c  compute gridded fields from their vertical EOFS
      call Sproject (q,qvmode,qmatrix,ni,nj,nk,ktrunc,nk,1,1,.t.)
c
      return
      end
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE QNqs2qv (qspecs,qvmode,work,qfac,
     &                    ni,nj,itrunc,jtrunc,ktrunc)
      real qspecs(itrunc,jtrunc,ktrunc), qvmode(ni,nj,ktrunc),
     &     work(ni,nj,2)
c
c  compute spectral q from q in vertical EOF space on grid
c
      ni1=ni-1
      nj1=nj-1
      nn=ni*nj
c
      call Nrsetc (work,nn,0.)
      do kl=1,ktrunc
        do kj=1,jtrunc
          do ki=1,itrunc
            work(ki,kj,1)=qspecs(ki,kj,kl)/qfac
          enddo
        enddo
        call Iftproj (qvmode(1,1,kl),work,work(1,1,2),ni,ni1,nj1,.true.)
      enddo
c
      return
      end
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE QNf2sA (qspecs,q,qinversT,qvmode,work,qnorm,
     &                   ni,nj,nk,itrunc,jtrunc,ktrunc)
c
c  Adjoint of routine = QNf2s
c  (fac accounts for factor in Iftproj not in adjoint of Ifttran)
c
      real qspecs(itrunc,jtrunc,ktrunc)
     &,    q(ni,nj,nk), qvmode(ni,nj,ktrunc), work(ni,nj,2)
     &,    qinversT(nk,nk)
c
      nijk=ni*nj*nk
      fac=4./((ni-2)*(nj-2))
      qfac=fac*qnorm
c
      call Nrsetc (q,nijk,0.)
      call QNqv2qsA (qspecs,qvmode,work,qfac,
     &              ni,nj,itrunc,jtrunc,ktrunc)
      call Sproject (q,qvmode,qinversT,ni,nj,nk,ktrunc,nk,1,1,.f.)
c
      return
      end
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE QNqv2qsA (qspecs,qvmode,work,qfac,
     &                    ni,nj,itrunc,jtrunc,ktrunc)
      real qspecs(itrunc,jtrunc,ktrunc), qvmode(ni,nj,ktrunc),
     &     work(ni,nj,2)
c
c  Adjoint of routine = QNqv2qs
c
      ni1=ni-1
      nj1=nj-1
      nij=ni*nj
c
      call Nrsetc (work,nij,0.)
      do kl=1,ktrunc
        do kj=1,jtrunc
          do ki=1,itrunc
            work(ki,kj,1)=qfac*qspecs(ki,kj,kl)
          enddo
        enddo
        call Iftproj (qvmode(1,1,kl),work,work(1,1,2),ni,ni1,nj1,.true.)
      enddo
c
      return
      end
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE QNs2fA (qspecs,q,qmatrixT,qvmode,work,qnorm,
     &                  ni,nj,nk,itrunc,jtrunc,ktrunc)
c
c  Adjoint of routine = QNs2f 
c  (fac accounts for factor in Iftproj not in adjoint of Ifttran)
c
      real qspecs(itrunc,jtrunc,ktrunc)
     &,    q(ni,nj,nk), qvmode(ni,nj,ktrunc), work(ni,nj,2)
     &,    qmatrixT(nk,nk)
c
      nij=ni*nj
      nijk=nij*nk
      nkmd=ktrunc
      nijkmd=nij*nkmd
      fac=4./((ni-2)*(nj-2))
      qfac=fac*qnorm
c
      call Nrsetc (qvmode,nijkmd,0.)
      call Sproject (qvmode,q,qmatrixT,ni,nj,ktrunc,nk,nk,1,1,.f.)
      call QNqs2qvA (qspecs,qvmode,work,qfac,
     &              ni,nj,itrunc,jtrunc,ktrunc)
c
      return
      end
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE QNqs2qvA (qspecs,qvmode,work,qfac,
     &                    ni,nj,itrunc,jtrunc,ktrunc)
      real qspecs(itrunc,jtrunc,ktrunc), qvmode(ni,nj,ktrunc),
     &     work(ni,nj,2)
c
c  Adjoint of routine = QNqs2qv 
c
      ni1=ni-1
      nj1=nj-1
c
      do kl=1,ktrunc
        call Ifttran (work,qvmode(1,1,kl),work(1,1,2),ni,ni1,nj1)
        do kj=1,jtrunc
          do ki=1,itrunc
            qspecs(ki,kj,kl)=work(ki,kj,1)/qfac
          enddo
        enddo
      enddo
c
      return
      end
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE QNdimv (xfull,q,ni,nj,nk,nfull)
c
c  Put in q format of dimensioned SVs
c
      real q(ni,nj,nk), xfull(nfull)
c
      ni2 = ni - 2
      nj2 = nj - 2
c
      nxu = 0
      nxv = nxu + ni2 * nj2 * nk
      nxt = nxv + ni2 * nj2 * nk
      nxp = nxt + (ni2-1) * (nj2-1) * nk
      nq  = nxp + (ni2-1) * (nj2-1)
c
      do k=1,nk
        do j=2,nj2
          do i=2,ni2
            nq=nq+1
            xfull(nq)=q(i,j,k)
          enddo
        enddo
      enddo
c
      return
      end









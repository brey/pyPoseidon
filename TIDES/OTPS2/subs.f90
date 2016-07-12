!=====================================================================
! AUTHORS:
!  Gary Egbert & Lana Erofeeva
!  College of Atmospheric and Oceanic Sciences
!  104 COAS Admin. Bldg.
!  Oregon State University
!  Corvallis, OR 97331-5503
!  
!  E-mail:  egbert@coas.oregonstate.edu                                      
!  Fax:     (541) 737-2064
!  Ph.:     (541) 737-2947                                        
!  http://volkov.oce.orst.edu/tides/
!
! COPYRIGHT: OREGON STATE UNIVERSITY, 2012
! (see the file COPYRIGHT for lisence agreement)
!=====================================================================
! Set of subroutines used by extract_HC and predict_tide
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine ptide(z1,cid,ncon,ind,lat,time_mjd,ntime, &
                         interp,zpred)
        implicit none
        include 'constit.h'
        integer ncon,ntime,k,ind(ncon),ierr
        complex z1(ncon)
        real lat,zpred(ntime)
        real*8 time_mjd(ntime),time,dh
        character*4 cid(ncon)
        integer SecondsPerDay
        real height
        real ww(17,8)
!        complex z1(ncon)
        complex, allocatable:: A(:)
! nodal arguments shoud be real*8
        real*8 pu(ncmx),pf(ncmx),dlat,t_mjd
        logical, intent(inout):: interp
!
        !if(interp)call mkw(interp,ind,ncon,ww)
        allocate(A(ncon))
	SecondsPerDay=86400
        dlat=lat
        ierr=0
        do k=1,ntime
	 call nodal(time_mjd(k),lat,pu,pf) ! changed dlat to lat
! to use phase shifts from constit.h time should be in seconds
! relatively Jan 1 1992 (48622mjd)       
         time=(time_mjd(k)-dble(48622))*dble(SecondsPerDay)
! .true. means NO solid Earth correction will be applied in make_a
         call make_a(.false.,ind,ncon,time,pu,pf,ww,A,.true.)
         zpred(k)=height(A,z1,ncon)
         if(interp)then
          dh=0.
          call infer_minor(z1,cid,ncon,time_mjd(k),dh,ierr)
          if(ierr.eq.-1)then
           write(*,*)'Not enough constituents for inference of'
           write(*,*)'minor constituents: IGNORED'
           interp=.false.
          endif
          zpred(k)=zpred(k)+dh ! add minor constituents
         endif
        enddo
        deallocate(A)
        return
        end
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Lana Erofeeva, June 2004
! Based on Richard Ray code perth2
! Return correction for the 16 minor constituents
! listed in subroutine
!         infer minor tides
!         -----------------
         subroutine infer_minor(zmaj,cid,ncon,time,dh,ierr)
         implicit none
         integer ncon,i,j,ni,k
         complex zmaj(ncon)       ! HC for GIVEN consituents
         character*4 cid(ncon)    ! GIVEN constituents
         real*8 time              ! time, mjd
         real*8 dh                ! OUTPUT: correction at given time
                                  ! for 16 minor constituents
         integer ierr             ! -1 if not enough constuents
                                  !    for inference 
         complex*16 zmin(18)
         complex*16 z8(8)
         real*8 hour,t1,t2,shpn(4),s,h,p,omega
         real*8 sinn,cosn,sin2n,cos2n
         real*8 u(18),f(18),arg(18)

         real*8, parameter:: pi=3.141592654D0
         real*8, parameter:: rad=pi/180.D0
         real*8, parameter:: PP=282.8D0 
         EQUIVALENCE (SHPN(1),S),(SHPN(2),H),(SHPN(3),P),(SHPN(4),omega) 

         character*4 cid8(8)      ! in order to correspnd RR coefficients
         data cid8/'q1  ','o1  ','p1  ','k1  ', &
                   'n2  ','m2  ','s2  ','k2  '/
! re-order to correspond to cid8
         ierr=0
         z8=0.
         ni=0
         do i=1,8
          do j=1,ncon
           if(cid(j).eq.cid8(i))then
              z8(i)=zmaj(j)
              if(i.ne.3.and.i.ne.8)ni=ni+1
           endif 
          enddo
         enddo

         if(ni.lt.6)then
          ierr=-1 ! not enough constituents for inference
          return 
         endif
!            
          zmin(1)  = 0.263 *z8(1) - 0.0252*z8(2)  !2Q1
          zmin(2) = 0.297 *z8(1) - 0.0264*z8(2)   !sigma1
          zmin(3) = 0.164 *z8(1) + 0.0048*z8(2)   !rho1 +
          zmin(4) = 0.0140*z8(2) + 0.0101*z8(4)   !M1
          zmin(5) = 0.0389*z8(2) + 0.0282*z8(4)   !M1
          zmin(6) = 0.0064*z8(2) + 0.0060*z8(4)   !chi1
          zmin(7) = 0.0030*z8(2) + 0.0171*z8(4)   !pi1
          zmin(8) =-0.0015*z8(2) + 0.0152*z8(4)   !phi1
          zmin(9) =-0.0065*z8(2) + 0.0155*z8(4)   !theta1
          zmin(10) =-0.0389*z8(2) + 0.0836*z8(4)  !J1 +
          zmin(11) =-0.0431*z8(2) + 0.0613*z8(4)  !OO1 +
          zmin(12) = 0.264 *z8(5) - 0.0253*z8(6)  !2N2 +
          zmin(13) = 0.298 *z8(5) - 0.0264*z8(6)  !mu2 +
          zmin(14) = 0.165 *z8(5) + 0.00487*z8(6) !nu2 +
          zmin(15) = 0.0040*z8(6) + 0.0074*z8(7)  !lambda2
          zmin(16) = 0.0131*z8(6) + 0.0326*z8(7)  !L2 +
          zmin(17) = 0.0033*z8(6) + 0.0082*z8(7)  !L2 +
          zmin(18) = 0.0585*z8(7)                 !t2 + 
!
          hour = (time - int(time))*24.D0
          t1 = 15.*hour
          t2 = 30.*hour
          call astrol( time, SHPN )
!
          arg(1) = t1 - 4.*S + H + 2.*P - 90.     ! 2Q1
          arg(2) = t1 - 4.*S + 3.*H - 90.         ! sigma1
          arg(3) = t1 - 3.*S + 3.*H - P - 90.     ! rho1
          arg(4) = t1 - S + H - P + 90.           ! M1
          arg(5) = t1 - S + H + P + 90.           ! M1
          arg(6) = t1 - S + 3.*H - P + 90.        ! chi1
          arg(7) = t1 - 2.*H + PP - 90.           ! pi1
          arg(8) = t1 + 3.*H + 90.                ! phi1
          arg(9) = t1 + S - H + P + 90.           ! theta1
          arg(10) = t1 + S + H - P + 90.          ! J1
          arg(11) = t1 + 2.*S + H + 90.           ! OO1
          arg(12) = t2 - 4.*S + 2.*H + 2.*P       ! 2N2
          arg(13) = t2 - 4.*S + 4.*H              ! mu2
          arg(14) = t2 - 3.*S + 4.*H - P          ! nu2
          arg(15) = t2 - S + P + 180.D0           ! lambda2
          arg(16) = t2 - S + 2.*H - P + 180.D0    ! L2
          arg(17) = t2 - S + 2.*H + P             ! L2
          arg(18) = t2 - H + PP                   ! t2
!
!     determine nodal corrections f and u
         sinn = SIN(omega*rad)
         cosn = COS(omega*rad)
         sin2n = SIN(2.*omega*rad)
         cos2n = COS(2.*omega*rad)
         F = 1.
         f(1) = sqrt((1.0 + 0.189*cosn - 0.0058*cos2n)**2 + &
                      (0.189*sinn - 0.0058*sin2n)**2)
         f(2) = f(1)
         f(3) = f(1)
         f(4) = sqrt((1.0 + 0.185*cosn)**2 + (0.185*sinn)**2)
         f(5) = sqrt((1.0 + 0.201*cosn)**2 + (0.201*sinn)**2)
         f(6) = sqrt((1.0 + 0.221*cosn)**2 + (0.221*sinn)**2)
         f(10) = sqrt((1.0 + 0.198*cosn)**2 + (0.198*sinn)**2)
         f(11) = sqrt((1.0 + 0.640*cosn + 0.134*cos2n)**2 + &
                      (0.640*sinn + 0.134*sin2n)**2 )
         f(12) = sqrt((1.0 - 0.0373*cosn)**2 + (0.0373*sinn)**2)
         f(13) = f(12)
         f(14) = f(12)
         f(16) = f(12)
         f(17) = sqrt((1.0 + 0.441*cosn)**2 + (0.441*sinn)**2)
!
         U = 0.
         u(1) = atan2(0.189*sinn - 0.0058*sin2n, &
                      1.0 + 0.189*cosn - 0.0058*sin2n)/rad
         u(2) = u(1)
         u(3) = u(1)
         u(4) = atan2( 0.185*sinn, 1.0 + 0.185*cosn)/rad
         u(5) = atan2(-0.201*sinn, 1.0 + 0.201*cosn)/rad
         u(6) = atan2(-0.221*sinn, 1.0 + 0.221*cosn)/rad
         u(10) = atan2(-0.198*sinn, 1.0 + 0.198*cosn)/rad
         u(11) = atan2(-0.640*sinn - 0.134*sin2n, &
                       1.0 + 0.640*cosn + 0.134*cos2n)/rad
         u(12) = atan2(-0.0373*sinn, 1.0 - 0.0373*cosn)/rad
         u(13) = u(12)
         u(14) = u(12)
         u(16) = u(12)
         u(17) = atan2(-0.441*sinn, 1.0 + 0.441*cosn)/rad

!     sum over all tides
!     ------------------
        dh = 0.D0
        do i=1,18
         dh = dh + dreal(zmin(i))*f(i)*cos((arg(i)+u(i))*rad)- &
                   dimag(zmin(i))*f(i)*sin((arg(i)+u(i))*rad)
        enddo
        return
        end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! define constituent indices of c_id in c_id_mod
! constituent names are assumed to be all in LOWER case
      subroutine def_con_ind(c_id,ncon,c_id_mod,nc,cind)
      implicit none
      include 'constit.h'
      character*4 c_id(ncmx),c_id_mod(ncmx)
      integer nc,ncon,cind(ncon),ic1,ic2
!
      cind=0     
      do ic1=1,ncon
       do ic2=1,nc
        if(c_id(ic1).eq.c_id_mod(ic2))cind(ic1)=ic2
       enddo
       if(cind(ic1)==0)then
        !write(*,*)'WARNING: No constituent ',c_id(ic1),'in the model'
       endif
      enddo
!
      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rd_inp(modname,lltname,zuv,c_id,ncon, &
                        APRI,Geo,outname,interp)
      implicit none
      include 'constit.h'
      character*80 modname,lltname,outname,tmp,rmCom
      character*1 zuv
      logical APRI,Geo,interp 
      character*4 c_id(ncmx),con
      integer k,ic,ncon
!
      interp=.false.
      read(*,'(a80)',err=1)modname
      modname=rmCom(modname)
      if(trim(modname).eq.'')then
       write(*,*)'Model name can not be empty.'
       write(*,*)'A tidal model is needed for OTPS to work.'
       stop
      endif
      read(*,'(a80)',err=1)lltname
      lltname=rmCom(lltname)
      read(*,'(a1)',err=1)zuv
      read(*,'(a80)',err=1)tmp
      tmp=rmCom(tmp)
      ic=0
      k=1
      do while(k.gt.0)
       ncon=ic
       k=index(tmp,',')
       ic=ic+1
       c_id(ic)=tmp(1:max(1,k-1))
       tmp=tmp(k+1:80)
      enddo
      if(trim(tmp).ne.'')ncon=ic
      c_id(ic)=tmp(1:80)
      read(*,'(a80)')tmp
      tmp=rmCom(tmp)
      APRI=.false.
      if(trim(tmp).eq.'AP')APRI=.true.
      geo=.false.
      read(*,'(a80)',err=1)tmp ! geo/oce
      tmp=rmCom(tmp)
      if(tmp(1:3).eq.'geo')geo=.true.
      read(*,*,err=1)k
      if(k.eq.1)interp=.true.
      read(*,'(a80)',err=1)outname
      outname=rmCom(outname)
      return
1     write(*,*)'Input file of WRONG FORMAT: see setup.inp'
      write(*,*)'for right format. Recovering setup.inp ... done'
      open(unit=111,file='dumsetup',status='old')
      open(unit=112,file='setup.inp',status='unknown')
3     read (111,'(a80)',end=2)tmp
      write(112,'(a)')trim(tmp)
      go to 3
2     close(111)
      close(112)      
      write(*,*)'Please edit ''setup.inp'' for your settings'
      write(*,*)'and re-run extract_HC or predict_tide'
      stop
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      character*80 function rmCom(str) ! removes everything after '!' 
                                       ! in a string 
      character str*(*)
      integer k 
      k=index(str,'!')
      rmCom=str
      if(k>0)rmCom=str(1:k-1)
      return
      end
!____________________________________________________________________
!
      function ipshft(i,ish,n)
!      periodic shift maps i to i+ish  , mod n;
!        (always between 1 and n;  never 0 )
      integer i,ipshft,n,ish
      ipshft = mod(i+ish+n-1,n)+1
      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! make_a:: computes A matrix elements for one data point 
	subroutine make_a(interp,ind,nc,time,pu,pf,w,A,l_sal) 
        implicit none
	include 'constit.h'
        integer, parameter:: ncon = 17 ! for case of using weights.h
!
        logical interp,l_sal
        real omega(ncmx),phase(ncmx)
	integer ind(nc)
	real w(ncon,8)
	real*8 pu(ncmx),pf(ncmx)
        real*8 time
        complex A(nc),c(ncmx)
        integer i,j,nc
! If l_sal=.true. - NO solid Earth correction is applied
! USING beta_SE coefficients
!
        if(.not.interp)then
          do j=1,ncmx
           omega(j)=omega_d(j)
           phase(j)=phase_mkB(j)
          enddo
  	  do j=1,nc
            i=ind(j)
            if(i.ne.0)then
 	    c(j) = cmplx( pf(i)*cos(omega(i)*time+phase(i)+pu(i)), &
      	                  pf(i)*sin(omega(i)*time+phase(i)+pu(i)))
            endif
	  enddo
!  remove solid Earth tide
          if(.not.l_sal)then
           do j=1,nc
            A(j)=0. 
            if(ind(j).ne.0)A(j)=c(j)*beta_SE(ind(j))
           enddo
          else
           do j=1,nc 
            A(j)=c(j)*1.
           enddo
          endif
        else
! this is the case when w from weights.h is used -> see
! comments to ../include/weights.h
!----------------------------------------------------------c        
          omega(1:ncon)=omega_d(1:ncon)
          phase(1:ncon)=phase_mkB(1:ncon)
!
  	  do i=1,ncon 
 	    c(i) = cmplx( pf(i)*cos(omega(i)*time+phase(i)+pu(i)), &
      	                  pf(i)*sin(omega(i)*time+phase(i)+pu(i)))
	  enddo
          A=cmplx(0,0)
!
! ind(j)=0 means the constituent is excluded
          do i=1,ncon
            do j=1,nc 
             if(ind(j).ne.0)A(j)=A(j)+c(i)*beta_SE(i)*w(i,ind(j))
            enddo
          enddo
!---------------------------------------------------------c

        endif      
!
        return
        end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      

	real function height(A,P,NC)
! returns height from model array of complex constituents
        complex A(nc),P(nc)
        real sum
        integer i,nc
        if(nc.eq.0)then
          height=0.
          return
        endif
        sum=0.
! height(i)=sum_of_real(A(i)*P(i))
        do i=1,nc
         sum=sum+real(P(i))*real(A(i))-aimag(P(i))*aimag(A(i))
        enddo
        height=sum
        return
        end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rd_inp_JS(setup_file_unit,modname,zuv,c_id,ncon, &
                        APRI,Geo,interp)
!
      implicit none
      include 'constit.h'
      integer setup_file_unit
      character*80 modname,outname,tmp,rmCom
      character*80 lltname
      character*1 zuv
      logical APRI,Geo,interp 
      character*4 c_id(ncmx),con
      integer k,ic,ncon
!
      interp=.false.
      read(setup_file_unit,'(a80)',err=1)modname
      modname=rmCom(modname)
      read(setup_file_unit,'(a80)',err=1)lltname
      lltname=rmCom(lltname)
      read(setup_file_unit,'(a1)',err=1)zuv
      read(setup_file_unit,'(a80)',err=1)tmp
      tmp=rmCom(tmp)
      write(*,*)tmp
      ic=0
      k=1
      do while(k.gt.0)
       ncon=ic
       k=index(tmp,',')
       ic=ic+1
       c_id(ic)=tmp(1:k-1)
       tmp=tmp(k+1:80)
      enddo
      if(trim(tmp).ne.'')ncon=ic
      c_id(ic)=tmp(1:80)
      read(setup_file_unit,'(a80)')tmp
      tmp=rmCom(tmp)
      APRI=.false.
      if(trim(tmp).eq.'AP')APRI=.true.
      geo=.false.
      read(setup_file_unit,'(a80)',err=1)tmp ! geo/oce
      tmp=rmCom(tmp)
      if(tmp(1:3).eq.'geo')geo=.true.
      read(setup_file_unit,*,err=1)k
      if(k.eq.1)interp=.true.
      read(setup_file_unit,'(a80)',err=1)outname
      outname=rmCom(outname)
      return
1     write(*,*)'Input file of WRONG FORMAT: see setup.inp'
      write(*,*)'for right format. Recovering setup.inp ... done'
      open(unit=111,file='dumsetup',status='old')
      open(unit=112,file='setup.inp',status='unknown')
3     read (111,'(a80)',end=2)tmp
      write(112,'(a)')trim(tmp)
      go to 3
2     close(111)
      close(112)      
      write(*,*)'Please edit ''setup.inp'' for your settings'
      write(*,*)'and re-run tide_js'
      stop
      end
!  ARGUMENTS and ASTROL subroutines SUPPLIED by RICHARD RAY, March 1999
!  see arguments.txt: corrections Apr 2, 2003
!  attached to OTIS by Lana Erofeeva (subroutine nodal.f)
!  
!  March 2011 - replaced old subroutine "arguments" with
!               new Richard Ray's "fu_nodal"
!               corresponding adjustments made to nodal
!
	subroutine nodal(dtime,lat,pu,pf)
        implicit none
        include 'constit.h'
!
        real*8 rad
        parameter(rad=3.141592654358979/180.D0)
        integer jndex(ncmx),i
        real*4 t1,lat
        real*8 latitude,pu(ncmx),pf(ncmx)
        real*8 shpn(5),s,h,p,omega,dtime,pprime
        equivalence (shpn(1),s),(shpn(2),h),(shpn(3),p),(shpn(4),omega), &
                    (shpn(5),pprime)
!
        latitude=lat
!  dtime - should be MJD, t1 is given in secs relatively jan 1 1985 (altimetry)
        !!! dtime=t1/tsec+2556.+46066. !!!! already in mjd
        !call astrol(dtime,shpn)
        call astro5(dtime,shpn)
        do i=1,ncmx
         call fu_nodal( omega, p, constid(i), 1, pf(i), pu(i) ) 
         pu(i)=pu(i)*rad
         !write(*,*)constid(i),pf(i),pu(i)
        enddo
        return
        end 
!======================================================================
      SUBROUTINE ASTROL( TIME, SHPN )
!
!  Computes the basic astronomical mean longitudes  s, h, p, N.
!  Note N is not N', i.e. N is decreasing with time.
!  These formulae are for the period 1990 - 2010, and were derived
!  by David Cartwright (personal comm., Nov. 1990).
!  TIME is UTC in decimal MJD.
!  All longitudes returned in degrees.
!  R. D. Ray    Dec. 1990
!
!  Non-vectorized version.
!
!       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      double precision circle,shpn,t,time
      DIMENSION  SHPN(4)
      PARAMETER  (CIRCLE=360.0D0)
!
      T = TIME - 51544.4993D0
!
!     mean longitude of moon
!     ----------------------
      SHPN(1) = 218.3164D0 + 13.17639648D0 * T
!
!     mean longitude of sun
!     ---------------------
      SHPN(2) = 280.4661D0 +  0.98564736D0 * T
!
!     mean longitude of lunar perigee
!     -------------------------------
      SHPN(3) =  83.3535D0 +  0.11140353D0 * T
!
!     mean longitude of ascending lunar node
!     --------------------------------------
      SHPN(4) = 125.0445D0 -  0.05295377D0 * T

      SHPN(1) = MOD(SHPN(1),CIRCLE)
      SHPN(2) = MOD(SHPN(2),CIRCLE)
      SHPN(3) = MOD(SHPN(3),CIRCLE)
      SHPN(4) = MOD(SHPN(4),CIRCLE)

      IF (SHPN(1).LT.0.D0) SHPN(1) = SHPN(1) + CIRCLE
      IF (SHPN(2).LT.0.D0) SHPN(2) = SHPN(2) + CIRCLE
      IF (SHPN(3).LT.0.D0) SHPN(3) = SHPN(3) + CIRCLE
      IF (SHPN(4).LT.0.D0) SHPN(4) = SHPN(4) + CIRCLE
      RETURN
      END
!-------------------------------------------------------------------
      recursive subroutine fu_nodal( omega, p, constituent, n, f, u )
!
!  Computes tidal nodal (& perigee) corrections "f" and "u" 
!  for n tidal constituents,
!  given the mean longitudes p and omega.
!
!  Language: Fortran 90.
!
!-------------------------------------------------------
!  Input:
!     omega - mean longitude of lunar node, in degrees.
!     p     - mean longitude of lunar perigee.
!     consituent - array of constituent names (char*4)
!     n  - length of arrays {f,u,constituent}.
!  Output:
!     f  - modulation factor for given tide(s).
!     u  - phase correction for given tide, in degrees.
!-------------------------------------------------------
!
!  Note: omega = -N' (the 5th Doodson variable), decreasing in time.
!  This is not a very efficient routine, but it needn't be
!  called very often, so is ok.
!
!  R Ray  21 August 2003
!  Revised 9/15/08.
!  Revised 3/28/11 - fixed Mt, MSt.
!
!
      implicit none
      integer          n
      double precision time, f(*), u(*)
      character*4      constituent(*)
!
      double precision s,h,p,omega,pp
      double precision term1,term2,ftmp(3),utmp(3)
      character*4      ctmp(3)
      integer          i
      real         cosn,cos2n,sinn,sin2n,sin3n,sinpn,cospn
      real         sinp,cosp,sin2p,cos2p
      real         two,three
      parameter   (two=2.0,three=3.0)

      double precision pi, rad
      parameter       (pi=3.141592654d0, rad=pi/180.D0)

!
!     Various trig functions of astronomical longitudes
!     -------------------------------------------------
      sinn = sin(omega*rad)
      cosn = cos(omega*rad)
      sin2n = sin(two*omega*rad)
      cos2n = cos(two*omega*rad)
      sin3n = sin(three*omega*rad)
      sinp  = sin(p*rad)
      cosp  = cos(p*rad)
      sin2p = sin(two*p*rad)
      cos2p = cos(two*p*rad)
      sinpn = sin((p-omega)*rad)
      cospn = cos((p-omega)*rad)
!
!     Compute standard nodal corrections f and u 
!     ------------------------------------------
      do i=1,n

         select case ( constituent(i) )

            case ( 'mm  ','msm ' )
               term1 = -.0534*sin2p - .0219*sin((two*p-omega)*rad)
               term2 = 1.0 - .1308*cosn - .0534*cos2p  &
                       - .0219*cos((two*p-omega)*rad)

            case ( 'mf  ','msq ','msp ','mq  ' )
               term1 = -0.04324*sin2p - 0.41465*sinn - 0.03873*sin2n
               term2 = 1.0 +0.04324*cos2p +0.41465*cosn +0.03873*cos2n

            case ( 'msf ' )     !  This is linear tide, not compound.
               term1 = 0.137*sinn
               term2 = 1.0

            case ( 'mt  ' )
               term1 = -0.018*sin2p - 0.4145*sinn - 0.040*sin2n
               term2 = 1.0 +0.018*cos2p +0.4145*cosn +0.040*cos2n

            case ( 'mst ' )
               term1 = -0.380*sin2p - 0.413*sinn - 0.037*sin2n
               term2 = 1.0 +0.380*cos2p +0.413*cosn +0.037*cos2n

            case ( 'o1  ','so3 ' )
               term1 = 0.1886*sinn - 0.0058*sin2n - 0.0065*sin2p
               term2 = 1.0 + 0.1886*cosn - 0.0058*cos2n - 0.0065*cos2p

            case ( '2q1 ','q1  ','rho1','sig1' )
               term1 = 0.1886*sinn 
               term2 = 1.0 + 0.1886*cosn 
            case ( 'tau1' )
               term1 = 0.219*sinn 
               term2 = 1.0 - 0.219*cosn 

            case ( 'm1  ' )        ! This assumes m1 argument includes p.
               term1 = -0.2294*sinn - 0.3594*sin2p &
                         - 0.0664*sin((two*p-omega)*rad)
               term2 = 1.0 + 0.1722*cosn + 0.3594*cos2p &
                         + 0.0664*cos((two*p-omega)*rad)
            case ( 'chi1' )
               term1 = -0.221*sinn 
               term2 = 1.0 + 0.221*cosn 
            case ( 'p1  ' )
               term1 = -0.0112*sinn 
               term2 = 1.0 - 0.0112*cosn 
            case ( 'k1  ','sk3 ' )
               term1 = -0.1554*sinn + 0.0031*sin2n
               term2 = 1.0 + 0.1158*cosn - 0.0028*cos2n
            case ( 'j1  ','the1' )
               term1 = -0.227*sinn
               term2 = 1.0 + 0.169*cosn
            case ( 'oo1 ' )
               term1 = -0.640*sinn - 0.134*sin2n - 0.150*sin2p
               term2 = 1.0 + 0.640*cosn + 0.134*cos2n + 0.150*cos2p

            case ( 'm2  ','2n2 ','mu2 ','n2  ','nu2 ','lam2','2sm2', &
                   'ms4 ','eps2' )
               term1 = -0.03731*sinn + 0.00052*sin2n
               term2 = 1.0 - 0.03731*cosn + 0.00052*cos2n

            case ( 'l2  ' )
               term1 = -0.250*sin2p - 0.110*sin((two*p-omega)*rad) &
                         - 0.037*sinn
               term2 = 1.0 - 0.250*cos2p - 0.110*cos((two*p-omega)*rad) &
                         - 0.037*cosn
            case ( 'k2  ' )
               term1 = -0.3108*sinn - 0.0324*sin2n
               term2 = 1.0 + 0.2853*cosn + 0.0324*cos2n

            case ( 'eta2' )
               term1 = -0.436*sinn
               term2 = 1.0 + 0.436*cosn

            case ( 'm3  ','o3  ','F3  ' )   ! Linear 3rd-deg terms
               term1 = -0.056*sinn
               term2 = 1.0 - 0.056*cosn
            case ( 'j3  ' )
               term1 = -0.43*sinn 
               term2 = 1.0 + 0.43*cosn 

            case default
               term1 = 0.
               term2 = 1.

         end select

         f(i) = sqrt(term1**2 + term2**2)
         u(i) = atan( term1/term2 )/rad


!        Following tides are all compound &  use recursion
!        -------------------------------------------------

         select case ( constituent(i) )

            case ( 'm4  ','mn4 ','mns2' )
               ctmp(1)='m2  '
               call fu_nodal( omega, p, ctmp, 1, ftmp, utmp )
               f(i) = ftmp(1)**2
               u(i) = 2.0*utmp(1)

            case ( 'mk3 ' )
               ctmp(1:2) = (/'m2  ','k1  '/)
               call fu_nodal( omega, p, ctmp, 2, ftmp, utmp )
               f(i) = ftmp(1)*ftmp(2)
               u(i) = utmp(1) + utmp(2)

            case ( 'mk4 ' )
               ctmp(1:2) = (/'m2  ','k2  '/)
               call fu_nodal( omega, p, ctmp, 2, ftmp, utmp )
               f(i) = ftmp(1)*ftmp(2)
               u(i) = utmp(1) + utmp(2)

            case ( 'm6  ' )
               ctmp(1)='m2  '
               call fu_nodal( omega, p, ctmp, 1, ftmp, utmp )
               f(i) = ftmp(1)**3
               u(i) = 3.0*utmp(1)

            case ( 'm8  ' )
               ctmp(1)='m2  '
               call fu_nodal( omega, p, ctmp, 1, ftmp, utmp )
               f(i) = ftmp(1)**4
               u(i) = 4.0*utmp(1)

            case ( 'mfDW' )   !  special test of old Doodson-Warburg formula
               f(i) = 1.043 + 0.414*cosn
               u(i) = -23.7*sinn + 2.7*sin2n - 0.4*sin3n 

         end select

      end do

      return
      end
!----------------------------------------------------------------------------
      SUBROUTINE ASTRO5( TIME, SHPNP )
!
!---------------------------------------------------------------------
!  Computes the 5 basic astronomical mean longitudes  s, h, p, N, p'.
!
!  Note N is not N', i.e. N is decreasing with time.
!
!  TIME is UTC in decimal Modified Julian Day (MJD).
!  All longitudes returned in degrees.
!
!  R. D. Ray, NASA/GSFC   August 2003
!
!  Most of the formulae for mean longitudes are extracted from 
!  Jean Meeus, Astronomical Algorithms, 2nd ed., 1998.  
!  Page numbers below refer to this book.
!
!  Note: This routine uses TIME in UT and does not distinguish between
!    the subtle differences of UTC, UT1, etc.  This is more than adequate
!   for the calculation of these arguments, especially in tidal studies.
!
!  Revised 4/15/10 - minor change to modulo function eliminates some extra code.
!---------------------------------------------------------------------
!
      IMPLICIT NONE
      DOUBLE PRECISION TIME, SHPNP(5)

      DOUBLE PRECISION TJD,T,CIRCLE,DEL,TJLAST,D
      PARAMETER       (CIRCLE=360.0D0)
      DOUBLE PRECISION DELTAT
      EXTERNAL         DELTAT
      SAVE             DEL,TJLAST
      DATA TJLAST/-1.d0/

!     Convert to Julian Day and to Ephemeris Time
!     -------------------------------------------
      TJD = TIME + 2400000.5D0
      IF (ABS(TJD-TJLAST).GT.100.D0) THEN
         DEL = DELTAT( TJD )/86400.D0
         TJLAST = TJD
      ENDIF
      TJD = TJD + DEL

!     Compute time argument in centuries relative to J2000
!     ----------------------------------------------------
      T = ( TJD - 2451545.d0 )/36525.d0
!
!
!     mean longitude of moon (p.338)
!     ------------------------------
      SHPNP(1) = (((-1.53388d-8*T + 1.855835d-6)*T - 1.5786d-3)*T + &
                  481267.88123421d0)*T + 218.3164477d0
!
!     mean elongation of moon (p.338)
!     -------------------------------
      D = (((-8.8445d-9*T + 1.83195d-6)*T - 1.8819d-3)*T + &
                445267.1114034d0)*T + 297.8501921d0
!
!     mean longitude of sun
!     ---------------------
      SHPNP(2) = SHPNP(1) - D
!
!     mean longitude of lunar perigee (p.343)
!     ---------------------------------------
      SHPNP(3) = ((-1.249172d-5*T - 1.032d-2)*T + 4069.0137287d0)*T + &
                  83.3532465d0
!
!     mean longitude of ascending lunar node (p.144)
!     ----------------------------------------------
      SHPNP(4) = ((2.22222d-6*T + 2.0708d-3)*T - 1934.136261d0)*T + &
                  125.04452d0
!
!     mean longitude of solar perigee (Simon et al., 1994)
!     ----------------------------------------------------
      SHPNP(5) = 282.94d0 + 1.7192d0 * T

      SHPNP = MODULO( SHPNP, CIRCLE )
      RETURN
      END
!-------------------------------------------------------------------------!
      double precision function deltat( tjd )

!..author: F X Timmes, University of Chicago

!..Slightly revised by R Ray, GSFC, Aug 2003.  Also updated table.
!  Updated periodically the tables.  -RDR

!..this routine computes the difference between 
!..universal time and dynamical time.

!..input:
!..tjd  = UT julian day number

!..output:
!  deltat = ET - UT in seconds


      implicit none
      save

!..declare the pass
      double precision tjd,tjde,secdif


!..for storing the table
      integer          tabstart,tabend,tabsize
      parameter       (tabstart = 1620, &
                       tabend   = 2015, &
                       tabsize  = tabend - tabstart + 1)
      double precision dt(tabsize),years(tabsize)


!..local variables
      integer          i,j,iy,im,id,iflag,mp,iat
      parameter        (mp = 2)
      double precision hr,b,xx,dy

!..parameter mp above determines the order of the interpolatant
!..when interpolating in the table. mp=2=linear mp=3=quadratic and so on.
!..don't be too greedy about the order of the interpolant since
!..for many years the data is flat, and trying to fit anything other
!..that a line will produce unwanted oscillations. thus i choose mp=2.
   


!..these tables of observed and extrapolated data are
!..are from the us naval observatory ftp://maia.usno.navy.mil/ser7/
!..years 1620 to 1710
      data  (dt(i), i=1,90) / &
          124.00d0, 119.00d0, 115.00d0, 110.00d0, 106.00d0, 102.00d0, &
           98.00d0,  95.00d0,  91.00d0,  88.00d0,  85.00d0,  82.00d0, &
           79.00d0,  77.00d0,  74.00d0,  72.00d0,  70.00d0,  67.00d0, & 
           65.00d0,  63.00d0,  62.00d0,  60.00d0,  58.00d0,  57.00d0, &
           55.00d0,  54.00d0,  53.00d0,  51.00d0,  50.00d0,  49.00d0, &
           48.00d0,  47.00d0,  46.00d0,  45.00d0,  44.00d0,  43.00d0, &
           42.00d0,  41.00d0,  40.00d0,  38.00d0,  37.00d0,  36.00d0, &
           35.00d0,  34.00d0,  33.00d0,  32.00d0,  31.00d0,  30.00d0, &
           28.00d0,  27.00d0,  26.00d0,  25.00d0,  24.00d0,  23.00d0, &
           22.00d0,  21.00d0,  20.00d0,  19.00d0,  18.00d0,  17.00d0, &
           16.00d0,  15.00d0,  14.00d0,  14.00d0,  13.00d0,  12.00d0, &
           12.00d0,  11.00d0,  11.00d0,  10.00d0,  10.00d0,  10.00d0, &
            9.00d0,   9.00d0,   9.00d0,   9.00d0,   9.00d0,   9.00d0, &
            9.00d0,   9.00d0,   9.00d0,   9.00d0,   9.00d0,   9.00d0, &
            9.00d0,   9.00d0,   9.00d0,   9.00d0,  10.00d0,  10.00d0/


!..years 1711 to 1799
      data  (dt(i), i=91, 180) / &
           10.00d0,  10.00d0,  10.00d0,  10.00d0,  10.00d0,  10.00d0,  &
           10.00d0,  11.00d0,  11.00d0,  11.00d0,  11.00d0,  11.00d0,  &
           11.00d0,  11.00d0,  11.00d0,  11.00d0,  11.00d0,  11.00d0,  &
           11.00d0,  11.00d0,  11.00d0,  11.00d0,  11.00d0,  11.00d0,  &
           12.00d0,  12.00d0,  12.00d0,  12.00d0,  12.00d0,  12.00d0,  &
           12.00d0,  12.00d0,  12.00d0,  12.00d0,  13.00d0,  13.00d0,  &
           13.00d0,  13.00d0,  13.00d0,  13.00d0,  13.00d0,  14.00d0,  &
           14.00d0,  14.00d0,  14.00d0,  14.00d0,  14.00d0,  14.00d0,  &
           15.00d0,  15.00d0,  15.00d0,  15.00d0,  15.00d0,  15.00d0,  &
           15.00d0,  16.00d0,  16.00d0,  16.00d0,  16.00d0,  16.00d0,  &
           16.00d0,  16.00d0,  16.00d0,  16.00d0,  16.00d0,  17.00d0,  & 
           17.00d0,  17.00d0,  17.00d0,  17.00d0,  17.00d0,  17.00d0,  &
           17.00d0,  17.00d0,  17.00d0,  17.00d0,  17.00d0,  17.00d0,  &
           17.00d0,  17.00d0,  17.00d0,  17.00d0,  16.00d0,  16.00d0,  &
           16.00d0,  16.00d0,  15.00d0,  15.00d0,  14.00d0,  14.00d0/


!..years 1800 to 1890
      data  (dt(i), i=181, 270) / &
           13.70d0,  13.40d0,  13.10d0,  12.90d0,  12.70d0,  12.60d0, &
           12.50d0,  12.50d0,  12.50d0,  12.50d0,  12.50d0,  12.50d0, &
           12.50d0,  12.50d0,  12.50d0,  12.50d0,  12.50d0,  12.40d0, & 
           12.30d0,  12.20d0,  12.00d0,  11.70d0,  11.40d0,  11.10d0, & 
           10.60d0,  10.20d0,   9.60d0,   9.10d0,   8.60d0,   8.00d0, & 
            7.50d0,   7.00d0,   6.60d0,   6.30d0,   6.00d0,   5.80d0, &
            5.70d0,   5.60d0,   5.60d0,   5.60d0,   5.70d0,   5.80d0, &
            5.90d0,   6.10d0,   6.20d0,   6.30d0,   6.50d0,   6.60d0, & 
            6.80d0,   6.90d0,   7.10d0,   7.20d0,   7.30d0,   7.40d0, &
            7.50d0,   7.60d0,   7.70d0,   7.70d0,   7.80d0,   7.80d0, &
            7.88d0,   7.82d0,   7.54d0,   6.97d0,   6.40d0,   6.02d0, &
            5.41d0,   4.10d0,   2.92d0,   1.82d0,   1.61d0,   0.10d0, &
           -1.02d0,  -1.28d0,  -2.69d0,  -3.24d0,  -3.64d0,  -4.54d0, &
           -4.71d0,  -5.11d0,  -5.40d0,  -5.42d0,  -5.20d0,  -5.46d0, &
           -5.46d0,  -5.79d0,  -5.63d0,  -5.64d0,  -5.80d0,  -5.66d0/


!..years 1890 to 1979
      data  (dt(i), i=271, 360) / &
           -5.87d0,  -6.01d0,  -6.19d0,  -6.64d0,  -6.44d0,  -6.47d0, & 
           -6.09d0,  -5.76d0,  -4.66d0,  -3.74d0,  -2.72d0,  -1.54d0, & 
           -0.2d0,   1.24d0,   2.64d0,   3.86d0,   5.37d0,   6.14d0, & 
            7.75d0,   9.13d0,  10.46d0,  11.53d0,  13.36d0,  14.65d0, & 
           16.01d0,  17.20d0,  18.24d0,  19.06d0,  20.25d0,  20.95d0, & 
           21.16d0,  22.25d0,  22.41d0,  23.03d0,  23.49d0,  23.62d0, &  
           23.86d0,  24.49d0,  24.34d0,  24.08d0,  24.02d0,  24.00d0, &  
           23.87d0,  23.95d0,  23.86d0,  23.93d0,  23.73d0,  23.92d0, &  
           23.96d0,  24.02d0,  24.33d0,  24.83d0,  25.30d0,  25.70d0, &  
           26.24d0,  26.77d0,  27.28d0,  27.78d0,  28.25d0,  28.71d0, & 
           29.15d0,  29.57d0,  29.97d0,  30.36d0,  30.72d0,  31.07d0, &  
           31.35d0,  31.68d0,  32.18d0,  32.68d0,  33.15d0,  33.59d0, &  
           34.00d0,  34.47d0,  35.03d0,  35.73d0,  36.54d0,  37.43d0, &  
           38.29d0,  39.20d0,  40.18d0,  41.17d0,  42.23d0,  43.37d0, &  
           44.49d0,  45.48d0,  46.46d0,  47.52d0,  48.53d0,  49.59d0/

!..years 1980 to 2013 
      data  (dt(i), i=361, 396) / &
           50.54d0,  51.38d0,  52.17d0,  52.96d0,  53.79d0,  54.34d0,  &
           54.87d0,  55.32d0,  55.82d0,  56.30d0,  56.86d0,  57.57d0,  &
           58.31d0,  59.12d0,  59.98d0,  60.78d0,  61.63d0,  62.29d0,  &
           62.97d0,  63.47d0,  63.83d0,  64.09d0,  64.30d0,  64.47d0,  &
           64.57d0,  64.69d0,  64.85d0,  65.15d0,  65.46d0,  65.78d0,  & !  last entry here is 2009
           66.07d0,  66.32d0, &      ! last entry here is 2011
           66.6d0,   66.9d0,   67.2d0,   67.5d0/
!  NOTE:  Values thru 2011 ok, 2012 & above are extrapolated.  -RDR     

!..first time flag
      data iflag/0/



!..if this is the first time, initialize 
!..put the julain date at the first of the year in the array
!..years so that we can interpolate/extrapolate directly on
!..the given ut julian date

      if (iflag .eq. 0) then
       iflag = 1
  
       do i=1,tabsize
        j = tabstart + (i-1)
        call juldat2(j,1,1,0.0d0,xx)
        years(i) = xx
       end do

      endif


!..convert the given ut julian date to a ut calander date
      call caldat2(tjd,iy,im,id,hr)


!..if we are outside the table on the low end
!..use the stephenson and morrison expression 948 to 1600,
!..and the borkowski formula for earlier years

      if (iy .lt. tabstart) then
       if (iy .ge. 948) then
        b      = 0.01d0 * float(iy - 2000)
        secdif = b * (b * 23.58d0 + 100.3d0) + 101.6d0
       else
        b      = 0.01d0 * float(iy -2000) + 3.75d0
        secdif = 35.0d0 * b * b + 40.0d0
       end if


!..if we are outside the table on the high end
!..use a linear extrapolation into the future

      else if (iy .gt. tabend) then
       b      = float(iy - tabend)
       secdif = dt(tabsize) + b*(dt(tabsize) - dt(tabsize-1))



!..otherwise we are in the table
!..get the table location and interpolate

      else
       iat = iy - tabstart + 1
       iat = max(1, min(iat - mp/2 + 1, tabsize - mp + 1))
       call polint(years(iat),dt(iat),mp,tjd,secdif,dy)


!..the astronomical almanac table is corrected by adding the expression
!..      -0.000091 (ndot + 26)(year-1955)^2  seconds
!..to entries prior to 1955 (page K8), where ndot is the secular tidal 
!..term in the mean motion of the moon. entries after 1955 are referred 
!..to atomic time standards and are not affected by errors in lunar 
!..or planetary theory.  a value of ndot = -25.8 arcsec per century squared 
!..is the value used in jpl's de403 ephemeris, the earlier de200 ephemeris
!..used the value -23.8946. note for years below the table (less than 1620)
!..the time difference is not adjusted for small improvements in the 
!..current estimate of ndot because the formulas were derived from 
!..studies of ancient eclipses and other historical information, whose 
!..interpretation depends only partly on ndot.
!..here we make the ndot correction.

       if (iy .lt. 1955) then
        b = float(iy - 1955)
        secdif = secdif - 0.000091d0 * (-25.8d0 + 26.0d0)*b*b
       end if
      end if



!..add the difference to the ut julian date to get the dynamical julian date
!     tjde = tjd + secdif/86400.0d0

      deltat = secdif

      return
      end




      subroutine polint(xa,ya,n,x,y,dy)
      implicit none
      save

!..given arrays xa and ya of length n and a value x, this routine returns a 
!..value y and an error estimate dy. if p(x) is the polynomial of degree n-1
!..such that ya = p(xa) ya then the returned value is y = p(x) 

!..input:
!..xa(1:n) = array of x values
!..ya(1:n) = array of y values
!..n       = order of interpolant, 2=linear, 3=quadratic ...
!..x       = x value where interpolation is desired

!..output:
!..y       = interpolated y value
!..dy      = error esimate


!..declare
      integer          n,nmax,ns,i,m
      parameter        (nmax=10)
      double precision xa(n),ya(n),x,y,dy,c(nmax),d(nmax),dif,dift, &
                       ho,hp,w,den

!..find the index ns of the closest table entry; initialize the c and d tables
      ns  = 1
      dif = abs(x - xa(1))
      do i=1,n
       dift = abs(x - xa(i))
       if (dift .lt. dif) then
        ns  = i
        dif = dift
       end if
       c(i)  = ya(i)
       d(i)  = ya(i)
      enddo

!..first guess for y
      y = ya(ns)

!..for each column of the table, loop over the c's and d's and update them
      ns = ns - 1
      do m=1,n-1
       do i=1,n-m
        ho   = xa(i) - x
        hp   = xa(i+m) - x
        w    = c(i+1) - d(i)
        den  = ho - hp
        if (den .eq. 0.0) stop ' 2 xa entries are the same in polint'
        den  = w/den
        d(i) = hp * den
        c(i) = ho * den
       enddo

!..after each column is completed, decide which correction c or d, to add
!..to the accumulating value of y, that is, which path to take in the table
!..by forking up or down. ns is updated as we go to keep track of where we
!..are. the last dy added is the error indicator.
       if (2*ns .lt. n-m) then
        dy = c(ns+1)
       else
        dy = d(ns)
        ns = ns - 1
       end if
       y = y + dy
      enddo
      return
      end


!..The next two routines are alternatives to the juldat and caldat routines
!..in the novas package, since the novas routines do not take care of the
!..switch from julian to gregorian calanders on 15oct1582. Hence,
!..the novas routines are of no use with the long jpl ephemris de406.

!..Input calendar dates 1582-oct-15 and after are taken to be expressed 
!..in the gregorian calendar system. Prior dates are assumed to be in the 
!..julian calendar system.

!..Historically, not all regions switched calendars at the same time 
!..(or even in the same century). Thus, the user must be aware of
!..which calendar was in effect for a particular historical record. 
!..It should not be assumed this system's calendar automatically
!..correlates with a date from an arbitrary historical document. 

!..Here is the progression near the calendar switch point: 
!
!       calendar type    calendar date   julian day number
!       -------------    -------------   -----------------
!        julian           1582-oct-03        2299158.5
!        julian           1582-oct-04        2299159.5 --->
!         (skipped)      "1582-oct-05"                    |
!         (skipped)      "1582-oct-06"                    |
!         (skipped)      "1582-oct-07"                    |
!         (skipped)      "1582-oct-08"                    |
!         (skipped)      "1582-oct-09"                    |
!         (skipped)      "1582-oct-10"                    |
!         (skipped)      "1582-oct-11"                    |
!         (skipped)      "1582-oct-12"                    |
!         (skipped)      "1582-oct-13"                    |
!         (skipped)      "1582-oct-14"                    |
!        gregorian        1582-oct-15        2299160.5 <---
!        gregorian        1582-oct-16        2299161.5
!        gregorian        1582-oct-17        2299162.5


!..in this system there are zero and negative years. 
!..the progression is as follows: 
!
!    julian day number       labeling-convention
!      (jan 1 00:00)       BC/AD      arithmetical 
!    -----------------     -----      ------------
!    1720327.5              3BC           -2
!    1720692.5              2BC           -1
!    1721057.5              1BC            0
!    1721423.5              1AD            1
!    1721788.5              2AD            2


!  Author: F. X. Timmes, U. Chicago



      subroutine juldat2(iy,im,id,rh,tjd)
      implicit none
      save

!..converts a calander date to a julian date.
!..input time value can be in any ut-like time scale (utc, ut1, tt, etc.) 
!..and the output julian date will have same basis.

!..the astronomical calendar is used. thus the year before 1 AD is 0, 
!..the year before that is 1 BC. the change to the gregorian calander 
!..on oct 15, 1582 is accounted for.

!..input:
!..iy = integer year
!..im = integer month
!..id = integer day
!..rh = number of hours past midnight

!..output:
!..tjd = julian date


!..declare
      integer          im,iy,id,igreg,i,jd
      parameter        (igreg = 15+31*(10+12*1582))
      double precision tjd,xy,xm,xa,xb,rh

      if (im .gt. 2) then
       xy = iy
       xm = im + 1
      else
       xy = iy - 1
       xm = im + 13
      end if
      i = id + 31 * (im + 12 * iy)
      if (i .ge. igreg) then
       xa = int(0.01d0 * xy)
       xb = 2.0d0 - xa + int(0.25d0 * xa)
      else
       xb = 0.0d0
      end if
      jd  = int(365.25d0*(xy + 4716.0d0)) + int(30.6001d0*xm) + id
      tjd = dble(jd) + xb - 1524.5d0 + rh/24.d0
      return
      end






      subroutine caldat2(tjd,iy,im,id,rh)
      implicit none
      save

!..converts a julian date to a calander date.
!..input time value can be in any ut-like time scale (utc, ut1, tt, etc.) 
!..and the output calander date will have same basis.

!..the astronomical calendar is used. thus the year before 1 ad is 0, 
!..the year before that is 1 bc. the change to the gregorian calander 
!..on oct 15, 1582 is accounted for.

!..input:
!..tjd = julian date
!..iy = integer year
!..im = integer month
!..id = integer day
!..rh = number of hours past midnight

!..output:
!..iy = integer year
!..im = integer month
!..id = integer day
!..rh = number of hours past midnight


!..declare
      integer          id,im,iy,igreg
      parameter        (igreg = 2299161)
      double precision tjd,rh,x1,z,f,x2,xa,xb,xc,xd,xe,rd,c1,c2,c3
      parameter        (c1 = 1.0d0/36524.25d0, &
                        c2 = 1.0d0/365.25d0, &
                        c3 = 1.0d0/30.6001d0) 

      x1 = tjd + 0.5d0
      z  = int(x1)
      f  = x1 - z

      if (x1 .ge. igreg) then
       x2 = int((x1-1867216.25d0)*c1)
       xa = z + 1.0d0 +x2 - int(0.25d0 * x2)
      else
       xa = z
      end if
      xb = xa + 1524.0d0
      xc = int((xb - 122.1d0)*c2)
      xd = int(365.25d0*xc)
      xe = int((xb-xd)*c3)

      rd = xb - xd - int(30.6001d0*xe) + f
      id = rd
      rh = 24.0d0*(rd - dble(id)) 
      im = xe - 1
      if (im .gt. 12) im = im - 12
      iy = xc - 4715
      if (im .gt. 2) iy = iy - 1
      return
      end
!------------------------------------------------------------------------------
         subroutine def_cid(nc0,cid)
         implicit none
         include 'derived_types.h'
         include 'constit.h'
         integer(kind=4) nc0,ic,jc,k
         type (constituents) cid
!
         k=1
         do ic = 1,nc0
          do jc = 1,ncmx             
            if(cid%cons(ic).eq.constid(jc)) go to 4
          enddo
          write(*,*) 'WARNING: Constituent ID ',cid%cons(ic), &
                     ' is not allowed'
          write(*,*) '         in constit_f90.h'
          cid%idc(k)=0
          go to 5      
4         continue
          cid%idc(k)=jc
5         k=k+1
         enddo
         return
         end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      character*10 function deblank(str)
      character*10 str
      integer k
1     k=index(str,' ')
      if(k.gt.0)then
       str(k:k)='0'
       go to 1
      endif
      deblank=str
      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc Lana Erofeeva, Feb 2002 cccccccccccccccccccccccccccccccccccccc
!cc Debugged for consistency with CALDAT for any dates mjd>=0
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
      subroutine date_mjd(mm,id,iyyy,mjd)
      implicit none
! converts date to mjd
! INPUT:  id - day, mm - month, iyyy - year
! OUTPUT: mjd>0 - modified julian days
! date>=11.17.1858 corresponds to mjd=0
      integer dpm(12),days,i,nleap,k
      integer mm,id,iyyy,mjd
      data dpm/31,28,31,30,31,30,31,31,30,31,30,31/
      mjd=0
! NO earlier dates
      if(iyyy.lt.1858)iyyy=1858
      if(iyyy.eq.1858.and.mm.gt.11)mm=11
      if(iyyy.eq.1858.and.mm.eq.11.and.id.gt.17)id=17
!
      days=0 
      do i=1,mm-1
       days=days+dpm(i)
       if(i.eq.2.and.int(iyyy/4)*4.eq.iyyy)days=days+1         
      enddo
      days=days+id-321
! leap day correction
      do k=1900,iyyy,100
       if(iyyy.eq.k.and.mm.gt.2)days=days-1
      enddo
      do k=2000,iyyy,400
       if(iyyy.eq.k.and.mm.gt.2)days=days+1
      enddo
! EACH 4th year is leap year
      nleap=int((iyyy-1-1860)*0.25)            
      if(iyyy.gt.1860)nleap=nleap+1
! EXCEPT
      do k=1900,iyyy-1,100
       if(iyyy.gt.k)nleap=nleap-1
       if(iyyy.eq.k.and.mm.gt.2)days=days-1
      enddo
! BUT each in the row 2000:400:... IS LEAP year again
      do k=2000,iyyy-1,400
       if(iyyy.gt.k)nleap=nleap+1
       if(iyyy.eq.k.and.mm.gt.2)days=days+1
      enddo
      mjd=365*(iyyy-1858)+nleap &
          +days   
      return
      end
!-------------------------------------
      SUBROUTINE CALDAT (JULIAN,MM,ID,IYYY)                             
!                                                                       
!   ROUTINE CONVERTS JULIAN DAY TO MONTH, DAY, & YEAR.                      
!   THIS CODE IS LIFTED FROM THE BOOK:                                      
!   W.H. PRESS ET AL., NUMERICAL RECIPES, CAMBRIDGE UNIV. PRESS, 1986.      
!   THE ONLY MODIFICATION IS THAT REAL ARITHMETIC IS DONE IN R*8.           
!                                                                           
!   TO CONVERT MODIFIED JULIAN DAY, CALL THIS ROUTINE WITH                 
!     JULIAN = MJD + 2400001                                               
!                                                                          
      PARAMETER (IGREG=2299161)                                           
      IF (JULIAN.GE.IGREG) THEN                                            
         JALPHA=INT((DBLE(JULIAN-1867216)-0.25D0)/36524.25D0)              
         JA=JULIAN+1+JALPHA-INT(0.25D0*DBLE(JALPHA))                       
      ELSE                                                                 
         JA=JULIAN                                                         
      ENDIF                                                                
      JB=JA+1524                                                           
      JC=INT(6680.D0+(DBLE(JB-2439870)-122.1D0)/365.25D0)                  
      JD=365*JC+INT(0.25D0*DBLE(JC))                                       
      JE=INT(DBLE(JB-JD)/30.6001D0)                                        
      ID=JB-JD-INT(30.6001D0*DBLE(JE))                                     
      MM=JE-1                                                              
      IF (MM.GT.12) MM=MM-12                                               
      IYYY=JC-4715                                                         
      IF (MM.GT.2) IYYY=IYYY-1                                             
      IF (IYYY.LE.0) IYYY=IYYY-1                                           
      RETURN                                                               
      END           
!--------------------------------------------------------------------------------------
      subroutine def_model(lat,lon,theta_lim,phi_lim,&
                           n,m,pmask,mid)
      real*4 lat,lon,lon1,theta_lim(2),phi_lim(2),dx,dy
      integer*4 i,j,n,m,k,l,kt,mid,pmask(n,m),tmp(3,3),median
!
      mid=0
      dx=(phi_lim(2)-phi_lim(1))/n
      dy=(theta_lim(2)-theta_lim(1))/m
      lon1=lon
      if(lon1.lt.phi_lim(1))lon1=lon1+360
      if(lon1.gt.phi_lim(2))lon1=lon1-360
      if(lon1.lt.phi_lim(1).or.lon1.gt.phi_lim(2))return
      j=int((lat-theta_lim(1))/dx+0.5)
      i=int((lon1-phi_lim(1))/dy+0.5)
      j=max(j,1);j=min(j,m)
      i=max(i,1);i=min(i,n)
      !write(*,*)lat,phi_lim,lon,lon1,i,j
      tmp=0
      do k=i-1,i+1
       kt=k
       if(kt<0)kt=kt+n
       if(kt>n)kt=kt-n
       do l=max(j-1,1),min(j+1,m)
        tmp(k-i+2,l-j+2)=pmask(kt,l)
       enddo
      enddo
      mid=maxval(tmp)
      return
      end
!--------------------------------------------------------------------------------------
      subroutine def_modelC(lat,lon,theta_lim,phi_lim,&
                           n,m,pmask,grid,mid,nodeID)
      implicit none
      include 'derived_types.h'
      type (c_grid) grid
      real*4 lat,lon,lon1,theta_lim(2),phi_lim(2),dx,dy,hmax
      integer*4 i,j,n,m,k,l,kt,mid,pmask(n,m),tmp(3,3),iz,it
      integer*4 nodeID
!
      mid=0
      dx=(phi_lim(2)-phi_lim(1))/n
      dy=(theta_lim(2)-theta_lim(1))/m
      lon1=lon
      if(lon1.lt.phi_lim(1))lon1=lon1+360
      if(lon1.gt.phi_lim(2))lon1=lon1-360
      if(lon1.lt.phi_lim(1).or.lon1.gt.phi_lim(2))return
      j=int((lat-theta_lim(1))/dx+0.5)
      i=int((lon1-phi_lim(1))/dy+0.5)
      j=max(j,1);j=min(j,m)
      i=max(i,1);i=min(i,n)
      !write(*,*)lat,phi_lim,lon,lon1,i,j
      tmp=0
      nodeID=1 ! sea
      iz=0
      hmax=0
      do k=i-1,i+1
       kt=k
       if(kt<0)kt=kt+n
       if(kt>n)kt=kt-n
       do l=max(j-1,1),min(j+1,m)
        tmp(k-i+2,l-j+2)=pmask(kt,l)
        if(grid%mz(kt,l).eq.0)nodeID=-1 ! coast
        if(grid%hz(kt,l).gt.hmax)hmax=grid%hz(kt,l)     
        iz=iz+grid%mz(kt,l)
       enddo
      enddo
      !if(hmax.lt.100)nodeID=-1
      mid=maxval(tmp)
      it=0
      do k=1,3
       do l=1,3
        if(tmp(k,l).eq.mid)it=it+1
       enddo
      enddo
!
      if(mid.ne.0.and.minval(tmp).eq.0.and. &
         iz.eq.9.and.nodeID.gt.0)nodeID=it+1 ! boundary btw local & global sol
                                             ! in deep ocean
      return
      end
!-----------------------------------------------------------------------------------------
      integer*4 function median(a,n,m)
      integer*4 a(n,m),n,m,i,j,c
      integer*4, allocatable:: b(:)
      allocate(b(n*m))
      k=1
      do i=1,n
       do j=1,m
         b(k)=a(i,j)
         k=k+1
       enddo
      enddo
! sort b
1     continue
      do i=2,n*m
       if(b(i).lt.b(i-1))then
        c=b(i);b(i)=b(i-1);b(i-1)=c
        go to 1
       endif
      enddo
      median=b(n*m/2)
      deallocate(b)  
      return
      end
!----------------------------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
      subroutine getDims(cgrid,grid)
      implicit none
      include 'derived_types.h'
      integer(kind=4) j
      real(kind=4)eps
      character cgrid*(*)
      type (c_grid) grid
      character*2 c2
!
!cc   READ IN GRID FILE HEADER
!
       grid%NPG=.false.
       j=index(cgrid,'/',back=.true.)
       !write(*,*)'Grid file:',trim(cgrid)
       if(cgrid(j+1:j+8).eq.'grid_NPG')then
         grid%NPG=.true.
         write(*,*)'ROTATED coordinate system (NPG)'
       endif
!
      open(unit=1,file=trim(cgrid),form='unformatted',status='old',&
            err=1)
      read(1,err=1) grid%n,grid%m,grid%theta_lim,grid%phi_lim,grid%dt,grid%nob
      close(1)
!cc Check for global
      eps=0.25*(grid%phi_lim(2)-grid%phi_lim(1))/grid%n
      grid%ll_km =.false.
      j=index(cgrid,'km')
      if(j.gt.0)grid%dt=-abs(grid%dt) ! this is for VERY obsolete grids
      if(grid%dt.lt.0)then ! grid is stereographic
       grid%dt=abs(grid%dt)
       grid%km_area_id=int(abs(grid%dt)+0.5) ! this will be in future, not for old grids
       grid%ll_km=.true.
       grid%x_lim=grid%phi_lim
       grid%y_lim=grid%theta_lim
       call what_km_area(grid%n,grid%m,grid%theta_lim,grid%phi_lim,grid%km_area_id)
       call def_ll_lims(grid) ! define grid limits in lats, lons
      endif
      if(abs(grid%phi_lim(2)-grid%phi_lim(1)-360.).lt.eps.and..not.grid%ll_km)then
        grid%lglobe=.true.
        !write(*,*)'Grid is global'
      else
        grid%lglobe=.false.
        !write(*,*)'Grid is NOT global'
      endif
      return
1     cgrid=trim(cgrid)
      write(*,*)'File NOT found or wrong format:'
      write(*,*)'*',trim(cgrid),'*'
      grid%n=0
      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! reads in all grid arrays ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid_in(cgrid,grid,rc_msk)
      implicit none
      include 'derived_types.h'
      integer(kind=4) j,k,n_x2,n,m,nob
      real(kind=4)dt,theta_lim(2),phi_lim(2)
      logical, optional:: rc_msk
      character cgrid*(*)
      type (c_grid),intent(inout):: grid
!
      if(.not.present(rc_msk))then
        rc_msk=.false.
      endif      
      if(.not.grid%allocated)stop 911
      open(unit=1,file=trim(cgrid),form='unformatted',status='old')
! pass header: header has been read before grid was allocated
      read(1)n,m,theta_lim,phi_lim,dt,nob
!     read open boundary conditions
      grid%iob(:,:)=0
      if(nob.gt.0)then
        read(1) ((grid%iob(k,j),k=1,2),j=1,nob)
	do k=1,nob
          if (grid%iob(1,k).lt.1.or.grid%iob(1,k).gt.n) go to 1
	  if (grid%iob(2,k).lt.1.or.grid%iob(2,k).gt.m) go to 1
	enddo
      else
        read(1)
      endif
!     read bathymetry for z nodes
      read(1) grid%hz
!     read in z node mask
      read(1) grid%mz
      close(1)
!
      if(.not.rc_msk)then      ! for global atlas zero mask is everywhere
                               ! where global is patched with local models 
       grid%hz=grid%hz*grid%mz
      else ! recover mask for atlas
       grid%mz=1
       where(grid%hz.eq.0)grid%mz=0
      endif
!
! this check used to be in repx,fwd_fac,rlc: does not belong there!
      if(grid%lglobe)then
        n_x2=0
        do j=1,grid%m
         if(grid%mz(grid%n,j).gt.0)n_x2=n_x2+1
        enddo
        if(n_x2.eq.0)then
         write(*,*)'Right grid column is masked off or has OB nodes'
         write(*,*)'while grid seems to be global.'
         write(*,*)'Check grid and OB!'
         stop
        endif
      endif
      return
1     write(*,*)'Error in OB definition'
      write(*,*)'OB index is out of grid size'
      stop
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine getDims_model(fname,model)
      implicit none
      include 'derived_types.h'
      integer(kind=4) j,n,m
      real(kind=4)eps,eps1
      character*80 fname
      type (z_model) model
!
! read in model header
!
      model%NPG=.false.
      j=index(fname,'NPG')
      if(j.gt.0)then
         model%NPG=.true.
         write(*,*)'Model ', trim(fname),' is in ROTATED coordinate system (NPG)'
      endif
!
      open(unit=1,file=trim(fname),form='unformatted',status='old',&
            err=1)
      read(1,err=1) model%n,model%m,model%nc,model%theta_lim,model%phi_lim
      close(1)
!
! Check for km case. This is NOT perfect recognition.
! Should put km_area_id somewhere in the model files.
! May be add last short record? Shall decide later.
!
      n=model%n;m=model%m
      model%ll_km =.false.
      j=index(fname,'km')
      eps=model%theta_lim(2)-model%theta_lim(1)
      eps1=model%phi_lim(2)-model%phi_lim(1)
      model%km_area_id=0
      if(j.gt.0.or.eps.gt.180+2*eps/m.or.eps1.gt.360+2*eps1/n.or. &
         model%theta_lim(1).lt.-90-2*eps/m .or. model%theta_lim(2).gt.90+2*eps/m .or. &
         model%phi_lim(1)  .lt.-180-2*eps1/n .or. model%phi_lim(2).gt.360+2*eps1/n)then
        model%ll_km=.true.
        write(*,*)'Model is on stereographic grid'
        call what_km_area(model%n,model%m,model%theta_lim,model%phi_lim,model%km_area_id)
        model%x_lim=model%phi_lim
        model%y_lim=model%theta_lim
        call def_ll_lims0(model%x_lim,model%y_lim,model%km_area_id, &
                          model%phi_lim,model%theta_lim)
        write(*,*)'Approximate lat limits:',model%theta_lim
        write(*,*)'Approximate lon limits:',model%phi_lim
      endif
!cc Check for global
      eps=0.25*(model%phi_lim(2)-model%phi_lim(1))/model%n
      if(abs(model%phi_lim(2)-model%phi_lim(1)-360.).lt.eps)then
        model%lglobe=.true.
        !write(*,*)'Model is global'
      else
        model%lglobe=.false.
        !write(*,*)'Model is NOT global'
      endif
      return
1     model%n=0
      write(*,*)'File ',trim(fname),' NOT found or wrong format'
      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine what_km_area(ng,mg,th_grid,ph_grid,km_area_id)
       implicit none
       integer(kind=4) ng,mg,n0,m0,km_area_id
       real(kind=4) dx,dy,th_grid(2),ph_grid(2),x1,y1
       real(kind=8) dlat,dlon,tx,ty,SLAT,SLON
       character*80 ll_km_sub
!
       open(unit=26,file='subname',form='formatted', &
            status='old')
       read(26,*)ll_km_sub
       close(26)
       km_area_id=0
       if(trim(ll_km_sub).eq.'xy_ll')km_area_id=1          ! Arctic
       if(trim(ll_km_sub).eq.'xy_ll_S')km_area_id=2        ! Antarctic
       if(trim(ll_km_sub).eq.'mapxy,SLON=+70')km_area_id=3 ! Amery Ice Shelf
       if(trim(ll_km_sub).eq.'mapxy,SLON=-70')km_area_id=4 ! Weddel/Ross Sea or CATs
       if(trim(ll_km_sub).eq.'mapxy,SLON=150')km_area_id=5 ! Mertz area of Anatractic
       if(km_area_id.eq.0)then
         write(*,*)'Can NOT identify converting function!'
         write(*,*)'Make sure enter function names EXACTLY AS WRITTEN BELOW'
         write(*,*)'in 4th line of Model_* file:'
         write(*,*)'KNOWN FUNCTIONS:'
         write(*,*)'1. Arctic                     xy_ll'
         write(*,*)'2. Antarctic                  xy_ll_S'
         write(*,*)'3. Amery Ice Shelf            mapxy,SLON=+70'
         write(*,*)'4. Weddel/Ross Sea or CATs    mapxy,SLON=-70'
         write(*,*)'5. Mertz area of Antarctic    mapxy,SLON=150'
         stop
       endif
111    write(*,*)'AREA IDENTIFIED:'
       if(km_area_id.eq.1)write(*,*)'1. Arctic (xy_ll,ll_xy)'
       if(km_area_id.eq.2)write(*,*)'2. Antarctic (xy_ll_S,ll_xy_S)'
       if(km_area_id.eq.3)write(*,*)'3. Amery Ice Shelf (mapxy, mapll, SLON=+70)'
       if(km_area_id.eq.4)write(*,*)'4. Weddel/Ross Sea or CATs (mapxy, mapll, SLON=-70)'
       if(km_area_id.eq.5)write(*,*)'5. Mertz area of Antarctic (mapxy, mapll, SLON=150)'
       return
       end
!--------------------------------------------------------------------------------
	subroutine def_ll_lims(grid)
        implicit none
        include 'derived_types.h'
        integer(kind=4)k
        real(kind=8)x(4),y(4),lon(4),lat(4)
        type (c_grid)grid
!
! 4 corner points of stereographic grid
        x(1)=grid%x_lim(1);x(2)=grid%x_lim(2);x(3)=grid%x_lim(2);x(4)=grid%x_lim(1)
        y(1)=grid%y_lim(1);y(2)=grid%y_lim(1);y(3)=grid%y_lim(2);y(4)=grid%y_lim(2)
!
        do k=1,4
          call convert_xy_ll(grid%km_area_id,x(k),y(k),lon(k),lat(k))
        enddo
        if(grid%x_lim(1)*grid%x_lim(2).lt.0.and. &
           grid%y_lim(1)*grid%y_lim(2).lt.0)then ! grid is over a pole
          grid%phi_lim(1)=0
          grid%phi_lim(2)=360
          if(lat(1).lt.0)then
           grid%theta_lim(2)=maxval(lat)
           grid%theta_lim(1)=-90.
          else
           grid%theta_lim(1)=minval(lat)
           grid%theta_lim(2)=90.
          endif
        else
         grid%phi_lim=lon(1);grid%theta_lim=lat(1)
         do k=1,4
          if(lon(k).lt.grid%phi_lim(1))grid%phi_lim(1)=lon(k)
          if(lon(k).gt.grid%phi_lim(2))grid%phi_lim(2)=lon(k)
          if(lat(k).lt.grid%theta_lim(1))grid%theta_lim(1)=lat(k)
          if(lat(k).gt.grid%theta_lim(2))grid%theta_lim(2)=lat(k)
         enddo
         if(grid%phi_lim(1).lt.-180.)grid%phi_lim=grid%phi_lim+360
        endif
        return
        end
! ---------------------------------------------------------------------------------
        logical function not_in_lims(lat,lon,grid)
        implicit none
        include 'derived_types.h'
        type (c_grid), intent(in):: grid
	real(kind=4)lat,lon,theta_lim(2),phi_lim(2)
        real(kind=8) xlon,xlat,x,y
        logical not_in_lims_km
        theta_lim=grid%theta_lim
        phi_lim=grid%phi_lim
        not_in_lims=.false.
        if(lat.gt.theta_lim(2).or.lat.lt.theta_lim(1))not_in_lims=.true.
        if(lon.lt.phi_lim(1).and.lon+360.gt.phi_lim(2))not_in_lims=.true.
        if(lon.gt.phi_lim(2).and.lon-360.lt.phi_lim(1))not_in_lims=.true.
        if(lon.lt.phi_lim(1).and..not.not_in_lims)lon=lon+360.
        if(lon.gt.phi_lim(2).and..not.not_in_lims)lon=lon-360.
        if(not_in_lims)return
        if(grid%ll_km)then
         xlon=lon;xlat=lat
         call convert_ll_xy(grid%km_area_id,xlon,xlat,x,y)
         not_in_lims=not_in_lims_km(x,y,grid%x_lim,grid%y_lim)
        endif
        return
        end
! ---------------------------------------------------------------------------------
        logical function not_in_lims_km(x,y,x_lim,y_lim)
        implicit none
	real(kind=4)x_lim(2),y_lim(2)
        real(kind=8)x,y
        not_in_lims_km=.false.
        if(x.gt.x_lim(2).or.x.lt.x_lim(1))not_in_lims_km=.true.
        if(y.gt.y_lim(2).or.y.lt.y_lim(1))not_in_lims_km=.true.
        return
        end
! ----------------------------------------------------------------------!
	subroutine def_ll_lims0(x_lim,y_lim,km_area_id, &
                                phi_lim,theta_lim)
        implicit none
        real(kind=4)phi_lim(2),theta_lim(2),x_lim(2),y_lim(2)
        integer(kind=4)km_area_id,k
        real(kind=8)x(4),y(4),lon(4),lat(4)
!
        x(1)=x_lim(1);x(2)=x_lim(2);x(3)=x_lim(2);x(4)=x_lim(1)
        y(1)=y_lim(1);y(2)=y_lim(1);y(3)=y_lim(2);y(4)=y_lim(2)
!
        do k=1,4
          call convert_xy_ll(km_area_id,x(k),y(k),lon(k),lat(k))
        enddo
        if(x_lim(1)*x_lim(2).lt.0.and.y_lim(1)*y_lim(2).lt.0)then ! grid is over a pole
          phi_lim(1)=0
          phi_lim(2)=360
          if(lat(1).lt.0)then
           theta_lim(2)=maxval(lat)
           theta_lim(1)=-90.
          else
           theta_lim(1)=minval(lat)
           theta_lim(2)=90.
          endif
        else
         phi_lim=lon(1);theta_lim=lat(1)
         do k=1,4
          if(lon(k).lt.phi_lim(1))phi_lim(1)=lon(k)
          if(lon(k).gt.phi_lim(2))phi_lim(2)=lon(k)
          if(lat(k).lt.theta_lim(1))theta_lim(1)=lat(k)
          if(lat(k).gt.theta_lim(2))theta_lim(2)=lat(k)
         enddo
         if(phi_lim(1).lt.-180.)phi_lim=phi_lim+360
        endif
        return
        end
! ----------------------------------------------------------------------------
        subroutine convert_xy_ll(km_area_id,x,y,lon,lat)
        integer(kind=4) km_area_id
        real(kind=8)x,y,lon,lat,SLAT,SLON
        character*1 HEMI
!
        HEMI='S'
              selectcase(km_area_id)
              case(1)
                call xy_ll(x,y,lon,lat)
              case(2)
                call xy_ll_S(x,y,lon,lat)
              case(3)
                SLAT=71.;SLON=70.
                call mapxy(1,1,x,y,lon,lat,SLAT,SLON,HEMI)
              case(4)
                SLAT=71.;SLON=-70.
                call mapxy(1,1,x,y,lon,lat,SLAT,SLON,HEMI)
              case(5)
                SLAT=71.;SLON=150.
                call mapxy(1,1,x,y,lon,lat,SLAT,SLON,HEMI)
              endselect
!
        return
        end
! ----------------------------------------------------------------------------
        subroutine convert_ll_xy(km_area_id,lon,lat,x,y)
        integer(kind=4) km_area_id
        real(kind=8)x,y,lon,lat,SLAT,SLON
        character*1 HEMI
!
        HEMI='S'
              selectcase(km_area_id)
              case(1)
                call ll_xy(lon,lat,x,y)
              case(2)
                call ll_xy_S(lon,lat,x,y)
              case(3)
                SLAT=71.;SLON=70.
                call mapll(1,1,lon,lat,x,y,SLAT,SLON,HEMI)
              case(4)
                SLAT=71.;SLON=-70.
                call mapll(1,1,lon,lat,x,y,SLAT,SLON,HEMI)
              case(5)
                SLAT=71.;SLON=150.
                call mapll(1,1,lon,lat,x,y,SLAT,SLON,HEMI)
              endselect
!
        return
        end
! translate x,y to lat,lon
! Uniform grid (in km) is centered in lat=90 (North Pole)
! Y oriented along lon=0
! Lana Erofeeva, May 2002
!
      subroutine xy_ll(x,y,lon,lat)
      real*8 y,x,lat,lon,pi
!
      pi=3.14159265358979
! 
      lat=90.-sqrt(x**2+y**2)/111.7
      lon=atan2(y,x)*180./pi
      if(lon.lt.0.)lon=lon+360.
      if(x.eq.0.and.y.eq.0.)lon=0.
      return
      end
!______________________________________________________________________
!
! translate lat,lon to x,y (km)
! Uniform grid (in km) is centered in lat=90 (North Pole)
! Y oriented along lon=0
! Lana Erofeeva, May 2002
      subroutine ll_xy(lon,lat,x,y)
      real*8 y,x,lat,lon,pi
!
      pi=3.14159265358979
! 
      x=(90.-lat)*111.7*cos(lon/180.*pi)
      y=(90.-lat)*111.7*sin(lon/180.*pi)
      return
      end
!
!----------------------------------------------------------------------
! translate x,y to lat,lon Antarctic
! Uniform grid (in km) is centered in lat=-90 (South Pole)
! Y oriented along lon=0
! Lana Erofeeva, May 2002
!
      subroutine xy_ll_S(x,y,lon,lat)
      real*8 y,x,lat,lon,pi
!
      pi=3.14159265358979
! 
      lat=-90.+sqrt(x**2+y**2)/111.7
      lon=-atan2(y,x)*180./pi+90.
      if(lon.gt.360.)lon=lon-360.
      if(lon.lt.  0.)lon=lon+360.
      return
      end
!______________________________________________________________________
!
! translate lat,lon to x,y (km) Antarctic
! Uniform grid (in km) is centered in lat=-90 (South Pole)
! Y oriented along lon=0
! Lana Erofeeva, May 2002
      subroutine ll_xy_S(lon,lat,x,y)
      real*8 y,x,lat,lon,pi
!
      pi=3.14159265358979
! 
      x=-(90.+lat)*111.7*cos((90+lon)/180.*pi)
      y= (90.+lat)*111.7*sin((90+lon)/180.*pi)
      return
      end
!___________________________________________________________________
!
	subroutine mapxy(n,m,x,y,lon,lat,SLAT,SLON,HEMI)
! (X,Y)km->(lon,lat)deg
        implicit none
        integer*4 n,m
        real*8 x(n,m),y(n,m),lat(n,m),lon(n,m),SLAT,SLON,SLAT1
        real*8 CDR,E2,E,pi,RE,SGN,CM,SL
        real*8 a1,a2,a3,a4,a5
        real*8, allocatable:: RHO(:,:),T(:,:),CHI(:,:)
        character*1 HEMI
!
        allocate(RHO(n,m),T(n,m),CHI(n,m))
!
        CDR= 57.29577951
        E2 = 6.694379852e-3           ! Eccentricity squared
        E  = sqrt(E2)
        pi = 3.14159265358979
        RE = 6378.273
!                                                                         
!*************************************************************************   
        if(HEMI.eq.'S'.or.HEMI.eq.'s')then
          SGN=-1
        else
          SGN= 1
        endif
        SLAT1=abs(SLAT)
        SL  = SLAT1/CDR
        RHO = sqrt(X**2+Y**2)
        if(n.eq.1.and.m.eq.1.and.RHO(1,1).lt.0.1)then
! Don't calculate if ONE POINT is on the equator
            lat=90.*SGN
            lon=0.0
            return
        else
           CM=cos(SL)/sqrt(1.0-E2*(sin(SL)**2)) 
           T=tan((pi/4.0)-(SL/2))/((1.0-E*sin(SL))/ &
                 (1.0+E*sin(SL)))**(E/2.0)
           if(abs(SLAT1-90.).lt.1.e-5)then
             T=RHO*sqrt((1.+E)**(1+E)*(1-E)**(1-E))/2/RE
           else
             T=RHO*T/(RE*CM)
           endif
           a1 =  5*E2**2 / 24
           a2 =    E2**3 / 12
           a3 =  7*E2**2 / 48
           a4 = 29*E2**3 /240
           a5 =  7*E2**3 /120
           CHI= (pi/2)-2*atan(T);
           lat= CHI+((E2/2) + a1 + a2)*sin(2*CHI)+(a3 + a4)*sin(4*CHI)+ &
                                  a5*sin(6*CHI)
           lat=SGN*lat*CDR
           lon= -(atan2(-X,Y)*CDR)+SLON
        endif
        deallocate(RHO,T,CHI)
        return
        end
!________________________________________________________________________
!
	subroutine mapll(n,m,lon,lat,x,y,SLAT,SLON,HEMI)
! (lon,lat)deg->(X,Y)km
        implicit none
        integer*4 n,m
        real*8 x(n,m),y(n,m),lat(n,m),lon(n,m),SLAT,SLON
        real*8 CDR,E2,E,pi,RE,SGN,RHO,SL,TC,MC
        character*1 HEMI
        real*8, allocatable:: T(:,:),lat1(:,:),lon1(:,:)
        allocate(T(n,m),lat1(n,m),lon1(n,m))
        CDR= 57.29577951
        E2 = 6.694379852e-3           ! Eccentricity squared
        E  = sqrt(E2)
        pi = 3.14159265358979
        RE = 6378.273
!                                                                         
!*************************************************************************   
        if(HEMI.eq.'S'.or.HEMI.eq.'s')then
          SGN=-1
        else
          SGN= 1
        endif
! This only works for Southern hemisphere lat/lon
!*************************************************************************
        if(abs(SLAT).eq.90)then
           RHO=2*RE/((1+E)**(1+E)*(1-E)**(1-E))**(E/2)
        else
           SL  = abs(SLAT)/CDR
           TC  = tan(pi/4-SL/2)/((1-E*sin(SL))/(1+E*sin(SL)))**(E/2)
           MC  = cos(SL)/sqrt(1-E2*(sin(SL)**2))
           RHO = RE*MC/TC;
        endif
        !lat1 = abs(lat)/CDR
         lat1= -lat/CDR
        T   = tan(pi/4-lat1/2)/((1-E*sin(lat1))/(1+E*sin(lat1)))**(E/2)
        lon1 =-(lon-SLON)/CDR
        x   =-RHO*T*sin(lon1)
        y   = RHO*T*cos(lon1)
        deallocate(T,lat1,lon1)
        return
        end
!-----------------------------------------------------------------ccc
! 
! loadModel:: reads binary z model file with the following format
!             fortran unformatted binary:
!
!       rec 1 (header):n,m,nc,theta_min,theta_max,phi_min,phi_max,
!                      const_1,const_2,...const_nc
! where const_j - constituent id char*4                     
!       rec 2:  1st constituent elevations (n by m complex)
!       rec 3:  2nd constituent elevations
!         .
!       rec nc+1: constituent n! elevations
!
!
	subroutine loadModel_z(model,fname)
        implicit none
        include 'derived_types.h'
        include 'constit.h'
        character fname*(*)
        type (z_model)model
        integer(kind=4) ic,j
        real(kind=4) xlim(2),ylim(2),dx
!
!--------------------------------------------------------------ccc
! reading the model
	open(unit=1,file=fname,form='unformatted',status='old',err=1)
	read(1,err=2)model%n,model%m,model%nc,ylim, &
                     xlim,(model%cid%cons(j),j=1,model%nc)
        model%lglobe=.false.
        dx=(xlim(2)-xlim(1))/model%n
        if(abs(xlim(2)-xlim(1)-360).le.dx)model%lglobe=.true. 
        do ic=1,model%nc
         read(1)model%z(ic,1:model%n,1:model%m)
        enddo
        model%mz=0
        where(model%z(1,:,:).ne.0)model%mz=1
        close(1)
        if(model%ll_km)then
          model%x_lim=xlim
          model%y_lim=ylim
        else
          model%phi_lim=xlim
          model%theta_lim=ylim
        endif
        return
1       write(*,*)'Model file ',trim(fname),' not found'
        stop
2       write(*,*)'Wrong format of model file ',trim(fname)
        stop
        end
!-----------------------------------------------------------------ccc
! 
! loadModel_uv:: reads binary uv model file with the following format
!                fortran unformatted binary:
!
!       rec 1 (header):n,m,nc,theta_min,theta_max,phi_min,phi_max,
!                      const_1,const_2,...const_nc
! where const_j - constituent id char*4                     
!       rec 2:  1st constituent transports (2 by n by m complex)
!       rec 3:  2nd constituent transports
!       ...     
!       rec nc+1: constituent nc transports (m2/s)
!
!
	subroutine loadModel_uv(model,fname)
        implicit none
        include 'derived_types.h'
        include 'constit.h'
        character fname*(*)
        complex*8, allocatable:: uv(:,:,:)
        integer(kind=4) ic,j
        real(kind=4) xlim(2),ylim(2),dx
        type (uv_model)model
!
!------------------------------------------------------------ccc
! reading the model
 	open(unit=1,file=fname,form='unformatted',status='old',err=1)
	read(1,err=2)model%n,model%m,model%nc,ylim, &
                     xlim,(model%cid%cons(j),j=1,model%nc)
        model%lglobe=.false.
        dx=(xlim(2)-xlim(1))/model%n
        if(abs(xlim(2)-xlim(1)-360).le.dx)model%lglobe=.true. 
        allocate(uv(2,model%n,model%m))
        do ic=1,model%nc
         read(1)uv
         model%u(ic,:,:)=uv(1,:,:)
         model%v(ic,:,:)=uv(2,:,:)
        enddo
        model%mu=0
        where(model%u(1,:,:).ne.0)model%mu=1
        model%mv=0
        where(model%v(1,:,:).ne.0)model%mv=1
        close(1)
        if(model%ll_km)then
          model%x_lim=xlim
          model%y_lim=ylim
        else
          model%phi_lim=xlim
          model%theta_lim=ylim
        endif
         deallocate(uv)
        return
1       write(*,*)'Model file ',trim(fname),' not found'
        stop
2       write(*,*)'Wrong format of model file ',trim(fname)
        stop
        end
!---------------------------------------------------------------------
        subroutine interp_driver_z(model,aloc,nrec,P,y)
! 
! interpolating z model on all data sites
! all info about the model is in model, grid - in grid,
! locations - aloc, result is in P(nc,nrec)
! y='y' - mask off land
        implicit none
        include 'derived_types.h'
        include 'constit.h'
        integer(kind=4)nrec,nloc,nc,nz
        type (z_model), intent(in):: model
        type (LocDict) aloc(nrec)
        type (c_grid), pointer:: mgrid
        type (sparsedf), pointer:: sdf
        complex(kind=4), intent(out):: P(model%nc,nrec)
        integer(kind=4) k,iGr,jGr,iloc,mGr,irec,ic
        complex(kind=8), allocatable:: Sc16(:)
        character*1 y
!
        allocate(Sc16(model%nc))
!
        P=0.
        allocate(mgrid)
        mgrid%NPG=model%NPG;mgrid%ll_km=model%ll_km;
        mgrid%lglobe=model%lglobe ! BUG Fix, Apr 2011
        mgrid%km_area_id=model%km_area_id
        mgrid%n=model%n;mgrid%m=model%m
        mgrid%theta_lim=model%theta_lim
        mgrid%phi_lim=model%phi_lim
        mgrid%x_lim=model%x_lim
        mgrid%y_lim=model%y_lim
        allocate(mgrid%mz(mgrid%n,mgrid%m))
        mgrid%mz=model%mz
        do k=1,nrec
         allocate(sdf)
         nloc=1
         sdf%ind_loc=k
         sdf%nloc=nloc
         allocate(sdf%iC(nloc,4),sdf%jC(nloc,4),sdf%ww(nloc,4))
         sdf%allocated=.true.
         call def_sparsedf(mgrid,aloc(k),sdf,'z')
         nz=0
         if(sum(sdf%iC).eq.0)then ! for out of grid data sites
           P(:,k)=0.
         else
         do iloc=1,nloc
           do ic=1,model%nc
            Sc16(ic)= model%z(ic,sdf%iC(iloc,1),sdf%jC(iloc,1))*sdf%ww(iloc,1)+ &
                      model%z(ic,sdf%iC(iloc,2),sdf%jC(iloc,2))*sdf%ww(iloc,2)+ &
                      model%z(ic,sdf%iC(iloc,3),sdf%jC(iloc,3))*sdf%ww(iloc,3)+ &
                      model%z(ic,sdf%iC(iloc,4),sdf%jC(iloc,4))*sdf%ww(iloc,4)
           enddo
           if(sum(Sc16).ne.0.)then
             P(:,k)=P(:,k)+Sc16
           else
             nz=nz+1
           endif
         enddo !iloc
         endif
         P(:,k)=P(:,k)/nloc
         if(nz.eq.nloc.and.y.eq.'y')then
            !write(*,*)minval(aloc(k)%lat),maxval(aloc(k)%lat),mgrid%theta_lim
            !write(*,*)minval(aloc(k)%lon),maxval(aloc(k)%lon),mgrid%phi_lim
            !write(*,*)nz,nloc
            !read(*,*)
            aloc(k)%land=.true.
         endif
         if(sum(P(:,k)).eq.0.and.y.eq.'y')aloc(k)%land=.true.
         deallocate(sdf%iC,sdf%jC,sdf%ww)     ! This is IMPORTANT
         deallocate(sdf)                      ! to deallocate sdf%iC etc.!
        enddo !k over data records
        deallocate(Sc16,mgrid%mz);deallocate(mgrid);
        return
        end
!---------------------------------------------------------------------
        subroutine interp_driver_d(mgrid,aloc,nrec,d,y)
! 
! interpolating grid depth on all data sites
! locations - aloc, result is in d(nrec)
! y='y' - mask off land
        implicit none
        include 'derived_types.h'
        include 'constit.h'
        integer(kind=4)nrec,nloc,nc,nz
        type (LocDict),intent(inout):: aloc(nrec)
        type (c_grid), intent(in):: mgrid
        type (sparsedf), pointer:: sdf
        real(kind=4), intent(out):: d(nrec)
        real(kind=8) d1
        integer(kind=4) k,iGr,jGr,iloc,mGr,irec,ic
        character*1 y
!
        d=0.
        do k=1,nrec
         allocate(sdf)
         nloc=1
         sdf%ind_loc=k
         sdf%nloc=nloc
         allocate(sdf%iC(nloc,4),sdf%jC(nloc,4),sdf%ww(nloc,4))
         sdf%allocated=.true.
         call def_sparsedf(mgrid,aloc(k),sdf,'z')
         nz=0
         if(sum(sdf%iC).gt.0)then ! for out of grid data sites
         do iloc=1,nloc
                 d1 = mgrid%hz(sdf%iC(iloc,1),sdf%jC(iloc,1))*sdf%ww(iloc,1)+ &
                      mgrid%hz(sdf%iC(iloc,2),sdf%jC(iloc,2))*sdf%ww(iloc,2)+ &
                      mgrid%hz(sdf%iC(iloc,3),sdf%jC(iloc,3))*sdf%ww(iloc,3)+ &
                      mgrid%hz(sdf%iC(iloc,4),sdf%jC(iloc,4))*sdf%ww(iloc,4)
           if(d1.gt.0)then
             d(k)=d(k)+d1
           else
             nz=nz+1
           endif
         enddo !iloc
         endif
         d(k)=d(k)/nloc
         if(d(k).eq.0.and.y.eq.'y')aloc(k)%land=.true.
         deallocate(sdf%iC,sdf%jC,sdf%ww) ! This is IMPORTANT
         deallocate(sdf)                  ! to deallocate sdf%iC etc.!
        enddo !k over data 
        return
        end
!----------------------------------------------------------------------------- !
        subroutine interp_driver_uv(model,aloc,grid,nrec,Puv,y)
!
! interpolating uv model on all data sites: to CURRENTS(m/s)
! all info about the model is in model, grid - in grid,
! locations - aloc, result is in Puv(nc,nrec)
! NOTE: grid is ALWAYS the model grid!
        implicit none
        include 'derived_types.h'
        include 'constit.h'
        integer(kind=4)nrec,nloc,ic
        type (uv_model)model
        type (LocDict) aloc(nrec)
        type (c_grid), intent(in):: grid
        type (c_grid), pointer:: mgrid
        type (sparsedf), pointer:: sdfz,sdfu,sdfv
        complex*8 Puv(model%nc,nrec)
        integer(kind=4) k,iGr,jGr,iloc,mGr,nz
        complex*16, allocatable:: Pu(:),Pv(:)
        complex*16 depth
        character*1 node,y ! if y=='T' -> do transports & mask land
                           !                otherwise do velocities
                           ! if y=='y' -> mask land
!
        allocate(mgrid)
      
! if model is on a different grid compared to the given
        !if(model%n.ne.grid%n.or.model%m.ne.grid%m.or. &
        !    model%theta_lim(1).ne.grid%theta_lim(1).or.model%theta_lim(2).ne.grid%theta_lim(2).or.&
        !    model%phi_lim(1).ne.grid%phi_lim(1).or.model%phi_lim(2).ne.grid%phi_lim(2))then
        !    write(*,*)''
        !    write(*,*)'FATAL ERROR: Wrong call of interp_driver_uv'
        !    write(*,*)'Interpolated model grid and default grid are not the same'
        !    write(*,*)'Currents can not be interpolated properly'
        !    write(*,*)model%n,grid%n,model%m,grid%m
        !    write(*,*)model%theta_lim(1),grid%theta_lim(1),model%theta_lim(2),grid%theta_lim(2)
        !    write(*,*)model%phi_lim(1),grid%phi_lim(1),model%phi_lim(2),grid%phi_lim(2)
        !    stop            
        ! endif

         mgrid%NPG=grid%NPG;mgrid%ll_km=grid%ll_km;
         mgrid%lglobe=grid%lglobe ! BUG Fix, Apr 2011
         mgrid%km_area_id=grid%km_area_id
         mgrid%n=grid%n;mgrid%m=grid%m
         mgrid%theta_lim=grid%theta_lim
         mgrid%phi_lim=grid%phi_lim
         mgrid%x_lim=grid%x_lim
         mgrid%y_lim=grid%y_lim
         allocate(mgrid%mz(mgrid%n,mgrid%m),mgrid%hz(mgrid%n,mgrid%m))
         mgrid%mz=grid%mz
         mgrid%hz=grid%hz
         allocate(Pu(model%nc),Pv(model%nc))
!
        Puv=0.
        do k=1,nrec
         Pu=0;Pv=0;
         allocate(sdfz,sdfu,sdfv)
         nloc=1
         sdfz%ind_loc=k
         sdfz%nloc=nloc
         allocate(sdfz%iC(nloc,4),sdfz%jC(nloc,4),sdfz%ww(nloc,4))
         sdfz%allocated=.true.
         if(aloc(k)%theta.eq.1.and.aloc(k)%phi.eq.0)then
           node='u'
         elseif(aloc(k)%theta.eq.0.and.aloc(k)%phi.eq.1)then
           node='v'
         else
           node='r'
         endif
         call def_sparsedf(mgrid,aloc(k),sdfz,'z')
         if(node.eq.'u'.or.node.eq.'r')then
          !allocate(sdfu)
          sdfu%ind_loc=k
          sdfu%nloc=nloc
          allocate(sdfu%iC(nloc,4),sdfu%jC(nloc,4),sdfu%ww(nloc,4))
          sdfu%allocated=.true.
          mgrid%mz=model%mu 
          call def_sparsedf(mgrid,aloc(k),sdfu,'u')
         endif
         if(node.eq.'v'.or.node.eq.'r')then
          !allocate(sdfv)
          sdfv%ind_loc=k
          sdfv%nloc=nloc
          allocate(sdfv%iC(nloc,4),sdfv%jC(nloc,4),sdfv%ww(nloc,4))
          sdfv%allocated=.true.
          mgrid%mz=model%mv
          call def_sparsedf(mgrid,aloc(k),sdfv,'v')
         endif
!
         nz=0
         do iloc=1,nloc
          if(y.eq.'T')then
           depth=1
          else
           depth = mgrid%hz(sdfz%iC(iloc,1),sdfz%jC(iloc,1))*sdfz%ww(iloc,1)+ &
                   mgrid%hz(sdfz%iC(iloc,2),sdfz%jC(iloc,2))*sdfz%ww(iloc,2)+ &
                   mgrid%hz(sdfz%iC(iloc,3),sdfz%jC(iloc,3))*sdfz%ww(iloc,3)+ &
                   mgrid%hz(sdfz%iC(iloc,4),sdfz%jC(iloc,4))*sdfz%ww(iloc,4)
          endif
          if(node.eq.'u'.or.node.eq.'r')then
            do ic=1,model%nc
             Pu(ic)= model%u(ic,sdfu%iC(iloc,1),sdfu%jC(iloc,1))*sdfu%ww(iloc,1)+ &
                     model%u(ic,sdfu%iC(iloc,2),sdfu%jC(iloc,2))*sdfu%ww(iloc,2)+ &
                     model%u(ic,sdfu%iC(iloc,3),sdfu%jC(iloc,3))*sdfu%ww(iloc,3)+ &
                     model%u(ic,sdfu%iC(iloc,4),sdfu%jC(iloc,4))*sdfu%ww(iloc,4)
            enddo
          endif
          if(node.eq.'v'.or.node.eq.'r')then
            do ic=1,model%nc
             Pv(ic)= model%v(ic,sdfv%iC(iloc,1),sdfv%jC(iloc,1))*sdfv%ww(iloc,1)+ &
                     model%v(ic,sdfv%iC(iloc,2),sdfv%jC(iloc,2))*sdfv%ww(iloc,2)+ &
                     model%v(ic,sdfv%iC(iloc,3),sdfv%jC(iloc,3))*sdfv%ww(iloc,3)+ &
                     model%v(ic,sdfv%iC(iloc,4),sdfv%jC(iloc,4))*sdfv%ww(iloc,4)
            enddo
          endif
          if(depth.ne.0.)then
            Puv(:,k)=Puv(:,k)+(Pu*aloc(k)%theta+Pv*aloc(k)%phi)/depth
          else
            nz=nz+1
          endif 
         enddo !iloc
         Puv(:,k)=Puv(:,k)/nloc
         if(abs(sum(Puv(:,k))).lt.1e-7)nz=nloc
         if(nz.eq.nloc.and.(y.eq.'y'.or.y.eq.'T'))aloc(k)%land=.true.
         deallocate(sdfz%iC,sdfz%jC,sdfz%ww);deallocate(sdfz);
         if(node.eq.'u'.or.node.eq.'r')deallocate(sdfu%iC,sdfu%jC,sdfu%ww)
         deallocate(sdfu)
         if(node.eq.'v'.or.node.eq.'r')deallocate(sdfv%iC,sdfv%jC,sdfv%ww)
         deallocate(sdfv)
        enddo !k over data records
        deallocate(Pu,Pv,mgrid%mz,mgrid%hz)
        deallocate(mgrid)
        return
        end
! -------------------------------------------------------------------------------- !
        subroutine def_sparsedf(grid,aloc,sdf,node)
!  Input:  grid, data locations "aloc", node ("z","u","v")
!  Output: sparse data functional "sdf"
        implicit none
        include 'derived_types.h'
        type (c_grid) grid
        type (LocDict) aloc
        type (sparsedf) sdf
        character*1 node
        real (kind=4)step_lat,step_lon,lat1,lon1
        real (kind=4)lat0,lon0,gth_lim(2),gph_lim(2)
        real (kind=8) dlat,dlon,tx,ty,dx,dy,xi,xj,x,y,w00,w01,w10,w11,wtot, &
                      SLAT,SLON,sm,tlat,tlon,blat,blon
        integer (kind=4) iGr,jGr,iloc,mGr,nloc,i0,j0,i1,j1,ipshft
        character*1 HEMI
!
        if(.not.sdf%allocated)then
         write(*,*)'Wrong call of def_sparsedf:'
         write(*,*)'sdf is not allocated'
         stop
        endif
        if(node.eq.'U')node='u'
        if(node.eq.'V')node='v'
        if(node.eq.'Z')node='z'
        if(node.ne.'u'.and.node.ne.'v'.and.node.ne.'z')then
          write(*,*)'Wrong call: def_sparsedf, no ',node,' nodes in C-grid'
          stop
        endif
        if(node.eq.'z')sdf%itype=1
        if(node.eq.'u')sdf%itype=2
        if(node.eq.'v')sdf%itype=3
        if(aloc%theta.ne.0.and.aloc%phi.ne.0)sdf%itype=4
!
 ! in majority of cases nloc=1, step*=0. Only needed for integral cases, like GRACE
        nloc=sdf%nloc
        dy = (grid%theta_lim(2)-grid%theta_lim(1))/grid%m
        dx = (grid%phi_lim(2)-grid%phi_lim(1))/grid%n
        !call find_step_ll(aloc,nloc,dx,dy,step_lat,step_lon,lat1,lon1,mGr)
        step_lat=0;step_lon=0;lat1=aloc%lat(1);lon1=aloc%lon(1)
        if(grid%ll_km)then
         gth_lim=grid%y_lim
         gph_lim=grid%x_lim
         dy = (grid%y_lim(2)-grid%y_lim(1))/grid%m
         dx = (grid%x_lim(2)-grid%x_lim(1))/grid%n
        else
         gth_lim=grid%theta_lim
         gph_lim=grid%phi_lim
         !dy = (grid%theta_lim(2)-grid%theta_lim(1))/grid%m
         !dx = (grid%phi_lim(2)-grid%phi_lim(1))/grid%n
        endif

        iGr=1;jGr=1
! these arrays should be already allocated in calling routine !!!
        sdf%iC=0;sdf%jC=0;sdf%ww=0;
        do iloc=1,nloc
           if(aloc%lat(3)*aloc%lat(4)*aloc%lon(3)*aloc%lon(4).ne.0) then
! GRACE: integral over rectangular: just passing each cell
              lon0=lon1+iGr*step_lon-step_lon/2
              lat0=lat1+jGr*step_lat-step_lat/2
              !write(*,*)iloc,iGr,jGr,lat0,lon0              
              jGr=jGr+1
             if(jGr.gt.mGr)then
              iGr=iGr+1
              jGr=1
              !read(*,*)
             endif 
           else ! everything else (here nloc>1 only for Tav case)
            lon0=lon1+step_lon*(iloc-1)
            lat0=lat1+step_lat*(iloc-1)
           endif
! take care of rotated coordinates & stereographic grids
           tlat=lat0;tlon=lon0
           !if(grid%NPG)then
           ! call BTRb(tlon,tlat,blon,blat,1,1)
           ! tlat=blat
           ! tlon=blon
           !endif
           if(grid%ll_km)then
              dlat=tlat
              dlon=tlon
              call convert_ll_xy(grid%km_area_id,dlon,dlat,tx,ty)
              tlon=tx
              tlat=ty
           else ! adjust lon convention if needed
            if(tlon.lt.gph_lim(1))tlon=tlon+360.
            if(tlon.gt.gph_lim(2))tlon=tlon-360.
           endif
            if(node.eq.'z')then
!      find grid coordinates for lower left hand z node
!      among the square (of z nodes which surround the obs point)
!
!                 |                     
!                \|/
!     -> -z---u---z---u---z--
!         |       |       |     x is the location given by theta,phi
!         v       v       v     the four surrounding elevation
!         | x     |       |     nodes are used for the bilinear
!     -> -z---u---z---u---z-    spline interpolation
!         |      /|\      |
!                 |
            xi = (tlon-gph_lim(1))/dx+.5
            xj = (tlat-gth_lim(1))/dy+.5
         elseif(node.eq.'u')then
!       interior point: measurement of current vector
!                       direction is given as unit vector (th,ph)
!                       (th,ph)=(1.,0.)->u
!                       (th,ph)=(0.,1.)->v 

!      find grid coordinates for lower left hand z node
!      among the square (of z nodes which surround the obs point)
!               
!            
!             u---z---u--
!             |       |         x is the location given by theta,phi
!             |   v---|---v     the 8 surrounding u,v 
!             |   | x |   |     nodes are used for the bilinear
!             u---z---u---z-    spline interpolation
!                 |       |     Current vector direction is given
!                 v-------v     with a unit vector (th,ph):
!                               (1.,0)->EW, (0.,1.)->NS
            xi = (tlon-gph_lim(1))/dx+1.0
            xj = (tlat-gth_lim(1))/dy+.5
         elseif(node.eq.'v')then
            xi = (tlon-gph_lim(1))/dx+.5
            xj = (tlat-gth_lim(1))/dy+1.0
         endif
         if(xi.lt.1.and.grid%lglobe) xi = grid%n+xi
!
         i0 = int(xi)
         x = xi-i0
         j0 = int(xj)
         y = xj-j0
!      check to see if calculated indices are in range
         if((i0.gt.grid%n).or.(i0.lt.1).or.(j0.gt.grid%m).or.(j0.lt.1) )then
          sdf%ic=1;sdf%jc=1
          !write(*,*)
          !write(*,*)aloc%lon(1),tlon,lon0
          !write(*,*)aloc%lat(1),tlat,lat0
          !write(*,*)i0,grid%n,j0,grid%m
          return
         endif
!
!        compute weights for bilinear spline interpolation; only
!        use ocean nodes (for which mask mask is = 1)
         j1 = ipshft(j0,1,grid%m)
         i1 = ipshft(i0,1,grid%n)
         sm=grid%mz(i0,j0)+grid%mz(i0,j1)+grid%mz(i1,j0)+grid%mz(i1,j1)
         if(sm.gt.0)then
          w00 = (1.-x)*(1.-y)*grid%mz(i0,j0)
          w01 = (1.-x)*y*grid%mz(i0,j1)
          w10 = x*(1.-y)*grid%mz(i1,j0)
          w11 = x*y*grid%mz(i1,j1)
          wtot = w00+w01+w10+w11
          if(wtot.eq.0)return
!
          sdf%ww(iloc,1) = w00/wtot
          sdf%ww(iloc,2) = w01/wtot
          sdf%ww(iloc,3) = w10/wtot
          sdf%ww(iloc,4) = w11/wtot
         endif
!
         sdf%iC(iloc,1)=i0;sdf%iC(iloc,2)=i0;sdf%iC(iloc,3)=i1;sdf%iC(iloc,4)=i1;
         sdf%jC(iloc,1)=j0;sdf%jC(iloc,2)=j1;sdf%jC(iloc,3)=j0;sdf%jC(iloc,4)=j1;
        enddo !iloc
        return
        end
!-----------------------------------------------------------------------------
      subroutine sdf_set(rloc,nrep,grid,MASK,sdfz,sdfu,sdfv)
      implicit none
      include 'derived_types.h'
      type (LocDict), intent(in):: rloc(nrep)
      type (c_grid), intent(in):: grid
      type (c_grid), pointer:: mgrid
      type (masks), intent(in):: MASK
      type (sparsedf), intent(inout):: sdfz(nrep),sdfu(nrep),sdfv(nrep)
      integer(kind=4)i,nrep
!
      allocate(mgrid)
      mgrid%NPG=grid%NPG;mgrid%ll_km=grid%ll_km;
      mgrid%km_area_id=grid%km_area_id
      mgrid%n=grid%n;mgrid%m=grid%m
      mgrid%theta_lim=grid%theta_lim
      mgrid%phi_lim=grid%phi_lim
      mgrid%x_lim=grid%x_lim
      mgrid%y_lim=grid%y_lim
      allocate(mgrid%mz(mgrid%n,mgrid%m))
      do i=1,nrep
         call def_sparsedf(grid,rloc(i),sdfz(i),'z')
         mgrid%mz=MASK%mu
         call def_sparsedf(mgrid,rloc(i),sdfu(i),'u')
         mgrid%mz=MASK%mv
         call def_sparsedf(mgrid,rloc(i),sdfv(i),'v')
      enddo
      deallocate(mgrid)
      return
      end
!----------------------------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! reads in all grid arrays ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine atlas_mask_in(gname,pmask,n,m,atlas,nmod,lname,ln)
      implicit none
      character*80 gname
      integer*4 n,m,k,nob,pmask(n,m),nmod,eof,n1,m1,nz,ln
      real*4 theta_lim(2),phi_lim(2),dt,ll_loc(4)
      integer(kind=4), allocatable:: i(:),j(:)
      real*4, allocatable:: d(:)
      logical atlas
      character*20 lname(ln)
!
      atlas=.false.
      open(unit=1,file=trim(gname),form='unformatted',status='old')
      call pass_record_grid(1,0)
      read(1,end=1,err=1)pmask
      atlas=.true.
! now define # of local solutions
      eof=0;k=1;nmod=1
      do while(eof.eq.0)
       if(ln.gt.1)k=nmod
       read(1,end=2,err=2)n1,m1,nz,ll_loc,lname(k)
       allocate(d(nz),i(nz),j(nz))
       read(1)i,j
       read(1)d
       deallocate(d,i,j)
       if(eof.eq.0)nmod=nmod+1
      enddo
!
2     nmod=nmod-1
      close(1)
      return
1     pmask=1
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
! reading dimensions of local grid (with ID==id) from atlas grid
      subroutine getDims_l(cgrid,id,grid,nz,iunit)
      implicit none
      include 'derived_types.h'
      integer(kind=4) j,id,nz,iunit
      character cgrid*(*)
      type (c_grid) grid
!
      open(unit=iunit,file=trim(cgrid),form='unformatted',status='old')
! pass records for global grid (including pmask)
      call pass_record_grid(iunit,1)
      do j=1,id-1
       call pass_record_lgrid(iunit)
      enddo
      call rd_lgrid_header(grid,iunit,nz)
      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
! reading dimensions of local model (with ID==id) from atlas
      subroutine getDims_z_l(fname,id,nz,iunit,model)
      implicit none
      include 'derived_types.h'
      integer(kind=4) j,id,nz,iunit
      character fname*(*)
      type (z_model) model
!
      open(unit=iunit,file=trim(fname),form='unformatted',status='old')
! pass records for global grid
      call pass_record_z(iunit)
      do j=1,id-1
       call pass_record_lz(iunit)
      enddo
      call rd_lz_header(model,iunit,nz)
      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
! reading dimensions of local model (with ID==id) from atlas
      subroutine getDims_u_l(fname,id,nu,nv,iunit,model)
      implicit none
      include 'derived_types.h'
      integer(kind=4) j,id,nu,nv,iunit
      character fname*(*)
      type (uv_model) model
!
      open(unit=iunit,file=trim(fname),form='unformatted',status='old')
! pass records for global grid 
      call pass_record_u(iunit)
      do j=1,id-1
       call pass_record_lu(iunit)
      enddo
      call rd_lu_header(model,iunit,nu,nv)
      return
      end
! reading local grid from atlas grid
! (correct record is already set by getDims_l, close iunit file when done)
      subroutine grid_in_l(iunit,grid,nz)
      implicit none
      include 'derived_types.h'
      integer(kind=4)nz,iunit,k
      integer(kind=4),allocatable:: i(:),j(:)
      type (c_grid) grid
      real(4),allocatable::depth(:)
!
      if(.not.grid%allocated)then
       write(*,*)'wrong usage og grid_in_l: gris is not allocated'
       stop
      endif
      allocate(depth(nz),i(nz),j(nz))
!
      read(iunit)i,j
      read(iunit)depth
      grid%mz=0;grid%hz=0
      do k=1,nz
       grid%mz(i(k),j(k))=1;grid%hz(i(k),j(k))=depth(k);
      enddo
      deallocate(depth,i,j)
      close(iunit)
      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rd_lgrid_header(grid,iunit,nz)
      implicit none
      include 'derived_types.h'
      integer(kind=4) iunit,nz
      real(kind=4) eps
      type (c_grid) grid
      character*20 lname
!
      read(iunit,err=1) grid%n,grid%m,nz,grid%theta_lim,grid%phi_lim,lname
      grid%dt=1;grid%nob=0;
      grid%ll_km =.false.
!cc Check for global (i.e. CATs)
      eps=0.25*(grid%phi_lim(2)-grid%phi_lim(1))/grid%n
      if(abs(grid%phi_lim(2)-grid%phi_lim(1)-360.).lt.eps)then
        grid%lglobe=.true.
      else
        grid%lglobe=.false.
      endif 
      return
1     grid%n=0
      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rd_lz_header(model,iunit,nz)
      implicit none
      include 'derived_types.h'
      integer(kind=4) iunit,nz
      real(kind=4) eps
      type (z_model) model
      character*20 lname
!
      read(iunit,err=1) model%n,model%m,model%nc,nz,model%theta_lim,model%phi_lim
      model%ll_km=.false.
!cc Check for global (i.e. CATs)
      eps=0.25*(model%phi_lim(2)-model%phi_lim(1))/model%n
      if(abs(model%phi_lim(2)-model%phi_lim(1)-360.).lt.eps)then
        model%lglobe=.true.
      else
        model%lglobe=.false.
      endif 
      return
1     model%n=0
      return
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rd_lu_header(model,iunit,nu,nv)
      implicit none
      include 'derived_types.h'
      integer(kind=4) iunit,nu,nv
      real(kind=4) eps
      type (uv_model) model
      character*20 lname
!
      read(iunit,err=1) model%n,model%m,model%nc,nu,nv,model%theta_lim,model%phi_lim
      model%ll_km=.false.
!cc Check for global (i.e. CATs)
      eps=0.25*(model%phi_lim(2)-model%phi_lim(1))/model%n
      if(abs(model%phi_lim(2)-model%phi_lim(1)-360.).lt.eps)then
        model%lglobe=.true.
      else
        model%lglobe=.false.
      endif 
      return
1     model%n=0
      return
      end
!--------------------------------------------------------------------------------------------------
     subroutine pass_record_grid(iunit,imask)
     implicit none
     integer(kind=4) iunit,n,m,nob,k,j,imask
     real(kind=4)theta_lim(2),phi_lim(2),dt
     integer(kind=4),allocatable:: mtmp(:,:),iob(:,:)
     real(kind=4), allocatable:: hz(:,:),pmask(:,:)
! pass header: header has been read before grid was allocated
      read(iunit)n,m,theta_lim,phi_lim,dt,nob
      allocate(mtmp(n,m),hz(n,m),iob(2,nob))
!     read open boundary conditions
      if(nob.gt.0)then
        read(iunit) ((iob(k,j),k=1,2),j=1,nob)
      else
        read(iunit)
      endif
!     read bathymetry for z nodes
      read(iunit) hz
!     read in z node mask
      read(iunit) mtmp
      if(imask.eq.1)then
       allocate(pmask(n,m))
       read(iunit) pmask
       deallocate(pmask)
      endif
      deallocate(mtmp,iob,hz)
      return
      end
!--------------------------------------------------------------------------------------------------
     subroutine pass_record_z(iunit)
     implicit none
     integer(kind=4) iunit,n,m,nc,k
     complex, allocatable:: z(:,:)
!
     read(iunit) n,m,nc ! pass the rest of the header
     allocate(z(n,m))
     do k=1,nc
      read(iunit)z
     enddo
     deallocate(z)
     return
     end
!--------------------------------------------------------------------------------------------------
     subroutine pass_record_u(iunit)
     implicit none
     integer(kind=4) iunit,n,m,nc,k
     complex, allocatable::uv(:,:,:)
!
     read(iunit) n,m,nc ! pass the rest of the header
     allocate(uv(2,n,m))
     do k=1,nc
      read(iunit)uv
     enddo
     deallocate(uv)
     return
     end
!--------------------------------------------------------------------------------------------------
     subroutine pass_record_lgrid(iunit)
     implicit none
     integer(kind=4)n,m,nz,iunit
     integer(kind=4),allocatable:: i(:),j(:)
     real(kind=4), allocatable:: d(:)
     read(iunit)n,m,nz ! pass the rest of header
     allocate(i(nz),j(nz),d(nz))
     read(iunit)i,j
     read(iunit)d
     deallocate(i,j,d)
     return
     end
!--------------------------------------------------------------------------------------------------
     subroutine pass_record_lz(iunit)
     implicit none
     integer(kind=4)n,m,nc,nz,iunit,ic
     integer(kind=4),allocatable:: i(:),j(:)
     complex, allocatable:: z(:)
     read(iunit)n,m,nc,nz ! pass the rest of header
     allocate(i(nz),j(nz),z(nz))
     read(iunit)i,j
     do ic=1,nc
      read(iunit)z
     enddo
     deallocate(i,j,z)
     return
     end
!--------------------------------------------------------------------------------------------------
     subroutine pass_record_lu(iunit)
     implicit none
     integer(kind=4)n,m,nc,nu,nv,iunit,ic
     integer(kind=4),allocatable:: iu(:),ju(:),iv(:),jv(:)
     complex, allocatable:: u(:),v(:)
     read(iunit)n,m,nc,nu,nv ! pass the rest of header
     allocate(iu(nu),ju(nu),u(nu))
     allocate(iv(nv),jv(nv),v(nv))
     read(iunit)iu,ju
     read(iunit)iv,jv
     do ic=1,nc
      read(iunit)u
      read(iunit)v
     enddo
     deallocate(iu,ju,iv,jv,u,v)
     return
     end
!-------------------------------------------------------------
	subroutine loadModel_z_l(iunit,model,nz)
        implicit none
        include 'derived_types.h'
        include 'constit.h'
        type (z_model)model
        integer(kind=4) ic,nz,k,iunit,i,j
        real(kind=4) xlim(2),ylim(2)
        integer(kind=4), allocatable:: iz(:),jz(:)
        complex(kind=4), allocatable::z1(:)
        character*20 lname
!
!--------------------------------------------------------------
! reading local model (header has been already read)
        backspace(iunit) ! one record back
	read(iunit)model%n,model%m,model%nc,nz,ylim, &
                     xlim,(model%cid%cons(j),j=1,model%nc),lname
        !write(*,*)lname
        !write(*,*)model%n,model%m,model%nc,nz,ylim,xlim
        !write(*,*)(model%cid%cons(j),j=1,model%nc)
        allocate(iz(nz),jz(nz),z1(nz))
        read(iunit)iz,jz
        model%z=0
        do ic=1,model%nc
         read(iunit)z1
         !write(*,*)ic,maxval(abs(z1))
         do k=1,nz
          model%z(ic,iz(k),jz(k))=z1(k)
         enddo
        enddo
        model%mz=0
        where(model%z(1,:,:).ne.0)model%mz=1
        close(iunit)
        if(model%ll_km)then
          model%x_lim=xlim
          model%y_lim=ylim
        else
          model%phi_lim=xlim
          model%theta_lim=ylim
        endif
        return
        end
!-------------------------------------------------------------
	subroutine loadModel_uv_l(iunit,model,nu,nv)
        implicit none
        include 'derived_types.h'
        include 'constit.h'
        type (uv_model)model
        integer(kind=4) ic,nu,nv,k,iunit,j
        real(kind=4) xlim(2),ylim(2)
        integer(kind=4), allocatable:: iu(:),ju(:),iv(:),jv(:)
        complex, allocatable::u1(:),v1(:)
        character*20 lname
!
!--------------------------------------------------------------
! reading local model
        backspace(iunit) ! one record back
	read(iunit)model%n,model%m,model%nc,nu,nv,ylim, &
                     xlim,(model%cid%cons(j),j=1,model%nc),lname
        !write(*,*)lname
        !write(*,*)model%n,model%m,model%nc,nu,nv,ylim,xlim
        !write(*,*)(model%cid%cons(j),j=1,model%nc)
        allocate(iu(nu),ju(nu),iv(nv),jv(nv),u1(nu),v1(nv))
        read(iunit)iu,ju
        read(iunit)iv,jv
        model%u=0;model%v=0
        do ic=1,model%nc
         read(iunit)u1
         read(iunit)v1
         do k=1,nu
          model%u(ic,iu(k),ju(k))=u1(k)
         enddo
         do k=1,nv
          model%v(ic,iv(k),jv(k))=v1(k)
         enddo
        enddo
        model%mu=0
        where(model%u(1,:,:).ne.0)model%mu=1
        model%mv=0
        where(model%v(1,:,:).ne.0)model%mv=1
        close(iunit)
        if(model%ll_km)then
          model%x_lim=xlim
          model%y_lim=ylim
        else
          model%phi_lim=xlim
          model%theta_lim=ylim
        endif
        return
        end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rd_lmsetup(modname,lth_lim,lph_lim,c_id,ncon,outname)
      implicit none
      include 'constit.h'
      character*80 modname,outname(4),tmp,rmCom,dirp,sfx
      character*4 c_id(ncmx),con
      real*4 lth_lim(2),lph_lim(2)
      integer k,ic,ncon
!
      read(*,'(a80)',err=1)modname
      modname=rmCom(modname)
      if(trim(modname).eq.'')then
       write(*,*)'1st line in setup file can not be empty'
       stop
      endif
      read(*,*,err=1)lth_lim
      read(*,*,err=1)lph_lim
!
      read(*,'(a80)',err=1)tmp
      tmp=rmCom(tmp)
      ic=0
      k=1
      do while(k.gt.0)
       ncon=ic
       k=index(tmp,',')
       ic=ic+1
       c_id(ic)=tmp(1:max(1,k-1))
       tmp=tmp(k+1:80)
      enddo
      if(trim(tmp).ne.'')ncon=ic
      c_id(ic)=tmp(1:80)
!
      read(*,'(a80)',err=1)dirp
      dirp=rmCom(dirp)
      read(*,'(a80)',err=1)sfx
      sfx=rmCom(sfx)
      outname(1)=trim(dirp)//'/grid_'//trim(sfx)
      outname(2)=trim(dirp)//'/h_'//trim(sfx)
      outname(3)=trim(dirp)//'/UV_'//trim(sfx)
      outname(4)='Model_'//trim(sfx)      
      return
1     write(*,*)'Input file of WRONG FORMAT: see lmsetup.example'
      write(*,*)'for right format'
      stop
      end

!----------------------------------------------------------------------------- !
        subroutine interp_driver_uvq(model,aloc,nrec,Puv,node,y)
!
! interpolating uv model on all data sites for U/V transports only
! quick version of interp_driver_uv (one kind of u/v data)
! all info about the model is in model, grid - in grid,
! locations - aloc, result is in Puv(nc,nrec)
! NOTE: grid is ALWAYS the model grid!
        implicit none
        include 'derived_types.h'
        include 'constit.h'
        integer(kind=4)nrec,nloc,ic
        type (uv_model)model
        type (LocDict) aloc(nrec)
        type (c_grid), pointer:: mgrid
        type (sparsedf), pointer:: sdfu
        complex(kind=4), intent(out):: Puv(model%nc,nrec)
        complex(kind=8), allocatable:: Pu(:)
        integer(kind=4) k,iGr,jGr,iloc,mGr,nz
        character*1 node,y ! 
                           ! if y=='y' -> mask land
!
        Puv=0
        allocate(Pu(model%nc))
        allocate(mgrid)
        mgrid%NPG=model%NPG;mgrid%ll_km=model%ll_km;
        mgrid%lglobe=model%lglobe ! BUG Fix, Apr 2011
        mgrid%km_area_id=model%km_area_id
        mgrid%n=model%n;mgrid%m=model%m
        mgrid%theta_lim=model%theta_lim
        mgrid%phi_lim=model%phi_lim
        mgrid%x_lim=model%x_lim
        mgrid%y_lim=model%y_lim
        allocate(mgrid%mz(mgrid%n,mgrid%m),mgrid%hz(mgrid%n,mgrid%m))
        if(node.eq.'U')node='u'
        if(node.eq.'V')node='v'
        if(node.eq.'u')mgrid%mz=model%mu
        if(node.eq.'v')mgrid%mz=model%mv
!
        do k=1,nrec
         allocate(sdfu)
         nloc=1
         sdfu%ind_loc=k
         sdfu%nloc=nloc
         allocate(sdfu%iC(nloc,4),sdfu%jC(nloc,4),sdfu%ww(nloc,4))
         sdfu%allocated=.true.
         call def_sparsedf(mgrid,aloc(k),sdfu,node)
!
         nz=0
         if(sum(sdfu%iC).ne.0)then ! for in-grid data sites
          do iloc=1,nloc
            do ic=1,model%nc
             if(node.eq.'u')then
             Pu(ic)= model%u(ic,sdfu%iC(iloc,1),sdfu%jC(iloc,1))*sdfu%ww(iloc,1)+ &
                     model%u(ic,sdfu%iC(iloc,2),sdfu%jC(iloc,2))*sdfu%ww(iloc,2)+ &
                     model%u(ic,sdfu%iC(iloc,3),sdfu%jC(iloc,3))*sdfu%ww(iloc,3)+ &
                     model%u(ic,sdfu%iC(iloc,4),sdfu%jC(iloc,4))*sdfu%ww(iloc,4)
             else
             Pu(ic)= model%v(ic,sdfu%iC(iloc,1),sdfu%jC(iloc,1))*sdfu%ww(iloc,1)+ &
                     model%v(ic,sdfu%iC(iloc,2),sdfu%jC(iloc,2))*sdfu%ww(iloc,2)+ &
                     model%v(ic,sdfu%iC(iloc,3),sdfu%jC(iloc,3))*sdfu%ww(iloc,3)+ &
                     model%v(ic,sdfu%iC(iloc,4),sdfu%jC(iloc,4))*sdfu%ww(iloc,4)
             endif
            enddo
            if(sum(Pu).ne.0.)then
             Puv(:,k)=Puv(:,k)+Pu(:)
            else
             nz=nz+1
            endif 
          enddo !iloc
         endif
         Puv(:,k)=Puv(:,k)/nloc
         if(nz.eq.nloc)aloc(k)%land=.true.
         deallocate(sdfu%iC,sdfu%jC,sdfu%ww)
         deallocate(sdfu)
        enddo !k over data records
        deallocate(Pu,mgrid%mz,mgrid%hz)
        deallocate(mgrid)
        return
        end
!------------------------------------------------------------------------------------
        subroutine write_grid(dloc,ndat,depth,outname,nl,ml,th_lim,ph_lim)
        implicit none
        include 'derived_types.h'
        type (LocDict), intent(in):: dloc(ndat)
        real(kind=4) depth(ndat),th_lim(2),ph_lim(2)
        integer(kind=4) nl,ml,idum,i,j,k,ndat
        integer(kind=4), allocatable:: mz(:,:)
        real(kind=4), allocatable:: hz(:,:)
        character*80 outname
! writing grid in OTIS binary format
         allocate(mz(nl,ml),hz(nl,ml))
         idum=0
         k=1
         do i=1,nl
          do j=1,ml
            hz(i,j)=depth(k)
            k=k+1
          enddo
         enddo
         open(unit=11,file = trim(outname), &
            form='unformatted',status='unknown')
         write(11) nl,ml,th_lim,ph_lim,12.,idum
         write(11)idum
         write(11)hz
         mz=0
         where(hz.gt.0)mz=1
         write(11) ((mz(i,j)   ,i=1,nl),j=1,ml)
         close(11)
         deallocate(hz,mz)        
         return
         end
!------------------------------------------------------------------------------------
        subroutine write_z(dloc,ndat,z1,depth,outname,nl,ml,th_lim,ph_lim, &
                   nc,ncon,cind,c_id)
        implicit none
        include 'derived_types.h'
        type (LocDict), intent(in):: dloc(ndat)
        real(kind=4), intent(in):: th_lim(2),ph_lim(2),depth(ndat)
        complex(kind=4) z1(nc,ndat)
        integer(kind=4) nl,ml,i,j,k,ic,ndat,ncon,nc,cind(ncon)
        complex(kind=4), allocatable:: z(:,:,:)
        character*80 outname
        character*4 c_id(ncon)
! writing grid in OTIS binary format
         allocate(z(ncon,nl,ml))
         z=0.
         do ic=1,ncon
          k=1
          do i=1,nl
           do j=1,ml
             if(depth(k).ne.0)z(ic,i,j)=z1(cind(ic),k)
             k=k+1
           enddo
          enddo
         enddo
         open(unit=11,file = trim(outname), &
            form='unformatted',status='unknown')
         write(11) nl,ml,ncon,th_lim,ph_lim,c_id
         do ic=1,ncon
          write(11)z(ic,:,:)
         enddo
         close(11)
         deallocate(z)        
         return
         end
!------------------------------------------------------------------------------------
        subroutine write_uv(dloc,ndat,u1,v1,depth,outname,nl,ml, &
                           th_lim,ph_lim,nc,ncon,cind,c_id)
        implicit none
        include 'derived_types.h'
        type (LocDict), intent(in):: dloc(ndat)
        real(kind=4), intent(in):: th_lim(2),ph_lim(2),depth(ndat)
        complex(kind=4) u1(nc,ndat),v1(nc,ndat)
        integer(kind=4) nl,ml,i,j,k,ic,ndat,ncon,nc,cind(ncon),i1,j1,m3(3,3)
        integer(kind=4) nh,niter,iter
        complex(kind=4), allocatable:: u(:,:,:),v(:,:,:),uv(:,:,:),us(:)
        integer(kind=4), allocatable:: mz(:,:),mu(:,:),mv(:,:)
        complex(kind=4) tmp(3,3)
        character*80 outname
        character*4 c_id(ncon)
! writing grid in OTIS binary format
         allocate(u(ncon,nl,ml),v(ncon,nl,ml),uv(2,nl,ml),us(ncon))
         allocate(mz(nl,ml),mu(nl,ml),mv(nl,ml))
         mz=0;mu=0;mv=0
         k=1
         do i=1,nl
          do j=1,ml
           if(depth(k).ne.0)mz(i,j)=1
           k=k+1
          enddo
         enddo
         call Muv(mz,mu,mv,nl,ml)
!
         u=0;v=0
         do ic=1,ncon
          k=1
          do i=1,nl
           do j=1,ml
             if(mu(i,j).ne.0)u(ic,i,j)=u1(cind(ic),k)
             if(mv(i,j).ne.0)v(ic,i,j)=v1(cind(ic),k)
             k=k+1
           enddo
          enddo
         enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          nh=1;niter=20;iter=1
          do while(nh.gt.0.and.iter.lt.niter)
          nh=0
          do i=2,nl-1
           do j=2,ml-1
             if(mu(i,j).ne.0.and.u(1,i,j).eq.0)then
                do ic=1,ncon
                 m3=0;tmp=u(ic,i-1:i+1,j-1:j+1)
                 where(tmp.ne.0)m3=1
                 us(ic)=sum(tmp)/max(sum(m3),1)
                enddo
                u(:,i,j)=us
                nh=1
             endif
             if(mv(i,j).ne.0.and.v(1,i,j).eq.0)then
                do ic=1,ncon
                 m3=0;tmp=v(ic,i-1:i+1,j-1:j+1)
                 where(tmp.ne.0)m3=1
                 us(ic)=sum(tmp)/max(sum(m3),1)
                enddo
                v(:,i,j)=us
                nh=1
             endif
           enddo
          enddo
          iter=iter+1
          enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         open(unit=11,file = trim(outname), &
            form='unformatted',status='unknown')
         write(11) nl,ml,ncon,th_lim,ph_lim,c_id
         do ic=1,ncon
          uv(1,:,:)=u(ic,:,:)
          uv(2,:,:)=v(ic,:,:)
          write(11)uv
         enddo
         close(11)
         deallocate(u,v,uv,mz,mu,mv,us)        
         return
         end
!-------------------------------------------------------------------------
      subroutine Muv(mz,mu,mv,n,m)
      implicit none
! mz->mu,mv for global
      integer(kind=4) n,m,mz(n,m),mu(n,m),mv(n,m)

       mu(:,:) = mz(:,:)*cshift(mz(:,:),dim=1,shift=-1)
       mv(:,:) = mz(:,:)*cshift(mz(:,:),dim=2,shift=-1)
!
       mv(:,1)=mz(:,1)
       return
       end

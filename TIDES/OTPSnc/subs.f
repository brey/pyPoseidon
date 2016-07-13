c=====================================================================
c AUTHORS:
c  Gary Egbert & Lana Erofeeva
c  College of Atmospheric and Oceanic Sciences
c  104 COAS Admin. Bldg.
c  Oregon State University
c  Corvallis, OR 97331-5503
c  
c  E-mail:  egbert@coas.oregonstate.edu                                      
c  Fax:     (541) 737-2064
c  Ph.:     (541) 737-2947                                        
c  http://volkov.oce.orst.edu/tides/
c
c COPYRIGHT: OREGON STATE UNIVERSITY, 2012
c (see the file COPYRIGHT for lisence agreement)
c=====================================================================
c Set of subroutines used by extract_HC and predict_tide
c 
      subroutine rd_mod_body_nc(modname,zuv,cind,nc,ncon,
     *                          mlon,mlat,zRe,zIm,n,m,mask)
      implicit none
      include 'netcdf.inc'
      character*80 modname,hname,uname,gname,fname,rmCom
      character*1 zuv
      integer ncon,n,m,ncid
      real*8 mlon(n,m),mlat(n,m)
      real*8 zRe(ncon,n,m),zIm(ncon,n,m)
      real*8, allocatable:: rtmp1(:,:,:),rtmp2(:,:,:)
      integer mask(n,m),cind(ncon)
      integer nc,ic,k,i,j,funit,status
      character*5 ctmp
c
      uname=' '
      gname=' '
      funit=111
      open(unit=funit,file=modname,status='old')
      read(funit,'(a80)')hname
      hname=rmCom(hname)
      read(funit,'(a80)',end=1)uname
      uname=rmCom(uname)
      read(funit,'(a80)',end=1)gname
      gname=rmCom(gname)
1     close(funit)
      fname=hname
      if(zuv.ne.'z')then
       fname=uname
      endif
      allocate(rtmp1(nc,n,m),rtmp2(nc,n,m))
c
      status = nf_open(trim(fname), nf_nowrite, ncid)
      if(zuv.eq.'z')then
       call read2d_double_nc(ncid,'lon_z',n,m,mlon)
       call read2d_double_nc(ncid,'lat_z',n,m,mlat)
       call read3d_double_nc(ncid,'hRe  ',nc,n,m,rtmp1)
       call read3d_double_nc(ncid,'hIm  ',nc,n,m,rtmp2)
      elseif(zuv.eq.'u'.or.zuv.eq.'U')then
       call read2d_double_nc(ncid,'lon_u',n,m,mlon)
       call read2d_double_nc(ncid,'lat_u',n,m,mlat)
       call read3d_double_nc(ncid,'URe  ',nc,n,m,rtmp1)
       call read3d_double_nc(ncid,'UIm  ',nc,n,m,rtmp2)
      else ! 'v' or 'V'
       call read2d_double_nc(ncid,'lon_v',n,m,mlon)
       call read2d_double_nc(ncid,'lat_v',n,m,mlat)
       call read3d_double_nc(ncid,'VRe  ',nc,n,m,rtmp1)
       call read3d_double_nc(ncid,'VIm  ',nc,n,m,rtmp2)
      endif
      mask=0
      where(rtmp1(1,:,:).ne.0.)mask=1
      zRe=0;zIm=0 ! Bug fix, Dec 2012
      do ic=1,ncon
! Dec 2012 inserted -> if(cind(ic).ne.0), Lana
        if(cind(ic).ne.0)zRe(ic,:,:)=rtmp1(cind(ic),:,:) ! Bug fix Dec 2012 
        if(cind(ic).ne.0)zIm(ic,:,:)=rtmp2(cind(ic),:,:) ! Bug fix Dec 2012
      enddo
      status=nf_close(ncid)
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c define constituent indices in a model for given set of constituent
c constituent names are assumed to be all in LOWER case
      subroutine def_con_ind(c_id,ncon,c_id_mod,nc,cind)
      implicit none
      include 'constit.h'
      character*4 c_id(ncmx),c_id_mod(ncmx)
      integer nc,ncon,cind(ncon),ic1,ic2
c
      cind=0     
      do ic1=1,ncon
       do ic2=1,nc
        if(c_id(ic1).eq.c_id_mod(ic2))cind(ic1)=ic2
       enddo
      enddo
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rd_inp(modname,lltname,zuv,c_id,ncon,
     *                  APRI,Geo,outname,interp)
      implicit none
      include 'constit.h'
      character*80 modname,lltname,outname,tmp,rmCom
      character*1 zuv
      logical APRI,Geo,interp 
      character*4 c_id(ncmx),con
      integer k,ic,ncon
c
      interp=.false.
      read(*,'(a80)',err=1)modname
      modname=rmCom(modname)
      if(trim(modname).eq.'')modname='model.list'
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rd_mod_header_nc(modname,zuv,n,m,th_lim,ph_lim,nc,
     *                         c_id,ll_km)
      implicit none
      include 'netcdf.inc'
      include 'constit.h'
      real*4 th_lim(2),ph_lim(2)
      real*8 dt,stx,sty
      real*8, allocatable:: x(:,:),y(:,:)
      character*4 c_id(ncmx)
      character*80 modname
      character*80 hname,uname,gname,fname,rmCom
      logical ll_km
      character*1 zuv
      integer k,nc,n,m,funit
      integer status,ncid,dimid,len,varid,ncid1
      character*3 dname(4)
      data dname/'nx ','ny ','nc ','nct'/
c
      ll_km=.false.
      uname=' '
      gname=' ' 
      funit=111
      open(unit=funit,file=modname,status='old',err=1)
      read(funit,'(a80)',err=1)hname ! ALWAYS there
      hname=rmCom(hname)
      read(funit,'(a80)',end=5)uname
      uname=rmCom(uname)
      read(funit,'(a80)',end=5)gname
5     close(funit)
      fname=hname
      if(zuv.ne.'z')fname=uname
      if(fname.eq.' ')then
       write(*,*)'No name for model file is given in '
       write(*,*) trim(modname)
       stop
      endif
      if(gname.eq.' '.and.zuv.ne.'z')then
       write(*,*)'No name for bathymetry grid file is given in '
       write(*,*) trim(modname)
       stop
      endif
c
      status = nf_open(gname, nf_nowrite, ncid)
      if(status.ne.0)then
c grid is not given in Model_* file. Apply
c default: grid is on uniform lat/lon grid (like *load file)
       ll_km=.false.
      else
c check is model is on grid unform in km cccccccccccc
       status=nf_inq_varid(ncid,'dt',varid)
       status=nf_get_var_double(ncid,varid,dt)
       if(dt.le.0)ll_km=.true.
      endif
      status=nf_close(ncid)
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c read model file header
      status = nf_open(fname, nf_nowrite, ncid)
      if(status.ne.0)then
       write(*,*)'No file ',trim(fname),' found'
       stop
      endif
cccc READ DIMENSIONS
      do k=1,3
       status=nf_inq_dimid(ncid,trim(dname(k)),dimid)
       status=nf_inq_dimlen(ncid,dimid,len)
       if(k.eq.1)n=len
       if(k.eq.2)m=len
       if(k.eq.3)nc=len
      enddo
      !write(*,*)'n,m,nc:',n,m,nc
      allocate(x(n,m),y(n,m))
      if(zuv.eq.'z')then
       if(ll_km)then
        status = nf_open(gname, nf_nowrite, ncid1)
        call read2d_double_nc(ncid1,'x_z  ',n,m,x)
        call read2d_double_nc(ncid1,'y_z  ',n,m,y)
        status=nf_close(ncid1)
       else
        call read2d_double_nc(ncid,'lon_z',n,m,x)
        call read2d_double_nc(ncid,'lat_z',n,m,y)
       endif
       stx=x(2,1)-x(1,1)
       sty=y(1,2)-y(1,1)
       th_lim(1)=y(1,1)-sty/2
       th_lim(2)=y(1,m)+sty/2
       ph_lim(1)=x(1,1)-stx/2
       ph_lim(2)=x(n,1)+stx/2
      else
       if(ll_km)then
        status = nf_open(gname, nf_nowrite, ncid1)
        call read2d_double_nc(ncid1,'x_u  ',n,m,x)
        call read2d_double_nc(ncid1,'y_u  ',n,m,y)
        status=nf_close(ncid1)
       else
        call read2d_double_nc(ncid,'lon_u',n,m,x)
        call read2d_double_nc(ncid,'lat_u',n,m,y)
       endif
       stx=x(2,1)-x(1,1)
       sty=y(1,2)-y(1,1)
       th_lim(1)=y(1,1)-sty/2
       th_lim(2)=y(1,m)+sty/2
       ph_lim(1)=x(1,1)
       ph_lim(2)=x(n,1)+stx
      endif
      deallocate(x,y)
      status=nf_inq_varid(ncid,'con',varid)
      status=nf_get_var_text(ncid,varid,c_id)
      return
1     write(*,*)'File ',trim(modname),' does NOT exist or empty:'
      write(*,*)'Check line 1 of setup.inp'
      stop
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      character*80 function rmCom(str) ! removes everything after '!' 
                                       ! in a string 
      character*80 str
      integer k 
      k=index(str,'!')
      rmCom=str
      if(k>0)rmCom=str(1:k-1)
      return
      end
C______________________________________________________________________
c
        subroutine interp_r8(uv,nt,n,m,mz,lon,lat,
     *                          xlat,xlon,uv1,ierr)
c
ccc     Interpolates REAL*8 array uv(nt,n,m) at point xlat,xlon
ccc     given latitude lat(n,m) and longitude lon(n,m) of grid
ccc
        implicit none
        integer ierr,n,m,mz(n,m),k,nt
        real*8 uv(nt,n,m),lon(n,m),lat(n,m)
        real*8 uv1(nt),xlonc
        real*8 ww(0:1,0:1),dlat,dlon
        real*8 xlon,xlat,minlon,maxlon,minlat,maxlat
        integer iw(0:1),jw(0:1)
c
        ierr = 0
c check on xlon longitude convention
        xlonc=xlon
        minlon=minval(lon)
        maxlon=maxval(lon)
        minlat=minval(lat)
        maxlat=maxval(lat)
        if(xlonc.lt.minlon)xlonc=xlonc+360.
        if(xlonc.gt.maxlon)xlonc=xlonc-360.
        if(xlonc.lt.minlon.or.xlonc.gt.maxlon)then
          ierr=-1
          return
        endif
        if(xlat.lt.minlat.or.xlat.gt.maxlat)then
          ierr=-1
          return
        endif
c
        call BSI_weights_r8(xlat,xlonc,lat,lon,
     *                      mz,n,m,ww,iw,jw)
        if(sum(ww).lt.0.01) then
          ierr = -2
          do k = 1,nt
            uv1(k) = 0.
          enddo
        else
         ierr = 0
         do k = 1,nt
            uv1(k) = uv(k,iw(0),jw(0))*ww(0,0)+
     &               uv(k,iw(1),jw(0))*ww(1,0)+
     &               uv(k,iw(0),jw(1))*ww(0,1)+
     &               uv(k,iw(1),jw(1))*ww(1,1)
         enddo   ! k
        endif
c
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c bilinear spline interpolation weights
         subroutine BSI_weights_r8(theta,phi,lat,lon,
     *                             mask,n,m,ww,iw,jw)
         implicit none
         integer n,m,mask(n,m)          ! grid dimensions and mask
         integer i,j
         real*8 lat(n,m),lon(n,m),theta,phi
         real*8 xi,xj,x,y
         integer i0,j0,i1,j1 
         real*8 ww(0:1,0:1)
         integer iw(0:1),jw(0:1)
         integer ipshft,sm
         real*8 w00,w01,w10,w11,wtot
c
         do i=1,n-1
          do j=1,m-1
           if(lat(i,j).le.theta.and.theta.le.lat(i,j+1).and.
     *        lon(i,j).le.phi.and.phi.le.lon(i+1,j))then
              i0=i
              j0=j
              x=(phi-lon(i,j))/(lon(i+1,j)-lon(i,j))
              y=(theta-lat(i,j))/(lat(i,j+1)-lat(i,j))
              go to 1
           endif
          enddo
         enddo
         do j=1,m-1
           if(lat(n,j).le.theta.and.theta.le.lat(n,j+1).and.
     *        lon(n,j).le.phi.and.phi.le.lon(1,j)+360.)then
              i0=n
              j0=j
              x=(phi-lon(n,j))/(lon(1,j)+360-lon(n,j))
              y=(theta-lat(n,j))/(lat(n,j+1)-lat(n,j))
              go to 1
           endif
          enddo
          ww=0.
          return
c
c        compute weights for bilinear spline interpolation; only
c        use ocean nodes (for which mask mask is = 1)
1        j1 = ipshft(j0,1,m)
         i1 = ipshft(i0,1,n)
         ww=0.
         sm=mask(i0,j0)+mask(i0,j1)+mask(i1,j0)+mask(i1,j1)
         if(sm.gt.0)then
          w00 = (1.-x)*(1.-y)*mask(i0,j0)
          w01 = (1.-x)*y*mask(i0,j1)
          w10 = x*(1.-y)*mask(i1,j0)
          w11 = x*y*mask(i1,j1)
          wtot = w00+w01+w10+w11
          if(wtot.eq.0)then
            ww=0
            return 
          endif
c
          ww(0,0) = w00/wtot
          ww(0,1) = w01/wtot
          ww(1,0) = w10/wtot
          ww(1,1) = w11/wtot
         endif
c
         iw(0)=i0
         iw(1)=i1
         jw(0)=j0
         jw(1)=j1
c 
         return
         end
c____________________________________________________________________
c
      function ipshft(i,ish,n)
c      periodic shift maps i to i+ish  , mod n;
c        (always between 1 and n;  never 0 )
      integer i,ipshft,n,ish
      ipshft = mod(i+ish+n-1,n)+1
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c make_a:: computes A matrix elements for one data point 
	subroutine make_a(interp,ind,nc,time,pu,pf,w,A,l_sal) 
        implicit none
	include 'constit.h'
        integer, parameter:: ncon = 17 ! for case of using weights.h
c
        logical interp,l_sal
        real omega(ncmx),phase(ncmx)
	integer ind(nc)
	real w(ncon,8)
	real*8 pu(ncmx),pf(ncmx)
        real*8 time
        complex A(nc),c(ncmx)
        integer i,j,nc
c If l_sal=.true. - NO solid Earth correction is applied
c USING beta_SE coefficients
c
        if(.not.interp)then
          do j=1,ncmx
           omega(j)=omega_d(j)
           phase(j)=phase_mkB(j)
          enddo
  	  do j=1,nc
            i=ind(j)
            if(i.ne.0)then
 	    c(j) = cmplx( pf(i)*cos(omega(i)*time+phase(i)+pu(i)),
     *	                  pf(i)*sin(omega(i)*time+phase(i)+pu(i)))
            endif
	  enddo
c  remove solid Earth tide
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
c this is the case when w from weights.h is used -> see
c comments to ../include/weights.h
c----------------------------------------------------------c        
          omega(1:ncon)=omega_d(1:ncon)
          phase(1:ncon)=phase_mkB(1:ncon)
c
  	  do i=1,ncon 
 	    c(i) = cmplx( pf(i)*cos(omega(i)*time+phase(i)+pu(i)),
     *	                  pf(i)*sin(omega(i)*time+phase(i)+pu(i)))
	  enddo
          A=cmplx(0,0)
c
c ind(j)=0 means the constituent is excluded
          do i=1,ncon
            do j=1,nc 
             if(ind(j).ne.0)A(j)=A(j)+c(i)*beta_SE(i)*w(i,ind(j))
            enddo
          enddo
c---------------------------------------------------------c

        endif      
c
        return
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
        subroutine mkw(interp,ind,nc,wr)
        include 'weights.h'
        real wr(17,8)
        logical interp
        integer i,j,nc,ind(nc)
        wr=w
        if(interp)return
c
        do j=1,nc
         if(ind(j).ne.0)wr(ind(j),:)=0.
        enddo
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	real function height(A,P,NC)
c returns height from model array of complex constituents
        complex A(nc),P(nc)
        real sum
        integer i,nc
        if(nc.eq.0)then
          height=0.
          return
        endif
        sum=0.
c height(i)=sum_of_real(A(i)*P(i))
        do i=1,nc
         sum=sum+real(P(i))*real(A(i))-aimag(P(i))*aimag(A(i))
        enddo
        height=sum
        return
        end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c argUMENTS and ASTROL subroutines SUPPLIED by RICHARD RAY, March 1999
c attached to OTIS by Lana Erofeeva (subroutine nodal.f)
c NOTE - "no1" in constit.h corresponds to "M1" in arguments
	subroutine nodal(dtime,latitude,pu,pf)
        implicit none
        include 'constit.h'
c
        integer index(ncmx),i
        real*8 latitude,pu(ncmx),pf(ncmx)
        real*8 dtime,arg(53),f(53),u(53),pi
        data pi/3.14159265358979/
c index gives correspondence between constit.h and Richard's subroutines
c constit.h:       M2,S2,K1,O1,N2,P1,K2,q1,2N2,mu2,nu2,L2,t2,
c                  J1,M1(no1),OO1,rho1,Mf,Mm,SSA,M4,
c                  MS4,MN4,M6,M8,MK3,S6,2SM2,2MK3
c
	data index/30,35,19,12,27,17,37,10,25,26,28,33,34,
     *             23,14,24,11,5,3,2,45,46,44,50,0,42,51,40,0/
c f, u - same as pf, pu in old nodal.f; arg is NOT needed; 
c dtime - MJD
        call arguments(dtime,arg,f,u)
        pu=0
        pf=1
        do i=1,ncmx
c u is returned by "arguments" in degrees
         if(index(i).gt.0)then 
          pu(i)=u(index(i))*pi/180.
          pf(i)=f(index(i))
         endif
         !write(*,*)pu(i),pf(i)
        enddo
        return
        end 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine arguments( time1, arg, f, u )
      implicit none
*
*   Kernel routine for subroutine hat53.    Calculate tidal arguments.
*
      real*8 time1, arg(*), f(*), u(*)
      real*8 shpn(4),s,h,p,omega,pp,hour,t1,t2
      real*8 tmp1,tmp2,temp1,temp2
      real*8 cosn,cos2n,sinn,sin2n,sin3n
      real*8 zero,one,two,three,four,five
      real*8 fiften,thirty,ninety
      real*8 pi, rad
      parameter       (pi=3.141592654d0, rad=pi/180.D0)
      parameter   (zero=0.d0, one=1.d0)
      parameter   (two=2.d0, three=3.d0, four=4.d0, five=5.d0)
      parameter   (fiften=15.d0, thirty=30.d0, ninety=90.d0)
      parameter   (pp=282.94) ! solar perigee at epoch 2000.
      equivalence (shpn(1),s),(shpn(2),h),(shpn(3),p),(shpn(4),omega)
*
*     Determine equilibrium arguments
*     -------------------------------
      call astrol( time1, shpn )
      hour = (time1 - int(time1))*24.d0
      t1 = fiften*hour
      t2 = thirty*hour
      arg( 1) = h - pp                                  ! Sa
      arg( 2) = two*h                                   ! Ssa
      arg( 3) = s - p                                   ! Mm
      arg( 4) = two*s - two*h                           ! MSf
      arg( 5) = two*s                                   ! Mf
      arg( 6) = three*s - p                             ! Mt
      arg( 7) = t1 - five*s + three*h + p - ninety      ! alpha1
      arg( 8) = t1 - four*s + h + two*p - ninety        ! 2Q1
      arg( 9) = t1 - four*s + three*h - ninety          ! sigma1
      arg(10) = t1 - three*s + h + p - ninety           ! q1
      arg(11) = t1 - three*s + three*h - p - ninety     ! rho1
      arg(12) = t1 - two*s + h - ninety                 ! o1
      arg(13) = t1 - two*s + three*h + ninety           ! tau1
      arg(14) = t1 - s + h + ninety                     ! M1
      arg(15) = t1 - s + three*h - p + ninety           ! chi1
      arg(16) = t1 - two*h + pp - ninety                ! pi1
      arg(17) = t1 - h - ninety                         ! p1
      arg(18) = t1 + ninety                             ! s1
      arg(19) = t1 + h + ninety                         ! k1
      arg(20) = t1 + two*h - pp + ninety                ! psi1
      arg(21) = t1 + three*h + ninety                   ! phi1
      arg(22) = t1 + s - h + p + ninety                 ! theta1
      arg(23) = t1 + s + h - p + ninety                 ! J1
      arg(24) = t1 + two*s + h + ninety                 ! OO1
      arg(25) = t2 - four*s + two*h + two*p             ! 2N2
      arg(26) = t2 - four*s + four*h                    ! mu2
      arg(27) = t2 - three*s + two*h + p                ! n2
      arg(28) = t2 - three*s + four*h - p               ! nu2
      arg(29) = t2 - two*s + h + pp                     ! M2a
      arg(30) = t2 - two*s + two*h                      ! M2
      arg(31) = t2 - two*s + three*h - pp               ! M2b
      arg(32) = t2 - s + p + 180.d0                     ! lambda2
      arg(33) = t2 - s + two*h - p + 180.d0             ! L2
      arg(34) = t2 - h + pp                             ! t2
      arg(35) = t2                                      ! S2
      arg(36) = t2 + h - pp + 180.d0                    ! R2
      arg(37) = t2 + two*h                              ! K2
      arg(38) = t2 + s + two*h - pp                     ! eta2
      arg(39) = t2 - five*s + 4.0*h + p                 ! MNS2
      arg(40) = t2 + two*s - two*h                      ! 2SM2
      arg(41) = 1.5*arg(30)                             ! M3
      arg(42) = arg(19) + arg(30)                       ! MK3
      arg(43) = three*t1                                ! S3
      arg(44) = arg(27) + arg(30)                       ! MN4
      arg(45) = two*arg(30)                             ! M4
      arg(46) = arg(30) + arg(35)                       ! MS4
      arg(47) = arg(30) + arg(37)                       ! MK4
      arg(48) = four*t1                                 ! S4
      arg(49) = five*t1                                 ! S5
      arg(50) = three*arg(30)                           ! M6
      arg(51) = three*t2                                ! S6
      arg(52) = 7.0*t1                                  ! S7
      arg(53) = four*t2                                 ! S8
*
*     determine nodal corrections f and u 
*     -----------------------------------
      sinn = sin(omega*rad)
      cosn = cos(omega*rad)
      sin2n = sin(two*omega*rad)
      cos2n = cos(two*omega*rad)
      sin3n = sin(three*omega*rad)
      f( 1) = one                                     ! Sa
      f( 2) = one                                     ! Ssa
      f( 3) = one - 0.130*cosn                        ! Mm
      f( 4) = one                                     ! MSf
      f( 5) = 1.043 + 0.414*cosn                      ! Mf
      f( 6) = sqrt((one+.203*cosn+.040*cos2n)**2 + 
     *              (.203*sinn+.040*sin2n)**2)        ! Mt

      f( 7) = one                                     ! alpha1
      f( 8) = sqrt((1.+.188*cosn)**2+(.188*sinn)**2)  ! 2Q1
      f( 9) = f(8)                                    ! sigma1
      f(10) = f(8)                                    ! q1
      f(11) = f(8)                                    ! rho1
      f(12) = sqrt((1.0+0.189*cosn-0.0058*cos2n)**2 +
     *             (0.189*sinn-0.0058*sin2n)**2)      ! O1
      f(13) = one                                     ! tau1
ccc   tmp1  = 2.*cos(p*rad)+.4*cos((p-omega)*rad)
ccc   tmp2  = sin(p*rad)+.2*sin((p-omega)*rad)         ! Doodson's
      tmp1  = 1.36*cos(p*rad)+.267*cos((p-omega)*rad)  ! Ray's
      tmp2  = 0.64*sin(p*rad)+.135*sin((p-omega)*rad)
      f(14) = sqrt(tmp1**2 + tmp2**2)                 ! M1
      f(15) = sqrt((1.+.221*cosn)**2+(.221*sinn)**2)  ! chi1
      f(16) = one                                     ! pi1
      f(17) = one                                     ! P1
      f(18) = one                                     ! S1
      f(19) = sqrt((1.+.1158*cosn-.0029*cos2n)**2 + 
     *             (.1554*sinn-.0029*sin2n)**2)       ! K1
      f(20) = one                                     ! psi1
      f(21) = one                                     ! phi1
      f(22) = one                                     ! theta1
      f(23) = sqrt((1.+.169*cosn)**2+(.227*sinn)**2)  ! J1
      f(24) = sqrt((1.0+0.640*cosn+0.134*cos2n)**2 +
     *             (0.640*sinn+0.134*sin2n)**2 )      ! OO1
      f(25) = sqrt((1.-.03731*cosn+.00052*cos2n)**2 +
     *             (.03731*sinn-.00052*sin2n)**2)     ! 2N2
      f(26) = f(25)                                   ! mu2
      f(27) = f(25)                                   ! N2
      f(28) = f(25)                                   ! nu2
      f(29) = one                                     ! M2a
      f(30) = f(25)                                   ! M2
      f(31) = one                                     ! M2b
      f(32) = one                                     ! lambda2
      temp1 = 1.-0.25*cos(two*p*rad)
     *        -0.11*cos((two*p-omega)*rad)-0.04*cosn
      temp2 = 0.25*sin(two*p)+0.11*sin((two*p-omega)*rad)
     *        + 0.04*sinn
      f(33) = sqrt(temp1**2 + temp2**2)               ! L2
      f(34) = one                                     ! t2
      f(35) = one                                     ! S2
      f(36) = one                                     ! R2
      f(37) = sqrt((1.+.2852*cosn+.0324*cos2n)**2 +
     *             (.3108*sinn+.0324*sin2n)**2)       ! K2
      f(38) = sqrt((1.+.436*cosn)**2+(.436*sinn)**2)  ! eta2
      f(39) = f(30)**2                                ! MNS2
      f(40) = f(30)                                   ! 2SM2
      f(41) = one   ! wrong                           ! M3
      f(42) = f(19)*f(30)                             ! MK3
      f(43) = one                                     ! S3
      f(44) = f(30)**2                                ! MN4
      f(45) = f(44)                                   ! M4
      f(46) = f(44)                                   ! MS4
      f(47) = f(30)*f(37)                             ! MK4
      f(48) = one                                     ! S4
      f(49) = one                                     ! S5
      f(50) = f(30)**3                                ! M6
      f(51) = one                                     ! S6
      f(52) = one                                     ! S7
      f(53) = one                                     ! S8

         u( 1) = zero                                    ! Sa
         u( 2) = zero                                    ! Ssa
         u( 3) = zero                                    ! Mm
         u( 4) = zero                                    ! MSf
         u( 5) = -23.7*sinn + 2.7*sin2n - 0.4*sin3n      ! Mf
         u( 6) = atan(-(.203*sinn+.040*sin2n)/
     *                 (one+.203*cosn+.040*cos2n))/rad   ! Mt
         u( 7) = zero                                    ! alpha1
         u( 8) = atan(.189*sinn/(1.+.189*cosn))/rad      ! 2Q1
         u( 9) = u(8)                                    ! sigma1
         u(10) = u(8)                                    ! q1
         u(11) = u(8)                                    ! rho1
         u(12) = 10.8*sinn - 1.3*sin2n + 0.2*sin3n       ! O1
         u(13) = zero                                    ! tau1
         u(14) = atan2(tmp2,tmp1)/rad                    ! M1
         u(15) = atan(-.221*sinn/(1.+.221*cosn))/rad     ! chi1
         u(16) = zero                                    ! pi1
         u(17) = zero                                    ! P1
         u(18) = zero                                    ! S1
         u(19) = atan((-.1554*sinn+.0029*sin2n)/
     *                (1.+.1158*cosn-.0029*cos2n))/rad   ! K1
         u(20) = zero                                    ! psi1
         u(21) = zero                                    ! phi1
         u(22) = zero                                    ! theta1
         u(23) = atan(-.227*sinn/(1.+.169*cosn))/rad     ! J1
         u(24) = atan(-(.640*sinn+.134*sin2n)/
     *                (1.+.640*cosn+.134*cos2n))/rad     ! OO1
         u(25) = atan((-.03731*sinn+.00052*sin2n)/ 
     *                (1.-.03731*cosn+.00052*cos2n))/rad ! 2N2
         u(26) = u(25)                                   ! mu2
         u(27) = u(25)                                   ! N2
         u(28) = u(25)                                   ! nu2
         u(29) = zero                                    ! M2a
         u(30) = u(25)                                   ! M2
         u(31) = zero                                    ! M2b
         u(32) = zero                                    ! lambda2
         u(33) = atan(-temp2/temp1)/rad                  ! L2
         u(34) = zero                                    ! t2
         u(35) = zero                                    ! S2
         u(36) = zero                                    ! R2
         u(37) = atan(-(.3108*sinn+.0324*sin2n)/ 
     *                (1.+.2852*cosn+.0324*cos2n))/rad   ! K2
         u(38) = atan(-.436*sinn/(1.+.436*cosn))/rad     ! eta2
         u(39) = u(30)*two                               ! MNS2
         u(40) = u(30)                                   ! 2SM2
         u(41) = 1.5d0*u(30)                             ! M3
         u(42) = u(30) + u(19)                           ! MK3
         u(43) = zero                                    ! S3
         u(44) = u(30)*two                               ! MN4
         u(45) = u(44)                                   ! M4
         u(46) = u(30)                                   ! MS4
         u(47) = u(30)+u(37)                             ! MK4
         u(48) = zero                                    ! S4
         u(49) = zero                                    ! S5
         u(50) = u(30)*three                             ! M6
         u(51) = zero                                    ! S6
         u(52) = zero                                    ! S7
         u(53) = zero                                    ! S8

      return
      end


*======================================================================
      SUBROUTINE ASTROL( time, SHPN )
*
*  Computes the basic astronomical mean longitudes  s, h, p, N.
*  Note N is not N', i.e. N is decreasing with time.
*  These formulae are for the period 1990 - 2010, and were derived
*  by David Cartwright (personal comm., Nov. 1990).
*  time is UTC in decimal MJD.
*  All longitudes returned in degrees.
*  R. D. Ray    Dec. 1990
*
*  Non-vectorized version.
*
c      IMPLICIT real*8 (A-H,O-Z)
      real*8 circle,shpn,t,time
      DIMENSION  SHPN(4)
      PARAMETER  (CIRCLE=360.0D0)
*
      T = time - 51544.4993D0
*
*     mean longitude of moon
*     ----------------------
      SHPN(1) = 218.3164D0 + 13.17639648D0 * T
*
*     mean longitude of sun
*     ---------------------
      SHPN(2) = 280.4661D0 +  0.98564736D0 * T
*
*     mean longitude of lunar perigee
*     -------------------------------
      SHPN(3) =  83.3535D0 +  0.11140353D0 * T
*
*     mean longitude of ascending lunar node
*     --------------------------------------
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
c-------------------------------------
      SUBROUTINE CALDAT (JULIAN,MM,ID,IYYY)                             
C                                                                       
C   ROUTINE CONVERTS JULIAN DAY TO MONTH, DAY, & YEAR.                      
C   THIS CODE IS LIFTED FROM THE BOOK:                                      
C   W.H. PRESS ET AL., NUMERICAL RECIPES, CAMBRIDGE UNIV. PRESS, 1986.      
C   THE ONLY MODIFICATION IS THAT REAL ARITHMETIC IS DONE IN R*8.           
C                                                                           
C   TO CONVERT MODIFIED JULIAN DAY, CALL THIS ROUTINE WITH                 
C     JULIAN = MJD + 2400001                                               
C                                                                          
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc Lana Erofeeva, Feb 2002 cccccccccccccccccccccccccccccccccccccc
ccc Debugged for consistency with CALDAT for any dates mjd>=0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
      subroutine date_mjd(mm,id,iyyy,mjd)
      implicit none
c converts date to mjd
c INPUT:  id - day, mm - month, iyyy - year
c OUTPUT: mjd>0 - modified julian days
c date>=11.17.1858 corresponds to mjd=0
      integer dpm(12),days,i,nleap,k
      integer mm,id,iyyy,mjd
      data dpm/31,28,31,30,31,30,31,31,30,31,30,31/
      mjd=0
c NO earlier dates
      if(iyyy.lt.1858)iyyy=1858
      if(iyyy.eq.1858.and.mm.gt.11)mm=11
      if(iyyy.eq.1858.and.mm.eq.11.and.id.gt.17)id=17
c
      days=0 
      do i=1,mm-1
       days=days+dpm(i)
       if(i.eq.2.and.int(iyyy/4)*4.eq.iyyy)days=days+1         
      enddo
      days=days+id-321
c leap day correction
      do k=1900,iyyy,100
       if(iyyy.eq.k.and.mm.gt.2)days=days-1
      enddo
      do k=2000,iyyy,400
       if(iyyy.eq.k.and.mm.gt.2)days=days+1
      enddo
c EACH 4th year is leap year
      nleap=int((iyyy-1-1860)*0.25)            
      if(iyyy.gt.1860)nleap=nleap+1
c EXCEPT
      do k=1900,iyyy-1,100
       if(iyyy.gt.k)nleap=nleap-1
       if(iyyy.eq.k.and.mm.gt.2)days=days-1
      enddo
c BUT each in the row 2000:400:... IS LEAP year again
      do k=2000,iyyy-1,400
       if(iyyy.gt.k)nleap=nleap+1
       if(iyyy.eq.k.and.mm.gt.2)days=days+1
      enddo
      mjd=365*(iyyy-1858)+nleap
     *    +days   
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine ptide(z1,cid,ncon,ind,lat,time_mjd,ntime,
     *                   interp,zpred)
        implicit none
        include 'constit.h'
        include 'weights.h'
        integer ncon,ntime,k,ind(ncon),ierr
        complex z1(ncon)
        real lat,zpred(ntime)
        real*8 time_mjd(ntime),time,dh
        character*4 cid(ncon)
        integer SecondsPerDay
        real height
        real ww(17,8)
c        complex z1(ncon)
        complex, allocatable:: A(:)
c nodal arguments shoud be real*8
        real*8 pu(ncmx),pf(ncmx),dlat,t_mjd
        logical interp
c
        if(interp)call mkw(interp,ind,ncon,ww)
        allocate(A(ncon))
	SecondsPerDay=86400
        dlat=lat
        ierr=0
        dh=0. 
        do k=1,ntime
	 call nodal(time_mjd(k),dlat,pu,pf)
c to use phase shifts from constit.h time should be in seconds
c relatively Jan 1 1992 (48622mjd)       
         time=(time_mjd(k)-dble(48622))*dble(SecondsPerDay)
c .true. means NO solid Earth correction will be applied in make_a
         call make_a(.false.,ind,ncon,time,pu,pf,ww,A,.true.)
         zpred(k)=height(A,z1,ncon)
         if(interp)call infer_minor(z1,cid,ncon,time_mjd(k),dh,ierr)
         if(ierr.eq.-1)then
          write(*,*)'Not enough constituents for inference of'
          write(*,*)'minor constituents: IGNORED'
          interp=.false.
         endif
         zpred(k)=zpred(k)+dh ! add minor constituents
        enddo
        deallocate(A)
        return
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         subroutine def_cid(nc0,cid,ind)
         implicit none
         include 'constit.h'
         integer nc0,ic,jc,ind(nc0),k
         character*4 cid(ncmx)
c
         k=1
         do ic = 1,nc0
          do jc = 1,ncmx             
            if(cid(ic).eq.constid(jc)) go to 4
          enddo
          write(*,*) 'WARNING: Constituent ID ',cid(ic),
     *               ' is not allowed'
          ind(k)=0
          k=k+1
4         continue
          ind(k)=jc
          k=k+1
         enddo
         return
         end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Lana Erofeeva, June 2004
c Based on Richard Ray code perth2
c Return correction for the 16 minor constituents
c listed in subroutine
c         infer minor tides
c         -----------------
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
         data cid8/'q1  ','o1  ','p1  ','k1  ',
     *             'n2  ','m2  ','s2  ','k2  '/
c re-order to correspond to cid8
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
c            
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
c
          hour = (time - int(time))*24.D0
          t1 = 15.*hour
          t2 = 30.*hour
          call astrol( time, SHPN )
c
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
c
c     determine nodal corrections f and u
         sinn = SIN(omega*rad)
         cosn = COS(omega*rad)
         sin2n = SIN(2.*omega*rad)
         cos2n = COS(2.*omega*rad)
         F = 1.
         f(1) = sqrt((1.0 + 0.189*cosn - 0.0058*cos2n)**2 +
     *                (0.189*sinn - 0.0058*sin2n)**2)
         f(2) = f(1)
         f(3) = f(1)
         f(4) = sqrt((1.0 + 0.185*cosn)**2 + (0.185*sinn)**2)
         f(5) = sqrt((1.0 + 0.201*cosn)**2 + (0.201*sinn)**2)
         f(6) = sqrt((1.0 + 0.221*cosn)**2 + (0.221*sinn)**2)
         f(10) = sqrt((1.0 + 0.198*cosn)**2 + (0.198*sinn)**2)
         f(11) = sqrt((1.0 + 0.640*cosn + 0.134*cos2n)**2 +
     *                (0.640*sinn + 0.134*sin2n)**2 )
         f(12) = sqrt((1.0 - 0.0373*cosn)**2 + (0.0373*sinn)**2)
         f(13) = f(12)
         f(14) = f(12)
         f(16) = f(12)
         f(17) = sqrt((1.0 + 0.441*cosn)**2 + (0.441*sinn)**2)
c
         U = 0.
         u(1) = atan2(0.189*sinn - 0.0058*sin2n,
     *                1.0 + 0.189*cosn - 0.0058*sin2n)/rad
         u(2) = u(1)
         u(3) = u(1)
         u(4) = atan2( 0.185*sinn, 1.0 + 0.185*cosn)/rad
         u(5) = atan2(-0.201*sinn, 1.0 + 0.201*cosn)/rad
         u(6) = atan2(-0.221*sinn, 1.0 + 0.221*cosn)/rad
         u(10) = atan2(-0.198*sinn, 1.0 + 0.198*cosn)/rad
         u(11) = atan2(-0.640*sinn - 0.134*sin2n,
     *                 1.0 + 0.640*cosn + 0.134*cos2n)/rad
         u(12) = atan2(-0.0373*sinn, 1.0 - 0.0373*cosn)/rad
         u(13) = u(12)
         u(14) = u(12)
         u(16) = u(12)
         u(17) = atan2(-0.441*sinn, 1.0 + 0.441*cosn)/rad

c     sum over all tides
c     ------------------
        dh = 0.D0
        do i=1,18
         dh = dh + dreal(zmin(i))*f(i)*cos((arg(i)+u(i))*rad)-
     *             dimag(zmin(i))*f(i)*sin((arg(i)+u(i))*rad)
        enddo
        return
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rd_inp_JS(setup_file_unit,modname,zuv,c_id,ncon,
     *                  APRI,Geo,interp)
c
      implicit none
      include 'constit.h'
      integer setup_file_unit
      character*80 modname,outname,tmp,rmCom
      character*80 lltname
      character*1 zuv
      logical APRI,Geo,interp 
      character*4 c_id(ncmx),con
      integer k,ic,ncon
c
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine read2d_int_nc(ncid,varname,nx,ny,iar)
      implicit none
      include 'netcdf.inc'
      character*5 varname
      integer*4 ncid,status,varid,nx,ny,iar(nx,ny),k
      integer*4, allocatable:: itmp(:,:)
c
      allocate(itmp(ny,nx))
      status=nf_inq_varid(ncid,trim(varname),varid)
      if(status.ne.0)go to 1
      status=nf_get_var_int(ncid,varid,itmp)
      if(status.ne.0)go to 2
      do k=1,nx
       iar(k,:)=itmp(:,k)
      enddo
      deallocate(itmp)
      return
1     write(*,*)'No variable ',trim(varname),' found'
      stop
2     write(*,*)'Error while reading integer 2d array ',
     *          trim(varname)
      stop
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine read2d_double_nc(ncid,varname,nx,ny,dar)
      implicit none
      include 'netcdf.inc'
      character*5 varname
      integer*4 ncid,status,varid,nx,ny,k
      real*8 dar(nx,ny)
      real*8, allocatable:: dtmp(:,:)
c
      allocate(dtmp(ny,nx))
      status=nf_inq_varid(ncid,trim(varname),varid)
      if(status.ne.0)go to 1
      status=nf_get_var_double(ncid,varid,dtmp)
      if(status.ne.0)go to 2
      do k=1,nx
       dar(k,:)=dtmp(:,k)
      enddo
      deallocate(dtmp)
      return
1     write(*,*)'No variable ',trim(varname),' found'
      stop
2     write(*,*)'Error while reading double 2d array ',
     *          trim(varname)
      stop
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine read3d_double_nc(ncid,varname,nc,nx,ny,dar)
      implicit none
      include 'netcdf.inc'
      character*5 varname
      integer*4 ncid,status,varid,nx,ny,nc,k,j
      real*8 dar(nc,nx,ny)
      real*8, allocatable:: dtmp(:,:,:)
c
      allocate(dtmp(ny,nx,nc))
      status=nf_inq_varid(ncid,trim(varname),varid)
      if(status.ne.0)go to 1
      status=nf_get_var_double(ncid,varid,dtmp)
      if(status.ne.0)go to 2
      do j=1,nc
       do k=1,nx
        dar(j,k,:)=dtmp(:,k,j)
       enddo
      enddo
      deallocate(dtmp)
      return
1     write(*,*)'No variable ',trim(varname),' found'
      stop
2     write(*,*)'Error while reading double 3d array ',
     *          trim(varname)
      stop
      end

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
      program extract_HC
ccc   LANA: 2008 remake of extract_HC(2004) for netcdf 
ccc   reads OTIS netcdf model file
ccc   (elevations OR transports), reads a list of locations 
ccc   and outputs ASCII file with the complex amplitudes/amp,phases of
ccc   elevations/transports/currents interpolated at the locations
ccc  
      implicit none
      include 'netcdf.inc'
      include 'constit.h'
      integer ncid,dimid,status
      real*8, allocatable:: zlRe(:,:,:),zlIm(:,:,:)
      real*8, allocatable:: mlon(:,:),mlat(:,:),slon(:,:),slat(:,:)
      real*8, allocatable:: glon(:,:),glat(:,:)
      real*8, allocatable:: zRe(:,:,:),zIm(:,:,:),depth(:,:)
      complex*16, allocatable:: z1(:),zl1(:)
      real*8, allocatable:: lat(:),lon(:),lon0(:),dtmp(:)
      real*8 latp,lonp,d1
      real th_lim(2),ph_lim(2),dum,lth_lim(2),lph_lim(2)
      integer, allocatable:: mask(:,:),cind(:),lcind(:),
     *                       mz(:,:),mzl(:,:)
c
      character*4 c_id(ncmx),c_id_mod(ncmx),lc_id(ncmx)
      character*80 modname,lltname,outname,gname,ctmp,lname
      character*2000 fmt,ftmp
      character*80 rmCom
      character*1 zuv,c1
      character*80 xy_ll_sub
      logical APRI,geo,interp_micon,ll_km,ltmp
      integer ncon,nc,n,m,ndat,k,ierr,ierr1,ic,n0,m0
      integer ncl,nl,ml,nmod,imod,ibl
c
      ll_km=.false.
      nmod=1
      ibl=0
      lname='DATA/Model_load'
      call rd_inp(modname,lltname,zuv,c_id,ncon,APRI,geo,
     *            outname,interp_micon)
      if(trim(modname).eq.'model.list')then
        nmod=0
        open(unit=12,file='model.list',status='old',err=11)
        write(*,*)'Models to extract HC from:'
13      read(12,'(a80)',end=12)ctmp
        write(*,*)trim(ctmp)
        nmod=nmod+1
        go to 13
12      rewind(12)
        read(12,'(a80)')modname
      endif
      write(*,*)
      write(*,*)'Lat/Lon file:',trim(lltname)
      if(ncon.gt.0)write(*,*)'Constituents to include: ',c_id(1:ncon)
      if(geo.and.zuv.eq.'z')then
         write(*,*)'Extract GEOCENTRIC tide HC'
      else
         write(*,*)'Extract OCEAN tide HC'
      endif
c
      open(unit=11,file=outname,status='unknown')
      do imod=1,nmod
      write(*,*)
      if(imod.gt.1)read(12,'(a80)')modname

      call rd_mod_header_nc(modname,zuv,n,m,
     *          th_lim,ph_lim,nc,c_id_mod,ll_km)

      write(*,*)'Model:        ',trim(modname(12:80))
      write(11,*)'Model:        ',trim(modname(12:80))
      if(ll_km)then
       write(*,*)'Model is on grid uniform in km'
       write(*,*)'x(km) limits:   ',th_lim
       write(*,*)'y(km) limits:   ',ph_lim
      else
       write(*,*)'Model is on grid uniform in lat,lon'
       write(*,*)'Lat limits:   ',th_lim
       write(*,*)'Lon limits:   ',ph_lim
      endif
      write(*,*)'Constituents: ',c_id_mod(1:nc)
c
      if(zuv.eq.'z')write(*,*)'Output elevations (m)'
      if(zuv.eq.'U')write(*,*)'Output WE transport (m^2/s)'
      if(zuv.eq.'u')write(*,*)'Output WE velocity  (cm/s)'
      if(zuv.eq.'V')write(*,*)'Output SN transport (m^2/s)'
      if(zuv.eq.'v')write(*,*)'Output SN velocity  (cm/s)'
c
      if(zuv.eq.'z')write(11,*)'Elevations (m)'
      if(zuv.eq.'U')write(11,*)'WE transport (m^2/s)'
      if(zuv.eq.'u')write(11,*)'WE velocity  (cm/s)'
      if(zuv.eq.'V')write(11,*)'SN transport (m^2/s)'
      if(zuv.eq.'v')write(11,*)'SN velocity  (cm/s)'
c
      if(ncon.eq.0)then
       ibl=1
       ncon=nc
       c_id=c_id_mod
       write(*,*)'Constituents to include: ',c_id(1:ncon)
      endif
c
      allocate(cind(ncon))
      call def_con_ind(c_id,ncon,c_id_mod,nc,cind)
c 
      ndat=0
      open(unit=1,file=lltname,status='old',err=1)
3     read(1,*,end=2)dum,dum
      ndat=ndat+1
      go to 3
2     rewind(1)
      allocate(lat(ndat),lon(ndat),lon0(ndat))
      do k=1,ndat
        read(1,*)lat(k),lon(k)
        lon0(k)=lon(k)
      enddo ! k
      close(1)
c
      allocate(zRe(ncon,n,m),zIm(ncon,n,m),z1(ncon),mask(n,m),
     *         mlon(n,m),mlat(n,m),dtmp(ncon))
      write(*,'(a,$)')' Reading model...'
      call rd_mod_body_nc(modname,zuv,cind,nc,ncon,
     *                    mlon,mlat,zRe,zIm,n,m,mask)
      write(*,*)'done'
      if(zuv.eq.'z'.and.geo)then
       ltmp=.false.
       call rd_mod_header_nc(lname,'z',nl,ml,
     *                    lth_lim,lph_lim,ncl,lc_id,ltmp)
       !write(*,*)ncl,lc_id(1:ncl)
       allocate(lcind(ncon))
       call def_con_ind(c_id,ncon,lc_id,ncl,lcind)
       !write(*,*)ncon,c_id(1:ncon)
       !write(*,*)lcind
       allocate(zlRe(ncon,nl,ml),zlIm(ncon,nl,ml),zl1(ncon),
     *          mzl(nl,ml),slon(nl,ml),slat(nl,ml))
       write(*,*)
       write(*,'(a,$)')' Reading file to apply load corrections...'
       call rd_mod_body_nc(lname,'z',lcind,ncl,ncon,
     *                     slon,slat,zlRe,zlIm,nl,ml,mzl)
       write(*,*)'done'
      endif
c
      if(zuv.eq.'u'.or.zuv.eq.'v')then
       allocate(depth(n,m),mz(n,m),glon(n,m),glat(n,m))
c currents case: need to read grid
       open(unit=1,file=modname,status='old')
       read(1,*)gname
       read(1,*)gname
       read(1,'(a80)')gname
       gname=rmCom(gname)
       close(1)
       status = nf_open(gname, nf_nowrite, ncid)
       if(status.ne.0)then
        write(*,*)'No grid file ',trim(gname),' found'
        stop
       endif
cccc READ GRID DIMENSIONS
       status=nf_inq_dimid(ncid,'nx',dimid)
       status=nf_inq_dimlen(ncid,dimid,n0)
       status=nf_inq_dimid(ncid,'ny',dimid)
       status=nf_inq_dimlen(ncid,dimid,m0)
       if(n0.ne.n.or.m0.ne.m)then
         write(*,*)'Wrong grid in ',trim(gname)
         write(*,*)'Grid size and model size are different!'
         stop
       endif
       call read2d_double_nc(ncid,'hz   ',n,m,depth)
       call read2d_int_nc(ncid,'mz   ',n,m,mz)
       call read2d_double_nc(ncid,'lat_z',n,m,glat)
       call read2d_double_nc(ncid,'lon_z',n,m,glon)

c       if(zuv.eq.'u')then
c          call read2d_double_nc(ncid,'hu   ',n,m,depth)
c          call read2d_int_nc(ncid,'mu   ',n,m,mz)
c       else
c           call read2d_double_nc(ncid,'hv   ',n,m,depth)
c           call read2d_int_nc(ncid,'mv   ',n,m,mz)
c       endif
       status=nf_close(ncid)
      endif
c mz now whatever we need, i.e. mu or mv
C
c output format
      fmt='(f9.3,f9.3,'
      c1=','
      do ic=1,ncon
       if(ic.eq.ncon)c1=')'
       if(APRI)then
        fmt=trim(fmt)//'f8.3,f8.1'//c1
       else
        fmt=trim(fmt)//'f8.3,f8.3'//c1
       endif
      enddo
c
      ftmp='   Lat     Lon  |'
      do ic=1,ncon
       if(APRI)then
        ftmp=trim(ftmp)//'   '//trim(c_id(ic))//'_amp  '//
     *                   trim(c_id(ic))//'_ph'
       else
        ftmp=trim(ftmp)//'   '//trim(c_id(ic))//'_Re   '//
     *                   trim(c_id(ic))//'_Im'
       endif
      enddo
      write(11,*)trim(ftmp)
c
      latp=0.
      lonp=0.
      do k=1,ndat
       if(lat(k).ne.latp.or.lon(k).ne.lonp)then
        call interp_r8(zRe,ncon,n,m,mask,mlon,mlat,
     *              lat(k),lon(k),dtmp,ierr)
        z1=dtmp
        call interp_r8(zIm,ncon,n,m,mask,mlon,mlat,
     *              lat(k),lon(k),dtmp,ierr)
        z1=cmplx(real(z1),real(dtmp))
        if(ierr.eq.0)then
         if(zuv.eq.'u'.or.zuv.eq.'v')then
           call interp_r8(depth,1,n,m,mz,glon,glat,
     *                 lat(k),lon(k),d1,ierr1)
           z1=z1/amax1(real(d1),1.)*100. ! currents cm/s
         elseif(zuv.eq.'z'.and.geo)then
           call interp_r8(zlRe,ncon,nl,ml,mzl,slon,slat,
     *                 lat(k),lon(k),dtmp,ierr1)
           zl1=dtmp
           call interp_r8(zlIm,ncon,nl,ml,mzl,slon,slat,
     *                 lat(k),lon(k),dtmp,ierr1)
           zl1=cmplx(real(zl1),real(dtmp))
           z1=z1+zl1     ! apply load correction to get geocentric tide
         endif  
        endif
        if(ierr.eq.0)then
         if(APRI)then
          write(11,fmt)lat(k),lon0(k),(abs(z1(ic)),atan2(-imag(z1(ic)),
     *                 real(z1(ic)))*180/3.141593,ic=1,ncon)
         else
           write(11,fmt)lat(k),lon0(k),
     *              (real(z1(ic)),imag(z1(ic)),ic=1,ncon)
         endif
        else
          write(11,'(1x,f8.3,f8.3,a70)')lat(k),lon(k),
     * '************* Site is out of model grid OR land ***************'
        endif
        latp=lat(k)
        lonp=lon(k)
       endif  
      enddo
      deallocate(zRe,zIm,z1,mask,cind,lat,lon,lon0)
      deallocate(mlon,mlat,dtmp)
      if(zuv.eq.'u'.or.zuv.eq.'v')deallocate(depth,mz,glon,glat)
      if(ibl.eq.1.and.geo.and.zuv.eq.'z')then
         ncon=0
         deallocate(zlRe,zlIm,zl1,mzl,lcind,slon,slat)
      endif
      enddo ! imod
      close(11)
      write(*,*)'Results are in ',trim(outname)
      if(zuv.eq.'z'.and.geo.and.ibl.eq.0)
     *   deallocate(zlRe,zlIm,zl1,mzl,lcind)
      stop
1     write(*,*)'Lat Lon file ',trim(lltname),' not found'
      write(*,*)'Check setup file setup, line 2.'
      stop
4     write(*,*)'Grid file ',trim(gname),' not found'
      write(*,*)'Check file ',trim(modname)
      stop
11    write(*,*)'File ''model.list'' was NOT found...'
      write(*,*)'TO CREATE please do:'
      write(*,*)'ls -1 DATA/Model_*>model.list'
      stop
      end

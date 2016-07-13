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
      program predict_tide
ccc   LANA, 2008 (remake of 2004 for netcdf format)
ccc   modified March 2011 to optimize for obtaining
ccc   time series at open boundaries
ccc
ccc   reads OTIS complex model file in netcdf format
ccc   (elevations OR transports), reads a list of locations,
ccc   reads list of times  
ccc   and outputs ASCII file with the tidal predictions of tidal 
ccc   elevations/transports/currents at the locations and times
ccc   
      implicit none
      include 'netcdf.inc'
      include 'constit.h'
      integer ncid,dimid,status
      real*8, allocatable:: zRe(:,:,:),zIm(:,:,:)
      real*8, allocatable:: zlRe(:,:,:),zlIm(:,:,:),depth(:,:)
      real*8, allocatable:: mlon(:,:),mlat(:,:),slon(:,:),slat(:,:)
      real*8, allocatable:: mlonu(:,:),mlatu(:,:),mlonv(:,:),mlatv(:,:)
      real*8, allocatable:: lat(:),lon(:),lon0(:),dtmp(:)
      real, allocatable:: zpred(:),upred(:),vpred(:)
      real*8, allocatable:: uRe(:,:,:),uIm(:,:,:),
     *                      vRe(:,:,:),vIm(:,:,:)
      complex, allocatable:: z1(:),zl1(:),u1(:),v1(:)
      real*8, allocatable:: time_mjd(:)
      real*8 d1
      real th_lim(2),ph_lim(2),dum,lth_lim(2),lph_lim(2)
      integer, allocatable:: mask(:,:),cind(:),lcind(:),ccind(:),
     *                       mz(:,:),mzl(:,:)
      integer, allocatable:: masku(:,:),maskv(:,:)
c
      character*4 c_id(ncmx),c_id_mod(ncmx),lc_id(ncmx),tcon(ncmx)
      character*80 modname,lltname,outname,gname,ctmp,lname
      character*2000 fmt
      character*80 rmCom,tfname
      character*1 zuv,c1,c2
      character*80 arg
      character*10 cdate
      character*8 ctime
      character*10 deblank 
      logical APRI,geo,interp_micon,ll_km,ltmp
      integer ncon,nc,n,m,ndat,i,j,k,k1,ierr,ierr1,ic,n0,m0
      integer ncl,nl,ml,nmod,imod,ibl,ntime,idum,mjd,julian
      integer yyyy1,mm1,dd1,iargc,narg,it
      integer, allocatable:: yyyy(:),mm(:),dd(:),hh(:),
     *                       mi(:),ss(:)
c
      narg=iargc()
      if(narg.gt.0)then
        call getarg(1,arg)
        read(arg(3:80),'(a80)')tfname
        ntime=0
        open(unit=1,file=tfname,status='old',err=16)
18      read(1,*,end=17)idum,idum,idum,idum,idum,idum
        ntime=ntime+1
        go to 18
17      rewind(1)
        ntime=ntime-1
        allocate(yyyy(ntime),mm(ntime),dd(ntime),
     *         hh(ntime),mi(ntime),ss(ntime),time_mjd(ntime))
        call read_time(ntime,1,2,yyyy,mm,dd,hh,mi,ss,time_mjd)
        close(1)
      endif       
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
        write(*,*)'Models to predict tides from:'
13      read(12,'(a80)',end=12)ctmp
        write(*,*)trim(ctmp)
        nmod=nmod+1
        go to 13
12      rewind(12)
        read(12,'(a80)')modname
      endif
      write(*,*)
      write(*,*)'Lat/Lon/Time file:',trim(lltname)
      if(ncon.gt.0)write(*,*)'Constituents to include: ',c_id(1:ncon)
      if(zuv.eq.'z')then
       if(geo)then
         write(*,*)'Predict GEOCENTRIC tide'
       else
         write(*,*)'Predict OCEAN tide'
       endif
      endif
      if(interp_micon)write(*,*)'Interpolate minor constituents'
c
      if(narg.eq.0)then ! default usage
c read times
       ntime=0
       open(unit=1,file=lltname,status='old',err=6)
8      read(1,*,end=7)dum,dum,idum,idum,idum,idum,idum,idum
       ntime=ntime+1
       go to 8
7      rewind(1)
       allocate(yyyy(ntime),mm(ntime),dd(ntime),
     *          hh(ntime),mi(ntime),ss(ntime),time_mjd(ntime))
       call read_time(ntime,1,1,yyyy,mm,dd,hh,mi,ss,time_mjd)
       close(1)
      endif
c
      open(unit=11,file=outname,status='unknown')
      do imod=1,nmod
      write(*,*)
      if(imod.gt.1)read(12,'(a80)')modname
      call rd_mod_header_nc(modname,zuv,n,m,
     *          th_lim,ph_lim,nc,c_id_mod,ll_km)
      write(*,*)'Model:        ',trim(modname(12:80))
      write(11,'(60a1)')('-',i=1,60)
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
      if(zuv.eq.'z')then
        write(*,*)'Predict elevations (m)'
      else
        write(*,*)'Predict transport (m^2/s) and currents (cm/s)'
      endif
c
      k1=1
      tcon=''  
      if(ncon.eq.0)then
       ibl=1
       ncon=nc
       c_id=c_id_mod
      else
c check if all required constituents are in the model
       do ic=1,ncon
        do k=1,nc
         if(c_id(ic).eq.c_id_mod(k))then
          tcon(k1)=c_id(ic)
          k1=k1+1
          go to 14
         endif
        enddo
        write(*,*)'Constituent ',c_id(ic), ' is NOT in the model'
14      continue
       enddo 
       ncon=k1-1
       c_id=tcon
      endif
 
      write(*,*)'Constituents to include: ',c_id(1:ncon)
      write(11,*)'Constituents included: ',c_id(1:ncon)
c
      allocate(cind(ncon),ccind(ncon))
      call def_con_ind(c_id,ncon,c_id_mod,nc,cind)
c find corresponding indices in constit.h
      call def_cid(ncon,c_id,ccind)
c
      if(narg.eq.0)then
       ndat=ntime
      else
       ndat=0
       open(unit=1,file=lltname,status='old',err=6)
 20    read(1,*,end=19)dum,dum
       ndat=ndat+1
       go to 20
 19    close(1)
      endif
c
      open(unit=1,file=lltname,status='old',err=6)
      allocate(lat(ndat),lon(ndat),lon0(ndat))
      if(zuv.eq.'z')then
       allocate(zpred(ndat))
      else
       allocate(upred(ndat),vpred(ndat))
      endif
      do k=1,ndat
       read(1,*)lat(k),lon(k) ! ignore rest of record
       lon0(k)=lon(k)
      enddo
      close(1)
c
      allocate(mlon(n,m),mlat(n,m),dtmp(ncon))
      if(zuv.eq.'z')then      
       allocate(zRe(ncon,n,m),zIm(ncon,n,m),z1(ncon),mask(n,m))
       if(geo) then
        ltmp=.false.
        call rd_mod_header_nc(lname,'z',nl,ml,
     *                    lth_lim,lph_lim,ncl,lc_id,ltmp)
        allocate(lcind(ncon))
        call def_con_ind(c_id,ncon,lc_id,ncl,lcind)
        allocate(zlRe(ncon,nl,ml),zlIm(ncon,nl,ml),zl1(ncon),
     *           mzl(nl,ml),slon(nl,ml),slat(nl,ml))
        write(*,*)
        write(*,'(a,$)')' Reading file to apply load corrections...'
        call rd_mod_body_nc(lname,'z',lcind,ncl,ncon,
     *                     slon,slat,zlRe,zlIm,nl,ml,mzl)
        write(*,*)'done'
       endif
      else
       allocate(uRe(ncon,n,m),uIm(ncon,n,m),vRe(ncon,n,m),
     *          vIm(ncon,n,m),u1(ncon),v1(ncon),
     *          masku(n,m),maskv(n,m))
       allocate(mlonu(n,m),mlatu(n,m),
     *          mlonv(n,m),mlatv(n,m))
      endif
c
      write(*,'(a,$)')' Reading model...'
      if(zuv.eq.'z')then
       call rd_mod_body_nc(modname,'z',cind,nc,ncon,
     *                    mlon,mlat,zRe,zIm,n,m,mask)
      else
       call rd_mod_body_nc(modname,'u',cind,nc,ncon,
     *                    mlonu,mlatu,uRe,uIm,n,m,masku)
       call rd_mod_body_nc(modname,'v',cind,nc,ncon,
     *                    mlonv,mlatv,vRe,vIm,n,m,maskv)
      endif
      write(*,*)'done'
c
c      if(zuv.ne.'z')then
       allocate(depth(n,m),mz(n,m))
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
       if(zuv.ne.'z')then
        call read2d_double_nc(ncid,'lon_z',n,m,mlon)
        call read2d_double_nc(ncid,'lat_z',n,m,mlat)
       endif
       status=nf_close(ncid)
c      endif
      ctmp='    Lat       Lon        mm.dd.yyyy hh:mm:ss'
      write(11,*)'' 
      if(zuv.eq.'z')then
        write(11,*)trim(ctmp),'     z(m)   Depth(m)'
      else
        write(11,*)trim(trim(ctmp)//'   U(m^2/s)  V(m^2/s)'//
     *                        '   u(cm/s)   v(cm/s) Depth(m)')
      endif
      write(11,*)''
c
      do k=1,ndat
        if(zuv.eq.'z')then
          call interp_r8(zRe,ncon,n,m,mask,mlon,mlat,
     *                   lat(k),lon(k),dtmp,ierr)
          z1=dtmp
          call interp_r8(zIm,ncon,n,m,mask,mlon,mlat,
     *              lat(k),lon(k),dtmp,ierr)
          z1=cmplx(real(z1),real(dtmp))
         else
          call interp_r8(uRe,ncon,n,m,masku,mlonu,mlatu,
     *                   lat(k),lon(k),dtmp,ierr)
          u1=dtmp
          call interp_r8(uIm,ncon,n,m,masku,mlonu,mlatu,
     *              lat(k),lon(k),dtmp,ierr)
          u1=cmplx(real(u1),real(dtmp))
          call interp_r8(vRe,ncon,n,m,maskv,mlonv,mlatv,
     *                   lat(k),lon(k),dtmp,ierr)
          v1=dtmp
          call interp_r8(vIm,ncon,n,m,maskv,mlonv,mlatv,
     *              lat(k),lon(k),dtmp,ierr)
          v1=cmplx(real(v1),real(dtmp))
        endif
        if(ierr.eq.0)then
          call interp_r8(depth,1,n,m,mz,mlon,mlat,
     *              lat(k),lon(k),d1,ierr)
         if(zuv.eq.'z'.and.geo)then    
           call interp_r8(zlRe,ncon,nl,ml,mzl,slon,slat,
     *                 lat(k),lon(k),dtmp,ierr1)
           zl1=dtmp
           call interp_r8(zlIm,ncon,nl,ml,mzl,slon,slat,
     *                 lat(k),lon(k),dtmp,ierr1)
           zl1=cmplx(real(zl1),real(dtmp))
           z1=z1+zl1     ! apply load correction to get geocentric tide
         endif  
c predict tide
         if(narg.eq.0)then ! default usage
          if(zuv.eq.'z')then
           call ptide(z1,c_id,ncon,ccind,lat(k),time_mjd(k),1,
     *               interp_micon,zpred(k))
          else
           call ptide(u1,c_id,ncon,ccind,lat(k),time_mjd(k),1,
     *               interp_micon,upred(k))
           call ptide(v1,c_id,ncon,ccind,lat(k),time_mjd(k),1,
     *               interp_micon,vpred(k))
          endif
c
           write(cdate,'(i2,a1,i2,a1,i2)')hh(k),':',mi(k),':',ss(k)
           ctime=deblank(cdate)
           write(cdate,'(i2,a1,i2,a1,i4)')mm(k),'.',dd(k),'.',yyyy(k)
           cdate=deblank(cdate)
           if(zuv.eq.'z')then
            write(11,'(1x,f10.4,f10.4,5x,a10,1x,a8,f10.3,f10.3)')
     *       lat(k),lon0(k),cdate,ctime,zpred(k),real(d1)
           else
            write(11,'(1x,f10.4,f10.4,5x,a10,1x,a8,5(f10.3))')
     *       lat(k),lon0(k),cdate,ctime,upred(k),vpred(k),
     *       upred(k)/real(d1)*100,vpred(k)/real(d1)*100,real(d1)
           endif
c OB usage March 2011
         else
          do it=1,ntime
           if(zuv.eq.'z')then
            call ptide(z1,c_id,ncon,ccind,lat(k),time_mjd(it),1,
     *                interp_micon,zpred(k))
           else
            call ptide(u1,c_id,ncon,ccind,lat(k),time_mjd(it),1,
     *                interp_micon,upred(k))
            call ptide(v1,c_id,ncon,ccind,lat(k),time_mjd(it),1,
     *               interp_micon,vpred(k))
           endif
c
           write(cdate,'(i2,a1,i2,a1,i2)')hh(it),':',mi(it),':',ss(it)
           ctime=deblank(cdate)
           write(cdate,'(i2,a1,i2,a1,i4)')mm(it),'.',dd(it),'.',yyyy(it)
           cdate=deblank(cdate)
           if(it.eq.1)then
            write(11,'(1x,f10.4,f10.4)')lat(k),lon0(k)
           endif
           if(zuv.eq.'z')then
            write(11,'(26x,a10,1x,a8,f10.3,f10.3)')
     *            cdate,ctime,zpred(k),real(d1)
           else
            write(11,'(26x,a10,1x,a8,5(f10.3))')
     *            cdate,ctime,upred(k),vpred(k),
     *       upred(k)/real(d1)*100,vpred(k)/real(d1)*100,real(d1)
           endif
          enddo
         endif
        else
          write(11,'(1x,f10.4,f10.4,a)')lat(k),lon(k),
     * '***** Site is out of model grid OR land *****'
        endif  
      enddo
      deallocate(cind,ccind,lat,lon,lon0,dtmp,mz,mlon,mlat)
      if(zuv.eq.'z')then
       deallocate(zRe,zIm,z1,mask,zpred)
      else
       deallocate(uRe,uIm,vRe,vIm,u1,v1,masku,maskv,upred,vpred)
      endif
      if(zuv.eq.'z'.and.ibl.eq.1.and.geo)then
         ncon=0
         deallocate(zlRe,zlIm,zl1,mzl,lcind,slon,slat)
      endif
      enddo ! imod
      close(11)
      close(12)
      write(*,*)'Results are in ',trim(outname)
      stop
1     write(*,*)'Lat lon file ',trim(lltname),' not found'
      write(*,*)'Check setup file, line 2.'
      stop
4     write(*,*)'Grid file ',trim(gname),' not found'
      write(*,*)'Check file ',trim(modname),', line 3'
      stop
6     write(*,*)'File ',trim(lltname),' not found'
      write(*,*)'Check setup file, line 7.'
      stop
16    write(*,*)'File ',trim(tfname),' not found'
      write(*,*)'Check spelling in command line'
      call usage()
      stop
11    write(*,*)'File ''model.list'' was NOT found...'
      write(*,*)'TO CREATE please do:'
      write(*,*)'ls -1 DATA/Model_*>model.list'
      stop
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine usage()
      write(*,*)'Usage:'
      write(*,*)'predict_tide [-t<time_file>]<setup.inp'
      write(*,*)'Default: lat_lon_time file is used for input'
      write(*,*)'          if option -t is given, then'
      write(*,*)'          lats/lons are read from lat_lon file'
      write(*,*)'          set in setup.inp, snd times are read'
      write(*,*)'          from <time_file>.'
      write(*,*)'          Use, when need output for time series'
      write(*,*)'          for open boundaries, i.e. times are'
      write(*,*)'          the same in all nodes'
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine read_time(ntime,iunit,iopt,yyyy,mm,dd,hh,
     *                     mi,ss,time_mjd)
      implicit none
      integer k,iunit,iopt,julian,mjd,mm1,dd1,yyyy1
      integer ntime,yyyy(ntime),mm(ntime),dd(ntime)
      integer hh(ntime),mi(ntime),ss(ntime)
      real dum
      real*8 time_mjd(ntime)
      character*10 cdate,deblank
c
      do k=1,ntime
       if(iopt.eq.1)then
        read(iunit,*)dum,dum,yyyy(k),mm(k),dd(k),hh(k),mi(k),ss(k)
       else
        read(iunit,*)yyyy(k),mm(k),dd(k),hh(k),mi(k),ss(k)
       endif
c convert to mjd
       call date_mjd(mm(k),dd(k),yyyy(k),mjd)
c check if exists such a date
        julian=mjd+2400001        
        call CALDAT (JULIAN,MM1,DD1,YYYY1)
        if(mm(k).ne.mm1.or.dd(k).ne.dd1.or.yyyy(k).ne.yyyy1)then
         write(cdate,'(i2,a1,i2,a1,i4)')mm(k),'.',dd(k),'.',yyyy(k)
         cdate=deblank(cdate)
         write(*,*)'Wrong date in (lat_lon)_time file:',cdate
         stop
        endif
        time_mjd(k)=dble(mjd)+dble(hh(k))/24.D0+
     *              dble(mi(k))/(24.D0*60.D0)  +
     *              dble(ss(k))/(24.D0*60.D0*60.D0)
      enddo
      return
      end

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
      program extract_HC
!cc   LANA, 2012 remake to provide for "atlas" models 
!cc   reads OTIS standard binary complex model file
!cc   (elevations OR transports), reads a list of locations 
!cc   and outputs ASCII file with the complex amplitudes/amp,phases of
!cc   elevations/transports/currents interpolated at the locations
! 
      implicit none
      include 'derived_types.h'
      include 'constit.h'
      type (c_grid), pointer:: grid0,grid(:)
      type (constituents) cid,cid0
      type (z_model), pointer:: Mod_z0,Mod_z(:),Mod_load
      type (uv_model), pointer:: Mod_uv0,Mod_uv(:)
      type (LocDict), pointer:: dloc(:)
      integer*4 nob,nlines,nz,nu,nv
!
      complex*8, allocatable:: z1(:,:),zl1(:,:),u1(:,:),r(:,:),r0(:,:)
      complex*8, allocatable:: z2(:,:),u2(:,:)
      real*4, allocatable:: phase(:)
      real*4 lon,lat,x,y,latp,lonp
      real*8 dlat,dlon
      real*4 th_lim(2),ph_lim(2)
      integer*4, allocatable:: cind(:),lcind(:),pmask(:,:),mid(:),id(:),mid1(:)
      integer*4, allocatable:: cind2(:)
!
      character*4 c_id(ncmx),c_id_mod(ncmx),lc_id(ncmx),c_id_mod1(ncmx)
      character*80 modname,lltname,outname,ctmp(1,3),lname,rmCom
      character*80 hname,uname,gname,ll_km_sub
      character*2000 fmt
      character*20, allocatable:: mnames(:)
      character*1 zuv,c1
      logical APRI,geo,atlas,ll_km,interp_micon
      integer*4 ncon,nc,n,m,ndat,i,k,l,nmod,nmodi,midp,ln
      integer(kind=4) ic,ncl,nl,ml,gid,n1,m1,nc1,nob1
!
      atlas=.false.
      lname='DATA/load_file'
      call rd_inp(modname,lltname,zuv,c_id,ncon,APRI,geo, &
                  outname,interp_micon)
      open(unit=11,file=outname,status='unknown')
      write(*,*)trim(modname)
      write(11,*)trim(modname)
      write(*,*)'Lat/Lon file:',trim(lltname)
      if(ncon.gt.0)write(*,*)'Constituents to include: ',c_id(1:ncon)
      if(geo.and.zuv.eq.'z')then
         write(*,*)'Extract GEOCENTRIC tide HC'
      else
         write(*,*)'Extract OCEAN tide HC'
      endif
! read Model_* file
       ll_km_sub=''
       open(unit=1,file=modname,status='old')
       read(1,'(a80)')hname
       read(1,'(a80)')uname
       read(1,'(a80)')gname
       read(1,*,end=13,err=13)ll_km_sub
13     hname=rmCom(hname)
       uname=rmCom(uname)
       gname=rmCom(gname)
       close(1)
       if(trim(ll_km_sub).ne.'')then
        open(unit=1,file='subname',status='unknown');
        write(1,*)trim(ll_km_sub)
        close(1)
       endif
! get model dimensions and limits
       allocate(grid0)
       call getDims(gname,grid0)
       th_lim=grid0%theta_lim;ph_lim=grid0%phi_lim;
       n=grid0%n;m=grid0%m;
       nob=max(1,grid0%nob);ll_km=grid0%ll_km;
       gid=grid0%km_area_id
       allocate(grid0%hz(n,m),grid0%mz(n,m),grid0%iob(2,nob))
       grid0%allocated=.true.
       call grid_in(gname,grid0,.true.)
!
      allocate(pmask(n,m),mnames(1))
      call atlas_mask_in(gname,pmask,n,m,atlas,nmod,mnames,1)
      write(*,'(a,$)')' Model type:'
      if(atlas)then
       write(*,*)'atlas (combined of several models)'
       write(*,*)'Models in atlas:',nmod
       deallocate(mnames);allocate(mnames(nmod))
       ln=nmod
       call atlas_mask_in(gname,pmask,n,m,atlas,nmod,mnames,ln)
      else
       write(*,*)'single model'
      endif
6     if(ll_km)then
         write(*,*)'Approxmate area limits:' 
      endif
      write(*,*)'Lat limits:   ',th_lim
      write(*,*)'Lon limits:   ',ph_lim
! Find nc
      allocate(Mod_z0)
      call getDims_model(hname,Mod_z0)
      nc=Mod_z0%nc
!
      if(zuv.eq.'z')then
         write(*,*)'Output elevations (m)'
         write(11,*)'Elevations (m)'
         allocate(Mod_z0%z(nc,n,m),Mod_z0%mz(n,m),Mod_z0%cid)
         allocate(Mod_z0%cid%cons(nc),Mod_z0%cid%idc(nc),Mod_z0%cid%group(nc))
         call loadModel_z(Mod_z0,hname)
         if(ll_km)then
          Mod_z0%theta_lim=th_lim
          Mod_z0%phi_lim=ph_lim
         endif
         call def_cid(nc,Mod_z0%cid) ! result in Mod_z0%cid%idc
         c_id_mod(1:nc)=Mod_z0%cid%cons
      else
         deallocate(Mod_z0)
         if(zuv.eq.'U')write(*,*)'Output WE transport (m^2/s)'
         if(zuv.eq.'u')write(*,*)'Output WE velocity  (cm/s)'
         if(zuv.eq.'V')write(*,*)'Output SN transport (m^2/s)'
         if(zuv.eq.'v')write(*,*)'Output SN velocity  (cm/s)'
         if(zuv.eq.'U')write(11,*)'WE transport (m^2/s)'
         if(zuv.eq.'u')write(11,*)'WE velocity  (cm/s)'
         if(zuv.eq.'V')write(11,*)'SN transport (m^2/s)'
         if(zuv.eq.'v')write(11,*)'SN velocity  (cm/s)'
         allocate(Mod_uv0)
         allocate(Mod_uv0%u(nc,n,m),Mod_uv0%v(nc,n,m), &
                  Mod_uv0%mu(n,m),Mod_uv0%mv(n,m),Mod_uv0%cid) 
         allocate(Mod_uv0%cid%cons(nc),Mod_uv0%cid%idc(nc),Mod_uv0%cid%group(nc))
         Mod_uv0%ll_km=ll_km
         call loadModel_uv(Mod_uv0,uname)
         if(ll_km)then
          Mod_uv0%theta_lim=th_lim
          Mod_uv0%phi_lim=ph_lim
         endif
         call def_cid(nc,Mod_uv0%cid) ! result in Mod_z0%cid%idc
         c_id_mod(1:nc)=Mod_uv0%cid%cons
      endif

      write(*,*)'Constituents: ',c_id_mod(1:nc)
!
      if(ncon.eq.0)then
       ncon=nc
       c_id=c_id_mod
      endif
      write(*,*)'Constituents to include: ',c_id(1:ncon)
!
      allocate(cind(ncon))
      call def_con_ind(c_id,ncon,c_id_mod,nc,cind)
!
      ndat=0;nlines=0;
      latp=0;lonp=0;
      open(unit=1,file=lltname,status='old',err=1)
3     read(1,*,end=2)lat,lon
      nlines=nlines+1
      if(lat.ne.latp.or.lon.ne.lonp)then
        ndat=ndat+1
        latp=lat;lonp=lon
      else
       !write(*,*)'Info: line is doubled in ',trim(lltname),':',lat,lon
      endif
      go to 3
2     rewind(1)
      allocate(dloc(ndat))
      allocate(mid(ndat),id(ndat),mid1(ndat))
      mid=1
      latp=0;lonp=0;k=0;nmodi=0
      !write(*,*)'max pmask:',maxval(pmask)
      do l=1,nlines
       read(1,*)lat,lon
       if(latp.ne.lat.or.lonp.ne.lon)then
        k=k+1
        dloc(k)%lat(1)=lat;dloc(k)%lon(1)=lon
        dloc(k)%theta=0;dloc(k)%phi=0;dloc(k)%land=.false.
        if(atlas)then
           call def_model(lat,lon,th_lim,ph_lim,n,m,pmask,mid(k))
           if(mid(k).gt.0)then
            do i=1,k-1
             if(mid(k).eq.mid(i))go to 7
            enddo
            nmodi=nmodi+1
            id(nmodi)=mid(k)
7           continue
           endif
        endif
        x=0;y=0;latp=lat;lonp=lon;dlat=lat;dlon=lon
        if(ll_km)call convert_ll_xy(gid,dlat,dlon,x,y)
        dloc(k)%x(1)=x;dloc(k)%y(1)=y        
        !write(*,*)dloc(k)%lat(1),dloc(k)%lon(1)
       endif
      enddo ! l
      close(1)
      mid1=mid ! save mid in mid1, since below mid is redefined
!
      if(zuv.eq.'z'.and.geo)then
! Read loading & self attraction
        write(*,'(a,$)')' Reading Ocean load&self attraction file...'
        allocate(Mod_load)
        call getDims_model(lname,Mod_load)
        ncl=Mod_load%nc;nl=Mod_load%n;ml=Mod_load%m
        allocate(Mod_load%z(ncl,nl,ml),Mod_load%mz(nl,ml),Mod_load%cid)
        allocate(Mod_load%cid%cons(ncl),Mod_load%cid%idc(ncl),Mod_load%cid%group(ncl))
        call loadModel_z(Mod_load,lname)
        call def_cid(ncl,Mod_load%cid) ! result in Mod_z0%cid%idc
        lc_id(1:ncl)=Mod_load%cid%cons
        allocate(lcind(ncon))
        call def_con_ind(c_id,ncon,lc_id,ncl,lcind)
        allocate(zl1(ncl,ndat))
        write(*,*)'done'
      endif
!
      if(nmodi.gt.0)then ! atlas case
       allocate(grid(nmodi))
       do k=1,nmodi
        write(*,*)'Loading local solution grid ',mnames(id(k))
        call getDims_l(gname,id(k),grid(k),nz,1)
        if(grid(k)%n.eq.0)then ! error reading local grid
          grid(k)%n=1;grid(k)%m=1;grid(k)%nob=1;
        endif
        n1=grid(k)%n;m1=grid(k)%m;
        nob1=max(1,grid(k)%nob);
        allocate(grid(k)%hz(n1,m1),grid(k)%mz(n1,m1),grid(k)%iob(2,nob1))
        grid(k)%allocated=.true.
        if(n1*m1.gt.1)call grid_in_l(1,grid(k),nz)
       enddo
       if(zuv.eq.'z')then
         allocate(Mod_z(nmodi))
         do k=1,nmodi
           write(*,*)'Loading local solution elevations ',mnames(id(k))
           call getDims_z_l(hname,id(k),nz,8,Mod_z(k))
           if(Mod_z(k)%n.eq.0)then ! Model not found
            where(mid.eq.id(k))mid=0
            id(k)=0;Mod_z(k)%n=1;Mod_z(k)%m=1;Mod_z(k)%nc=1;     
           endif
           nc1=Mod_z(k)%nc;n1=Mod_z(k)%n;m1=Mod_z(k)%m;
           allocate(Mod_z(k)%z(nc1,n1,m1),Mod_z(k)%mz(n1,m1),Mod_z(k)%cid)
           allocate(Mod_z(k)%cid%cons(nc1),Mod_z(k)%cid%idc(nc1),Mod_z(k)%cid%group(nc1))
           if(id(k).gt.0)then
            call loadModel_z_l(8,Mod_z(k),nz)
            call def_cid(nc1,Mod_z(k)%cid) ! result in Mod_z(k)%cid%idc
           endif
         enddo
        else
         allocate(Mod_uv(nmodi))
         do k=1,nmodi
          write(*,*)'Loading local solution transports ',mnames(id(k))
          call getDims_u_l(uname,id(k),nu,nv,8,Mod_uv(k))
          if(Mod_uv(k)%n.eq.0)then ! Model not found
            where(mid.eq.id(k))mid=0
            id(k)=0;Mod_uv(k)%n=1;Mod_uv(k)%m=1;Mod_uv(k)%nc=1;     
          endif
          nc1=Mod_uv(k)%nc;n1=Mod_uv(k)%n;m1=Mod_uv(k)%m;
          allocate(Mod_uv(k)%u(nc1,n1,m1),Mod_uv(k)%v(nc1,n1,m1), &
                   Mod_uv(k)%mu(n1,m1),Mod_uv(k)%mv(n1,m1),Mod_uv(k)%cid) 
          allocate(Mod_uv(k)%cid%cons(nc1),Mod_uv(k)%cid%idc(nc1),Mod_uv(k)%cid%group(nc1))
          if(id(k).gt.0)then
           call loadModel_uv_l(8,Mod_uv(k),nu,nv)
           call def_cid(nc1,Mod_uv(k)%cid) ! result in Mod_uv(k)%cid%idc
          endif
         enddo
        endif
     endif
!
8      allocate(phase(ncon))
! output format
      fmt='(f9.4,x,f9.4,'
      c1=','
      do ic=1,ncon
       if(ic.eq.ncon)c1=')'
       if(APRI)then
        fmt=trim(fmt)//'f8.3,f8.1'//c1
       else
        fmt=trim(fmt)//'f8.3,f8.3'//c1
       endif
      enddo
!
      if(APRI)then
        write(11,*)'  Lat     Lon       ', &
                    (trim(c_id(ic)),'_amp  ', &
                     trim(c_id(ic)),'_ph   ',ic=1,ncon)
      else
        write(11,*)'  Lat       Lon     ', &
                   (trim(c_id(ic)),'_Re   ', &
                    trim(c_id(ic)),'_Im   ',ic=1,ncon)
      endif
!
      allocate(z1(nc,ndat),u1(nc,ndat),r(ncon,ndat),r0(nc,ndat))
      select case(zuv)
      case ('z')
       z1=0
       call interp_driver_z(Mod_z0,dloc,ndat,z1,'y')
       r0=z1
!
       !write(*,*)mid1
       do k=1,nmodi
        if(id(k).gt.0)then
         allocate(z2(Mod_z(k)%nc,ndat),cind2(nc))
         do l=1,ndat
          dloc(l)%land=.false.
         enddo
         call interp_driver_z(Mod_z(k),dloc,ndat,z2,'y')
         !do l=1,ndat
         ! if(dloc(l)%land)write(*,*)k,l,dloc(l)%lat(1),dloc(l)%lon(1)
         !enddo
         call def_con_ind(c_id_mod,nc,Mod_z(k)%cid%cons,Mod_z(k)%nc,cind2)
         do l=1,ndat
          do ic=1,nc
           if(.not.dloc(l)%land)then
            if(cind2(ic).gt.0)then
             if(z2(cind2(ic),l).ne.0)then
                z1(ic,l)=z2(cind2(ic),l)
                r0(ic,l)=0
             endif
            endif
           else
            if(z1(ic,l).ne.0.and.mid1(l).ne.id(k))dloc(l)%land=.false.       
           endif
          enddo
         enddo
         deallocate(z2,cind2)
        endif        
       enddo
       do ic=1,ncon
         if(cind(ic).ne.0)r(ic,:)=z1(cind(ic),:)  
       enddo
       if(geo)then
        zl1=0
        call interp_driver_z(Mod_load,dloc,ndat,zl1,'n')
        do ic=1,ncon
         if(lcind(ic).ne.0)r(ic,:)=r(ic,:)+zl1(lcind(ic),:)
        enddo
       endif
      case ('u','v','U','V')
       u1=0.
       if(zuv.eq.'u'.or.zuv.eq.'U')then
        dloc(:)%theta=1;dloc(:)%phi=0
       endif
       if(zuv.eq.'v'.or.zuv.eq.'V')then
        dloc(:)%theta=0;dloc(:)%phi=1;
       endif
       if(zuv.eq.'u'.or.zuv.eq.'v')call interp_driver_uv(Mod_uv0,dloc,grid0,ndat,u1,'y')
       if(zuv.eq.'U'.or.zuv.eq.'V')call interp_driver_uv(Mod_uv0,dloc,grid0,ndat,u1,'T')
       if(zuv.eq.'u'.or.zuv.eq.'v')u1=u1*100
       r0=u1
       do k=1,nmodi
        l=grid(k)%n*grid(k)%m
        if(id(k).gt.0.and.l.gt.1)then
         allocate(u2(Mod_uv(k)%nc,ndat),cind2(nc))
         do l=1,ndat
          dloc(l)%land=.false.
         enddo
         if(zuv.eq.'u'.or.zuv.eq.'v')call interp_driver_uv(Mod_uv(k),dloc,grid(k),ndat,u2,'y')
         if(zuv.eq.'U'.or.zuv.eq.'V')call interp_driver_uv(Mod_uv(k),dloc,grid(k),ndat,u2,'T')
         if(zuv.eq.'u'.or.zuv.eq.'v')u2=u2*100
         call def_con_ind(c_id_mod,nc,Mod_uv(k)%cid%cons,Mod_uv(k)%nc,cind2)
         do l=1,ndat
          do ic=1,nc
           if(.not.dloc(l)%land)then
            if(cind2(ic).gt.0)then
             if(u2(cind2(ic),l).ne.0)then
              u1(ic,l)=u2(cind2(ic),l)
              r0(ic,l)=0.
             endif
            endif
           else
            if(u1(ic,l).ne.0.and.mid1(l).ne.id(k))dloc(l)%land=.false.
           endif
          enddo
         enddo
         deallocate(u2,cind2)
        endif        
       enddo
       do ic=1,ncon
        if(cind(ic).ne.0)r(ic,:)=u1(cind(ic),:)
       enddo
      endselect
!
      midp=-1
      do k=1,ndat
       if(.not.dloc(k)%land)then
        if(atlas)then
         if(mid(k).ne.midp)then
          midp=mid(k)
          if(mid(k).gt.0.and.r0(1,k).eq.0)then
           write(11,*)'HC extracted from local solution :',mnames(mid(k))
          else
           write(11,*)'HC extracted from global solution'
          endif
         endif
        endif
        if(APRI)then
          phase=atan2(-imag(r(:,k)),real(r(:,k)))*180/3.141593
          write(11,fmt)dloc(k)%lat(1),dloc(k)%lon(1), &
                (abs(r(ic,k)),phase(ic),ic=1,ncon)
        else
           write(11,fmt)dloc(k)%lat(1),dloc(k)%lon(1), &
                    (real(r(ic,k)),imag(r(ic,k)),ic=1,ncon)
        endif
       else
          write(11,'(1x,f8.3,f8.3,a70)')dloc(k)%lat(1),dloc(k)%lon(1), &
       '************* Site is out of model grid OR land ***************'
       endif
      enddo
      deallocate(z1,u1,dloc,grid0,r,r0)
      if(geo.and.zuv.eq.'z')deallocate(Mod_load,zl1)
      if(zuv.eq.'z')then
        deallocate(Mod_z0)
      else
        deallocate(Mod_uv0)
      endif
      close(11)
      write(*,*)'Results are in ',trim(outname)
      stop
1     write(*,*)'Lat Lon file ',trim(lltname),' not found'
      write(*,*)'Check setup file, line 2.'
      stop
4     write(*,*)'Grid file ',trim(gname),' not found'
      write(*,*)'Check file ',trim(modname)
      stop
      end

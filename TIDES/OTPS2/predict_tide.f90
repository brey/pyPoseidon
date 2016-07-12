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
      program predict_tide
!cc   LANA, 2012 remake to provide for "atlas" models 
!cc   modified March 2011 to optimize for obtaining
!cc   time series at open boundaries
!cc
!cc   reads OTIS standard binary complex model file
!cc   (elevations OR transports), reads a list of locations,
!cc   reads list of times  
!cc   and outputs ASCII file with the tidal predictions of tidal 
!cc   elevations/transports/currents at the locations and times
!cc   
      implicit none
      include 'derived_types.h'
      include 'constit.h'
      type (c_grid), pointer:: grid0,grid(:)
      type (constituents) cid,cid0
      type (z_model), pointer:: Mod_z0,Mod_z(:),Mod_load
      type (uv_model), pointer:: Mod_uv0,Mod_uv(:)
      type (LocDict), pointer:: dloc(:)
!
      complex, allocatable:: z1(:,:),zl1(:,:),r(:,:),z2(:,:)
      complex, allocatable:: u1(:,:),v1(:,:),ut1(:,:),vt1(:,:)
      complex, allocatable:: ru(:,:),rv(:,:),rut(:,:),rvt(:,:)
      complex, allocatable:: u2(:,:),v2(:,:),ut2(:,:),vt2(:,:)
      real*4, allocatable:: zpred(:),upred(:),vpred(:),d(:),d2(:)
      real*8, allocatable:: time_mjd(:),t1(:)
      real*4 th_lim(2),ph_lim(2),d1,dum,latp,lonp,lat,lon,x,y
      real*8 dlat,dlon
      integer*4, allocatable:: pmask(:,:),cind(:),lcind(:),ccind(:),&
                             mid(:),id(:),cind2(:),tid(:)
!
      character*4 c_id(ncmx),c_id_mod(ncmx),lc_id(ncmx),tcon(ncmx)
      character*4 c_id_mod1(ncmx)
      character*80 modname,lltname,outname,ctmp(1,3)
      character*80 lname,rmCom,tfname,arg
      character*80 gname,hname,uname,deblank,ll_km_sub
      character*20, allocatable:: mnames(:)
      character*2000 fmt
      character*1 zuv,c1
      character*10 cdate
      character*8 ctime
      logical APRI,geo,atlas,ll_km,interp_micon
      integer*4 ncon,nc,n,m,ndat,i,j,k,l,ic,nmod,nmodi,midp,nlines,nz
      integer*4 ncl,nl,ml,imod,ibl,ntime,idum,mjd,julian,narg,nu,nv,ln
      integer*4 yyyy1,mm1,dd1,iargc,gid,n1,m1,nc1,nob,nob1,it,nt
      integer*4, allocatable:: yyyy(:),mm(:),dd(:),hh(:), &
                             mi(:),ss(:)
!
      atlas=.false.
      lname='DATA/load_file'
      call rd_inp(modname,lltname,zuv,c_id,ncon,APRI,geo, &
                  outname,interp_micon)
      open(unit=11,file=outname,status='unknown')
      write(*,*)trim(modname)
      write(11,*)trim(modname)
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
      narg=iargc()
      if(narg.gt.0)then
!       read-in time file if provided
        call getarg(1,arg)
        read(arg(3:80),'(a80)')tfname
        ntime=0
        open(unit=1,file=tfname,status='old',err=16)
18      read(1,*,end=17)idum,idum,idum,idum,idum,idum
        ntime=ntime+1
        go to 18
17      rewind(1)
        ntime=ntime-1
        allocate(yyyy(ntime),mm(ntime),dd(ntime), &
               hh(ntime),mi(ntime),ss(ntime),time_mjd(ntime))
        call read_time(ntime,1,2,yyyy,mm,dd,hh,mi,ss,time_mjd)
        close(1)
      else
! default usage: read times from lat_lon_time
       ntime=0
       open(unit=1,file=lltname,status='old',err=6)
8      read(1,*,end=7)dum,dum,idum,idum,idum,idum,idum,idum
       ntime=ntime+1
       go to 8
7      rewind(1)
       allocate(yyyy(ntime),mm(ntime),dd(ntime), &
                hh(ntime),mi(ntime),ss(ntime),time_mjd(ntime))
       call read_time(ntime,1,1,yyyy,mm,dd,hh,mi,ss,time_mjd)
       close(1)
      endif     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
        write(*,*)'Grid is in km, area identified:'
        select case (gid)
          case(1)
           write(*,*)'Arctic (xy_ll,ll_xy)'
          case(2)
           write(*,*)'Antarctic (xy_ll_S,ll_xy_S)'
          case(3)
           write(*,*)'Amery Ice Shelf (mapxy, mapll, SLON=+70)'
          case(4)
           write(*,*)'Weddel/Ross Sea or CATs (mapxy, mapll, SLON=-70)'
          case(5)
           write(*,*)'Mertz area of Antarctic (mapxy, mapll, SLON=150)'
         endselect
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
         write(*,*)'Predict elevations (m)'
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
         write(*,*)'Predict transport (m^2/s) and currents (cm/s)'
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
      write(11,*)'Constituents included: ',c_id(1:ncon)
!
      allocate(cind(ncon),ccind(ncon))
      call def_con_ind(c_id,ncon,c_id_mod,nc,cind)
! check this (ccind needed for ptide!)
      do ic=1,ncon
       if(zuv.eq.'z')then
        ccind(ic)=Mod_z0%cid%idc(cind(ic))
       else
        ccind(ic)=Mod_uv0%cid%idc(cind(ic))
       endif
      enddo
!
! read lats, lons from lat_lon_time file (times are already read)
      ndat=0;nlines=0;
      latp=0;lonp=0;
      open(unit=1,file=lltname,status='old',err=1)
3     read(1,*,end=2)lat,lon
      nlines=nlines+1
      if(lat.ne.latp.or.lon.ne.lonp)then
        ndat=ndat+1
        latp=lat;lonp=lon
      endif
      go to 3
2     rewind(1)
      allocate(dloc(ndat))
      allocate(mid(ndat),id(ndat),tid(nlines))
      mid=1
      latp=0;lonp=0;k=0;nmodi=0
      do l=1,nlines
       read(1,*)lat,lon
       if(latp.ne.lat.or.lonp.ne.lon)then
        k=k+1
        dloc(k)%lat(1)=lat;dloc(k)%lon(1)=lon
        dloc(k)%theta=0;dloc(k)%phi=0;dloc(k)%land=.false.
        if(atlas)then
          call def_model(lat,lon,th_lim,&
                 ph_lim,n,m,pmask,mid(k))
          if(mid(k).gt.0)then
           do i=1,k-1
            if(mid(k).eq.mid(i))go to 9
           enddo
           nmodi=nmodi+1
           id(nmodi)=mid(k)
9          continue
         endif
        endif
        x=0;y=0;dlat=lat;dlon=lon;
        if(ll_km)call convert_ll_xy(gid,dlat,dlon,x,y)
        dloc(k)%x(1)=x;dloc(k)%y(1)=y
        latp=lat;lonp=lon
        !write(*,*)dloc(k)%lat(1),dloc(k)%lon(1)
       endif
       tid(l)=k ! index showing what data time belongs (default usage)
      enddo ! l
      close(1)
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
      if(nmodi.gt.0)then
       allocate(grid(nmodi))
       do k=1,nmodi
        write(*,*)'Loading local solution grid ',trim(mnames(id(k)))
        call getDims_l(gname,id(k),grid(k),nz,1)
        if(grid(k)%n.eq.0)then ! Model grid not found
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
      ctmp(1,1)='    Lat       Lon        mm.dd.yyyy hh:mm:ss'
      write(11,*)'' 
      if(zuv.eq.'z')then
        write(11,*)trim(ctmp(1,1)),'     z(m)   Depth(m)'
      else
        write(11,*)trim(ctmp(1,1)),'   U(m^2/s)  V(m^2/s)', &
                              '   u(cm/s)   v(cm/s) Depth(m)'
      endif
      write(11,*)''
! Interpolate model on locations
      allocate(z1(nc,ndat),r(ncon,ndat), &
               u1(nc,ndat),v1(nc,ndat),ut1(nc,ndat),vt1(nc,ndat),&
               ru(ncon,ndat),rv(ncon,ndat),rut(ncon,ndat),rvt(ncon,ndat),d(ndat))
      call interp_driver_d(grid0,dloc,ndat,d,'y')
!
      if(atlas)then
       allocate(d2(ndat))
       do k=1,nmodi
         do l=1,ndat
          dloc(l)%land=.false.
         enddo
         call interp_driver_d(grid(k),dloc,ndat,d2,'y')
         do l=1,ndat
          if(d2(l).ne.0)d(l)=d2(l)
         enddo
       enddo
       deallocate(d2)
      endif
!
      select case(zuv)
      case ('z')
       z1=0
       call interp_driver_z(Mod_z0,dloc,ndat,z1,'y')
       do k=1,nmodi
        if(id(k).gt.0)then
         allocate(z2(Mod_z(k)%nc,ndat),cind2(nc))
         do l=1,ndat
          dloc(l)%land=.false.
         enddo
         call interp_driver_z(Mod_z(k),dloc,ndat,z2,'y')
         call def_con_ind(c_id_mod,nc,Mod_z(k)%cid%cons,Mod_z(k)%nc,cind2)
         do l=1,ndat
          do ic=1,nc
           if(.not.dloc(l)%land)then
            if(cind2(ic).gt.0)then
             !write(*,*)z1(ic,l)
             if(z2(cind2(ic),l).ne.0)z1(ic,l)=z2(cind2(ic),l)
             !write(*,*)z1(ic,l)
            endif
           else
            if(z1(ic,l).ne.0)dloc(l)%land=.false.
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
        call interp_driver_z(Mod_load,dloc,ndat,zl1,'n')
        do ic=1,ncon
         if(lcind(ic).ne.0)r(ic,:)=r(ic,:)+zl1(lcind(ic),:)
        enddo
       endif
      case ('u','v','U','V')       
       dloc(:)%theta=1;dloc(:)%phi=0;dloc%land=.false.
       call interp_driver_uv(Mod_uv0,dloc,grid0,ndat,u1,'y')
       call interp_driver_uv(Mod_uv0,dloc,grid0,ndat,ut1,'T')
       dloc(:)%theta=0;dloc(:)%phi=1;dloc%land=.false.
       call interp_driver_uv(Mod_uv0,dloc,grid0,ndat,v1,'y')
       call interp_driver_uv(Mod_uv0,dloc,grid0,ndat,vt1,'T')
       u1=u1*100;v1=v1*100
       do k=1,nmodi
        l=grid(k)%n*grid(k)%m
        if(id(k).gt.0.and.l.gt.1)then
         allocate(u2(Mod_uv(k)%nc,ndat),v2(Mod_uv(k)%nc,ndat),&
                  ut2(Mod_uv(k)%nc,ndat),vt2(Mod_uv(k)%nc,ndat),cind2(nc))
         do l=1,ndat
          dloc(l)%land=.false.
         enddo
         dloc(:)%theta=1;dloc(:)%phi=0
         call interp_driver_uv(Mod_uv(k),dloc,grid(k),ndat,u2,'y')
         call interp_driver_uv(Mod_uv(k),dloc,grid(k),ndat,ut2,'T')
         dloc(:)%theta=0;dloc(:)%phi=1;
         call interp_driver_uv(Mod_uv(k),dloc,grid(k),ndat,v2,'y')
         call interp_driver_uv(Mod_uv(k),dloc,grid(k),ndat,vt2,'T')
         u2=u2*100;v2=v2*100
         call def_con_ind(c_id_mod,nc,Mod_uv(k)%cid%cons,Mod_uv(k)%nc,cind2)
         do l=1,ndat
          do ic=1,nc
           if(.not.dloc(l)%land)then
            if(cind2(ic).gt.0)then
             if(u2(cind2(ic),l).ne.0)u1(ic,l)=u2(cind2(ic),l)
             if(v2(cind2(ic),l).ne.0)v1(ic,l)=v2(cind2(ic),l)
             if(ut2(cind2(ic),l).ne.0)ut1(ic,l)=ut2(cind2(ic),l)
             if(vt2(cind2(ic),l).ne.0)vt1(ic,l)=vt2(cind2(ic),l)
            endif
           else
            if(u1(ic,l).ne.0.and.v1(ic,l).ne.0)dloc(l)%land=.false.
           endif
          enddo
         enddo
         deallocate(u2,v2,ut2,vt2,cind2)
        endif        
       enddo
       do ic=1,ncon
        if(cind(ic).ne.0)then
          ru(ic,:)=u1(cind(ic),:)
          rv(ic,:)=v1(cind(ic),:)
          rut(ic,:)=ut1(cind(ic),:)
          rvt(ic,:)=vt1(cind(ic),:)
        endif
       enddo
      endselect
! predict tide
      nz=max(nlines,ntime);it=0;midp=-1
      allocate(zpred(nz),upred(nz),vpred(nz),t1(nz))
      do k=1,ndat
       if(.not.dloc(k)%land)then
        lat=dloc(k)%lat(1);lon=dloc(k)%lon(1)
        if(narg.eq.0)then ! default usage
          if(atlas)then
           if(mid(k).ne.midp)then
            midp=mid(k)
            if(mid(k).gt.0)then
             write(11,*)'Predicted from local solution :',mnames(mid(k))
            else
             write(11,*)'Predicted from global solution'
            endif
           endif
          endif
          j=0;t1=0;it=0
          do l=1,nlines
           if(tid(l).eq.k)then
            it=l-1
            go to 11
           endif
          enddo
11        do l=1,nlines
           if(tid(l).eq.k)then
            j=j+1;t1(j)=time_mjd(l)
           endif
          enddo
          nt=j;
          if(zuv.eq.'z')then          
           call ptide(r(:,k),c_id,ncon,ccind,lat,t1,nt, &
                     interp_micon,zpred)
          else
           call ptide(rut(:,k),c_id,ncon,ccind,lat,t1,nt, &
                     interp_micon,upred)
           call ptide(rvt(:,k),c_id,ncon,ccind,lat,t1,nt, &
                     interp_micon,vpred)
          endif
!
          do l=1,nt
           j=it+l
           write(cdate,'(i2,a1,i2,a1,i2)')hh(j),':',mi(j),':',ss(j)
           ctime=deblank(cdate)
           write(cdate,'(i2,a1,i2,a1,i4)')mm(j),'.',dd(j),'.',yyyy(j)
           cdate=deblank(cdate)
           if(zuv.eq.'z')then
            write(11,'(1x,f10.4,f10.4,5x,a10,1x,a8,f10.3,f10.3)') &
                 lat,lon,cdate,ctime,zpred(l),d(k)
           else
              write(11,'(1x,f10.4,f10.4,5x,a10,1x,a8,5(f10.3))') &
                  lat,lon,cdate,ctime,upred(l),vpred(l), &
                  upred(l)/d(k)*100,vpred(l)/d(k)*100,d(k)
           endif
          enddo
! OB usage March 2011
         else
           if(zuv.eq.'z')then
            call ptide(r(:,k),c_id,ncon,ccind,lat,time_mjd,ntime, &
                      interp_micon,zpred)
           else
            call ptide(rut(:,k),c_id,ncon,ccind,lat,time_mjd,ntime, &
                      interp_micon,upred)
            call ptide(rvt(:,k),c_id,ncon,ccind,lat,time_mjd,ntime, &
                     interp_micon,vpred)
           endif
!
           do l=1,ntime
            write(cdate,'(i2,a1,i2,a1,i2)')hh(l),':',mi(l),':',ss(l)
            ctime=deblank(cdate)
            write(cdate,'(i2,a1,i2,a1,i4)')mm(l),'.',dd(l),'.',yyyy(l)
            cdate=deblank(cdate)
            if(l.eq.1)then
             if(atlas)then
              if(mid(k).gt.0)then
               write(11,*)'Predicted from local solution :',mnames(mid(k))
              else
               write(11,*)'Predicted from global solution'
              endif
             endif
             write(11,'(1x,f10.4,f10.4)')lat,lon
            endif
            if(zuv.eq.'z')then
             write(11,'(26x,a10,1x,a8,f10.3,f10.3)') &
                  cdate,ctime,zpred(l),d(k)
            else
             write(11,'(26x,a10,1x,a8,5(f10.3))') &
                  cdate,ctime,upred(l),vpred(l), &
             upred(l)/d(k)*100,vpred(l)/d(k)*100,d(k)
            endif
           enddo
         endif
        else
          lat=dloc(k)%lat(1);lon=dloc(k)%lon(1)
          write(11,'(1x,f10.4,f10.4,a)')lat,lon, &
       '***** Site is out of model grid OR land *****'
        endif  
      enddo
      deallocate(dloc,grid0)
      if(geo)deallocate(Mod_load,zl1)
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
16    write(*,*)'File ',trim(tfname),' not found'
      write(*,*)'Check spelling in command line'
      call usage()
      stop
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine read_time(ntime,iunit,iopt,yyyy,mm,dd,hh, &
                           mi,ss,time_mjd)
      implicit none
      integer k,iunit,iopt,julian,mjd,mm1,dd1,yyyy1
      integer ntime,yyyy(ntime),mm(ntime),dd(ntime)
      integer hh(ntime),mi(ntime),ss(ntime)
      real dum
      real*8 time_mjd(ntime)
      character*10 cdate,deblank
!
      do k=1,ntime
       if(iopt.eq.1)then
        read(iunit,*)dum,dum,yyyy(k),mm(k),dd(k),hh(k),mi(k),ss(k)
       else
        read(iunit,*)yyyy(k),mm(k),dd(k),hh(k),mi(k),ss(k)
       endif
! convert to mjd
       call date_mjd(mm(k),dd(k),yyyy(k),mjd)
! check if exists such a date
        julian=mjd+2400001        
        call CALDAT (JULIAN,MM1,DD1,YYYY1)
        if(mm(k).ne.mm1.or.dd(k).ne.dd1.or.yyyy(k).ne.yyyy1)then
         write(cdate,'(i2,a1,i2,a1,i4)')mm(k),'.',dd(k),'.',yyyy(k)
         cdate=deblank(cdate)
         write(*,*)'Wrong date in (lat_lon)_time file:',cdate
         stop
        endif
        time_mjd(k)=dble(mjd)+dble(hh(k))/24.D0+ &
                    dble(mi(k))/(24.D0*60.D0)  + &
                    dble(ss(k))/(24.D0*60.D0*60.D0)
      enddo
      return
      end

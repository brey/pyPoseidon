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
      program extract_local_model
!cc   Extracts local model at 1/30 resolution for a given rectangular area 
!cc   from tpxo8-atlas-compact files (elevations, transports and bathymetry)
!cc   and  outputs the model in standard OTIS binary format
! 
      implicit none
      include 'derived_types.h'
      include 'constit.h'
      type (c_grid), pointer:: grid0,grid(:)
      type (constituents) cid,cid0
      type (z_model), pointer:: Mod_z0,Mod_z(:)
      type (uv_model), pointer:: Mod_uv0,Mod_uv(:)
      type (LocDict), pointer:: dloc(:)
      integer*4 nob,nlines,nz,nu,nv,iz
      integer*4 i1,i2,j1,j2
      logical atlas,ll_km
!
      complex*8, allocatable:: z1(:,:),u1(:,:),z0(:,:),u0(:,:),v0(:,:)
      complex*8, allocatable:: z2(:,:),u2(:,:),v1(:,:),v2(:,:)
      real*4 lon,lat,x,y,latp,lonp,lat1,lat2,lon1,lon2
      real*8 dlat,dlon,qq
      real*4 th_lim(2),ph_lim(2),lth_lim(2),lph_lim(2)
      integer*4, allocatable:: cind(:),lcind(:),pmask(:,:),mid(:),id(:),mid1(:)
      integer*4, allocatable:: cind2(:),nodeID(:)
      integer(kind=4), allocatable:: mu(:,:),mv(:,:)
      real*4, allocatable:: depth(:),d2(:)
!
      character*4 c_id(ncmx),c_id_mod(ncmx),lc_id(ncmx),c_id_mod1(ncmx)
      character*80 modname,outname(4),ctmp(1,3),lname,rmCom
      character*80 hname,uname,gname
      character*20, allocatable:: mnames(:)
      character*1 zuv(3),c1
      integer*4 ncon,nc,n,m,ndat,i,k,l,nmod,nmodi,midp,ln
      integer(kind=4) ic,ncl,nl,ml,gid,n1,m1,nc1,nob1,ii,jj
      data zuv/'z','U','V'/
!
      atlas=.false.
      call rd_lmsetup(modname,lth_lim,lph_lim,c_id,ncon,outname)
      write(*,*)trim(modname)
      if(ncon.gt.0)write(*,*)'Constituents to include: ',c_id(1:ncon)
! read Model_* file
       open(unit=1,file=modname,status='old')
       read(1,'(a80)')hname
       read(1,'(a80)')uname
       read(1,'(a80)')gname
       hname=rmCom(hname)
       uname=rmCom(uname)
       gname=rmCom(gname)
       close(1)
! get model dimensions and limits
       allocate(grid0)
       call getDims(gname,grid0)
       th_lim=grid0%theta_lim;ph_lim=grid0%phi_lim;
       n=grid0%n;m=grid0%m;
       nob=max(1,grid0%nob);ll_km=grid0%ll_km;
       gid=grid0%km_area_id
       allocate(grid0%hz(n,m),grid0%mz(n,m),grid0%iob(2,nob))
       grid0%allocated=.true.
       call grid_in(gname,grid0,.false.)
       allocate(mu(n,m),mv(n,m))
       call Muv(grid0%mz,mu,mv,n,m)
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
      write(*,*)'Lat limits:   ',th_lim
      write(*,*)'Lon limits:   ',ph_lim
! Find nc
      allocate(Mod_z0)
      call getDims_model(hname,Mod_z0)
      nc=Mod_z0%nc
! Global elevations
      write(*,'(a,$)')' Reading global elevations...'
      allocate(Mod_z0%z(nc,n,m),Mod_z0%mz(n,m),Mod_z0%cid)
      allocate(Mod_z0%cid%cons(nc),Mod_z0%cid%idc(nc),Mod_z0%cid%group(nc))
      call loadModel_z(Mod_z0,hname)
      call def_cid(nc,Mod_z0%cid) ! result in Mod_z0%cid%idc
      c_id_mod(1:nc)=Mod_z0%cid%cons
      write(*,*)'done'
! Global transports
      write(*,'(a,$)')' Reading global transports...'
      allocate(Mod_uv0)
      allocate(Mod_uv0%u(nc,n,m),Mod_uv0%v(nc,n,m), &
               Mod_uv0%mu(n,m),Mod_uv0%mv(n,m),Mod_uv0%cid) 
      allocate(Mod_uv0%cid%cons(nc),Mod_uv0%cid%idc(nc),Mod_uv0%cid%group(nc))
      call loadModel_uv(Mod_uv0,uname)
      write(*,*)'done'

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

! Define set of lats/lons to extract from tpxo8-atlas-compact
! align for 1/30 degree
      lat1=int(lth_lim(1)*30)/30.
      lat2=int(lth_lim(2)*30+0.5)/30.
      lon1=int(lph_lim(1)*30)/30.
      lon2=int(lph_lim(2)*30+0.5)/30.
! correct limits for C-grid
      lth_lim(1)=lat1-1./60.
      lth_lim(2)=lat2+1./60.
      lph_lim(1)=lon1-1./60.
      lph_lim(2)=lon2+1./60.
!
      nl=(lph_lim(2)-lph_lim(1))*30
      ml=(lth_lim(2)-lth_lim(1))*30       
      ndat=nl*ml
!
      write(*,'(a,$)')' Defining local grid...'
      allocate(dloc(ndat))
      allocate(mid(ndat),id(ndat),mid1(ndat),nodeID(ndat))
! zero dloc
      dloc(:)%theta=0;dloc(:)%phi=0;dloc(:)%land=.false.
      do k=1,4
        dloc(:)%lat(k)=0;dloc(:)%lon(k)=0;
      enddo
!
      k=1
      nmodi=0
      do ii=1,nl
       lon=lon1+(ii-1.)/30.
       do jj=1,ml
         lat=lat1+(jj-1.)/30.
         dloc(k)%lat(1)=lat
         dloc(k)%lon(1)=lon
         if(atlas)then
! mid(k)<0, if coastal node on global grid
           call def_modelC(lat,lon,th_lim,ph_lim,n,m,pmask,grid0,mid(k),nodeID(k))
           if(mid(k).gt.0)then
            do i=1,k-1
             if(mid(k).eq.mid(i))go to 7
            enddo
            nmodi=nmodi+1
            id(nmodi)=mid(k)
7           continue
           endif
         endif
         k=k+1
         if(k.gt.ndat+1)then
          write(*,*)ii,nl,jj,ml,k,ndat
          write(*,*)lat,lth_lim
          write(*,*)lon,lph_lim
          stop 911
         endif
       enddo
      enddo

      write(*,*)'done'
      write(*,*)'Local models in grid:',nmodi
      !write(*,*)mnames(id(1:nmodi))      
      write(*,'(a,$)')' Interpolating global bathymetry...'
      allocate(z1(nc,ndat),u1(nc,ndat),v1(nc,ndat),depth(ndat))
      allocate(z0(nc,ndat),u0(nc,ndat),v0(nc,ndat))

!
      !depth=nodeID
      !outname(1)='DATA/grid_nodeID'
      !call write_grid(dloc,ndat,depth,outname(1),nl,ml,lth_lim,lph_lim)      
      !stop
!
      call interp_driver_d(grid0,dloc,ndat,depth,'y')
      write(*,*)'done' 
!
      write(*,'(a,$)')' Interpolating global elevations...'
      do ic=1,nc
       where(Mod_z0%z(ic,:,:).eq.0)Mod_z0%mz=0
       where(Mod_uv0%u(ic,:,:).eq.0)Mod_uv0%mu=0
       where(Mod_uv0%v(ic,:,:).eq.0)Mod_uv0%mv=0
      enddo   
      call interp_driver_z(Mod_z0,dloc,ndat,z0,'y')
      where(grid0%mz.eq.0)Mod_z0%mz(:,:)=0
      do ic=1,nc
       where(grid0%mz.eq.0)Mod_z0%z(ic,:,:)=0
      enddo
      call interp_driver_z(Mod_z0,dloc,ndat,z1,'y')
      write(*,*)'done'
      write(*,'(a,$)')' Interpolating global U transports...'
      call interp_driver_uvq(Mod_uv0,dloc,ndat,u0,'u','y') 
      where(mu.eq.0)Mod_uv0%mu(:,:)=0
      do ic=1,nc
       where(mu.eq.0)Mod_uv0%u(ic,:,:)=0
      enddo
      call interp_driver_uvq(Mod_uv0,dloc,ndat,u1,'u','y')    
      write(*,*)'done' 
      write(*,'(a,$)')' Interpolating global V transports...'
      call interp_driver_uvq(Mod_uv0,dloc,ndat,v0,'v','y')
      where(mv.eq.0)Mod_uv0%mv(:,:)=0
      do ic=1,nc
       where(mv.eq.0)Mod_uv0%v(ic,:,:)=0
      enddo
      call interp_driver_uvq(Mod_uv0,dloc,ndat,v1,'v','y')    
      write(*,*)'done'
      if(nmodi.gt.0)then
!-------------------------------------------------------------------------
       write(*,*)'Extracting local model data...'
       allocate(grid(nmodi),d2(ndat))
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
! interpolate depth
        call interp_driver_d(grid(k),dloc,ndat,d2,'y')
        do i=1,ndat
         if(mid(i).eq.id(k))then
          if(d2(i).ne.0)then
           depth(i)=d2(i)
          else
           if(nodeID(i).lt.0.and.k.lt.2)depth(i)=d2(i) ! if coastal - prefer local zeros
          endif
         endif
        enddo      
       enddo
!--------------------------------------------------------------------------------------
       !go to 111
       allocate(Mod_z(nmodi))
       do k=1,nmodi
           write(*,*)'Loading local solution elevations ',mnames(id(k))
           call getDims_z_l(hname,id(k),nz,8,Mod_z(k))
           if(Mod_z(k)%n.eq.0)then ! Model not found
            !where(abs(mid).eq.id(k))mid=0       
            id(k)=0;Mod_z(k)%n=1;Mod_z(k)%m=1;Mod_z(k)%nc=1;     
           endif
           nc1=Mod_z(k)%nc;n1=Mod_z(k)%n;m1=Mod_z(k)%m;
           allocate(Mod_z(k)%z(nc1,n1,m1),Mod_z(k)%mz(n1,m1),Mod_z(k)%cid)
           allocate(Mod_z(k)%cid%cons(nc1),Mod_z(k)%cid%idc(nc1),Mod_z(k)%cid%group(nc1))
           if(id(k).gt.0)then
            call loadModel_z_l(8,Mod_z(k),nz)
            call def_cid(nc1,Mod_z(k)%cid) ! result in Mod_z(k)%cid%idc
            write(*,*)'Constituents:',Mod_z(k)%cid%cons
           endif
!
        if(id(k).gt.0)then
         allocate(z2(Mod_z(k)%nc,ndat),cind2(nc))
         call interp_driver_z(Mod_z(k),dloc,ndat,z2,'y')
         call def_con_ind(c_id_mod,nc,Mod_z(k)%cid%cons,Mod_z(k)%nc,cind2)
!
         do l=1,ndat
          do ic=1,nc
           if(mid(l).eq.id(k))then
            if(cind2(ic).gt.0)then
             if(z2(cind2(ic),l).ne.0)then
               if(nodeID(l).le.1)then
                z1(ic,l)=z2(cind2(ic),l)
               else
                qq=(nodeID(l)-1.)/9.
                !write(*,*)z1(ic,l),z2(cind2(ic),l)
                z1(ic,l)=(1.-qq)*z1(ic,l)+qq*z2(cind2(ic),l)
                !write(*,*)z1(ic,l)
               endif
             else
                if(nodeID(l).lt.0.and.k.lt.2)z1(ic,l)=0 ! if coastal - prefer local zeros
             endif
            else !if icin2(ic)==0 - no such constituent in local solution
             z1(ic,l)=z0(ic,l)  
            endif   
           endif
          enddo
         enddo
         deallocate(z2,cind2)
        endif        
!
        enddo
!----------------------------------------------------------------------------------------
        allocate(Mod_uv(nmodi))
        do k=1,nmodi
          write(*,*)'Loading local solution transports ',mnames(id(k))
          call getDims_u_l(uname,id(k),nu,nv,8,Mod_uv(k))
          if(Mod_uv(k)%n.eq.0)then ! Model not found
            !where(abs(mid).eq.id(k))mid=0
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
!
        if(id(k).gt.0)then
         allocate(u2(Mod_uv(k)%nc,ndat),v2(Mod_uv(k)%nc,ndat),cind2(nc))
         call interp_driver_uvq(Mod_uv(k),dloc,ndat,u2,'u','y')
         call interp_driver_uvq(Mod_uv(k),dloc,ndat,v2,'v','y')
         call def_con_ind(c_id_mod,nc,Mod_uv(k)%cid%cons,Mod_uv(k)%nc,cind2)
!
         do l=1,ndat
          do ic=1,nc
           if(mid(l).eq.id(k))then
            if(cind2(ic).gt.0)then
             if(u2(cind2(ic),l).ne.0)then
               if(nodeID(l).le.1)then
                u1(ic,l)=u2(cind2(ic),l)
               else
                qq=(nodeID(l)-1.)/9.
                u1(ic,l)=(1.-qq)*u1(ic,l)+qq*u2(cind2(ic),l)
               endif
             else
                if(nodeID(l).lt.0.and.k.lt.2)u1(ic,l)=0 ! if coastal - prefer local zeros
             endif
             if(v2(cind2(ic),l).ne.0)then
               if(nodeID(l).le.1)then
                v1(ic,l)=v2(cind2(ic),l)
               else
                qq=(nodeID(l)-1.)/9.
                v1(ic,l)=(1.-qq)*v1(ic,l)+qq*v2(cind2(ic),l)
               endif
             else
                if(nodeID(l).lt.0.and.k.lt.2)v1(ic,l)=0 ! if coastal - prefer local zeros
             endif
            else !if cind2(ic)==0 - no such constituent in local solution
             u1(ic,l)=u0(ic,l)     
             v1(ic,l)=v0(ic,l)     
            endif     
           endif
          enddo
         enddo
         deallocate(u2,v2,cind2)
        endif        
!
        enddo
!---------------------------------------------------------
      endif
! writing output files
111      write(*,'(a,$)')' Writing output files...'
      call write_grid(dloc,ndat,depth,outname(1),nl,ml,lth_lim,lph_lim)
     
      call write_z(dloc,ndat,z1,depth,outname(2),nl,ml,lth_lim,lph_lim, &
                   nc,ncon,cind,c_id)
      call write_uv(dloc,ndat,u1,v1,depth,outname(3),nl,ml,lth_lim,lph_lim, &
                   nc,ncon,cind,c_id)
      open(unit=1,file=outname(4),status='unknown')
      write(1,*)outname(2)
      write(1,*)outname(3)
      write(1,*)outname(1)
      close(1)      
      write(*,*)'done'
      stop
!
      deallocate(z1,u1,v1,depth,dloc,grid0)
      deallocate(Mod_z0)
      deallocate(Mod_uv0)
      stop
      end

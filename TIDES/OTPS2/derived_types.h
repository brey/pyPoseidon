	type datavec
		integer:: Ndat, nc, itype              ! itype=1,2,3,4 for z,u,v,radial velocity
                real(kind=8), pointer:: data(:),time(:),err(:) ! time series
                character*6 type_tag                   ! just nice to have
                integer(kind=4):: ind_loc              ! points to data location in location dictionary
		logical:: ha, subprior, apsal, correct ! Are data HA? Can prior be subtracted,
	        logical:: xodif,itide                  ! sal applied, correction model
                                                       ! be subtracted from the data?
                                                       ! if this x/o difference or internal tide "combined" HC?
	end type datavec
!
        type sparsedf                                     ! sparse data functional on C-grid
		integer(kind=4):: nloc,itype              ! itype=1,2,3,4 as above
		integer(kind=4), pointer:: iC(:,:),jC(:,:)! C-grid i,j indices, surrounding
                                                          ! measurement points (allocated (nloc,4))
                real (kind=8),pointer:: ww(:,:)           ! weights for bi-linear spline 
                                                          ! interpolation, ordered as iC, jC (allocated (nloc,4))
                integer (kind=4):: ind_loc                ! points to data location in location dictionary
                type (c_grid),  pointer::   grid          ! pointer to the parental c-grid
                logical:: allocated
        end type sparsedf
!
! theta=1, phi=0 for u
! theta=0, phi=1 for v
! 0<theta, phi<2*pi for radial velocity
! theta=0, phi=0 for elevation
! 
        type LocDict
                real(kind=8) lat(4),lon(4),theta,phi
                real(kind=8) x(4),y(4)              ! coordinates in km corresponding to lat,lon
                logical land
        end type LocDict
!
        type c_grid
               integer(kind=4)::n,m,nob,km_area_id        ! km_area_id - work out later:
               logical:: lglobe,ll_km,NPG,allocated       ! could be kept in dt<0, i.e. 
               real(kind=4)theta_lim(2),phi_lim(2),dt     ! if dt<0; km_area_id=int(abs(dt));endif
               real(kind=4)x_lim(2),y_lim(2)              ! if grid is in km, these are x_lim,y_lim
               integer(kind=4), pointer :: mz(:,:),iob(:,:)
               real(kind=4), pointer :: hz(:,:)
        end type c_grid
!
        type z_model
               integer(kind=4)::n,m,nc,km_area_id         ! same as in c_grid
               logical:: lglobe,ll_km,NPG                 ! same as in c_grid
               real(kind=4)theta_lim(2),phi_lim(2)        !
               real(kind=4)x_lim(2),y_lim(2)              ! if model is in km, these are xy_lims
               integer(kind=4), pointer:: mz(:,:)         !
               complex(kind=4), pointer :: z(:,:,:)       ! elevations, z(nc,n,m)
               type (constituents), pointer :: cid        ! constituents
        end type z_model
!
        type uv_model
               integer(kind=4)::n,m,nc,km_area_id         ! km_area_id - same as in c_grid
               logical:: lglobe,ll_km,NPG                 ! same as in c_grid
               real(kind=4)theta_lim(2),phi_lim(2)        !
               real(kind=4)x_lim(2),y_lim(2)              ! if model is in km, these are xy_lims
               integer(kind=4), pointer:: mu(:,:),mv(:,:) !
               complex(kind=4), pointer :: u(:,:,:)       ! transports EW: u(:,:,:)
               complex(kind=4), pointer :: v(:,:,:)       ! transports NS: v(:,:,:)
               type (constituents), pointer :: cid        ! constituents
        end type uv_model
!
	type constituents
		integer(kind=4):: nc
		character(len=4),pointer:: cons(:)
                integer(kind=4), pointer:: idc(:)   ! index in constit.h
                integer(kind=4), pointer:: group(:) ! constituent group id in 
                                                    ! ../include/constit_f90.h
	end type constituents
!
        type masks
                 integer(kind=4), pointer:: mu(:,:),mv(:,:),mz(:,:)
                 integer(kind=4), pointer:: mbu(:,:),mbv(:,:)
                 integer(kind=4), pointer:: nu(:,:),nv(:,:)  
        end type masks

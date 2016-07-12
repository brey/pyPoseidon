!cc   This file contains the standard parameters which define the
!cc   amplitudes, frequencies, etc. for the primary tidal constituents
!cc   Currently knows about 20 tidal constituents:
!cc      NOW KNOWS about 21: M4 also
      integer, parameter :: ncmx = 29
      character*4 constid(ncmx)
      data constid &
                  /'m2  ','s2  ','k1  ','o1  ', &
                   'n2  ','p1  ','k2  ','q1  ', &
                   '2n2 ','mu2 ','nu2 ','l2  ', &
                   't2  ','j1  ','m1  ','oo1 ', &
                   'rho1','mf  ','mm  ','ssa ', &
                   'm4  ','ms4 ','mn4 ','m6  ', &
                   'm8  ','mk3 ','s6  ','2sm2', &
                   '2mk3'/

!    FOR EACH POSSIBLE CONSTIUENT, these parameters are given:
!    alpha = correction factor for first order load tides
!    amp = amplitude of equilibrium tide in m
!    ph = Currently set to zero ...   phases for
!             each constituent are referred to the time
!             when the phase of the forcing for that
!             constituent is zero on the Greenich meridian.)

!    omega = angular frequency of constituent, in radians
      real alpha_d(ncmx),ph_d(ncmx),amp_d(ncmx),omega_d(ncmx) &
            ,phase_mkB(ncmx),beta_SE(ncmx)
      integer ispec_d(ncmx)

!     Tidal parameters taken from Rodney's constituent.h, 2/23/96:
!     (except for ispec).
      data ispec_d/ &
          2,2,1,1, &
          2,1,2,1, &
          2,2,2,2, &
          2,1,1,1, &
          1,0,0,0, &
          0,0,0,0, &
          0,0,0,0, &
          0/
!cc     note: for now I am just leaving ispec for M4 set to 0 (ispec
!cc     is only used to define forcing in atgf, and this is always  0
!cc     for M4)

      data alpha_d/ &
          0.693,0.693,0.736,0.695,&
          0.693,0.706,0.693,0.695,&
          0.693,0.693,0.693,0.693,&
          0.693,0.695,0.695,0.695,&
          0.695,0.693,0.693,0.693,&
          0.693,0.693,0.693,0.693,&
          0.693,0.693,0.693,0.693,&
          0.693/

      data omega_d/ &
          1.405189e-04,1.454441e-04,7.292117e-05,6.759774e-05,&
          1.378797e-04,7.252295e-05,1.458423e-04,6.495854e-05,&
          1.352405e-04,1.355937e-04,1.382329e-04,1.431581e-04,&
          1.452450e-04,7.556036e-05,7.028195e-05,7.824458e-05,&
          6.531174e-05,0.053234e-04,0.026392e-04,0.003982e-04,&
          2.810377e-04,2.859630e-04,2.783984e-04,4.215566e-04,&
          5.620755e-04,2.134402e-04,4.363323e-04,1.503693e-04,&
          2.081166e-04/

      data ph_d/29*0.0/

      data amp_d/ &
          0.242334,0.112743,0.141565,0.100661, &
          0.046397,0.046848,0.030684,0.019273, &
          0.006141,0.007408,0.008811,0.006931, &
          0.006608,0.007915,0.007915,0.004338, &
          0.003661,0.042041,0.022191,0.019567, &
!cc       amplitude for M4 etc. is zero
          0.,0.,0.,0., &
          0.,0.,0.,0., &
          0./
 
! Astronomical arguments, obtained with Richard Ray's
! "arguments" and "astrol", for Jan 1, 1992, 00:00 Greenwich time
! Corrected July 12, 2000  
       data phase_mkB/ &
          1.731557546,0.000000000,0.173003674,1.558553872, &
          6.050721243,6.110181633,3.487600001,5.877717569, &
          4.086699633,3.463115091,5.427136701,0.553986502, &
          0.052841931,2.137025284,2.436575100,1.929046130, &
          5.254133027,1.756042456,1.964021610,3.487600001, &
          3.463115091,1.731557546,1.499093481,5.194672637, &
          6.926230184,1.904561220,0.000000000,4.551627762, &
          3.809122439/
! I am putting 0 for ms2,mn4 etc. for now: correct later
! Now this correction is done using the SAL file (h_TPXO3_90-90.load)
! I replace beta_SE with units for now (on case we decide to switch back
! to old version) and comment the old numbers - this way I do NOT change
! anything in subroutines  
! This was in weights.h before - placed here not to mix with w!
! to remove solid Earth tide multily by beta:
       data beta_SE/ &
          0.9540,0.9540,0.9400,0.9400, &
          0.9540,0.9400,0.9540,0.9400, &
          0.9540,0.9540,0.9540,0.9540, &
          0.9540,0.9400,0.9400,0.9400, &
          0.9400,0.9400,0.9400,0.9400, &
!ccc      for M4 just using value for semi-diurnals (no good reason!)
          0.9540,0.9540,0.9540,0.954, &
          0.9540,0.9540,0.9540,0.954, &
          0.9540/
!       data beta_SE/29*1./


ccc   This file contains the standard parameters which define the
ccc   amplitudes, frequencies, etc. for the primary tidal constituents
ccc   Currently knows about 20 tidal constituents:
ccc      NOW KNOWS about 21: M4 also
      integer, parameter :: ncmx = 29
      character*4 constid(ncmx)
      data constid
     1            /'m2  ','s2  ','k1  ','o1  ',
     2             'n2  ','p1  ','k2  ','q1  ',
     3             '2n2 ','mu2 ','nu2 ','l2  ',
     4             't2  ','j1  ','m1  ','oo1 ',
     5             'rho1','mf  ','mm  ','ssa ',
     6             'm4  ','ms4 ','mn4 ','m6  ',
     7             'm8  ','mk3 ','s6  ','2sm2',
     8             '2mk3'/

c    FOR EACH POSSIBLE CONSTIUENT, these parameters are given:
c    alpha = correction factor for first order load tides
c    amp = amplitude of equilibrium tide in m
c    ph = Currently set to zero ...   phases for
c             each constituent are referred to the time
c             when the phase of the forcing for that
c             constituent is zero on the Greenich meridian.)

c    omega = angular frequency of constituent, in radians
      real alpha_d(ncmx),ph_d(ncmx),amp_d(ncmx),omega_d(ncmx)
     *      ,phase_mkB(ncmx),beta_SE(ncmx)
      integer ispec_d(ncmx)

c     Tidal parameters taken from Rodney's constituent.h, 2/23/96:
c     (except for ispec).
      data ispec_d/
     1    2,2,1,1,
     2    2,1,2,1,
     3    2,2,2,2,
     4    2,1,1,1,
     5    1,0,0,0,
     6    0,0,0,0,
     7    0,0,0,0,
     8    0/
ccc     note: for now I am just leaving ispec for M4 set to 0 (ispec
ccc     is only used to define forcing in atgf, and this is always  0
ccc     for M4)

      data alpha_d/
     1    0.693,0.693,0.736,0.695,
     2    0.693,0.706,0.693,0.695,
     3    0.693,0.693,0.693,0.693,
     4    0.693,0.695,0.695,0.695,
     5    0.695,0.693,0.693,0.693,
     6    0.693,0.693,0.693,0.693,
     7    0.693,0.693,0.693,0.693,
     8    0.693/

      data omega_d/
     1    1.405189e-04,1.454441e-04,7.292117e-05,6.759774e-05,
     2    1.378797e-04,7.252295e-05,1.458423e-04,6.495854e-05,
     3    1.352405e-04,1.355937e-04,1.382329e-04,1.431581e-04,
     4    1.452450e-04,7.556036e-05,7.028195e-05,7.824458e-05,
     5    6.531174e-05,0.053234e-04,0.026392e-04,0.003982e-04,
     6    2.810377e-04,2.859630e-04,2.783984e-04,4.215566e-04,
     7    5.620755e-04,2.134402e-04,4.363323e-04,1.503693e-04,
     8    2.081166e-04/

      data ph_d/29*0.0/

      data amp_d/
     1    0.242334,0.112743,0.141565,0.100661,
     2    0.046397,0.046848,0.030684,0.019273,
     3    0.006141,0.007408,0.008811,0.006931,
     4    0.006608,0.007915,0.007915,0.004338,
     5    0.003661,0.042041,0.022191,0.019567,
ccc       amplitude for M4 etc. is zero
     6    0.,0.,0.,0.,
     7    0.,0.,0.,0.,
     8    0./
 
C Astronomical arguments, obtained with Richard Ray's
c "arguments" and "astrol", for Jan 1, 1992, 00:00 Greenwich time
c Corrected July 12, 2000  
       data phase_mkB/
     1    1.731557546,0.000000000,0.173003674,1.558553872,
     2    6.050721243,6.110181633,3.487600001,5.877717569,
     3    4.086699633,3.463115091,5.427136701,0.553986502,
     4    0.052841931,2.137025284,2.436575100,1.929046130,
     5    5.254133027,1.756042456,1.964021610,3.487600001,
     6    3.463115091,1.731557546,1.499093481,5.194672637,
     7    6.926230184,1.904561220,0.000000000,4.551627762,
     8    3.809122439/
c I am putting 0 for ms2,mn4 etc. for now: correct later
c Now this correction is done using the SAL file (h_TPXO3_90-90.load)
c I replace beta_SE with units for now (on case we decide to switch back
c to old version) and comment the old numbers - this way I do NOT change
c anything in subroutines  
c This was in weights.h before - placed here not to mix with w!
c to remove solid Earth tide multily by beta:
       data beta_SE/
     1    0.9540,0.9540,0.9400,0.9400,
     2    0.9540,0.9400,0.9540,0.9400,
     3    0.9540,0.9540,0.9540,0.9540,
     4    0.9540,0.9400,0.9400,0.9400,
     5    0.9400,0.9400,0.9400,0.9400,
cccc      for M4 just using value for semi-diurnals (no good reason!)
     7    0.9540,0.9540,0.9540,0.954,
     6    0.9540,0.9540,0.9540,0.954,
     8    0.9540/
c       data beta_SE/29*1./


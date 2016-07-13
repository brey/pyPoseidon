c Gary's approach to minor constituents interpolation
c DISABLED in this package
c This is comments from old weights.h
c/* 
c * weights.h:: original file from Gary for NLP=4 
c *
c * $Id: weights.h,v 1.4 1996/02/13 22:54:26 rodney Exp $
c */
c Redone by Lana 1998.02.23 for a new FORTRAN version
c This is good only for 8 constituents in order:
c 'm2  ','s2  ','k1  ','o1  ','n2  ','p1  ','k2  ','q1  '
c 
c The same order is supported in constit.h - so
c do not need to care about the corresponding indices
        real w(17,8),beta(17)
        data w(1,:)/1.0,  .00,  .00,  .00,  .00,  .00,  .00,  .00/
        data w(2,:)/0.0, 1.00,  .00,  .00,  .00,  .00,  .00,  .00/
        data w(3,:)/0.0,  .00,  1.0,  .00,  .00,  .00,  .00,  .00/
        data w(4,:)/0.0,  .00,  .00, 1.00,  .00,  .00,  .00,  .00/
        data w(5,:)/0.0,  .00,  .00,  .00, 1.00,  .00,  .00,  .00/
        data w(6,:)/0.0,  .00,  .00,  .00,  .00, 1.00,  .00,  .00/
        data w(7,:)/0.0,  .00,  .00,  .00,  .00,  .00, 1.00,  .00/
        data w(8,:)/0.0,  .00,  .00,  .00,  .00,  .00,  .00, 1.00/
        data w(9,:)/-0.0379, .0,.00,  .00,  .30859 ,0.0, .03289,.0/
        data w(10,:)/-0.03961,.0,.00,  .00,  .34380, 0.0, .03436,.0/
        data w(11,:)/.00696,  .0,.00,  .00,  .15719, 0.0, -.00547,.0/
        data w(12,:)/.02884,  .0,.00,  .00, -.05036, 0.0,  .07424,.0/
        data w(13,:)/.00854,  .0,.00,  .00, -.01913, 0.0,  .17685,.0/
        data w(14,:)/.0,  .0, -.00571, .11234, .0, .05285, .0, -.26257/
        data w(15,:)/.0,  .0,  .00749, .07474, .0, .03904, .0, -.12959/
        data w(16,:)/.0,  .0, -.03748, .12419, .0, .05843, .0, -.29027/
        data w(17,:)/.0,  .0,  .00842, .01002, .0,-.03064, .0,  .15028/
c








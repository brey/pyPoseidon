%This script creates a bathymetry file that can be used in Dflow?fm in the
%refinement method.

load diepten16


lon16=-180:1/16:180;
lat16=-90:1/16:90;
%The new longitude will be from ?181 till 181 so the refinement is also
%good on the sides.
lon=-181:2/16:181;
lat=-90:2/16:90;

[lon,lat]=ndgrid(lon,lat);

%Interpolate the data to a vector version of a 1/16 degree grid.
d=interp2(lon16,lat16,diepten16,reshape(lon,4174577,1),reshape(lat,4174577,1));

%Create a 3?column array with respectively the lon, lat and depth values.
A=[reshape(lon,4174577,1) flipud(reshape(lat,4174577,1)) d];

for i=0:1440
%Find all coordinates with ?181<lon<?180 and make their depth values equal to
%the values of 179<lon<180.
b1=i*2897+1;
b2=(i+1)*2897-15;
e1=i*2897+8;
e2=(i+1)*2897-8;
A(b1:e1,3)=flipud(A(b2:e2,3));

%Find all coordinates with 180<lon<181 and make their depth values equal to
%the values of ?180<lon<?179.
b3=(i+1)*2897-7;
b4=i*2897+9;
e3=(i+1)*2897;
e4=i*2897+16;
A(b3:e3,3)=flipud(A(b4:e4,3));
end
%Save the values as a sample file so it can be used in Dflow?FM
save diepten.xyz -ascii A

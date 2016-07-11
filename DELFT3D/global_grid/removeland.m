function [lon,lat] = removeland(lon,lat)
 %This function removes the land from the grid.

 %Load the depth file create with readgebco.
 load diepten16

 lon16=[-180:1/16:180];
 lat16=[-90:1/16:90];
 %Interpolate the depths to the grid.
 d=interp2(lon16, lat16, diepten16, lon, lat);
 %Remove all latitudes and longitudes of the grid with height bigger than 0.
 lon(d>0)=NaN;
 lat(d>0)=NaN;
 
 end

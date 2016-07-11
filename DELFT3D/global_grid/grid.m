clear all; close all; clc;
%This script creates a standard lat?lon grid and transforms is to a grid
%with 'North? and South Pole' on place of your choice.

%Create standard grid with resolution delta.
delta=1/2;
slon = -180:delta:180;
slat = -90:delta:90;
[slon,slat] = ndgrid(slon,slat);

%Move all points to their new location with the new 'North Pole' on Canada
%and the 'South Pole' still on the South Pole but shifted.
[nlon, nlat, map]=changepole(deg2rad(80),deg2rad(-85),deg2rad(-100), ...
deg2rad(65),deg2rad(slon),deg2rad(slat));

deglon=rad2deg(nlon);
deglat=rad2deg(nlat);

%Remove all points that are on land from the grid.
[lon,lat]=removeland(deglon,deglat);

%Save the grid.
wlgrid('write','NPonNA',lon,lat,'spherical');

wlgrid('write','NPonNA1',lon,flipud(lat),'spherical');

%Determine the maximum and minimum local enlargement.
mapfactor=[min(map) max(map)];

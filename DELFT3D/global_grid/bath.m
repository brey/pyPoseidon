%This script reads the dataset from GEBCO and interpolates it to a grid
%with a resolution of 1/16 degree.
clear all; close all; clc;
tic
i=1;
%Read GEBCO part by part so computer doesn't give a memory error
while i<234
start=(i-1)*10^6+1;
count=10^6;
diepten(start:i*count)=ncread('GEBCO1/GRIDONE_1D.nc','z',start,count);
i=i+1;
end
if i==234
start=(i-1)*10^6+1;
diepten(start:233312401)=ncread('GEBCO1/GRIDONE_1D.nc','z',start,Inf);
end
toc

E=reshape(diepten,21601,10801)';

toc

x=0:1/60:360;
y=-90:1/60:90;
toc

xn=0:1/16:360;
yn=-90:1/16:90;

[X,Y] = meshgrid(x,y);
[Xn,Yn] = meshgrid(xn,yn);

clear x; clear y; clear xn; clear yn; clear diepten, clear count, clear start
%Interpolate the GEBCO data into a 1/16 degree grid
diepten16=interp2(X,Y,E,Xn,Yn);
toc

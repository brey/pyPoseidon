
function [nlon,nlat,mapfactor]=changepole(ZPlonnew,ZPlatnew,NPlonnew,NPlatnew,lon,lat)
 %This function rotates the standard grid into the new grid with the shifted
 %'North?' and 'South' Pole.

 %Calculate parameters.
 Chi_p=0.5*fc(NPlonnew-ZPlonnew,pi/2-ZPlatnew,pi/2-NPlatnew);
 alpha=ft(NPlonnew-ZPlonnew,pi/2-ZPlatnew,pi/2-NPlatnew);
 Qlon=NPlonnew-ft(alpha,Chi_p,pi/2-NPlatnew);
 Qlat=fc(alpha,pi/2-NPlatnew,Chi_p);
 Beta=ft(alpha,pi/2-NPlatnew,Chi_p);
 r_p=tan(Chi_p/2);

 %Transform the old points to the new points.

 Chi_c=fc(lon,pi/2-lat,pi/2);
 theta_c=ft(lon,pi/2-lat,pi/2);

 Chi_s=2*atan(r_p*tan(0.5*Chi_c));
 theta_s=theta_c+Beta;

 nlon=Qlon+ft(theta_s,Chi_s,Qlat);
 nlat=pi/2-fc(theta_s,Chi_s,Qlat);

 %Make sure the final grid has longitude between ?180 and 180 degrees.
 nlon(nlon<(-pi))=nlon(nlon<(-pi))+2*pi;
 nlon(nlon>pi)=nlon(nlon>pi)-2*pi;

 mapfactor=(r_p*sec(Chi_c/2).^2)./(1+r_p^2*tan(Chi_c/2).^2);

 end
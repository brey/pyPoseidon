lat1=25+04.70/60;lon1=- 80 - 13.89/60;
lat2=25+02.42/60;lon2=- 80 - 15.50/60;
lat3=24+58.26/60;lon3=- 80 - 19.06/60;
lat=[lat1,lat2,lat3];lon=[lon1,lon2,lon3]; 
%lat=[ 34.70182, 34.58622];
%lon=[-72.93232,-73.29314]+360;	
dt=1/24;
%dt=1/24/6; % days - 10 minutes 
%lat=45.8892;lon=360-123.9608;
%
%d1=datenum(2013,09,20,00,0,0);
%d2=datenum(2013,10,30,00,0,0);
d1=datenum(2012,12,12,00,0,0);
d2=datenum(2013,02,13,00,0,0);
for l=1:length(lat)
  fid=fopen(['llt_Flo' int2str(l)],'w');
  for t=d1:dt:d2
   d=datevec(t);d(6)=round(d(6));d(5)=d(5)+round(d(6)/60);d(6)=0;
  fprintf(fid,'%10.4f %10.4f %8d %4d %4d %4d %4d %4d\n',lat(l),lon(l),d);
 end
 fclose(fid)
end

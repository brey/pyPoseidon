lat=44.5567;
lon=-125.1500;
dt=1/24;
%
d1=datenum(2000,06,21,00,0,0);
d2=datenum(2000,06,30,00,0,0);
fid=fopen('lat_lon_time_Anne','w');
for t=d1:dt:d2
 fprintf(fid,'%10.4f %10.4f %8d %4d %4d %4d %4d %4d\n',lat,lon,datevec(t));
end
fclose(fid);

function [lat_new,lon_new]=find_water(Model,lat,lon,R);

lon_new=NaN*ones(size(lon));
lat_new=NaN*ones(size(lat));

% R is search radius in km

% Get model bathy grid
[long,latg,h]=tmd_get_bathy(Model);
loc=find(h==0); h(loc)=NaN; clear loc;
mask=isnan(h);
dx=mean(diff(long));
dy=mean(diff(latg));
[ny nx]=size(h);

% Find nearest ocean point for each GZ point
for i=1:length(lon);
    
    nx_search=round(R/(abs(dx*111.7*cos(pi*lon(i)/180))));
    ny_search=round(R/(dy*111.7));
    
    %disp([num2str(i) '/' num2str(length(lon))]);
    ix=round((lon(i)-min(long))/dx)+1;
    iy=round((lat(i)-min(latg))/dy)+1;
    ix1=max([ix-nx_search 1]); ix2=min([ix+nx_search nx]);
    iy1=max([iy-ny_search 1]); iy2=min([iy+ny_search ny]);
    distmin=1e34; indx=NaN; indy=NaN;
    for iy=iy1:iy2
        for ix=ix1:ix2
            dist=distance(lat(i),lon(i),latg(iy),long(ix));
            if(dist<distmin & mask(iy,ix)==0);
                indy=iy; indx=ix;
                distmin=dist;
            end
        end
    end
    %disp([ix1 ix2 iy1 iy2 nx ny nx_search ny_search])
    if(~isnan(indx));
        lat_new(i)=latg(indy);
        lon_new(i)=long(indx);
    else;
        lat_new(i)=NaN; lon_new(i)=NaN;
    end
end
return
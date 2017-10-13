% This script reads the geos-chem variables and regrids them to match the hsrl flights
% output is a file: geoschem_curtain_<date_f#>_sub.mat
% The output contains the same variable names as were loaded from the datafile.

hsrlpath = '/gpfs_backup/meskhidze_data/DISCOVER_AQ/HSRL2/';

d = dir([hsrlpath,'*.h5']);
frmt = '201\d\d\d\d\d_R\d_F\d';
nf = numel(d);

load([pwd,filesep,'lat_long_NA_HighRes_grid.mat']);
LAT = LAT(:,1);
LONG = LONG(1,:);
for i1 = 1:nf
fprintf('Working on %d of %d\n',i1,nf)
	hlat = h5read([hsrlpath,d(i1).name],'/ApplanixIMU/gps_lat');
	hlon = h5read([hsrlpath,d(i1).name],'/ApplanixIMU/gps_lon');
	halt = h5read([hsrlpath,d(i1).name],'/DataProducts/Altitude');
	htime = h5read([hsrlpath,d(i1).name],'/ApplanixIMU/gps_time');
	
	hshrtname = regexp(d(i1).name,frmt,'match');
	hshrtname = hshrtname{:};
	hdate = regexp(hshrtname,'201\d\d\d\d\d','match');
	hdate = hdate{:};
	hsrldata.lat = hlat;
	hsrldata.lon = hlon;
	hsrldata.time = htime;
	hsrldata.z = halt;

  % Load the GEOSChem data:

	load(['../cluster_variables/newemissions/vars_',hdate,'.mat'])
	
	vars = fieldnames(optics);
	
	idx = [1,2,3,4,5,8];
	nvars = numel(idx);
        
	gcdata.lat = LAT;
	gcdata.lon = LONG;
	gcdata.alts = optics.z_centers;
	gcdata.time = 1:24;
	
	for i2 = 1:nvars
	
		if i2 == 1
		  out = grid2hsrl(gcdata,hsrldata,optics.(vars{idx(i2)}));
 		  index = out.index;
		  gridded.(vars{idx(i2)}) = out.data;
		  gridded.alts = out.alts;
		  gridded.latlon = out.coords;
		else
		  temp = optics.(vars{idx(i2)});
		  gridded.(vars{idx(i2)}) = temp(index);
		end
	end

save(['../curtains/newemissions/geoschem_curtain_',hshrtname,'.mat'],'gridded','index','-v7.3')

end 
			
clear,clc,exit

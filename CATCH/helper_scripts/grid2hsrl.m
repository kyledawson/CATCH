function OUT = grid2hsrl(gcdata,hsrldata,gridthisdata)
% This function is intendend to take an input GEOS dataset and regrid it to
% a matching HSRL dataset.
% INPUTS:  gcdata is a structure with lat, lon, time vectors and Z matrix.
%          hsrldata is a structure with lat, lon, altitude, and time.
%          gridthisdata is a matrix same size as gchem data to be gridded
%               to HSRL data size.
% OUTPUTS: out - a matrix same size as hsrldata but gridded for
%                   gcdata.  
% DATE: 8-Dec-2015
% AUTHOR: Kyle Dawson
% UPDATES:
%

% Getting field names:
fngc = fieldnames(gcdata);
latidx = strncmpi(fngc,'lat',3);
lonidx = strncmpi(fngc,'lon',3);
timeidx = strncmpi(fngc,'time',4);
zidx = ~latidx & ~lonidx & ~timeidx;

gctime = gcdata.(fngc{timeidx});
gclat = gcdata.(fngc{latidx});
gclon = gcdata.(fngc{lonidx});
gcz = gcdata.(fngc{zidx});

fnh = fieldnames(hsrldata);
hlatidx = strncmpi(fnh,'lat',3);
hlonidx = strncmpi(fnh,'lon',3);
htimeidx = strncmpi(fnh,'time',4);
haltidx = ~hlatidx & ~hlonidx & ~htimeidx;

htime = hsrldata.(fnh{htimeidx});
haltc = hsrldata.(fnh{haltidx});
hlat = hsrldata.(fnh{hlatidx});
hlon = hsrldata.(fnh{hlonidx});

% Set the GC vectors to common grid:
s = size(gcz);
zvec = zeros(size(1:s(2)));
[t,~,y,x] = ndgrid(gctime,zvec,gclat,gclon);

z = gcz;  clear gcz

dz = mode(diff(haltc));
halte = haltc-0.5*dz;
halte(end+1) = halte(end)+dz;

gctypegrid = zeros(s(2),length(htime));
gcgrididx = gctypegrid;
alts = gctypegrid;
altidx = 1:s(2);

for i1 = 1:length(hlat)
    tclose = closest(gctime,ceil(htime(i1)));
    latclose = closest(gclat,hlat(i1));
    lonclose = closest(gclon,hlon(i1));
    
    gclatclose(i1,1) = gclat(latclose);
    gclonclose(i1,1) = gclon(lonclose);
    
    alts(:,i1) = z(tclose,:,latclose,lonclose);

    % TO DO: ADD THE ALT PART
    if nargin > 2
        gctypegrid(:,i1) = squeeze(gridthisdata(tclose,:,latclose,lonclose));
    end

    % Assign the index numbers
    for i2 = altidx
        gcgrididx(i2,i1) = sub2ind(size(z),tclose,i2,latclose,lonclose);
    end
       
end

if nargin > 2
    OUT.data = gctypegrid;
end
OUT.index = gcgrididx;
OUT.alts = alts;
OUT.coords = [gclatclose,gclonclose];







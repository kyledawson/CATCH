% Kyle Dawson
% Matt sent a .dat file for GChem OPs that I need to read in so we can
% calculate the AODs on an hourly basis. Here goes:
% INPUTS:  filename, string will full filepath.
%          UseColette, logical - if true, uses Colette's alpha data.
% OUTPUTS: headers:  Optical property data headers
%          data:     Optical property data
%          info:     info on the RH and dust masses and types from the
%                    file.
% EXAMPLE: [headers,data,info] = read_OpticalProps_GeosChem
%           *** Note this example runs on my machine because the default
%           looks in my directory. ***  If running on another machine, pass
%           full filename as <filename> and execute:
%           [headers,data,info] = read_OpticalProps_GeosChem(<filename>)
% UPDATES:
% 

function [headers,data,info] = read_OpticalProps_GeosChem(filename,UseColette)

p = inputParser;
addRequired(p,'filename',@(x) ischar(x))
addRequired(p,'UseColette',@(x) islogical(x))
parse(p,filename,UseColette);
r = p.Results;
filename = r.filename;
usecolette = r.UseColette;

if ~isempty(filename)
    [fpath,~,~] = fileparts(filename);
    if ~isempty(fpath)
        fullfile = filename;
    else
        fpath = '/Volumes/cloud/GEOSChem/';
        filename = 'jv_spec_aod.dat';
        fullfile = [fpath,filename];
        if ~logical(exist(fullfile,'file'));
            error('Please pass the full filename including the path!')
        end
    end
else
    fpath = '/Volumes/cloud/GEOSChem/';
    filename = 'jv_spec_aod.dat';
    fullfile = [fpath,filename];
    if ~logical(exist(fullfile,'file'));
        error('Please pass the full filename including the path!')
    end
end

delim = ' ';
frmt=repmat('%s',1,12);
Nheader = 3;
fid = fopen(fullfile,'r');
headers = textscan(fid,frmt,1,'delimiter',delim,...
    'MultipleDelimsAsOne',1,'HeaderLines',Nheader,...
    'collectoutput',1);
headers = headers{1};
fclose(fid);

idx = getIDX(fullfile,Nheader,delim);

[data,info] = getDATA(fullfile,idx);



% NOW To Parse all this info data into useful stuff:
cellcontents = @(x) isempty(x);

% GET RH:
rhvec = regexp(info,'RH=[\d]+','match');
norhvecidx = cellfun(cellcontents,rhvec);
RH = [];

% GET M = N+Ni
refrac = regexp(info,'n@550=[\d]+[-.\d]*[\S]','match');
norefrac = cellfun(cellcontents,refrac);

% GET TYPE:
temptype = regexpi(info,['Mdust|Trop sulphate|Black C|'...
        'Organic C|Sea Salt (a\w+[)]|Sea Salt (c\w+[)]'],'match');
    
% GET THE size dist params:
% Radius:
sizedist = regexp(info,'r[_g]*=[\d]*[.]+[\d]+','match');
nosizedist = cellfun(cellcontents,sizedist);
% Sigma:
sizedist2 = regexp(info,'sigma=[\d]*[.]+[\d]+','match');

% Get the Dust Masses:
masses = regexp(info,'Mdust [\d]+[.]+[\d]+','match');

for i1 = 1:length(norhvecidx)
    
    % ASSIGN RH and dust masses:
    if norhvecidx(i1)
        rh = NaN;
        masstemp = masses{i1}{:};
        M(i1,1) = str2double(masstemp(7:end));
    else
        rh = rhvec{i1}{:};
        rh = str2double(rh(end-1:end));
        M(i1,1) = NaN;
    end
    
    RH = [RH;rh];
    
    % ASSIGN TYPE:
    
    TYPE(i1,1) = temptype{i1};
    
    % ASSIGN INDEX OF REFRAC:
     
    if norefrac(i1)
        Nr(i1,1) = NaN;
        Ni(i1,1) = NaN;
    else
        n = refrac{i1}{:};
        n = str2double(n(7:end));
        Nr(i1,1) = real(n);
        Ni(i1,1) = imag(n);
    end
    
    % ASSIGN R and Sigma:
    if nosizedist(i1)
        R(i1,1) = NaN;
        sig(i1,1) = NaN;
    else
        Rtemp = sizedist{i1}{:};
        Rtemp = regexp(Rtemp,'[\d]*[.]*\d+','match');
        sigtemp = sizedist2{i1}{:};
        
        R(i1,1) = str2double(Rtemp{:});
        sig(i1,1) = str2double(sigtemp(7:end));
    end
    
end
info0 = info; % bs added
clearvars info % 

info.raw = info0;
info.Type = TYPE;
info.Dust_mass = M;
info.RH = RH;
info.N_real = Nr;
info.N_imaginary = Ni;
info.Size_Distribution = [R,sig];

% Add density from: 
% http://www.atmos.colostate.edu/~heald/docs/GEOS_Chem_optics_description.pdf

densityDUST = [2.5;2.5;2.5;2.5;2.65;2.65;2.65];
densitySULFATE = repmat(1.7,7,1);
densityBC = repmat(1.0,7,1);
densityOC = repmat(1.8,7,1);
densitySSA = repmat(2.2,7,1);
densitySSC = repmat(2.2,7,1);

info.Density = [densityDUST;densitySULFATE;densityBC;densityOC;...
    densitySSA;densitySSC];

info.Alpha = 0.75.*data(:,2)./data(:,3)./info.Density;
if usecolette
    calpha = getColetteAlpha;
    alphatemp = info.Alpha;
    alphatemp(length(densityDUST)+1:end) = calpha;
    info.Alpha = alphatemp;
end
info.Units = {'Density (g/cm-3)','Alpha (m2/g)','Qext (m2/m2)','Reff (um)',...
    'Dust Mass (g)'};

% CALLBACKS
function idx = getIDX(fullfile,Nheader,delim)
    
    fid2 = fopen(fullfile,'r');
    
    idx=textscan(fid2,'%d','HeaderLines',Nheader,...
        'CommentStyle',delim,'Delimiter',delim,...
        'CollectOutput',1);
    idx = idx{1};
    
    fclose(fid2);
    
end

function [data,info] = getDATA(fullfile,idx)
    data = []; info = data;
    fid3 = fopen(fullfile,'r');
    cont = 1;
    while cont
        tline = fgetl(fid3);
        cont = tline~=-1;
        if cont
        xx = str2double(tline(1:2));
        isinfoline = any(xx==idx);
        if isinfoline
            info = [info;{tline}];
            tline2 = fgetl(fid3);
            dtemp = str2num(tline2);
            data = [data;dtemp];
        end
        end
    end
    
    fclose(fid3);
end

function coletteAlpha = getColetteAlpha()
    % from: 
    % http://www.atmos.colostate.edu/~heald/docs/GEOS_Chem_optics_
    % description.pdf
    rhvals = [0 50 70 80 90 95 99];
    so4alpha = [2.2 7.0 10.0 13.2 20.8 34.8 inf];
    bcalpha = [8.0 8.0 8.0 9.7 11.2 11.9 16.9];
    ocalpha = [2.8 4.8 6.1 7.6 11.4 17.9 45.5];
    ssaalpha = [2.4 9.4 13.8 18.6 32.4 57.8 inf];
    sscalpha = [0.9 2.3 2.8 3.4 4.7 6.9 18.7];
    coletteAlpha = [so4alpha';bcalpha';ocalpha';...
        ssaalpha';sscalpha'];
        
end
end
% EOF
            
    
    
        